#include "../include/pjtree.h"
#include "../include/JetCorrector.h"
#include "../include/JetUncertainty.h"

#include "../git/config/include/configurer.h"

#include "../git/foliage/include/eggen.h"
#include "../git/foliage/include/electrons.h"
#include "../git/foliage/include/event.h"
#include "../git/foliage/include/jets.h"
#include "../git/foliage/include/photons.h"
#include "../git/foliage/include/triggers.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/history.h"

#include "../git/tricks-and-treats/include/train.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;

float weight_for(std::vector<int32_t> const& divisions,
                 std::vector<float> const& weights, float value) {
    int64_t index = -1;
    for (auto edge : divisions)
        if (value > edge)
            ++index;

    return weights[index];
}

int extract(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto array_size = conf->get<int64_t>("array_size");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");
    auto mc_branches = conf->get<bool>("mc_branches");
    auto hlt_branches = conf->get<bool>("hlt_branches");
    auto apply_residual = conf->get<bool>("apply_residual");

    auto selections = conf->get<std::vector<std::string>>("selections");
    auto paths = conf->get<std::vector<std::string>>("paths");
    auto skim = conf->get<std::vector<std::string>>("skim");
    auto algo = conf->get<std::string>("algo");
    auto jecs = conf->get<std::vector<std::string>>("jecs");
    auto jeu = conf->get<std::string>("jeu");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto residual = conf->get<std::string>("residual");

    auto pthat = conf->get<std::vector<int32_t>>("pthat");
    auto pthatw = conf->get<std::vector<float>>("pthatw");
    auto vzw = conf->get<std::vector<float>>("vzw");

    auto ihf = std::make_shared<interval>(dhf);

    /* load forest */
    auto forest = new train(files);
    auto chain_evt = forest->attach("hiEvtAnalyzer/HiTree", true);
    auto chain_sel = forest->attach("skimanalysis/HltTree", true);
    auto chain_eg = forest->attach("ggHiNtuplizerGED/EventTree", true);
    auto chain_jet = forest->attach((algo + "/t").data(), !algo.empty());
    auto chain_hlt = forest->attach("hltanalysis/HltTree", hlt_branches);

    (*forest)();

    auto tree_evt = harvest<event>(chain_evt, mc_branches);
    auto tree_sel = harvest<triggers>(chain_sel, selections);
    auto tree_egg = harvest<eggen>(chain_eg, mc_branches);
    auto tree_ele = harvest<electrons>(chain_eg);
    auto tree_pho = harvest<photons>(chain_eg);
    auto tree_jet = harvest<jets>(chain_jet, mc_branches, array_size);
    auto tree_hlt = harvest<triggers>(chain_hlt, paths);

    /* setup output tree */
    TTree::SetMaxTreeSize(1000000000000LL);

    TFile* fout = new TFile(output, "recreate");
    TTree* tout = new TTree("pj", "photon-jet");
    auto tree_pj = new pjtree(tout, mc_branches, hlt_branches);

    /* weights, corrections */
    auto JEC = new JetCorrector(jecs);
    auto JEU = new JetUncertainty(jeu);

    TF1* fweight = new TF1("fweight", "(gaus(0))/(gaus(3))");
    if (mc_branches) { fweight->SetParameters(
        vzw[0], vzw[1], vzw[2], vzw[3], vzw[4], vzw[5]); }

    TF1** fres = nullptr;
    if (apply_residual) {
        TFile* fh = new TFile(residual.data(), "read");
        auto hres = new history(fh, tag + "_es_dhf_f_pt");

        fres = new TF1*[hres->size()];
        hres->apply([&](TH1* h, int64_t index) {
            fres[index] = h->GetFunction(
                ("f_es_dhf_f_pt_"s + std::to_string(index)).data()); });
    }

    /* boom! */
    int64_t nentries = forest->count();
    if (max_entries) nentries = std::min(nentries, max_entries);
    for (int64_t i = 0; i < nentries; ++i) {
        tree_pj->clear();

        if (i % 10000 == 0)
            printf("entry: %li/%li\n", i, nentries);

        forest->get(i);

        if (!selections.empty()) {
            bool pass_selection = false;
            for (auto const& path : selections)
                if (tree_sel->accept(path) == 1)
                    pass_selection = true;

            if (!pass_selection) { continue; }
        }

        if (!skim.empty()) {
            if (tree_pho->nPho < 1) { continue; }

            bool pass_skim = false;
            for (auto const& path : skim)
                if (tree_hlt->accept(path) == 1)
                    pass_skim = true;

            if (!pass_skim) { continue; }
        }

        tree_pj->copy(tree_evt);
        tree_pj->copy(tree_egg);
        tree_pj->copy(tree_ele);
        tree_pj->copy(tree_pho);
        tree_pj->copy(tree_jet);
        tree_pj->copy(tree_hlt);

        if (!heavyion) {
            tree_pj->hiBin = 0;
            tree_pj->hiHF = 0;
            tree_pj->Ncoll = 1000;
        }

        tree_pj->weight = mc_branches
            ? tree_pj->Ncoll / 1000.f
                * fweight->Eval(tree_pj->vz)
                * weight_for(pthat, pthatw, tree_pj->pthat)
            : 1.f;

        auto hf_x = ihf->index_for(tree_pj->hiHF);

        /* apply jet energy corrections and evaluate uncertainties */
        for (int64_t j = 0; j < tree_pj->nref; ++j) {
            JEC->SetJetPT((*tree_pj->rawpt)[j]);
            JEC->SetJetEta((*tree_pj->jteta)[j]);
            JEC->SetJetPhi((*tree_pj->jtphi)[j]);

            float corr = JEC->GetCorrectedPT();
            float cres = (apply_residual) ? fres[hf_x]->Eval(corr) : 1.f;
            cres = cres > 0.8 ? cres : 0.8;
            (*tree_pj->jtpt)[j] = corr / cres;

            auto unc = JEU->GetUncertainty();
            tree_pj->jtpt_up->push_back(corr / cres * (1. - unc.first));
            tree_pj->jtpt_down->push_back(corr / cres * (1. + unc.second));
        }

        tout->Fill();
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();
    
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return extract(argv[1], argv[2]);

    return 0;
}
