#include "../include/pjtree.h"

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

#include "../include/JetCorrector.h"
#include "../include/JetUncertainty.h"
#include "../include/phoERegression.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TTree.h"

#include <algorithm>
#include <array>
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

float jer(std::vector<float> const& csn, float pt) {
    return std::sqrt((csn[1] - csn[4]) / pt + (csn[2] - csn[5]) / (pt * pt));
}

int regulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto array_size = conf->get<int64_t>("array_size");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");
    auto mc_branches = conf->get<bool>("mc_branches");
    auto hlt_branches = conf->get<bool>("hlt_branches");
    auto apply_weights = conf->get<bool>("apply_weights");
    auto apply_residual = conf->get<bool>("apply_residual");
    auto active = conf->get<std::vector<bool>>("active");

    auto selections = conf->get<std::vector<std::string>>("selections");
    auto paths = conf->get<std::vector<std::string>>("paths");
    auto skim = conf->get<std::vector<std::string>>("skim");
    auto algo = conf->get<std::string>("algo");
    auto jecs = conf->get<std::vector<std::string>>("jecs");
    auto jeu = conf->get<std::string>("jeu");
    auto direction = conf->get<bool>("direction");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto residual = conf->get<std::string>("residual");
    auto csn = conf->get<std::vector<float>>("csn");
    auto xmls = conf->get<std::vector<std::string>>("xmls");

    auto pthat = conf->get<std::vector<int32_t>>("pthat");
    auto pthatw = conf->get<std::vector<float>>("pthatw");
    auto vzw = conf->get<std::vector<float>>("vzw");

    auto ihf = new interval(dhf);

    for (auto& v : csn) { v = v * v; }

    /* load forest */
    auto forest = new train(files);
    auto chain_evt = forest->attach("hiEvtAnalyzer/HiTree", true);
    auto chain_sel = forest->attach("skimanalysis/HltTree", true);
    auto chain_eg = forest->attach("ggHiNtuplizerGED/EventTree", true);
    auto chain_jet = forest->attach((algo + "/t").data(), !algo.empty());
    auto chain_hlt = forest->attach("hltanalysis/HltTree", hlt_branches);

    (*forest)();

    auto tevt = harvest<event>(chain_evt, mc_branches);
    auto tsel = harvest<triggers>(chain_sel, selections);
    auto tegg = harvest<eggen>(chain_eg, mc_branches);
    auto tele = harvest<electrons>(chain_eg);
    auto tpho = harvest<photons>(chain_eg);
    auto tjet = harvest<jets>(chain_jet, mc_branches, array_size);
    auto thlt = harvest<triggers>(chain_hlt, paths);

    std::array<bool, tt::ntt> flags = { tevt, tegg, tpho, tele, tjet, thlt };
    std::transform(flags.begin(), flags.end(), active.begin(), flags.begin(),
                   [](bool a, bool b) -> bool { return a && b; });

    /* setup output tree */
    TTree::SetMaxTreeSize(1000000000000LL);

    TFile* fout = new TFile(output, "recreate");
    TTree* tout = new TTree("pj", "photon-jet");
    auto tree_pj = new pjtree(tout, mc_branches, hlt_branches, flags);

    /* weights, corrections */
    auto JEC = new JetCorrector(jecs);
    auto JEU = new JetUncertainty(jeu);

    TF1* fweight = new TF1("fweight", "(gaus(0))/(gaus(3))");
    if (mc_branches && apply_weights) { fweight->SetParameters(
        vzw[0], vzw[1], vzw[2], vzw[3], vzw[4], vzw[5]); }

    TF1** fres = nullptr;
    if (apply_residual) {
        TFile* fh = new TFile(residual.data(), "read");
        auto hres = new history<TH1F>(fh, tag + "_scale_s_dhf_f_pt");

        fres = new TF1*[hres->size()];
        hres->apply([&](TH1* h, int64_t index) {
            fres[index] = h->GetFunction(
                ("f_s_dhf_f_pt_"s + std::to_string(index)).data()); });
    }

    auto regr = new phoERegression();
    if (!xmls.empty()) {
        regr->initialiseReaderEB(xmls[0]);
        regr->initialiseReaderEE(xmls[1]);
    }

    std::vector<float> regr_variables(17, 0);
    auto fill_regr_variables = [&](int64_t index) {
        regr_variables[0] = (*tree_pj->phoSCRawE)[index];
        regr_variables[1] = (*tree_pj->phoSCEta)[index];
        regr_variables[2] = (*tree_pj->phoSCPhi)[index];
        regr_variables[3] = (*tree_pj->phoSCEtaWidth)[index];
        regr_variables[4] = (*tree_pj->phoSCPhiWidth)[index];
        regr_variables[5] = (*tree_pj->phoE3x3_2012)[index];
        regr_variables[6] = (*tree_pj->phoMaxEnergyXtal_2012)[index];
        regr_variables[7] = (*tree_pj->phoE2nd_2012)[index];
        regr_variables[8] = (*tree_pj->phoELeft_2012)[index];
        regr_variables[9] = (*tree_pj->phoERight_2012)[index];
        regr_variables[10] = (*tree_pj->phoETop_2012)[index];
        regr_variables[11] = (*tree_pj->phoEBottom_2012)[index];
        regr_variables[12] = (*tree_pj->phoSigmaIEtaIEta_2012)[index];
        regr_variables[13] = (*tree_pj->phoSigmaIEtaIPhi_2012)[index];
        regr_variables[14] = (*tree_pj->phoSigmaIPhiIPhi_2012)[index];
        regr_variables[15] = tree_pj->rho;
        regr_variables[16] = (*tree_pj->phoESEn)[index];
    };

    auto rng = new TRandom3(144);

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
                if (tsel->accept(path) == 1)
                    pass_selection = true;

            if (!pass_selection) { continue; }
        }

        if (!skim.empty()) {
            if (tpho->nPho < 1) { continue; }

            bool pass_skim = false;
            for (auto const& path : skim)
                if (thlt->accept(path) == 1)
                    pass_skim = true;

            if (!pass_skim) { continue; }
        }

        tree_pj->copy(tevt, tegg, tpho, tele, tjet, thlt);

        if (!heavyion) {
            tree_pj->hiBin = 0;
            tree_pj->hiHF = 0;
            tree_pj->Ncoll = 1000;
        }

        tree_pj->weight = (mc_branches && apply_weights)
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
            corr = corr / (apply_residual ? fres[hf_x]->Eval(corr) : 1.f);

            if (!jeu.empty()) {
                auto unc = JEU->GetUncertainty();
                corr *= direction ? (1. + unc.second) : (1. - unc.first);
            }

            if (!csn.empty()) { corr *= rng->Gaus(1., jer(csn, corr)); }

            (*tree_pj->jtpt)[j] = corr;
        }

        /* apply photon energy corrections */
        if (!xmls.empty()) {
            for (int64_t j = 0; j < tree_pj->nPho; ++j) {
                fill_regr_variables(j);
                (*tree_pj->phoEt)[j] = regr->getCorrectedPt(
                    regr_variables,
                    (*tree_pj->phoEt)[j],
                    (*tree_pj->phoEta)[j],
                    (*tree_pj->phoSCEta)[j]);
            }
        }

        tout->Fill();
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return regulate(argv[1], argv[2]);

    return 0;
}
