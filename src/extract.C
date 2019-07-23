#include "../include/pjtree.h"
#include "../include/JetCorrector.h"

#include "../git/config/include/configurer.h"

#include "../git/foliage/include/event.h"
#include "../git/foliage/include/jets.h"
#include "../git/foliage/include/photons.h"
#include "../git/foliage/include/triggers.h"

#include "../git/tricks-and-treats/include/train.h"

#include "TFile.h"
#include "TTree.h"

#include <string>
#include <vector>

int extract(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto centrality = conf->get<bool>("centrality");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto mc_branches = conf->get<bool>("mc_branches");
    auto hlt_branches = conf->get<bool>("hlt_branches");
    auto selections = conf->get<std::vector<std::string>>("selections");
    auto paths = conf->get<std::vector<std::string>>("paths");
    auto skim = conf->get<std::vector<std::string>>("skim");
    auto jet_algo = conf->get<std::string>("jet_algo");
    auto array_size = conf->get<int64_t>("array_size");
    auto jec_file = conf->get<std::string>("jec_file");

    auto forest = new train(files);
    auto chain_evt = forest->attach("hiEvtAnalyzer/HiTree", true);
    auto chain_sel = forest->attach("skimanalysis/HltTree", true);
    auto chain_eg = forest->attach("ggHiNtuplizerGED/EventTree", true);
    auto chain_jet = forest->attach((jet_algo + "/t").data(), true);
    auto chain_hlt = forest->attach("hltanalysis/HltTree", hlt_branches);

    (*forest)();

    auto tree_evt = new event(chain_evt, mc_branches);
    auto tree_sel = new triggers(chain_sel, selections);
    auto tree_egg = new eggen(chain_eg, mc_branches);
    auto tree_ele = new electrons(chain_eg);
    auto tree_pho = new photons(chain_eg);
    auto tree_jet = new jets(chain_jet, mc_branches, array_size);
    auto tree_hlt = new triggers(chain_hlt, paths);

    TTree::SetMaxTreeSize(1000000000000LL);

    TFile* fout = new TFile(output, "recreate");
    TTree* tout = new TTree("pj", "photon-jet");
    auto tree_pj = new pjtree(tout, mc_branches, hlt_branches);

    auto JEC = new JetCorrector(jec_file);

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

        if (!centrality) {
            tree_pj->hiBin = 0;
            tree_pj->hiHF = 0;
            tree_pj->Ncoll = 1000;
        }

        tree_pj->weight = mc_branches ? tree_pj->Ncoll / 1000.f : 1.f;

        /* apply l2 jet energy corrections */
        for (int64_t j = 0; j < tree_pj->nref; ++j) {
            JEC->SetJetPT((*tree_pj->rawpt)[j]);
            JEC->SetJetEta((*tree_pj->jteta)[j]);

            (*tree_pj->jtpt)[j] = JEC->GetCorrectedPT();
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
