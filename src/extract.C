#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include <string>
#include <vector>

#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

#include "../git/foliage/include/foliage.h"
#include "../git/foliage/include/jets.h"
#include "../git/foliage/include/photons.h"
#include "../git/foliage/include/tracks.h"

#include "../git/tricks-and-treats/include/train.h"

int extract(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto mc_branches = conf->get<bool>("mc_branches");
    auto hlt_branches = conf->get<bool>("hlt_branches");
    auto paths = conf->get<std::vector<std::string>>("paths");
    auto jet_algo = conf->get<std::string>("jet_algo");
    auto array_size = conf->get<int64_t>("array_size");

    auto forest = new train(files);
    auto chain_eg = forest->attach("ggHiNtuplizerGED/EventTree", true);
    auto chain_jet = forest->attach((jet_algo + "/t").data(), true);
    auto chain_trk = forest->attach("ppTrack/trackTree", true);
    auto chain_hlt = forest->attach("hltanalysis/HltTree", hlt_branches);

    (*forest)();

    auto tree_eg = new photons(chain_eg, mc_branches);
    auto tree_jet = new jets(chain_jet, mc_branches, array_size);
    auto tree_trk = new tracks(chain_trk, array_size);
    auto tree_hlt = new triggers(chain_hlt, paths);

    TTree::SetMaxTreeSize(1000000000000LL);

    TFile* fout = new TFile(output, "recreate");
    TTree* tout = new TTree("pj", "photon-jet");
    auto tree_pj = new pjtree(tout, mc_branches, hlt_branches);

    int64_t nentries = forest->count();
    if (max_entries) nentries = std::min(nentries, max_entries);
    for (int64_t i = 0; i < nentries; ++i) {
        tree_pj->clear();

        if (i % 10000 == 0)
            printf("entry: %li\n", i);

        forest->get(i);

        tree_pj->copy(tree_eg);
        tree_pj->copy(tree_jet);
        tree_pj->copy(tree_trk);
        tree_pj->copy(tree_hlt);

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
