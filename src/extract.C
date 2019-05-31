#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include <string>
#include <vector>

#include "../include/defines.h"
#include "../include/jettree.h"
#include "../include/photontree.h"
#include "../include/pjtree.h"
#include "../include/tracktree.h"

#include "../git/config/include/configurer.h"

static constexpr float mindr2 = 0.15 * 0.15;

int extract(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto mc_branches = conf->get<bool>("mc_branches");
    auto hlt_branches = conf->get<bool>("hlt_branches");

    TChain* chain_eg = new TChain("ggHiNtuplizerGED/EventTree");
    TChain* chain_jet = new TChain("ak4PFJetAnalyzer/t");
    TChain* chain_trk = new TChain("ppTrack/trackTree");
    TChain* chain_hlt = hlt_branches ? new TChain("hltanalysis/HltTree") : 0;

    chain_eg->SetBranchStatus("*", 0);
    chain_jet->SetBranchStatus("*", 0);
    chain_trk->SetBranchStatus("*", 0);
    if (hlt_branches)
        chain_hlt->SetBranchStatus("*", 0);

    for (auto const& f : files) {
        chain_eg->Add(f.data());
        chain_jet->Add(f.data());
        chain_trk->Add(f.data());
        if (hlt_branches)
            chain_hlt->Add(f.data());
    }

    auto tree_eg = new photontree(chain_eg, mc_branches);
    auto tree_jet = new jettree(chain_jet, mc_branches);
    auto tree_trk = new tracktree(chain_trk);

    TFile* fout = new TFile(output, "recreate");
    TTree* tout = new TTree("pj", "photon-jet");

    auto tree_pj = new pjtree(tout, mc_branches);

    int64_t nentries = chain_eg->GetEntries();
    if (max_entries) nentries = std::min(nentries, max_entries);

    printf("entries: %li\n", nentries);

    for (int64_t i = 0; i < nentries; ++i) {
        tree_pj->clear();

        chain_eg->GetEntry(i);
        chain_jet->GetEntry(i);
        chain_trk->GetEntry(i);
        if (hlt_branches)
            chain_hlt->GetEntry(i);

        if (i % 10000 == 0)
            printf("entry: %li\n", i);

        tree_pj->copy(tree_eg);
        tree_pj->copy(tree_jet);
        tree_pj->copy(tree_trk);

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
