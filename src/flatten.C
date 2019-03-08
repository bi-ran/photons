#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

int flatten(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto mc_branches = conf->get<bool>("mc_branches");

    auto photon_pt_min = conf->get<float>("photon_pt_min");
    auto photon_eta_abs = conf->get<float>("photon_eta_abs");
    auto jet_pt_min = conf->get<float>("jet_pt_min");
    auto jet_eta_abs = conf->get<float>("jet_eta_abs");
    auto trk_pt_min = conf->get<float>("trk_pt_min");
    auto trk_eta_abs = conf->get<float>("trk_eta_abs");

    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("pyjamas");
    auto pjt = new pjtree(t, mc_branches);

    TFile* fout = new TFile(output, "recreate");

    TH1F* h_photon_pt = new TH1F("h_photon_pt", "", 100, 0, 200);
    TH1F* h_jet_pt = new TH1F("h_jet_pt", "", 100, 0, 200);
    TH1F* h_trk_pt = new TH1F("h_trk_pt", "", 100, 0, 200);

    int64_t nentries = std::min((int64_t)t->GetEntries(), max_entries);
    for (int64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        for (int j = 0; j < pjt->nPho; ++j) {
            if ((*pjt->phoEt)[j] > photon_pt_min
                    && std::abs((*pjt->phoEta)[j]) < photon_eta_abs) {
                h_photon_pt->Fill((*pjt->phoEt)[j]);
            }
        }

        for (int j = 0; j < pjt->nref; ++j) {
            if ((*pjt->jtpt)[j] > jet_pt_min
                    && std::abs((*pjt->jteta)[j]) < jet_eta_abs) {
                h_jet_pt->Fill((*pjt->jtpt)[j]);
            }
        }

        for (int j = 0; j < pjt->nTrk; ++j) {
            if ((*pjt->trkPt)[j] > trk_pt_min
                    && std::abs((*pjt->trkEta)[j]) < trk_eta_abs) {
                h_trk_pt->Fill((*pjt->trkPt)[j]);
            }
        }
    }

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return flatten(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
