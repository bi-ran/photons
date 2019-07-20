#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <memory>
#include <string>
#include <vector>

using namespace std::literals::string_literals;

int invigilate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto tag = conf->get<std::string>("tag");

    auto pt_min = conf->get<float>("pt_min");
    auto eta_max = conf->get<float>("eta_max");

    auto res = conf->get<std::vector<float>>("es_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto deta = conf->get<std::vector<float>>("eta_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");

    /* exclude most peripheral events */
    auto hf_min = dhf.front();

    /* prepare histograms */
    auto ies = std::make_shared<interval>("energy scale"s,
        static_cast<int64_t>(res[0]), res[1], res[2]);

    auto ipt = std::make_shared<interval>(dpt);
    auto ieta = std::make_shared<interval>(deta);
    auto ihf = std::make_shared<interval>(dhf);

    auto mptetahf = std::make_shared<multival>(dpt, deta, dhf);

    auto scale = std::make_unique<history>("scale"s, "counts", ies, mptetahf);

    /* load input */
    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("pj");
    auto p = new pjtree(true, t);

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* fill histograms */
    auto nentries = static_cast<int64_t>(t->GetEntries());
    for (int64_t i = 0; i < nentries; ++i) {
        if (i % 100000 == 0) { printf("%li/%li\n", i, nentries); }

        t->GetEntry(i);

        if (p->hiHF <= hf_min) { continue; }

        for (int64_t j = 0; j < p->nref; ++j) {
            if ((*p->subid)[j] > 0) { continue; }

            auto reco_pt = (*p->jtpt)[j];
            if (reco_pt <= pt_min) { continue; }

            auto reco_eta = (*p->jteta)[j];
            if (std::abs(reco_eta) > eta_max) { continue; }

            auto gen_pt = (*p->refpt)[j];
            if (gen_pt < 0) { continue; }

            (*scale)[v{reco_pt, reco_eta, p->hiHF}]->Fill(
                reco_pt / gen_pt, p->weight);
        }
    }

    /* save output */
    TFile* fout = new TFile(output, "recreate");

    scale->save(tag);

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return invigilate(argv[1], argv[2]);

    return 0;
}
