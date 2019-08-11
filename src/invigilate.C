#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/memory.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <memory>
#include <string>
#include <vector>

using namespace std::literals::string_literals;

static bool in_hem_failure_region(float eta, float phi) {
    return (eta < -1.242 && -1.72 < phi && phi < -0.72);
}

int invigilate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");

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

    auto scale = std::make_unique<memory>("scale"s, "counts", ies, mptetahf);

    /* load input */
    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("pj");
    auto p = new pjtree(true, false, t, { 1, 0, 1, 0, 1, 0 });

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* fill histograms */
    auto nentries = static_cast<int64_t>(t->GetEntries());
    for (int64_t i = 0; i < nentries; ++i) {
        if (i % 100000 == 0) { printf("%li/%li\n", i, nentries); }

        t->GetEntry(i);

        if (p->hiHF <= hf_min) { continue; }

        std::vector<int64_t> exclusion;

        for (int64_t j = 0; j < p->nMC; ++j) {
            auto pid = (*p->mcPID)[j];
            auto mpid = (*p->mcMomPID)[j];
            if (pid != 22 || (std::abs(mpid) > 22 && mpid != -999)) { continue; }

            /* gen isolation requirement */
            float isolation = (*p->mcCalIsoDR04)[j];
            if (isolation > 5.) { continue; }

            exclusion.push_back(j);
        }

        for (int64_t j = 0; j < p->nref; ++j) {
            if ((*p->subid)[j] > 0) { continue; }

            auto gen_pt = (*p->refpt)[j];
            if (gen_pt < pt_min) { continue; }

            auto gen_eta = (*p->refeta)[j];
            if (std::abs(gen_eta) >= eta_max) { continue; }

            bool match = false;
            for (auto const& index : exclusion) {
                auto photon_phi = convert_radian((*p->mcPhi)[index]);
                auto gen_phi = convert_radian((*p->genphi)[j]);

                auto deta = (*p->mcEta)[index] - (*p->geneta)[j];
                auto dphi = revert_radian(photon_phi - gen_phi);
                auto dr2 = deta * deta + dphi * dphi;

                if (dr2 < 0.01) { match = true; break; }
            }

            if (match == true) { continue; }

            if (heavyion && in_hem_failure_region(gen_eta, (*p->refphi)[j]))
                continue;

            auto reco_pt = (*p->jtpt)[j];

            (*scale)[v{gen_pt, gen_eta, p->hiHF}]->Fill(
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
