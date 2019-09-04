#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/memory.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"
#include "../git/tricks-and-treats/include/trunk.h"

#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;

static bool in_hem_failure_region(float eta, float phi) {
    return (eta < -1.242 && -1.72 < phi && phi < -0.72);
}

static float dr2(float eta1, float eta2, float phi1, float phi2) {
    auto deta = eta1 - eta2;
    auto dphi = revert_radian(convert_radian(phi1) - convert_radian(phi2));

    return deta * deta + dphi * dphi;
}

template <typename T>
static int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int fabulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");

    auto pt_min = conf->get<float>("pt_min");
    auto eta_max = conf->get<float>("eta_max");

    auto res = conf->get<std::vector<float>>("es_range");
    auto rdr = conf->get<std::vector<float>>("dr_range");
    auto rde = conf->get<std::vector<float>>("de_range");
    auto rdp = conf->get<std::vector<float>>("dp_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto deta = conf->get<std::vector<float>>("eta_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");

    /* exclude most peripheral events */
    auto hf_min = dhf.front();

    /* prepare histograms */
    auto ies = new interval("energy scale"s, (int64_t)res[0], res[1], res[2]);
    auto idr = new interval("#deltar^{2}"s, (int64_t)rdr[0], rdr[1], rdr[2]);
    auto ide = new interval("#delta#eta"s, (int64_t)rde[0], rde[1], rde[2]);
    auto idp = new interval("#delta#phi"s, (int64_t)rdp[0], rdp[1], rdp[2]);

    auto mptetahf = new multival(dpt, deta, dhf);

    auto scale = new memory<TH1F>("scale"s, "counts", ies, mptetahf);
    auto angle = new memory<TH1F>("angle"s, "counts", idr, mptetahf);
    auto eta = new memory<TH1F>("eta"s, "counts", ide, mptetahf);
    auto phi = new memory<TH1F>("phi"s, "counts", idp, mptetahf);

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load input */
    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("pj");
    auto p = new pjtree(true, false, t, { 1, 1, 1, 0, 1, 0 });

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
                if (dr2((*p->mcEta)[index], (*p->geneta)[j],
                        (*p->mcPhi)[index], (*p->genphi)[j]) < 0.01) {
                    match = true; break; }
            }

            if (match == true) { continue; }

            if (heavyion && in_hem_failure_region(gen_eta, (*p->refphi)[j]))
                continue;

            auto index = mptetahf->index_for(v{gen_pt, gen_eta, p->hiHF});

            auto deta = (*p->jteta)[j] - (*p->refeta)[j];
            auto dphi = revert_radian(convert_radian((*p->jtphi)[j])
                - convert_radian((*p->refphi)[j]));

            (*scale)[index]->Fill((*p->jtpt)[j] / gen_pt, p->weight);
            (*angle)[index]->Fill(sgn((*p->refphi)[j])
                * (deta * deta + dphi * dphi), p->weight);

            (*eta)[index]->Fill(deta, p->weight);
            (*phi)[index]->Fill(dphi, p->weight);
        }
    }

    /* save output */
    in(output, [&]() {
        scale->save(tag);
        angle->save(tag);

        eta->save(tag);
        phi->save(tag);
    });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return fabulate(argv[1], argv[2]);

    return 0;
}
