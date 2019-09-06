#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/memory.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"
#include "../git/tricks-and-treats/include/trunk.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

static bool in_hem_failure_region(float eta, float phi) {
    return (eta < -1.242 && -1.72 < phi && phi < -0.72);
}

static float dr2(float eta1, float eta2, float phi1, float phi2) {
    auto deta = eta1 - eta2;
    auto dphi = revert_radian(convert_radian(phi1) - convert_radian(phi2));

    return deta * deta + dphi * dphi;
}

template <typename T>
static int sgn(T val) { return (T(0) < val) - (val < T(0)); }

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

    auto rar = conf->get<std::vector<float>>("ar_range");
    auto rag = conf->get<std::vector<float>>("ag_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto deta = conf->get<std::vector<float>>("eta_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");

    /* exclude most peripheral events */
    auto hf_min = dhf.front();

    /* prepare histograms */
    auto mptetahf = new multival(dpt, deta, dhf);

    auto mes = new multival("energy scale"s, res[0], res[1], res[2]);
    auto mdr = new multival("#deltar^{2}"s, rdr[0], rdr[1], rdr[2]);
    auto mde = new multival("#delta#eta"s, rde[0], rde[1], rde[2]);
    auto mdp = new multival("#delta#phi"s, rdp[0], rdp[1], rdp[2]);

    auto iar = new interval("#deltaj"s, rar[0], rar[1], rar[2]);
    auto iag = new interval("#deltaj"s, rag[0], rag[1], rag[2]);

    auto marag = new multival(*iar, *iag);

    auto fes = std::bind(&multival::book<TH1F>, mes, _1, _2);
    auto fdr = std::bind(&multival::book<TH1F>, mdr, _1, _2);
    auto fde = std::bind(&multival::book<TH1F>, mde, _1, _2);
    auto fdp = std::bind(&multival::book<TH1F>, mdp, _1, _2);

    auto farag = std::bind(&multival::book<TH2F>, marag, _1, _2);

    auto scale = new memory<TH1F>("scale"s, "counts", fes, mptetahf);
    auto angle = new memory<TH1F>("angle"s, "counts", fdr, mptetahf);
    auto eta = new memory<TH1F>("eta"s, "counts", fde, mptetahf);
    auto phi = new memory<TH1F>("phi"s, "counts", fdp, mptetahf);

    auto axis = new memory<TH2F>("axis"s, "counts", farag, mptetahf);

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

        std::unordered_map<float, int64_t> genid;
        for (int64_t j = 0; j < p->ngen; ++j)
            genid[(*p->genpt)[j]] = j;

        for (int64_t j = 0; j < p->nref; ++j) {
            if ((*p->subid)[j] > 0) { continue; }

            auto gen_pt = (*p->refpt)[j];
            if (gen_pt < pt_min) { continue; }

            auto gen_eta = (*p->refeta)[j];
            if (std::abs(gen_eta) >= eta_max) { continue; }

            bool match = false;
            for (auto const& index : exclusion) {
                if (dr2((*p->mcEta)[index], (*p->refeta)[j],
                        (*p->mcPhi)[index], (*p->refphi)[j]) < 0.01) {
                    match = true; break; }
            }

            if (match == true) { continue; }

            if (heavyion && in_hem_failure_region(gen_eta, (*p->refphi)[j]))
                continue;

            auto index = mptetahf->index_for(v{gen_pt, gen_eta, p->hiHF});

            (*scale)[index]->Fill((*p->jtpt)[j] / gen_pt, p->weight);

            auto deta = (*p->jteta)[j] - (*p->refeta)[j];
            auto dphi = revert_radian(convert_radian((*p->jtphi)[j])
                - convert_radian((*p->refphi)[j]));

            (*eta)[index]->Fill(deta, p->weight);
            (*phi)[index]->Fill(dphi, p->weight);

            (*angle)[index]->Fill(sgn((*p->refphi)[j])
                * (deta * deta + dphi * dphi), p->weight);

            auto id = genid[gen_pt];
            auto rdr2 = dr2((*p->jteta)[j], (*p->WTAeta)[j],
                            (*p->jtphi)[j], (*p->WTAphi)[j]);
            auto gdr2 = dr2((*p->refeta)[j], (*p->WTAgeneta)[id],
                            (*p->refphi)[j], (*p->WTAgenphi)[id]);

            (*axis)[index]->Fill(std::sqrt(rdr2), std::sqrt(gdr2), p->weight);
        }
    }

    /* save output */
    in(output, [&]() {
        scale->save(tag);
        angle->save(tag);
        eta->save(tag);
        phi->save(tag);

        axis->save(tag);
    });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return fabulate(argv[1], argv[2]);

    return 0;
}
