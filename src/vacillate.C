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

int vacillate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");
    auto eta_max = conf->get<float>("eta_max");

    auto rdrr = conf->get<std::vector<float>>("drr_range");
    auto rdrg = conf->get<std::vector<float>>("drg_range");
    auto rptr = conf->get<std::vector<float>>("ptr_range");
    auto rptg = conf->get<std::vector<float>>("ptg_range");

    auto dhf = conf->get<std::vector<float>>("hf_diff");

    /* exclude most peripheral events */
    auto hf_min = dhf.front();

    /* prepare histograms */
    auto ihf = new interval(dhf);
    auto mhf = new multival(dhf);

    auto mcdr = new multival(rdrr, rdrg);
    auto mcpt = new multival(rptr, rptg);

    auto mr = new multival(rdrr, rptr);
    auto mg = new multival(rdrg, rptg);

    auto fcdr = std::bind(&multival::book<TH2F>, mcdr, _1, _2, _3);
    auto fcpt = std::bind(&multival::book<TH2F>, mcpt, _1, _2, _3);

    auto fc = [&](int64_t, std::string const& name, std::string const& label) {
        return new TH2F(name.data(), (";reco;gen;"s + label).data(),
            mr->size(), 0, mr->size(), mg->size(), 0, mg->size()); };

    auto cdr = new memory<TH2F>("cdr"s, "counts", fcdr, mhf);
    auto cpt = new memory<TH2F>("cpt"s, "counts", fcpt, mhf);

    auto c = new memory<TH2F>("c"s, "counts", fc, mhf);

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

        auto hf_x = ihf->index_for(p->hiHF);

        for (int64_t j = 0; j < p->nref; ++j) {
            if ((*p->subid)[j] > 0) { continue; }

            auto gen_pt = (*p->refpt)[j];
            auto gen_eta = (*p->refeta)[j];
            auto gen_phi = (*p->refphi)[j];

            if (std::abs(gen_eta) >= eta_max) { continue; }

            bool match = false;
            for (auto const& index : exclusion) {
                if (dr2((*p->mcEta)[index], gen_eta,
                        (*p->mcPhi)[index], gen_phi) < 0.01) {
                    match = true; break; }
            }

            if (match == true) { continue; }

            if (heavyion && in_hem_failure_region(gen_eta, gen_phi))
                continue;

            auto reco_pt = (*p->jtpt)[j];
            auto reco_eta = (*p->jteta)[j];
            auto reco_phi = (*p->jtphi)[j];

            auto id = genid[gen_pt];
            auto gdr = std::sqrt(dr2(gen_eta, (*p->WTAgeneta)[id],
                                     gen_phi, (*p->WTAgenphi)[id]));

            if (reco_pt <= rptr.front() || reco_pt >= rptr.back()) {
                (*c)[hf_x]->Fill(-1, mg->index_for(v{gdr, gen_pt}), p->weight);
                continue;
            }

            auto rdr = std::sqrt(dr2(reco_eta, (*p->WTAeta)[j],
                                     reco_phi, (*p->WTAphi)[j]));

            (*cdr)[hf_x]->Fill(rdr, gdr, p->weight);
            (*cpt)[hf_x]->Fill(reco_pt, gen_pt, p->weight);
            (*c)[hf_x]->Fill(mr->index_for(v{rdr, reco_pt}),
                             mg->index_for(v{gdr, gen_pt}),
                             p->weight);
        }
    }

    /* save output */
    in(output, [&]() {
        cdr->save(tag);
        cpt->save(tag);
        c->save(tag);
    });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return vacillate(argv[1], argv[2]);

    return 0;
}