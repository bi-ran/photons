#include "../include/pjtree.h"
#include "../include/specifics.h"

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

    auto start = conf->get<int64_t>("start");
    auto end = conf->get<int64_t>("end");

    auto heavyion = conf->get<bool>("heavyion");
    auto jet_eta_max = conf->get<float>("jet_eta_max");
    auto photon_pt_min = conf->get<float>("photon_pt_min");
    auto photon_eta_max = conf->get<float>("photon_eta_max");
    auto hovere_max = conf->get<float>("hovere_max");
    auto see_min = conf->get<float>("see_min");
    auto see_max = conf->get<float>("see_max");
    auto iso_max = conf->get<float>("iso_max");

    auto rdrr = conf->get<std::vector<float>>("drr_range");
    auto rdrg = conf->get<std::vector<float>>("drg_range");
    auto rptr = conf->get<std::vector<float>>("ptr_range");
    auto rptg = conf->get<std::vector<float>>("ptg_range");

    auto dhf = conf->get<std::vector<float>>("hf_diff");

    /* exclude most peripheral events */
    auto hf_min = dhf.front();

    /* prepare histograms */
    auto incl = new interval(""s, 1, 0.f, 9999.f);
    auto ihf = new interval(dhf);

    auto mcdr = new multival(rdrr, rdrg);
    auto mcpt = new multival(rptr, rptg);

    auto mr = new multival(rdrr, rptr);
    auto mg = new multival(rdrg, rptg);

    auto fn = std::bind(&interval::book<TH1F>, incl, _1, _2, _3);
    auto fcdr = std::bind(&multival::book<TH2F>, mcdr, _1, _2, _3);
    auto fcpt = std::bind(&multival::book<TH2F>, mcpt, _1, _2, _3);

    auto fr = [&](int64_t, std::string const& name, std::string const& label) {
        return new TH1F(name.data(), (";reco;"s + label).data(),
            mr->size(), 0, mr->size()); };

    auto fg = [&](int64_t, std::string const& name, std::string const& label) {
        return new TH1F(name.data(), (";gen;"s + label).data(),
            mg->size(), 0, mg->size()); };

    auto fc = [&](int64_t, std::string const& name, std::string const& label) {
        return new TH2F(name.data(), (";reco;gen;"s + label).data(),
            mr->size(), 0, mr->size(), mg->size(), 0, mg->size()); };

    auto n = new history<TH1F>("n"s, "events", fn, ihf->size());
    auto r = new history<TH1F>("r"s, "counts", fr, ihf->size());
    auto g = new history<TH1F>("g"s, "counts", fg, ihf->size());
    auto cdr = new history<TH2F>("cdr"s, "counts", fcdr, ihf->size());
    auto cpt = new history<TH2F>("cpt"s, "counts", fcpt, ihf->size());
    auto c = new history<TH2F>("c"s, "counts", fc, ihf->size());

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load input */
    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("pj");
    auto p = new pjtree(true, false, t, { 1, 1, 1, 0, 1, 0 });

    /* fill histograms */
    if (!end) { end = t->GetEntries(); }
    for (int64_t i = start; i < end; ++i) {
        if (i % 100000 == 0) { printf("%li/%li\n", i, end); }

        t->GetEntry(i);

        if (p->hiHF <= hf_min) { continue; }

        int64_t leading = -1;
        for (int64_t j = 0; j < p->nPho; ++j) {
            if ((*p->phoEt)[j] <= photon_pt_min) { continue; }
            if (std::abs((*p->phoSCEta)[j]) >= photon_eta_max) { continue; }
            if ((*p->phoHoverE)[j] > hovere_max) { continue; }

            leading = j;
            break;
        }

        /* require leading photon */
        if (leading < 0) { continue; }

        if ((*p->phoSigmaIEtaIEta_2012)[leading] > see_max
                || (*p->phoSigmaIEtaIEta_2012)[leading] < see_min)
            continue;

        /* hem failure region exclusion */
        if (heavyion && within_hem_failure_region(p, leading)) { continue; }

        /* require match to gen */
        auto gen_index = (*p->pho_genMatchedIndex)[leading];
        if (gen_index == -1) { continue; }

        /* isolation requirement */
        float isolation = (*p->pho_ecalClusterIsoR3)[leading]
            + (*p->pho_hcalRechitIsoR3)[leading]
            + (*p->pho_trackIsoR3PtCut20)[leading];
        if (isolation > iso_max) { continue; }

        /* photon axis */
        auto photon_eta = (*p->phoEta)[leading];
        auto photon_phi = convert_radian((*p->phoPhi)[leading]);

        /* electron rejection */
        bool electron = false;
        for (int64_t j = 0; j < p->nEle; ++j) {
            if (std::abs((*p->eleSCEta)[j]) > 1.4442) { continue; }

            auto deta = photon_eta - (*p->eleEta)[j];
            if (deta > 0.1) { continue; }

            auto ele_phi = convert_radian((*p->elePhi)[j]);
            auto dphi = revert_radian(photon_phi - ele_phi);
            auto dr2 = deta * deta + dphi * dphi;

            if (dr2 < 0.01 && passes_electron_id<
                        det::barrel, wp::loose, pjtree
                    >(p, j, heavyion)) {
                electron = true; break; }
        }

        if (electron) { continue; }

        /* fill event weight */
        auto hf_x = ihf->index_for(p->hiHF);
        (*n)[hf_x]->Fill(1., p->weight);

        /* map reco jet to gen jet */
        std::unordered_map<float, int64_t> genid;
        for (int64_t j = 0; j < p->ngen; ++j)
            genid[(*p->genpt)[j]] = j;

        for (int64_t j = 0; j < p->nref; ++j) {
            if ((*p->subid)[j] > 0) { continue; }

            auto gen_pt = (*p->refpt)[j];
            auto gen_eta = (*p->refeta)[j];
            auto gen_phi = (*p->refphi)[j];

            if (gen_pt < rptg.front()) { continue; }
            if (std::abs(gen_eta) >= jet_eta_max) { continue; }

            auto reco_pt = (*p->jtpt)[j];
            auto reco_eta = (*p->jteta)[j];
            auto reco_phi = (*p->jtphi)[j];

            if (heavyion && in_hem_failure_region(reco_eta, reco_phi))
                continue;

            /* require back-to-back jets */
            if (std::abs(photon_phi - convert_radian(reco_phi)) < 0.875_pi)
                continue;

            auto id = genid[gen_pt];
            auto gdr = std::sqrt(dr2(gen_eta, (*p->WTAgeneta)[id],
                                     gen_phi, (*p->WTAgenphi)[id]));
            auto g_x = mg->index_for(v{gdr, gen_pt});

            (*g)[hf_x]->Fill(g_x, p->weight);

            if (reco_pt > rptr.front() && reco_pt < rptr.back()) {
                auto rdr = std::sqrt(dr2(reco_eta, (*p->WTAeta)[j],
                                         reco_phi, (*p->WTAphi)[j]));
                auto r_x = mr->index_for(v{rdr, reco_pt});

                (*r)[hf_x]->Fill(r_x, p->weight);

                (*cdr)[hf_x]->Fill(rdr, gdr, p->weight);
                (*cpt)[hf_x]->Fill(reco_pt, gen_pt, p->weight);
                (*c)[hf_x]->Fill(r_x, g_x, p->weight);
            } else {
                /* missed */
                (*cpt)[hf_x]->Fill(-1, gen_pt, p->weight);
                (*c)[hf_x]->Fill(-1, g_x, p->weight);
            }
        }
    }

    r->divide(*n);
    g->divide(*n);

    /* save output */
    in(output, [&]() {
        n->save(tag);
        r->save(tag);
        g->save(tag);

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
