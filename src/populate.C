#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"

#include "../include/pjtree.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define FP_TH1_FILL (int (TH1::*)(double))&TH1::Fill

using namespace std::literals::string_literals;
using namespace std::placeholders;

template <typename... T>
void scale(double factor, std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->scale(factor), 0)... };
}

void fill_axes(pjtree* pjt, float jet_pt_min, float jet_eta_abs,
               double photon_pt, int64_t photon_phi,
               int64_t pt_x, int64_t hf_x, int64_t pthf_x,
               std::shared_ptr<interval>& ix,
               std::unique_ptr<history>& nevt,
               std::unique_ptr<history>& pjet_es_f_dphi,
               std::unique_ptr<history>& pjet_wta_f_dphi,
               std::unique_ptr<history>& pjet_f_x,
               std::unique_ptr<history>& pjet_f_ddr) {
    (*nevt)[pthf_x]->Fill(1.);

    for (int64_t j = 0; j < pjt->nref; ++j) {
        auto jet_pt = (*pjt->jtpt)[j];
        if (jet_pt <= jet_pt_min) { continue; }

        auto jet_eta = (*pjt->jteta)[j];
        if (std::abs(jet_eta) >= jet_eta_abs) { continue; }

        auto jet_phi = convert_radian((*pjt->jtphi)[j]);

        auto jet_wta_eta = (*pjt->WTAeta)[j];
        auto jet_wta_phi = convert_radian((*pjt->WTAphi)[j]);

        auto photon_jet_dphi = std::abs(photon_phi - jet_phi);
        auto photon_wta_dphi = std::abs(photon_phi - jet_wta_phi); 

        (*pjet_es_f_dphi)[pthf_x]->Fill(photon_jet_dphi);
        (*pjet_wta_f_dphi)[pthf_x]->Fill(photon_wta_dphi);

        /* require back-to-back jets */
        if (photon_jet_dphi < 0.875_pi) { continue; }

        double pjet_x = jet_pt / photon_pt;
        (*pjet_f_x)[pthf_x]->Fill(pjet_x);

        /* calculate index */
        auto x_x = ix->index_for(pjet_x);

        double jt_deta = jet_eta - jet_wta_eta;
        double jt_dphi = revert_radian(jet_phi - jet_wta_phi);
        double jt_dr = jt_deta * jt_deta + jt_dphi * jt_dphi;

        (*pjet_f_ddr)[x{pt_x, hf_x, x_x}]->Fill(std::sqrt(jt_dr));
    }
}

static int64_t within_hem_failure_region(pjtree* t, int64_t index) {
    return ((*t->phoSCEta)[index] < -1.3
        && (*t->phoSCPhi)[index] < -0.9
        && (*t->phoSCPhi)[index] > -1.6);
}

static int64_t passes_basic_selections(pjtree* t, int64_t index) {
    return (*t->eleMissHits)[index] <= 1 && (*t->eleIP3D)[index] < 0.03;
}

static int64_t passes_looseid_barrel(pjtree* t, int64_t index, bool pp) {
    if (!passes_basic_selections(t, index)) { return 0; }
    if (!(std::abs((*t->eleSCEta)[index]) < 1.442)) { return 0; }

    if (!pp) {
        if (t->hiBin < 60) {
            return (*t->eleHoverEBc)[index] < 0.1616
                && (*t->eleSigmaIEtaIEta_2012)[index] < 0.0135
                && std::abs((*t->eledEtaSeedAtVtx)[index]) < 0.0038
                && std::abs((*t->eledPhiAtVtx)[index]) < 0.0376
                && std::abs((*t->eleEoverPInv)[index]) < 0.0177;
        }

        return (*t->eleHoverEBc)[index] < 0.1268
            && (*t->eleSigmaIEtaIEta_2012)[index] < 0.0107
            && std::abs((*t->eledEtaSeedAtVtx)[index]) < 0.0035
            && std::abs((*t->eledPhiAtVtx)[index]) < 0.0327
            && std::abs((*t->eleEoverPInv)[index]) < 0.0774;
    } else {
        return (*t->eleHoverE)[index] < 0.02711
            && (*t->eleSigmaIEtaIEta_2012)[index] < 0.01016
            && std::abs((*t->eledEtaSeedAtVtx)[index]) < 0.00316
            && std::abs((*t->eledPhiAtVtx)[index]) < 0.03937
            && std::abs((*t->eleEoverPInv)[index]) < 0.05304;
    }
}

int populate(char const* config, char const* output) {
    printf("load config options..\n");

    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto mix = conf->get<std::string>("mix");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto frequency = conf->get<int64_t>("frequency");
    auto events_to_mix = conf->get<int64_t>("events_to_mix");

    auto type = conf->get<std::string>("type");

    /* options */
    auto pp = conf->get<bool>("pp");
    auto generator_isolation = conf->get<bool>("generator_isolation");
    auto electron_rejection = conf->get<bool>("electron_rejection");

    /* selections */
    auto const photon_pt_min = conf->get<float>("photon_pt_min");
    auto const photon_eta_abs = conf->get<float>("photon_eta_abs");
    auto const jet_pt_min = conf->get<float>("jet_pt_min");
    auto const jet_eta_abs = conf->get<float>("jet_eta_abs");
    auto const hovere_max = conf->get<float>("hovere_max");
    auto const see_min = conf->get<float>("see_min");
    auto const see_max = conf->get<float>("see_max");
    auto const iso_max = conf->get<float>("iso_max");

    auto rppt = conf->get<std::vector<float>>("ppt_range");
    auto rjpt = conf->get<std::vector<float>>("jpt_range");
    auto rdphi = conf->get<std::vector<float>>("dphi_range");
    auto rx = conf->get<std::vector<float>>("x_range");
    auto rdr = conf->get<std::vector<float>>("dr_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dx = conf->get<std::vector<float>>("x_diff");

    /* convert to integral angle units (cast to double) */
    convert_in_place_pi(rdphi);

    /* default values for config options */
    frequency = frequency ? frequency : 10000;

    /* exclude most peripheral events */
    auto hf_min = dhf.front();

    printf("prepare histograms..\n");

    auto incl = std::make_shared<interval>(1, 0.f, 9999.f);
    auto ipt = std::make_shared<interval>(dpt);
    auto ihf = std::make_shared<interval>(dhf);
    auto ix = std::make_shared<interval>(dx);

    auto mincl = std::make_shared<multival>(*incl);
    auto mpthf = std::make_shared<multival>(dpt, dhf);
    auto mpthfx = std::make_shared<multival>(dpt, dhf, dx);

    auto nevt = std::make_unique<history>("nevt"s, "", incl, mpthf);
    auto nmix = std::make_unique<history>("nmix"s, "", incl, mpthf);

    auto photon_f_pt = std::make_unique<history>("photon_f_pt"s,
        "dN/dp_{T}^{#gamma}", "p_{T}^{#gamma}", rppt, mincl);

    auto pjet_es_f_dphi = std::make_unique<history>("pjet_es_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}", rdphi, mpthf);
    auto pjet_wta_f_dphi = std::make_unique<history>("pjet_wta_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}", rdphi, mpthf);
    auto pjet_f_x = std::make_unique<history>("pjet_f_x"s,
        "dN/dx^{#gammaj}", "x^{#gammaj}", rx, mpthf);

    auto mix_pjet_es_f_dphi = std::make_unique<history>("mix_pjet_es_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}", rdphi, mpthf);
    auto mix_pjet_wta_f_dphi = std::make_unique<history>("mix_pjet_wta_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}", rdphi, mpthf);
    auto mix_pjet_f_x = std::make_unique<history>("mix_pjet_f_x"s,
        "dN/dx^{#gammaj}", "x^{#gammaj}", rx, mpthf);

    auto pjet_f_ddr = std::make_unique<history>("pjet_f_ddr"s,
        "dN/d#Deltar^{jj}", "#Deltar^{jj}", rdr, mpthfx);
    auto mix_pjet_f_ddr = std::make_unique<history>("mix_pjet_f_ddr",
        "dN/d#Deltar^{jj}", "#Deltar^{jj}", rdr, mpthfx);

    printf("iterate..\n");

    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("pj");
    auto pjt = new pjtree(false, t);

    TFile* fm = new TFile(mix.data(), "read");
    TTree* tm = (TTree*)fm->Get("pj");
    auto pjtm = new pjtree(false, tm);

    TFile* fout = new TFile(output, "recreate");

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    int64_t nentries = static_cast<int64_t>(t->GetEntries());
    if (max_entries) { nentries = std::min(nentries, max_entries); }
    int64_t mentries = static_cast<int64_t>(tm->GetEntries());
    for (int64_t i = 0, m = 0; i < nentries; ++i) {
        if (i % frequency == 0) { printf("entry: %li/%li\n", i, nentries); }

        t->GetEntry(i);

        if (pjt->hiHF <= hf_min) { continue; }

        int64_t leading = -1;
        for (int64_t j = 0; j < pjt->nPho; ++j) {
            if ((*pjt->phoEt)[j] <= photon_pt_min) { continue; }
            if (std::abs((*pjt->phoSCEta)[j]) >= photon_eta_abs) { continue; }
            if ((*pjt->phoHoverE)[j] > hovere_max) { continue; }

            leading = j;
            break;
        }

        if ((*pjt->phoSigmaIEtaIEta_2012)[leading] > see_max
                || (*pjt->phoSigmaIEtaIEta_2012)[leading] < see_min)
            continue;

        /* require leading photon */
        if (leading < 0) { continue; }

        /* hem failure region exclusion */
        if (!pp && within_hem_failure_region(pjt, leading)) { continue; }

        /* isolation requirement */
        if (generator_isolation) {
            float gen_index = (*pjt->pho_genMatchedIndex)[leading];
            if (gen_index == -1) { continue; }

            float isolation = (*pjt->mcCalIsoDR04)[gen_index];
            if (isolation > iso_max) { continue; }
        } else {
            float isolation = (*pjt->pho_ecalClusterIsoR3)[leading]
                + (*pjt->pho_hcalRechitIsoR3)[leading]
                + (*pjt->pho_trackIsoR3PtCut20)[leading];
            if (isolation > iso_max) { continue; }
        }

        /* leading photon axis */
        auto photon_eta = (*pjt->phoEta)[leading];
        auto photon_phi = convert_radian((*pjt->phoPhi)[leading]);

        /* electron rejection */
        if (electron_rejection) {
            bool electron_match = false;
            for (int64_t j = 0; j < pjt->nEle; ++j) {
                auto deta = photon_eta - (*pjt->eleEta)[j];
                if (deta > 0.1) { continue; }

                auto ele_phi = convert_radian((*pjt->elePhi)[j]);
                auto dphi = revert_radian(photon_phi - ele_phi);
                auto dr2 = deta * deta + dphi * dphi;

                if (dr2 < 0.01 && passes_looseid_barrel(pjt, j, pp)) {
                    electron_match = true;
                    break;
                }
            }

            if (electron_match) { continue; }
        }

        double photon_pt = (*pjt->phoEt)[leading];
        auto pt_x = ipt->index_for(photon_pt);

        (*photon_f_pt)(0, FP_TH1_FILL, photon_pt);

        double hf = pjt->hiHF;
        auto hf_x = ihf->index_for(hf);

        auto pthf_x = mpthf->index_for(x{pt_x, hf_x});

        fill_axes(pjt, jet_pt_min, jet_eta_abs,
                  photon_pt, photon_phi,
                  pt_x, hf_x, pthf_x, ix,
                  nevt, pjet_es_f_dphi, pjet_wta_f_dphi,
                  pjet_f_x, pjet_f_ddr);

        /* mixing events in minimum bias */
        for (int64_t k = 0; k < events_to_mix; m = (m + 1) % mentries) {
            tm->GetEntry(m);

            /* hf within +/- 10% */
            if (std::abs(pjtm->hiHF / pjt->hiHF - 1.) > 0.1) { continue; }

            fill_axes(pjtm, jet_pt_min, jet_eta_abs,
                      photon_pt, photon_phi,
                      pt_x, hf_x, pthf_x, ix,
                      nmix, mix_pjet_es_f_dphi, mix_pjet_wta_f_dphi,
                      mix_pjet_f_x, mix_pjet_f_ddr);

            ++k;
        }
    }

    /* normalise histograms */
    if (events_to_mix > 0)
        scale(1. / events_to_mix,
            mix_pjet_es_f_dphi,
            mix_pjet_wta_f_dphi,
            mix_pjet_f_ddr,
            mix_pjet_f_x);

    /* subtract histograms */
    auto sub_pjet_es_f_dphi = new history(*pjet_es_f_dphi, "sub");
    auto sub_pjet_wta_f_dphi = new history(*pjet_wta_f_dphi, "sub");
    auto sub_pjet_f_x = new history(*pjet_f_x, "sub");
    auto sub_pjet_f_ddr = new history(*pjet_f_ddr, "sub");

    *sub_pjet_es_f_dphi -= *mix_pjet_es_f_dphi;
    *sub_pjet_wta_f_dphi -= *mix_pjet_wta_f_dphi;
    *sub_pjet_f_x -= *mix_pjet_f_x;
    *sub_pjet_f_ddr -= *mix_pjet_f_ddr;

    /* normalise by number of photons (events) */
    sub_pjet_es_f_dphi->divide(*nevt);
    sub_pjet_wta_f_dphi->divide(*nevt);
    sub_pjet_f_x->divide(*nevt);
    sub_pjet_f_ddr->divide(*nevt);

    /* save histograms */
    nevt->save(type);
    nmix->save(type);

    photon_f_pt->save(type);

    pjet_es_f_dphi->save(type);
    pjet_wta_f_dphi->save(type);
    pjet_f_x->save(type);
    pjet_f_ddr->save(type);

    mix_pjet_es_f_dphi->save(type);
    mix_pjet_wta_f_dphi->save(type);
    mix_pjet_f_x->save(type);
    mix_pjet_f_ddr->save(type);

    sub_pjet_es_f_dphi->save(type);
    sub_pjet_wta_f_dphi->save(type);
    sub_pjet_f_x->save(type);
    sub_pjet_f_ddr->save(type);

    fout->Close();

    printf("destroying objects..\n");

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return populate(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
