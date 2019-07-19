#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"

#include "../include/lambdas.h"
#include "../include/pjtree.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLatex.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#define FP_TH1_FILL (int (TH1::*)(double))&TH1::Fill
#define FP_TH1_FILLW (int (TH1::*)(double, double))&TH1::Fill
#define FP_TH1_GETBC (double (TH1::*)(int) const)&TH1::GetBinContent
#define FP_TH1_SETBC (void (TH1::*)(int, double))&TH1::SetBinContent
#define FP_TH1_DRAW (void (TH1::*)(char const*))&TH1::Draw

using namespace std::literals::string_literals;
using namespace std::placeholders;

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
        if ((*pjt->jtpt)[j] < jet_pt_min) { continue; }
        if (std::abs((*pjt->jteta)[j]) > jet_eta_abs) { continue; }

        auto jet_eta = (*pjt->jteta)[j];
        auto jet_phi = convert_radian((*pjt->jtphi)[j]);

        auto jet_wta_eta = (*pjt->WTAeta)[j];
        auto jet_wta_phi = convert_radian((*pjt->WTAphi)[j]);

        auto photon_jet_dphi = std::abs(photon_phi - jet_phi);
        auto photon_wta_dphi = std::abs(photon_phi - jet_wta_phi); 

        (*pjet_es_f_dphi)[pthf_x]->Fill(photon_jet_dphi);
        (*pjet_wta_f_dphi)[pthf_x]->Fill(photon_wta_dphi);

        /* require back-to-back jets */
        if (photon_jet_dphi < 0.875_pi) { continue; }

        double jet_pt = (*pjt->jtpt)[j];
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

template <typename... T>
void normalise(std::unique_ptr<history>& norm,
               std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->divide(*norm), 0)... };
}

template <typename... T>
void scale(double factor, std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->scale(factor), 0)... };
}

template <typename... T>
void scale_bin_width(std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        obj->Scale(1., "width"); }), 0)... };
}

template <typename... T>
void normalise_to_unity(std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        obj->Scale(1. / obj->Integral(), "width"); }), 0)... };
}

int populate(char const* config, char const* output) {
    printf("load config options..\n");

    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto mix = conf->get<std::string>("mix");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto frequency = conf->get<int64_t>("frequency");

    auto type = conf->get<std::string>("type");

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

    auto events_to_mix = conf->get<int64_t>("events_to_mix");

    /* convert to integral angle units (cast to double) */
    convert_in_place_pi(rdphi);

    /* default values for config options */
    frequency = frequency ? frequency : 10000;
    events_to_mix = std::max(100L, events_to_mix);

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

        int64_t leading = -1;
        for (int64_t j = 0; j < pjt->nPho; ++j) {
            if ((*pjt->phoEt)[j] < photon_pt_min) { continue; }
            if (std::abs((*pjt->phoEta)[j]) > photon_eta_abs) { continue; }
            if ((*pjt->phoHoverE)[j] > hovere_max) { continue; }

            leading = j;
            break;
        }

        if ((*pjt->phoSigmaIEtaIEta_2012)[leading] > see_max
                || (*pjt->phoSigmaIEtaIEta_2012)[leading] < see_min)
            continue;

        /* require leading photon */
        if (leading < 0) { continue; }

        /* isolation requirement */
        float isolation = (*pjt->pho_ecalClusterIsoR4)[leading]
            + (*pjt->pho_hcalRechitIsoR4)[leading]
            + (*pjt->pho_trackIsoR4PtCut20)[leading];
        if (isolation > iso_max) { continue; }

        double photon_pt = (*pjt->phoEt)[leading];
        auto pt_x = ipt->index_for(photon_pt);

        (*photon_f_pt)(0, FP_TH1_FILL, photon_pt);

        /* set (phi) axis of leading photon */
        auto photon_phi = convert_radian((*pjt->phoPhi)[leading]);

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
    scale(1. / events_to_mix,
        mix_pjet_es_f_dphi, mix_pjet_wta_f_dphi,
        mix_pjet_f_ddr, mix_pjet_f_x);

    /* integrate histograms */
    auto nevt_d_pt = nevt->sum(1);
    auto nevt_d_hf = nevt->sum(0);

    auto pjet_es_f_dphi_d_pt = pjet_es_f_dphi->sum(1);
    auto pjet_wta_f_dphi_d_pt = pjet_wta_f_dphi->sum(1);
    auto pjet_f_x_d_pt = pjet_f_x->sum(1);
    auto pjet_f_x_d_hf = pjet_f_x->sum(0);

    auto pjet_f_ddr_d_pt = pjet_f_ddr->sum(1, 1);
    auto pjet_f_ddr_d_hf = pjet_f_ddr->sum(0, 1);

    /* mixed events */
    auto mix_pjet_es_f_dphi_d_pt = mix_pjet_es_f_dphi->sum(1);
    auto mix_pjet_wta_f_dphi_d_pt = mix_pjet_wta_f_dphi->sum(1);
    auto mix_pjet_f_x_d_pt = mix_pjet_f_x->sum(1);
    auto mix_pjet_f_x_d_hf = mix_pjet_f_x->sum(0);

    auto mix_pjet_f_ddr_d_pt = mix_pjet_f_ddr->sum(1, 1);
    auto mix_pjet_f_ddr_d_hf = mix_pjet_f_ddr->sum(0, 1);

    /* subtract histograms */
    /* *pjet_es_f_dphi -= *mix_pjet_es_f_dphi; */
    /* *pjet_wta_f_dphi -= *mix_pjet_wta_f_dphi; */
    /* *pjet_f_ddr -= *mix_pjet_f_ddr; */
    /* *pjet_f_x -= *mix_pjet_f_x; */

    /* normalise to number of photons (events) */
    normalise(nevt_d_pt,
        pjet_es_f_dphi_d_pt, mix_pjet_es_f_dphi_d_pt,
        pjet_wta_f_dphi_d_pt, mix_pjet_wta_f_dphi_d_pt,
        pjet_f_ddr_d_pt, mix_pjet_f_ddr_d_pt,
        pjet_f_x_d_pt, mix_pjet_f_x_d_pt);

    normalise(nevt_d_hf,
        pjet_f_ddr_d_hf, mix_pjet_f_ddr_d_hf,
        pjet_f_x_d_hf, mix_pjet_f_x_d_hf);

    /* scale by bin width */
    scale_bin_width(
        pjet_f_ddr_d_pt, mix_pjet_f_ddr_d_pt,
        pjet_f_ddr_d_hf, mix_pjet_f_ddr_d_hf,
        pjet_f_x_d_pt, mix_pjet_f_x_d_pt,
        pjet_f_x_d_hf, mix_pjet_f_x_d_hf);

    /* normalise to unity */

    printf("painting..\n");

    auto photon_pt_selection = [&](int64_t index) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{#gamma} < %.0f",
            (*ipt)[index - 1], (*ipt)[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, buffer);
    };

    auto hf_selection = [&](int64_t index) {
        char buffer[128] = { '\0' };
        std::vector<int32_t> bins = { 90, 50, 30, 10, 0 };
        sprintf(buffer, "%i - %i%%", bins[index], bins[index - 1]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, buffer);
    };

    auto hb = new pencil();
    hb->category("axis", "escheme", "wta");
    hb->category("type", "raw", "mix");

    /* hb->alias("raw", "sub."); */

    auto system = "PbPb #sqrt{s_{NN}} = 5.02 TeV"s;

    auto c1 = new paper("dphi_d_pt", hb);
    apply_default_style(c1, system, 0., 0.8);
    c1->accessory(photon_pt_selection);

    for (int64_t i = 0; i < ipt->size() - 1; ++i) {
        c1->add((*pjet_es_f_dphi_d_pt)[i], "escheme", "raw");
        c1->stack((*pjet_wta_f_dphi_d_pt)[i], "wta", "raw");
        c1->stack((*mix_pjet_es_f_dphi_d_pt)[i], "escheme", "mix");
        c1->stack((*mix_pjet_wta_f_dphi_d_pt)[i], "wta", "mix");
    }

    auto c2 = new paper("ddr_d_pt", hb);
    apply_default_style(c2, system, 0., 12.);
    c2->accessory(photon_pt_selection);

    for (int64_t i = 0; i < ipt->size() - 1; ++i) {
        c2->add((*pjet_f_ddr_d_pt)[i], "raw");
        c2->stack((*mix_pjet_f_ddr_d_pt)[i], "mix");
    }

    auto c3 = new paper("ddr_d_hf", hb);
    apply_default_style(c3, system, 0., 12.);
    c3->accessory(hf_selection);

    for (int64_t i = 1; i < ihf->size(); ++i) {
        c3->add((*pjet_f_ddr_d_hf)[i], "raw");
        c3->stack((*mix_pjet_f_ddr_d_hf)[i], "mix");
    }

    auto c4 = new paper("x_d_pt", hb);
    apply_default_style(c4, system, 0., 1.5);
    c4->accessory(photon_pt_selection);

    for (int64_t i = 0; i < ipt->size() - 1; ++i) {
        c4->add((*pjet_f_x_d_pt)[i], "raw");
        c4->stack((*mix_pjet_f_x_d_pt)[i], "mix");
    }

    auto c5 = new paper("x_d_hf", hb);
    apply_default_style(c5, system, 0., 1.5);
    c5->accessory(hf_selection);

    for (int64_t i = 1; i < ihf->size(); ++i) {
        c5->add((*pjet_f_x_d_hf)[i], "raw");
        c5->stack((*mix_pjet_f_x_d_hf)[i], "mix");
    }

    hb->set_binary("type");
    hb->sketch();

    c1->draw("pdf");
    c2->draw("pdf");
    c3->draw("pdf");
    c4->draw("pdf");
    c5->draw("pdf");

    pjet_es_f_dphi->save(type);
    pjet_wta_f_dphi->save(type);
    pjet_f_ddr->save(type);
    pjet_f_x->save(type);

    pjet_es_f_dphi_d_pt->save(type);
    pjet_wta_f_dphi_d_pt->save(type);
    pjet_f_ddr_d_pt->save(type);
    pjet_f_ddr_d_hf->save(type);
    pjet_f_x_d_pt->save(type);
    pjet_f_x_d_hf->save(type);

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
