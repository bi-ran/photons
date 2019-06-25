#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLatex.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"

#include "../include/lambdas.h"
#include "../include/pjtree.h"

#define FP_TH1_FILL (int (TH1::*)(double))&TH1::Fill
#define FP_TH1_FILLW (int (TH1::*)(double, double))&TH1::Fill
#define FP_TH1_GETBC (double (TH1::*)(int) const)&TH1::GetBinContent
#define FP_TH1_SETBC (void (TH1::*)(int, double))&TH1::SetBinContent
#define FP_TH1_DRAW (void (TH1::*)(char const*))&TH1::Draw

using namespace std::literals::string_literals;
using namespace std::placeholders;

int64_t leading_track_index(pjtree* pjt, float jet_eta, int64_t jet_phi) {
    int64_t index = -1;
    float max_track_pt = -1.f;

    for (int64_t j = 0; j < pjt->nTrk; ++j) {
        auto track_pt = (*pjt->trkPt)[j];
        if (track_pt < 1.f) { continue; }

        auto track_eta = (*pjt->trkEta)[j];
        if (std::abs(track_eta) > 2.f) { continue; }

        /* track quality selections */
        if (!(*pjt->highPurity)[j]) { continue; }
        if ((*pjt->trkPtError)[j] / (*pjt->trkPt)[j] > 0.1) { continue; }
        if ((*pjt->trkDxy1)[j] / (*pjt->trkDxyError1)[j] > 3.) { continue; }
        if ((*pjt->trkDz1)[j] / (*pjt->trkDzError1)[j] > 3.) { continue; }

        auto track_phi = convert_radian((*pjt->trkPhi)[j]);

        float deta = jet_eta - track_eta;
        float dphi = revert_radian(jet_phi - track_phi);
        float dr2 = deta * deta + dphi * dphi;

        if (track_pt > max_track_pt && dr2 < 0.16) {
            max_track_pt = track_pt;
            index = j;
        }
    }

    return index;
}

void fill_axes(pjtree* pjt, float jet_pt_min, float jet_eta_abs,
               double photon_leading_pt, int64_t photon_phi,
               int64_t photon_pt_x,
               std::shared_ptr<multival>& mdiff,
               std::unique_ptr<history>& nevt,
               std::unique_ptr<history>& pjet_f_x,
               std::unique_ptr<history>& pjet_es_f_dphi,
               std::unique_ptr<history>& pjet_wta_f_dphi,
               std::unique_ptr<history>& pjet_f_ddr) {
    (*nevt)[photon_pt_x]->Fill(1.);

    for (int64_t j = 0; j < pjt->nref; ++j) {
        if ((*pjt->jtpt)[j] < jet_pt_min) { continue; }
        if (std::abs((*pjt->jteta)[j]) > jet_eta_abs) { continue; }

        auto jet_eta = (*pjt->jteta)[j];
        auto jet_phi = convert_radian((*pjt->jtphi)[j]);

        /* find leading track */
        auto track_index = leading_track_index(pjt, jet_eta, jet_phi);

        if (track_index < 0) { return; }

        auto track_eta = (*pjt->trkEta)[track_index];
        auto track_phi = convert_radian((*pjt->trkPhi)[track_index]);

        double jet_pt = (*pjt->jtpt)[j];
        double pjet_x = jet_pt / photon_leading_pt;

        /* calculate index for jet pt, x */
        int64_t index = mdiff->index_for(v{photon_leading_pt, jet_pt, pjet_x});

        auto photon_jet_dphi = std::abs(photon_phi - jet_phi);
        auto photon_track_dphi = std::abs(photon_phi - track_phi); 

        (*pjet_es_f_dphi)[index]->Fill(photon_jet_dphi);
        (*pjet_wta_f_dphi)[index]->Fill(photon_track_dphi);

        /* require back-to-back jets */
        if (photon_jet_dphi < 0.875_pi) { continue; }

        (*pjet_f_x)[photon_pt_x]->Fill(pjet_x);

        double jt_deta = jet_eta - track_eta;
        double jt_dphi = revert_radian(jet_phi - track_phi);
        double jt_dr = jt_deta * jt_deta + jt_dphi * jt_dphi;

        (*pjet_f_ddr)[index]->Fill(std::sqrt(jt_dr));
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

int diffaxis(char const* config, char const* output) {
    printf("load config options..\n");

    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto mix = conf->get<std::string>("mix");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto mc_branches = conf->get<bool>("mc_branches");

    auto freq = conf->get<int64_t>("frequency");

    auto photon_pt_min = conf->get<float>("photon_pt_min");
    auto photon_eta_abs = conf->get<float>("photon_eta_abs");
    auto jet_pt_min = conf->get<float>("jet_pt_min");
    auto jet_eta_abs = conf->get<float>("jet_eta_abs");

    auto rppt = conf->get<std::vector<float>>("ppt_range");
    auto rjpt = conf->get<std::vector<float>>("jpt_range");
    auto rdphi = conf->get<std::vector<float>>("dphi_range");
    auto rx = conf->get<std::vector<float>>("x_range");
    auto rdr = conf->get<std::vector<float>>("dr_range");

    auto dppt = conf->get<std::vector<float>>("ppt_diff");
    auto djpt = conf->get<std::vector<float>>("jpt_diff");
    auto dx = conf->get<std::vector<float>>("x_diff");

    auto events_to_mix = conf->get<int64_t>("events_to_mix");

    /* convert to integral angle units (cast to double) */
    convert_in_place_pi(rdphi);

    /* default values for config options */
    freq = freq ? freq : 10000;
    events_to_mix = std::max(100L, events_to_mix);

    printf("prepare histograms..\n");

    auto incl = std::make_shared<interval>(0.f, 9999.f);
    auto ippt = std::make_shared<interval>(dppt);
    auto ijpt = std::make_shared<interval>(djpt);
    auto ix = std::make_shared<interval>(dx);

    auto mincl = std::make_shared<multival>(*incl);
    auto mppt = std::make_shared<multival>(dppt);
    auto mdiff = std::make_shared<multival>(dppt, djpt, dx);

    auto nevt = std::make_unique<history>("nevt"s, "", incl, mppt);
    auto nmix = std::make_unique<history>("nmix"s, "", incl, mppt);

    auto photon_f_pt = std::make_unique<history>("photon_f_pt"s,
        "dN/dp_{T}^{#gamma}", "p_{T}^{#gamma}", rppt, mincl);

    auto pjet_f_x = std::make_unique<history>("pjet_f_x"s,
        "dN/dx^{#gammaj}", "x^{#gammaj}", rx, mppt);
    auto mix_pjet_f_x = std::make_unique<history>("mix_pjet_f_x"s,
        "dN/dx^{#gammaj}", "x^{#gammaj}", rx, mppt);

    auto pjet_es_f_dphi = std::make_unique<history>("pjet_es_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}", rdphi, mdiff);
    auto pjet_wta_f_dphi = std::make_unique<history>("pjet_wta_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}", rdphi, mdiff);
    auto pjet_f_ddr = std::make_unique<history>("pjet_f_ddr"s,
        "dN/d#Deltar^{jj}", "#Deltar^{jj}", rdr, mdiff);

    auto mix_pjet_es_f_dphi = std::make_unique<history>("mix_pjet_es_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}", rdphi, mdiff);
    auto mix_pjet_wta_f_dphi = std::make_unique<history>("mix_pjet_wta_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}", rdphi, mdiff);
    auto mix_pjet_f_ddr = std::make_unique<history>("mix_pjet_f_ddr",
        "dN/d#Deltar^{jj}", "#Deltar^{jj}", rdr, mdiff);

    printf("iterate..\n");

    TFile* f = new TFile(input.data(), "read");
    TTree* t = (TTree*)f->Get("pj");
    auto pjt = new pjtree(mc_branches, t);

    TFile* fm = new TFile(mix.data(), "read");
    TTree* tm = (TTree*)fm->Get("pj");
    auto pjtm = new pjtree(mc_branches, tm);

    TFile* fout = new TFile(output, "recreate");

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    int64_t nentries = static_cast<int64_t>(t->GetEntries());
    if (max_entries) { nentries = std::min(nentries, max_entries); }
    int64_t mentries = static_cast<int64_t>(tm->GetEntries());
    for (int64_t i = 0, m = 0; i < nentries; ++i) {
        if (i % freq == 0) { printf("entry: %li/%li\n", i, nentries); }

        t->GetEntry(i);

        /* event selections */

        int64_t photon_leading = -1;
        for (int64_t j = 0; j < pjt->nPho; ++j) {
            if ((*pjt->phoEt)[j] < photon_pt_min) { continue; }
            if (std::abs((*pjt->phoEta)[j]) > photon_eta_abs) { continue; }
            if ((*pjt->phoHoverE)[j] > 0.1) { continue; }
            if ((*pjt->phoSigmaIEtaIEta_2012)[j] > 0.01) { continue; }

            /* isolation requirement */
            float isolation = (*pjt->pho_ecalClusterIsoR4)[j]
                + (*pjt->pho_hcalRechitIsoR4)[j]
                + (*pjt->pho_trackIsoR4PtCut20)[j];
            if (isolation > 1.f) { continue; }

            photon_leading = j;
            break;
        }

        /* require leading photon */
        if (photon_leading < 0) { continue; }

        double photon_leading_pt = (*pjt->phoEt)[photon_leading];
        int64_t photon_pt_x = ippt->index_for(photon_leading_pt);

        (*photon_f_pt)(0, FP_TH1_FILL, photon_leading_pt);

        /* set (phi) axis of leading photon */
        auto photon_phi = convert_radian((*pjt->phoPhi)[photon_leading]);

        fill_axes(pjt, jet_pt_min, jet_eta_abs,
                  photon_leading_pt, photon_phi,
                  photon_pt_x, mdiff, nevt,
                  pjet_f_x,
                  pjet_es_f_dphi, pjet_wta_f_dphi,
                  pjet_f_ddr);

        /* mixing events in minimum bias */
        for (int64_t k = 0; k < events_to_mix; ++k) {
            tm->GetEntry(m);

            fill_axes(pjtm, jet_pt_min, jet_eta_abs,
                      photon_leading_pt, photon_phi,
                      photon_pt_x, mdiff, nmix,
                      mix_pjet_f_x,
                      mix_pjet_es_f_dphi, mix_pjet_wta_f_dphi,
                      mix_pjet_f_ddr);

            m = (m + 1) % mentries;
        }
    }

    /* normalise histograms */
    scale(1. / events_to_mix, mix_pjet_es_f_dphi, mix_pjet_wta_f_dphi,
        mix_pjet_f_ddr, mix_pjet_f_x);

    /* integrate histograms */
    auto pjet_es_f_dphi_d_ppt = pjet_es_f_dphi->sum(1, 1);
    auto pjet_wta_f_dphi_d_ppt = pjet_wta_f_dphi->sum(1, 1);
    auto pjet_f_ddr_d_ppt = pjet_f_ddr->sum(1, 1);
    auto pjet_f_ddr_d_x = pjet_f_ddr->sum(0, 0);

    /* mixed events */
    auto mix_pjet_es_f_dphi_d_ppt = mix_pjet_es_f_dphi->sum(1, 1);
    auto mix_pjet_wta_f_dphi_d_ppt = mix_pjet_wta_f_dphi->sum(1, 1);
    auto mix_pjet_f_ddr_d_ppt = mix_pjet_f_ddr->sum(1, 1);
    auto mix_pjet_f_ddr_d_x = mix_pjet_f_ddr->sum(0, 0);

    /* subtract histograms */
    /* *pjet_f_x -= *mix_pjet_f_x; */

    /* *pjet_es_f_dphi_d_ppt -= *mix_pjet_es_f_dphi_d_ppt; */
    /* *pjet_wta_f_dphi_d_ppt -= *mix_pjet_wta_f_dphi_d_ppt; */
    /* *pjet_f_ddr_d_ppt -= *mix_pjet_f_ddr_d_ppt; */

    /* normalise to number of photons (events) */
    normalise(nevt,
        pjet_es_f_dphi_d_ppt, pjet_wta_f_dphi_d_ppt,
        mix_pjet_es_f_dphi_d_ppt, mix_pjet_wta_f_dphi_d_ppt,
        pjet_f_ddr_d_ppt, mix_pjet_f_ddr_d_ppt,
        pjet_f_x, mix_pjet_f_x);

    /* scale by bin width */
    scale_bin_width(pjet_f_ddr_d_ppt, mix_pjet_f_ddr_d_ppt,
        pjet_f_ddr_d_x, mix_pjet_f_ddr_d_x, pjet_f_x, mix_pjet_f_x);

    /* normalise to unity */
    normalise_to_unity(pjet_f_ddr_d_x, mix_pjet_f_ddr_d_x);

    printf("painting..\n");

    auto photon_pt_selection = [&](int64_t index) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{#gamma} < %.0f",
            (*ijpt)[index - 1], (*ijpt)[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, buffer);
    };

    auto x_selection = [&](int64_t index) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.1f < x^{#gammaj} < %.1f",
            (*ix)[index - 1], (*ix)[index]);

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

    auto c1 = new paper("dphi_d_ppt", hb);
    apply_default_style(c1, system, 0., 0.8);
    c1->accessory(photon_pt_selection);

    for (int64_t i = 0; i < ippt->size() - 1; ++i) {
        c1->add((*pjet_es_f_dphi_d_ppt)[i], "escheme", "raw");
        c1->stack((*pjet_wta_f_dphi_d_ppt)[i], "wta", "raw");
        c1->stack((*mix_pjet_es_f_dphi_d_ppt)[i], "escheme", "mix");
        c1->stack((*mix_pjet_wta_f_dphi_d_ppt)[i], "wta", "mix");
    }

    auto c2 = new paper("ddr_d_ppt", hb);
    apply_default_style(c2, system, 0., 12.);
    c2->accessory(photon_pt_selection);

    for (int64_t i = 0; i < ippt->size() - 1; ++i) {
        c2->add((*pjet_f_ddr_d_ppt)[i], "raw");
        c2->stack((*mix_pjet_f_ddr_d_ppt)[i], "mix");
    }

    auto c3 = new paper("x_d_ppt", hb);
    apply_default_style(c3, system, 0., 1.2);
    c3->accessory(photon_pt_selection);

    for (int64_t i = 0; i < ippt->size() - 1; ++i) {
        c3->add((*pjet_f_x)[i], "raw");
        c3->stack((*mix_pjet_f_x)[i], "mix");
    }

    auto c4 = new paper("ddr_d_x", hb);
    apply_default_style(c4, system, 0., 20.);
    c4->accessory(x_selection);

    for (int64_t i = 0; i < ix->size() - 1; ++i) {
        c4->add((*pjet_f_ddr_d_x)[i], "raw");
        c4->stack((*mix_pjet_f_ddr_d_x)[i], "mix");
    }

    hb->set_binary("type");
    hb->sketch();

    c1->draw("pdf");
    c2->draw("pdf");
    c3->draw("pdf");
    c4->draw("pdf");

    /* fout->Write("", TObject::kOverwrite); */

    (void)fout;

    printf("destroying objects..\n");

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return diffaxis(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
