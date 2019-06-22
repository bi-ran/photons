#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLatex.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../include/pjtree.h"

#include "../include/integral_angles.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../include/lambdas.h"

#define FP_TH1_FILL (int (TH1::*)(double))&TH1::Fill
#define FP_TH1_FILLW (int (TH1::*)(double, double))&TH1::Fill
#define FP_TH1_GETBC (double (TH1::*)(int) const)&TH1::GetBinContent
#define FP_TH1_SETBC (void (TH1::*)(int, double))&TH1::SetBinContent
#define FP_TH1_DRAW (void (TH1::*)(char const*))&TH1::Draw

using namespace std::literals::string_literals;
using namespace std::placeholders;

using x = std::initializer_list<int64_t>;
using v = std::initializer_list<double>;

void fill_tracks(pjtree* pjt, float trk_pt_min, float trk_eta_abs,
                 double photon_leading_pt, int64_t photon_phi, int64_t photon_pt_x,
                 std::shared_ptr<interval>& idphi,
                 std::shared_ptr<multival>& mdphi,
                 std::unique_ptr<history>& ntrk,
                 std::unique_ptr<history>& sumpt,
                 std::unique_ptr<history>& ntrk_f_pt,
                 std::unique_ptr<history>& sumpt_f_pt,
                 std::unique_ptr<history>& trk_f_dphi,
                 std::unique_ptr<history>& trk_f_pt,
                 std::unique_ptr<history>& evt_f_ntrk,
                 std::unique_ptr<history>& evt_f_sumpt) {
    for (int64_t j = 0; j < pjt->nTrk; ++j) {
        if ((*pjt->trkPt)[j] < trk_pt_min) { continue; }
        if (std::abs((*pjt->trkEta)[j]) > trk_eta_abs) { continue; }

        /* track quality selections */
        if (!(*pjt->highPurity)[j]) { continue; }
        if ((*pjt->trkPtError)[j] / (*pjt->trkPt)[j] > 0.1) { continue; }
        if ((*pjt->trkDxy1)[j] / (*pjt->trkDxyError1)[j] > 3.) { continue; }
        if ((*pjt->trkDz1)[j] / (*pjt->trkDzError1)[j] > 3.) { continue; }

        double trk_pt = (*pjt->trkPt)[j];
        auto trk_phi = convert_radians((*pjt->trkPhi)[j]);
        auto photon_trk_dphi = std::abs(photon_phi - trk_phi);
        int64_t dphi_x = idphi->index_for(photon_trk_dphi);

        (*ntrk)(dphi_x, FP_TH1_FILLW, 1., 1.);
        (*sumpt)(dphi_x, FP_TH1_FILLW, 1., trk_pt);

        (*trk_f_dphi)(photon_pt_x, FP_TH1_FILL,
                      static_cast<double>(photon_trk_dphi));
        (*trk_f_pt)(x{photon_pt_x, dphi_x}, FP_TH1_FILL, trk_pt);
    }

    double norm = 2. * trk_eta_abs * 2. * M_PI / 3.;

    ntrk->scale(1. / norm);
    sumpt->scale(1. / norm);

    for (int64_t j = 0; j < mdphi->size(); ++j) {
        double evt_ntrk = (*ntrk)(j, FP_TH1_GETBC, 1);
        double evt_sumpt = (*sumpt)(j, FP_TH1_GETBC, 1);

        (*evt_f_ntrk)(x{photon_pt_x, j}, FP_TH1_FILL, evt_ntrk);
        (*evt_f_sumpt)(x{photon_pt_x, j}, FP_TH1_FILL, evt_sumpt);

        (*ntrk_f_pt)(j, FP_TH1_FILLW, photon_leading_pt, evt_ntrk);
        (*sumpt_f_pt)(j, FP_TH1_FILLW, photon_leading_pt, evt_sumpt);
    }
}

void fill_jets(pjtree* pjt, float jet_pt_min, float jet_eta_abs,
               double photon_leading_pt, int64_t photon_phi, int64_t photon_pt_x,
               std::unique_ptr<history>& ntrk,
               std::unique_ptr<history>& sumpt,
               std::shared_ptr<interval>& intrk,
               std::shared_ptr<interval>& isumpt,
               std::unique_ptr<history>& nevt,
               std::unique_ptr<history>& photon_f_pt,
               std::unique_ptr<history>& pjet_f_dphi,
               std::unique_ptr<history>& pjet_f_jetpt,
               std::unique_ptr<history>& pjet_f_x) {
    int64_t near_ntrk_x = intrk->index_for((*ntrk)(0, FP_TH1_GETBC, 1));
    int64_t perp_ntrk_x = intrk->index_for((*ntrk)(1, FP_TH1_GETBC, 1));

    int64_t near_sumpt_x = isumpt->index_for((*sumpt)(0, FP_TH1_GETBC, 1));
    int64_t perp_sumpt_x = isumpt->index_for((*sumpt)(1, FP_TH1_GETBC, 1));

    auto index = nevt->index_for(x{photon_pt_x, near_ntrk_x, perp_ntrk_x,
                                   near_sumpt_x, perp_sumpt_x});

    (*nevt)(index, FP_TH1_FILL, 1.);
    (*photon_f_pt)[x{near_ntrk_x, perp_ntrk_x,
                     near_sumpt_x, perp_sumpt_x}]->Fill(photon_leading_pt);

    for (int64_t j = 0; j < pjt->nref; ++j) {
        if ((*pjt->jtpt)[j] < jet_pt_min) { continue; }
        if (std::abs((*pjt->jteta)[j]) > jet_eta_abs) { continue; }

        auto jet_phi = convert_radians((*pjt->jtphi)[j]);
        auto photon_jet_dphi = std::abs(photon_phi - jet_phi);

        (*pjet_f_dphi)(index, FP_TH1_FILL, static_cast<double>(photon_jet_dphi));

        /* require back-to-back jets */
        if (photon_jet_dphi < 0.875_pi) { continue; }

        double jet_pt = (*pjt->jtpt)[j];

        (*pjet_f_jetpt)(index, FP_TH1_FILL, jet_pt);
        (*pjet_f_x)(index, FP_TH1_FILL, jet_pt / photon_leading_pt);
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

int flatten(char const* config, char const* output) {
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
    auto trk_pt_min = conf->get<float>("trk_pt_min");
    auto trk_eta_abs = conf->get<float>("trk_eta_abs");

    auto rpt = conf->get<std::vector<float>>("pt_range");
    auto rtrkpt = conf->get<std::vector<float>>("trkpt_range");
    auto rdphi = conf->get<std::vector<float>>("dphi_range");
    auto rx = conf->get<std::vector<float>>("x_range");
    auto rntrk = conf->get<std::vector<float>>("ntrk_range");
    auto rsumpt = conf->get<std::vector<float>>("sumpt_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto ddphi = conf->get<std::vector<float>>("dphi_diff");
    auto dntrk = conf->get<std::vector<float>>("ntrk_diff");
    auto dsumpt = conf->get<std::vector<float>>("sumpt_diff");

    auto events_to_mix = conf->get<int64_t>("events_to_mix");

    /* convert to integral angle units (cast to double) */
    for (auto& i : rdphi) { i = convert_pis(i); }
    for (auto& i : ddphi) { i = convert_degrees(i); }

    /* default values for config options */
    freq = freq ? freq : 10000;
    events_to_mix = std::max(100L, events_to_mix);

    printf("prepare histograms..\n");

    auto incl = std::make_shared<interval>(0.f, 9999.f);
    auto ipt = std::make_shared<interval>(dpt);
    auto idphi = std::make_shared<interval>(ddphi);
    auto intrk = std::make_shared<interval>(dntrk);
    auto isumpt = std::make_shared<interval>(dsumpt);

    auto mincl = std::make_shared<multival>(*incl);
    auto mpt = std::make_shared<multival>(dpt);
    auto mdphi = std::make_shared<multival>(ddphi);
    auto mptdphi = std::make_shared<multival>(dpt, ddphi);
    auto mntrk2sumpt2 = std::make_shared<multival>(
        dntrk, dntrk, dsumpt, dsumpt);
    auto mptntrk2sumpt2 = std::make_shared<multival>(
        dpt, dntrk, dntrk, dsumpt, dsumpt);

    auto photon_f_pt = std::make_unique<history>("photon_f_pt"s,
        "dN/dp_{T}^{#gamma}", "p_{T}^{#gamma}", rpt, mntrk2sumpt2);
    auto mix_photon_f_pt = std::make_unique<history>("mix_photon_f_pt"s,
        "dN/dp_{T}^{#gamma}", "p_{T}^{#gamma}", rpt, mntrk2sumpt2);

    auto ntrk = std::make_unique<history>("ntrk"s,
        "N^{h^{#pm}}", incl, mdphi);
    auto sumpt = std::make_unique<history>("sumpt"s,
        "#sump_{T}^{h^{#pm}}", incl, mdphi);
    auto mix_ntrk = std::make_unique<history>("mix_ntrk"s,
        "N^{h^{#pm}}", incl, mdphi);
    auto mix_sumpt = std::make_unique<history>("mix_sumpt"s,
        "#sump_{T}^{h^{#pm}}", incl, mdphi);

    auto ntrk_f_pt = std::make_unique<history>("ntrk_f_pt"s,
        "dN^{h^{#pm}}/d#etad#phi", "p_{T}^{#gamma}", rpt, mdphi);
    auto sumpt_f_pt = std::make_unique<history>("sumpt_f_pt"s,
        "d#Sigmap_{T}^{h^{#pm}}/d#etad#phi", "p_{T}^{#gamma}", rpt, mdphi);
    auto mix_ntrk_f_pt = std::make_unique<history>("mix_ntrk_f_pt"s,
        "dN^{h^{#pm}}/d#etad#phi", "p_{T}^{#gamma}", rpt, mdphi);
    auto mix_sumpt_f_pt = std::make_unique<history>("mix_sumpt_f_pt"s,
        "d#Sigmap_{T}^{h^{#pm}}/d#etad#phi", "p_{T}^{#gamma}", rpt, mdphi);

    auto evt_f_ntrk = std::make_unique<history>("evt_f_ntrk"s,
        "dN/dN_{h^{#pm}}", "N^{h^{#pm}}", rntrk, mptdphi);
    auto evt_f_sumpt = std::make_unique<history>("evt_f_sumpt"s,
        "dN/d#Sigmap_{T}^{h^{#pm}}", "#Sigmap_{T}^{h^{#pm}}", rsumpt, mptdphi);
    auto mix_evt_f_ntrk = std::make_unique<history>("mix_evt_f_ntrk"s,
        "dN/dN_{h^{#pm}}", "N^{h^{#pm}}", rntrk, mptdphi);
    auto mix_evt_f_sumpt = std::make_unique<history>("mix_evt_f_sumpt"s,
        "dN/d#Sigmap_{T}^{h^{#pm}}", "#Sigmap_{T}^{h^{#pm}}", rsumpt, mptdphi);

    auto trk_f_dphi = std::make_unique<history>("trk_f_dphi"s,
        "dN/d#Delta#phi^{#gammah}", "#Delta#phi^{#gammah}", rdphi, mpt);
    auto trk_f_pt = std::make_unique<history>("trk_f_pt"s,
        "dN/dp_{T}^{h^{#pm}}", "p_{T}^{h^{#pm}}", rtrkpt, mptdphi);
    auto mix_trk_f_dphi = std::make_unique<history>("mix_trk_f_dphi"s,
        "dN/d#Delta#phi^{#gammah}", "#Delta#phi^{#gammah}", rdphi, mpt);
    auto mix_trk_f_pt = std::make_unique<history>("mix_trk_f_pt"s,
        "dN/dp_{T}^{h^{#pm}}", "p_{T}^{h^{#pm}}", rtrkpt, mptdphi);

    auto nevt = std::make_unique<history>("nevt"s, "", incl, mptntrk2sumpt2);
    auto nmix = std::make_unique<history>("nmix"s, "", incl, mptntrk2sumpt2);

    auto pjet_f_dphi = std::make_unique<history>("pjet_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}",
        rdphi, mptntrk2sumpt2);
    auto pjet_f_jetpt = std::make_unique<history>("pjet_f_jetpt"s,
        "dN/dp_{T}^{j}", "p_{T}^{j}", rpt, mptntrk2sumpt2);
    auto pjet_f_x = std::make_unique<history>("pjet_f_x"s,
        "dN/dx^{#gammaj}", "x^{#gammaj}", rx, mptntrk2sumpt2);

    auto mix_pjet_f_dphi = std::make_unique<history>("mix_pjet_f_dphi"s,
        "dN/d#Delta#phi^{#gammaj}", "#Delta#phi^{#gammaj}",
        rdphi, mptntrk2sumpt2);
    auto mix_pjet_f_jetpt = std::make_unique<history>("mix_pjet_f_jetpt"s,
        "dN/dp_{T}^{j}", "p_{T}^{j}", rpt, mptntrk2sumpt2);
    auto mix_pjet_f_x = std::make_unique<history>("mix_pjet_f_x"s,
        "dN/dx^{#gammaj}", "x^{#gammaj}", rx, mptntrk2sumpt2);

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
        int64_t photon_pt_x = ipt->index_for(photon_leading_pt);

        /* set (phi) axis of leading photon */
        auto photon_phi = convert_radians((*pjt->phoPhi)[photon_leading]);

        /* clear counts (sums) */
        for (int64_t j = 0; j < mdphi->size(); ++j) {
            (*ntrk)(j, FP_TH1_SETBC, 1, 0.);
            (*sumpt)(j, FP_TH1_SETBC, 1, 0.);
            (*mix_ntrk)(j, FP_TH1_SETBC, 1, 0.);
            (*mix_sumpt)(j, FP_TH1_SETBC, 1, 0.);
        }

        fill_tracks(pjt, trk_pt_min, trk_eta_abs,
                    photon_leading_pt, photon_phi, photon_pt_x,
                    idphi, mdphi,
                    ntrk, sumpt,
                    ntrk_f_pt, sumpt_f_pt,
                    trk_f_dphi, trk_f_pt,
                    evt_f_ntrk, evt_f_sumpt);

        fill_jets(pjt, jet_pt_min, jet_eta_abs,
                  photon_leading_pt, photon_phi, photon_pt_x,
                  ntrk, sumpt, intrk, isumpt,
                  nevt, photon_f_pt,
                  pjet_f_dphi, pjet_f_jetpt, pjet_f_x);

        /* mixing events in minimum bias */
        for (int64_t k = 0; k < events_to_mix; ++k) {
            tm->GetEntry(m);

            fill_tracks(pjtm, trk_pt_min, trk_eta_abs,
                        photon_leading_pt, photon_phi, photon_pt_x,
                        idphi, mdphi,
                        mix_ntrk, mix_sumpt,
                        mix_ntrk_f_pt, mix_sumpt_f_pt,
                        mix_trk_f_dphi, mix_trk_f_pt,
                        mix_evt_f_ntrk, mix_evt_f_sumpt);

            fill_jets(pjtm, jet_pt_min, jet_eta_abs,
                      photon_leading_pt, photon_phi, photon_pt_x,
                      mix_ntrk, mix_sumpt, intrk, isumpt,
                      nmix, mix_photon_f_pt,
                      mix_pjet_f_dphi, mix_pjet_f_jetpt, mix_pjet_f_x);

            m = (m + 1) % mentries;
        }
    }

    /* normalise histograms */
    scale(1. / events_to_mix, mix_pjet_f_dphi, mix_pjet_f_jetpt, mix_pjet_f_x,
        mix_ntrk_f_pt, mix_sumpt_f_pt, mix_evt_f_ntrk, mix_evt_f_sumpt,
        mix_photon_f_pt);

    /* integrate histograms */
    /* photon (event) count */
    auto nevt_d_photon_pt = nevt->sum(1, 1, 1, 1);
    auto nevt_d_np_ntrk_np_sumpt = nevt->sum(0);

    auto nevt_d_np_ntrk = nevt_d_np_ntrk_np_sumpt->sum(2, 2);
    auto nevt_d_perp_ntrk = nevt_d_np_ntrk->sum(0);
    auto nevt_d_near_ntrk = nevt_d_np_ntrk->sum(1);

    auto nevt_d_np_sumpt = nevt_d_np_ntrk_np_sumpt->sum(0, 0);
    auto nevt_d_perp_sumpt = nevt_d_np_sumpt->sum(0);
    auto nevt_d_near_sumpt = nevt_d_np_sumpt->sum(1);

    /* photon pt spectra */
    auto photon_f_pt_incl = photon_f_pt->sum(0, 0, 0, 0);

    auto photon_f_pt_d_np_ntrk = photon_f_pt->sum(2, 2);
    auto photon_f_pt_d_perp_ntrk = photon_f_pt_d_np_ntrk->sum(0);
    auto photon_f_pt_d_near_ntrk = photon_f_pt_d_np_ntrk->sum(1);

    auto photon_f_pt_d_np_sumpt = photon_f_pt->sum(0, 0);
    auto photon_f_pt_d_perp_sumpt = photon_f_pt_d_np_sumpt->sum(0);
    auto photon_f_pt_d_near_sumpt = photon_f_pt_d_np_sumpt->sum(1);

    /* photon-jet momentum imbalance as function of ntrk */
    auto pjet_f_x_d_np_ntrk = pjet_f_x->sum(0, 2, 2);

    auto pjet_f_x_d_perp_ntrk = pjet_f_x_d_np_ntrk->sum(0);
    auto pjet_f_x_d_near_ntrk = pjet_f_x_d_np_ntrk->sum(1);

    /* mixed events */
    auto mix_pjet_f_x_d_np_ntrk = mix_pjet_f_x->sum(0, 2, 2);

    auto mix_pjet_f_x_d_perp_ntrk = mix_pjet_f_x_d_np_ntrk->sum(0);
    auto mix_pjet_f_x_d_near_ntrk = mix_pjet_f_x_d_np_ntrk->sum(1);

    /* photon-jet momentum imbalance as function of sumpt */
    auto pjet_f_x_d_np_sumpt = pjet_f_x->sum(0, 0, 0);

    auto pjet_f_x_d_perp_sumpt = pjet_f_x_d_np_sumpt->sum(0);
    auto pjet_f_x_d_near_sumpt = pjet_f_x_d_np_sumpt->sum(1);

    /* mixed events */
    auto mix_pjet_f_x_d_np_sumpt = mix_pjet_f_x->sum(0, 0, 0);

    auto mix_pjet_f_x_d_perp_sumpt = mix_pjet_f_x_d_np_sumpt->sum(0);
    auto mix_pjet_f_x_d_near_sumpt = mix_pjet_f_x_d_np_sumpt->sum(1);

    /* subtract histograms */
    /* *pjet_f_x_d_perp_ntrk -= *mix_pjet_f_x_d_perp_ntrk; */
    /* *pjet_f_x_d_near_ntrk -= *mix_pjet_f_x_d_near_ntrk; */
    /* *pjet_f_x_d_perp_sumpt -= *mix_pjet_f_x_d_perp_sumpt; */
    /* *pjet_f_x_d_near_sumpt -= *mix_pjet_f_x_d_near_sumpt; */

    ntrk_f_pt->divide((*photon_f_pt_incl)[0]);
    sumpt_f_pt->divide((*photon_f_pt_incl)[0]);
    mix_ntrk_f_pt->divide((*photon_f_pt_incl)[0]);
    mix_sumpt_f_pt->divide((*photon_f_pt_incl)[0]);

    evt_f_ntrk->divide((*nevt_d_photon_pt), x{1});
    evt_f_sumpt->divide((*nevt_d_photon_pt), x{1});
    mix_evt_f_ntrk->divide((*nevt_d_photon_pt), x{1});
    mix_evt_f_sumpt->divide((*nevt_d_photon_pt), x{1});

    /* normalise to number of photons (events) */
    normalise(nevt_d_perp_ntrk, photon_f_pt_d_perp_ntrk,
        pjet_f_x_d_perp_ntrk, mix_pjet_f_x_d_perp_ntrk);
    normalise(nevt_d_near_ntrk, photon_f_pt_d_near_ntrk,
        pjet_f_x_d_near_ntrk, mix_pjet_f_x_d_near_ntrk);
    normalise(nevt_d_perp_sumpt, photon_f_pt_d_perp_sumpt,
        pjet_f_x_d_perp_sumpt, mix_pjet_f_x_d_perp_sumpt);
    normalise(nevt_d_near_sumpt, photon_f_pt_d_near_sumpt,
        pjet_f_x_d_near_sumpt, mix_pjet_f_x_d_near_sumpt);

    /* scale by bin width */
    scale_bin_width(
        photon_f_pt_d_near_ntrk, photon_f_pt_d_perp_ntrk,
        photon_f_pt_d_near_sumpt, photon_f_pt_d_perp_sumpt,
        pjet_f_x_d_near_ntrk, mix_pjet_f_x_d_near_ntrk,
        pjet_f_x_d_perp_ntrk, mix_pjet_f_x_d_perp_ntrk,
        pjet_f_x_d_near_sumpt, mix_pjet_f_x_d_near_sumpt,
        pjet_f_x_d_perp_sumpt, mix_pjet_f_x_d_perp_sumpt,
        evt_f_ntrk, evt_f_sumpt, mix_evt_f_ntrk, mix_evt_f_sumpt,
        photon_f_pt_incl);

    printf("painting..\n");

    auto ntrk_selection = [&](int64_t index) {
        auto text = "N^{h^{#pm}}"s;
        if (index != 1)
            text = std::to_string((*intrk)[index - 1]) + " < "s + text;
        if (index != intrk->size())
            text = text + " < "s + std::to_string((*intrk)[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, text.data());
    };

    auto sumpt_selection = [&](int64_t index) {
        auto text = "#Sigmap_{T}^{h^{#pm}}"s;
        if (index != 1)
            text = std::to_string((*isumpt)[index - 1]) + " < "s + text;
        if (index != isumpt->size())
            text = text + " < "s + std::to_string((*isumpt)[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, text.data());
    };

    auto photon_pt_selection = [&](int64_t index) {
        auto text = std::to_string((*ipt)[index - 1]) + " < p_{T}^{#gamma} < "s
            + std::to_string((*ipt)[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, text.data());
    };

    auto hb = new pencil();
    hb->category("dphi", "perp", "near");
    hb->category("type", "raw", "mix");

    hb->alias("perp", "#pi/3 < #||{#Delta#phi^{#gammah}} < 2#pi/3");
    hb->alias("near", "0 < #||{#Delta#phi^{#gammah}} < #pi/3");
    /* hb->alias("raw", "sub."); */

    auto system = "pPb #sqrt{s_{NN}} = 8.16 TeV"s;

    auto c1 = new paper("pjet_x_d_sumpt", hb);
    apply_default_style(c1, system, 0., 1.2);
    c1->accessory(sumpt_selection);

    for (int64_t i = 0; i < isumpt->size(); ++i) {
        c1->add((*pjet_f_x_d_perp_sumpt)[i], "perp", "raw");
        c1->stack((*pjet_f_x_d_near_sumpt)[i], "near", "raw");
        c1->stack((*mix_pjet_f_x_d_perp_sumpt)[i], "perp", "mix");
        c1->stack((*mix_pjet_f_x_d_near_sumpt)[i], "near", "mix");
    }

    auto c2 = new paper("pjet_x_d_ntrk", hb);
    apply_default_style(c2, system, 0., 1.2);
    c2->accessory(ntrk_selection);

    for (int64_t i = 0; i < intrk->size(); ++i) {
        c2->add((*pjet_f_x_d_perp_ntrk)[i], "perp", "raw");
        c2->stack((*pjet_f_x_d_near_ntrk)[i], "near", "raw");
        c2->stack((*mix_pjet_f_x_d_perp_ntrk)[i], "perp", "mix");
        c2->stack((*mix_pjet_f_x_d_near_ntrk)[i], "near", "mix");
    }

    auto c3 = new paper("density_f_pt", hb);
    apply_default_style(c3, system, 0., 10.);
    c3->divide(2, 1);

    c3->add((*ntrk_f_pt)[0], "near", "raw");
    c3->stack((*ntrk_f_pt)[1], "perp", "raw");
    c3->stack((*mix_ntrk_f_pt)[0], "near", "mix");
    c3->stack((*mix_ntrk_f_pt)[1], "perp", "mix");

    c3->add((*sumpt_f_pt)[0], "near", "raw");
    c3->stack((*sumpt_f_pt)[1], "perp", "raw");
    c3->stack((*mix_sumpt_f_pt)[0], "near", "mix");
    c3->stack((*mix_sumpt_f_pt)[1], "perp", "mix");

    auto c4 = new paper("ntrk_d_pt", hb);
    apply_default_style(c4, system, 0., 0.4);
    c4->accessory(photon_pt_selection);

    for (int64_t i = 0; i < ipt->size() - 1; ++i) {
        c4->add((*evt_f_ntrk)[x{i, 0}], "near", "raw");
        c4->stack((*evt_f_ntrk)[x{i, 1}], "perp", "raw");
        c4->stack((*mix_evt_f_ntrk)[x{i, 0}], "near", "mix");
        c4->stack((*mix_evt_f_ntrk)[x{i, 1}], "perp", "mix");
    }

    auto c5 = new paper("sumpt_d_pt", hb);
    apply_default_style(c5, system, 0., 0.4);
    c5->accessory(photon_pt_selection);

    for (int64_t i = 0; i < ipt->size() - 1; ++i) {
        c5->add((*evt_f_sumpt)[x{i, 0}], "near", "raw");
        c5->stack((*evt_f_sumpt)[x{i, 1}], "perp", "raw");
        c5->stack((*mix_evt_f_sumpt)[x{i, 0}], "near", "mix");
        c5->stack((*mix_evt_f_sumpt)[x{i, 1}], "perp", "mix");
    }

    auto c6 = new paper("pt_d_sumpt", hb);
    apply_default_style(c6, system, 0., 0.12);
    c6->accessory(sumpt_selection);

    for (int64_t i = 0; i < isumpt->size(); ++i) {
        c6->add((*photon_f_pt_d_perp_sumpt)[i], "perp");
        c6->stack((*photon_f_pt_d_near_sumpt)[i], "near");
    }

    auto c7 = new paper("pt_d_ntrk", hb);
    apply_default_style(c7, system, 0., 0.12);
    c7->accessory(ntrk_selection);

    for (int64_t i = 0; i < intrk->size(); ++i) {
        c7->add((*photon_f_pt_d_perp_ntrk)[i], "perp");
        c7->stack((*photon_f_pt_d_near_ntrk)[i], "near");
    }

    hb->set_binary("type");
    hb->sketch();

    c1->draw("pdf");
    c2->draw("pdf");
    c3->draw("pdf");
    c4->draw("pdf");
    c5->draw("pdf");
    c6->draw("pdf");
    c7->draw("pdf");

    /* fout->Write("", TObject::kOverwrite); */

    (void)fout;

    printf("destroying objects..\n");

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return flatten(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
