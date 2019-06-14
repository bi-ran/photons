#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"

#include <memory>
#include <string>
#include <vector>

#include "../include/differential_histograms.h"
#include "../include/integral_angles.h"
#include "../include/interval.h"
#include "../include/multival.h"

#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

#define FP_TH1_FILL (int (TH1::*)(double))&TH1::Fill
#define FP_TH1_FILLW (int (TH1::*)(double, double))&TH1::Fill
#define FP_TH1_GETBC (double (TH1::*)(int) const)&TH1::GetBinContent
#define FP_TH1_SETBC (void (TH1::*)(int, double))&TH1::SetBinContent
#define FP_TH1_DRAW (void (TH1::*)(char const*))&TH1::Draw

using namespace std::literals::string_literals;

using dhist = differential_histograms;

void fill_tracks(pjtree* pjt, float trk_pt_min, float trk_eta_abs,
                 int64_t photon_phi, int64_t photon_pt_x,
                 std::shared_ptr<interval>& idphi,
                 std::shared_ptr<multival>& mdphi,
                 std::unique_ptr<differential_histograms>& ntrk,
                 std::unique_ptr<differential_histograms>& sumpt,
                 std::unique_ptr<differential_histograms>& trk_f_dphi,
                 std::unique_ptr<differential_histograms>& trk_f_pt,
                 std::unique_ptr<differential_histograms>& evt_f_ntrk,
                 std::unique_ptr<differential_histograms>& evt_f_sumpt) {
    for (int64_t j = 0; j < pjt->nTrk; ++j) {
        if ((*pjt->trkPt)[j] < trk_pt_min) { continue; }
        if (std::abs((*pjt->trkEta)[j]) > trk_eta_abs) { continue; }

        /* track quality selections */

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

    for (int64_t j = 0; j < mdphi->size(); ++j) {
        double evt_ntrk = (*ntrk)(j, FP_TH1_GETBC, 1);
        double evt_sumpt = (*sumpt)(j, FP_TH1_GETBC, 1);

        (*evt_f_ntrk)(x{photon_pt_x, j}, FP_TH1_FILL, evt_ntrk);
        (*evt_f_sumpt)(x{photon_pt_x, j}, FP_TH1_FILL, evt_sumpt);
    }
}

void fill_jets(pjtree* pjt, float jet_pt_min, float jet_eta_abs,
               double photon_leading_pt, int64_t photon_phi, int64_t photon_pt_x,
               std::unique_ptr<differential_histograms>& ntrk,
               std::unique_ptr<differential_histograms>& sumpt,
               std::shared_ptr<interval>& intrk,
               std::shared_ptr<interval>& isumpt,
               double norm,
               std::unique_ptr<differential_histograms>& nevt,
               std::unique_ptr<differential_histograms>& pjet_f_dphi,
               std::unique_ptr<differential_histograms>& pjet_f_jetpt,
               std::unique_ptr<differential_histograms>& pjet_f_x) {
    int64_t near_ntrk_x = intrk->index_for((*ntrk)(0, FP_TH1_GETBC, 1) / norm);
    int64_t perp_ntrk_x = intrk->index_for((*ntrk)(1, FP_TH1_GETBC, 1) / norm);

    int64_t near_sumpt_x = isumpt->index_for((*sumpt)(0, FP_TH1_GETBC, 1) / norm);
    int64_t perp_sumpt_x = isumpt->index_for((*sumpt)(1, FP_TH1_GETBC, 1) / norm);

    auto index = nevt->index_for(
        x{photon_pt_x, near_ntrk_x, perp_ntrk_x, near_sumpt_x, perp_sumpt_x});

    (*nevt)(index, FP_TH1_FILL, 1.);

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

void normalise(std::unique_ptr<differential_histograms>& norm,
               std::unique_ptr<differential_histograms>& last) {
    /* assume object, normalisation have equal shapes */
    for (int64_t j = 0; j < norm->size(); ++j) {
        auto count = (*norm)[j]->GetBinContent(1);
        if (count != 0) { (*last)[j]->Scale(1. / count); }
    }
}

template <typename... T>
void normalise(std::unique_ptr<differential_histograms>& norm,
               std::unique_ptr<differential_histograms>& first,
               std::unique_ptr<T>&... others) {
    normalise(norm, first);
    normalise(norm, others...);
}

int flatten(char const* config, char const* output) {
    printf("load config options..\n");

    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto mix = conf->get<std::string>("mix");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto mc_branches = conf->get<bool>("mc_branches");

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

    /* convert to integral angle units (cast to double) */
    for (auto& i : rdphi) { i = convert_degrees(i); }
    for (auto& i : ddphi) { i = convert_degrees(i); }

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
    auto mptntrk2sumpt2 = std::make_shared<multival>(
        dpt, dntrk, dntrk, dsumpt, dsumpt);

    auto photon_f_pt = std::make_unique<dhist>("photon_f_pt"s, rpt, mincl);

    auto ntrk = std::make_unique<dhist>("ntrk"s, incl, mdphi);
    auto sumpt = std::make_unique<dhist>("sumpt"s, incl, mdphi);
    auto mix_ntrk = std::make_unique<dhist>("ntrk_mix"s, incl, mdphi);
    auto mix_sumpt = std::make_unique<dhist>("sumpt_mix"s, incl, mdphi);

    auto evt_f_ntrk = std::make_unique<dhist>("evt_f_ntrk"s, rntrk, mptdphi);
    auto evt_f_sumpt = std::make_unique<dhist>("evt_f_sumpt"s, rsumpt, mptdphi);
    auto mix_evt_f_ntrk = std::make_unique<dhist>("mix_evt_f_ntrk"s, rntrk, mptdphi);
    auto mix_evt_f_sumpt = std::make_unique<dhist>("mix_evt_f_sumpt"s, rsumpt, mptdphi);

    auto trk_f_dphi = std::make_unique<dhist>("trk_f_dphi"s, rdphi, mpt);
    auto trk_f_pt = std::make_unique<dhist>("trk_f_pt"s, rtrkpt, mptdphi);
    auto mix_trk_f_dphi = std::make_unique<dhist>("mix_trk_f_dphi"s, rdphi, mpt);
    auto mix_trk_f_pt = std::make_unique<dhist>("mix_trk_f_pt"s, rtrkpt, mptdphi);

    auto nevt = std::make_unique<dhist>("nevt"s, incl, mptntrk2sumpt2);
    auto nmix = std::make_unique<dhist>("nmix"s, incl, mptntrk2sumpt2);

    auto pjet_f_dphi = std::make_unique<dhist>(
        "pjet_f_dphi"s, rdphi, mptntrk2sumpt2);
    auto pjet_f_jetpt = std::make_unique<dhist>(
        "pjet_f_jetpt"s, rpt, mptntrk2sumpt2);
    auto pjet_f_x = std::make_unique<dhist>(
        "pjet_f_x"s, rx, mptntrk2sumpt2);

    auto mix_pjet_f_dphi = std::make_unique<dhist>(
        "mix_pjet_f_dphi"s, rdphi, mptntrk2sumpt2);
    auto mix_pjet_f_jetpt = std::make_unique<dhist>(
        "mix_pjet_f_jetpt"s, rpt, mptntrk2sumpt2);
    auto mix_pjet_f_x = std::make_unique<dhist>(
        "mix_pjet_f_x"s, rx, mptntrk2sumpt2);

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
        if (i % 1000 == 0) { printf("entry: %li/%li\n", i, nentries); }

        t->GetEntry(i);

        /* event selections */

        int64_t photon_leading = -1;
        for (int64_t j = 0; j < pjt->nPho; ++j) {
            if ((*pjt->phoEt)[j] < photon_pt_min) { continue; }
            if (std::abs((*pjt->phoEta)[j]) > photon_eta_abs) { continue; }
            if ((*pjt->phoHoverE)[j] > 0.1) { continue; }
            if ((*pjt->phoSigmaIEtaIEta_2012)[j] > 0.01) { continue; }

            /* isolation requirement */

            photon_leading = j;
            break;
        }

        /* require leading photon */
        if (photon_leading < 0) { continue; }

        double photon_leading_pt = (*pjt->phoEt)[photon_leading];
        int64_t photon_pt_x = ipt->index_for(photon_leading_pt);

        (*photon_f_pt)(0, FP_TH1_FILL, photon_leading_pt);

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
                    photon_phi, photon_pt_x,
                    idphi, mdphi,
                    ntrk, sumpt,
                    trk_f_dphi, trk_f_pt,
                    evt_f_ntrk, evt_f_sumpt);

        double norm = 2. * trk_eta_abs * 2. * M_PI / 3.;

        fill_jets(pjt, jet_pt_min, jet_eta_abs,
                  photon_leading_pt, photon_phi, photon_pt_x,
                  ntrk, sumpt, intrk, isumpt, norm,
                  nevt, pjet_f_dphi, pjet_f_jetpt, pjet_f_x);

        /* mixing events in minimum bias */
        for (int64_t k = 0; k < 100; ++k) {
            tm->GetEntry(m);

            fill_tracks(pjtm, trk_pt_min, trk_eta_abs,
                        photon_phi, photon_pt_x,
                        idphi, mdphi,
                        mix_ntrk, mix_sumpt,
                        mix_trk_f_dphi, mix_trk_f_pt,
                        mix_evt_f_ntrk, mix_evt_f_sumpt);

            fill_jets(pjtm, jet_pt_min, jet_eta_abs,
                      photon_leading_pt, photon_phi, photon_pt_x,
                      mix_ntrk, mix_sumpt, intrk, isumpt, norm,
                      nmix, mix_pjet_f_dphi, mix_pjet_f_jetpt, mix_pjet_f_x);

            m = (m + 1) % mentries;
        }
    }

    /* normalise histograms */
    normalise(nmix, mix_pjet_f_dphi, mix_pjet_f_jetpt, mix_pjet_f_x);

    /* integrate histograms */
    /* photon (event) count */
    auto nevt_d_near_perp_sumpt = nevt->sum(0)->sum(0)->sum(0);

    auto nevt_d_perp_sumpt = nevt_d_near_perp_sumpt->sum(0);
    auto nevt_d_near_sumpt = nevt_d_near_perp_sumpt->sum(1);

    /* photon-jet momentum imbalance as function of perp sumpt */
    auto pjet_f_x_d_near_perp_sumpt = pjet_f_x
        ->sum(0)    /* photon pt */
        ->sum(0)    /* near_ntrk */
        ->sum(0);   /* perp_ntrk */

    auto pjet_f_x_d_perp_sumpt = pjet_f_x_d_near_perp_sumpt
        ->sum(0);   /* near_sumpt */
    auto pjet_f_x_d_near_sumpt = pjet_f_x_d_near_perp_sumpt
        ->sum(1);   /* perp_sumpt */

    /* mixed events */
    auto mix_pjet_f_x_d_near_perp_sumpt = mix_pjet_f_x
        ->sum(0)    /* photon pt */
        ->sum(0)    /* near_ntrk */
        ->sum(0);   /* perp_ntrk */

    auto mix_pjet_f_x_d_perp_sumpt = mix_pjet_f_x_d_near_perp_sumpt
        ->sum(0);   /* near_sumpt */
    auto mix_pjet_f_x_d_near_sumpt = mix_pjet_f_x_d_near_perp_sumpt
        ->sum(1);   /* perp_sumpt */

    /* subtract histograms */
    *pjet_f_x_d_perp_sumpt -= *mix_pjet_f_x_d_perp_sumpt;
    *pjet_f_x_d_near_sumpt -= *mix_pjet_f_x_d_near_sumpt;

    /* normalise to number of photons (events) */
    normalise(nevt_d_perp_sumpt, pjet_f_x_d_perp_sumpt);
    normalise(nevt_d_near_sumpt, pjet_f_x_d_near_sumpt);
    normalise(nevt_d_perp_sumpt, mix_pjet_f_x_d_perp_sumpt);
    normalise(nevt_d_near_sumpt, mix_pjet_f_x_d_near_sumpt);

    printf("painting..\n");

    TCanvas* c1 = new TCanvas("c1", "", 800, 1200);

    c1->Divide(2, 3);

    for (int64_t i = 0; i < isumpt->size(); ++i) {
        c1->cd(i + 1);

        (*pjet_f_x_d_perp_sumpt)[i]->SetStats(0);
        (*pjet_f_x_d_perp_sumpt)[i]->SetMarkerStyle(21);
        (*pjet_f_x_d_perp_sumpt)[i]->SetAxisRange(0, 0.12, "Y");
        (*pjet_f_x_d_perp_sumpt)(i, FP_TH1_DRAW, "p e");

        (*pjet_f_x_d_near_sumpt)[i]->SetStats(0);
        (*pjet_f_x_d_near_sumpt)[i]->SetLineColor(2);
        (*pjet_f_x_d_near_sumpt)[i]->SetMarkerColor(2);
        (*pjet_f_x_d_near_sumpt)[i]->SetMarkerStyle(20);
        (*pjet_f_x_d_near_sumpt)(i, FP_TH1_DRAW, "same p e");

        (*mix_pjet_f_x_d_perp_sumpt)[i]->SetStats(0);
        (*mix_pjet_f_x_d_perp_sumpt)[i]->SetMarkerStyle(25);
        (*mix_pjet_f_x_d_perp_sumpt)(i, FP_TH1_DRAW, "same p e");

        (*mix_pjet_f_x_d_near_sumpt)[i]->SetStats(0);
        (*mix_pjet_f_x_d_near_sumpt)[i]->SetLineColor(2);
        (*mix_pjet_f_x_d_near_sumpt)[i]->SetMarkerColor(2);
        (*mix_pjet_f_x_d_near_sumpt)[i]->SetMarkerStyle(24);
        (*mix_pjet_f_x_d_near_sumpt)(i, FP_TH1_DRAW, "same p e");
    }

    c1->SaveAs("c1.pdf");

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
