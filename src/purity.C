#include "../include/lambdas.h"
#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TTree.h"

#include <string>
#include <tuple>
#include <vector>

using namespace std::literals::string_literals;

static bool in_hem_failure_region(float eta, float phi) {
    return (-3. < eta && eta < -1.3 && -1.57 < phi && phi < -0.87);
}

void fill_data(std::unique_ptr<history>& see_iso,
               std::unique_ptr<history>& see_noniso,
               std::shared_ptr<multival>& mpthf,
               TTree* t, pjtree* p,
               float pt_min, float eta_max, float hovere_max,
               float iso_max, float noniso_min) {
    printf("fill data\n");

    auto nentries = static_cast<int64_t>(t->GetEntries());
    for (int64_t i = 0; i < nentries; ++i) {
        if (i % 100000 == 0) { printf("%li/%li\n", i, nentries); }

        t->GetEntry(i);

        /* event selections */

        int64_t leading = -1;
        for (int64_t j = 0; j < p->nPho; ++j) {
            if ((*p->phoEt)[j] < pt_min) { continue; }
            if (std::abs((*p->phoSCEta)[j]) > eta_max) { continue; }
            if ((*p->phoHoverE)[j] > hovere_max) { continue; }

            if (in_hem_failure_region((*p->phoSCEta)[j], (*p->phoSCPhi)[j]))
                continue;

            leading = j;
            break;
        }

        /* require leading photon */
        if (leading < 0) { continue; }

        /* isolation requirement */
        float isolation = (*p->pho_ecalClusterIsoR4)[leading]
            + (*p->pho_hcalRechitIsoR4)[leading]
            + (*p->pho_trackIsoR4PtCut20)[leading];

        if (isolation > iso_max && isolation < noniso_min) { continue; }

        auto const& see = isolation > iso_max ? see_noniso : see_iso;
        int64_t index = mpthf->index_for(v{(*p->phoEt)[leading], p->hiHF});
        (*see)[index]->Fill((*p->phoSigmaIEtaIEta_2012)[leading], p->weight);
    }

    printf("\n");
}

void fill_signal(std::unique_ptr<history>& see,
                 std::shared_ptr<multival>& mpthf,
                 TTree* t, pjtree* p,
                 float pt_min, float eta_max, float hovere_max) {
    printf("fill signal\n");

    auto nentries = static_cast<int64_t>(t->GetEntries());
    for (int64_t i = 0; i < nentries; ++i) {
        if (i % 100000 == 0) { printf("%li/%li\n", i, nentries); }

        t->GetEntry(i);

        /* event selections */

        int64_t leading = -1;
        for (int64_t j = 0; j < p->nPho; ++j) {
            if ((*p->phoEt)[j] < pt_min) { continue; }
            if (std::abs((*p->phoSCEta)[j]) > eta_max) { continue; }
            if ((*p->phoHoverE)[j] > hovere_max) { continue; }

            leading = j;
            break;
        }

        /* require leading photon */
        if (leading < 0) { continue; }

        /* require gen-matching */
        int64_t gen_index = (*p->pho_genMatchedIndex)[leading];
        if (gen_index == -1) { continue; }

        auto pid = (*p->mcPID)[gen_index];
        auto mpid = (*p->mcMomPID)[gen_index];
        if (pid != 22 || (std::abs(mpid) > 22 && mpid != -999)) { continue; }

        /* gen isolation requirement */
        float isolation = (*p->mcCalIsoDR04)[gen_index];
        if (isolation > 5.) { continue; }

        int64_t index = mpthf->index_for(v{(*p->phoEt)[leading], p->hiHF});
        (*see)[index]->Fill((*p->phoSigmaIEtaIEta_2012)[leading], p->weight);
    }

    printf("\n");
}

auto fit_templates(TH1F* hdata, TH1F* hsig, TH1F* hbkg) {
    auto tag = "_"s + hdata->GetName();

    TH1F* tdata = (TH1F*)hdata->Clone(("t"s + tag).data());
    TH1F* tsig = (TH1F*)hsig->Clone(("t_s"s + tag).data());
    TH1F* tbkg = (TH1F*)hbkg->Clone(("t_b"s + tag).data());

    tsig->Scale(1. / tsig->Integral());
    tbkg->Scale(1. / tbkg->Integral());

    auto evaluate = [&](double* x, double* p) {
        float nsig = tsig->GetBinContent(tsig->FindBin(x[0]));
        float nbkg = tbkg->GetBinContent(tbkg->FindBin(x[0]));

        return p[0] * (nsig * p[1] + nbkg * (1 - p[1]));
    };

    auto range_low = tdata->GetBinLowEdge(1);
    auto range_high = tdata->GetBinLowEdge(tdata->GetNbinsX() + 1);

    TF1* f = new TF1(("f"s + tag).data(), evaluate, range_low, range_high, 2);
    f->SetParameters(tdata->Integral(), 0.8);
    f->SetParLimits(1, 0., 1.);

    tdata->Fit(("f"s + tag).data(), "L0Q", "", range_low, range_high);
    tdata->Fit(("f"s + tag).data(), "L0Q", "", range_low, range_high);
    tdata->Fit(("f"s + tag).data(), "L0QM", "", range_low, range_high);

    auto p0 = f->GetParameter(0);
    auto p1 = f->GetParameter(1);

    auto p0_err = f->GetParError(0);
    auto p1_err = f->GetParError(1);

    auto chisq = f->GetChisquare();
    auto ndof = f->GetNDF();

    return std::make_tuple(p0, p1, p0_err, p1_err, chisq, ndof);
}

int purity(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto data = conf->get<std::string>("data");
    auto signal = conf->get<std::string>("signal");

    auto pt_min = conf->get<float>("pt_min");
    auto eta_max = conf->get<float>("eta_max");
    auto hovere_max = conf->get<float>("hovere_max");
    auto iso_max = conf->get<float>("iso_max");
    auto noniso_min = conf->get<float>("noniso_min");

    auto see_nbins = conf->get<int64_t>("see_nbins");
    auto see_min = conf->get<float>("see_min");
    auto see_max = conf->get<float>("see_max");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");

    auto rsee = std::make_shared<interval>(
        see_nbins, see_min, see_max, "#sigma_{#eta#eta}");

    auto ipt = std::make_shared<interval>(dpt);
    auto ihf = std::make_shared<interval>(dhf);

    auto mpthf = std::make_shared<multival>(dpt, dhf);

    auto see_data = std::make_unique<history>("sigma_ieta_ieta_data"s,
        "counts", rsee, mpthf);
    auto see_sig = std::make_unique<history>("sigma_ieta_ieta_sig"s,
        "counts", rsee, mpthf);
    auto see_bkg = std::make_unique<history>("sigma_ieta_ieta_bkg"s,
        "counts", rsee, mpthf);

    TFile* fd = new TFile(data.data(), "read");
    TTree* td = (TTree*)fd->Get("pj");
    auto pd = new pjtree(false, td);

    TFile* fs = new TFile(signal.data(), "read");
    TTree* ts = (TTree*)fs->Get("pj");
    auto ps = new pjtree(true, ts);

    TFile* fout = new TFile(output, "recreate");

    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    fill_data(see_data, see_bkg, mpthf, td, pd, pt_min, eta_max, hovere_max,
              iso_max, noniso_min);

    fill_signal(see_sig, mpthf, ts, ps, pt_min, eta_max, hovere_max);

    auto hb = new pencil();
    hb->category("type", "data", "sig", "bkg");

    hb->alias("sig", "PYTHIA8 + HYDJET");
    hb->alias("bkg", "noniso. data");

    auto formatter = [](TH1* obj) {
        obj->SetStats(0);
        obj->SetMarkerSize(0.84);
        obj->GetXaxis()->CenterTitle();
        obj->GetYaxis()->CenterTitle();
    };

    auto info_text = [&](int64_t index) {
        TLatex* text = new TLatex();
        text->SetTextFont(43);
        text->SetTextSize(12);

        auto pt_x = (index - 1) % ipt->size();
        auto hf_x = (index - 1) / ipt->size();

        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{#gamma} < %0.f",
            (*ipt)[pt_x], (*ipt)[pt_x + 1]);

        text->DrawLatexNDC(0.54, 0.67, buffer);

        std::vector<std::string> bins = { "90"s, "50"s, "30"s, "10"s, "0"s };
        sprintf(buffer, "%s - %s%%", bins[hf_x + 1].data(), bins[hf_x].data());

        text->DrawLatexNDC(0.54, 0.63, buffer);
    };

    auto c1 = new paper("purity", hb);
    apply_default_style(c1, "PbPb #sqrt{s_{NN}} = 5.02 TeV"s, 0., 1.);
    c1->format(formatter);
    c1->accessory(info_text);
    c1->divide(ipt->size(), -1);

    printf("fit templates\n");

    std::vector<float> purities = { 0 };

    for (int64_t i = ipt->size(); i < mpthf->size(); ++i) {
        auto res = fit_templates((*see_data)[i], (*see_sig)[i], (*see_bkg)[i]);

        auto tag = "p_"s + (*see_data)[i]->GetName();
        auto pfit = (TH1F*)(*see_sig)[i]->Clone((tag + "f").data());
        auto pbkg = (TH1F*)(*see_bkg)[i]->Clone((tag + "b").data());

        auto entries = std::get<0>(res);
        auto fraction = std::get<1>(res);

        pfit->Scale(entries * fraction / pfit->Integral());
        pbkg->Scale(entries * (1. - fraction) / pbkg->Integral());

        pfit->Add(pbkg);

        c1->add((*see_data)[i], "data");
        c1->stack(pfit, "sig");
        c1->stack(pbkg, "bkg");

        c1->adjust(pfit, "hist f", "lf");
        c1->adjust(pbkg, "hist f", "lf");

        auto ntot = pfit->Integral(1, pfit->FindBin(0.01));
        auto nbkg = pbkg->Integral(1, pbkg->FindBin(0.01));

        printf("purity: %.3f\n", 1. - nbkg / ntot);
        purities.push_back(1. - nbkg / ntot);
    }

    auto purity_text = [&](int64_t index) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "purity: %.3f", purities[index]);

        TLatex* text = new TLatex();
        text->SetTextFont(43);
        text->SetTextSize(12);
        text->DrawLatexNDC(0.54, 0.56, buffer);
    };

    c1->accessory(purity_text);

    auto sig_style = [](TH1* h) {
        h->SetLineColor(kPink);
        h->SetFillColor(kOrange + 7);
        h->SetFillStyle(3004);
    };

    auto bkg_style = [](TH1* h) {
        h->SetLineColor(kGreen + 4);
        h->SetFillColor(kGreen + 1);
        h->SetFillStyle(3001);
    };

    hb->style("sig", sig_style);
    hb->style("bkg", bkg_style);
    hb->sketch();

    c1->draw("pdf");

    fout->Write("", TObject::kOverwrite);
    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return purity(argv[1], argv[2]);

    return 0;
}
