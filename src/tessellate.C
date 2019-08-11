#include "../include/lambdas.h"
#include "../include/pjtree.h"
#include "../include/specifics.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/memory.h"

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

void fill_data(std::unique_ptr<memory>& see_iso,
               std::unique_ptr<memory>& see_noniso,
               std::shared_ptr<multival>& mpthf,
               TTree* t, pjtree* p, bool heavyion,
               float pt_min, float eta_max, float hovere_max,
               float iso_max, float noniso_min, float noniso_max,
               float hf_min) {
    printf("fill data\n");

    auto nentries = static_cast<int64_t>(t->GetEntries());
    for (int64_t i = 0; i < nentries; ++i) {
        if (i % 100000 == 0) { printf("%li/%li\n", i, nentries); }

        t->GetEntry(i);

        if (p->hiHF <= hf_min) { continue; }

        int64_t leading = -1;
        for (int64_t j = 0; j < p->nPho; ++j) {
            if ((*p->phoEt)[j] < pt_min) { continue; }
            if (std::abs((*p->phoSCEta)[j]) > eta_max) { continue; }
            if ((*p->phoHoverE)[j] > hovere_max) { continue; }

            if (heavyion && within_hem_failure_region(p, j)) { continue; }

            leading = j;
            break;
        }

        /* require leading photon */
        if (leading < 0) { continue; }

        /* isolation requirement */
        float isolation = (*p->pho_ecalClusterIsoR3)[leading]
            + (*p->pho_hcalRechitIsoR3)[leading]
            + (*p->pho_trackIsoR3PtCut20)[leading];

        if ((isolation > iso_max && isolation < noniso_min)
            || isolation > noniso_max) { continue; }

        auto const& see = isolation > iso_max ? see_noniso : see_iso;
        int64_t index = mpthf->index_for(v{(*p->phoEt)[leading], p->hiHF});
        (*see)[index]->Fill((*p->phoSigmaIEtaIEta_2012)[leading], p->weight);
    }

    printf("\n");
}

void fill_signal(std::unique_ptr<memory>& see,
                 std::shared_ptr<multival>& mpthf,
                 TTree* t, pjtree* p, bool heavyion,
                 float pt_min, float eta_max, float hovere_max,
                 float hf_min) {
    printf("fill signal\n");

    auto nentries = static_cast<int64_t>(t->GetEntries());
    for (int64_t i = 0; i < nentries; ++i) {
        if (i % 100000 == 0) { printf("%li/%li\n", i, nentries); }

        t->GetEntry(i);

        if (p->hiHF <= hf_min) { continue; }

        int64_t leading = -1;
        for (int64_t j = 0; j < p->nPho; ++j) {
            if ((*p->phoEt)[j] < pt_min) { continue; }
            if (std::abs((*p->phoSCEta)[j]) > eta_max) { continue; }
            if ((*p->phoHoverE)[j] > hovere_max) { continue; }

            if (heavyion && within_hem_failure_region(p, j)) { continue; }

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

auto fit_templates(TH1F* hdata, TH1F* hsig, TH1F* hbkg,
                   std::vector<float> const& range) {
    auto stub = "_"s + hdata->GetName();

    TH1F* tdata = (TH1F*)hdata->Clone(("t"s + stub).data());
    TH1F* tsig = (TH1F*)hsig->Clone(("t_s"s + stub).data());
    TH1F* tbkg = (TH1F*)hbkg->Clone(("t_b"s + stub).data());

    tsig->Scale(1. / tsig->Integral());
    tbkg->Scale(1. / tbkg->Integral());

    auto evaluate = [&](double* x, double* p) {
        float nsig = tsig->GetBinContent(tsig->FindBin(x[0]));
        float nbkg = tbkg->GetBinContent(tbkg->FindBin(x[0]));

        return p[0] * (nsig * p[1] + nbkg * (1 - p[1]));
    };

    TF1* f = new TF1(("f"s + stub).data(), evaluate, range[0], range[1], 2);
    f->SetParameters(tdata->Integral(), 0.8);
    f->SetParLimits(1, 0., 1.);

    tdata->Fit(("f"s + stub).data(), "L0Q", "", range[0], range[1]);
    tdata->Fit(("f"s + stub).data(), "L0Q", "", range[0], range[1]);
    tdata->Fit(("f"s + stub).data(), "L0QM", "", range[0], range[1]);

    auto p0 = f->GetParameter(0);
    auto p1 = f->GetParameter(1);

    auto p0_err = f->GetParError(0);
    auto p1_err = f->GetParError(1);

    auto chisq = f->GetChisquare();
    auto ndof = f->GetNDF();

    return std::make_tuple(p0, p1, p0_err, p1_err, chisq, ndof);
}

int tessellate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto data = conf->get<std::string>("data");
    auto signal = conf->get<std::string>("signal");
    auto heavyion = conf->get<bool>("heavyion");
    auto tag = conf->get<std::string>("tag");

    auto pt_min = conf->get<float>("pt_min");
    auto eta_max = conf->get<float>("eta_max");
    auto hovere_max = conf->get<float>("hovere_max");
    auto iso_max = conf->get<float>("iso_max");
    auto noniso_min = conf->get<float>("noniso_min");
    auto noniso_max = conf->get<float>("noniso_max");
    auto see_max = conf->get<float>("see_max");

    auto see_nbins = conf->get<int64_t>("see_nbins");
    auto see_low = conf->get<float>("see_low");
    auto see_high = conf->get<float>("see_high");
    auto rfit = conf->get<std::vector<float>>("rfit");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    /* exclude most peripheral events */
    auto hf_min = dhf.front();

    auto rsee = std::make_shared<interval>("#sigma_{#eta#eta}",
        see_nbins, see_low, see_high);

    auto ipt = std::make_shared<interval>(dpt);
    auto ihf = std::make_shared<interval>(dhf);

    auto mpthf = std::make_shared<multival>(dpt, dhf);

    auto see_data = std::make_unique<memory>("sigma_ieta_ieta_data"s,
        "counts", rsee, mpthf);
    auto see_sig = std::make_unique<memory>("sigma_ieta_ieta_sig"s,
        "counts", rsee, mpthf);
    auto see_bkg = std::make_unique<memory>("sigma_ieta_ieta_bkg"s,
        "counts", rsee, mpthf);

    TFile* fd = new TFile(data.data(), "read");
    TTree* td = (TTree*)fd->Get("pj");
    auto pd = new pjtree(false, false, td, { 1, 0, 1, 0, 0, 0 });

    TFile* fs = new TFile(signal.data(), "read");
    TTree* ts = (TTree*)fs->Get("pj");
    auto ps = new pjtree(true, false, ts, { 1, 1, 1, 0, 0, 0 });

    TFile* fout = new TFile(output, "recreate");

    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    fill_data(see_data, see_bkg, mpthf, td, pd, heavyion,
              pt_min, eta_max, hovere_max, iso_max, noniso_min, noniso_max,
              hf_min);

    fill_signal(see_sig, mpthf, ts, ps, heavyion,
                pt_min, eta_max, hovere_max,
                hf_min);

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

        sprintf(buffer, "%i - %i%%", dcent[hf_x + 1], dcent[hf_x]);
        text->DrawLatexNDC(0.54, 0.63, buffer);
    };

    auto c1 = new paper(tag + "_purity", hb);
    apply_default_style(c1, "PbPb #sqrt{s_{NN}} = 5.02 TeV"s, 0., 1.);
    c1->format(formatter);
    c1->accessory(info_text);
    c1->divide(ipt->size(), -1);

    printf("fit templates\n");

    std::vector<float> purities(mpthf->size(), 1.);
    for (int64_t i = 0; i < mpthf->size(); ++i) {
        auto res = fit_templates((*see_data)[i], (*see_sig)[i], (*see_bkg)[i],
                                 rfit);

        auto stub = "p_"s + (*see_data)[i]->GetName();
        auto pfit = (TH1F*)(*see_sig)[i]->Clone((stub + "f").data());
        auto pbkg = (TH1F*)(*see_bkg)[i]->Clone((stub + "b").data());

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

        auto ntot = pfit->Integral(1, pfit->FindBin(see_max));
        auto nbkg = pbkg->Integral(1, pbkg->FindBin(see_max));

        printf("purity: %.3f\n", 1. - nbkg / ntot);
        purities[i] = 1. - nbkg / ntot;
    }

    auto purity_text = [&](int64_t index) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "purity: %.3f", purities[index - 1]);

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

    /* save purities in history format */
    auto incl = std::make_shared<interval>(1, 0., 1.);
    auto purity = std::make_unique<memory>("pthf"s, "purity"s, incl, mpthf);

    for (int64_t i = 0; i < purity->size(); ++i)
        (*purity)[i]->SetBinContent(1, purities[i]);

    purity->save(tag);

    see_data->save(tag);
    see_sig->save(tag);
    see_bkg->save(tag);

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return tessellate(argv[1], argv[2]);

    return 0;
}
