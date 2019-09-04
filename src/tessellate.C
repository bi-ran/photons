#include "../include/lambdas.h"
#include "../include/pjtree.h"
#include "../include/specifics.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/memory.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/trunk.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TTree.h"

#include <string>
#include <tuple>
#include <vector>

using namespace std::literals::string_literals;

void fill_data(memory* see_iso, memory* see_noniso, multival* mpthf,
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

void fill_signal(memory* see, multival* mpthf,
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
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");

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

    auto rsee = new interval("#sigma_{#eta#eta}", see_nbins, see_low, see_high);

    auto incl = new interval(1, 0., 1.);
    auto ipt = new interval(dpt);

    auto mpthf = new multival(dpt, dhf);

    auto see_data = new memory("sigma_ieta_ieta_data"s, "counts", rsee, mpthf);
    auto see_sig = new memory("sigma_ieta_ieta_sig"s, "counts", rsee, mpthf);
    auto see_bkg = new memory("sigma_ieta_ieta_bkg"s, "counts", rsee, mpthf);

    auto purity = new memory("pthf"s, "purity"s, incl, mpthf);

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load inputs */
    TFile* fd = new TFile(data.data(), "read");
    TTree* td = (TTree*)fd->Get("pj");
    auto pd = new pjtree(false, false, td, { 1, 0, 1, 0, 0, 0 });

    TFile* fs = new TFile(signal.data(), "read");
    TTree* ts = (TTree*)fs->Get("pj");
    auto ps = new pjtree(true, false, ts, { 1, 1, 1, 0, 0, 0 });

    fill_data(see_data, see_bkg, mpthf, td, pd, heavyion,
              pt_min, eta_max, hovere_max, iso_max, noniso_min, noniso_max,
              hf_min);

    fill_signal(see_sig, mpthf, ts, ps, heavyion,
                pt_min, eta_max, hovere_max,
                hf_min);

    auto hb = new pencil();
    hb->category("type", "data", "sig", "bkg");

    hb->alias("sig", "PYTHIA8");
    hb->alias("bkg", "noniso. data");

    hb->style("sig", [](TH1* h) {
        h->SetLineColor(kPink);
        h->SetFillColor(kOrange + 7);
        h->SetFillStyle(3004);
    });

    hb->style("bkg", [](TH1* h) {
        h->SetLineColor(kGreen + 4);
        h->SetFillColor(kGreen + 1);
        h->SetFillStyle(3001);
    });

    std::function<void(int64_t, float)> pt_info = [&](int64_t x, float pos) {
        info_text(x, pos, "%.0f < p_{T}^{#gamma} < %.0f", dpt, false); };

    std::function<void(int64_t, float)> hf_info = [&](int64_t x, float pos) {
        info_text(x, pos, "%i - %i%%", dcent, true); };

    auto pthf_info = [&](int64_t index) {
        stack_text(index, 0.75, 0.04, mpthf, pt_info, hf_info); };

    auto purity_info = [&](int64_t index) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "purity: %.3f",
            (*purity)[index - 1]->GetBinContent(1));

        TLatex* text = new TLatex();
        text->SetTextFont(43);
        text->SetTextSize(12);
        text->DrawLatexNDC(0.54, 0.56, buffer);
    };

    auto c1 = new paper(tag + "_purity", hb);
    apply_style(c1, system + " #sqrt{s_{NN}} = 5.02 TeV"s);
    c1->accessory(pthf_info);
    c1->accessory(purity_info);
    c1->divide(ipt->size(), -1);

    printf("fit templates\n");

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

        (*purity)[i]->SetBinContent(1, 1. - nbkg / ntot);

        printf("purity: %.3f\n", (*purity)[i]->GetBinContent(1));
    }

    hb->sketch();
    c1->draw("pdf");

    /* save purities */
    in(output, [&]() {
        purity->save(tag);

        see_data->save(tag);
        see_sig->save(tag);
        see_bkg->save(tag);
    });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return tessellate(argv[1], argv[2]);

    return 0;
}
