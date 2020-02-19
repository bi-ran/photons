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
using namespace std::placeholders;

void fill_data(memory<TH1F>* see_iso, memory<TH1F>* see_noniso,
               multival* mpthf, TTree* t, pjtree* p, bool heavyion,
               float pt_min, float eta_max, float hovere_max, float hf_min,
               float iso_max, float noniso_min, float noniso_max) {
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

void fill_signal(memory<TH1F>* see, memory<TH1F>* sfrac,
                 multival* mpthf, TTree* t, pjtree* p, bool heavyion,
                 float pt_min, float eta_max, float hovere_max, float hf_min,
                 float iso_max, float noniso_min, float noniso_max) {
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

        /* isolation requirement */
        float recoiso = (*p->pho_ecalClusterIsoR3)[leading]
            + (*p->pho_hcalRechitIsoR3)[leading]
            + (*p->pho_trackIsoR3PtCut20)[leading];

        if ((recoiso > iso_max && recoiso < noniso_min)
            || recoiso > noniso_max) { continue; }

        (*sfrac)[index]->Fill(recoiso > iso_max ? 1.5 : 0.5, p->weight);
    }

    printf("\n");
}

auto fit_templates(TH1F* hdata, TH1F* hsig, TH1F* hbkg,
                   TH1F* hpurity, TH1F* hsfrac,
                   std::vector<float> const& range) {
    auto stub = "_"s + hdata->GetName();

    TH1F* tdata = (TH1F*)hdata->Clone(("t"s + stub).data());
    TH1F* tsig = (TH1F*)hsig->Clone(("t_s"s + stub).data());
    TH1F* tbkg = (TH1F*)hbkg->Clone(("t_b"s + stub).data());

    tsig->Scale(1. / tsig->Integral());

    auto signorm = tbkg->Integral() * hpurity->GetBinContent(1)
        * hsfrac->GetBinContent(2) / hsfrac->GetBinContent(1);

    tbkg->Add(tsig, -signorm);
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

    return std::make_tuple(f, tsig, tbkg);
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

    auto rsee = conf->get<std::vector<float>>("rsee");
    auto rfit = conf->get<std::vector<float>>("rfit");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    /* exclude most peripheral events */
    auto hf_min = dhf.front();

    /* prepare histograms */
    auto mpthf = new multival(dpt, dhf);

    auto incl = new interval(""s, 1, 0., 1.);
    auto itwo = new interval(""s, 2, 0., 2.);
    auto isee = new interval("#sigma_{#eta#eta}"s, rsee[0], rsee[1], rsee[2]);

    auto fincl = std::bind(&interval::book<TH1F>, incl, _1, _2, _3);
    auto ftwo = std::bind(&interval::book<TH1F>, itwo, _1, _2, _3);
    auto fsee = std::bind(&interval::book<TH1F>, isee, _1, _2, _3);

    auto ipt = new interval(dpt);

    auto see_data = new memory<TH1F>("see_data"s, "counts", fsee, mpthf);
    auto see_sig = new memory<TH1F>("see_sig"s, "counts", fsee, mpthf);
    auto see_bkg = new memory<TH1F>("see_bkg"s, "counts", fsee, mpthf);

    auto purity = new memory<TH1F>("pthf"s, "purity"s, fincl, mpthf);
    auto sfrac = new memory<TH1F>("sfrac"s, "counts"s, ftwo, mpthf);

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
              pt_min, eta_max, hovere_max, hf_min,
              iso_max, noniso_min, noniso_max);

    fill_signal(see_sig, sfrac, mpthf, ts, ps, heavyion,
                pt_min, eta_max, hovere_max, hf_min,
                iso_max, noniso_min, noniso_max);

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
        std::vector<float> iters = {0};

        TF1* f = nullptr;
        TH1F* pfit = nullptr;
        TH1F* pbkg = nullptr;

        do {
            auto res = fit_templates(
                (*see_data)[i], (*see_sig)[i], (*see_bkg)[i],
                (*purity)[i], (*sfrac)[i], rfit);

            f = std::get<0>(res);
            pfit = std::get<1>(res);
            pbkg = std::get<2>(res);

            auto frac = f->GetParameter(1);
            auto chi2 = f->GetChisquare();
            auto ndof = f->GetNDF();

            auto nsig = pfit->Integral(1, pfit->FindBin(see_max))
                * frac / pfit->Integral();
            auto nbkg = pbkg->Integral(1, pbkg->FindBin(see_max))
                * (1. - frac) / pbkg->Integral();

            auto val = nsig / (nsig + nbkg);
            (*purity)[i]->SetBinContent(1, val);
            printf("  purity: %.3f, chi2/ndof: %.0f/%i\n", val, chi2, ndof);

            iters.push_back(val);
        } while (std::abs(iters.end()[-1] - iters.end()[-2]) > 0.001);

        printf("iterations: %zu\n", iters.size());

        auto entries = f->GetParameter(0);
        auto fraction = f->GetParameter(1);

        pfit->Scale(entries * fraction / pfit->Integral());
        pbkg->Scale(entries * (1. - fraction) / pbkg->Integral());

        pfit->Add(pbkg);

        c1->add((*see_data)[i], "data");
        c1->stack(pfit, "sig");
        c1->stack(pbkg, "bkg");

        c1->adjust(pfit, "hist f", "lf");
        c1->adjust(pbkg, "hist f", "lf");
    }

    hb->sketch();
    c1->draw("pdf");

    /* save purities */
    in(output, [&]() {
        purity->save(tag);
        sfrac->save(tag);

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
