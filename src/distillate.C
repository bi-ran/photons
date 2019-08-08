#include "../include/lambdas.h"

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

#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

int distillate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");

    auto rpt = conf->get<std::vector<float>>("pt_range");
    auto reta = conf->get<std::vector<float>>("eta_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto deta = conf->get<std::vector<float>>("eta_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    auto remove = conf->get<std::vector<int64_t>>("remove");
    auto csn = conf->get<std::vector<float>>("csn");

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load input */
    TFile* f = new TFile(input.data(), "read");
    auto scale = new history(f, tag + "_scale");

    /* prepare histograms */
    auto ival = std::make_shared<interval>(1, 0., 1.);

    auto ipt = std::make_shared<interval>(dpt);
    auto ieta = std::make_shared<interval>(deta);
    auto ihf = std::make_shared<interval>(dhf);

    auto hf_shape = x{ ihf->size() };
    auto pthf_shape = x{ ipt->size(), ihf->size() };
    auto etahf_shape = x{ ieta->size(), ihf->size() };

    /* fully differential (pt, eta, hf) */
    auto es = new history("es"s, "", ival, scale->shape());
    auto er = new history("er"s, "", ival, scale->shape());

    auto es_f_pt = std::make_unique<history>("es_f_pt"s,
        "reco p_{T}/gen p_{T}", "jet p_{T}", rpt, etahf_shape);
    auto er_f_pt = std::make_unique<history>("er_f_pt"s,
        "#sigma(p_{T})/p_{T}", "jet p_{T}", rpt, etahf_shape);

    /* differential in pt, hf */
    auto scale_dpthf = scale->sum(1);

    auto es_dpthf = new history("es_dpthf", "", ival, pthf_shape);
    auto er_dpthf = new history("er_dpthf", "", ival, pthf_shape);

    auto es_dhf_f_pt = std::make_unique<history>("es_dhf_f_pt"s,
        "reco p_{T}/gen p_{T}", "jet p_{T}", rpt, hf_shape);
    auto er_dhf_f_pt = std::make_unique<history>("er_dhf_f_pt"s,
        "#sigma(p_{T})/p_{T}", "jet p_{T}", rpt, hf_shape);

    /* differential in eta, hf */
    std::vector<int64_t> resize = {ipt->size() - 1, ieta->size(), ihf->size()};
    auto scale_detahf = scale->shrink("valid", resize, remove)->sum(0);

    auto es_detahf = new history("es_detahf", "", ival, etahf_shape);
    auto er_detahf = new history("er_detahf", "", ival, etahf_shape);

    auto es_dhf_f_eta = std::make_unique<history>("es_dhf_f_eta"s,
        "reco p_{T}/gen p_{T}", "jet #eta", reta, hf_shape);
    auto er_dhf_f_eta = std::make_unique<history>("er_dhf_f_eta"s,
        "#sigma(p_{T})/p_{T}", "jet #eta", reta, hf_shape);

    /* load fitting parameters */
    auto fl = new std::vector<float>*[ihf->size()];
    auto fh = new std::vector<float>*[ihf->size()];

    auto flp = new std::vector<float>[ihf->size()];
    auto fhp = new std::vector<float>[ihf->size()];
    auto fle = new std::vector<float>[ihf->size()];
    auto fhe = new std::vector<float>[ihf->size()];

    for (int64_t i = 0; i < ihf->size(); ++i) {
        auto hf_str = std::to_string(i);

        fl[i] = new std::vector<float>[ieta->size()];
        fh[i] = new std::vector<float>[ieta->size()];

        for (int64_t j = 0; j < ieta->size(); ++j) {
            auto eta_str = std::to_string(j);
            fl[i][j] = conf->get<std::vector<float>>(
                "fl_"s + hf_str + "_"s + eta_str);
            fh[i][j] = conf->get<std::vector<float>>(
                "fh_"s + hf_str + "_"s + eta_str);
        }

        flp[i] = conf->get<std::vector<float>>("flp_"s + hf_str);
        fhp[i] = conf->get<std::vector<float>>("fhp_"s + hf_str);
        fle[i] = conf->get<std::vector<float>>("fle_"s + hf_str);
        fhe[i] = conf->get<std::vector<float>>("fhe_"s + hf_str);
    }

    auto set_csn = [&](TF1* f, int64_t hf_x) {
        if (!heavyion || csn.empty()) {
            f->SetParameters(0.08, 0.32, 0.);
            return;
        }

        f->SetParameters(csn[0], csn[1], csn[2]);
        if (hf_x > 0) {
            f->FixParameter(0, csn[0]);
            f->FixParameter(1, csn[1]);
        }
    };

    /* info text */
    auto eta_info = [&](int64_t, int64_t eta_x, int64_t offset) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.1f < #eta < %.1f",
            (*ieta)[eta_x], (*ieta)[eta_x + 1]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75 - offset * 0.04, buffer);
    };

    auto pt_hf_selection = [&](int64_t index) {
        auto pt_x = scale_dpthf->indices_for(index - 1)[0];
        auto hf_x = scale_dpthf->indices_for(index - 1)[1];

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);

        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{jet} < %.0f",
            (*ipt)[pt_x], (*ipt)[pt_x + 1]);
        l->DrawLatexNDC(0.135, 0.75, buffer);

        sprintf(buffer, "%i - %i%%", dcent[hf_x + 1], dcent[hf_x]);
        l->DrawLatexNDC(0.135, 0.71, buffer);
    };

    auto eta_hf_selection = [&](int64_t index) {
        auto eta_x = scale_detahf->indices_for(index - 1)[0];
        auto hf_x = scale_detahf->indices_for(index - 1)[1];

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);

        char buffer[128] = { '\0' };
        sprintf(buffer, "%.1f < #eta < %.1f",
            (*ieta)[eta_x], (*ieta)[eta_x + 1]);
        l->DrawLatexNDC(0.135, 0.75, buffer);

        sprintf(buffer, "%i - %i%%", dcent[hf_x + 1], dcent[hf_x]);
        l->DrawLatexNDC(0.135, 0.71, buffer);
    };

    auto hf_selection = [&](int64_t index) {
        auto hf_x = (index - 1) % ihf->size();

        char buffer[128] = { '\0' };
        sprintf(buffer, "%i - %i%%", dcent[hf_x + 1], dcent[hf_x]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, buffer);
    };

    auto system_info = system + " #sqrt{s_{NN}} = 5.02 TeV";

    /* draw plots */
    auto hb = new pencil();
    hb->category("sample", "mc");

    hb->alias("mc", "AllQCDPhoton");

    auto c1 = new paper(tag + "_dpthf_jesr_fits", hb);
    apply_default_style(c1, system_info, 0., 1.);
    c1->format(simple_formatter);
    c1->accessory(pt_hf_selection);
    c1->divide(ipt->size(), -1);

    /* fit scale and resolution */
    scale_dpthf->apply([&](TH1* h, int64_t index) {
        auto indices = scale_dpthf->indices_for(index);
        auto pt_x = indices[0];
        auto hf_x = indices[1];

        auto label = "f_scale_dpthf_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", 0, 2);
        h->Fit(label.data(), "WLMQ", "", flp[hf_x][pt_x], fhp[hf_x][pt_x]);

        (*es_dpthf)[index]->SetBinContent(1, f->GetParameter(1));
        (*es_dpthf)[index]->SetBinError(1, f->GetParError(1));
        (*er_dpthf)[index]->SetBinContent(1, f->GetParameter(2));
        (*er_dpthf)[index]->SetBinError(1, f->GetParError(2));

        ++pt_x;

        (*es_dhf_f_pt)[hf_x]->SetBinContent(pt_x, f->GetParameter(1));
        (*es_dhf_f_pt)[hf_x]->SetBinError(pt_x, f->GetParError(1));
        (*er_dhf_f_pt)[hf_x]->SetBinContent(pt_x, f->GetParameter(2));
        (*er_dhf_f_pt)[hf_x]->SetBinError(pt_x, f->GetParError(2));

        c1->add(h, "mc");
    });

    auto c2 = new paper(tag + "_dhf_f_pt_jesr", hb);
    apply_default_style(c2, system_info, 0., 1.);
    c2->format(simple_formatter);
    c2->accessory(hf_selection);
    c2->divide(ihf->size(), -1);
    c2->set(paper::flags::logx);

    es_dhf_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(0.8, 1.5, "Y");

        auto label = "f_es_dhf_f_pt_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "[0]+[1]/x+[2]/(x*x)",
            rpt.front(), rpt.back());
        f->SetParameters(1.1, 1.2, 4.8);
        h->Fit(label.data(), "MEQ", "", 30, rpt.back());

        c2->add(h, "mc");
    });

    er_dhf_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(0., 1., "Y");

        auto label = "f_er_dhf_f_pt_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",
            rpt.front(), rpt.back());
        set_csn(f, index);
        h->Fit(label.data(), "MEQ", "", 30, rpt.back());

        csn[1] = f->GetParameter(1);
        csn[2] = f->GetParameter(2);

        printf("%i - %i%%: %.3f, %.3f, %.3f\n",
            dcent[index + 1], dcent[index], csn[0], csn[1], csn[2]);

        c2->add(h, "mc");
    });

    auto c3 = new paper(tag + "_detahf_jesr_fits", hb);
    apply_default_style(c3, system_info, 0., 1.);
    c3->format(simple_formatter);
    c3->accessory(eta_hf_selection);
    c3->divide(ieta->size(), -1);

    /* fit scale and resolution */
    scale_detahf->apply([&](TH1* h, int64_t index) {
        auto indices = scale_detahf->indices_for(index);
        auto eta_x = indices[0];
        auto hf_x = indices[1];

        auto label = "f_scale_detahf_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", 0, 2);
        h->Fit(label.data(), "WLMQ", "", fle[hf_x][eta_x], fhe[hf_x][eta_x]);

        (*es_detahf)[index]->SetBinContent(1, f->GetParameter(1));
        (*es_detahf)[index]->SetBinError(1, f->GetParError(1));
        (*er_detahf)[index]->SetBinContent(1, f->GetParameter(2));
        (*er_detahf)[index]->SetBinError(1, f->GetParError(2));

        ++eta_x;

        (*es_dhf_f_eta)[hf_x]->SetBinContent(eta_x, f->GetParameter(1));
        (*es_dhf_f_eta)[hf_x]->SetBinError(eta_x, f->GetParError(1));
        (*er_dhf_f_eta)[hf_x]->SetBinContent(eta_x, f->GetParameter(2));
        (*er_dhf_f_eta)[hf_x]->SetBinError(eta_x, f->GetParError(2));

        c3->add(h, "mc");
    });

    auto c4 = new paper(tag + "_dhf_f_eta_jesr", hb);
    apply_default_style(c4, system_info, 0., 1.);
    c4->format(simple_formatter);
    c4->accessory(hf_selection);
    c4->divide(ihf->size(), -1);

    es_dhf_f_eta->apply([&](TH1* h) {
        h->SetAxisRange(0.8, 1.5, "Y");
        c4->add(h, "mc"); });

    er_dhf_f_eta->apply([&](TH1* h) {
        h->SetAxisRange(0., 1., "Y");
        c4->add(h, "mc"); });

    auto c5 = std::vector<paper*>(ieta->size());
    auto c6 = std::vector<paper*>(ieta->size());

    for (int64_t i = 0; i < ieta->size(); ++i) {
        c5[i] = new paper(tag + "_jesr_fits_s" + std::to_string(i), hb);
        apply_default_style(c5[i], system_info, 0., 1.);
        c5[i]->format(simple_formatter);
        c5[i]->accessory(pt_hf_selection);
        c5[i]->accessory(std::bind(eta_info, _1, i, 2));
        c5[i]->divide(ipt->size(), -1);

        c6[i] = new paper(tag + "_f_pt_jesr_s" + std::to_string(i), hb);
        apply_default_style(c6[i], system_info, 0., 1.);
        c6[i]->format(simple_formatter);
        c6[i]->accessory(hf_selection);
        c6[i]->accessory(std::bind(eta_info, _1, i, 1));
        c6[i]->divide(ihf->size(), -1);
        c6[i]->set(paper::flags::logx);
    }

    /* fit scale and resolution */
    scale->apply([&](TH1* h, int64_t index) {
        auto indices = scale->indices_for(index);
        auto pt_x = indices[0];
        auto eta_x = indices[1];
        auto hf_x = indices[2];

        auto label = "f_scale_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", 0, 2);
        h->Fit(label.data(), "WLMQ", "",
            fl[hf_x][eta_x][pt_x], fh[hf_x][eta_x][pt_x]);

        (*es)[index]->SetBinContent(1, f->GetParameter(1));
        (*es)[index]->SetBinError(1, f->GetParError(1));
        (*er)[index]->SetBinContent(1, f->GetParameter(2));
        (*er)[index]->SetBinError(1, f->GetParError(2));

        ++pt_x;

        (*es_f_pt)[x{eta_x, hf_x}]->SetBinContent(pt_x, f->GetParameter(1));
        (*es_f_pt)[x{eta_x, hf_x}]->SetBinError(pt_x, f->GetParError(1));
        (*er_f_pt)[x{eta_x, hf_x}]->SetBinContent(pt_x, f->GetParameter(2));
        (*er_f_pt)[x{eta_x, hf_x}]->SetBinError(pt_x, f->GetParError(2));

        c5[eta_x]->add(h, "mc");
    });

    es_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(0.8, 1.5, "Y");

        auto label = "f_es_f_pt_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "[0]+[1]/x+[2]/(x*x)");
        f->SetParameters(1.1, 1.2, 4.8);
        h->Fit(label.data(), "MEQ", "", 30, rpt.back());

        auto eta_x = es_f_pt->indices_for(index)[0];
        c6[eta_x]->add(h, "mc");
    });

    er_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(0., 1., "Y");

        auto label = "f_er_f_pt_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))");
        set_csn(f, index);
        h->Fit(label.data(), "MEQ", "", 30, rpt.back());

        csn[1] = f->GetParameter(1);
        csn[2] = f->GetParameter(2);

        auto eta_x = er_f_pt->indices_for(index)[0];
        c6[eta_x]->add(h, "mc");
    });

    hb->sketch();

    for (auto const& p : { c1, c2, c3, c4 })
        p->draw("pdf");

    for (auto const& c : { c5, c6 })
        for (auto p : c) { p->draw("pdf"); }

    /* save output */
    TFile* fout = new TFile(output, "recreate");

    scale_dpthf->save(tag);
    scale_detahf->save(tag);

    es->save(tag);
    er->save(tag);
    es_f_pt->save(tag);
    er_f_pt->save(tag);

    es_dpthf->save(tag);
    er_dpthf->save(tag);
    es_dhf_f_pt->save(tag);
    er_dhf_f_pt->save(tag);

    es_detahf->save(tag);
    er_detahf->save(tag);
    es_dhf_f_eta->save(tag);
    er_dhf_f_eta->save(tag);

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return distillate(argv[1], argv[2]);

    return 0;
}
