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

int distillate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto rpt = conf->get<std::vector<float>>("pt_range");
    auto reta = conf->get<std::vector<float>>("eta_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto deta = conf->get<std::vector<float>>("eta_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

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

    auto mhf = std::make_shared<multival>(dhf);

    /* differential in pt, hf */
    auto scale_d_pthf = scale->sum(1);

    auto es_d_pt = new history("es_d_pt"s, "", ival, scale_d_pthf->shape());
    auto er_d_pt = new history("er_d_pt"s, "", ival, scale_d_pthf->shape());

    auto es_f_pt = std::make_unique<history>("es_f_pt"s,
        "reco p_{T}/gen p_{T}", "jet p_{T}", rpt, mhf);
    auto er_f_pt = std::make_unique<history>("er_f_pt"s,
        "#sigma(p_{T})/p_{T}", "jet p_{T}", rpt, mhf);

    /* differential in eta, hf */

    /* remove first pt interval - hardcoded for [25, 30] */
    std::vector<int64_t> resize = {ipt->size() - 1, ieta->size(), ihf->size()};
    auto scale_d_etahf = scale->shrink("valid", resize, {1, 0, 0})->sum(0);

    auto es_d_eta = new history("es_d_eta"s, "", ival, scale_d_etahf->shape());
    auto er_d_eta = new history("er_d_eta"s, "", ival, scale_d_etahf->shape());

    auto es_f_eta = std::make_unique<history>("es_f_eta"s,
        "reco p_{T}/gen p_{T}", "jet #eta", reta, mhf);
    auto er_f_eta = std::make_unique<history>("er_f_eta"s,
        "#sigma(p_{T})/p_{T}", "jet #eta", reta, mhf);

    /* load fitting parameters */
    auto flp = new std::vector<float>[ihf->size()];
    auto fhp = new std::vector<float>[ihf->size()];
    auto fle = new std::vector<float>[ihf->size()];
    auto fhe = new std::vector<float>[ihf->size()];
    for (int64_t i = 0; i < ihf->size(); ++i) {
        auto index_str = std::to_string(i);
        flp[i] = conf->get<std::vector<float>>("flp_"s + index_str);
        fhp[i] = conf->get<std::vector<float>>("fhp_"s + index_str);
        fle[i] = conf->get<std::vector<float>>("fle_"s + index_str);
        fhe[i] = conf->get<std::vector<float>>("fhe_"s + index_str);
    }

    /* info text */
    auto pt_hf_selection = [&](int64_t index) {
        auto pt_x = scale_d_pthf->indices_for(index - 1)[0];
        auto hf_x = scale_d_pthf->indices_for(index - 1)[1];

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

    auto hf_selection = [&](int64_t index) {
        auto hf_x = (index - 1) % ihf->size();

        char buffer[128] = { '\0' };
        sprintf(buffer, "%i - %i%%", dcent[hf_x + 1], dcent[hf_x]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, buffer);
    };

    /* draw plots */
    auto hb = new pencil();
    hb->category("sample", "mc");

    hb->alias("mc", "AllQCDPhoton");

    auto c1 = new paper(tag + "_pt_jesr_fits", hb);
    apply_default_style(c1, system + " #sqrt{s_{NN}} = 5.02 TeV" , 0., 1.);
    c1->format(simple_formatter);
    c1->accessory(pt_hf_selection);
    c1->divide(ipt->size(), -1);

    /* fit scale and resolution */
    scale_d_pthf->apply([&](TH1* h, int64_t index) {
        auto indices = scale_d_pthf->indices_for(index);
        auto pt_x = indices[0];
        auto hf_x = indices[1];

        auto label = "f_pt_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", 0, 2);
        h->Fit(label.data(), "WLMQ", "", flp[hf_x][pt_x], fhp[hf_x][pt_x]);

        (*es_d_pt)[index]->SetBinContent(1, f->GetParameter(1));
        (*es_d_pt)[index]->SetBinError(1, f->GetParError(1));
        (*er_d_pt)[index]->SetBinContent(1, f->GetParameter(2));
        (*er_d_pt)[index]->SetBinError(1, f->GetParError(2));

        ++pt_x;

        (*es_f_pt)[hf_x]->SetBinContent(pt_x, f->GetParameter(1));
        (*es_f_pt)[hf_x]->SetBinError(pt_x, f->GetParError(1));
        (*er_f_pt)[hf_x]->SetBinContent(pt_x, f->GetParameter(2));
        (*er_f_pt)[hf_x]->SetBinError(pt_x, f->GetParError(2));

        c1->add(h, "mc");
    });

    auto c2 = new paper(tag + "_pt_jesr", hb);
    apply_default_style(c2, system + " #sqrt{s_{NN}} = 5.02 TeV" , 0., 1.);
    c2->format(simple_formatter);
    c2->accessory(hf_selection);
    c2->divide(ihf->size(), -1);
    c2->set(paper::flags::logx);

    es_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(0.8, 1.5, "Y");

        auto label = "f_es_pt_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "[0]+[1]/x+[2]/(x*x)",
            rpt.front(), rpt.back());
        f->SetParameters(1.1, 1.2, 4.8);
        h->Fit(label.data(), "MEQ", "", 30, rpt.back());

        c2->add(h, "mc");
    });

    er_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(0., 1., "Y");

        auto label = "f_er_pt_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",
            rpt.front(), rpt.back());
        f->SetParameters(0.1, 1.2, 4.8);
        h->Fit(label.data(), "MEQ", "", 30, rpt.back());

        c2->add(h, "mc");
    });

    auto c3 = new paper(tag + "_eta_jesr_fits", hb);
    apply_default_style(c3, system + " #sqrt{s_{NN}} = 5.02 TeV" , 0., 1.);
    c3->format(simple_formatter);
    c3->divide(ieta->size(), -1);

    /* fit scale and resolution */
    scale_d_etahf->apply([&](TH1* h, int64_t index) {
        auto indices = scale_d_etahf->indices_for(index);
        auto eta_x = indices[0];
        auto hf_x = indices[1];

        auto label = "f_eta_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", 0, 2);
        h->Fit(label.data(), "WLMQ", "", fle[hf_x][eta_x], fhe[hf_x][eta_x]);

        (*es_d_eta)[index]->SetBinContent(1, f->GetParameter(1));
        (*es_d_eta)[index]->SetBinError(1, f->GetParError(1));
        (*er_d_eta)[index]->SetBinContent(1, f->GetParameter(2));
        (*er_d_eta)[index]->SetBinError(1, f->GetParError(2));

        ++eta_x;

        (*es_f_eta)[hf_x]->SetBinContent(eta_x, f->GetParameter(1));
        (*es_f_eta)[hf_x]->SetBinError(eta_x, f->GetParError(1));
        (*er_f_eta)[hf_x]->SetBinContent(eta_x, f->GetParameter(2));
        (*er_f_eta)[hf_x]->SetBinError(eta_x, f->GetParError(2));

        c3->add(h, "mc");
    });

    auto c4 = new paper(tag + "_eta_jesr", hb);
    apply_default_style(c4, system + " #sqrt{s_{NN}} = 5.02 TeV" , 0., 1.);
    c4->format(simple_formatter);
    c4->accessory(hf_selection);
    c4->divide(ihf->size(), -1);

    es_f_eta->apply([&](TH1* h) {
        h->SetAxisRange(0.8, 1.5, "Y");
        c4->add(h, "mc"); });

    er_f_eta->apply([&](TH1* h) {
        h->SetAxisRange(0., 1., "Y");
        c4->add(h, "mc"); });

    hb->sketch();

    c1->draw("pdf");
    c2->draw("pdf");
    c3->draw("pdf");
    c4->draw("pdf");

    /* save output */
    TFile* fout = new TFile(output, "recreate");

    scale_d_pthf->save(tag);
    scale_d_etahf->save(tag);

    es_d_pt->save(tag);
    er_d_pt->save(tag);
    es_f_pt->save(tag);
    er_f_pt->save(tag);

    es_d_eta->save(tag);
    er_d_eta->save(tag);
    es_f_eta->save(tag);
    er_f_eta->save(tag);

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return distillate(argv[1], argv[2]);

    return 0;
}
