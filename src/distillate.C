#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/trunk.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

int distillate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto object = conf->get<std::string>("object");
    auto label = conf->get<std::string>("label");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");
    auto fit = conf->get<bool>("fit");

    auto rpt = conf->get<std::vector<float>>("pt_range");
    auto reta = conf->get<std::vector<float>>("eta_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto deta = conf->get<std::vector<float>>("eta_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    auto remove = conf->get<std::vector<int64_t>>("remove");
    auto csn = conf->get<std::vector<float>>("csn");

    auto s_range = conf->get<std::vector<float>>("s_range");
    auto s_lines = conf->get<std::vector<float>>("s_lines");
    auto r_range = conf->get<std::vector<float>>("r_range");

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load input */
    TFile* f = new TFile(input.data(), "read");
    auto obj = new history(f, tag + "_" + object);

    /* prepare histograms */
    auto ival = std::make_shared<interval>(1, 0., 1.);

    auto ipt = std::make_shared<interval>(dpt);
    auto ieta = std::make_shared<interval>(deta);
    auto ihf = std::make_shared<interval>(dhf);

    auto hf_shape = x{ ihf->size() };
    auto pthf_shape = x{ ipt->size(), ihf->size() };
    auto etahf_shape = x{ ieta->size(), ihf->size() };

    auto title = "#sigma("s + label + ")";

    /* fully differential (pt, eta, hf) */
    auto s = new history("s"s, "", ival, obj->shape());
    auto r = new history("r"s, "", ival, obj->shape());

    auto s_f_pt = std::make_unique<history>("s_f_pt"s,
        label.data(), "jet p_{T}", rpt, etahf_shape);
    auto r_f_pt = std::make_unique<history>("r_f_pt"s,
        title.data(), "jet p_{T}", rpt, etahf_shape);

    /* differential in pt, hf */
    auto obj_dpthf = obj->sum(1);

    auto s_dpthf = new history("s_dpthf", "", ival, pthf_shape);
    auto r_dpthf = new history("r_dpthf", "", ival, pthf_shape);

    auto s_dhf_f_pt = std::make_unique<history>("s_dhf_f_pt"s,
        label.data(), "jet p_{T}", rpt, hf_shape);
    auto r_dhf_f_pt = std::make_unique<history>("r_dhf_f_pt"s,
        title.data(), "jet p_{T}", rpt, hf_shape);

    /* differential in eta, hf */
    std::vector<int64_t> resize = {ipt->size() - 1, ieta->size(), ihf->size()};
    auto obj_detahf = obj->shrink("valid", resize, remove)->sum(0);

    auto s_detahf = new history("s_detahf", "", ival, etahf_shape);
    auto r_detahf = new history("r_detahf", "", ival, etahf_shape);

    auto s_dhf_f_eta = std::make_unique<history>("s_dhf_f_eta"s,
        label.data(), "jet #eta", reta, hf_shape);
    auto r_dhf_f_eta = std::make_unique<history>("r_dhf_f_eta"s,
        title.data(), "jet #eta", reta, hf_shape);

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
        } else {
            f->SetParameters(csn[0], csn[1], csn[2]);
            if (hf_x > 0) {
                f->FixParameter(0, csn[0]);
                f->FixParameter(1, csn[1]);
            }
        }
    };

    /* info text */
    auto eta_info = [&](int64_t eta_x, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.1f < #eta < %.1f",
            (*ieta)[eta_x - 1], (*ieta)[eta_x]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto pt_info = [&](int64_t pt_x, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{jet} < %.0f",
            (*ipt)[pt_x - 1], (*ipt)[pt_x]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto hf_info = [&](int64_t hf_x, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%i - %i%%", dcent[hf_x], dcent[hf_x - 1]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto pthf_info = [&](int64_t index) {
        auto indices = obj_dpthf->indices_for(index - 1);

        pt_info(indices[0] + 1, 0.75);
        hf_info(indices[1] + 1, 0.71);
    };

    auto etahf_info = [&](int64_t index) {
        auto indices = obj_detahf->indices_for(index - 1);

        eta_info(indices[0] + 1, 0.75);
        hf_info(indices[1] + 1, 0.71);
    };

    auto tag_object = tag + "_" + object;
    auto system_info = system + " #sqrt{s_{NN}} = 5.02 TeV";

    /* draw plots */
    auto hb = new pencil();
    hb->category("sample", "mc");

    hb->alias("mc", "AllQCDPhoton");

    auto c1 = new paper(tag_object + "_dpthf_sr_fits", hb);
    apply_style(c1, system_info);
    c1->accessory(pthf_info);
    c1->divide(ipt->size(), -1);

    /* fit obj and resolution */
    obj_dpthf->apply([&](TH1* h, int64_t index) {
        auto indices = obj_dpthf->indices_for(index);
        auto pt_x = indices[0];
        auto hf_x = indices[1];

        auto label = "f_obj_dpthf_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", 0, 2);
        h->Fit(label.data(), "WLMQ", "", flp[hf_x][pt_x], fhp[hf_x][pt_x]);

        (*s_dpthf)[index]->SetBinContent(1, f->GetParameter(1));
        (*s_dpthf)[index]->SetBinError(1, f->GetParError(1));
        (*r_dpthf)[index]->SetBinContent(1, f->GetParameter(2));
        (*r_dpthf)[index]->SetBinError(1, f->GetParError(2));

        ++pt_x;

        (*s_dhf_f_pt)[hf_x]->SetBinContent(pt_x, f->GetParameter(1));
        (*s_dhf_f_pt)[hf_x]->SetBinError(pt_x, f->GetParError(1));
        (*r_dhf_f_pt)[hf_x]->SetBinContent(pt_x, f->GetParameter(2));
        (*r_dhf_f_pt)[hf_x]->SetBinError(pt_x, f->GetParError(2));

        c1->add(h, "mc");
    });

    auto c2 = new paper(tag_object + "_dhf_f_pt_s", hb);
    apply_style(c2, system_info);
    c2->accessory(std::bind(hf_info, _1, 0.75));
    c2->divide(ihf->size(), -1);
    c2->set(paper::flags::logx);

    s_dhf_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(s_range[0], s_range[1], "Y");

        if (fit) {
            auto label = "f_s_dhf_f_pt_"s + std::to_string(index);
            TF1* f = new TF1(label.data(), "[0]+[1]/x+[2]/(x*x)");
            f->SetParameters(1.1, 1.2, 4.8);
            h->Fit(label.data(), "MEQ", "", 30, rpt.back());
        }

        c2->add(h, "mc");
    });

    auto c3 = new paper(tag_object + "_dhf_f_pt_r", hb);
    apply_style(c3, system_info);
    c3->accessory(std::bind(hf_info, _1, 0.75));
    c3->divide(ihf->size(), -1);
    c3->set(paper::flags::logx);

    r_dhf_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(r_range[0], r_range[1], "Y");

        if (fit) {
            auto label = "f_r_dhf_f_pt_"s + std::to_string(index);
            TF1* f = new TF1(label.data(),
                "sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))");
            set_csn(f, index);
            h->Fit(label.data(), "MEQ", "", 30, rpt.back());

            csn[0] = f->GetParameter(0);
            csn[1] = f->GetParameter(1);
            csn[2] = f->GetParameter(2);

            printf("%i - %i%%: %.3f, %.3f, %.3f\n",
                dcent[index + 1], dcent[index], csn[0], csn[1], csn[2]);
        }

        c3->add(h, "mc");
    });

    auto c4 = new paper(tag_object + "_detahf_sr_fits", hb);
    apply_style(c4, system_info);
    c4->accessory(etahf_info);
    c4->divide(ieta->size(), -1);

    /* fit mean and resolution */
    obj_detahf->apply([&](TH1* h, int64_t index) {
        auto indices = obj_detahf->indices_for(index);
        auto eta_x = indices[0];
        auto hf_x = indices[1];

        auto label = "f_obj_detahf_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", 0, 2);
        h->Fit(label.data(), "WLMQ", "", fle[hf_x][eta_x], fhe[hf_x][eta_x]);

        (*s_detahf)[index]->SetBinContent(1, f->GetParameter(1));
        (*s_detahf)[index]->SetBinError(1, f->GetParError(1));
        (*r_detahf)[index]->SetBinContent(1, f->GetParameter(2));
        (*r_detahf)[index]->SetBinError(1, f->GetParError(2));

        ++eta_x;

        (*s_dhf_f_eta)[hf_x]->SetBinContent(eta_x, f->GetParameter(1));
        (*s_dhf_f_eta)[hf_x]->SetBinError(eta_x, f->GetParError(1));
        (*r_dhf_f_eta)[hf_x]->SetBinContent(eta_x, f->GetParameter(2));
        (*r_dhf_f_eta)[hf_x]->SetBinError(eta_x, f->GetParError(2));

        c4->add(h, "mc");
    });

    auto c5 = new paper(tag_object + "_dhf_f_eta_s", hb);
    apply_style(c5, system_info);
    c5->accessory(std::bind(hf_info, _1, 0.75));
    c5->divide(ihf->size(), -1);

    s_dhf_f_eta->apply([&](TH1* h) {
        h->SetAxisRange(s_range[0], s_range[1], "Y");
        c5->add(h, "mc"); });

    auto c6 = new paper(tag_object + "_dhf_f_eta_r", hb);
    apply_style(c6, system_info);
    c6->accessory(std::bind(hf_info, _1, 0.75));
    c6->divide(ihf->size(), -1);

    r_dhf_f_eta->apply([&](TH1* h) {
        h->SetAxisRange(r_range[0], r_range[1], "Y");
        c6->add(h, "mc"); });

    auto c7 = std::vector<paper*>(ieta->size());
    auto c8 = std::vector<paper*>(ieta->size());
    auto c9 = std::vector<paper*>(ieta->size());

    for (int64_t i = 0; i < ieta->size(); ++i) {
        c7[i] = new paper(tag_object + "_sr_fits_s" + std::to_string(i), hb);
        apply_style(c7[i], system_info);
        c7[i]->accessory(pthf_info);
        c7[i]->ornaments(std::bind(eta_info, i + 1, 0.67));
        c7[i]->divide(ipt->size(), -1);

        c8[i] = new paper(tag_object + "_f_pt_s_s" + std::to_string(i), hb);
        apply_style(c8[i], system_info);
        c8[i]->accessory(std::bind(hf_info, _1, 0.75));
        c8[i]->ornaments(std::bind(eta_info, i + 1, 0.71));
        c8[i]->divide(ihf->size(), -1);
        c8[i]->set(paper::flags::logx);

        c9[i] = new paper(tag_object + "_f_pt_r_s" + std::to_string(i), hb);
        apply_style(c9[i], system_info);
        c9[i]->accessory(std::bind(hf_info, _1, 0.75));
        c9[i]->ornaments(std::bind(eta_info, i + 1, 0.71));
        c9[i]->divide(ihf->size(), -1);
        c9[i]->set(paper::flags::logx);
    }

    /* fit mean and resolution */
    obj->apply([&](TH1* h, int64_t index) {
        auto indices = obj->indices_for(index);
        auto pt_x = indices[0];
        auto eta_x = indices[1];
        auto hf_x = indices[2];

        auto label = "f_obj_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", 0, 2);
        h->Fit(label.data(), "WLMQ", "",
            fl[hf_x][eta_x][pt_x], fh[hf_x][eta_x][pt_x]);

        (*s)[index]->SetBinContent(1, f->GetParameter(1));
        (*s)[index]->SetBinError(1, f->GetParError(1));
        (*r)[index]->SetBinContent(1, f->GetParameter(2));
        (*r)[index]->SetBinError(1, f->GetParError(2));

        ++pt_x;

        (*s_f_pt)[x{eta_x, hf_x}]->SetBinContent(pt_x, f->GetParameter(1));
        (*s_f_pt)[x{eta_x, hf_x}]->SetBinError(pt_x, f->GetParError(1));
        (*r_f_pt)[x{eta_x, hf_x}]->SetBinContent(pt_x, f->GetParameter(2));
        (*r_f_pt)[x{eta_x, hf_x}]->SetBinError(pt_x, f->GetParError(2));

        c7[eta_x]->add(h, "mc");
    });

    s_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(s_range[0], s_range[1], "Y");

        if (fit) {
            auto label = "f_s_f_pt_"s + std::to_string(index);
            TF1* f = new TF1(label.data(), "[0]+[1]/x+[2]/(x*x)");
            f->SetParameters(1.1, 1.2, 4.8);
            h->Fit(label.data(), "MEQ", "", 30, rpt.back());
        }

        auto eta_x = s_f_pt->indices_for(index)[0];
        c8[eta_x]->add(h, "mc");
    });

    r_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(r_range[0], r_range[1], "Y");

        if (fit) {
            auto label = "f_r_f_pt_"s + std::to_string(index);
            TF1* f = new TF1(label.data(),
                "sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))");
            set_csn(f, index);
            h->Fit(label.data(), "MEQ", "", 30, rpt.back());

            csn[1] = f->GetParameter(1);
            csn[2] = f->GetParameter(2);
        }

        auto eta_x = r_f_pt->indices_for(index)[0];
        c9[eta_x]->add(h, "mc");
    });

    hb->sketch();

    for (auto const& p : { c1, c2, c3, c4, c5, c6 })
        p->draw("pdf");
    for (auto const& c : { c7, c8, c9 })
        for (auto p : c)
            p->draw("pdf");

    /* save output */
    in(output, [&]() {
        obj_dpthf->save(tag_object);
        obj_detahf->save(tag_object);

        s->save(tag_object);
        r->save(tag_object);
        s_f_pt->save(tag_object);
        r_f_pt->save(tag_object);

        s_dpthf->save(tag_object);
        r_dpthf->save(tag_object);
        s_dhf_f_pt->save(tag_object);
        r_dhf_f_pt->save(tag_object);

        s_detahf->save(tag_object);
        r_detahf->save(tag_object);
        s_dhf_f_eta->save(tag_object);
        r_dhf_f_eta->save(tag_object);
    });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return distillate(argv[1], argv[2]);

    return 0;
}
