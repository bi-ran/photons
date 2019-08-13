#include "../include/lambdas.h"
#include "../include/specifics.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"

#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"

#include <memory>
#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

template <typename... T>
void scale_bin_width(std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        obj->Scale(1., "width"); }), 0)... };
}

template <typename... T>
void normalise_to_unity(std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        obj->Scale(1. / obj->Integral("width")); }), 0)... };
}

template <typename... T>
void scale_ia_bin_width(std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        for (int64_t i = 1; i <= obj->GetNbinsX(); ++i) {
            auto width = revert_radian(obj->GetBinWidth(i));
            obj->SetBinContent(i, obj->GetBinContent(i) / width);
        }
    }), 0)... };
}

template <typename... T>
void normalise_ia_to_unity(std::unique_ptr<T>&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        auto width = revert_radian(obj->GetBinLowEdge(1)
            - obj->GetBinLowEdge(obj->GetNbinsX() + 1));
        obj->Scale(1. / (obj->Integral() * std::abs(width)));
    }), 0)... };
}

int accumulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto label = conf->get<std::string>("label");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto rppt = conf->get<std::vector<float>>("ppt_range");
    auto rdphi = conf->get<std::vector<float>>("dphi_range");
    auto rx = conf->get<std::vector<float>>("x_range");
    auto rdr = conf->get<std::vector<float>>("dr_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dx = conf->get<std::vector<float>>("x_diff");

    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    /* convert to integral angle units (cast to double) */
    convert_in_place_pi(rdphi);

    auto ipt = std::make_shared<interval>(dpt);
    auto ihf = std::make_shared<interval>(dhf);
    auto ix = std::make_shared<interval>(dx);

    /* open input files */
    TFile* f = new TFile(input.data(), "read");

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    TFile* fout = new TFile(output, "recreate");

    /* load histograms */
    auto nevt = std::make_unique<history>(
        f, label + "_raw_nevt"s);

    auto pjet_es_f_dphi = std::make_unique<history>(
        f, label + "_raw_sub_pjet_es_f_dphi"s);
    auto pjet_wta_f_dphi = std::make_unique<history>(
        f, label + "_raw_sub_pjet_wta_f_dphi"s);
    auto pjet_f_x = std::make_unique<history>(
        f, label + "_raw_sub_pjet_f_x"s);
    auto pjet_f_ddr = std::make_unique<history>(
        f, label + "_raw_sub_pjet_f_ddr"s);

    /* rescale by number of signal photons (events) */
    pjet_es_f_dphi->multiply(*nevt);
    pjet_wta_f_dphi->multiply(*nevt);
    pjet_f_x->multiply(*nevt);
    pjet_f_ddr->multiply(*nevt);

    /* discard overflow photon pt bin */
    auto discard = [](std::unique_ptr<history>& h, int64_t axis) {
        auto shape = h->shape();
        shape[axis] = shape[axis] - 1;
        h = h->shrink("s", shape, std::vector<int64_t>(h->dims(), 0));
    };

    discard(nevt, 0);
    discard(pjet_es_f_dphi, 0);
    discard(pjet_wta_f_dphi, 0);
    discard(pjet_f_x, 0);
    discard(pjet_f_ddr, 0);

    /* integrate histograms */
    auto nevt_d_pt = nevt->sum(1);
    auto nevt_d_hf = nevt->sum(0);

    auto pjet_f_ddr_d_pthf = pjet_f_ddr->sum(2);

    auto pjet_es_f_dphi_d_pt = pjet_es_f_dphi->sum(1);
    auto pjet_wta_f_dphi_d_pt = pjet_wta_f_dphi->sum(1);
    auto pjet_f_x_d_pt = pjet_f_x->sum(1);
    auto pjet_f_x_d_hf = pjet_f_x->sum(0);
    auto pjet_f_ddr_d_pt = pjet_f_ddr_d_pthf->sum(1);
    auto pjet_f_ddr_d_hf = pjet_f_ddr_d_pthf->sum(0);

    /* normalise by number of signal photons (events) */
    pjet_es_f_dphi->divide(*nevt);
    pjet_wta_f_dphi->divide(*nevt);
    pjet_f_ddr_d_pthf->divide(*nevt);
    pjet_f_x->divide(*nevt);

    pjet_es_f_dphi_d_pt->divide(*nevt_d_pt);
    pjet_wta_f_dphi_d_pt->divide(*nevt_d_pt);
    pjet_f_x_d_pt->divide(*nevt_d_pt);
    pjet_f_x_d_hf->divide(*nevt_d_hf);
    pjet_f_ddr_d_pt->divide(*nevt_d_pt);
    pjet_f_ddr_d_hf->divide(*nevt_d_hf);

    /* scale by bin width */
    scale_bin_width(
        pjet_f_x,
        pjet_f_ddr_d_pthf,
        pjet_f_x_d_pt,
        pjet_f_x_d_hf,
        pjet_f_ddr_d_pt,
        pjet_f_ddr_d_hf);

    scale_ia_bin_width(
        pjet_es_f_dphi,
        pjet_wta_f_dphi,
        pjet_es_f_dphi_d_pt,
        pjet_wta_f_dphi_d_pt);

    /* normalise to unity */
    normalise_to_unity(
        pjet_f_ddr_d_pthf,
        pjet_f_ddr_d_pt,
        pjet_f_ddr_d_hf);

    normalise_ia_to_unity(
        pjet_es_f_dphi,
        pjet_wta_f_dphi,
        pjet_es_f_dphi_d_pt,
        pjet_wta_f_dphi_d_pt);

    /* save histograms */
    nevt->save(tag);

    pjet_es_f_dphi->save(tag);
    pjet_wta_f_dphi->save(tag);
    pjet_f_x->save(tag);
    pjet_f_ddr_d_pthf->save(tag);

    pjet_es_f_dphi_d_pt->save(tag);
    pjet_wta_f_dphi_d_pt->save(tag);
    pjet_f_x_d_pt->save(tag);
    pjet_f_x_d_hf->save(tag);
    pjet_f_ddr_d_pt->save(tag);
    pjet_f_ddr_d_hf->save(tag);

    /* draw plots */
    printf("painting..\n");

    auto pt_selection = [&](int64_t index, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{#gamma} < %.0f",
            (*ipt)[index - 1], (*ipt)[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto hf_selection = [&](int64_t index, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%i - %i%%", dcent[index], dcent[index - 1]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto info_text = [&](int64_t index) {
        auto indices = nevt->indices_for(index - 1);
        auto pt_x = indices[0];
        auto hf_x = indices[1];

        pt_selection(pt_x + 1, 0.75);
        hf_selection(hf_x + 1, 0.71);
    };

    auto hb = new pencil();
    hb->category("system", "PbPb", "pp");
    hb->category("axis", "na", "wta", "es");

    hb->alias("na", "");
    hb->alias("es", "E-Scheme");
    hb->alias("wta", "WTA");

    hb->ditto("es", "na");

    auto collisions = system + " #sqrt{s_{NN}} = 5.02 TeV"s;

    auto c1 = new paper(tag + "_dphi_d_pt", hb);
    apply_default_style(c1, collisions, -0.04, 0.24);
    c1->accessory(std::bind(pt_selection, _1, 0.75));
    c1->accessory(std::bind(line_at, _1, 0.f, rdphi[0], rdphi[1]));
    c1->divide(-1, 1);

    nevt_d_pt->apply([&](TH1*, int64_t index) {
        c1->add((*pjet_es_f_dphi_d_pt)[index], system, "es");
        c1->stack((*pjet_wta_f_dphi_d_pt)[index], system, "wta");
    });

    auto c2 = new paper(tag + "_ddr_d_pt", hb);
    apply_default_style(c2, collisions, -2., 27.);
    c2->accessory(std::bind(pt_selection, _1, 0.75));
    c2->accessory(std::bind(line_at, _1, 0.f, rdr[0], rdr[1]));
    c2->divide(-1, 1);

    pjet_f_ddr_d_pt->apply([&](TH1* h) { c2->add(h, system); });

    auto c3 = new paper(tag + "_ddr_d_hf", hb);
    apply_default_style(c3, collisions, -2., 27.);
    c3->accessory(std::bind(hf_selection, _1, 0.75));
    c3->accessory(std::bind(line_at, _1, 0.f, rdr[0], rdr[1]));
    c3->divide(-1, 1);

    pjet_f_ddr_d_hf->apply([&](TH1* h) { c3->add(h, system); });

    auto c4 = new paper(tag + "_x_d_pt", hb);
    apply_default_style(c4, collisions, -0.1, 1.5);
    c4->accessory(std::bind(pt_selection, _1, 0.75));
    c4->accessory(std::bind(line_at, _1, 0.f, rx[0], rx[1]));
    c4->divide(-1, 1);

    pjet_f_x_d_pt->apply([&](TH1* h) { c4->add(h, system); });

    auto c5 = new paper(tag + "_x_d_hf", hb);
    apply_default_style(c5, collisions, -0.1, 1.5);
    c5->accessory(std::bind(hf_selection, _1, 0.75));
    c5->accessory(std::bind(line_at, _1, 0.f, rx[0], rx[1]));
    c5->divide(-1, 1);

    pjet_f_x_d_hf->apply([&](TH1* h) { c5->add(h, system); });

    auto c6 = new paper(tag + "_dphi_d_pthf", hb);
    apply_default_style(c6, collisions, -0.04, 0.24);
    c6->accessory(info_text);
    c6->accessory(std::bind(line_at, _1, 0.f, rdphi[0], rdphi[1]));
    c6->divide(-1, ihf->size());

    nevt->apply([&](TH1*, int64_t index) {
        c6->add((*pjet_es_f_dphi)[index], system, "es");
        c6->stack((*pjet_wta_f_dphi)[index], system, "wta");
    });

    auto c7 = new paper(tag + "_x_d_pthf", hb);
    apply_default_style(c7, collisions, -0.1, 1.5);
    c7->accessory(info_text);
    c7->accessory(std::bind(line_at, _1, 0.f, rx[0], rx[1]));
    c7->divide(-1, ihf->size());

    pjet_f_x->apply([&](TH1* h) { c7->add(h, system); });

    auto c8 = new paper(tag + "_ddr_d_pthf", hb);
    apply_default_style(c8, collisions, -2., 27.);
    c8->accessory(info_text);
    c8->accessory(std::bind(line_at, _1, 0.f, rdr[0], rdr[1]));
    c8->divide(-1, ihf->size());

    pjet_f_ddr_d_pthf->apply([&](TH1* h) { c8->add(h, system); });

    hb->set_binary("system");
    hb->sketch();

    for (auto c : { c1, c2, c3, c4, c5, c6, c7, c8 })
        c->draw("pdf");

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return accumulate(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
