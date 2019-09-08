#include "../include/lambdas.h"
#include "../include/specifics.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"
#include "../git/tricks-and-treats/include/trunk.h"
#include "../git/tricks-and-treats/include/zip.h"

#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

template <typename... T>
void normalise_to_unity(T*&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        obj->Scale(1. / obj->Integral("width")); }), 0)... };
}

template <typename... T>
void normalise_ia_to_unity(T*&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        auto width = revert_radian(obj->GetBinLowEdge(1)
            - obj->GetBinLowEdge(obj->GetNbinsX() + 1));
        obj->Scale(1. / (obj->Integral() * std::abs(width)));
    }), 0)... };
}

template <typename... T>
void title(std::function<void(TH1*)> f, T*&... args) {
    (void)(int [sizeof...(T)]) { (args->apply(f), 0)... };
}

int accumulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto label = conf->get<std::string>("label");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto rjpt = conf->get<std::vector<float>>("jpt_range");
    auto rdphi = conf->get<std::vector<float>>("dphi_range");
    auto rx = conf->get<std::vector<float>>("x_range");
    auto rdr = conf->get<std::vector<float>>("dr_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");

    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    /* convert to integral angle units (cast to double) */
    convert_in_place_pi(rdphi);

    auto ihf = new interval(dhf);

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* open input files */
    TFile* f = new TFile(input.data(), "read");

    /* load histograms */
    auto nevt = new history<TH1F>(f, label + "_raw_nevt"s);

    auto pjet_es_f_dphi = new history<TH1F>(
        f, label + "_raw_sub_pjet_es_f_dphi"s);
    auto pjet_wta_f_dphi = new history<TH1F>(
        f, label + "_raw_sub_pjet_wta_f_dphi"s);
    auto pjet_f_x = new history<TH1F>(
        f, label + "_raw_sub_pjet_f_x"s);
    auto pjet_f_ddr = new history<TH1F>(
        f, label + "_raw_sub_pjet_f_ddr"s);
    auto pjet_f_jpt = new history<TH1F>(
        f, label + "_raw_sub_pjet_f_jpt"s);

    /* rescale by number of signal photons (events) */
    pjet_es_f_dphi->multiply(*nevt);
    pjet_wta_f_dphi->multiply(*nevt);
    pjet_f_x->multiply(*nevt);
    pjet_f_ddr->multiply(*nevt);
    pjet_f_jpt->multiply(*nevt);

    /* discard overflow photon pt bin */
    auto discard = [](history<TH1F>*& h, int64_t axis) {
        auto shape = h->shape();
        shape[axis] = shape[axis] - 1;
        h = h->shrink("s", shape, std::vector<int64_t>(h->dims(), 0));
    };

    discard(nevt, 0);
    discard(pjet_es_f_dphi, 0);
    discard(pjet_wta_f_dphi, 0);
    discard(pjet_f_x, 0);
    discard(pjet_f_ddr, 0);
    discard(pjet_f_jpt, 0);

    /* integrate histograms */
    auto nevt_d_pt = nevt->sum(1);
    auto nevt_d_hf = nevt->sum(0);

    auto pjet_es_f_dphi_d_pt = pjet_es_f_dphi->sum(1);
    auto pjet_es_f_dphi_d_hf = pjet_es_f_dphi->sum(0);
    auto pjet_wta_f_dphi_d_pt = pjet_wta_f_dphi->sum(1);
    auto pjet_wta_f_dphi_d_hf = pjet_wta_f_dphi->sum(0);
    auto pjet_f_x_d_pt = pjet_f_x->sum(1);
    auto pjet_f_x_d_hf = pjet_f_x->sum(0);
    auto pjet_f_ddr_d_pt = pjet_f_ddr->sum(1);
    auto pjet_f_ddr_d_hf = pjet_f_ddr->sum(0);
    auto pjet_f_jpt_d_pt = pjet_f_jpt->sum(1);
    auto pjet_f_jpt_d_hf = pjet_f_jpt->sum(0);

    /* normalise by number of signal photons (events) */
    pjet_es_f_dphi->divide(*nevt);
    pjet_wta_f_dphi->divide(*nevt);
    pjet_f_x->divide(*nevt);
    pjet_f_ddr->divide(*nevt);
    pjet_f_jpt->divide(*nevt);

    pjet_es_f_dphi_d_pt->divide(*nevt_d_pt);
    pjet_es_f_dphi_d_hf->divide(*nevt_d_hf);
    pjet_wta_f_dphi_d_pt->divide(*nevt_d_pt);
    pjet_wta_f_dphi_d_hf->divide(*nevt_d_hf);
    pjet_f_x_d_pt->divide(*nevt_d_pt);
    pjet_f_x_d_hf->divide(*nevt_d_hf);
    pjet_f_ddr_d_pt->divide(*nevt_d_pt);
    pjet_f_ddr_d_hf->divide(*nevt_d_hf);
    pjet_f_jpt_d_pt->divide(*nevt_d_pt);
    pjet_f_jpt_d_hf->divide(*nevt_d_hf);

    /* normalise to unity */
    normalise_to_unity(
        pjet_f_ddr,
        pjet_f_ddr_d_pt,
        pjet_f_ddr_d_hf);

    normalise_ia_to_unity(
        pjet_es_f_dphi,
        pjet_wta_f_dphi,
        pjet_es_f_dphi_d_pt,
        pjet_es_f_dphi_d_hf,
        pjet_wta_f_dphi_d_pt,
        pjet_wta_f_dphi_d_hf);

    title(std::bind(rename_axis, _1, "1/N^{#gammaj}dN/d#deltaj"),
        pjet_f_ddr,
        pjet_f_ddr_d_pt,
        pjet_f_ddr_d_hf);

    title(std::bind(rename_axis, _1, "1/N^{#gammaj}dN/d#Delta#phi^{#gammaj}"),
        pjet_es_f_dphi,
        pjet_wta_f_dphi,
        pjet_es_f_dphi_d_pt,
        pjet_es_f_dphi_d_hf,
        pjet_wta_f_dphi_d_pt,
        pjet_wta_f_dphi_d_hf);

    /* save histograms */
    in(output, [&]() {
        nevt->save(tag);

        pjet_es_f_dphi->save(tag);
        pjet_wta_f_dphi->save(tag);
        pjet_f_x->save(tag);
        pjet_f_ddr->save(tag);
        pjet_f_jpt->save(tag);

        pjet_es_f_dphi_d_pt->save(tag);
        pjet_es_f_dphi_d_hf->save(tag);
        pjet_wta_f_dphi_d_pt->save(tag);
        pjet_wta_f_dphi_d_hf->save(tag);
        pjet_f_x_d_pt->save(tag);
        pjet_f_x_d_hf->save(tag);
        pjet_f_ddr_d_pt->save(tag);
        pjet_f_ddr_d_hf->save(tag);
        pjet_f_jpt_d_pt->save(tag);
        pjet_f_jpt_d_hf->save(tag);
    });

    /* draw plots */
    printf("painting..\n");

    auto redraw_dphi_axis = [&](TH1* h, int64_t) {
        transform_axis(h, [](int64_t val) -> float {
            return std::abs(revert_radian(val)); }); };

    std::function<void(int64_t, float)> pt_info = [&](int64_t x, float pos) {
        info_text(x, pos, "%.0f < p_{T}^{#gamma} < %.0f", dpt, false); };

    std::function<void(int64_t, float)> hf_info = [&](int64_t x, float pos) {
        info_text(x, pos, "%i - %i%%", dcent, true); };

    auto pthf_info = [&](int64_t index) {
        stack_text(index, 0.75, 0.04, nevt, pt_info, hf_info); };

    auto hb = new pencil();
    hb->category("system", "PbPb", "pp");
    hb->category("axis", "na", "wta", "es");

    hb->alias("na", "");
    hb->alias("es", "E-Scheme");
    hb->alias("wta", "WTA");

    hb->ditto("es", "na");

    auto collisions = system + " #sqrt{s_{NN}} = 5.02 TeV"s;

    auto suffixes = { "d_pthf"s, "d_pt"s, "d_hf"s };
    auto texts = std::vector<std::function<void(int64_t)>> {
        pthf_info, std::bind(pt_info, _1, 0.75), std::bind(hf_info, _1, 0.75) };

    std::vector<paper*> c1(3, nullptr);
    std::vector<paper*> c2(3, nullptr);
    std::vector<paper*> c3(3, nullptr);
    std::vector<paper*> c4(3, nullptr);

    zip([&](paper*& c, int64_t rows, std::string const& suffix,
            std::function<void(int64_t)> text) {
        c = new paper(tag + "_dphi_" + suffix, hb);
        c->divide(-1, rows);
        c->accessory(text);

        apply_style(c, collisions, -0.04, 0.24);
        c->accessory(std::bind(line_at, _1, 0.f, rdphi[0], rdphi[1]));
        c->jewellery(redraw_dphi_axis);
    }, c1, x{ ihf->size(), 1L, 1L }, suffixes, texts);

    nevt->apply([&](TH1*, int64_t index) {
        c1[0]->add((*pjet_es_f_dphi)[index], system, "es");
        c1[0]->stack((*pjet_wta_f_dphi)[index], system, "wta");
    });

    nevt_d_pt->apply([&](TH1*, int64_t index) {
        c1[1]->add((*pjet_es_f_dphi_d_pt)[index], system, "es");
        c1[1]->stack((*pjet_wta_f_dphi_d_pt)[index], system, "wta");
    });

    nevt_d_hf->apply([&](TH1*, int64_t index) {
        c1[2]->add((*pjet_es_f_dphi_d_hf)[index], system, "es");
        c1[2]->stack((*pjet_wta_f_dphi_d_hf)[index], system, "wta");
    });

    zip([&](paper*& c, int64_t rows, std::string const& suffix,
            std::function<void(int64_t)> text) {
        c = new paper(tag + "_x_" + suffix, hb);
        c->divide(-1, rows);
        c->accessory(text);

        apply_style(c, collisions, -0.1, 1.5);
        c->accessory(std::bind(line_at, _1, 0.f, rx[0], rx[1]));
    }, c2, x{ ihf->size(), 1L, 1L }, suffixes, texts);

    pjet_f_x->apply([&](TH1* h) { c2[0]->add(h, system); });
    pjet_f_x_d_pt->apply([&](TH1* h) { c2[1]->add(h, system); });
    pjet_f_x_d_hf->apply([&](TH1* h) { c2[2]->add(h, system); });

    zip([&](paper*& c, int64_t rows, std::string const& suffix,
            std::function<void(int64_t)> text) {
        c = new paper(tag + "_ddr_" + suffix, hb);
        c->divide(-1, rows);
        c->accessory(text);

        apply_style(c, collisions, -2., 27.);
        c->accessory(std::bind(line_at, _1, 0.f, rdr[0], rdr[1]));
    }, c3, x{ ihf->size(), 1L, 1L }, suffixes, texts);

    pjet_f_ddr->apply([&](TH1* h) { c3[0]->add(h, system); });
    pjet_f_ddr_d_pt->apply([&](TH1* h) { c3[1]->add(h, system); });
    pjet_f_ddr_d_hf->apply([&](TH1* h) { c3[2]->add(h, system); });

    zip([&](paper*& c, int64_t rows, std::string const& suffix,
            std::function<void(int64_t)> text) {
        c = new paper(tag + "_jpt_" + suffix, hb);
        c->divide(-1, rows);
        c->accessory(text);

        apply_style(c, collisions, -0.001, 0.02);
        c->accessory(std::bind(line_at, _1, 0.f, rjpt[0], rjpt[1]));
    }, c4, x{ ihf->size(), 1L, 1L }, suffixes, texts);

    pjet_f_jpt->apply([&](TH1* h) { c4[0]->add(h, system); });
    pjet_f_jpt_d_pt->apply([&](TH1* h) { c4[1]->add(h, system); });
    pjet_f_jpt_d_hf->apply([&](TH1* h) { c4[2]->add(h, system); });

    hb->set_binary("system");
    hb->sketch();

    for (auto const& c : { c1, c2, c3, c4 })
        for (auto p : c) { p->draw("pdf"); }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return accumulate(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
