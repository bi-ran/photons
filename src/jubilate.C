#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/overflow_angles.h"
#include "../git/tricks-and-treats/include/trunk.h"

#include "TFile.h"
#include "TH1.h"

using namespace std::literals::string_literals;
using namespace std::placeholders;

template <typename... T>
void scale_bin_width(T&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        obj->Scale(1., "width"); }), 0)... };
}

template <typename... T>
void normalise_to_unity(T&... args) {
    (void)(int [sizeof...(T)]) { (args->apply([](TH1* obj) {
        obj->Scale(1. / obj->Integral("width")); }), 0)... };
}

int jubilate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto rdphi = conf->get<std::vector<float>>("dphi_range");
    auto rx = conf->get<std::vector<float>>("x_range");
    auto rdr = conf->get<std::vector<float>>("dr_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");

    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    /* convert to integral angle units (cast to double) */
    convert_in_place_pi(rdphi);

    auto ipt = std::make_shared<interval>(dpt);
    auto ihf = std::make_shared<interval>(dhf);

    /* load history objects */
    TFile* f = new TFile(input.data(), "read");

    TH1::SetDefaultSumw2();

    auto nevt = std::make_unique<history>(f, "raw_nevt");

    auto pjet_es_f_dphi = std::make_unique<history>(f, "raw_pjet_es_f_dphi");
    auto pjet_wta_f_dphi = std::make_unique<history>(f, "raw_pjet_wta_f_dphi");
    auto pjet_f_x = std::make_unique<history>(f, "raw_pjet_f_x");
    auto pjet_f_ddr = std::make_unique<history>(f, "raw_pjet_f_ddr");

    auto mix_pjet_es_f_dphi = std::make_unique<history>(
        f, "raw_mix_pjet_es_f_dphi");
    auto mix_pjet_wta_f_dphi = std::make_unique<history>(
        f, "raw_mix_pjet_wta_f_dphi");
    auto mix_pjet_f_x = std::make_unique<history>(
        f, "raw_mix_pjet_f_x");
    auto mix_pjet_f_ddr = std::make_unique<history>(
        f, "raw_mix_pjet_f_ddr");

    /* shrink to remove overflow photon pt bin */
    auto shape = nevt->shape();
    shape[0] = shape[0] - 1;

    auto wrap = [&](std::unique_ptr<history>& h) {
        h = h->shrink("s", shape, { 0, 0 }); };

    wrap(nevt);
    wrap(pjet_es_f_dphi);
    wrap(pjet_wta_f_dphi);
    wrap(pjet_f_x);
    wrap(pjet_f_ddr);
    wrap(mix_pjet_es_f_dphi);
    wrap(mix_pjet_wta_f_dphi);
    wrap(mix_pjet_f_x);
    wrap(mix_pjet_f_ddr);

    /* scale by bin width */
    scale_bin_width(
        pjet_f_x,
        pjet_f_ddr,
        mix_pjet_f_x,
        mix_pjet_f_ddr);

    /* draw figures */
    auto redraw_dphi_axis = [&](TH1* h, int64_t) {
        transform_axis(h, [](int64_t val) -> float {
            return std::abs(revert_radian(val)); }); };

    auto info_text = [&](int64_t index) {
        auto indices = nevt->indices_for(index - 1);
        auto pt_x = indices[0];
        auto hf_x = indices[1];

        char buffer[128] = { '\0' };

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);

        sprintf(buffer, "%.0f < p_{T}^{#gamma} < %.0f",
            (*ipt)[pt_x], (*ipt)[pt_x + 1]);
        l->DrawLatexNDC(0.135, 0.75, buffer);

        sprintf(buffer, "%i - %i%%", dcent[hf_x + 1], dcent[hf_x]);
        l->DrawLatexNDC(0.135, 0.71, buffer);
    };

    auto hb = new pencil();
    hb->category("system", "pp", "PbPb");
    hb->category("type", "raw", "mix");

    hb->set_binary("type");

    auto collisions = system + " #sqrt{s_{NN}} = 5.02 TeV"s;

    auto c1 = new paper(tag + "_mixing_dphi_es_d_pthf", hb);
    apply_style(c1, collisions, -0.04, 0.4);
    c1->accessory(std::bind(line_at, _1, 0.f, rdphi[0], rdphi[1]));
    c1->accessory(info_text);
    c1->jewellery(redraw_dphi_axis);
    c1->divide(-1 , ihf->size());

    auto c2 = new paper(tag + "_mixing_dphi_wta_d_pthf", hb);
    apply_style(c2, collisions, -0.04, 0.4);
    c2->accessory(std::bind(line_at, _1, 0.f, rdphi[0], rdphi[1]));
    c2->accessory(info_text);
    c2->jewellery(redraw_dphi_axis);
    c2->divide(-1 , ihf->size());

    auto c3 = new paper(tag + "_mixing_x_d_pthf", hb);
    apply_style(c3, collisions, -0.1, 2.0);
    c3->accessory(std::bind(line_at, _1, 0.f, rx[0], rx[1]));
    c3->accessory(info_text);
    c3->divide(-1 , ihf->size());

    auto c4 = new paper(tag + "_mixing_ddr_d_pthf", hb);
    apply_style(c4, collisions, -1., 24.);
    c4->accessory(std::bind(line_at, _1, 0.f, rdr[0], rdr[1]));
    c4->accessory(info_text);
    c4->divide(-1 , ihf->size());

    for (int64_t i = 0; i < nevt->size(); ++i) {
        c1->add((*pjet_es_f_dphi)[i], system, "raw");
        c1->stack((*mix_pjet_es_f_dphi)[i], system, "mix");

        c2->add((*pjet_wta_f_dphi)[i], system, "raw");
        c2->stack((*mix_pjet_wta_f_dphi)[i], system, "mix");

        c3->add((*pjet_f_x)[i], system, "raw");
        c3->stack((*mix_pjet_f_x)[i], system, "mix");

        c4->add((*pjet_f_ddr)[i], system, "raw");
        c4->stack((*mix_pjet_f_ddr)[i], system, "mix");
    }

    hb->sketch();

    for (auto const& c : { c1, c2, c3, c4 })
        c->draw("pdf");

    /* save output */
    in(output, []() {});

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return jubilate(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
