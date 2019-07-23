#include "../include/lambdas.h"

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
        obj->Scale(1. / obj->Integral(), "width"); }), 0)... };
}

int accumulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto label = conf->get<std::string>("label");
    auto tag = conf->get<std::string>("tag");
    auto system = conf->get<std::string>("system");

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
    auto nevt = new history(f, label + "_raw_nevt"s);

    auto pjet_es_f_dphi = new history(f, label + "_raw_sub_pjet_es_f_dphi"s);
    auto pjet_wta_f_dphi = new history(f, label + "_raw_sub_pjet_wta_f_dphi"s);
    auto pjet_f_x = new history(f, label + "_raw_sub_pjet_f_x"s);
    auto pjet_f_ddr = new history(f, label + "_raw_sub_pjet_f_ddr"s);

    /* rescale by number of signal photons (events) */
    pjet_es_f_dphi->multiply(*nevt);
    pjet_wta_f_dphi->multiply(*nevt);
    pjet_f_x->multiply(*nevt);
    pjet_f_ddr->multiply(*nevt);

    /* integrate histograms */
    auto nevt_d_pt = nevt->sum(1);
    auto nevt_d_hf = nevt->sum(0);

    auto pjet_es_f_dphi_d_pt = pjet_es_f_dphi->sum(1);
    auto pjet_wta_f_dphi_d_pt = pjet_wta_f_dphi->sum(1);
    auto pjet_f_x_d_pt = pjet_f_x->sum(1);
    auto pjet_f_x_d_hf = pjet_f_x->sum(0);
    auto pjet_f_ddr_d_pt = pjet_f_ddr->sum(1, 1);
    auto pjet_f_ddr_d_hf = pjet_f_ddr->sum(0, 1);

    /* normalise by number of signal photons (events) */
    pjet_es_f_dphi_d_pt->divide(*nevt_d_pt);
    pjet_wta_f_dphi_d_pt->divide(*nevt_d_pt);
    pjet_f_x_d_pt->divide(*nevt_d_pt);
    pjet_f_x_d_hf->divide(*nevt_d_hf);
    pjet_f_ddr_d_pt->divide(*nevt_d_pt);
    pjet_f_ddr_d_hf->divide(*nevt_d_hf);

    /* scale by bin width */
    scale_bin_width(
        pjet_f_ddr_d_pt,
        pjet_f_ddr_d_hf,
        pjet_f_x_d_pt,
        pjet_f_x_d_hf);

    /* normalise to unity */

    /* save histograms */
    pjet_es_f_dphi_d_pt->save(tag);
    pjet_wta_f_dphi_d_pt->save(tag);
    pjet_f_x_d_pt->save(tag);
    pjet_f_x_d_hf->save(tag);
    pjet_f_ddr_d_pt->save(tag);
    pjet_f_ddr_d_hf->save(tag);

    /* draw plots */
    printf("painting..\n");

    auto photon_pt_selection = [&](int64_t index) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{#gamma} < %.0f",
            (*ipt)[index - 1], (*ipt)[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, buffer);
    };

    auto hf_selection = [&](int64_t index) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%i - %i%%", dcent[index], dcent[index - 1]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, 0.75, buffer);
    };

    auto line_at_unity = [&](int64_t, float low, float high) {
        TLine* l1 = new TLine(low, 0., high, 0.);
        l1->SetLineStyle(7);
        l1->Draw();
    };

    auto hb = new pencil();
    hb->category("system", "PbPb", "pp");
    hb->category("axis", "es", "wta", "na");

    hb->ditto("na", "es");

    hb->alias("es", "E-Scheme");
    hb->alias("wta", "WTA");
    hb->alias("na", "");

    auto system_tag = system + " #sqrt{s_{NN}} = 5.02 TeV"s;

    auto c1 = new paper("dphi_d_pt", hb);
    apply_default_style(c1, system_tag, -0.08, 0.6);
    c1->accessory(photon_pt_selection);
    c1->accessory(std::bind(line_at_unity, _1, rdphi[0], rdphi[1]));

    for (int64_t i = 0; i < ipt->size() - 1; ++i) {
        c1->add((*pjet_es_f_dphi_d_pt)[i], system, "es");
        c1->stack((*pjet_wta_f_dphi_d_pt)[i], system, "wta");
    }

    auto c2 = new paper("ddr_d_pt", hb);
    apply_default_style(c2, system_tag, -1., 12.);
    c2->accessory(photon_pt_selection);
    c2->accessory(std::bind(line_at_unity, _1, rdr[0], rdr[1]));

    for (int64_t i = 0; i < ipt->size() - 1; ++i)
        c2->add((*pjet_f_ddr_d_pt)[i], system, "na");

    auto c3 = new paper("ddr_d_hf", hb);
    apply_default_style(c3, system_tag, -1., 12.);
    c3->accessory(hf_selection);
    c3->accessory(std::bind(line_at_unity, _1, rdr[0], rdr[1]));

    for (int64_t i = 0; i < ihf->size(); ++i)
        c3->add((*pjet_f_ddr_d_hf)[i], system, "na");

    auto c4 = new paper("x_d_pt", hb);
    apply_default_style(c4, system_tag, -0.1, 1.2);
    c4->accessory(photon_pt_selection);
    c4->accessory(std::bind(line_at_unity, _1, rx[0], rx[1]));

    for (int64_t i = 0; i < ipt->size() - 1; ++i)
        c4->add((*pjet_f_x_d_pt)[i], system, "na");

    auto c5 = new paper("x_d_hf", hb);
    apply_default_style(c5, system_tag, -0.1, 1.2);
    c5->accessory(hf_selection);
    c5->accessory(std::bind(line_at_unity, _1, rx[0], rx[1]));

    for (int64_t i = 0; i < ihf->size(); ++i)
        c5->add((*pjet_f_x_d_hf)[i], system, "na");

    hb->set_binary("type");
    hb->sketch();

    c1->draw("pdf");
    c2->draw("pdf");
    c3->draw("pdf");
    c4->draw("pdf");
    c5->draw("pdf");

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return accumulate(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
