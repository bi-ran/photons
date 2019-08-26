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

#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLine.h"

#include <memory>
#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

static auto const red = TColor::GetColor("#f2777a");
static auto const blue = TColor::GetColor("#6699cc");

int congratulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto inputs = conf->get<std::vector<std::string>>("inputs");
    auto figures = conf->get<std::vector<std::string>>("figures");
    auto xmins = conf->get<std::vector<float>>("xmin");
    auto xmaxs = conf->get<std::vector<float>>("xmax");
    auto ymins = conf->get<std::vector<float>>("ymin");
    auto ymaxs = conf->get<std::vector<float>>("ymax");
    auto oflows = conf->get<std::vector<bool>>("oflow");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");

    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    auto ipt = std::make_shared<interval>(dpt);
    auto ihf = std::make_shared<interval>(dhf);

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* open input files */
    std::vector<TFile*> files(inputs.size(), nullptr);
    zip([&](auto& file, auto const& input) {
        file = new TFile(input.data(), "read");
    }, files, inputs);

    /* load histograms */
    std::vector<std::string> tags = {
        "aa"s,
        "pp"s,
        "pp_smear_50_90"s,
        "pp_smear_30_50"s,
        "pp_smear_10_30"s,
        "pp_smear_0_10"s
    };

    std::vector<std::string> base_stubs(6);
    std::vector<std::string> syst_stubs(6);

    zip([&](auto& base, auto& syst, auto const& tag) {
        base = tag + "_base_" + tag + "_nominal_s_pure_raw_sub_";
        syst = tag + "_total_base_" + tag + "_nominal_s_pure_raw_sub_";
    }, base_stubs, syst_stubs, tags);

    /* prepare plots */
    auto hb = new pencil();
    hb->category("system", "pp", "aa", "ss");

    hb->alias("aa", "PbPb");
    hb->alias("ss", "pp (smeared)");

    auto decorator = [](std::string const& system) {
        TLatex* cms = new TLatex();
        cms->SetTextFont(62);
        cms->SetTextSize(0.084);
        cms->SetTextAlign(13);
        cms->DrawLatexNDC(0.135, 0.87, "CMS");

        TLatex* prel = new TLatex();
        prel->SetTextFont(52);
        prel->SetTextSize(0.056);
        prel->SetTextAlign(13);
        prel->DrawLatexNDC(0.135, 0.81, "Preliminary");

        TLatex* com = new TLatex();
        com->SetTextFont(42);
        com->SetTextSize(0.056);
        com->SetTextAlign(11);
        com->DrawLatexNDC(0.11, 0.92, "#sqrt{s_{NN}} = 5.02 TeV");

        TLatex* info = new TLatex();
        info->SetTextFont(42);
        info->SetTextSize(0.056);
        info->SetTextAlign(31);
        info->DrawLatexNDC(0.89, 0.92, system.data());
    };

    auto pt_info = [&](int64_t index, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{#gamma} < %.0f",
            (*ipt)[index - 1], (*ipt)[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(20);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto hf_info = [&](int64_t index, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%i - %i%%", dcent[index], dcent[index - 1]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(20);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto pp_text = [&](int64_t index) {
        pt_info(index, 0.67);
    };

    auto aa_text = [&](int64_t index) {
        auto pt_x = (index - 1) % (ipt->size() - 1);
        auto hf_x = (index - 1) / (ipt->size() - 1);

        pt_info(pt_x + 1, 0.67);
        hf_info(hf_x + 1, 0.60);
    };

    zip([&](auto const& figure, auto xmin, auto xmax, auto ymin, auto ymax,
            auto integral) {
        /* get histograms */
        std::vector<history*> hists(6, nullptr);
        std::vector<history*> systs(6, nullptr);

        zip([&](auto& hist, auto& syst, auto const file,
                auto const& base_stub, auto const& syst_stub) {
            hist = new history(file, base_stub + figure);
            syst = new history(file, syst_stub + figure);
        }, hists, systs, files, base_stubs, syst_stubs);

        /* link histograms, uncertainties */
        std::unordered_map<TH1*, TH1*> links;
        zip([&](auto hist, auto syst) {
            hist->apply([&](TH1* h, int64_t index) {
                links[h] = (*syst)[index]; });
        }, hists, systs);

        std::unordered_map<TH1*, int32_t> colours;
        hists[0]->apply([&](TH1* h) { colours[h] = 1; });

        /* uncertainty box */
        auto box = [&](TH1* h, int64_t) {
            TGraph* gr = new TGraph();
            gr->SetFillStyle(1001);
            gr->SetFillColorAlpha(colours[h] ? red : blue, 0.48);

            for (int i = 1; i <= h->GetNbinsX(); ++i) {
                if (h->GetBinError(i) == 0) continue;

                double x = h->GetBinCenter(i);
                double width = h->GetBinWidth(i);
                double val = h->GetBinContent(i);
                double err = links[h]->GetBinContent(i);

                gr->SetPoint(0, x - (width / 2), val - err);
                gr->SetPoint(1, x + (width / 2), val - err);
                gr->SetPoint(2, x + (width / 2), val + err);
                gr->SetPoint(3, x - (width / 2), val + err);

                gr->DrawClone("f");
            }
        };

        /* minor adjustments */
        if (integral) { xmin = convert_pi(xmin); xmax = convert_pi(xmax); }

        /* prepare papers */
        auto p = new paper("results_pp_" + figure, hb);
        apply_style(p, "", ymin, ymax);
        p->decorate(std::bind(decorator, "pp 320 pb^{-1}"));
        p->accessory(std::bind(line_at, _1, 0.f, xmin, xmax));
        p->accessory(pp_text);
        p->jewellery(box);
        p->divide(-1, 1);

        auto a = new paper("results_aa_" + figure, hb);
        apply_style(a, "", ymin, ymax);
        a->decorate(std::bind(decorator, "PbPb 1.6 nb^{-1}"));
        a->accessory(std::bind(line_at, _1, 0.f, xmin, xmax));
        a->accessory(aa_text);
        a->jewellery(box);
        a->divide(-1, ihf->size());

        auto s = new paper("results_ss_" + figure, hb);
        apply_style(s, "", ymin, ymax);
        s->decorate(std::bind(decorator, "pp 320 pb^{-1}, PbPb 1.6 nb^{-1}"));
        s->accessory(std::bind(line_at, _1, 0.f, xmin, xmax));
        s->accessory(aa_text);
        s->jewellery(box);
        s->divide(-1, ihf->size());

        /* draw histograms with uncertainties */
        hists[0]->apply([&](TH1* h) { a->add(h, "aa"); s->add(h, "aa"); });
        hists[1]->apply([&](TH1* h) { p->add(h, "pp"); });

        for (int64_t i = 0; i < 4; ++i) {
            hists[i + 2]->apply([&](TH1* h, int64_t index) {
                s->stack((ipt->size() - 1) * i + index + 1, h, "ss");
            });
        }

        auto pp_style = [](TH1* h) {
            h->SetLineColor(1);
            h->SetMarkerStyle(25);
            h->SetMarkerSize(0.60);
        };

        auto aa_style = [](TH1* h) {
            h->SetLineColor(1);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.60);
        };

        hb->style("pp", pp_style);
        hb->style("aa", aa_style);
        hb->style("ss", pp_style);
        hb->sketch();

        p->draw("pdf");
        a->draw("pdf");
        s->draw("pdf");
    }, figures, xmins, xmaxs, ymins, ymaxs, oflows);

    in(output, []() {});

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return congratulate(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
