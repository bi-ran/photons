#ifndef LAMBDAS_H
#define LAMBDAS_H

#include <array>
#include <string>

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

#include "../git/paper-and-pencil/include/paper.h"

auto _for_bins = [](TH1* h, float (*f)(float)) {
    for (int64_t j = 1; j <= h->GetNbinsX(); ++j) {
        auto val = h->GetBinContent(j);
        h->SetBinContent(j, f(val));
    }
};

auto graph_formatter = [](TGraph* obj) {
    obj->SetMarkerSize(0.84);
    obj->GetXaxis()->CenterTitle();
    obj->GetYaxis()->CenterTitle();
};

auto simple_formatter = [](TH1* obj) {
    obj->SetStats(0);
    obj->SetMarkerSize(0.84);
    obj->GetXaxis()->CenterTitle();
    obj->GetYaxis()->CenterTitle();
};

auto default_formatter = [](TH1* obj, double min, double max) {
    obj->SetStats(0);
    obj->SetMarkerSize(0.84);
    obj->SetAxisRange(min, max, "Y");
    obj->GetXaxis()->CenterTitle();
    obj->GetYaxis()->CenterTitle();
};

auto default_decorator = [](std::string const& system) {
    TLatex* cms = new TLatex();
    cms->SetTextFont(62);
    cms->SetTextSize(0.048);
    cms->SetTextAlign(13);
    cms->DrawLatexNDC(0.135, 0.87, "CMS");

    TLatex* lwip = new TLatex();
    lwip->SetTextFont(52);
    lwip->SetTextSize(0.032);
    lwip->SetTextAlign(13);
    lwip->DrawLatexNDC(0.135, 0.83, "Work-in-progress");

    TLatex* info = new TLatex();
    info->SetTextFont(42);
    info->SetTextSize(0.032);
    info->SetTextAlign(31);
    info->DrawLatexNDC(0.89, 0.92, system.data());
};

auto coordinates = [](float x0, float x1, float y1, float dy) {
    return std::array<float, 4>({ x0, x1, y1, dy });
};

auto default_legend_style = [](TLegend* l, int font, float size) {
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextSize(size);
};

void apply_default_style(paper* p, std::string const& system,
                         double min, double max) {
    using namespace std::placeholders;

    p->format(graph_formatter);
    p->format(std::bind(default_formatter, _1, min, max));
    p->decorate(std::bind(default_decorator, system));
    p->legend(std::bind(coordinates, 0.45, 0.9, 0.87, 0.04));
    p->style(std::bind(default_legend_style, _1, 43, 12));
}

auto line_at = [&](int64_t, float val, float low, float high) {
    TLine* l1 = new TLine(low, val, high, val);
    l1->SetLineStyle(7);
    l1->Draw();
};

#endif /* LAMBDAS_H */
