#ifndef LAMBDAS_H
#define LAMBDAS_H

#include <array>
#include <string>

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"

class paper;

auto histogram_formatter = [](TH1* obj, double min, double max) {
    obj->SetStats(0);
    obj->SetMarkerSize(0.84);
    obj->SetAxisRange(min, max, "Y");
};

auto default_decorator = [](std::string const& system) {
    TLatex* cms = new TLatex();
    cms->SetTextFont(62);
    cms->SetTextSize(0.048);
    cms->SetTextAlign(13);
    cms->DrawLatexNDC(0.135, 0.87, "CMS");

    TLatex* prelim = new TLatex();
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.032);
    prelim->SetTextAlign(13);
    prelim->DrawLatexNDC(0.135, 0.83, "Preliminary");

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

void apply_default_style(paper* p, std::string const& system) {
    using namespace std::placeholders;

    p->decorate(std::bind(default_decorator, system));
    p->legend(std::bind(coordinates, 0.45, 0.9, 0.87, 0.04));
    p->style(std::bind(default_legend_style, _1, 43, 12));
}

#endif /* LAMBDAS_H */
