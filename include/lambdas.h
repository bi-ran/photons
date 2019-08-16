#ifndef LAMBDAS_H
#define LAMBDAS_H

#include <array>
#include <functional>
#include <string>

#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

auto _for_content = [](TH1* h, float (*f)(float)) {
    for (int64_t j = 1; j <= h->GetNbinsX(); ++j) {
        auto val = h->GetBinContent(j);
        h->SetBinContent(j, f(val));
    }
};

auto _for_content_index = [](TH1* h, float (*f)(float, int64_t)) {
    for (int64_t j = 1; j <= h->GetNbinsX(); ++j) {
        auto val = h->GetBinContent(j);
        h->SetBinContent(j, f(val, j));
    }
};

auto graph_formatter = [](TGraph* obj) {
    obj->SetMarkerSize(0.84);
    obj->GetXaxis()->CenterTitle();
    obj->GetYaxis()->CenterTitle();
};

auto hist_formatter = [](TH1* obj) {
    obj->SetStats(0);
    obj->SetMarkerSize(0.84);
    obj->GetXaxis()->CenterTitle();
    obj->GetYaxis()->CenterTitle();
};

auto default_formatter = [](TH1* obj, double min, double max) {
    hist_formatter(obj);
    obj->SetAxisRange(min, max, "Y");
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

auto line_at = [&](int64_t, float val, float low, float high) {
    TLine* l1 = new TLine(low, val, high, val);
    l1->SetLineStyle(7);
    l1->Draw();
};

template <typename T, typename U>
void apply_style(T p, std::string const& text, U formatter) {
    using namespace std::placeholders;

    p->format(formatter);
    p->format(graph_formatter);
    p->decorate(std::bind(default_decorator, text));
    p->legend(std::bind(coordinates, 0.45, 0.9, 0.87, 0.04));
    p->style(std::bind(default_legend_style, _1, 43, 12));
}

template <typename T>
void apply_style(T p, std::string const& text, double min, double max) {
    using namespace std::placeholders;

    apply_style(p, text, std::bind(default_formatter, _1, min, max));
}

template <typename T>
void apply_style(T p, std::string const& text) {
    apply_style(p, text, hist_formatter);
}

#endif /* LAMBDAS_H */
