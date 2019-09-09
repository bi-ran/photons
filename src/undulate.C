#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"
#include "../git/history/include/memory.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/trunk.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"

#include "TUnfoldDensity.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

template <typename T>
T* null(int64_t, std::string const&, std::string const&) {
    return nullptr;
}

template <typename T>
TH1F* fold(T* flat, multival const* m, int64_t axis) {
    auto name = std::string(flat->GetName()) + "_fold" + std::to_string(axis);
    auto hfold = m->axis(axis).book<TH1F>(0, name, "");

    for (int64_t i = 0; i < m->size(); ++i) {
        auto index = m->indices_for(i)[axis];
        hfold->SetBinContent(index + 1, hfold->GetBinContent(index + 1)
            + flat->GetBinContent(i + 1));
        hfold->SetBinError(index + 1, std::sqrt(hfold->GetBinError(index + 1)
            * hfold->GetBinError(index + 1) + flat->GetBinError(i + 1)
            * flat->GetBinError(i + 1)));
    }

    return hfold;
}

int undulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto victim = conf->get<std::string>("victim");
    auto label = conf->get<std::string>("label");

    auto rdrr = conf->get<std::vector<float>>("drr_range");
    auto rdrg = conf->get<std::vector<float>>("drg_range");
    auto rptr = conf->get<std::vector<float>>("ptr_range");
    auto rptg = conf->get<std::vector<float>>("ptg_range");

    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    auto mhf = new multival(dhf);

    auto mg = new multival(rdrg, rptg);

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load input (and victims) */
    TFile* fi = new TFile(input.data(), "read");
    auto matrices = new history<TH2F>(fi, tag + "_c");

    TFile* fv = new TFile(victim.data(), "read");
    auto victims = new history<TH1F>(fv, label);

    /* prepare objects for unfolding */
    auto factory = [&](int64_t i, std::string const&, std::string const&) {
        return new TUnfoldDensity((*matrices)[i], TUnfold::kHistMapOutputVert);
    };

    auto uf = new memory<TUnfoldDensity>("uf"s, "", factory, mhf);

    auto logtaux = new memory<TSpline>("logtaux"s, "", null<TSpline>, mhf);
    auto logtauy = new memory<TSpline>("logtauy"s, "", null<TSpline>, mhf);
    auto lcurve = new memory<TGraph>("lcurve"s, "", null<TGraph>, mhf);

    auto result = new memory<TH1>("result"s, "", null<TH1>, mhf);
    auto refold = new memory<TH1>("refold"s, "", null<TH1>, mhf);

    auto fold0 = new memory<TH1F>("fold0"s, "", null<TH1F>, mhf);
    auto fold1 = new memory<TH1F>("fold1"s, "", null<TH1F>, mhf);

    constexpr int points = 30;

    std::vector<int32_t> tau(mhf->size(), -1);

    /* info text */
    std::function<void(int64_t, float)> hf_info = [&](int64_t x, float pos) {
        info_text(x, pos, "%i - %i%%", dcent, true); };

    auto system_info = system + " #sqrt{s_{NN}} = 5.02 TeV";

    /* figures */
    auto hb = new pencil();
    hb->category("type", "data", "unfolded", "refolded", "optimal");

    auto c1 = new paper(tag + "_dpthf_matrices", hb);
    apply_style(c1, system_info);
    c1->accessory(std::bind(hf_info, _1, 0.75));
    c1->divide(1, -1);

    matrices->apply([&](TH2F* h) { c1->add(h); c1->adjust(h, "colz", ""); });

    auto c2 = new paper(tag + "_dpthf_unfold", hb);
    apply_style(c2, system_info, -0.005, 0.05);
    c2->accessory(std::bind(hf_info, _1, 0.75));
    c2->divide(1, -1);

    auto c3 = new paper(tag + "_dpthf_refold", hb);
    apply_style(c3, system_info, -0.005, 0.05);
    c3->accessory(std::bind(hf_info, _1, 0.75));
    c3->divide(1, -1);

    auto c4 = new paper(tag + "_dpthf_logtaux", hb);
    apply_style(c4, system_info);
    c4->accessory(std::bind(hf_info, _1, 0.75));
    c4->divide(1, -1);

    auto c5 = new paper(tag + "_dpthf_lcurve", hb);
    apply_style(c5, system_info);
    c5->accessory(std::bind(hf_info, _1, 0.75));
    c5->divide(1, -1);

    auto c6 = new paper(tag + "_dpthf_fold0", hb);
    apply_style(c6, system_info);
    c6->accessory(std::bind(hf_info, _1, 0.75));
    c6->divide(1, -1);

    auto c7 = new paper(tag + "_dpthf_fold1", hb);
    apply_style(c7, system_info);
    c7->accessory(std::bind(hf_info, _1, 0.75));
    c7->divide(1, -1);

    /* unfold */
    uf->apply([&](TUnfoldDensity* u, int64_t i) {
        if (u->SetInput((*victims)[i]) > 9999) {
            printf("  [!] error: set input\n"); exit(1);
        }

        tau[i] = u->ScanLcurve(points, 0., 0., &(*lcurve)[i],
                               &(*logtaux)[i], &(*logtauy)[i]);

        (*result)[i] = u->GetOutput("Unfolded");
        (*refold)[i] = u->GetFoldedOutput("FoldedBack");

        (*fold0)[i] = fold((*result)[i], mg, 0);
        (*fold1)[i] = fold((*result)[i], mg, 1);

        double t;
        double x;
        double y;

        (*logtaux)[i]->GetKnot(tau[i], t, x);
        (*logtauy)[i]->GetKnot(tau[i], t, y);

        auto logtau_opt = new TGraph(1, &t, &x);
        logtau_opt->SetMarkerStyle(21);
        auto lcurve_opt = new TGraph(1, &x, &y);
        lcurve_opt->SetMarkerStyle(21);

        auto hframe = frame((*logtaux)[i], (*lcurve)[i]->GetXaxis());

        c2->add((*result)[i], "unfolded");

        c3->add((*victims)[i], "data");
        c3->stack((*refold)[i], "refolded");

        c4->add(hframe);
        c4->stack((*logtaux)[i]);
        c4->adjust((*logtaux)[i], "l", "");
        c4->stack(logtau_opt, "optimal");
        c4->adjust(logtau_opt, "p", "");

        c5->add((*lcurve)[i]);
        c5->adjust((*lcurve)[i], "al", "");
        c5->stack(lcurve_opt, "optimal");
        c5->adjust(lcurve_opt, "p", "");

        c6->add((*fold0)[i]);

        c7->add((*fold1)[i]);
    });

    hb->sketch();
    for (auto c : { c1, c2, c3, c4, c5, c6, c7 })
        c->draw("pdf");

    /* save output */
    in(output, [&]() {
        matrices->save("");

        logtaux->save(tag);
        logtauy->save(tag);
        lcurve->save(tag);

        result->save(tag);
        refold->save(tag);
        fold0->save(tag);
        fold1->save(tag);
    });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return undulate(argv[1], argv[2]);

    return 0;
}
