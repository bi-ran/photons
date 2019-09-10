#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"
#include "../git/history/include/memory.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/trunk.h"
#include "../git/tricks-and-treats/include/zip.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

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
TH1F* fold(T* flat, multival const* m, int64_t axis,
           std::array<int64_t, 2> const& offset) {
    auto name = std::string(flat->GetName()) + "_fold" + std::to_string(axis);
    auto hfold = m->axis(axis).book<TH1F, 2>(0, name, "", offset);

    for (int64_t i = 0; i < m->size(); ++i) {
        auto index = m->indices_for(i)[axis] - offset[0];

        if (index < 0 || index >= hfold->GetNbinsX()) { continue; }

        hfold->SetBinContent(index + 1, hfold->GetBinContent(index + 1)
            + flat->GetBinContent(i + 1));
        hfold->SetBinError(index + 1, std::sqrt(hfold->GetBinError(index + 1)
            * hfold->GetBinError(index + 1) + flat->GetBinError(i + 1)
            * flat->GetBinError(i + 1)));
    }

    return hfold;
}

template <typename T>
TH2F* shade(T* flat, multival const* m, std::array<int64_t, 4> const& offset) {
    auto name = std::string(flat->GetName()) + "_shade";
    auto hshade = m->book<TH2F, 4>(0, name, "", offset);

    for (int64_t i = 0; i < m->size(); ++i) {
        auto indices = m->indices_for(i);
        indices[0] = indices[0] - offset[0];
        indices[1] = indices[1] - offset[2];

        if (indices[0] < 0 || indices[0] >= hshade->GetNbinsX()
                || indices[1] < 0 || indices[1] >= hshade->GetNbinsY()) {
            continue; }

        hshade->SetBinContent(indices[0] + 1, indices[1] + 1,
            flat->GetBinContent(i + 1));
        hshade->SetBinError(indices[0] + 1, indices[1] + 1,
            flat->GetBinError(i + 1));
    }

    return hshade;
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

    auto ihf = new interval(dhf);
    auto mhf = new multival(dhf);

    auto mr = new multival(rdrr, rptr);
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

    auto sresult = new memory<TH2F>("sresult"s, "", null<TH2F>, mhf);
    auto srefold = new memory<TH2F>("srefold"s, "", null<TH2F>, mhf);
    auto svictim = new memory<TH2F>("svictim"s, "", null<TH2F>, mhf);

    constexpr int points = 30;

    std::vector<int32_t> tau(mhf->size(), -1);

    /* info text */
    std::function<void(int64_t, float)> hf_info = [&](int64_t x, float pos) {
        info_text(x, pos, "%i - %i%%", dcent, true); };

    auto system_info = system + " #sqrt{s_{NN}} = 5.02 TeV";

    /* figures */
    auto hb = new pencil();
    hb->category("type", "data", "unfolded", "refolded", "optimal");

    std::vector<paper*> cs(10, nullptr);
    zip([&](paper*& c, std::string const& title) {
        c = new paper(tag + "_dhf_" + title, hb);
        apply_style(c, system_info);
        c->accessory(std::bind(hf_info, _1, 0.75));
        c->divide(ihf->size(), -1);
    }, cs, (std::initializer_list<std::string> const) {
        "matrices"s, "unfold"s, "refold"s, "logtaux"s, "lcurve"s,
        "fold0"s, "fold1"s, "sresult"s, "srefold"s, "svictim"s });

    matrices->apply([&](TH2F* h) {
        cs[0]->add(h);
        cs[0]->adjust(h, "colz", "");
    });

    /* unfold */
    uf->apply([&](TUnfoldDensity* u, int64_t i) {
        if (u->SetInput((*victims)[i]) > 9999) {
            printf("  [!] error: set input\n"); exit(1);
        }

        tau[i] = u->ScanLcurve(points, 0., 0., &(*lcurve)[i],
                               &(*logtaux)[i], &(*logtauy)[i]);

        (*result)[i] = u->GetOutput("Unfolded");
        (*refold)[i] = u->GetFoldedOutput("FoldedBack");

        (*fold0)[i] = fold((*result)[i], mg, 0, { 0, 0 });
        (*fold1)[i] = fold((*result)[i], mg, 1, { 1, 1 });

        (*sresult)[i] = shade((*result)[i], mg, { 0, 0, 1, 1 });
        (*srefold)[i] = shade((*refold)[i], mr, { 0, 0, 0, 0 });
        (*svictim)[i] = shade((*victims)[i], mr, { 0, 0, 0, 0 });

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

        cs[1]->add((*result)[i], "unfolded");

        cs[2]->add((*victims)[i], "data");
        cs[2]->stack((*refold)[i], "refolded");

        cs[3]->add(hframe);
        cs[3]->stack((*logtaux)[i]);
        cs[3]->adjust((*logtaux)[i], "l", "");
        cs[3]->stack(logtau_opt, "optimal");
        cs[3]->adjust(logtau_opt, "p", "");

        cs[4]->add((*lcurve)[i]);
        cs[4]->adjust((*lcurve)[i], "al", "");
        cs[4]->stack(lcurve_opt, "optimal");
        cs[4]->adjust(lcurve_opt, "p", "");

        cs[5]->add((*fold0)[i]);
        cs[6]->add((*fold1)[i]);

        cs[7]->add((*sresult)[i]);
        cs[7]->adjust((*sresult)[i], "colz", "");
        cs[8]->add((*srefold)[i]);
        cs[8]->adjust((*srefold)[i], "colz", "");
        cs[9]->add((*svictim)[i]);
        cs[9]->adjust((*svictim)[i], "colz", "");
    });

    hb->sketch();
    for (auto c : cs) { c->draw("pdf"); }

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

        sresult->save(tag);
        srefold->save(tag);
        svictim->save(tag);
    });

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return undulate(argv[1], argv[2]);

    return 0;
}
