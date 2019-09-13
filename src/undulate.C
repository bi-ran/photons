#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

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

template <typename T, int64_t N>
TH1F* fold(T* flat, multival const* m, int64_t axis,
           std::array<int64_t, N> const& offsets) {
    auto name = std::string(flat->GetName()) + "_fold" + std::to_string(axis);
    auto hfold = m->axis(axis).book<TH1F, 2>(0, name, "",
        { offsets[axis << 1], offsets[(axis << 1) + 1] });

    auto shape = m->shape();

    for (int64_t i = 0; i < m->size(); ++i) {
        auto indices = m->indices_for(i);

        bool flag = false;
        for (int64_t j = 0; j < m->dims(); ++j)
            if (indices[j] >= shape[j] - offsets[(j << 1) + 1]
                    || indices[j] < offsets[j << 1])
                flag = true;
        if (flag) { continue; }

        auto index = indices[axis] - offsets[axis << 1];

        hfold->SetBinContent(index + 1, hfold->GetBinContent(index + 1)
            + flat->GetBinContent(i + 1));
        hfold->SetBinError(index + 1, std::sqrt(hfold->GetBinError(index + 1)
            * hfold->GetBinError(index + 1) + flat->GetBinError(i + 1)
            * flat->GetBinError(i + 1)));
    }

    hfold->Scale(1., "width");

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

    hshade->Scale(1., "width");

    return hshade;
}

int undulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto victim = conf->get<std::string>("victim");
    auto label = conf->get<std::string>("label");

    auto refs = conf->get<std::vector<std::string>>("refs");
    auto extensions = conf->get<std::vector<int64_t>>("extensions");
    auto dimensions = conf->get<std::vector<int64_t>>("dimensions");

    auto points = conf->get<int64_t>("points");

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

    /* load input, victims, and references */
    TFile* fi = new TFile(input.data(), "read");
    auto matrices = new history<TH2F>(fi, tag + "_c");

    TFile* fv = new TFile(victim.data(), "read");
    auto victims = new history<TH1F>(fv, label);

    std::vector<history<TH1F>*> notes(refs.size(), nullptr);
    zip([&](history<TH1F>*& n, std::string const& ref) {
        n = new history<TH1F>(fv, ref); }, notes, refs);

    zip([&](int64_t extension, int64_t dimension) {
        matrices = matrices->extend("ext"s + std::to_string(extension),
            extension, dimension); }, extensions, dimensions);

    auto shape = victims->shape();

    /* prepare objects for unfolding */
    auto factory = [&](int64_t i, std::string const&, std::string const&) {
        return new TUnfoldDensity((*matrices)[i],
            TUnfold::EHistMap::kHistMapOutputVert,
            TUnfold::ERegMode::kRegModeCurvature,
            TUnfold::EConstraint::kEConstraintNone,
            TUnfoldDensity::EDensityMode::kDensityModeNone
        );
    };

    auto uf = new history<TUnfoldDensity>("uf"s, "", factory, shape);

    auto bias = new history<TH1>("bias"s, "", null<TH1>, shape);
    auto ematrix = new history<TH2>("ematrix"s, "", null<TH2>, shape);
    auto lmatrix = new history<TH2>("lmatrix"s, "", null<TH2>, shape);
    auto logtaux = new history<TSpline>("logtaux"s, "", null<TSpline>, shape);
    auto logtauy = new history<TSpline>("logtauy"s, "", null<TSpline>, shape);
    auto lcurve = new history<TGraph>("lcurve"s, "", null<TGraph>, shape);

    auto shaded = new history<TH2F>(label + "_shade", "", null<TH2F>, shape);
    auto side0 = new history<TH1F>(label + "_side0", "", null<TH1F>, shape);
    auto side1 = new history<TH1F>(label + "_side1", "", null<TH1F>, shape);

    auto result = new history<TH1>("result"s, "", null<TH1>, shape);
    auto refold = new history<TH1>("refold"s, "", null<TH1>, shape);

    auto sresult = new history<TH2F>("sresult"s, "", null<TH2F>, shape);
    auto srefold = new history<TH2F>("srefold"s, "", null<TH2F>, shape);

    auto fold0 = new history<TH1F>("fold0"s, "", null<TH1F>, shape);
    auto fold1 = new history<TH1F>("fold1"s, "", null<TH1F>, shape);

    std::vector<int32_t> tau(mhf->size(), -1);

    /* info text */
    std::function<void(int64_t, float)> hf_info = [&](int64_t x, float pos) {
        info_text(x, pos, "%i - %i%%", dcent, true); };

    auto system_info = system + " #sqrt{s_{NN}} = 5.02 TeV";

    /* figures */
    auto hb = new pencil();
    hb->category("type", "data", "unfolded", "refolded", "optimal", "bias");
    hb->category("bins", "gen", "reco");

    hb->alias("gen", "");
    hb->alias("reco", "");

    std::vector<paper*> cs(13, nullptr);
    zip([&](paper*& c, std::string const& title) {
        c = new paper(tag + "_dhf_" + title, hb);
        apply_style(c, system_info);
        c->accessory(std::bind(hf_info, _1, 0.75));
        c->divide(-1, ihf->size());
    }, cs, (std::initializer_list<std::string> const) {
        "matrices"s, "bias"s, "ematrix"s, "lmatrix"s, "logtaux"s, "lcurve"s,
        "unfold"s, "refold"s, "sresult"s, "srefold"s, "fold0"s, "fold1"s,
        "shaded"s });

    cs[10]->format(std::bind(default_formatter, _1, -2., 27.));
    cs[11]->format(std::bind(default_formatter, _1, -0.001, 0.02));

    /* unfold */
    uf->apply([&](TUnfoldDensity* u, int64_t i) {
        if (u->SetInput((*victims)[i]) > 9999) {
            printf("  [!] error: set input\n"); exit(1);
        }

        tau[i] = u->ScanLcurve(points, 0., 0., &(*lcurve)[i],
                               &(*logtaux)[i], &(*logtauy)[i]);

        (*result)[i] = u->GetOutput("Unfolded");
        (*refold)[i] = u->GetFoldedOutput("FoldedBack");

        (*sresult)[i] = shade((*result)[i], mg, { 0, 0, 2, 1 });
        (*srefold)[i] = shade((*refold)[i], mr, { 0, 0, 0, 0 });

        (*fold0)[i] = fold<TH1, 4>((*result)[i], mg, 0, { 0, 0, 2, 1 });
        (*fold1)[i] = fold<TH1, 4>((*result)[i], mg, 1, { 0, 0, 2, 1 });

        (*bias)[i] = u->GetBias(nullptr);
        (*ematrix)[i] = u->GetEmatrixInput(nullptr);
        (*lmatrix)[i] = u->GetL(nullptr);

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

        /* input folds */
        (*shaded)[i] = shade((*victims)[i], mr, { 0, 0, 0, 0 });
        (*side0)[i] = fold<TH1, 4>((*victims)[i], mr, 0, { 0, 0, 0, 0 });
        (*side1)[i] = fold<TH1, 4>((*victims)[i], mr, 1, { 0, 0, 0, 0 });

        /* normalise to unity */
        (*fold0)[i]->Scale(1. / (*fold0)[i]->Integral("width"));
        (*side0)[i]->Scale(1. / (*side0)[i]->Integral("width"));

        /* figures */
        cs[0]->add((*matrices)[i]);
        cs[0]->adjust((*matrices)[i], "colz", "");

        cs[1]->add((*bias)[i], "bias");

        cs[2]->add((*ematrix)[i]);
        cs[2]->adjust((*ematrix)[i], "colz", "");

        cs[3]->add((*lmatrix)[i]);
        cs[3]->adjust((*lmatrix)[i], "colz", "");

        cs[4]->add(hframe);
        cs[4]->stack((*logtaux)[i]);
        cs[4]->adjust((*logtaux)[i], "l", "");
        cs[4]->stack(logtau_opt, "optimal");
        cs[4]->adjust(logtau_opt, "p", "");

        cs[5]->add((*lcurve)[i]);
        cs[5]->adjust((*lcurve)[i], "al", "");
        cs[5]->stack(lcurve_opt, "optimal");
        cs[5]->adjust(lcurve_opt, "p", "");

        cs[6]->add((*result)[i], "unfolded");

        cs[7]->add((*victims)[i], "data");
        cs[7]->stack((*refold)[i], "refolded");

        cs[8]->add((*sresult)[i]);
        cs[8]->adjust((*sresult)[i], "colz", "");
        cs[9]->add((*srefold)[i]);
        cs[9]->adjust((*srefold)[i], "colz", "");

        cs[10]->add((*notes[0])[i], "data", "gen");
        cs[10]->stack((*side0)[i], "data", "reco");
        cs[10]->stack((*fold0)[i], "unfolded", "gen");
        cs[11]->add((*notes[1])[i], "data", "gen");
        cs[11]->stack((*side1)[i], "data", "reco");
        cs[11]->stack((*fold1)[i], "unfolded", "gen");

        cs[12]->add((*shaded)[i]);
        cs[12]->adjust((*shaded)[i], "colz", "");
    });

    hb->set_binary("bins");
    hb->sketch();

    for (auto c : cs)
        c->draw("pdf");

    bias->rename();
    ematrix->rename();
    lmatrix->rename();
    logtaux->rename();
    logtauy->rename();
    lcurve->rename();

    result->rename();
    refold->rename();
    sresult->rename();
    srefold->rename();
    fold0->rename();
    fold1->rename();

    /* save output */
    in(output, [&]() {
        matrices->save("");
        victims->save("");
        shaded->save("");
        side0->save("");
        side1->save("");

        bias->save(tag);
        ematrix->save(tag);
        lmatrix->save(tag);
        logtaux->save(tag);
        logtauy->save(tag);
        lcurve->save(tag);

        result->save(tag);
        refold->save(tag);
        sresult->save(tag);
        srefold->save(tag);
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
