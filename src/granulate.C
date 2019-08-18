#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/history.h"

#include "../git/tricks-and-treats/include/zip.h"

#include "TFile.h"
#include "TH1.h"

#include <functional>
#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

int granulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto file = conf->get<std::string>("file");
    auto label = conf->get<std::string>("label");
    auto tag = conf->get<std::string>("tag");

    auto refs = conf->get<std::vector<std::string>>("refs");
    auto vars = conf->get<std::vector<std::string>>("vars");
    auto lrefs = conf->get<std::vector<std::string>>("lrefs");
    auto lvars = conf->get<std::vector<std::string>>("lvars");
    auto values = conf->get<std::vector<float>>("values");
    auto figures = conf->get<std::vector<std::string>>("figures");

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load input files */
    TFile* f = new TFile(file.data(), "read");

    std::vector<TFile*> frefs(refs.size(), nullptr);
    std::vector<TFile*> fvars(vars.size(), nullptr);
    zip([&](auto& fref, auto& fvar, auto const& ref, auto const& var) {
        fref = new TFile(ref.data(), "read");
        fvar = new TFile(var.data(), "read");
    }, frefs, fvars, refs, vars);

    /* prepare output */
    TFile* fout = new TFile(output, "recreate");

    /* calculate variations */
    zip([&](auto const& figure) {
        auto stub = "_"s + figure;

        auto base = new history(f, tag + "_"s + label + stub);

        std::vector<history*> references(refs.size(), nullptr);
        std::vector<history*> variations(vars.size(), nullptr);

        zip([&](auto& ref, auto& var, auto fref, auto fvar,
                auto const& lref, auto const& lvar, auto value) {
            ref = new history(fref, tag + "_"s + lref + stub);
            var = new history(fvar, tag + "_"s + lvar + stub);

            var->apply([&](TH1* h, int64_t index) {
                h->Divide((*ref)[index]); });

            /* scale uncertainties */
            if (value != 0) { var->apply([&](TH1* h) {
                for_contents([&](std::array<double, 1> val) -> float {
                    return 1. + value * (val[0] - 1.); }, h); }); }

            /* apply uncertainty to base */
            var->apply([&](TH1* h, int64_t index) {
                h->Multiply((*base)[index]); });

            /* save histograms */
            var->save(tag + "_mod");
        }, references, variations, frefs, fvars, lrefs, lvars, values);
    }, figures);

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return granulate(argv[1], argv[2]);

    return 0;
}
