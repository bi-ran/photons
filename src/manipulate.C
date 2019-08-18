#include "../git/config/include/configurer.h"

#include "../git/history/include/history.h"

#include "../git/tricks-and-treats/include/zip.h"

#include <string>
#include <unordered_map>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TNamed.h"

using namespace std::literals::string_literals;

int manipulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto purity = conf->get<std::string>("purity");
    auto type = conf->get<std::string>("type");
    auto inputs = conf->get<std::vector<std::string>>("inputs");
    auto groups = conf->get<std::vector<std::string>>("groups");
    auto labels = conf->get<std::vector<std::string>>("labels");
    auto tag = conf->get<std::string>("tag");

    /* manage meory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load input */
    TFile* fp = new TFile(purity.data(), "read");

    auto purities = new history(fp, type);
    auto impurities = new history(*purities, "im"s);

    purities->apply([](TH1* h) {
        auto purity = h->GetBinContent(1);
        h->SetBinContent(1, 1. / purity);
    });
    impurities->apply([](TH1* h) {
        auto purity = h->GetBinContent(1);
        h->SetBinContent(1, 1. - 1. / purity);
    });

    std::vector<TFile*> files(inputs.size(), nullptr);
    zip([&](auto& file, auto const& input) {
        file = new TFile(input.data(), "read");
    }, files, inputs);

    /* prepare output */
    TFile* fout = new TFile(output, "recreate");

    /* load histograms and perform purity subtraction */
    for (auto const& label : labels) {
        std::vector<history*> histograms(groups.size(), nullptr);
        zip([&](auto const file, auto const& group, auto& hist) {
            auto name = group + "_"s + label;
            hist = new history(file, name);

            auto factors = group.find("raw") != std::string::npos
                ? purities : impurities;

            /* scale by appropriate (im)purity values */
            hist->multiply(*factors);
        }, files, groups, histograms);

        /* subtract scaled histograms */
        histograms[0]->add(*histograms[1], 1);
        histograms[0]->save(tag);
    }

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return manipulate(argv[1], argv[2]);

    printf("usage: %s [config] [output]\n", argv[0]);
    return 1;
}
