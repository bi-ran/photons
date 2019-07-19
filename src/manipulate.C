#include "../git/config/include/configurer.h"

#include "../git/history/include/history.h"

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
    auto plabel = conf->get<std::string>("plabel");
    auto inputs = conf->get<std::vector<std::string>>("inputs");
    auto groups = conf->get<std::vector<std::string>>("groups");
    auto labels = conf->get<std::vector<std::string>>("labels");
    auto dims = conf->get<std::vector<int32_t>>("dims");
    auto tag = conf->get<std::string>("tag");

    if (inputs.empty()) { return 1; }
    if (inputs.size() != groups.size()) { return 1; }
    auto size = static_cast<int64_t>(inputs.size());

    std::vector<TFile*> files(size, nullptr);
    for (int64_t i = 0; i < size; ++i)
        files[i] = new TFile(inputs[i].data(), "read");

    TFile* fp = new TFile(purity.data(), "read");

    TFile* fout = new TFile(output, "recreate");

    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load purity histograms */
    auto purities = new history(fp, plabel);
    auto impurities = new history(*purities, "im"s);

    purities->apply([](TH1* h) {
        auto purity = h->GetBinContent(1);
        h->SetBinContent(1, 1. / purity);
    });
    impurities->apply([](TH1* h) {
        auto purity = h->GetBinContent(1);
        h->SetBinContent(1, 1. - 1. / purity);
    });

    /* load histograms and perform purity subtraction */
    for (int64_t i = 0; i < static_cast<int64_t>(labels.size()); ++i) {
        auto const& label = labels[i];

        std::vector<history*> histograms;
        for (int64_t j = 0; j < size; ++j) {
            auto name = groups[j] + "_"s + label;
            histograms.emplace_back(new history(files[j], name));

            auto factors = groups[j].find("raw") != std::string::npos
                ? purities : impurities;

            /* scale by appropriate (im)purity values */
            auto pdim = factors->dims();
            if (dims[i] == pdim) {
                histograms.back()->multiply(*factors);
            } else if (dims[i] > pdim) {
                std::vector<int64_t> axes(dims[i] - pdim);
                std::iota(std::begin(axes), std::end(axes), pdim);
                histograms.back()->multiply(*factors, axes);
            } else {
                printf("error! dimension mismatch\n");
                return 1;
            }
        }

        /* subtract scaled histograms */
        auto result = histograms.front();
        for (int64_t j = 1; j < size; ++j)
            result->add(*histograms[j], 1);

        result->save(tag);
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
