#include "../include/lambdas.h"
#include "../include/pjtree.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TChain.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;

int emulate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto data = conf->get<std::string>("data");
    auto files = conf->get<std::vector<std::string>>("files");
    auto pthats = conf->get<std::vector<int32_t>>("pthats");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto rpthat = conf->get<std::vector<float>>("pthat_range");
    auto rvz = conf->get<std::vector<float>>("vz_range");

    auto ipthat = std::make_shared<interval>("pthat"s,
        static_cast<int64_t>(rpthat[0]), rpthat[1], rpthat[2]);
    auto ivz = std::make_shared<interval>("v_{z}"s,
        static_cast<int64_t>(rvz[0]), rvz[1], rvz[2]);

    /* manange memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* merged gen inputs */
    TChain* tcomb = new TChain("pj");
    for (auto const& f : files)
        tcomb->Add(f.data());

    /* calculate pthat weights */
    TChain* tbase = new TChain("pj");
    tbase->Add(files[0].data());

    auto count = static_cast<int64_t>(pthats.size());
    auto incl = std::make_shared<interval>(1L, 0., 1.);
    auto pthat = std::make_unique<history>("pthat"s, "", incl, count);

    pthats.push_back(999999);
    for (int64_t i = 0; i < count; ++i) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "pthat>=%i&&pthat<%i", pthats[i], pthats[i + 1]);
        float nbase = tbase->GetEntries(buffer);
        float ncomb = tcomb->GetEntries(buffer);

        auto weight = nbase / ncomb;
        (*pthat)[i]->SetBinContent(1, weight);
        printf("[%i, %i]: %f\n", pthats[i], pthats[i + 1], weight);
    }

    printf("\n");

    /* calculate vz weights */
    auto vz = std::make_unique<history>("vz"s, "", ivz, 3L);

    TChain* tdata = new TChain("pj");
    tdata->Add(data.data());

    tdata->Project((*vz)[0]->GetName(), "vz");
    tcomb->Project((*vz)[1]->GetName(), "vz");

    (*vz)[0]->Scale(1. / (*vz)[0]->Integral());
    (*vz)[1]->Scale(1. / (*vz)[1]->Integral());

    (*vz)[0]->Fit("gaus", "LMQ", "", -15, 15);
    (*vz)[1]->Fit("gaus", "LMQ", "", -15, 15);

    auto fdata = (*vz)[0]->GetFunction("gaus");
    printf("data: (%f, %f, %f)\n", fdata->GetParameter(0),
        fdata->GetParameter(1), fdata->GetParameter(2));

    auto fcomb = (*vz)[1]->GetFunction("gaus");
    printf("gen: (%f, %f, %f)\n", fcomb->GetParameter(0),
        fcomb->GetParameter(1), fcomb->GetParameter(2));

    (*vz)[2]->Divide((*vz)[0], (*vz)[1]);

    TF1* fweight = new TF1("fweight", "(gaus(0))/(gaus(3))", -15, 15);
    fweight->SetParameters(
        fdata->GetParameter(0), fdata->GetParameter(1), fdata->GetParameter(2),
        fcomb->GetParameter(0), fcomb->GetParameter(1), fcomb->GetParameter(2)
    );

    /* draw vertex distributions */
    auto hb = new pencil();
    hb->category("system", "data", "mc");

    auto system_tag = system + " #sqrt{s_{NN}} = 5.02 TeV"s;

    auto c1 = new paper(tag + "_vz", hb);
    apply_style(c1, system_tag, 0., 0.04);
    c1->jewellery([&](TH1* h, int64_t index) {
        if (index == 3) {
            fweight->Draw("same");
            h->SetAxisRange(0.4, 1.6, "Y");
        }
    });

    c1->divide(3, -1);

    c1->add((*vz)[0], "data");
    c1->add((*vz)[1], "mc");
    c1->add((*vz)[2]);

    hb->sketch();
    c1->draw("pdf");

    /* save output */
    TFile* fout = new TFile(output, "recreate");

    pthat->save(tag);
    vz->save(tag);

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return emulate(argv[1], argv[2]);

    return 0;
}
