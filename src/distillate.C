#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;

int distillate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto input = conf->get<std::string>("input");
    auto tag = conf->get<std::string>("tag");

    auto rpt = conf->get<std::vector<float>>("pt_range");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto deta = conf->get<std::vector<float>>("eta_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    auto system = conf->get<std::string>("system");

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* load input */
    TFile* f = new TFile(input.data(), "read");
    auto scale = new history(f, tag);

    /* prepare histograms */
    auto ival = std::make_shared<interval>(1, 0., 1.);

    auto ipt = std::make_shared<interval>(dpt);
    auto ieta = std::make_shared<interval>(deta);
    auto ihf = std::make_shared<interval>(dhf);

    auto mhf = std::make_shared<multival>(dhf);

    auto scale_d_pthf = scale->sum(1);

    auto es = new history("es_"s + tag, "", ival, scale_d_pthf->shape());
    auto er = new history("er_"s + tag, "", ival, scale_d_pthf->shape());

    auto es_f_pt = std::make_unique<history>("es_d_eta_hf"s,
        "reco p_{T}/gen p_{T}", "jet p_{T}", rpt, mhf);
    auto er_f_pt = std::make_unique<history>("er_d_eta_hf"s,
        "#sigma_{p_{T}}/p_{T}", "jet p_{T}", rpt, mhf);

    /* load fitting parameters */
    auto flow = new std::vector<float>[ihf->size()];
    auto fhigh = new std::vector<float>[ihf->size()];
    for (int64_t i = 0; i < ihf->size(); ++i) {
        auto index_str = std::to_string(i);
        flow[i] = conf->get<std::vector<float>>("flow_"s + index_str);
        fhigh[i] = conf->get<std::vector<float>>("fhigh_"s + index_str);
    }

    /* draw plots */
    auto hb = new pencil();
    hb->category("centrality", "0", "1", "2", "3");

    hb->alias("0", "50 - 90%");
    hb->alias("1", "30 - 50%");
    hb->alias("2", "10 - 30%");
    hb->alias("3", "0 - 10%");

    auto c1 = new paper("jesr_fits_"s + tag, hb);
    apply_default_style(c1, system, 0., 1.);
    c1->format(simple_formatter);
    c1->divide(ipt->size(), -1);

    /* fit scale and resolution */
    scale_d_pthf->apply([&](TH1* h, int64_t index) {
        auto indices = scale_d_pthf->indices_for(index);
        auto pt_x = indices[0];
        auto hf_x = indices[1];

        auto label = "f_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "gaus", rpt.front(), rpt.back());
        h->Fit(label.data(), "WLMQ", "", flow[hf_x][pt_x], fhigh[hf_x][pt_x]);

        (*es)[index]->SetBinContent(1, f->GetParameter(1));
        (*es)[index]->SetBinError(1, f->GetParError(1));
        (*er)[index]->SetBinContent(1, f->GetParameter(2));
        (*er)[index]->SetBinError(1, f->GetParError(2));

        ++pt_x;

        (*es_f_pt)[hf_x]->SetBinContent(pt_x, f->GetParameter(1));
        (*es_f_pt)[hf_x]->SetBinError(pt_x, f->GetParError(1));
        (*er_f_pt)[hf_x]->SetBinContent(pt_x, f->GetParameter(2));
        (*er_f_pt)[hf_x]->SetBinError(pt_x, f->GetParError(2));

        c1->add(h, std::to_string(pt_x));
    });

    auto c2 = new paper("jesr_"s + tag, hb);
    apply_default_style(c2, system, 0., 1.);
    c2->format(simple_formatter);
    c2->divide(ihf->size(), -1);
    c2->set(paper::flags::logx);

    es_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(0.5, 1.5, "Y");

        auto label = "f_es_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "[0]+[1]/x+[2]/(x*x)",
            rpt.front(), rpt.back());
        f->SetParameters(1.1, 1.2, 4.8);
        h->Fit(label.data(), "MEQ", "", rpt.front(), rpt.back());

        c2->add(h, std::to_string(index));
    });

    er_f_pt->apply([&](TH1* h, int64_t index) {
        h->SetAxisRange(0., 1., "Y");

        auto label = "f_er_"s + std::to_string(index);
        TF1* f = new TF1(label.data(), "sqrt([0]*[0]+[1]*[1]/x+[2]*[2]/(x*x))",
            rpt.front(), rpt.back());
        f->SetParameters(0.1, 1.2, 4.8);
        h->Fit(label.data(), "MEQ", "", rpt.front(), rpt.back());

        c2->add(h, std::to_string(index));
    });

    hb->sketch();

    c1->draw("pdf");
    c2->draw("pdf");

    /* save output */
    TFile* fout = new TFile(output, "recreate");

    scale_d_pthf->save("pthf");

    es->save("jes");
    er->save("jer");

    es_f_pt->save("jes_f_pt");
    er_f_pt->save("jer_f_pt");

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return distillate(argv[1], argv[2]);

    return 0;
}
