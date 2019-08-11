#include "../include/lambdas.h"
#include "../include/specifics.h"

#include "../git/config/include/configurer.h"

#include "../git/foliage/include/event.h"
#include "../git/foliage/include/photons.h"
#include "../git/foliage/include/triggers.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/train.h"

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TTree.h"

#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

int speculate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto files = conf->get<std::vector<std::string>>("files");
    auto max_entries = conf->get<int64_t>("max_entries");
    auto system = conf->get<std::string>("system");
    auto tag = conf->get<std::string>("tag");

    auto heavyion = conf->get<bool>("heavyion");

    auto selections = conf->get<std::vector<std::string>>("selections");
    auto paths = conf->get<std::vector<std::string>>("paths");

    auto const eta_abs = conf->get<float>("eta_abs");
    auto const hovere_max = conf->get<float>("hovere_max");
    auto const see_min = conf->get<float>("see_min");
    auto const see_max = conf->get<float>("see_max");
    auto const iso_max = conf->get<float>("iso_max");

    auto rpt = conf->get<std::vector<float>>("pt_range");

    /* load forest */
    auto forest = new train(files);
    auto chain_evt = forest->attach("hiEvtAnalyzer/HiTree", true);
    auto chain_sel = forest->attach("skimanalysis/HltTree", true);
    auto chain_eg = forest->attach("ggHiNtuplizerGED/EventTree", true);
    auto chain_hlt = forest->attach("hltanalysis/HltTree", true);

    (*forest)();

    auto tevt = harvest<event>(chain_evt, false);
    auto tsel = harvest<triggers>(chain_sel, selections);
    auto tpho = harvest<photons>(chain_eg);
    auto thlt = harvest<triggers>(chain_hlt, paths);

    auto counts = new history("count", "counts", "photon p_{T}", rpt, 2);

    /* iterate */
    int64_t nentries = forest->count();
    if (max_entries) nentries = std::min(nentries, max_entries);
    for (int64_t i = 0; i < nentries; ++i) {
        if (i % 100000 == 0)
            printf("entry: %li/%li\n", i, nentries);

        forest->get(i);

        if (std::abs(tevt->vz) > 15) { continue; }

        if (!selections.empty()) {
            bool pass_selection = false;
            for (auto const& sel : selections)
                if (tsel->accept(sel) == 1)
                    pass_selection = true;

            if (!pass_selection) { continue; }
        }

        int64_t leading = -1;
        for (int64_t j = 0; j < tpho->nPho; ++j) {
            if (std::abs((*tpho->phoSCEta)[j]) >= eta_abs) { continue; }
            if ((*tpho->phoHoverE)[j] > hovere_max) { continue; }

            leading = j;
            break;
        }

        /* require leading photon */
        if (leading < 0) { continue; }

        if ((*tpho->phoSigmaIEtaIEta_2012)[leading] > see_max
                || (*tpho->phoSigmaIEtaIEta_2012)[leading] < see_min)
            continue;

        /* hem failure region exclusion */
        if (heavyion && within_hem_failure_region(tpho, leading)) { continue; }

        /* isolation requirement */
        float isolation = (*tpho->pho_ecalClusterIsoR3)[leading]
            + (*tpho->pho_hcalRechitIsoR3)[leading]
            + (*tpho->pho_trackIsoR3PtCut20)[leading];
        if (isolation > iso_max) { continue; }

        auto et = (*tpho->phoEt)[leading];

        (*counts)[0]->Fill(et);
        for (auto const& path : paths) {
            if (thlt->accept(path) == 1) {
                (*counts)[1]->Fill(et);
                break;
            }
        }
    }

    /* calculate efficiency */
    auto frame = (TH1F*)(*counts)[0]->Clone("frame");
    frame->GetYaxis()->SetTitle("efficiency");
    frame->Reset("MICES");

    auto eff = new TGraphAsymmErrors((*counts)[1], (*counts)[0],
        "c1=0.683 b(1,1) mode");

    /* draw efficiency */
    auto hb = new pencil();
    hb->category("sample", "pp", "aa");

    hb->alias("aa", "PbPb");

    auto c1 = new paper(tag + "_efficiency", hb);
    apply_default_style(c1, system + " #sqrt{s} = 5.02 TeV"s, 0., 1.2);
    c1->accessory(std::bind(line_at, _1, 1., rpt.front(), rpt.back()));

    c1->add(frame);
    c1->stack(eff, tag);

    hb->sketch();
    
    c1->draw("pdf");

    /* save output */
    TFile* fout = new TFile(output, "recreate");

    counts->save(tag);

    fout->Write("", TObject::kOverwrite);
    fout->Close();
    
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return speculate(argv[1], argv[2]);

    return 0;
}
