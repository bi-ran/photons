#ifndef SPECIFICS_H
#define SPECIFICS_H

template <typename T>
bool within_hem_failure_region(T* t, int64_t index) {
    return ((*t->phoSCEta)[index] < -1.3
        && (*t->phoSCPhi)[index] < -0.87
        && (*t->phoSCPhi)[index] > -1.57);
}

template <typename T>
bool passes_basic_electron_selections(T* t, int64_t index) {
    return (*t->eleMissHits)[index] <= 1 && (*t->eleIP3D)[index] < 0.03;
}

enum det { barrel, endcap };
enum wp { veto, loose, medium, tight, claustrophobic };

#include "pjtree.h"

template <det T, wp U>
bool passes_electron_id(pjtree* t, int64_t index, bool heavyion);

template <>
bool passes_electron_id<det::barrel, wp::loose>(
        pjtree* t, int64_t index, bool heavyion) {
    if (!passes_basic_electron_selections(t, index)) { return false; }
    if (!(std::abs((*t->eleSCEta)[index]) < 1.442)) { return false; }

    if (heavyion) {
        if (t->hiBin < 60) {
            return (*t->eleHoverEBc)[index] < 0.1616
                && (*t->eleSigmaIEtaIEta_2012)[index] < 0.0135
                && std::abs((*t->eledEtaSeedAtVtx)[index]) < 0.0038
                && std::abs((*t->eledPhiAtVtx)[index]) < 0.0376
                && std::abs((*t->eleEoverPInv)[index]) < 0.0177;
        }

        return (*t->eleHoverEBc)[index] < 0.1268
            && (*t->eleSigmaIEtaIEta_2012)[index] < 0.0107
            && std::abs((*t->eledEtaSeedAtVtx)[index]) < 0.0035
            && std::abs((*t->eledPhiAtVtx)[index]) < 0.0327
            && std::abs((*t->eleEoverPInv)[index]) < 0.0774;
    } else {
        return (*t->eleHoverE)[index] < 0.02711
            && (*t->eleSigmaIEtaIEta_2012)[index] < 0.01016
            && std::abs((*t->eledEtaSeedAtVtx)[index]) < 0.00316
            && std::abs((*t->eledPhiAtVtx)[index]) < 0.03937
            && std::abs((*t->eleEoverPInv)[index]) < 0.05304;
    }
}

#endif /* SPECIFICS_H */
