#ifndef SPECIFICS_H
#define SPECIFICS_H

template <typename T>
bool within_hem_failure_region(T* t, int64_t i) {
    return ((*t->phoSCEta)[i] < -1.3
        && (*t->phoSCPhi)[i] < -0.87
        && (*t->phoSCPhi)[i] > -1.57);
}

template <typename T>
bool passes_basic_electron_selections(T* t, int64_t i) {
    return (*t->eleMissHits)[i] <= 1 && (*t->eleIP3D)[i] < 0.03;
}

enum ip { incl, cent, peri, nip };
enum det { barrel, endcap, ndet };
enum wp { veto, loose, medium, tight, claustro, nwp };
enum var { hoe, see, deta, dphi, eop, nele,
        /* hoe, see, */ iso = 2, npho};

constexpr float ecuts[ip::nip][det::ndet][wp::nwp][var::nele] = {
    {    /* ip::incl */
        {    /* det::barrel */
            {    -1.,    -1.,    -1.,    -1.,    -1. }, /* wp::veto */
            {    -1.,    -1.,    -1.,    -1.,    -1. }, /* wp::loose */
            {    -1.,    -1.,    -1.,    -1.,    -1. }, /* wp::medium */
            {    -1.,    -1.,    -1.,    -1.,    -1. }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }, { /* det::endcap */
            {    -1.,    -1.,    -1.,    -1.,    -1. }, /* wp::veto */
            {    -1.,    -1.,    -1.,    -1.,    -1. }, /* wp::loose */
            {    -1.,    -1.,    -1.,    -1.,    -1. }, /* wp::medium */
            {    -1.,    -1.,    -1.,    -1.,    -1. }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }
    }, { /* ip::cent */
        {    /* det::barrel */
            { 0.2733, 0.0147, 0.0041, 0.0853, 0.0367 }, /* wp::veto */
            { 0.1616, 0.0135, 0.0038, 0.0376, 0.0177 }, /* wp::loose */
            { 0.1589, 0.0116, 0.0037, 0.0224, 0.0173 }, /* wp::medium */
            { 0.1459, 0.0104, 0.0029, 0.0206, 0.0105 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }, { /* det::endcap */
            { 0.1898, 0.0480, 0.0097, 0.2348, 0.0300 }, /* wp::veto */
            { 0.1317, 0.0466, 0.0063, 0.1186, 0.0201 }, /* wp::loose */
            { 0.1092, 0.0418, 0.0062, 0.0373, 0.0133 }, /* wp::medium */
            { 0.0925, 0.0358, 0.0051, 0.0266, 0.0065 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }
    }, { /* ip::peri */
        {    /* det::barrel */
            { 0.1814, 0.0113, 0.0037, 0.1280, 0.1065 }, /* wp::veto */
            { 0.1268, 0.0107, 0.0035, 0.0327, 0.0774 }, /* wp::loose */
            { 0.0311, 0.0101, 0.0033, 0.0210, 0.0701 }, /* wp::medium */
            { 0.0067, 0.0099, 0.0026, 0.0170, 0.0077 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }, { /* det::endcap */
            { 0.1138, 0.0376, 0.0074, 0.2085, 0.0237 }, /* wp::veto */
            { 0.0977, 0.0339, 0.0067, 0.0838, 0.0193 }, /* wp::loose */
            { 0.0810, 0.0316, 0.0051, 0.0384, 0.0192 }, /* wp::medium */
            { 0.0655, 0.0288, 0.0044, 0.0266, 0.0123 }, /* wp::tight */
            {    -1.,    -1.,    -1.,    -1.,    -1. }  /* wp::claustro */
        }
    }
};

template <det T, wp U, typename V>
bool passes_electron_id(V* t, int64_t i, bool heavyion) {
    if (!passes_basic_electron_selections(t, i)) { return false; }

    auto iptype = heavyion ? (t->hiBin < 60 ? ip::cent : ip::peri) : ip::incl;
    return (*t->eleHoverEBc)[i] < ecuts[iptype][T][U][var::hoe]
        && (*t->eleSigmaIEtaIEta_2012)[i] < ecuts[iptype][T][U][var::see]
        && std::abs((*t->eledEtaSeedAtVtx)[i]) < ecuts[iptype][T][U][var::deta]
        && std::abs((*t->eledPhiAtVtx)[i]) < ecuts[iptype][T][U][var::dphi]
        && std::abs((*t->eleEoverPInv)[i]) < ecuts[iptype][T][U][var::eop];
}

#endif /* SPECIFICS_H */
