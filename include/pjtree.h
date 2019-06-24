#ifndef PJTREE_H
#define PJTREE_H

#include "../git/foliage/include/foliage.h"
#include "../git/foliage/include/photons.h"
#include "../git/foliage/include/jets.h"
#include "../git/foliage/include/tracks.h"
#include "../git/foliage/include/triggers.h"

#include "TTree.h"

#include <vector>

#define B_VEC_JET_RECO(ACTION, ...)                                         \
    ACTION(sv<float>,       rawpt,                      ## __VA_ARGS__)     \
    ACTION(sv<float>,       jtpt,                       ## __VA_ARGS__)     \
    ACTION(sv<float>,       jteta,                      ## __VA_ARGS__)     \
    ACTION(sv<float>,       jtphi,                      ## __VA_ARGS__)     \

#define B_VEC_JET_GEN(ACTION, ...)                                          \
    ACTION(sv<float>,       genpt,                      ## __VA_ARGS__)     \
    ACTION(sv<float>,       geneta,                     ## __VA_ARGS__)     \
    ACTION(sv<float>,       genphi,                     ## __VA_ARGS__)     \

#define B_VEC_TRK_RECO(ACTION, ...)                                         \
    ACTION(sv<float>,       trkPt,                      ## __VA_ARGS__)     \
    ACTION(sv<float>,       trkPtError,                 ## __VA_ARGS__)     \
    ACTION(sv<uint8_t>,     trkNHit,                    ## __VA_ARGS__)     \
    ACTION(sv<uint8_t>,     trkNlayer,                  ## __VA_ARGS__)     \
    ACTION(sv<float>,       trkEta,                     ## __VA_ARGS__)     \
    ACTION(sv<float>,       trkPhi,                     ## __VA_ARGS__)     \
    ACTION(sv<int32_t>,     trkCharge,                  ## __VA_ARGS__)     \
    ACTION(sv<bool>,        highPurity,                 ## __VA_ARGS__)     \
    ACTION(sv<float>,       trkChi2,                    ## __VA_ARGS__)     \
    ACTION(sv<uint8_t>,     trkNdof,                    ## __VA_ARGS__)     \
    ACTION(sv<float>,       trkDxy1,                    ## __VA_ARGS__)     \
    ACTION(sv<float>,       trkDxyError1,               ## __VA_ARGS__)     \
    ACTION(sv<float>,       trkDz1,                     ## __VA_ARGS__)     \
    ACTION(sv<float>,       trkDzError1,                ## __VA_ARGS__)     \
    ACTION(sv<bool>,        trkFake,                    ## __VA_ARGS__)     \
    ACTION(sv<int32_t>,     pfType,                     ## __VA_ARGS__)     \
    ACTION(sv<float>,       pfCandPt,                   ## __VA_ARGS__)     \
    ACTION(sv<float>,       pfEcal,                     ## __VA_ARGS__)     \
    ACTION(sv<float>,       pfHcal,                     ## __VA_ARGS__)     \

#define B_VEC_TRG(ACTION, ...)                                              \
    ACTION(sv<int32_t>,     accepts,                    ## __VA_ARGS__)     \

class pjtree {
  public:
    pjtree(TTree* t, bool gen, bool hlt)
            : _gen(gen),
              _hlt(hlt) {
        B_VAL_PHO_RECO(SETMONE)
        B_VAL_JET_RECO(SETMONE)
        B_VAL_TRK_RECO(SETMONE)
        B_VEC_PHO_RECO(ALLOCOBJ)
        B_VEC_JET_RECO(ALLOCOBJ)
        B_VEC_TRK_RECO(ALLOCOBJ)

        if (_gen) {
            B_VAL_PHO_GEN(SETMONE)
            B_VAL_JET_GEN(SETMONE)
            B_VEC_PHO_GEN(ALLOCOBJ)
            B_VEC_JET_GEN(ALLOCOBJ)
        }

        if (_hlt) {
            B_VEC_TRG(ALLOCOBJ)
        }

        branch(t);
    }

    pjtree(TTree* t, bool gen)
        : pjtree(t, gen, false) { }

    pjtree(bool gen, bool hlt, TTree* t)
            : _gen(gen),
              _hlt(hlt) {
        B_VAL_PHO_RECO(SETZERO)
        B_VAL_JET_RECO(SETZERO)
        B_VAL_TRK_RECO(SETZERO)
        B_VEC_PHO_RECO(SETZERO)
        B_VEC_JET_RECO(SETZERO)
        B_VEC_TRK_RECO(SETZERO)

        if (_gen) {
            B_VAL_PHO_GEN(SETZERO)
            B_VAL_JET_GEN(SETZERO)
            B_VEC_PHO_GEN(SETZERO)
            B_VEC_JET_GEN(SETZERO)
        }

        if (_hlt) {
            B_VEC_TRG(SETZERO)
        }

        read(t);
    }

    pjtree(bool gen, TTree* t)
        : pjtree(gen, false, t) { }

    ~pjtree() = default;

    void clear() {
        B_VEC_PHO_RECO(CLEAROBJ)
        B_VEC_JET_RECO(CLEAROBJ)
        B_VEC_TRK_RECO(CLEAROBJ)

        if (_gen) {
            B_VEC_PHO_GEN(CLEAROBJ)
            B_VEC_JET_GEN(CLEAROBJ)
        }

        if (_hlt) {
            B_VEC_TRG(CLEAROBJ)
        }
    }

    void copy(photons* t) {
        B_VAL_PHO_RECO(COPYVAL, t)
        B_VEC_PHO_RECO(COPYOBJ, t)

        if (_gen) {
            B_VAL_PHO_GEN(COPYVAL, t)
            B_VEC_PHO_GEN(COPYOBJ, t)
        }
    }

    void copy(jets* t) {
        B_VAL_JET_RECO(COPYVAL, t)
        B_VEC_JET_RECO(COPYPTR, t, nref)

        if (_gen) {
            B_VAL_JET_GEN(COPYVAL, t)
            B_VEC_JET_GEN(COPYPTR, t, ngen)
        }
    }

    void copy(tracks* t) {
        B_VAL_TRK_RECO(COPYVAL, t)
        B_VEC_TRK_RECO(COPYPTR, t, nTrk)
    }

    void copy(triggers* t) {
        if (_hlt) {
            B_VEC_TRG(COPYPTR, t, t->size())
        }
    }

    B_VAL_PHO_RECO(DECLVAL)
    B_VEC_PHO_RECO(DECLPTR)
    B_VAL_JET_RECO(DECLVAL)
    B_VEC_JET_RECO(DECLPTR)
    B_VAL_TRK_RECO(DECLVAL)
    B_VEC_TRK_RECO(DECLPTR)
    B_VAL_PHO_GEN(DECLVAL)
    B_VEC_PHO_GEN(DECLPTR)
    B_VAL_JET_GEN(DECLVAL)
    B_VEC_JET_GEN(DECLPTR)
    B_VEC_TRG(DECLPTR)

  private:
    void branch(TTree* t) {
        B_VAL_PHO_RECO(BRANCHVAL, t)
        B_VAL_JET_RECO(BRANCHVAL, t)
        B_VAL_TRK_RECO(BRANCHVAL, t)
        B_VEC_PHO_RECO(BRANCHPTR, t)
        B_VEC_JET_RECO(BRANCHPTR, t)
        B_VEC_TRK_RECO(BRANCHPTR, t)

        if (_gen) {
            B_VAL_PHO_GEN(BRANCHVAL, t)
            B_VAL_JET_GEN(BRANCHVAL, t)
            B_VEC_PHO_GEN(BRANCHPTR, t)
            B_VEC_JET_GEN(BRANCHPTR, t)
        }

        if (_hlt) {
            B_VEC_TRG(BRANCHPTR, t)
        }
    }

    void read(TTree* t) {
        B_VAL_PHO_RECO(SETVALADDR, t)
        B_VAL_JET_RECO(SETVALADDR, t)
        B_VAL_TRK_RECO(SETVALADDR, t)
        B_VEC_PHO_RECO(SETVALADDR, t)
        B_VEC_JET_RECO(SETVALADDR, t)
        B_VEC_TRK_RECO(SETVALADDR, t)

        if (_gen) {
            B_VAL_PHO_GEN(SETVALADDR, t)
            B_VAL_JET_GEN(SETVALADDR, t)
            B_VEC_PHO_GEN(SETVALADDR, t)
            B_VEC_JET_GEN(SETVALADDR, t)
        }

        if (_hlt) {
            B_VEC_TRG(SETVALADDR, t)
        }
    }

    bool _gen;
    bool _hlt;
};

#endif /* PJTREE_H */
