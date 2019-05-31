#ifndef TRACKTREE_H
#define TRACKTREE_H

#include "TTree.h"

#include <vector>

#include "defines.h"

#define B_ALT_D(ACTION, ...)                                                \
    B_AVT_D(ACTION, ## __VA_ARGS__)                                         \
    B_AAT_D(ACTION, ## __VA_ARGS__)                                         \

#define B_AVT_D(ACTION, ...)                                                \
    ACTION(int32_t, nTrk, ## __VA_ARGS__)                                   \

#define B_AAT_D(ACTION, ...)                                                \
    ACTION(float, trkPt, ## __VA_ARGS__)                                    \
    ACTION(float, trkPtError, ## __VA_ARGS__)                               \
    ACTION(uint8_t, trkNHit, ## __VA_ARGS__)                                \
    ACTION(uint8_t, trkNlayer, ## __VA_ARGS__)                              \
    ACTION(float, trkEta, ## __VA_ARGS__)                                   \
    ACTION(float, trkPhi, ## __VA_ARGS__)                                   \
    ACTION(int32_t, trkCharge, ## __VA_ARGS__)                              \
    ACTION(bool, highPurity, ## __VA_ARGS__)                                \
    ACTION(float, trkChi2, ## __VA_ARGS__)                                  \
    ACTION(uint8_t, trkNdof, ## __VA_ARGS__)                                \
    ACTION(float, trkDxy1, ## __VA_ARGS__)                                  \
    ACTION(float, trkDxyError1, ## __VA_ARGS__)                             \
    ACTION(float, trkDz1, ## __VA_ARGS__)                                   \
    ACTION(float, trkDzError1, ## __VA_ARGS__)                              \
    ACTION(bool, trkFake, ## __VA_ARGS__)                                   \
    ACTION(int32_t, pfType, ## __VA_ARGS__)                                 \
    ACTION(float, pfCandPt, ## __VA_ARGS__)                                 \
    ACTION(float, pfEcal, ## __VA_ARGS__)                                   \
    ACTION(float, pfHcal, ## __VA_ARGS__)                                   \

class pjtree;

class tracktree {
    friend pjtree;

    public:
        tracktree() {
            B_AVT_D(ZERO)
            B_AAT_D(NEWARR, size)
        };

        tracktree(TTree* t) {
            read(t);
        }

        ~tracktree() = default;

        void read(TTree* t) {
            B_AVT_D(RREF, t)
            B_AAT_D(RVAR, t)
        };

    private:
        B_AVT_D(DECLARE)
        B_AAT_D(DECLPTR)
};

#endif  /* TRACKTREE_H */
