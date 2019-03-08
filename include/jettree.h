#ifndef JETTREE_H
#define JETTREE_H

#include "TTree.h"

#include <vector>

#include "defines.h"

#define B_ALJ_D(ACTION, ...)                                                \
    B_AVJ_D(ACTION, ## __VA_ARGS__)                                         \
    B_AAJ_D(ACTION, ## __VA_ARGS__)                                         \

#define B_AVJ_D(ACTION, ...)                                                \
    ACTION(int, nref, ## __VA_ARGS__)                                       \

#define B_AAJ_D(ACTION, ...)                                                \
    ACTION(float, rawpt, ## __VA_ARGS__)                                    \
    ACTION(float, jtpt, ## __VA_ARGS__)                                     \
    ACTION(float, jteta, ## __VA_ARGS__)                                    \
    ACTION(float, jtphi, ## __VA_ARGS__)                                    \

#define B_ALJ_M(ACTION, ...)                                                \
    B_AVJ_M(ACTION, ## __VA_ARGS__)                                         \
    B_AAJ_M(ACTION, ## __VA_ARGS__)                                         \

#define B_AVJ_M(ACTION, ...)                                                \
    ACTION(int, ngen, ## __VA_ARGS__)                                       \

#define B_AAJ_M(ACTION, ...)                                                \
    ACTION(float, genpt, ## __VA_ARGS__)                                    \
    ACTION(float, geneta, ## __VA_ARGS__)                                   \
    ACTION(float, genphi, ## __VA_ARGS__)                                   \

static constexpr int size = 1000;

class pjtree;

class jettree {
    friend pjtree;

    public:
        jettree() {
            this->mc_branches = 0;
            B_AVJ_D(ZERO)
            B_AVJ_M(ZERO)
            B_AAJ_D(NEWARR, size)
            B_AAJ_M(NEWARR, size)
        };

        jettree(TTree* t, bool mc_branches)
                : jettree() {
            this->mc_branches = mc_branches;
            read(t);
        }

        ~jettree() = default;

        void read(TTree* t) {
            B_AVJ_D(RREF, t)
            B_AAJ_D(RVAR, t)
            if (mc_branches) {
                B_AVJ_M(RREF, t)
                B_AAJ_M(RVAR, t) }
        };

    private:
        bool mc_branches;

        B_AVJ_D(DECLARE)
        B_AVJ_M(DECLARE)
        B_AAJ_D(DECLPTR)
        B_AAJ_M(DECLPTR)
};

#endif  /* JETTREE_H */
