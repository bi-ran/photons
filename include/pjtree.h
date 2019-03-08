#ifndef PJTREE_H
#define PJTREE_H

#include "TTree.h"

#include <vector>

#include "defines.h"

#define B_VAR_D(ACTION, ...)                                                \
    ACTION(UInt_t, run, ## __VA_ARGS__)                                     \
    ACTION(ULong64_t, event, ## __VA_ARGS__)                                \
    ACTION(UInt_t, lumis, ## __VA_ARGS__)                                   \
                                                                            \
    ACTION(Int_t, nPho, ## __VA_ARGS__)                                     \

#define B_VEC_D(ACTION, ...)                                                \
    ACTION(std::vector<float>, phoE, ## __VA_ARGS__)                        \
    ACTION(std::vector<float>, phoEt, ## __VA_ARGS__)                       \
    ACTION(std::vector<float>, phoEta, ## __VA_ARGS__)                      \
    ACTION(std::vector<float>, phoPhi, ## __VA_ARGS__)                      \
                                                                            \
    ACTION(std::vector<float>, phoSCE, ## __VA_ARGS__)                      \
    ACTION(std::vector<float>, phoSCRawE, ## __VA_ARGS__)                   \
    ACTION(std::vector<float>, phoSCEta, ## __VA_ARGS__)                    \
    ACTION(std::vector<float>, phoSCPhi, ## __VA_ARGS__)                    \
    ACTION(std::vector<float>, phoSCEtaWidth, ## __VA_ARGS__)               \
    ACTION(std::vector<float>, phoSCPhiWidth, ## __VA_ARGS__)               \
    ACTION(std::vector<float>, phoSCBrem, ## __VA_ARGS__)                   \
    ACTION(std::vector<int>, phoSCnHits, ## __VA_ARGS__)                    \
    ACTION(std::vector<uint32_t>, phoSCflags, ## __VA_ARGS__)               \
    ACTION(std::vector<int>, phoSCinClean, ## __VA_ARGS__)                  \
    ACTION(std::vector<int>, phoSCinUnClean, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, phoSCnBC, ## __VA_ARGS__)                    \
    ACTION(std::vector<float>, phoESEn, ## __VA_ARGS__)                     \
                                                                            \
    ACTION(std::vector<int>, phoIsPFPhoton, ## __VA_ARGS__)                 \
    ACTION(std::vector<int>, phoIsStandardPhoton, ## __VA_ARGS__)           \
    ACTION(std::vector<int>, phoHasPixelSeed, ## __VA_ARGS__)               \
    ACTION(std::vector<int>, phoHasConversionTracks, ## __VA_ARGS__)        \
                                                                            \
    ACTION(std::vector<float>, phoHoverE, ## __VA_ARGS__)                   \
    ACTION(std::vector<float>, phoR9_2012, ## __VA_ARGS__)                  \
    ACTION(std::vector<float>, phoSigmaEtaEta_2012, ## __VA_ARGS__)         \
    ACTION(std::vector<float>, phoSigmaIEtaIEta_2012, ## __VA_ARGS__)       \
    ACTION(std::vector<float>, phoE1x5_2012, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, phoE2x5_2012, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, phoE3x3_2012, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, phoE5x5_2012, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, phoR1x5_2012, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, phoR2x5_2012, ## __VA_ARGS__)                \
                                                                            \
    ACTION(std::vector<float>, pho_ecalClusterIsoR3, ## __VA_ARGS__)        \
    ACTION(std::vector<float>, pho_ecalClusterIsoR4, ## __VA_ARGS__)        \
    ACTION(std::vector<float>, pho_hcalRechitIsoR3, ## __VA_ARGS__)         \
    ACTION(std::vector<float>, pho_hcalRechitIsoR4, ## __VA_ARGS__)         \
    ACTION(std::vector<float>, pho_trackIsoR3PtCut20, ## __VA_ARGS__)       \
    ACTION(std::vector<float>, pho_trackIsoR4PtCut20, ## __VA_ARGS__)       \
                                                                            \
    ACTION(std::vector<float>, pho_swissCrx, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, pho_seedTime, ## __VA_ARGS__)                \

#define B_VAR_M(ACTION, ...)                                                \
    ACTION(Int_t, nMC, ## __VA_ARGS__)                                      \

#define B_VEC_M(ACTION, ...)                                                \
    ACTION(std::vector<float>, mcVtx_x, ## __VA_ARGS__)                     \
    ACTION(std::vector<float>, mcVtx_y, ## __VA_ARGS__)                     \
    ACTION(std::vector<float>, mcVtx_z, ## __VA_ARGS__)                     \
                                                                            \
    ACTION(std::vector<int>, mcPID, ## __VA_ARGS__)                         \
    ACTION(std::vector<int>, mcStatus, ## __VA_ARGS__)                      \
    ACTION(std::vector<float>, mcPt, ## __VA_ARGS__)                        \
    ACTION(std::vector<float>, mcEta, ## __VA_ARGS__)                       \
    ACTION(std::vector<float>, mcPhi, ## __VA_ARGS__)                       \
    ACTION(std::vector<float>, mcE, ## __VA_ARGS__)                         \
    ACTION(std::vector<float>, mcEt, ## __VA_ARGS__)                        \
    ACTION(std::vector<float>, mcMass, ## __VA_ARGS__)                      \
                                                                            \
    ACTION(std::vector<int>, mcParentage, ## __VA_ARGS__)                   \
    ACTION(std::vector<int>, mcMomPID, ## __VA_ARGS__)                      \
    ACTION(std::vector<float>, mcMomPt, ## __VA_ARGS__)                     \
    ACTION(std::vector<float>, mcMomEta, ## __VA_ARGS__)                    \
    ACTION(std::vector<float>, mcMomPhi, ## __VA_ARGS__)                    \
    ACTION(std::vector<float>, mcMomMass, ## __VA_ARGS__)                   \
    ACTION(std::vector<int>, mcGMomPID, ## __VA_ARGS__)                     \
                                                                            \
    ACTION(std::vector<float>, mcCalIsoDR03, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, mcCalIsoDR04, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, mcTrkIsoDR03, ## __VA_ARGS__)                \
    ACTION(std::vector<float>, mcTrkIsoDR04, ## __VA_ARGS__)                \

#define B_VAR_J(ACTION, ...)                                                \
    ACTION(int, nref, ## __VA_ARGS__)                                       \

#define B_ARR_J(ACTION, ...)                                                \
    ACTION(std::vector<float>, rawpt, ## __VA_ARGS__)                       \
    ACTION(std::vector<float>, jtpt, ## __VA_ARGS__)                        \
    ACTION(std::vector<float>, jteta, ## __VA_ARGS__)                       \
    ACTION(std::vector<float>, jtphi, ## __VA_ARGS__)                       \

#define B_VAR_G(ACTION, ...)                                                \
    ACTION(int, ngen, ## __VA_ARGS__)                                       \

#define B_ARR_G(ACTION, ...)                                                \
    ACTION(std::vector<float>, genpt, ## __VA_ARGS__)                       \
    ACTION(std::vector<float>, geneta, ## __VA_ARGS__)                      \
    ACTION(std::vector<float>, genphi, ## __VA_ARGS__)                      \

class photontree;
class jettree;

class pjtree {
    public:
        pjtree() {
            this->mc_branches = 0;
            B_VAR_D(INVALID)
            B_VAR_J(INVALID)
            B_VAR_M(INVALID)
            B_VAR_G(INVALID)
            B_VEC_D(NEWVEC)
            B_ARR_J(NEWVEC)
            B_VEC_M(NEWVEC)
            B_ARR_G(NEWVEC)
        };

        pjtree(TTree* t, bool mc_branches)
                : pjtree() {
            this->mc_branches = mc_branches;
            branch(t);
        };

        ~pjtree() = default;

        void branch(TTree* t) {
            B_VAR_D(BRNREF, t)
            B_VAR_J(BRNREF, t)
            B_VEC_D(BRNVAR, t)
            B_ARR_J(BRNVAR, t)
            if (mc_branches) {
                B_VAR_M(BRNREF, t)
                B_VAR_G(BRNREF, t)
                B_VEC_M(BRNVAR, t)
                B_ARR_G(BRNVAR, t) }
        };

        void clear() {
            B_VEC_D(CLEAR)
            B_ARR_J(CLEAR)
            if (mc_branches) {
                B_VEC_M(CLEAR)
                B_ARR_G(CLEAR) }
        };

        void copy(photontree* t) {
            B_VAR_D(VARCOPY, t)
            B_VEC_D(VECCOPY, t)
            if (mc_branches) {
                B_VAR_M(VARCOPY, t)
                B_VEC_M(VECCOPY, t) }
        };

        void copy(jettree* t) {
            B_VAR_J(VARCOPY, t)
            B_ARR_J(ARRCOPY, t, nref)
            if (mc_branches) {
                B_VAR_G(VARCOPY, t)
                B_ARR_G(ARRCOPY, t, nref) }
        };

    private:
        bool mc_branches;

        B_VAR_D(DECLARE)
        B_VAR_M(DECLARE)
        B_VEC_D(DECLPTR)
        B_VEC_M(DECLPTR)
        B_VAR_J(DECLARE)
        B_VAR_G(DECLARE)
        B_ARR_J(DECLPTR)
        B_ARR_G(DECLPTR)
};

#endif /* PJTREE_H */
