#ifndef PJTREE_H
#define PJTREE_H

#include "TTree.h"

#include <vector>

#include "defines.h"

#define BRANCHES(ACTION, ...)                                               \
    B_VAR_D(ACTION, ## __VA_ARGS__)                                         \
    B_VAR_M(ACTION, ## __VA_ARGS__)                                         \
    B_VEC_D(ACTION, ## __VA_ARGS__)                                         \
    B_VEC_M(ACTION, ## __VA_ARGS__)                                         \

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

class photontree;

class pjtree {
    public:
        pjtree() {
            this->mc_branches = 0;
            B_VAR_D(INVALID)
            B_VAR_M(INVALID)
        };

        pjtree(TTree* t, bool mc_branches)
                : pjtree() {
            this->mc_branches = mc_branches;
            branch(t);
        };

        ~pjtree() = default;

        void branch(TTree* t) {
            B_VAR_D(CREATE, t)
            B_VEC_D(CREATE, t)
            if (mc_branches) {
                B_VAR_M(CREATE, t)
                B_VEC_M(CREATE, t) }
        };

        void clear() {
            B_VEC_D(CLEAR)
            B_VEC_M(CLEAR)
        };

        void copy(photontree* t) {
            B_VAR_D(VARCOPY, t)
            B_VEC_D(VECCOPY, t)
            if (mc_branches) {
                B_VAR_M(VARCOPY, t)
                B_VEC_M(VECCOPY, t) }
        };

    private:
        bool mc_branches;

        BRANCHES(DECLARE)
};

#endif /* PJTREE_H */
