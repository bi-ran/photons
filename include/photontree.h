#ifndef PHOTONTREE_H
#define PHOTONTREE_H

#include "TTree.h"

#include <vector>

#include "defines.h"

#define B_AVE_D(ACTION, ...)                                                \
    ACTION(UInt_t, run, ## __VA_ARGS__)                                     \
    ACTION(ULong64_t, event, ## __VA_ARGS__)                                \
    ACTION(UInt_t, lumis, ## __VA_ARGS__)                                   \
                                                                            \
    ACTION(Int_t, nPho, ## __VA_ARGS__)                                     \

#define B_ARE_D(ACTION, ...)                                                \
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

#define B_AVE_M(ACTION, ...)                                                \
    ACTION(Int_t, nMC, ## __VA_ARGS__)                                      \

#define B_ARE_M(ACTION, ...)                                                \
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

class pjtree;

class photontree {
    friend pjtree;

    public:
        photontree() {
            this->mc_branches = 0;
            B_AVE_D(ZERO)
            B_AVE_M(ZERO)
            B_ARE_D(ZERO)
            B_ARE_M(ZERO)
        };

        photontree(TTree* t, bool mc_branches)
                : photontree() {
            this->mc_branches = mc_branches;
            read(t);
        }

        ~photontree() = default;

        void read(TTree* t) {
            B_AVE_D(RREF, t)
            B_ARE_D(RVAR, t)
            if (mc_branches) {
                B_AVE_M(RREF, t)
                B_ARE_M(RVAR, t) }
        };

    private:
        bool mc_branches;

        B_AVE_D(DECLARE)
        B_AVE_M(DECLARE)
        B_ARE_D(DECLPTR)
        B_ARE_M(DECLPTR)
};

#endif  /* PHOTONTREE_H */
