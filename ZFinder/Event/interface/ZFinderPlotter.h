#ifndef ZFINDER_ZFINDERPLOTTER_H_
#define ZFINDER_ZFINDERPLOTTER_H_

// Root
#include <TH1D.h>  // TH1D
#include <TH2.h>  // TH2D

// CMSSW
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ZFinder Code
#include "ZFinderEvent.h"  // ZFinderEvent


namespace zf {
    class ZFinderPlotter{
        public:
            // Constructor
            ZFinderPlotter(TFileDirectory& tdir, const bool USE_MC = false);

            // Add events
            void Fill(
                    const ZFinderEvent& zf_event,
                    const int first_electron = 0,
                    const int second_electron = 1,
                    const int first_muon = 0,
                    const int second_muon = 1,
                    const double EVENT_WEIGHT = 1.
                    );
            // Make PNGs
            void Print(const std::string& basename);

        protected:
            // Histograms
            TH1D* z0_mass_all_;
            TH1D* z0_mass_coarse_;
            TH1D* z0_mass_fine_;
            TH1D* z0_rapidity_;
            TH1D* z0_pt_;
            TH1D* e0_pt_;
            TH1D* e1_pt_;
            TH1D* e0_eta_;
            TH1D* e1_eta_;
            TH1D* e0_phi_;
            TH1D* e1_phi_;
            TH1D* e0_charge_;
            TH1D* e1_charge_;
            TH1D* phistar_;
            TH1D* pileup_;
            TH1D* nelectrons_;

            TH1D* jpsi0_mass_all_;
            TH1D* jpsi0_mass_coarse_;
            TH1D* jpsi0_mass_fine_;
            TH1D* jpsi0_rapidity_;
            TH1D* jpsi0_pt_;
            TH1D* jpsi0_distance_;
            TH1D* jpsi0_dist_err_;
            TH1D* jpsi0_chi2_;
            TH1D* jpsi0_tau_;
            TH1D* jpsi0_zpt_difference_;
            TH1D* mu0_pt_;
            TH1D* mu1_pt_;
            TH1D* mu0_eta_;
            TH1D* mu1_eta_;
            TH1D* mu0_phi_;
            TH1D* mu1_phi_;
            TH1D* mu0_charge_;
            TH1D* mu1_charge_;
            TH1D* nmuons_;
            TH1D* njpsis_;
            TH1D* jpsi_vtx_distance_z_vtx_x_;
            TH1D* jpsi_vtx_distance_z_vtx_y_;
            TH1D* jpsi_vtx_distance_z_vtx_z_;

            //TODO make this its own folder
            TH1D* jpsi_iso_mu0_;
            TH1D* jpsi_iso_sum_charged_hadron_pt_mu0_;
            TH1D* jpsi_iso_sum_charged_particle_pt_mu0_;
            TH1D* jpsi_iso_sum_neutral_hadron_et_mu0_;
            TH1D* jpsi_iso_sum_photon_et_mu0_;
            TH1D* jpsi_iso_sum_pileup_pt_mu0_;

            TH1D* jpsi_iso_mu1_;
            TH1D* jpsi_iso_sum_charged_hadron_pt_mu1_;
            TH1D* jpsi_iso_sum_charged_particle_pt_mu1_;
            TH1D* jpsi_iso_sum_neutral_hadron_et_mu1_;
            TH1D* jpsi_iso_sum_photon_et_mu1_;
            TH1D* jpsi_iso_sum_pileup_pt_mu1_;

            TH2* jpsi0_mass_vs_chi2_;

            // Use the MC or reco data
            const bool USE_MC_;

            // Plotting variables
            static const int X_SIZE = 1280;
            static const int Y_SIZE = 640;
    };
}  // namespace zf
#endif  // ZFINDER_ZFINDERPLOTTER_H_
