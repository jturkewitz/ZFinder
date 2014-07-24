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
            ZFinderPlotter(TFileDirectory& tdir, const bool USE_MC = false, const bool APPLY_MUON_MIN_PT = false, const bool APPLY_JPSI_MASS_WINDOW = false, const bool APPLY_VERTEX_Z_POS_WINDOW = false);

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
            // TODO why z0? there is only ever 1 Z in code
            TH1D* z0_mass_all_;
            TH1D* z0_mass_coarse_;
            TH1D* z0_mass_fine_;
            TH1D* z0_rapidity_;
            TH1D* z0_pt_;
            TH1D* z_vtx_prob_;

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
            TH2D* jpsi0_pt_vs_rap_;
            TH1D* jpsi0_distance_;
            TH1D* jpsi0_dist_err_;
            TH1D* jpsi0_chi2_;
            TH1D* jpsi0_distance_xy_;
            TH1D* jpsi0_dist_err_xy_;
            TH1D* jpsi0_chi2_xy_;
            TH1D* jpsi0_tau_xy_;
            TH1D* jpsi0_tau_xy_fine_;
            TH1D* jpsi0_tau_xy_very_fine_;

            TH1D* jpsi0_tau_xy_very_fine_ptUnder10_;
            TH1D* jpsi0_tau_xy_very_fine_pt10to15_;
            TH1D* jpsi0_tau_xy_very_fine_pt15to20_;
            TH1D* jpsi0_tau_xy_very_fine_pt20to25_;
            TH1D* jpsi0_tau_xy_very_fine_pt25to30_;
            TH1D* jpsi0_tau_xy_very_fine_ptAbove20_;
            TH1D* jpsi0_tau_xy_very_fine_ptAbove30_;

            TH1D* jpsi0_tau_xy_very_fine_rap0_0to0_3_;
            TH1D* jpsi0_tau_xy_very_fine_rap0_3to0_6_;
            TH1D* jpsi0_tau_xy_very_fine_rap0_6to0_9_;
            TH1D* jpsi0_tau_xy_very_fine_rap0_9to1_2_;
            TH1D* jpsi0_tau_xy_very_fine_rap1_2to1_5_;
            TH1D* jpsi0_tau_xy_very_fine_rap1_5to1_8_;
            TH1D* jpsi0_tau_xy_very_fine_rap1_8to2_1_;
            TH1D* jpsi0_tau_xy_very_fine_rap2_1to2_4_;

            TH1D* jpsi0_tau_xy_very_fine_rap0_0to0_9pt10to15_;
            TH1D* jpsi0_tau_xy_very_fine_rap0_9to1_2pt10to15_;
            TH1D* jpsi0_tau_xy_very_fine_rap1_2to2_1pt10to15_;

            TH1D* jpsi0_tau_z_;
            TH1D* jpsi0_tau_z_fine_;
            TH1D* jpsi0_tau_z_very_fine_;
            TH1D* jpsi0_zpt_difference_;

            TH1D* dimuon_vtx_prob_;
            TH1D* dimuon_delta_phi_;
            TH1D* dimuon_delta_eta_;
            TH1D* dimuon_deltaR_;

            TH1D* z_jpsi_delta_phi_;

            TH1D* mu0_pt_;
            TH1D* mu1_pt_;
            TH1D* mu0_eta_;
            TH1D* mu1_eta_;
            TH1D* mu0_phi_;
            TH1D* mu1_phi_;
            TH1D* mu0_charge_;
            TH1D* mu1_charge_;

            TH1D* jet_pt_;
            TH1D* jet_eta_;
            TH1D* jet_btag_discriminator_;
            TH1D* muon_jet_pt_;
            TH1D* muon_jet_pt_diff_z_pt_;
            TH1D* muon_jet_pt_diff_dimuon_pt_;
            TH1D* muon_jet_eta_;
            TH1D* muon_jet_btag_discriminator_;

            TH2D* muon_jet_pt_z_pt_;
            TH2D* muon_jet_pt_dimuon_pt_;
            TH2D* muon_jet_phi_z_phi_;
            TH2D* muon_jet_phi_dimuon_phi_;

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

            TH2D* jpsi0_mass_vs_chi2_;
            TH2D* jpsi0_tau_xy_vs_tau_z_;
            TH2D* jpsi0_tau_xy_vs_distance_z_;
            TH2D* jpsi0_tau_z_vs_distance_z_;

            TH1D* vtx_x_;
            TH1D* vtx_y_;
            TH1D* vtx_z_;

            TH1D* primary_vtx_x_;
            TH1D* primary_vtx_y_;
            TH1D* primary_vtx_z_;

            TH1D* z_vtx_x_;
            TH1D* z_vtx_y_;
            TH1D* z_vtx_z_;

            TH1D* dimuon_vtx_x_;
            TH1D* dimuon_vtx_y_;
            TH1D* dimuon_vtx_z_;

            TH2D* primary_vtx_x_vs_z_vtx_x_;
            TH2D* primary_vtx_y_vs_z_vtx_y_;
            TH2D* primary_vtx_z_vs_z_vtx_z_;

            TH1D* nmuons_;
            TH1D* njets_;
            TH1D* n_muonjets_;
            TH1D* njpsis_;

            // Use the MC or reco data
            const bool USE_MC_;
            const bool APPLY_MUON_MIN_PT_;
            const bool APPLY_JPSI_MASS_WINDOW_;
            const bool APPLY_VERTEX_Z_POS_WINDOW_;
            // Plotting variables
            static const int X_SIZE = 1280;
            static const int Y_SIZE = 640;
    };
}  // namespace zf
#endif  // ZFINDER_ZFINDERPLOTTER_H_
