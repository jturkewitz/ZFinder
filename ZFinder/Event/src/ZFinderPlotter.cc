#include "ZFinder/Event/interface/ZFinderPlotter.h"

// Standard Library
#include <string>  // string

// Root
#include <TCanvas.h>  // TCanvas

// ZFinder Code
#include "ZFinder/Event/interface/ZFinderElectron.h"  // ZFinderElectron


namespace zf {
    // Constructor
    ZFinderPlotter::ZFinderPlotter(TFileDirectory& tdir, const bool USE_MC) : USE_MC_(USE_MC) {
        /*
         * Initialize a set of histograms and associate them with a given TDirectory.
         */

        // Set up histograms
        // z0_mass_all_
        const std::string z0_mass_all_name = "Z0 Mass: All";
        z0_mass_all_ = tdir.make<TH1D>(z0_mass_all_name.c_str(), z0_mass_all_name.c_str(), 300, 0., 300.);
        z0_mass_all_->GetXaxis()->SetTitle("m_{ee} [GeV]");
        z0_mass_all_->GetYaxis()->SetTitle("Counts / GeV");

        // z0_mass_coarse_
        const std::string z0_mass_coarse_name = "Z0 Mass: Coarse";
        z0_mass_coarse_ = tdir.make<TH1D>(z0_mass_coarse_name.c_str(), z0_mass_coarse_name.c_str(), 100, 50., 150.);
        z0_mass_coarse_->GetXaxis()->SetTitle("m_{ee} [GeV]");
        z0_mass_coarse_->GetYaxis()->SetTitle("Counts / GeV");

        // z0_mass_fine_
        const std::string z0_mass_fine_name = "Z0 Mass: Fine";
        z0_mass_fine_ = tdir.make<TH1D>(z0_mass_fine_name.c_str(), z0_mass_fine_name.c_str(), 80, 80., 100.);
        z0_mass_fine_->GetXaxis()->SetTitle("m_{ee} [GeV]");
        z0_mass_fine_->GetYaxis()->SetTitle("Counts / 0.25 GeV");

        // z0_rapidity_
        const std::string z0_rapidity_name = "Z0 Rapidity";
        z0_rapidity_ = tdir.make<TH1D>(z0_rapidity_name.c_str(), z0_rapidity_name.c_str(), 100, -5., 5.);
        z0_rapidity_->GetXaxis()->SetTitle("Z_{Y}");
        z0_rapidity_->GetYaxis()->SetTitle("Counts");

        // z0_pt
        const std::string z0_pt_name = "Z0 p_{T}";
        z0_pt_ = tdir.make<TH1D>(z0_pt_name.c_str(), z0_pt_name.c_str(), 200, 0., 200.);
        z0_pt_->GetXaxis()->SetTitle("p_{T,Z}");
        z0_pt_->GetYaxis()->SetTitle("Counts / GeV");

        // e0_pt
        const std::string e0_pt_name = "p_{T,e_{0}}";
        e0_pt_ = tdir.make<TH1D>(e0_pt_name.c_str(), e0_pt_name.c_str(), 200, 0., 200.);
        e0_pt_->GetXaxis()->SetTitle("p_{T,e_{0}}");
        e0_pt_->GetYaxis()->SetTitle("Counts / GeV");

        // e1_pt
        const std::string e1_pt_name = "p_{T,e_{1}}";
        e1_pt_ = tdir.make<TH1D>(e1_pt_name.c_str(), e1_pt_name.c_str(), 200, 0., 200.);
        e1_pt_->GetXaxis()->SetTitle("p_{T,e_{0}}");
        e1_pt_->GetYaxis()->SetTitle("Counts / GeV");

        // e0_eta_
        const std::string e0_eta_name = "#eta_{e_{0}}";
        e0_eta_ = tdir.make<TH1D>(e0_eta_name.c_str(), e0_eta_name.c_str(), 50, -5., 5.);
        e0_eta_->GetXaxis()->SetTitle("#eta_{e_{0}}");
        e0_eta_->GetYaxis()->SetTitle("Counts");

        // e1_eta_
        const std::string e1_eta_name = "#eta_{e_{1}}";
        e1_eta_ = tdir.make<TH1D>(e1_eta_name.c_str(), e1_eta_name.c_str(), 50, -5., 5.);
        e1_eta_->GetXaxis()->SetTitle("#eta_{e_{1}}");
        e1_eta_->GetYaxis()->SetTitle("Counts");

        // e0_phi_
        const std::string e0_phi_name = "#phi_{e_{0}}";
        e0_phi_ = tdir.make<TH1D>(e0_phi_name.c_str(), e0_phi_name.c_str(), 60, -3.15, 3.15);
        e0_phi_->GetXaxis()->SetTitle("#phi_{e_{0}}");
        e0_phi_->GetYaxis()->SetTitle("Counts");

        // e1_phi_
        const std::string e1_phi_name = "#phi_{e_{1}}";
        e1_phi_ = tdir.make<TH1D>(e1_phi_name.c_str(), e1_phi_name.c_str(), 50, -3.15, 3.15);
        e1_phi_->GetXaxis()->SetTitle("#phi_{e_{1}}");
        e1_phi_->GetYaxis()->SetTitle("counts");

        // e0_charge_
        const std::string e0_charge_name = "charge_{e_{0}}";
        e0_charge_ = tdir.make<TH1D>(e0_charge_name.c_str(), e0_charge_name.c_str(), 60, -3.15, 3.15);
        e0_charge_->GetXaxis()->SetTitle("charge_{e_{0}}");
        e0_charge_->GetYaxis()->SetTitle("Counts");

        // e1_charge_
        const std::string e1_charge_name = "charge_{e_{1}}";
        e1_charge_ = tdir.make<TH1D>(e1_charge_name.c_str(), e1_charge_name.c_str(), 60, -3.15, 3.15);
        e1_charge_->GetXaxis()->SetTitle("charge_{e_{1}}");
        e1_charge_->GetYaxis()->SetTitle("counts");

        // phistar
        const std::string phistar_name = "#phi*";
        phistar_ = tdir.make<TH1D>(phistar_name.c_str(), phistar_name.c_str(), 100, 0., 1.);
        phistar_->GetXaxis()->SetTitle("#phi*");
        phistar_->GetYaxis()->SetTitle("Counts");

        // jpsi0_mass_all_
        const std::string jpsi0_mass_all_name = "jpsi0 Mass: All";
        jpsi0_mass_all_ = tdir.make<TH1D>(jpsi0_mass_all_name.c_str(), jpsi0_mass_all_name.c_str(), 300, 0., 300.);
        jpsi0_mass_all_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
        jpsi0_mass_all_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi0_mass_coarse_
        const std::string jpsi0_mass_coarse_name = "jpsi0 Mass: Coarse";
        jpsi0_mass_coarse_ = tdir.make<TH1D>(jpsi0_mass_coarse_name.c_str(), jpsi0_mass_coarse_name.c_str(), 100, 0.0, 10.0);
        jpsi0_mass_coarse_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
        jpsi0_mass_coarse_->GetYaxis()->SetTitle("Counts / 0.1 GeV");

        // jpsi0_mass_fine_
        const std::string jpsi0_mass_fine_name = "jpsi0 Mass: Fine";
        jpsi0_mass_fine_ = tdir.make<TH1D>(jpsi0_mass_fine_name.c_str(), jpsi0_mass_fine_name.c_str(), 100, 2.5, 3.5);
        jpsi0_mass_fine_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
        jpsi0_mass_fine_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

        // jpsi0_rapidity_
        const std::string jpsi0_rapidity_name = "jpsi0 Rapidity";
        jpsi0_rapidity_ = tdir.make<TH1D>(jpsi0_rapidity_name.c_str(), jpsi0_rapidity_name.c_str(), 100, -5., 5.);
        jpsi0_rapidity_->GetXaxis()->SetTitle("JPsi_{Y}");
        jpsi0_rapidity_->GetYaxis()->SetTitle("Counts");

        // jpsi0_pt
        const std::string jpsi0_pt_name = "jpsi0 p_{T}";
        jpsi0_pt_ = tdir.make<TH1D>(jpsi0_pt_name.c_str(), jpsi0_pt_name.c_str(), 200, 0., 200.);
        jpsi0_pt_->GetXaxis()->SetTitle("p_{T,JPsi}");
        jpsi0_pt_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_vtx_distance_z_vtx_x
        const std::string jpsi_vtx_distance_z_vtx_x_name = "jpsi_vtx_x - z_vtx_x";
        jpsi_vtx_distance_z_vtx_x_ = tdir.make<TH1D>(jpsi_vtx_distance_z_vtx_x_name.c_str(), jpsi_vtx_distance_z_vtx_x_name.c_str(), 1000, -5., 5.);
        jpsi_vtx_distance_z_vtx_x_->GetXaxis()->SetTitle("distance [cm]");
        jpsi_vtx_distance_z_vtx_x_->GetYaxis()->SetTitle("Counts / 0.01 cm");

        // jpsi_vtx_distance_z_vtx_y
        const std::string jpsi_vtx_distance_z_vtx_y_name = "jpsi_vtx_y - z_vtx_y";
        jpsi_vtx_distance_z_vtx_y_ = tdir.make<TH1D>(jpsi_vtx_distance_z_vtx_y_name.c_str(), jpsi_vtx_distance_z_vtx_y_name.c_str(), 1000, -5., 5.);
        jpsi_vtx_distance_z_vtx_y_->GetXaxis()->SetTitle("distance [cm]");
        jpsi_vtx_distance_z_vtx_y_->GetYaxis()->SetTitle("Counts / 0.01 cm");

        // jpsi_vtx_distance_z_vtx_z
        const std::string jpsi_vtx_distance_z_vtx_z_name = "jpsi_vtx_z - z_vtx_z";
        jpsi_vtx_distance_z_vtx_z_ = tdir.make<TH1D>(jpsi_vtx_distance_z_vtx_z_name.c_str(), jpsi_vtx_distance_z_vtx_z_name.c_str(), 1000, -5., 5.);
        jpsi_vtx_distance_z_vtx_z_->GetXaxis()->SetTitle("distance [cm]");
        jpsi_vtx_distance_z_vtx_z_->GetYaxis()->SetTitle("Counts / 0.01 cm");

        // jpsi0_distance
        const std::string jpsi0_distance_name = "jpsi0 distance";
        jpsi0_distance_ = tdir.make<TH1D>(jpsi0_distance_name.c_str(), jpsi0_distance_name.c_str(), 500, 0., 5.);
        jpsi0_distance_->GetXaxis()->SetTitle("distance [cm]");
        jpsi0_distance_->GetYaxis()->SetTitle("Counts / 0.01 cm");

        // jpsi0_dist_err
        const std::string jpsi0_dist_err_name = "jpsi0 dist_err";
        jpsi0_dist_err_ = tdir.make<TH1D>(jpsi0_dist_err_name.c_str(), jpsi0_dist_err_name.c_str(), 500, 0., 5.);
        jpsi0_dist_err_->GetXaxis()->SetTitle("dist_err [cm]");
        jpsi0_dist_err_->GetYaxis()->SetTitle("Counts / 0.01 cm");

        // jpsi0_chi2
        const std::string jpsi0_chi2_name = "jpsi0 chi2";
        jpsi0_chi2_ = tdir.make<TH1D>(jpsi0_chi2_name.c_str(), jpsi0_chi2_name.c_str(), 1000, 0., 100.);
        jpsi0_chi2_->GetXaxis()->SetTitle("chi2");
        jpsi0_chi2_->GetYaxis()->SetTitle("Counts / 0.1 ");

        // jpsi0_zpt_difference_
        const std::string jpsi0_zpt_difference__name = "zPt - jpsipT";
        jpsi0_zpt_difference_ = tdir.make<TH1D>(jpsi0_zpt_difference__name.c_str(), jpsi0_zpt_difference__name.c_str(), 200, -100., 100.);
        jpsi0_zpt_difference_->GetXaxis()->SetTitle("pT [GeV]");
        jpsi0_zpt_difference_->GetYaxis()->SetTitle("Counts / 1 GeV ");

        // jpsi0_tau
        const std::string jpsi0_tau_name = "jpsi0 tau";
        jpsi0_tau_ = tdir.make<TH1D>(jpsi0_tau_name.c_str(), jpsi0_tau_name.c_str(), 200, -100., 100.);
        jpsi0_tau_->GetXaxis()->SetTitle("tau [ps]");
        jpsi0_tau_->GetYaxis()->SetTitle("Counts / 1 ps ");

        // jpsi0_mass_vs_chi2
        const std::string jpsi0_mass_vs_chi2_name = "jpsi0 mass vs chi2";
        jpsi0_mass_vs_chi2_ = tdir.make<TH2D>(jpsi0_mass_vs_chi2_name.c_str(), jpsi0_mass_vs_chi2_name.c_str(), 100, 0., 10., 100, 0., 10.);
        jpsi0_mass_vs_chi2_->GetXaxis()->SetTitle("chi2");
        jpsi0_mass_vs_chi2_->GetYaxis()->SetTitle("mass");

        //TODO remove 0 from jpsi0
        // jpsi_iso_mu0
        const std::string jpsi_iso_mu0_name = "(sum charged had pt + neutral Et + photon Et ) / (mu0 pt)";
        jpsi_iso_mu0_ = tdir.make<TH1D>(jpsi_iso_mu0_name.c_str(), jpsi_iso_mu0_name.c_str(), 50, 0., 5.);
        jpsi_iso_mu0_->GetXaxis()->SetTitle("isolation (sum had + photon pt) / (mu0 pt) ");
        jpsi_iso_mu0_->GetYaxis()->SetTitle("Counts / 0.1");

        // jpsi_iso_sum_charged_hadron_pt_mu0
        const std::string jpsi_iso_sum_charged_hadron_pt_mu0_name = "charged_hadron_pt_mu0 p_{T}";
        jpsi_iso_sum_charged_hadron_pt_mu0_ = tdir.make<TH1D>(jpsi_iso_sum_charged_hadron_pt_mu0_name.c_str(), jpsi_iso_sum_charged_hadron_pt_mu0_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_charged_hadron_pt_mu0_->GetXaxis()->SetTitle("p_{T} charged hadrons");
        jpsi_iso_sum_charged_hadron_pt_mu0_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_sum_charged_particle_pt_mu0
        const std::string jpsi_iso_sum_charged_particle_pt_mu0_name = "charged_particle_pt_mu0 p_{T}";
        jpsi_iso_sum_charged_particle_pt_mu0_ = tdir.make<TH1D>(jpsi_iso_sum_charged_particle_pt_mu0_name.c_str(), jpsi_iso_sum_charged_particle_pt_mu0_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_charged_particle_pt_mu0_->GetXaxis()->SetTitle("p_{T} charged particles");
        jpsi_iso_sum_charged_particle_pt_mu0_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_sum_neutral_hadron_et_mu0
        const std::string jpsi_iso_sum_neutral_hadron_et_mu0_name = "neutral_hadron_et_mu0 E_{T}";
        jpsi_iso_sum_neutral_hadron_et_mu0_ = tdir.make<TH1D>(jpsi_iso_sum_neutral_hadron_et_mu0_name.c_str(), jpsi_iso_sum_neutral_hadron_et_mu0_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_neutral_hadron_et_mu0_->GetXaxis()->SetTitle("E_{T} neutral hadrons");
        jpsi_iso_sum_neutral_hadron_et_mu0_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_sum_photon_et_mu0
        const std::string jpsi_iso_sum_photon_et_mu0_name = "photon_et_mu0 E_{T}";
        jpsi_iso_sum_photon_et_mu0_ = tdir.make<TH1D>(jpsi_iso_sum_photon_et_mu0_name.c_str(), jpsi_iso_sum_photon_et_mu0_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_photon_et_mu0_->GetXaxis()->SetTitle("E_{T} photons");
        jpsi_iso_sum_photon_et_mu0_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_sum_pileup_pt_mu0
        const std::string jpsi_iso_sum_pileup_pt_mu0_name = "pileup_pt_mu0 p_{T}";
        jpsi_iso_sum_pileup_pt_mu0_ = tdir.make<TH1D>(jpsi_iso_sum_pileup_pt_mu0_name.c_str(), jpsi_iso_sum_pileup_pt_mu0_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_pileup_pt_mu0_->GetXaxis()->SetTitle("p_{T} pileup");
        jpsi_iso_sum_pileup_pt_mu0_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_mu1
        const std::string jpsi_iso_mu1_name = "(sum charged had pt + neutral Et + photon Et ) / (mu1 pt)";
        jpsi_iso_mu1_ = tdir.make<TH1D>(jpsi_iso_mu1_name.c_str(), jpsi_iso_mu1_name.c_str(), 50, 0., 5.);
        jpsi_iso_mu1_->GetXaxis()->SetTitle("isolation (sum had + photon pt) / (mu1 pt) ");
        jpsi_iso_mu1_->GetYaxis()->SetTitle("Counts / 0.1");

        // jpsi_iso_sum_charged_hadron_pt_mu1
        const std::string jpsi_iso_sum_charged_hadron_pt_mu1_name = "charged_hadron_pt_mu1 p_{T}";
        jpsi_iso_sum_charged_hadron_pt_mu1_ = tdir.make<TH1D>(jpsi_iso_sum_charged_hadron_pt_mu1_name.c_str(), jpsi_iso_sum_charged_hadron_pt_mu1_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_charged_hadron_pt_mu1_->GetXaxis()->SetTitle("p_{T} charged hadrons");
        jpsi_iso_sum_charged_hadron_pt_mu1_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_sum_charged_particle_pt_mu1
        const std::string jpsi_iso_sum_charged_particle_pt_mu1_name = "charged_particle_pt_mu1 p_{T}";
        jpsi_iso_sum_charged_particle_pt_mu1_ = tdir.make<TH1D>(jpsi_iso_sum_charged_particle_pt_mu1_name.c_str(), jpsi_iso_sum_charged_particle_pt_mu1_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_charged_particle_pt_mu1_->GetXaxis()->SetTitle("p_{T} charged particles");
        jpsi_iso_sum_charged_particle_pt_mu1_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_sum_neutral_hadron_et_mu1
        const std::string jpsi_iso_sum_neutral_hadron_et_mu1_name = "neutral_hadron_et_mu1 E_{T}";
        jpsi_iso_sum_neutral_hadron_et_mu1_ = tdir.make<TH1D>(jpsi_iso_sum_neutral_hadron_et_mu1_name.c_str(), jpsi_iso_sum_neutral_hadron_et_mu1_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_neutral_hadron_et_mu1_->GetXaxis()->SetTitle("E_{T} neutral hadrons");
        jpsi_iso_sum_neutral_hadron_et_mu1_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_sum_photon_et_mu1
        const std::string jpsi_iso_sum_photon_et_mu1_name = "photon_et_mu1 E_{T}";
        jpsi_iso_sum_photon_et_mu1_ = tdir.make<TH1D>(jpsi_iso_sum_photon_et_mu1_name.c_str(), jpsi_iso_sum_photon_et_mu1_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_photon_et_mu1_->GetXaxis()->SetTitle("E_{T} photons");
        jpsi_iso_sum_photon_et_mu1_->GetYaxis()->SetTitle("Counts / GeV");

        // jpsi_iso_sum_pileup_pt_mu1
        const std::string jpsi_iso_sum_pileup_pt_mu1_name = "pileup_pt_mu1 p_{T}";
        jpsi_iso_sum_pileup_pt_mu1_ = tdir.make<TH1D>(jpsi_iso_sum_pileup_pt_mu1_name.c_str(), jpsi_iso_sum_pileup_pt_mu1_name.c_str(), 200, 0., 200.);
        jpsi_iso_sum_pileup_pt_mu1_->GetXaxis()->SetTitle("p_{T} pileup");
        jpsi_iso_sum_pileup_pt_mu1_->GetYaxis()->SetTitle("Counts / GeV");

        // mu0_pt
        const std::string mu0_pt_name = "p_{T,mu_{0}}";
        mu0_pt_ = tdir.make<TH1D>(mu0_pt_name.c_str(), mu0_pt_name.c_str(), 200, 0., 200.);
        mu0_pt_->GetXaxis()->SetTitle("p_{T,mu_{0}}");
        mu0_pt_->GetYaxis()->SetTitle("Counts / GeV");

        // mu1_pt
        const std::string mu1_pt_name = "p_{T,mu_{1}}";
        mu1_pt_ = tdir.make<TH1D>(mu1_pt_name.c_str(), mu1_pt_name.c_str(), 200, 0., 200.);
        mu1_pt_->GetXaxis()->SetTitle("p_{T,mu_{0}}");
        mu1_pt_->GetYaxis()->SetTitle("Counts / GeV");

        // mu0_eta_
        const std::string mu0_eta_name = "#eta_{mu_{0}}";
        mu0_eta_ = tdir.make<TH1D>(mu0_eta_name.c_str(), mu0_eta_name.c_str(), 50, -5., 5.);
        mu0_eta_->GetXaxis()->SetTitle("#eta_{mu_{0}}");
        mu0_eta_->GetYaxis()->SetTitle("Counts");

        // mu1_eta_
        const std::string mu1_eta_name = "#eta_{mu_{1}}";
        mu1_eta_ = tdir.make<TH1D>(mu1_eta_name.c_str(), mu1_eta_name.c_str(), 50, -5., 5.);
        mu1_eta_->GetXaxis()->SetTitle("#eta_{mu_{1}}");
        mu1_eta_->GetYaxis()->SetTitle("Counts");

        // mu0_phi_
        const std::string mu0_phi_name = "#phi_{mu_{0}}";
        mu0_phi_ = tdir.make<TH1D>(mu0_phi_name.c_str(), mu0_phi_name.c_str(), 60, -3.15, 3.15);
        mu0_phi_->GetXaxis()->SetTitle("#phi_{mu_{0}}");
        mu0_phi_->GetYaxis()->SetTitle("Counts");

        // mu1_phi_
        const std::string mu1_phi_name = "#phi_{mu_{1}}";
        mu1_phi_ = tdir.make<TH1D>(mu1_phi_name.c_str(), mu1_phi_name.c_str(), 50, -3.15, 3.15);
        mu1_phi_->GetXaxis()->SetTitle("#phi_{mu_{1}}");
        mu1_phi_->GetYaxis()->SetTitle("counts");

        // mu0_charge_
        const std::string mu0_charge_name = "charge_{mu_{0}}";
        mu0_charge_ = tdir.make<TH1D>(mu0_charge_name.c_str(), mu0_charge_name.c_str(), 60, -3.15, 3.15);
        mu0_charge_->GetXaxis()->SetTitle("charge_{mu_{0}}");
        mu0_charge_->GetYaxis()->SetTitle("Counts");

        // mu1_charge_
        const std::string mu1_charge_name = "charge_{mu_{1}}";
        mu1_charge_ = tdir.make<TH1D>(mu1_charge_name.c_str(), mu1_charge_name.c_str(), 60, -3.15, 3.15);
        mu1_charge_->GetXaxis()->SetTitle("charge_{mu_{1}}");
        mu1_charge_->GetYaxis()->SetTitle("counts");

        // nmuons
        const std::string nmuons_name = "N_{mu}";
        nmuons_ = tdir.make<TH1D>(nmuons_name.c_str(), nmuons_name.c_str(), 10, 0., 10.);
        nmuons_->GetXaxis()->SetTitle("N_{mu}");
        nmuons_->GetYaxis()->SetTitle("Events");

        // njpsis
        const std::string njpsis_name = "N_{jpsi}";
        njpsis_ = tdir.make<TH1D>(njpsis_name.c_str(), njpsis_name.c_str(), 10, 0., 10.);
        njpsis_->GetXaxis()->SetTitle("N_{jpsi}");
        njpsis_->GetYaxis()->SetTitle("Events");

        // nelectrons
        const std::string nelectrons_name = "N_{e}";
        nelectrons_ = tdir.make<TH1D>(nelectrons_name.c_str(), nelectrons_name.c_str(), 10, 0., 10.);
        nelectrons_->GetXaxis()->SetTitle("N_{e}");
        nelectrons_->GetYaxis()->SetTitle("Events");

        // pileup
        const std::string pileup_name = "N_{Vertices}";
        pileup_ = tdir.make<TH1D>(pileup_name.c_str(), pileup_name.c_str(), 100, 0., 100.);
        pileup_->GetXaxis()->SetTitle("Number of Vertices");
        pileup_->GetYaxis()->SetTitle("Counts");
    }

    void ZFinderPlotter::Fill(
            const ZFinderEvent& zf_event,
            const int electron_0,
            const int electron_1,
            const int muon_0,
            const int muon_1,
            const double EVENT_WEIGHT
            ) {
        /*
         * Given a zf_event, fills all the histograms.
         *
         * electron_0 and electron_1 can be used to assign zf_event.eN to the given
         * number in the histogram. For example, assigning electron_0 = 1 will fill
         * the e0 histograms with data from zf_event.e1.
         */
        // Z Info
        if (!USE_MC_) {
            z0_mass_all_->Fill(zf_event.reco_z.m, EVENT_WEIGHT);
            z0_mass_coarse_->Fill(zf_event.reco_z.m, EVENT_WEIGHT);
            z0_mass_fine_->Fill(zf_event.reco_z.m, EVENT_WEIGHT);
            z0_rapidity_->Fill(zf_event.reco_z.y, EVENT_WEIGHT);
            z0_pt_->Fill(zf_event.reco_z.pt, EVENT_WEIGHT);
            phistar_->Fill(zf_event.reco_z.phistar, EVENT_WEIGHT);

            for (unsigned int i = 0; i < zf_event.reco_jpsi.m.size() ; ++i ) {
                jpsi0_mass_all_->Fill(zf_event.reco_jpsi.m.at(i), EVENT_WEIGHT);
                jpsi0_mass_coarse_->Fill(zf_event.reco_jpsi.m.at(i), EVENT_WEIGHT);
                jpsi0_mass_fine_->Fill(zf_event.reco_jpsi.m.at(i), EVENT_WEIGHT);
                jpsi0_rapidity_->Fill(zf_event.reco_jpsi.y.at(i), EVENT_WEIGHT);
                jpsi0_pt_->Fill(zf_event.reco_jpsi.pt.at(i), EVENT_WEIGHT);
                jpsi_vtx_distance_z_vtx_x_->Fill( zf_event.reco_jpsi.vtx_x.at(i) - zf_event.reco_z.vtx_x , EVENT_WEIGHT);
                jpsi_vtx_distance_z_vtx_y_->Fill( zf_event.reco_jpsi.vtx_y.at(i) - zf_event.reco_z.vtx_y , EVENT_WEIGHT);
                jpsi_vtx_distance_z_vtx_z_->Fill( zf_event.reco_jpsi.vtx_z.at(i) - zf_event.reco_z.vtx_z , EVENT_WEIGHT);
                jpsi0_distance_->Fill(zf_event.reco_jpsi.distance.at(i), EVENT_WEIGHT);
                jpsi0_dist_err_->Fill(zf_event.reco_jpsi.dist_err.at(i), EVENT_WEIGHT);
                jpsi0_chi2_->Fill(zf_event.reco_jpsi.chi2.at(i), EVENT_WEIGHT);
                jpsi0_zpt_difference_->Fill( zf_event.reco_z.pt - zf_event.reco_jpsi.pt.at(i), EVENT_WEIGHT);
                jpsi0_tau_->Fill(zf_event.reco_jpsi.tau.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
                jpsi0_mass_vs_chi2_->Fill(zf_event.reco_jpsi.m.at(i) , zf_event.reco_jpsi.chi2.at(i) );

                mu0_pt_->Fill(zf_event.mu0.at(i).pt(), EVENT_WEIGHT);
                mu0_eta_->Fill(zf_event.mu0.at(i).eta(), EVENT_WEIGHT);
                mu0_phi_->Fill(zf_event.mu0.at(i).phi(), EVENT_WEIGHT);
                mu0_charge_->Fill(zf_event.mu0.at(i).charge(), EVENT_WEIGHT);
                mu1_pt_->Fill(zf_event.mu1.at(i).pt(), EVENT_WEIGHT);
                mu1_eta_->Fill(zf_event.mu1.at(i).eta(), EVENT_WEIGHT);
                mu1_phi_->Fill(zf_event.mu1.at(i).phi(), EVENT_WEIGHT);
                mu1_charge_->Fill(zf_event.mu1.at(i).charge(), EVENT_WEIGHT);
            }
            if (zf_event.reco_jpsi.m.size() > 0) {
                jpsi_iso_mu0_->Fill(zf_event.reco_jpsi.iso_mu0 ) ;
                jpsi_iso_sum_charged_hadron_pt_mu0_->Fill(zf_event.reco_jpsi.iso_sum_charged_hadron_pt_mu0 ) ;
                jpsi_iso_sum_charged_particle_pt_mu0_->Fill(zf_event.reco_jpsi.iso_sum_charged_particle_pt_mu0 ) ;
                jpsi_iso_sum_neutral_hadron_et_mu0_->Fill(zf_event.reco_jpsi.iso_sum_neutral_hadron_et_mu0 ) ;
                jpsi_iso_sum_photon_et_mu0_->Fill(zf_event.reco_jpsi.iso_sum_photon_et_mu0 ) ;
                jpsi_iso_sum_pileup_pt_mu0_->Fill(zf_event.reco_jpsi.iso_sum_pileup_pt_mu0 ) ;

                jpsi_iso_mu1_->Fill(zf_event.reco_jpsi.iso_mu1 ) ;
                jpsi_iso_sum_charged_hadron_pt_mu1_->Fill(zf_event.reco_jpsi.iso_sum_charged_hadron_pt_mu1 ) ;
                jpsi_iso_sum_charged_particle_pt_mu1_->Fill(zf_event.reco_jpsi.iso_sum_charged_particle_pt_mu1 ) ;
                jpsi_iso_sum_neutral_hadron_et_mu1_->Fill(zf_event.reco_jpsi.iso_sum_neutral_hadron_et_mu1 ) ;
                jpsi_iso_sum_photon_et_mu1_->Fill(zf_event.reco_jpsi.iso_sum_photon_et_mu1 ) ;
                jpsi_iso_sum_pileup_pt_mu1_->Fill(zf_event.reco_jpsi.iso_sum_pileup_pt_mu1 ) ;
            }

            // Fill the histograms with the information from the approriate electron
            if (zf_event.e0 != NULL && zf_event.e1 != NULL){
                if (electron_0 == 0 && electron_1 == 1) {
                    e0_pt_->Fill(zf_event.e0->pt, EVENT_WEIGHT);
                    e0_eta_->Fill(zf_event.e0->eta, EVENT_WEIGHT);
                    e0_phi_->Fill(zf_event.e0->phi, EVENT_WEIGHT);
                    e0_charge_->Fill(zf_event.e0->charge, EVENT_WEIGHT);
                    e1_pt_->Fill(zf_event.e1->pt, EVENT_WEIGHT);
                    e1_eta_->Fill(zf_event.e1->eta, EVENT_WEIGHT);
                    e1_phi_->Fill(zf_event.e1->phi, EVENT_WEIGHT);
                    e1_charge_->Fill(zf_event.e1->charge, EVENT_WEIGHT);
                } else if (electron_0 == 1 && electron_1 == 0) {
                    e0_pt_->Fill(zf_event.e1->pt, EVENT_WEIGHT);
                    e0_eta_->Fill(zf_event.e1->eta, EVENT_WEIGHT);
                    e0_phi_->Fill(zf_event.e1->phi, EVENT_WEIGHT);
                    e0_charge_->Fill(zf_event.e1->charge, EVENT_WEIGHT);
                    e1_pt_->Fill(zf_event.e0->pt, EVENT_WEIGHT);
                    e1_eta_->Fill(zf_event.e0->eta, EVENT_WEIGHT);
                    e1_phi_->Fill(zf_event.e0->phi, EVENT_WEIGHT);
                    e1_charge_->Fill(zf_event.e0->charge, EVENT_WEIGHT);
                }
            }

            // Fill the histograms with the information from the approriate muon
            // TODO remove clutter from electron_0 vs electron_1 (never used and convoluted)
            //if (zf_event.mu0 != NULL && zf_event.mu1 != NULL){
            
            //TODO clean this up

            //if (zf_event.reco_jpsi.m.size() > 0){
            //    if (electron_0 == 0 && electron_1 == 1) {
            //        mu0_pt_->Fill(zf_event.mu0.pt(), EVENT_WEIGHT);
            //        mu0_eta_->Fill(zf_event.mu0.eta(), EVENT_WEIGHT);
            //        mu0_phi_->Fill(zf_event.mu0.phi(), EVENT_WEIGHT);
            //        mu0_charge_->Fill(zf_event.mu0.charge(), EVENT_WEIGHT);
            //        mu1_pt_->Fill(zf_event.mu1.pt(), EVENT_WEIGHT);
            //        mu1_eta_->Fill(zf_event.mu1.eta(), EVENT_WEIGHT);
            //        mu1_phi_->Fill(zf_event.mu1.phi(), EVENT_WEIGHT);
            //        mu1_charge_->Fill(zf_event.mu1.charge(), EVENT_WEIGHT);
            //    } else if (electron_0 == 1 && electron_1 == 0) {
            //        mu0_pt_->Fill(zf_event.mu1.pt(), EVENT_WEIGHT);
            //        mu0_eta_->Fill(zf_event.mu1.eta(), EVENT_WEIGHT);
            //        mu0_phi_->Fill(zf_event.mu1.phi(), EVENT_WEIGHT);
            //        mu0_charge_->Fill(zf_event.mu1.charge(), EVENT_WEIGHT);
            //        mu1_pt_->Fill(zf_event.mu0.pt(), EVENT_WEIGHT);
            //        mu1_eta_->Fill(zf_event.mu0.eta(), EVENT_WEIGHT);
            //        mu1_phi_->Fill(zf_event.mu0.phi(), EVENT_WEIGHT);
            //        mu1_charge_->Fill(zf_event.mu0.charge(), EVENT_WEIGHT);
            //    }
            //}

            // Event Info
            pileup_->Fill(zf_event.reco_vert.num, EVENT_WEIGHT);
            nelectrons_->Fill(zf_event.n_reco_electrons, EVENT_WEIGHT);
            nmuons_->Fill(zf_event.n_reco_muons, EVENT_WEIGHT);
            njpsis_->Fill(zf_event.reco_jpsi.m.size(), EVENT_WEIGHT); 
        } else if (USE_MC_ && !zf_event.is_real_data) {
            z0_mass_all_->Fill(zf_event.truth_z.m, EVENT_WEIGHT);
            z0_mass_coarse_->Fill(zf_event.truth_z.m, EVENT_WEIGHT);
            z0_mass_fine_->Fill(zf_event.truth_z.m, EVENT_WEIGHT);
            z0_rapidity_->Fill(zf_event.truth_z.y, EVENT_WEIGHT);
            z0_pt_->Fill(zf_event.truth_z.pt, EVENT_WEIGHT);
            phistar_->Fill(zf_event.truth_z.phistar, EVENT_WEIGHT);

            // Fill the histograms with the information from the approriate electron
            if (zf_event.e0_truth != NULL && zf_event.e1_truth != NULL){
                if (electron_0 == 0 && electron_1 == 1) {
                    e0_pt_->Fill(zf_event.e0_truth->pt, EVENT_WEIGHT);
                    e0_eta_->Fill(zf_event.e0_truth->eta, EVENT_WEIGHT);
                    e0_phi_->Fill(zf_event.e0_truth->phi, EVENT_WEIGHT);
                    e1_pt_->Fill(zf_event.e1_truth->pt, EVENT_WEIGHT);
                    e1_eta_->Fill(zf_event.e1_truth->eta, EVENT_WEIGHT);
                    e1_phi_->Fill(zf_event.e1_truth->phi, EVENT_WEIGHT);
                } else if (electron_0 == 1 && electron_1 == 0) {
                    e0_pt_->Fill(zf_event.e1_truth->pt, EVENT_WEIGHT);
                    e0_eta_->Fill(zf_event.e1_truth->eta, EVENT_WEIGHT);
                    e0_phi_->Fill(zf_event.e1_truth->phi, EVENT_WEIGHT);
                    e1_pt_->Fill(zf_event.e0_truth->pt, EVENT_WEIGHT);
                    e1_eta_->Fill(zf_event.e0_truth->eta, EVENT_WEIGHT);
                    e1_phi_->Fill(zf_event.e0_truth->phi, EVENT_WEIGHT);
                }
            }
            // Event Info
            pileup_->Fill(zf_event.truth_vert.num, EVENT_WEIGHT);
            nelectrons_->Fill(2, EVENT_WEIGHT);  // We only ever grab the two electrons from the Z
        }
    }

    void ZFinderPlotter::Print(const std::string& basename) {
        // Write all PNGs
        std::string z0_mass_all_Str = basename + "_z0_mass_all" ;
        TCanvas* z0_mass_all_C = new TCanvas(z0_mass_all_Str.c_str(), z0_mass_all_Str.c_str(), X_SIZE, Y_SIZE);
        z0_mass_all_->Draw();
        z0_mass_all_C->Print((z0_mass_all_Str+".png").c_str());

        std::string z0_mass_coarse_Str = basename + "_z0_mass_coarse" ;
        TCanvas* z0_mass_coarse_C = new TCanvas(z0_mass_coarse_Str.c_str(), z0_mass_coarse_Str.c_str(), X_SIZE, Y_SIZE);
        z0_mass_coarse_->Draw();
        z0_mass_coarse_C->Print((z0_mass_coarse_Str+".png").c_str());

        std::string z0_mass_fine_Str = basename + "_z0_mass_fine";
        TCanvas* z0_mass_fine_C = new TCanvas(z0_mass_fine_Str.c_str(), z0_mass_fine_Str.c_str(), X_SIZE, Y_SIZE);
        z0_mass_fine_->Draw();
        z0_mass_fine_C->Print((z0_mass_fine_Str+".png").c_str());

        std::string z0_rapidity_Str = basename + "_z0_rapidity";
        TCanvas* z0_rapidity_C = new TCanvas(z0_rapidity_Str.c_str(), z0_rapidity_Str.c_str(), X_SIZE, Y_SIZE);
        z0_rapidity_->Draw();
        z0_rapidity_C->Print((z0_rapidity_Str+".png").c_str());

        std::string z0_ptStr = basename + "_z0_pt";
        TCanvas* z0_ptC = new TCanvas(z0_ptStr.c_str(), z0_ptStr.c_str(), X_SIZE, Y_SIZE);
        z0_pt_->Draw();
        z0_ptC->Print((z0_ptStr+".png").c_str());

        std::string e0_ptStr = basename + "_e0_pt";
        TCanvas* e0_ptC = new TCanvas(e0_ptStr.c_str(), e0_ptStr.c_str(), X_SIZE, Y_SIZE);
        e0_pt_->Draw();
        e0_ptC->Print((e0_ptStr+".png").c_str());

        std::string e1_ptStr = basename + "_e1_pt";
        TCanvas* e1_ptC = new TCanvas(e1_ptStr.c_str(), e1_ptStr.c_str(), X_SIZE, Y_SIZE);
        e1_pt_->Draw();
        e1_ptC->Print((e1_ptStr+".png").c_str());

        std::string e0_eta_Str = basename + "_e0_eta";
        TCanvas* e0_eta_C = new TCanvas(e0_eta_Str.c_str(), e0_eta_Str.c_str(), X_SIZE, Y_SIZE);
        e0_eta_->Draw();
        e0_eta_C->Print((e0_eta_Str+".png").c_str());

        std::string e1_eta_Str = basename + "_e1_eta";
        TCanvas* e1_eta_C = new TCanvas(e1_eta_Str.c_str(), e1_eta_Str.c_str(), X_SIZE, Y_SIZE);
        e1_eta_->Draw();
        e1_eta_C->Print((e1_eta_Str+".png").c_str());

        std::string e0_phi_Str = basename + "_e0_phi";
        TCanvas* e0_phi_C = new TCanvas(e0_phi_Str.c_str(), e0_phi_Str.c_str(), X_SIZE, Y_SIZE);
        e0_phi_->Draw();
        e0_phi_C->Print((e0_phi_Str+".png").c_str());

        std::string e1_phi_Str = basename + "_e1_phi";
        TCanvas* e1_phi_C = new TCanvas(e1_phi_Str.c_str(), e1_phi_Str.c_str(), X_SIZE, Y_SIZE);
        e1_phi_->Draw();
        e1_phi_C->Print((e1_phi_Str+".png").c_str());

        std::string e0_charge_Str = basename + "_e0_charge";
        TCanvas* e0_charge_C = new TCanvas(e0_charge_Str.c_str(), e0_charge_Str.c_str(), X_SIZE, Y_SIZE);
        e0_charge_->Draw();
        e0_charge_C->Print((e0_charge_Str+".png").c_str());

        std::string e1_charge_Str = basename + "_e1_charge";
        TCanvas* e1_charge_C = new TCanvas(e1_charge_Str.c_str(), e1_charge_Str.c_str(), X_SIZE, Y_SIZE);
        e1_charge_->Draw();
        e1_charge_C->Print((e1_charge_Str+".png").c_str());

        std::string phistarStr = basename + "_phistar";
        TCanvas* phistarC = new TCanvas(phistarStr.c_str(), phistarStr.c_str(), X_SIZE, Y_SIZE);
        phistar_->Draw();
        phistarC->Print((phistarStr+".png").c_str());

        std::string nelectronsStr = basename + "_nelectrons";
        TCanvas* nelectronsC = new TCanvas(nelectronsStr.c_str(), nelectronsStr.c_str(), X_SIZE, Y_SIZE);
        nelectrons_->Draw();
        nelectronsC->Print((nelectronsStr+".png").c_str());

        std::string jpsi0_mass_all_Str = basename + "_jpsi0_mass_all" ;
        TCanvas* jpsi0_mass_all_C = new TCanvas(jpsi0_mass_all_Str.c_str(), jpsi0_mass_all_Str.c_str(), X_SIZE, Y_SIZE);
        jpsi0_mass_all_->Draw();
        jpsi0_mass_all_C->Print((jpsi0_mass_all_Str+".png").c_str());

        std::string jpsi0_mass_coarse_Str = basename + "_jpsi0_mass_coarse" ;
        TCanvas* jpsi0_mass_coarse_C = new TCanvas(jpsi0_mass_coarse_Str.c_str(), jpsi0_mass_coarse_Str.c_str(), X_SIZE, Y_SIZE);
        jpsi0_mass_coarse_->Draw();
        jpsi0_mass_coarse_C->Print((jpsi0_mass_coarse_Str+".png").c_str());

        std::string jpsi0_mass_fine_Str = basename + "_jpsi0_mass_fine";
        TCanvas* jpsi0_mass_fine_C = new TCanvas(jpsi0_mass_fine_Str.c_str(), jpsi0_mass_fine_Str.c_str(), X_SIZE, Y_SIZE);
        jpsi0_mass_fine_->Draw();
        jpsi0_mass_fine_C->Print((jpsi0_mass_fine_Str+".png").c_str());

        std::string jpsi0_rapidity_Str = basename + "_jpsi0_rapidity";
        TCanvas* jpsi0_rapidity_C = new TCanvas(jpsi0_rapidity_Str.c_str(), jpsi0_rapidity_Str.c_str(), X_SIZE, Y_SIZE);
        jpsi0_rapidity_->Draw();
        jpsi0_rapidity_C->Print((jpsi0_rapidity_Str+".png").c_str());

        std::string jpsi0_ptStr = basename + "_jpsi0_pt";
        TCanvas* jpsi0_ptC = new TCanvas(jpsi0_ptStr.c_str(), jpsi0_ptStr.c_str(), X_SIZE, Y_SIZE);
        jpsi0_pt_->Draw();
        jpsi0_ptC->Print((jpsi0_ptStr+".png").c_str());

        std::string mu0_ptStr = basename + "_mu0_pt";
        TCanvas* mu0_ptC = new TCanvas(mu0_ptStr.c_str(), mu0_ptStr.c_str(), X_SIZE, Y_SIZE);
        mu0_pt_->Draw();
        mu0_ptC->Print((mu0_ptStr+".png").c_str());

        std::string mu1_ptStr = basename + "_mu1_pt";
        TCanvas* mu1_ptC = new TCanvas(mu1_ptStr.c_str(), mu1_ptStr.c_str(), X_SIZE, Y_SIZE);
        mu1_pt_->Draw();
        mu1_ptC->Print((mu1_ptStr+".png").c_str());

        std::string mu0_eta_Str = basename + "_mu0_eta";
        TCanvas* mu0_eta_C = new TCanvas(mu0_eta_Str.c_str(), mu0_eta_Str.c_str(), X_SIZE, Y_SIZE);
        mu0_eta_->Draw();
        mu0_eta_C->Print((mu0_eta_Str+".png").c_str());

        std::string mu1_eta_Str = basename + "_mu1_eta";
        TCanvas* mu1_eta_C = new TCanvas(mu1_eta_Str.c_str(), mu1_eta_Str.c_str(), X_SIZE, Y_SIZE);
        mu1_eta_->Draw();
        mu1_eta_C->Print((mu1_eta_Str+".png").c_str());

        std::string mu0_phi_Str = basename + "_mu0_phi";
        TCanvas* mu0_phi_C = new TCanvas(mu0_phi_Str.c_str(), mu0_phi_Str.c_str(), X_SIZE, Y_SIZE);
        mu0_phi_->Draw();
        mu0_phi_C->Print((mu0_phi_Str+".png").c_str());

        std::string mu1_phi_Str = basename + "_mu1_phi";
        TCanvas* mu1_phi_C = new TCanvas(mu1_phi_Str.c_str(), mu1_phi_Str.c_str(), X_SIZE, Y_SIZE);
        mu1_phi_->Draw();
        mu1_phi_C->Print((mu1_phi_Str+".png").c_str());

        std::string mu0_charge_Str = basename + "_mu0_charge";
        TCanvas* mu0_charge_C = new TCanvas(mu0_charge_Str.c_str(), mu0_charge_Str.c_str(), X_SIZE, Y_SIZE);
        mu0_charge_->Draw();
        mu0_charge_C->Print((mu0_charge_Str+".png").c_str());

        std::string mu1_charge_Str = basename + "_mu1_charge";
        TCanvas* mu1_charge_C = new TCanvas(mu1_charge_Str.c_str(), mu1_charge_Str.c_str(), X_SIZE, Y_SIZE);
        mu1_charge_->Draw();
        mu1_charge_C->Print((mu1_charge_Str+".png").c_str());

        std::string nmuonsStr = basename + "_nmuons";
        TCanvas* nmuonsC = new TCanvas(nmuonsStr.c_str(), nmuonsStr.c_str(), X_SIZE, Y_SIZE);
        nmuons_->Draw();
        nmuonsC->Print((nmuonsStr+".png").c_str());

        std::string njpsisStr = basename + "_njpsis";
        TCanvas* njpsisC = new TCanvas(njpsisStr.c_str(), njpsisStr.c_str(), X_SIZE, Y_SIZE);
        njpsis_->Draw();
        njpsisC->Print((njpsisStr+".png").c_str());

        std::string pileupStr = basename + "_pileup";
        TCanvas* pileupC = new TCanvas(pileupStr.c_str(), pileupStr.c_str(), X_SIZE, Y_SIZE);
        pileup_->Draw();
        pileupC->Print((pileupStr+".png").c_str());
    }
}  // namespace zf
