#include "ZFinder/Event/interface/ZFinderPlotter.h"

// Standard Library
#include <string>  // string

// Root
#include <TCanvas.h>  // TCanvas

// ZFinder Code
#include "ZFinder/Event/interface/ZFinderElectron.h"  // ZFinderElectron
#include "ZFinder/Event/interface/ZFinderCuts.h"  // ZFinderCuts


namespace zf {
  // Constructor
  ZFinderPlotter::ZFinderPlotter(TFileDirectory& tdir, const bool USE_MC, const bool APPLY_MUON_MIN_PT, const bool APPLY_SOFT_MUONS, const bool APPLY_DIMUON_VTX_COMPATIBILITY, 
      const bool APPLY_JPSI_MASS_WINDOW , const bool APPLY_VERTEX_Z_POS_WINDOW)
    : USE_MC_(USE_MC), APPLY_MUON_MIN_PT_(APPLY_MUON_MIN_PT), APPLY_SOFT_MUONS_(APPLY_SOFT_MUONS), APPLY_DIMUON_VTX_COMPATIBILITY_(APPLY_DIMUON_VTX_COMPATIBILITY),
    APPLY_JPSI_MASS_WINDOW_(APPLY_JPSI_MASS_WINDOW), APPLY_VERTEX_Z_POS_WINDOW_(APPLY_VERTEX_Z_POS_WINDOW) {
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

    // z_vtx_prob_
    const std::string z_vtx_prob_name = "dielectron vertex probability";
    z_vtx_prob_ = tdir.make<TH1D>(z_vtx_prob_name.c_str(), z_vtx_prob_name.c_str(), 200, 0., 1.);
    z_vtx_prob_->GetXaxis()->SetTitle("Probability");
    z_vtx_prob_->GetYaxis()->SetTitle("Probability / 0.005 ");

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

    // jpsi_mass_all_
    const std::string jpsi_mass_all_name = "jpsi Mass: All";
    jpsi_mass_all_ = tdir.make<TH1D>(jpsi_mass_all_name.c_str(), jpsi_mass_all_name.c_str(), 300, 0., 300.);
    jpsi_mass_all_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
    jpsi_mass_all_->GetYaxis()->SetTitle("Counts / GeV");

    // jpsi_mass_coarse_
    const std::string jpsi_mass_coarse_name = "jpsi Mass: Coarse";
    jpsi_mass_coarse_ = tdir.make<TH1D>(jpsi_mass_coarse_name.c_str(), jpsi_mass_coarse_name.c_str(), 150, 0.0, 15.0);
    jpsi_mass_coarse_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
    jpsi_mass_coarse_->GetYaxis()->SetTitle("Counts / 0.1 GeV");

    // jpsi_mass_fine_
    const std::string jpsi_mass_fine_name = "jpsi Mass: Fine";
    jpsi_mass_fine_ = tdir.make<TH1D>(jpsi_mass_fine_name.c_str(), jpsi_mass_fine_name.c_str(), 100, 2.5, 3.5);
    jpsi_mass_fine_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
    jpsi_mass_fine_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

    // jpsi_rapidity_
    const std::string jpsi_rapidity_name = "jpsi Rapidity";
    jpsi_rapidity_ = tdir.make<TH1D>(jpsi_rapidity_name.c_str(), jpsi_rapidity_name.c_str(), 100, -5., 5.);
    jpsi_rapidity_->GetXaxis()->SetTitle("JPsi_{Y}");
    jpsi_rapidity_->GetYaxis()->SetTitle("Counts");

    // jpsi_pt
    const std::string jpsi_pt_name = "jpsi p_{T}";
    jpsi_pt_ = tdir.make<TH1D>(jpsi_pt_name.c_str(), jpsi_pt_name.c_str(), 200, 0., 200.);
    jpsi_pt_->GetXaxis()->SetTitle("p_{T,JPsi}");
    jpsi_pt_->GetYaxis()->SetTitle("Counts / GeV");

    // jpsi_pt
    const std::string jpsi_pt_vs_rap_name = "jpsi_pt_vs_rap";
    jpsi_pt_vs_rap_ = tdir.make<TH2D>(jpsi_pt_vs_rap_name.c_str(), jpsi_pt_vs_rap_name.c_str(), 100, -5., 5., 200, 0., 200.);
    jpsi_pt_vs_rap_->GetXaxis()->SetTitle("jpsi Rapidity");
    jpsi_pt_vs_rap_->GetYaxis()->SetTitle("jpsi Pt");

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
    jpsi_vtx_distance_z_vtx_z_ = tdir.make<TH1D>(jpsi_vtx_distance_z_vtx_z_name.c_str(), jpsi_vtx_distance_z_vtx_z_name.c_str(), 2000, -20., 20.);
    jpsi_vtx_distance_z_vtx_z_->GetXaxis()->SetTitle("distance [cm]");
    jpsi_vtx_distance_z_vtx_z_->GetYaxis()->SetTitle("Counts / 0.02 cm");

    // jpsi_distance
    const std::string jpsi_distance_name = "jpsi distance";
    jpsi_distance_ = tdir.make<TH1D>(jpsi_distance_name.c_str(), jpsi_distance_name.c_str(), 500, 0., 5.);
    jpsi_distance_->GetXaxis()->SetTitle("distance [cm]");
    jpsi_distance_->GetYaxis()->SetTitle("Counts / 0.01 cm");

    // jpsi_dist_err
    const std::string jpsi_dist_err_name = "jpsi dist_err";
    jpsi_dist_err_ = tdir.make<TH1D>(jpsi_dist_err_name.c_str(), jpsi_dist_err_name.c_str(), 500, 0., 5.);
    jpsi_dist_err_->GetXaxis()->SetTitle("dist_err [cm]");
    jpsi_dist_err_->GetYaxis()->SetTitle("Counts / 0.01 cm");

    // jpsi_chi2
    const std::string jpsi_chi2_name = "jpsi chi2";
    jpsi_chi2_ = tdir.make<TH1D>(jpsi_chi2_name.c_str(), jpsi_chi2_name.c_str(), 1000, 0., 100.);
    jpsi_chi2_->GetXaxis()->SetTitle("chi2");
    jpsi_chi2_->GetYaxis()->SetTitle("Counts / 0.1 ");

    // jpsi_distance_xy
    const std::string jpsi_distance_xy_name = "jpsi distance_xy";
    jpsi_distance_xy_ = tdir.make<TH1D>(jpsi_distance_xy_name.c_str(), jpsi_distance_xy_name.c_str(), 500, 0., 5.);
    jpsi_distance_xy_->GetXaxis()->SetTitle("distance_xy [cm]");
    jpsi_distance_xy_->GetYaxis()->SetTitle("Counts / 0.01 cm");

    // jpsi_dist_err_xy
    const std::string jpsi_dist_err_xy_name = "jpsi dist_err_xy";
    jpsi_dist_err_xy_ = tdir.make<TH1D>(jpsi_dist_err_xy_name.c_str(), jpsi_dist_err_xy_name.c_str(), 500, 0., 5.);
    jpsi_dist_err_xy_->GetXaxis()->SetTitle("dist_err_xy [cm]");
    jpsi_dist_err_xy_->GetYaxis()->SetTitle("Counts / 0.01 cm");

    // jpsi_chi2_xy
    const std::string jpsi_chi2_xy_name = "jpsi chi2_xy";
    jpsi_chi2_xy_ = tdir.make<TH1D>(jpsi_chi2_xy_name.c_str(), jpsi_chi2_xy_name.c_str(), 1000, 0., 100.);
    jpsi_chi2_xy_->GetXaxis()->SetTitle("chi2_xy");
    jpsi_chi2_xy_->GetYaxis()->SetTitle("Counts / 0.1 ");

    // jpsi_zpt_difference_
    const std::string jpsi_zpt_difference__name = "zPt - jpsipT";
    jpsi_zpt_difference_ = tdir.make<TH1D>(jpsi_zpt_difference__name.c_str(), jpsi_zpt_difference__name.c_str(), 200, -100., 100.);
    jpsi_zpt_difference_->GetXaxis()->SetTitle("pT [GeV]");
    jpsi_zpt_difference_->GetYaxis()->SetTitle("Counts / 1 GeV ");

    // jpsi_tau_xy
    const std::string jpsi_tau_xy_name = "jpsi tau_xy";
    jpsi_tau_xy_ = tdir.make<TH1D>(jpsi_tau_xy_name.c_str(), jpsi_tau_xy_name.c_str(), 200, -100., 100.);
    jpsi_tau_xy_->GetXaxis()->SetTitle("tau_xy [ps]");
    jpsi_tau_xy_->GetYaxis()->SetTitle("Counts / 1 ps ");

    // jpsi_tau_xy_fine
    const std::string jpsi_tau_xy_fine_name = "jpsi tau_xy_fine";
    jpsi_tau_xy_fine_ = tdir.make<TH1D>(jpsi_tau_xy_fine_name.c_str(), jpsi_tau_xy_fine_name.c_str(), 200, -10., 10.);
    jpsi_tau_xy_fine_->GetXaxis()->SetTitle("tau_xy_fine [ps]");
    jpsi_tau_xy_fine_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

    // jpsi_tau_xy_very_fine
    const std::string jpsi_tau_xy_very_fine_name = "jpsi_tau_xy_very_fine_all";
    jpsi_tau_xy_very_fine_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_name.c_str(), jpsi_tau_xy_very_fine_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_ptUnder10
    const std::string jpsi_tau_xy_very_fine_ptUnder10_name = "jpsi_tau_xy_very_fine_ptUnder10";
    jpsi_tau_xy_very_fine_ptUnder10_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_ptUnder10_name.c_str(), jpsi_tau_xy_very_fine_ptUnder10_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_ptUnder10_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_ptUnder10_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_pt10to15
    const std::string jpsi_tau_xy_very_fine_pt10to15_name = "jpsi_tau_xy_very_fine_pt_10_to_15";
    jpsi_tau_xy_very_fine_pt10to15_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_pt10to15_name.c_str(), jpsi_tau_xy_very_fine_pt10to15_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_pt10to15_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_pt10to15_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_pt15to20
    const std::string jpsi_tau_xy_very_fine_pt15to20_name = "jpsi_tau_xy_very_fine_pt_15_to_20";
    jpsi_tau_xy_very_fine_pt15to20_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_pt15to20_name.c_str(), jpsi_tau_xy_very_fine_pt15to20_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_pt15to20_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_pt15to20_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_pt20to25
    const std::string jpsi_tau_xy_very_fine_pt20to25_name = "jpsi_tau_xy_very_fine_pt_20_to_25";
    jpsi_tau_xy_very_fine_pt20to25_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_pt20to25_name.c_str(), jpsi_tau_xy_very_fine_pt20to25_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_pt20to25_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_pt20to25_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_pt25to30
    const std::string jpsi_tau_xy_very_fine_pt25to30_name = "jpsi_tau_xy_very_fine_pt_25_to_30";
    jpsi_tau_xy_very_fine_pt25to30_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_pt25to30_name.c_str(), jpsi_tau_xy_very_fine_pt25to30_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_pt25to30_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_pt25to30_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_ptAbove20
    const std::string jpsi_tau_xy_very_fine_ptAbove20_name = "jpsi_tau_xy_very_fine_ptAbove20";
    jpsi_tau_xy_very_fine_ptAbove20_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_ptAbove20_name.c_str(), jpsi_tau_xy_very_fine_ptAbove20_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_ptAbove20_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_ptAbove20_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_ptAbove30
    const std::string jpsi_tau_xy_very_fine_ptAbove30_name = "jpsi_tau_xy_very_fine_ptAbove30";
    jpsi_tau_xy_very_fine_ptAbove30_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_ptAbove30_name.c_str(), jpsi_tau_xy_very_fine_ptAbove30_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_ptAbove30_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_ptAbove30_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap0_0to0_3
    const std::string jpsi_tau_xy_very_fine_rap0_0to0_3_name = "jpsi_tau_xy_very_fine_rap_0.0_to_0.3";
    jpsi_tau_xy_very_fine_rap0_0to0_3_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_0to0_3_name.c_str(), jpsi_tau_xy_very_fine_rap0_0to0_3_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap0_0to0_3_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap0_0to0_3_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap0_3to0_6
    const std::string jpsi_tau_xy_very_fine_rap0_3to0_6_name = "jpsi_tau_xy_very_fine_rap_0.3_to_0.6";
    jpsi_tau_xy_very_fine_rap0_3to0_6_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_3to0_6_name.c_str(), jpsi_tau_xy_very_fine_rap0_3to0_6_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap0_3to0_6_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap0_3to0_6_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap0_6to0_9
    const std::string jpsi_tau_xy_very_fine_rap0_6to0_9_name = "jpsi_tau_xy_very_fine_rap_0.6_to_0.9";
    jpsi_tau_xy_very_fine_rap0_6to0_9_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_6to0_9_name.c_str(), jpsi_tau_xy_very_fine_rap0_6to0_9_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap0_6to0_9_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap0_6to0_9_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap0_9to1_2
    const std::string jpsi_tau_xy_very_fine_rap0_9to1_2_name = "jpsi_tau_xy_very_fine_rap_0.9_to_1.2";
    jpsi_tau_xy_very_fine_rap0_9to1_2_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_9to1_2_name.c_str(), jpsi_tau_xy_very_fine_rap0_9to1_2_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap0_9to1_2_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap0_9to1_2_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap1_2to1_5
    const std::string jpsi_tau_xy_very_fine_rap1_2to1_5_name = "jpsi_tau_xy_very_fine_rap_1.2_to_1.5";
    jpsi_tau_xy_very_fine_rap1_2to1_5_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap1_2to1_5_name.c_str(), jpsi_tau_xy_very_fine_rap1_2to1_5_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap1_2to1_5_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap1_2to1_5_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap1_5to1_8
    const std::string jpsi_tau_xy_very_fine_rap1_5to1_8_name = "jpsi_tau_xy_very_fine_rap_1.5_to_1.8";
    jpsi_tau_xy_very_fine_rap1_5to1_8_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap1_5to1_8_name.c_str(), jpsi_tau_xy_very_fine_rap1_5to1_8_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap1_5to1_8_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap1_5to1_8_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap1_8to2_1
    const std::string jpsi_tau_xy_very_fine_rap1_8to2_1_name = "jpsi_tau_xy_very_fine_rap_1.8_to_2.1";
    jpsi_tau_xy_very_fine_rap1_8to2_1_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap1_8to2_1_name.c_str(), jpsi_tau_xy_very_fine_rap1_8to2_1_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap1_8to2_1_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap1_8to2_1_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap2_1to2_4
    const std::string jpsi_tau_xy_very_fine_rap2_1to2_4_name = "jpsi_tau_xy_very_fine_rap_2.1_to_2.4";
    jpsi_tau_xy_very_fine_rap2_1to2_4_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap2_1to2_4_name.c_str(), jpsi_tau_xy_very_fine_rap2_1to2_4_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap2_1to2_4_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap2_1to2_4_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap0_0to0_9_pt10to15
    const std::string jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_name = "jpsi_tau_xy_very_fine_rap_0.0_to_0.9_pt10to15";
    jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_name.c_str(), jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap0_9to1_2_pt10to15
    const std::string jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_name = "jpsi_tau_xy_very_fine_rap_0.9_to_1.2_pt10to15";
    jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_name.c_str(), jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_rap1_2to2_1_pt10to15
    const std::string jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_name = "jpsi_tau_xy_very_fine_rap_1.2_to_2.1_pt10to15";
    jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_name.c_str(), jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_above_12_tracker_layers
    const std::string jpsi_tau_xy_very_fine_above_12_tracker_layers_name = "jpsi_tau_xy_very_fine_above_12_tracker_layers_";
    jpsi_tau_xy_very_fine_above_12_tracker_layers_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_above_12_tracker_layers_name.c_str(), jpsi_tau_xy_very_fine_above_12_tracker_layers_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_above_12_tracker_layers_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_above_12_tracker_layers_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_xy_very_fine_similar_pt_muons_
    const std::string jpsi_tau_xy_very_fine_similar_pt_muons_name = "jpsi_tau_xy_very_fine_similar_pt_muons_";
    jpsi_tau_xy_very_fine_similar_pt_muons_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_similar_pt_muons_name.c_str(), jpsi_tau_xy_very_fine_similar_pt_muons_name.c_str(), 2000, -10., 10.);
    jpsi_tau_xy_very_fine_similar_pt_muons_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
    jpsi_tau_xy_very_fine_similar_pt_muons_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_tau_z
    const std::string jpsi_tau_z_name = "jpsi tau_z";
    jpsi_tau_z_ = tdir.make<TH1D>(jpsi_tau_z_name.c_str(), jpsi_tau_z_name.c_str(), 200, -100., 100.);
    jpsi_tau_z_->GetXaxis()->SetTitle("tau_z [ps]");
    jpsi_tau_z_->GetYaxis()->SetTitle("Counts / 1 ps ");

    // jpsi_tau_z_fine
    const std::string jpsi_tau_z_fine_name = "jpsi tau_z_fine";
    jpsi_tau_z_fine_ = tdir.make<TH1D>(jpsi_tau_z_fine_name.c_str(), jpsi_tau_z_fine_name.c_str(), 200, -10., 10.);
    jpsi_tau_z_fine_->GetXaxis()->SetTitle("tau_z_fine [ps]");
    jpsi_tau_z_fine_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

    // jpsi_tau_z_very_fine
    const std::string jpsi_tau_z_very_fine_name = "jpsi tau_z_very_fine";
    jpsi_tau_z_very_fine_ = tdir.make<TH1D>(jpsi_tau_z_very_fine_name.c_str(), jpsi_tau_z_very_fine_name.c_str(), 2000, -10., 10.);
    jpsi_tau_z_very_fine_->GetXaxis()->SetTitle("tau_z_very_fine [ps]");
    jpsi_tau_z_very_fine_->GetYaxis()->SetTitle("Counts / 0.01 ps ");

    // jpsi_mass_vs_chi2
    const std::string jpsi_mass_vs_chi2_name = "jpsi mass vs chi2";
    jpsi_mass_vs_chi2_ = tdir.make<TH2D>(jpsi_mass_vs_chi2_name.c_str(), jpsi_mass_vs_chi2_name.c_str(), 100, 0., 10., 100, 0., 10.);
    jpsi_mass_vs_chi2_->GetXaxis()->SetTitle("chi2");
    jpsi_mass_vs_chi2_->GetYaxis()->SetTitle("mass");

    // jpsi_tau_xy_vs_tau_z
    const std::string jpsi_tau_xy_vs_tau_z_name = "jpsi tau_xy vs tau_z";
    jpsi_tau_xy_vs_tau_z_ = tdir.make<TH2D>(jpsi_tau_xy_vs_tau_z_name.c_str(), jpsi_tau_xy_vs_tau_z_name.c_str(), 200, -10., 10., 200, -10., 10.);
    jpsi_tau_xy_vs_tau_z_->GetXaxis()->SetTitle("tau_xy [ps]");
    jpsi_tau_xy_vs_tau_z_->GetYaxis()->SetTitle("tau_z [ps]");

    // jpsi_tau_xy_vs_distance_z
    const std::string jpsi_tau_xy_vs_distance_z_name = "jpsi tau_xy vs z vertex difference";
    jpsi_tau_xy_vs_distance_z_ = tdir.make<TH2D>(jpsi_tau_xy_vs_distance_z_name.c_str(), jpsi_tau_xy_vs_distance_z_name.c_str(), 200, -10., 10., 200, -10., 10.);
    jpsi_tau_xy_vs_distance_z_->GetXaxis()->SetTitle("tau_xy [ps]");
    jpsi_tau_xy_vs_distance_z_->GetYaxis()->SetTitle("dZ [cm]");

    // jpsi_tau_z_vs_distance_z
    const std::string jpsi_tau_z_vs_distance_z_name = "jpsi tau_z vs z vertex difference";
    jpsi_tau_z_vs_distance_z_ = tdir.make<TH2D>(jpsi_tau_z_vs_distance_z_name.c_str(), jpsi_tau_z_vs_distance_z_name.c_str(), 200, -10., 10., 200, -10., 10.);
    jpsi_tau_z_vs_distance_z_->GetXaxis()->SetTitle("tau_z [ps]");
    jpsi_tau_z_vs_distance_z_->GetYaxis()->SetTitle("dZ [cm]");

    //TODO remove 0 from jpsi
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

    // dimuon_vtx_prob_
    const std::string dimuon_vtx_prob_name = "dimuon vertex probability";
    dimuon_vtx_prob_ = tdir.make<TH1D>(dimuon_vtx_prob_name.c_str(), dimuon_vtx_prob_name.c_str(), 200, 0., 1.);
    dimuon_vtx_prob_->GetXaxis()->SetTitle("Probability");
    dimuon_vtx_prob_->GetYaxis()->SetTitle("Probability / 0.005 ");

    // dimuon_delta_phi_
    const std::string dimuon_delta_phi_name = "dimuon delta phi";
    dimuon_delta_phi_ = tdir.make<TH1D>(dimuon_delta_phi_name.c_str(), dimuon_delta_phi_name.c_str(), 35, 0.0, 3.5);
    dimuon_delta_phi_->GetXaxis()->SetTitle("Angle [rad]");
    dimuon_delta_phi_->GetYaxis()->SetTitle("delta Phi / 0.1");

    // dimuon_delta_eta_
    const std::string dimuon_delta_eta_name = "dimuon delta eta";
    dimuon_delta_eta_ = tdir.make<TH1D>(dimuon_delta_eta_name.c_str(), dimuon_delta_eta_name.c_str(), 50, 0.0, 5.);
    dimuon_delta_eta_->GetXaxis()->SetTitle("Eta");
    dimuon_delta_eta_->GetYaxis()->SetTitle("delta Eta / 0.1");

    // dimuon_deltaR_
    const std::string dimuon_deltaR_name = "dimuon deltaR";
    dimuon_deltaR_ = tdir.make<TH1D>(dimuon_deltaR_name.c_str(), dimuon_deltaR_name.c_str(), 100, 0.0, 10.);
    dimuon_deltaR_->GetXaxis()->SetTitle("deltaR");
    dimuon_deltaR_->GetYaxis()->SetTitle("deltaR / 0.1");

    // mu0_pt
    const std::string mu0_pt_name = "p_{T,mu_{0}}";
    mu0_pt_ = tdir.make<TH1D>(mu0_pt_name.c_str(), mu0_pt_name.c_str(), 400, 0., 200.);
    mu0_pt_->GetXaxis()->SetTitle("p_{T,mu_{0}}");
    mu0_pt_->GetYaxis()->SetTitle("Counts / 0.5 GeV");

    // mu1_pt
    const std::string mu1_pt_name = "p_{T,mu_{1}}";
    mu1_pt_ = tdir.make<TH1D>(mu1_pt_name.c_str(), mu1_pt_name.c_str(), 400, 0., 200.);
    mu1_pt_->GetXaxis()->SetTitle("p_{T,mu_{0}}");
    mu1_pt_->GetYaxis()->SetTitle("Counts / 0.5 GeV");

    // mu0_eta_
    const std::string mu0_eta_name = "#eta_{mu_{0}}";
    mu0_eta_ = tdir.make<TH1D>(mu0_eta_name.c_str(), mu0_eta_name.c_str(), 50, -5., 5.);
    mu0_eta_->GetXaxis()->SetTitle("#eta_{mu_{0}}");
    mu0_eta_->GetYaxis()->SetTitle("Counts / 0.2");

    // mu1_eta_
    const std::string mu1_eta_name = "#eta_{mu_{1}}";
    mu1_eta_ = tdir.make<TH1D>(mu1_eta_name.c_str(), mu1_eta_name.c_str(), 50, -5., 5.);
    mu1_eta_->GetXaxis()->SetTitle("#eta_{mu_{1}}");
    mu1_eta_->GetYaxis()->SetTitle("Counts / 0.2");

    // mu0_phi_
    const std::string mu0_phi_name = "#phi_{mu_{0}}";
    mu0_phi_ = tdir.make<TH1D>(mu0_phi_name.c_str(), mu0_phi_name.c_str(), 60, -3.15, 3.15);
    mu0_phi_->GetXaxis()->SetTitle("#phi_{mu_{0}}");
    mu0_phi_->GetYaxis()->SetTitle("Counts / 0.105");

    // mu1_phi_
    const std::string mu1_phi_name = "#phi_{mu_{1}}";
    mu1_phi_ = tdir.make<TH1D>(mu1_phi_name.c_str(), mu1_phi_name.c_str(), 50, -3.15, 3.15);
    mu1_phi_->GetXaxis()->SetTitle("#phi_{mu_{1}}");
    mu1_phi_->GetYaxis()->SetTitle("Counts / 0.105");

    // mu0_charge_
    const std::string mu0_charge_name = "charge_{mu_{0}}";
    mu0_charge_ = tdir.make<TH1D>(mu0_charge_name.c_str(), mu0_charge_name.c_str(), 60, -3.15, 3.15);
    mu0_charge_->GetXaxis()->SetTitle("charge_{mu_{0}}");
    mu0_charge_->GetYaxis()->SetTitle("Counts");

    // mu1_charge_
    const std::string mu1_charge_name = "charge_{mu_{1}}";
    mu1_charge_ = tdir.make<TH1D>(mu1_charge_name.c_str(), mu1_charge_name.c_str(), 60, -3.15, 3.15);
    mu1_charge_->GetXaxis()->SetTitle("charge_{mu_{1}}");
    mu1_charge_->GetYaxis()->SetTitle("Counts");

    // mu0_deltaR_truth_
    const std::string mu0_deltaR_truth_name = "deltaR_truth_mu0";
    mu0_deltaR_truth_ = tdir.make<TH1D>(mu0_deltaR_truth_name.c_str(), mu0_deltaR_truth_name.c_str(), 500, 0, 5);
    mu0_deltaR_truth_->GetXaxis()->SetTitle("deltaR_truth_{mu_{0}}");
    mu0_deltaR_truth_->GetYaxis()->SetTitle("Counts / 0.01");

    // mu1_deltaR_truth_
    const std::string mu1_deltaR_truth_name = "deltaR_truth_mu1";
    mu1_deltaR_truth_ = tdir.make<TH1D>(mu1_deltaR_truth_name.c_str(), mu1_deltaR_truth_name.c_str(), 500, 0, 5);
    mu1_deltaR_truth_->GetXaxis()->SetTitle("deltaR_truth_{mu_{1}}");
    mu1_deltaR_truth_->GetYaxis()->SetTitle("Counts / 0.01");

    // n_truth_matched_jpsi_muons
    const std::string n_truth_matched_jpsi_muons_name = "N_truth_matched_jpsi_muons";
    n_truth_matched_jpsi_muons_ = tdir.make<TH1D>(n_truth_matched_jpsi_muons_name.c_str(), n_truth_matched_jpsi_muons_name.c_str(), 3, 0, 3);
    n_truth_matched_jpsi_muons_->GetXaxis()->SetTitle("Number of Truth Matched Jpsi Muons");
    n_truth_matched_jpsi_muons_->GetYaxis()->SetTitle("Events");

    // mu0_tracker_layers_
    const std::string mu0_tracker_layers_name = "tracker_layers_{mu_{0}}";
    mu0_tracker_layers_ = tdir.make<TH1D>(mu0_tracker_layers_name.c_str(), mu0_tracker_layers_name.c_str(), 20, 0, 20);
    mu0_tracker_layers_->GetXaxis()->SetTitle("tracker_layers_{mu_{0}}");
    mu0_tracker_layers_->GetYaxis()->SetTitle("Counts");

    // mu1_tracker_layers_
    const std::string mu1_tracker_layers_name = "tracker_layers_{mu_{1}}";
    mu1_tracker_layers_ = tdir.make<TH1D>(mu1_tracker_layers_name.c_str(), mu1_tracker_layers_name.c_str(), 20, 0, 20);
    mu1_tracker_layers_->GetXaxis()->SetTitle("tracker_layers_{mu_{1}}");
    mu1_tracker_layers_->GetYaxis()->SetTitle("Counts");

    // jet_pt
    const std::string jet_pt_name = "p_{T,jet}";
    jet_pt_ = tdir.make<TH1D>(jet_pt_name.c_str(), jet_pt_name.c_str(), 200, 0., 200.);
    jet_pt_->GetXaxis()->SetTitle("p_{T,jet}");
    jet_pt_->GetYaxis()->SetTitle("Counts / GeV");

    // jet_eta_
    const std::string jet_eta_name = "#eta_{jet}";
    jet_eta_ = tdir.make<TH1D>(jet_eta_name.c_str(), jet_eta_name.c_str(), 50, -5., 5.);
    jet_eta_->GetXaxis()->SetTitle("#eta_{jet}");
    jet_eta_->GetYaxis()->SetTitle("Counts");

    // jet_btag_discriminator_
    const std::string jet_btag_discriminator_name = "jet btag discriminator";
    jet_btag_discriminator_ = tdir.make<TH1D>(jet_btag_discriminator_name.c_str(), jet_btag_discriminator_name.c_str(), 1400, -120., 20.);
    jet_btag_discriminator_->GetXaxis()->SetTitle("discriminator");
    jet_btag_discriminator_->GetYaxis()->SetTitle("Counts / 0.1");

    // muon_jet_pt
    const std::string muon_jet_pt_name = "p_{T,muon_jet}";
    muon_jet_pt_ = tdir.make<TH1D>(muon_jet_pt_name.c_str(), muon_jet_pt_name.c_str(), 200, 0., 200.);
    muon_jet_pt_->GetXaxis()->SetTitle("p_{T,muon_jet}");
    muon_jet_pt_->GetYaxis()->SetTitle("Counts / GeV");

    // muon_jet_pt_diff_z_pt_
    const std::string muon_jet_pt_diff_z_pt_name = "z p_{T} - muon jet p_{T}";
    muon_jet_pt_diff_z_pt_ = tdir.make<TH1D>(muon_jet_pt_diff_z_pt_name.c_str(), muon_jet_pt_diff_z_pt_name.c_str(), 400, -200., 200.);
    muon_jet_pt_diff_z_pt_->GetXaxis()->SetTitle("p_{T,muon_jet}");
    muon_jet_pt_diff_z_pt_->GetYaxis()->SetTitle("Counts / GeV");

    // muon_jet_pt_diff_dimuon_pt_
    const std::string muon_jet_pt_diff_dimuon_pt_name = "dimuon p_{T} - muon jet p_{T}";
    muon_jet_pt_diff_dimuon_pt_ = tdir.make<TH1D>(muon_jet_pt_diff_dimuon_pt_name.c_str(), muon_jet_pt_diff_dimuon_pt_name.c_str(), 400, -200., 200.);
    muon_jet_pt_diff_dimuon_pt_->GetXaxis()->SetTitle("p_{T,muon_jet}");
    muon_jet_pt_diff_dimuon_pt_->GetYaxis()->SetTitle("Counts / GeV");

    // muon_jet_pt_z_pt_
    const std::string muon_jet_pt_z_pt_name = "muon jet p_{T} vs z p_{T}";
    muon_jet_pt_z_pt_ = tdir.make<TH2D>(muon_jet_pt_z_pt_name.c_str(), muon_jet_pt_z_pt_name.c_str(), 200, 0., 200., 200, 0., 200.);
    muon_jet_pt_z_pt_->GetXaxis()->SetTitle("p_{T,z}");
    muon_jet_pt_z_pt_->GetYaxis()->SetTitle("p_{T,muon_jet}");

    // muon_jet_pt_dimuon_pt_
    const std::string muon_jet_pt_dimuon_pt_name = "muon jet p_{T} vs dimuon p_{T}";
    muon_jet_pt_dimuon_pt_ = tdir.make<TH2D>(muon_jet_pt_dimuon_pt_name.c_str(), muon_jet_pt_dimuon_pt_name.c_str(), 200, 0., 200., 200, 0., 200.);
    muon_jet_pt_dimuon_pt_->GetXaxis()->SetTitle("p_{T,dimuon}");
    muon_jet_pt_dimuon_pt_->GetYaxis()->SetTitle("p_{T,muon_jet}");

    // muon_jet_phi_z_phi_
    const std::string muon_jet_phi_z_phi_name = "muon jet phi vs z phi";
    muon_jet_phi_z_phi_ = tdir.make<TH2D>(muon_jet_phi_z_phi_name.c_str(), muon_jet_phi_z_phi_name.c_str(), 8, -4., 4., 8, -4., 4.);
    muon_jet_phi_z_phi_->GetXaxis()->SetTitle("phi z");
    muon_jet_phi_z_phi_->GetYaxis()->SetTitle("phi muon_jet");

    // muon_jet_phi_dimuon_phi_
    const std::string muon_jet_phi_dimuon_phi_name = "muon jet phi vs dimuon phi";
    muon_jet_phi_dimuon_phi_ = tdir.make<TH2D>(muon_jet_phi_dimuon_phi_name.c_str(), muon_jet_phi_dimuon_phi_name.c_str(), 8, -4., 4., 8, -4., 4.);
    muon_jet_phi_dimuon_phi_->GetXaxis()->SetTitle("phi dimuon");
    muon_jet_phi_dimuon_phi_->GetYaxis()->SetTitle("phi muon_jet");

    // muon_jet_eta_
    const std::string muon_jet_eta_name = "#eta_{muon_jet}";
    muon_jet_eta_ = tdir.make<TH1D>(muon_jet_eta_name.c_str(), muon_jet_eta_name.c_str(), 50, -5., 5.);
    muon_jet_eta_->GetXaxis()->SetTitle("#eta_{muon_jet}");
    muon_jet_eta_->GetYaxis()->SetTitle("Counts");

    // muon_jet_btag_discriminator_
    const std::string muon_jet_btag_discriminator_name = "muon jet btag discriminator";
    muon_jet_btag_discriminator_ = tdir.make<TH1D>(muon_jet_btag_discriminator_name.c_str(), muon_jet_btag_discriminator_name.c_str(), 1400, -120., 20.);
    muon_jet_btag_discriminator_->GetXaxis()->SetTitle("discriminator");
    muon_jet_btag_discriminator_->GetYaxis()->SetTitle("Counts / 0.1");

    // vtx_x_
    const std::string vtx_x_name = "vertex position x";
    vtx_x_ = tdir.make<TH1D>(vtx_x_name.c_str(), vtx_x_name.c_str(), 2000, -10., 10.);
    vtx_x_->GetXaxis()->SetTitle("Position [cm]");
    vtx_x_->GetYaxis()->SetTitle("Vertexes / 0.001 cm");

    // vtx_y_
    const std::string vtx_y_name = "vertex position y";
    vtx_y_ = tdir.make<TH1D>(vtx_y_name.c_str(), vtx_y_name.c_str(), 2000, -10., 10.);
    vtx_y_->GetXaxis()->SetTitle("Position [cm]");
    vtx_y_->GetYaxis()->SetTitle("Vertexes / 0.001 cm");

    // vtx_z_
    const std::string vtx_z_name = "vertex position z";
    vtx_z_ = tdir.make<TH1D>(vtx_z_name.c_str(), vtx_z_name.c_str(), 200, -10., 10.);
    vtx_z_->GetXaxis()->SetTitle("Position [cm]");
    vtx_z_->GetYaxis()->SetTitle("Vertexes / 0.01 cm");

    // primary_vtx_x_
    const std::string primary_vtx_x_name = "primary vertex position x";
    primary_vtx_x_ = tdir.make<TH1D>(primary_vtx_x_name.c_str(), primary_vtx_x_name.c_str(), 2000, -10., 10.);
    primary_vtx_x_->GetXaxis()->SetTitle("Position [cm]");
    primary_vtx_x_->GetYaxis()->SetTitle("Primary Vertex / 0.001 cm");

    // primary_vtx_y_
    const std::string primary_vtx_y_name = "primary vertex position y";
    primary_vtx_y_ = tdir.make<TH1D>(primary_vtx_y_name.c_str(), primary_vtx_y_name.c_str(), 2000, -10., 10.);
    primary_vtx_y_->GetXaxis()->SetTitle("Position [cm]");
    primary_vtx_y_->GetYaxis()->SetTitle("Primary Vertex / 0.001 cm");

    // primary_vtx_z_
    const std::string primary_vtx_z_name = "primary vertex position z";
    primary_vtx_z_ = tdir.make<TH1D>(primary_vtx_z_name.c_str(), primary_vtx_z_name.c_str(), 200, -10., 10.);
    primary_vtx_z_->GetXaxis()->SetTitle("Position [cm]");
    primary_vtx_z_->GetYaxis()->SetTitle("Primary Vertex / 0.01 cm");

    // z_vtx_x_
    const std::string z_vtx_x_name = "z vertex position x";
    z_vtx_x_ = tdir.make<TH1D>(z_vtx_x_name.c_str(), z_vtx_x_name.c_str(), 2000, -10., 10.);
    z_vtx_x_->GetXaxis()->SetTitle("Position [cm]");
    z_vtx_x_->GetYaxis()->SetTitle("z Vertex / 0.001 cm");

    // z_vtx_y_
    const std::string z_vtx_y_name = "z vertex position y";
    z_vtx_y_ = tdir.make<TH1D>(z_vtx_y_name.c_str(), z_vtx_y_name.c_str(), 2000, -10., 10.);
    z_vtx_y_->GetXaxis()->SetTitle("Position [cm]");
    z_vtx_y_->GetYaxis()->SetTitle("z Vertex / 0.001 cm");

    // z_vtx_z_
    const std::string z_vtx_z_name = "z vertex position z";
    z_vtx_z_ = tdir.make<TH1D>(z_vtx_z_name.c_str(), z_vtx_z_name.c_str(), 200, -10., 10.);
    z_vtx_z_->GetXaxis()->SetTitle("Position [cm]");
    z_vtx_z_->GetYaxis()->SetTitle("z Vertex / 0.01 cm");

    // primary_vtx_x_vs_z_vtx_x_
    const std::string primary_vtx_x_vs_z_vtx_x_name = "primary vertex x pos vs Z vertex x pos";
    primary_vtx_x_vs_z_vtx_x_ = tdir.make<TH2D>(primary_vtx_x_vs_z_vtx_x_name.c_str(), primary_vtx_x_vs_z_vtx_x_name.c_str(), 200, -1., 1., 200, -1., 1.);
    primary_vtx_x_vs_z_vtx_x_->GetXaxis()->SetTitle("Z Vertex x Pos [cm]");
    primary_vtx_x_vs_z_vtx_x_->GetYaxis()->SetTitle("Primary Vertex x Pos [cm]");

    // primary_vtx_y_vs_z_vtx_y_
    const std::string primary_vtx_y_vs_z_vtx_y_name = "primary vertex y pos vs Z vertex y pos";
    primary_vtx_y_vs_z_vtx_y_ = tdir.make<TH2D>(primary_vtx_y_vs_z_vtx_y_name.c_str(), primary_vtx_y_vs_z_vtx_y_name.c_str(), 200, -1., 1., 200, -1., 1.);
    primary_vtx_y_vs_z_vtx_y_->GetXaxis()->SetTitle("Z Vertex y Pos [cm]");
    primary_vtx_y_vs_z_vtx_y_->GetYaxis()->SetTitle("Primary Vertex y Pos [cm]");

    // primary_vtx_z_vs_z_vtx_z_
    const std::string primary_vtx_z_vs_z_vtx_z_name = "primary vertex z pos vs Z vertex z pos";
    primary_vtx_z_vs_z_vtx_z_ = tdir.make<TH2D>(primary_vtx_z_vs_z_vtx_z_name.c_str(), primary_vtx_z_vs_z_vtx_z_name.c_str(), 200, -10., 10., 200, -10., 10.);
    primary_vtx_z_vs_z_vtx_z_->GetXaxis()->SetTitle("Z Vertex z Pos [cm]");
    primary_vtx_z_vs_z_vtx_z_->GetYaxis()->SetTitle("Primary Vertex z Pos [cm]");

    // dimuon_vtx_x_
    const std::string dimuon_vtx_x_name = "dimuon vertex position x";
    dimuon_vtx_x_ = tdir.make<TH1D>(dimuon_vtx_x_name.c_str(), dimuon_vtx_x_name.c_str(), 2000, -10., 10.);
    dimuon_vtx_x_->GetXaxis()->SetTitle("Position [cm]");
    dimuon_vtx_x_->GetYaxis()->SetTitle("dimuon Vertex / 0.001 cm");

    // dimuon_vtx_y_
    const std::string dimuon_vtx_y_name = "dimuon vertex position y";
    dimuon_vtx_y_ = tdir.make<TH1D>(dimuon_vtx_y_name.c_str(), dimuon_vtx_y_name.c_str(), 2000, -10., 10.);
    dimuon_vtx_y_->GetXaxis()->SetTitle("Position [cm]");
    dimuon_vtx_y_->GetYaxis()->SetTitle("dimuon Vertex / 0.001 cm");

    // dimuon_vtx_z_
    const std::string dimuon_vtx_z_name = "dimuon vertex position z";
    dimuon_vtx_z_ = tdir.make<TH1D>(dimuon_vtx_z_name.c_str(), dimuon_vtx_z_name.c_str(), 200, -10., 10.);
    dimuon_vtx_z_->GetXaxis()->SetTitle("Position [cm]");
    dimuon_vtx_z_->GetYaxis()->SetTitle("dimuon Vertex / 0.01 cm");

    // z_jpsi_delta_phi_
    const std::string z_jpsi_delta_phi_name = "z jpsi delta phi";
    z_jpsi_delta_phi_ = tdir.make<TH1D>(z_jpsi_delta_phi_name.c_str(), z_jpsi_delta_phi_name.c_str(), 35, 0.0, 3.5);
    z_jpsi_delta_phi_->GetXaxis()->SetTitle("Angle [rad]");
    z_jpsi_delta_phi_->GetYaxis()->SetTitle("delta Phi / 0.1");

    // nmuons
    const std::string nmuons_name = "N_{mu}";
    nmuons_ = tdir.make<TH1D>(nmuons_name.c_str(), nmuons_name.c_str(), 10, 0., 10.);
    nmuons_->GetXaxis()->SetTitle("N_{mu}");
    nmuons_->GetYaxis()->SetTitle("Events");

    // njets
    const std::string njets_name = "N_{jet}";
    njets_ = tdir.make<TH1D>(njets_name.c_str(), njets_name.c_str(), 20, 0., 20.);
    njets_->GetXaxis()->SetTitle("N_{jet}");
    njets_->GetYaxis()->SetTitle("Events");

    // n_muonjets
    const std::string n_muonjets_name = "N_{muon_jet}";
    n_muonjets_ = tdir.make<TH1D>(n_muonjets_name.c_str(), n_muonjets_name.c_str(), 20, 0., 20.);
    n_muonjets_->GetXaxis()->SetTitle("N_{muon_jet}");
    n_muonjets_->GetYaxis()->SetTitle("Events");

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
      z_vtx_x_->Fill(zf_event.reco_z.vtx_x, EVENT_WEIGHT);
      z_vtx_y_->Fill(zf_event.reco_z.vtx_y, EVENT_WEIGHT);
      z_vtx_z_->Fill(zf_event.reco_z.vtx_z, EVENT_WEIGHT);
      z_vtx_prob_->Fill(zf_event.reco_z.vtx_prob, EVENT_WEIGHT);
      phistar_->Fill(zf_event.reco_z.phistar, EVENT_WEIGHT);

      primary_vtx_x_vs_z_vtx_x_->Fill(zf_event.reco_z.vtx_x, zf_event.reco_vert.primary_x );
      primary_vtx_y_vs_z_vtx_y_->Fill(zf_event.reco_z.vtx_y, zf_event.reco_vert.primary_y );
      primary_vtx_z_vs_z_vtx_z_->Fill(zf_event.reco_z.vtx_z, zf_event.reco_vert.primary_z );

      int n_jpsi = 0;
      for (unsigned int i = 0; i < zf_event.reco_jpsi.m.size() ; ++i ) {
        //TODO move to integer method (integer to represent cut level), use ! relative to ZFinder.cc
        if (APPLY_MUON_MIN_PT_ && (zf_event.mu0.at(i).pt() < MIN_MUON_PT || zf_event.mu1.at(i).pt() < MIN_MUON_PT) ) {
          continue;
        }
        if (APPLY_SOFT_MUONS_ && !(muon::isSoftMuon(zf_event.mu0.at(i), zf_event.reco_vert.primary_vert ) 
            && muon::isSoftMuon(zf_event.mu1.at(i), zf_event.reco_vert.primary_vert) ) ) {
          continue;
        }
        if (APPLY_JPSI_MASS_WINDOW_ && (zf_event.reco_jpsi.m.at(i) > MAX_JPSI_MASS || zf_event.reco_jpsi.m.at(i) < MIN_JPSI_MASS) ) {
          continue;
        }
        if (APPLY_VERTEX_Z_POS_WINDOW_ && fabs(zf_event.reco_vert.primary_z - zf_event.reco_jpsi.vtx_z.at(i)) > MAX_JPSI_VERTEX_Z_DISPLACEMENT ) {
          continue;
        }
        if (APPLY_DIMUON_VTX_COMPATIBILITY_ && !( zf_event.reco_jpsi.vtx_prob.at(i) >= MIN_VERTEX_PROB ) ) {
          continue;
        }
        n_jpsi++;

        jpsi_mass_all_->Fill(zf_event.reco_jpsi.m.at(i), EVENT_WEIGHT);
        jpsi_mass_coarse_->Fill(zf_event.reco_jpsi.m.at(i), EVENT_WEIGHT);
        jpsi_mass_fine_->Fill(zf_event.reco_jpsi.m.at(i), EVENT_WEIGHT);
        jpsi_rapidity_->Fill(zf_event.reco_jpsi.y.at(i), EVENT_WEIGHT);
        jpsi_pt_->Fill(zf_event.reco_jpsi.pt.at(i), EVENT_WEIGHT);
        jpsi_pt_vs_rap_->Fill(zf_event.reco_jpsi.y.at(i), zf_event.reco_jpsi.pt.at(i) , EVENT_WEIGHT);
        jpsi_vtx_distance_z_vtx_x_->Fill( zf_event.reco_jpsi.distance_x.at(i), EVENT_WEIGHT);
        jpsi_vtx_distance_z_vtx_y_->Fill( zf_event.reco_jpsi.distance_y.at(i), EVENT_WEIGHT);
        jpsi_vtx_distance_z_vtx_z_->Fill( zf_event.reco_jpsi.distance_z.at(i), EVENT_WEIGHT);
        jpsi_distance_->Fill(zf_event.reco_jpsi.distance.at(i), EVENT_WEIGHT);
        jpsi_dist_err_->Fill(zf_event.reco_jpsi.dist_err.at(i), EVENT_WEIGHT);
        jpsi_chi2_->Fill(zf_event.reco_jpsi.chi2.at(i), EVENT_WEIGHT);
        jpsi_distance_xy_->Fill(zf_event.reco_jpsi.distance_xy.at(i), EVENT_WEIGHT);
        jpsi_dist_err_xy_->Fill(zf_event.reco_jpsi.dist_err_xy.at(i), EVENT_WEIGHT);
        jpsi_chi2_xy_->Fill(zf_event.reco_jpsi.chi2_xy.at(i), EVENT_WEIGHT);
        jpsi_zpt_difference_->Fill( zf_event.reco_z.pt - zf_event.reco_jpsi.pt.at(i), EVENT_WEIGHT);
        jpsi_tau_xy_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        jpsi_tau_xy_fine_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        jpsi_tau_xy_very_fine_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps

        //pt
        if ( zf_event.reco_jpsi.pt.at(i) < 10.0 ) {
          jpsi_tau_xy_very_fine_ptUnder10_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( zf_event.reco_jpsi.pt.at(i) >= 10.0 && zf_event.reco_jpsi.pt.at(i) < 15 ) {
          jpsi_tau_xy_very_fine_pt10to15_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( zf_event.reco_jpsi.pt.at(i) >= 15.0 && zf_event.reco_jpsi.pt.at(i) < 20 ) {
          jpsi_tau_xy_very_fine_pt15to20_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( zf_event.reco_jpsi.pt.at(i) >= 20.0 && zf_event.reco_jpsi.pt.at(i) < 25 ) {
          jpsi_tau_xy_very_fine_pt20to25_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( zf_event.reco_jpsi.pt.at(i) >= 25.0 && zf_event.reco_jpsi.pt.at(i) < 30 ) {
          jpsi_tau_xy_very_fine_pt25to30_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( zf_event.reco_jpsi.pt.at(i) >= 20.0 ) {
          jpsi_tau_xy_very_fine_ptAbove20_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( zf_event.reco_jpsi.pt.at(i) >= 30.0 ) {
          jpsi_tau_xy_very_fine_ptAbove30_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }

        //rap
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 0.0 && fabs(zf_event.reco_jpsi.y.at(i)) < 0.3 ) {
          jpsi_tau_xy_very_fine_rap0_0to0_3_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 0.3 && fabs(zf_event.reco_jpsi.y.at(i)) < 0.6 ) {
          jpsi_tau_xy_very_fine_rap0_3to0_6_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 0.6 && fabs(zf_event.reco_jpsi.y.at(i)) < 0.9 ) {
          jpsi_tau_xy_very_fine_rap0_6to0_9_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 0.9 && fabs(zf_event.reco_jpsi.y.at(i)) < 1.2 ) {
          jpsi_tau_xy_very_fine_rap0_9to1_2_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 1.2 && fabs(zf_event.reco_jpsi.y.at(i)) < 1.5 ) {
          jpsi_tau_xy_very_fine_rap1_2to1_5_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 1.5 && fabs(zf_event.reco_jpsi.y.at(i)) < 1.8 ) {
          jpsi_tau_xy_very_fine_rap1_5to1_8_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 1.8 && fabs(zf_event.reco_jpsi.y.at(i)) < 2.1 ) {
          jpsi_tau_xy_very_fine_rap1_8to2_1_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 2.1 && fabs(zf_event.reco_jpsi.y.at(i)) < 2.4 ) {
          jpsi_tau_xy_very_fine_rap2_1to2_4_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }

        //rap and pt
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 0.0 && fabs(zf_event.reco_jpsi.y.at(i)) < 0.9 && (zf_event.reco_jpsi.pt.at(i) >= 10.0 && zf_event.reco_jpsi.pt.at(i) < 15 ) ) {
          jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 0.9 && fabs(zf_event.reco_jpsi.y.at(i)) < 1.2 && (zf_event.reco_jpsi.pt.at(i) >= 10.0 && zf_event.reco_jpsi.pt.at(i) < 15 ) ) {
          jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zf_event.reco_jpsi.y.at(i)) >= 1.2 && fabs(zf_event.reco_jpsi.y.at(i)) < 2.1 && (zf_event.reco_jpsi.pt.at(i) >= 10.0 && zf_event.reco_jpsi.pt.at(i) < 15 ) ) {
          jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }

        //similar pt muons
        //TODO do this in a better way, or maybe remove it
        if ( fabs( zf_event.mu0.at(i).pt() - zf_event.mu1.at(i).pt() ) < 2 ) {
          jpsi_tau_xy_very_fine_similar_pt_muons_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }

        //number of muon segments
        //TODO
        bool mu0ID = muon::isGoodMuon(zf_event.mu0.at(i), muon::TMOneStationTight);
        int mu0_layers = 0;
        if (mu0ID) {
          mu0_layers = zf_event.mu0.at(i).innerTrack()->hitPattern().trackerLayersWithMeasurement();
        }
        bool mu1ID = muon::isGoodMuon(zf_event.mu1.at(i), muon::TMOneStationTight);
        int mu1_layers = 0;
        if (mu1ID) {
          mu1_layers = zf_event.mu1.at(i).innerTrack()->hitPattern().trackerLayersWithMeasurement();
        }
        if ( mu0_layers > 12 && mu1_layers > 12 ) {
          jpsi_tau_xy_very_fine_above_12_tracker_layers_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        }

        jpsi_tau_z_->Fill(zf_event.reco_jpsi.tau_z.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        jpsi_tau_z_fine_->Fill(zf_event.reco_jpsi.tau_z.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        jpsi_tau_z_very_fine_->Fill(zf_event.reco_jpsi.tau_z.at(i) * 1000, EVENT_WEIGHT); // multiply by 1000 to go from ns to ps
        jpsi_mass_vs_chi2_->Fill(zf_event.reco_jpsi.m.at(i) , zf_event.reco_jpsi.chi2.at(i) );
        jpsi_tau_xy_vs_tau_z_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000 , zf_event.reco_jpsi.tau_z.at(i) * 1000 );
        jpsi_tau_xy_vs_distance_z_->Fill(zf_event.reco_jpsi.tau_xy.at(i) * 1000 , (zf_event.reco_jpsi.vtx_z.at(i) - zf_event.reco_z.vtx_z) );
        jpsi_tau_z_vs_distance_z_->Fill(zf_event.reco_jpsi.tau_z.at(i) * 1000 , (zf_event.reco_jpsi.vtx_z.at(i) - zf_event.reco_z.vtx_z) );

        dimuon_vtx_x_->Fill(zf_event.reco_jpsi.vtx_x.at(i), EVENT_WEIGHT);
        dimuon_vtx_y_->Fill(zf_event.reco_jpsi.vtx_y.at(i), EVENT_WEIGHT);
        dimuon_vtx_z_->Fill(zf_event.reco_jpsi.vtx_z.at(i), EVENT_WEIGHT);
        dimuon_vtx_prob_->Fill(zf_event.reco_jpsi.vtx_prob.at(i), EVENT_WEIGHT);

        z_jpsi_delta_phi_->Fill(zf_event.reco_jpsi.z_delta_phi.at(i), EVENT_WEIGHT);

        //TODO clean this up - deltaR vs delta_phi inconsistency in naming
        //TODO muon variables calculated in the ZFinderEvent part of the code?
        dimuon_delta_phi_->Fill(zf_event.reco_jpsi.muons_delta_phi.at(i), EVENT_WEIGHT);
        dimuon_deltaR_->Fill(zf_event.reco_jpsi.muons_deltaR.at(i), EVENT_WEIGHT);
        dimuon_delta_eta_->Fill(zf_event.reco_jpsi.muons_delta_eta.at(i), EVENT_WEIGHT);

        mu0_pt_->Fill(zf_event.mu0.at(i).pt(), EVENT_WEIGHT);
        mu0_eta_->Fill(zf_event.mu0.at(i).eta(), EVENT_WEIGHT);
        mu0_phi_->Fill(zf_event.mu0.at(i).phi(), EVENT_WEIGHT);
        mu0_charge_->Fill(zf_event.mu0.at(i).charge(), EVENT_WEIGHT);
        mu0_deltaR_truth_->Fill(zf_event.reco_jpsi_muon0.deltaR_truth.at(i), EVENT_WEIGHT);
        mu0_tracker_layers_->Fill(mu0_layers, EVENT_WEIGHT);

        mu1_pt_->Fill(zf_event.mu1.at(i).pt(), EVENT_WEIGHT);
        mu1_eta_->Fill(zf_event.mu1.at(i).eta(), EVENT_WEIGHT);
        mu1_phi_->Fill(zf_event.mu1.at(i).phi(), EVENT_WEIGHT);
        mu1_charge_->Fill(zf_event.mu1.at(i).charge(), EVENT_WEIGHT);
        mu1_deltaR_truth_->Fill(zf_event.reco_jpsi_muon1.deltaR_truth.at(i), EVENT_WEIGHT);
        mu1_tracker_layers_->Fill(mu1_layers, EVENT_WEIGHT);

        //TODO testing -------------------------------------
        int n_truth_matched_jpsi_muons = 0;
        if ( zf_event.reco_jpsi_muon0.deltaR_truth.at(i) <= MAX_DELTAR_TRUTH_MATCHED_JPSI_MUONS ) {
          n_truth_matched_jpsi_muons++;
        }
        if ( zf_event.reco_jpsi_muon1.deltaR_truth.at(i) <= MAX_DELTAR_TRUTH_MATCHED_JPSI_MUONS ) {
          n_truth_matched_jpsi_muons++;
        }
        n_truth_matched_jpsi_muons_->Fill (n_truth_matched_jpsi_muons, EVENT_WEIGHT);

        jpsi_iso_mu0_->Fill(zf_event.reco_jpsi.iso_mu0.at(i) ) ;
        jpsi_iso_sum_charged_hadron_pt_mu0_->Fill(zf_event.reco_jpsi.iso_sum_charged_hadron_pt_mu0.at(i) ) ;
        jpsi_iso_sum_charged_particle_pt_mu0_->Fill(zf_event.reco_jpsi.iso_sum_charged_particle_pt_mu0.at(i) ) ;
        jpsi_iso_sum_neutral_hadron_et_mu0_->Fill(zf_event.reco_jpsi.iso_sum_neutral_hadron_et_mu0.at(i) ) ;
        jpsi_iso_sum_photon_et_mu0_->Fill(zf_event.reco_jpsi.iso_sum_photon_et_mu0.at(i) ) ;
        jpsi_iso_sum_pileup_pt_mu0_->Fill(zf_event.reco_jpsi.iso_sum_pileup_pt_mu0.at(i) ) ;

        jpsi_iso_mu1_->Fill(zf_event.reco_jpsi.iso_mu1.at(i) ) ;
        jpsi_iso_sum_charged_hadron_pt_mu1_->Fill(zf_event.reco_jpsi.iso_sum_charged_hadron_pt_mu1.at(i) ) ;
        jpsi_iso_sum_charged_particle_pt_mu1_->Fill(zf_event.reco_jpsi.iso_sum_charged_particle_pt_mu1.at(i) ) ;
        jpsi_iso_sum_neutral_hadron_et_mu1_->Fill(zf_event.reco_jpsi.iso_sum_neutral_hadron_et_mu1.at(i) ) ;
        jpsi_iso_sum_photon_et_mu1_->Fill(zf_event.reco_jpsi.iso_sum_photon_et_mu1.at(i) ) ;
        jpsi_iso_sum_pileup_pt_mu1_->Fill(zf_event.reco_jpsi.iso_sum_pileup_pt_mu1.at(i) ) ;
      }
      for (unsigned int i = 0; i < zf_event.reco_jets.pt.size() ; ++i ) {
        jet_pt_->Fill(zf_event.reco_jets.pt.at(i), EVENT_WEIGHT);
        jet_eta_->Fill(zf_event.reco_jets.eta.at(i), EVENT_WEIGHT);
        jet_btag_discriminator_->Fill(zf_event.reco_jets.btag_discriminator.at(i), EVENT_WEIGHT);
      }
      for (unsigned int i = 0; i < zf_event.reco_muon_jets.pt.size() ; ++i ) {
        muon_jet_pt_->Fill(zf_event.reco_muon_jets.pt.at(i), EVENT_WEIGHT);
        muon_jet_pt_diff_z_pt_->Fill(zf_event.reco_z.pt - zf_event.reco_muon_jets.pt.at(i), EVENT_WEIGHT);
        muon_jet_pt_z_pt_->Fill(zf_event.reco_z.pt , zf_event.reco_muon_jets.pt.at(i), EVENT_WEIGHT);
        muon_jet_phi_z_phi_->Fill(zf_event.reco_z.phi , zf_event.reco_muon_jets.phi.at(i), EVENT_WEIGHT);
        for (unsigned int j = 0; j < zf_event.reco_jpsi.m.size() ; ++j ) {
          //TODO fix this, need more apply flags
          if (APPLY_MUON_MIN_PT_ && (zf_event.mu0.at(i).pt() < MIN_MUON_PT || zf_event.mu1.at(i).pt() < MIN_MUON_PT) ) {
            continue;
          }
          if (APPLY_SOFT_MUONS_ && !(muon::isSoftMuon(zf_event.mu0.at(i), zf_event.reco_vert.primary_vert ) 
                && muon::isSoftMuon(zf_event.mu1.at(i), zf_event.reco_vert.primary_vert) ) ) {
            continue;
          }
          if (APPLY_DIMUON_VTX_COMPATIBILITY_ && !( zf_event.reco_jpsi.vtx_prob.at(i) >= MIN_VERTEX_PROB ) ) {
            continue;
          }
          if (APPLY_VERTEX_Z_POS_WINDOW_ && fabs(zf_event.reco_vert.primary_z - zf_event.reco_jpsi.vtx_z.at(i)) > MAX_JPSI_VERTEX_Z_DISPLACEMENT ) {
            continue;
          }
          if (APPLY_JPSI_MASS_WINDOW_ && (zf_event.reco_jpsi.m.at(i) > MAX_JPSI_MASS || zf_event.reco_jpsi.m.at(i) < MIN_JPSI_MASS) ) {
            continue;
          }
          muon_jet_pt_diff_dimuon_pt_->Fill(zf_event.reco_jpsi.pt.at(j) - zf_event.reco_muon_jets.pt.at(i), EVENT_WEIGHT);
          muon_jet_pt_dimuon_pt_->Fill(zf_event.reco_jpsi.pt.at(j) , zf_event.reco_muon_jets.pt.at(i), EVENT_WEIGHT);
          muon_jet_phi_dimuon_phi_->Fill(zf_event.reco_jpsi.phi.at(j) , zf_event.reco_muon_jets.phi.at(i), EVENT_WEIGHT);
        }
        muon_jet_eta_->Fill(zf_event.reco_muon_jets.eta.at(i), EVENT_WEIGHT);
        muon_jet_btag_discriminator_->Fill(zf_event.reco_muon_jets.btag_discriminator.at(i), EVENT_WEIGHT);
      }
      for (unsigned int i = 0; i < zf_event.reco_vert.x.size() ; ++i ) {
        vtx_x_->Fill(zf_event.reco_vert.x.at(i), EVENT_WEIGHT);
        vtx_y_->Fill(zf_event.reco_vert.y.at(i), EVENT_WEIGHT);
        vtx_z_->Fill(zf_event.reco_vert.z.at(i), EVENT_WEIGHT);
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
      //}

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
      primary_vtx_x_->Fill(zf_event.reco_vert.primary_x, EVENT_WEIGHT);
      primary_vtx_y_->Fill(zf_event.reco_vert.primary_y, EVENT_WEIGHT);
      primary_vtx_z_->Fill(zf_event.reco_vert.primary_z, EVENT_WEIGHT);
      pileup_->Fill(zf_event.reco_vert.num, EVENT_WEIGHT);
      nelectrons_->Fill(zf_event.n_reco_electrons, EVENT_WEIGHT);
      nmuons_->Fill(zf_event.n_reco_muons, EVENT_WEIGHT);
      njets_->Fill(zf_event.n_reco_jets, EVENT_WEIGHT);
      n_muonjets_->Fill(zf_event.n_reco_muon_jets, EVENT_WEIGHT);
      njpsis_->Fill(n_jpsi, EVENT_WEIGHT); 
    } else if (USE_MC_ && !zf_event.is_real_data) {
      z0_mass_all_->Fill(zf_event.truth_z.m, zf_event.event_weight);
      z0_mass_coarse_->Fill(zf_event.truth_z.m, zf_event.event_weight);
      z0_mass_fine_->Fill(zf_event.truth_z.m, zf_event.event_weight);
      z0_rapidity_->Fill(zf_event.truth_z.y, zf_event.event_weight);
      z0_pt_->Fill(zf_event.truth_z.pt, zf_event.event_weight);
      phistar_->Fill(zf_event.truth_z.phistar, zf_event.event_weight);

      // Fill the histograms with the information from the approriate electron
      // TODO get rid of option to switch electrons, clutter and never used
      if (zf_event.e0_truth != NULL && zf_event.e1_truth != NULL){
        if (electron_0 == 0 && electron_1 == 1) {
          e0_pt_->Fill(zf_event.e0_truth->pt, zf_event.event_weight);
          e0_eta_->Fill(zf_event.e0_truth->eta, zf_event.event_weight);
          e0_phi_->Fill(zf_event.e0_truth->phi, zf_event.event_weight);
          e1_pt_->Fill(zf_event.e1_truth->pt, zf_event.event_weight);
          e1_eta_->Fill(zf_event.e1_truth->eta, zf_event.event_weight);
          e1_phi_->Fill(zf_event.e1_truth->phi, zf_event.event_weight);
        } else if (electron_0 == 1 && electron_1 == 0) {
          e0_pt_->Fill(zf_event.e1_truth->pt, zf_event.event_weight);
          e0_eta_->Fill(zf_event.e1_truth->eta, zf_event.event_weight);
          e0_phi_->Fill(zf_event.e1_truth->phi, zf_event.event_weight);
          e1_pt_->Fill(zf_event.e0_truth->pt, zf_event.event_weight);
          e1_eta_->Fill(zf_event.e0_truth->eta, zf_event.event_weight);
          e1_phi_->Fill(zf_event.e0_truth->phi, zf_event.event_weight);
        }
      }
      //jpsi plots
      int n_truth_jpsi = 0;
      for (unsigned int i = 0; i < zf_event.truth_jpsi.m.size() ; ++i ) {
        if (APPLY_JPSI_MASS_WINDOW_ && (zf_event.truth_jpsi.m.at(i) > MAX_JPSI_MASS || zf_event.truth_jpsi.m.at(i) < MIN_JPSI_MASS) ) {
          continue;
        }
        if (APPLY_VERTEX_Z_POS_WINDOW_ && fabs(zf_event.truth_jpsi.distance_z.at(i)) > MAX_JPSI_VERTEX_Z_DISPLACEMENT ) {
          continue;
        }
        if (APPLY_MUON_MIN_PT_ && (zf_event.jpsi_muon0.at(i)->pt() < MIN_MUON_PT || zf_event.jpsi_muon1.at(i)->pt() < MIN_MUON_PT) ) {
          continue;
        }
        
        n_truth_jpsi++;
        
        jpsi_mass_all_->Fill(zf_event.truth_jpsi.m.at(i), zf_event.event_weight);
        jpsi_mass_coarse_->Fill(zf_event.truth_jpsi.m.at(i), zf_event.event_weight);
        jpsi_mass_fine_->Fill(zf_event.truth_jpsi.m.at(i), zf_event.event_weight);
        jpsi_rapidity_->Fill(zf_event.truth_jpsi.y.at(i), zf_event.event_weight);
        jpsi_pt_->Fill(zf_event.truth_jpsi.pt.at(i), zf_event.event_weight);

        mu0_pt_->Fill(zf_event.jpsi_muon0.at(i)->pt(), zf_event.event_weight);
        mu0_eta_->Fill(zf_event.jpsi_muon0.at(i)->eta(), zf_event.event_weight);
        mu0_phi_->Fill(zf_event.jpsi_muon0.at(i)->phi(), zf_event.event_weight);
        mu0_charge_->Fill(zf_event.jpsi_muon0.at(i)->charge(), zf_event.event_weight);
        mu1_pt_->Fill(zf_event.jpsi_muon1.at(i)->pt(), zf_event.event_weight);
        mu1_eta_->Fill(zf_event.jpsi_muon1.at(i)->eta(), zf_event.event_weight);
        mu1_phi_->Fill(zf_event.jpsi_muon1.at(i)->phi(), zf_event.event_weight);
        mu1_charge_->Fill(zf_event.jpsi_muon1.at(i)->charge(), zf_event.event_weight);
      }
      // Event Info
      pileup_->Fill(zf_event.truth_vert.num, zf_event.event_weight);
      nelectrons_->Fill(2, zf_event.event_weight);  // We only ever grab the two electrons from the Z
      njpsis_->Fill(n_truth_jpsi, zf_event.event_weight); 
    }
  }

    void ZFinderPlotter::Print(const std::string& basename) {
      // Write all PNGs
      //std::string z0_mass_all_Str = basename + "_z0_mass_all" ;
      //TCanvas* z0_mass_all_C = new TCanvas(z0_mass_all_Str.c_str(), z0_mass_all_Str.c_str(), X_SIZE, Y_SIZE);
      //z0_mass_all_->Draw();
      //z0_mass_all_C->Print((z0_mass_all_Str+".png").c_str());

      //std::string z0_mass_coarse_Str = basename + "_z0_mass_coarse" ;
      //TCanvas* z0_mass_coarse_C = new TCanvas(z0_mass_coarse_Str.c_str(), z0_mass_coarse_Str.c_str(), X_SIZE, Y_SIZE);
      //z0_mass_coarse_->Draw();
      //z0_mass_coarse_C->Print((z0_mass_coarse_Str+".png").c_str());

      //std::string z0_mass_fine_Str = basename + "_z0_mass_fine";
      //TCanvas* z0_mass_fine_C = new TCanvas(z0_mass_fine_Str.c_str(), z0_mass_fine_Str.c_str(), X_SIZE, Y_SIZE);
      //z0_mass_fine_->Draw();
      //z0_mass_fine_C->Print((z0_mass_fine_Str+".png").c_str());

      //std::string z0_rapidity_Str = basename + "_z0_rapidity";
      //TCanvas* z0_rapidity_C = new TCanvas(z0_rapidity_Str.c_str(), z0_rapidity_Str.c_str(), X_SIZE, Y_SIZE);
      //z0_rapidity_->Draw();
      //z0_rapidity_C->Print((z0_rapidity_Str+".png").c_str());

      //std::string z0_ptStr = basename + "_z0_pt";
      //TCanvas* z0_ptC = new TCanvas(z0_ptStr.c_str(), z0_ptStr.c_str(), X_SIZE, Y_SIZE);
      //z0_pt_->Draw();
      //z0_ptC->Print((z0_ptStr+".png").c_str());

      //std::string e0_ptStr = basename + "_e0_pt";
      //TCanvas* e0_ptC = new TCanvas(e0_ptStr.c_str(), e0_ptStr.c_str(), X_SIZE, Y_SIZE);
      //e0_pt_->Draw();
      //e0_ptC->Print((e0_ptStr+".png").c_str());

      //std::string e1_ptStr = basename + "_e1_pt";
      //TCanvas* e1_ptC = new TCanvas(e1_ptStr.c_str(), e1_ptStr.c_str(), X_SIZE, Y_SIZE);
      //e1_pt_->Draw();
      //e1_ptC->Print((e1_ptStr+".png").c_str());

      //std::string e0_eta_Str = basename + "_e0_eta";
      //TCanvas* e0_eta_C = new TCanvas(e0_eta_Str.c_str(), e0_eta_Str.c_str(), X_SIZE, Y_SIZE);
      //e0_eta_->Draw();
      //e0_eta_C->Print((e0_eta_Str+".png").c_str());

      //std::string e1_eta_Str = basename + "_e1_eta";
      //TCanvas* e1_eta_C = new TCanvas(e1_eta_Str.c_str(), e1_eta_Str.c_str(), X_SIZE, Y_SIZE);
      //e1_eta_->Draw();
      //e1_eta_C->Print((e1_eta_Str+".png").c_str());

      //std::string e0_phi_Str = basename + "_e0_phi";
      //TCanvas* e0_phi_C = new TCanvas(e0_phi_Str.c_str(), e0_phi_Str.c_str(), X_SIZE, Y_SIZE);
      //e0_phi_->Draw();
      //e0_phi_C->Print((e0_phi_Str+".png").c_str());

      //std::string e1_phi_Str = basename + "_e1_phi";
      //TCanvas* e1_phi_C = new TCanvas(e1_phi_Str.c_str(), e1_phi_Str.c_str(), X_SIZE, Y_SIZE);
      //e1_phi_->Draw();
      //e1_phi_C->Print((e1_phi_Str+".png").c_str());

      //std::string e0_charge_Str = basename + "_e0_charge";
      //TCanvas* e0_charge_C = new TCanvas(e0_charge_Str.c_str(), e0_charge_Str.c_str(), X_SIZE, Y_SIZE);
      //e0_charge_->Draw();
      //e0_charge_C->Print((e0_charge_Str+".png").c_str());

      //std::string e1_charge_Str = basename + "_e1_charge";
      //TCanvas* e1_charge_C = new TCanvas(e1_charge_Str.c_str(), e1_charge_Str.c_str(), X_SIZE, Y_SIZE);
      //e1_charge_->Draw();
      //e1_charge_C->Print((e1_charge_Str+".png").c_str());

      //std::string phistarStr = basename + "_phistar";
      //TCanvas* phistarC = new TCanvas(phistarStr.c_str(), phistarStr.c_str(), X_SIZE, Y_SIZE);
      //phistar_->Draw();
      //phistarC->Print((phistarStr+".png").c_str());

      //std::string nelectronsStr = basename + "_nelectrons";
      //TCanvas* nelectronsC = new TCanvas(nelectronsStr.c_str(), nelectronsStr.c_str(), X_SIZE, Y_SIZE);
      //nelectrons_->Draw();
      //nelectronsC->Print((nelectronsStr+".png").c_str());

      //std::string jpsi_mass_all_Str = basename + "_jpsi_mass_all" ;
      //TCanvas* jpsi_mass_all_C = new TCanvas(jpsi_mass_all_Str.c_str(), jpsi_mass_all_Str.c_str(), X_SIZE, Y_SIZE);
      //jpsi_mass_all_->Draw();
      //jpsi_mass_all_C->Print((jpsi_mass_all_Str+".png").c_str());

      //std::string jpsi_mass_coarse_Str = basename + "_jpsi_mass_coarse" ;
      //TCanvas* jpsi_mass_coarse_C = new TCanvas(jpsi_mass_coarse_Str.c_str(), jpsi_mass_coarse_Str.c_str(), X_SIZE, Y_SIZE);
      //jpsi_mass_coarse_->Draw();
      //jpsi_mass_coarse_C->Print((jpsi_mass_coarse_Str+".png").c_str());

      //std::string jpsi_mass_fine_Str = basename + "_jpsi_mass_fine";
      //TCanvas* jpsi_mass_fine_C = new TCanvas(jpsi_mass_fine_Str.c_str(), jpsi_mass_fine_Str.c_str(), X_SIZE, Y_SIZE);
      //jpsi_mass_fine_->Draw();
      //jpsi_mass_fine_C->Print((jpsi_mass_fine_Str+".png").c_str());

      //std::string jpsi_rapidity_Str = basename + "_jpsi_rapidity";
      //TCanvas* jpsi_rapidity_C = new TCanvas(jpsi_rapidity_Str.c_str(), jpsi_rapidity_Str.c_str(), X_SIZE, Y_SIZE);
      //jpsi_rapidity_->Draw();
      //jpsi_rapidity_C->Print((jpsi_rapidity_Str+".png").c_str());

      //std::string jpsi_ptStr = basename + "_jpsi_pt";
      //TCanvas* jpsi_ptC = new TCanvas(jpsi_ptStr.c_str(), jpsi_ptStr.c_str(), X_SIZE, Y_SIZE);
      //jpsi_pt_->Draw();
      //jpsi_ptC->Print((jpsi_ptStr+".png").c_str());

      //std::string mu0_ptStr = basename + "_mu0_pt";
      //TCanvas* mu0_ptC = new TCanvas(mu0_ptStr.c_str(), mu0_ptStr.c_str(), X_SIZE, Y_SIZE);
      //mu0_pt_->Draw();
      //mu0_ptC->Print((mu0_ptStr+".png").c_str());

      //std::string mu1_ptStr = basename + "_mu1_pt";
      //TCanvas* mu1_ptC = new TCanvas(mu1_ptStr.c_str(), mu1_ptStr.c_str(), X_SIZE, Y_SIZE);
      //mu1_pt_->Draw();
      //mu1_ptC->Print((mu1_ptStr+".png").c_str());

      //std::string mu0_eta_Str = basename + "_mu0_eta";
      //TCanvas* mu0_eta_C = new TCanvas(mu0_eta_Str.c_str(), mu0_eta_Str.c_str(), X_SIZE, Y_SIZE);
      //mu0_eta_->Draw();
      //mu0_eta_C->Print((mu0_eta_Str+".png").c_str());

      //std::string mu1_eta_Str = basename + "_mu1_eta";
      //TCanvas* mu1_eta_C = new TCanvas(mu1_eta_Str.c_str(), mu1_eta_Str.c_str(), X_SIZE, Y_SIZE);
      //mu1_eta_->Draw();
      //mu1_eta_C->Print((mu1_eta_Str+".png").c_str());

      //std::string mu0_phi_Str = basename + "_mu0_phi";
      //TCanvas* mu0_phi_C = new TCanvas(mu0_phi_Str.c_str(), mu0_phi_Str.c_str(), X_SIZE, Y_SIZE);
      //mu0_phi_->Draw();
      //mu0_phi_C->Print((mu0_phi_Str+".png").c_str());

      //std::string mu1_phi_Str = basename + "_mu1_phi";
      //TCanvas* mu1_phi_C = new TCanvas(mu1_phi_Str.c_str(), mu1_phi_Str.c_str(), X_SIZE, Y_SIZE);
      //mu1_phi_->Draw();
      //mu1_phi_C->Print((mu1_phi_Str+".png").c_str());

      //std::string mu0_charge_Str = basename + "_mu0_charge";
      //TCanvas* mu0_charge_C = new TCanvas(mu0_charge_Str.c_str(), mu0_charge_Str.c_str(), X_SIZE, Y_SIZE);
      //mu0_charge_->Draw();
      //mu0_charge_C->Print((mu0_charge_Str+".png").c_str());

      //std::string mu1_charge_Str = basename + "_mu1_charge";
      //TCanvas* mu1_charge_C = new TCanvas(mu1_charge_Str.c_str(), mu1_charge_Str.c_str(), X_SIZE, Y_SIZE);
      //mu1_charge_->Draw();
      //mu1_charge_C->Print((mu1_charge_Str+".png").c_str());

      //std::string nmuonsStr = basename + "_nmuons";
      //TCanvas* nmuonsC = new TCanvas(nmuonsStr.c_str(), nmuonsStr.c_str(), X_SIZE, Y_SIZE);
      //nmuons_->Draw();
      //nmuonsC->Print((nmuonsStr+".png").c_str());

      //std::string njpsisStr = basename + "_njpsis";
      //TCanvas* njpsisC = new TCanvas(njpsisStr.c_str(), njpsisStr.c_str(), X_SIZE, Y_SIZE);
      //njpsis_->Draw();
      //njpsisC->Print((njpsisStr+".png").c_str());

      //std::string pileupStr = basename + "_pileup";
      //TCanvas* pileupC = new TCanvas(pileupStr.c_str(), pileupStr.c_str(), X_SIZE, Y_SIZE);
      //pileup_->Draw();
      //pileupC->Print((pileupStr+".png").c_str());
    }
  }  // namespace zf
