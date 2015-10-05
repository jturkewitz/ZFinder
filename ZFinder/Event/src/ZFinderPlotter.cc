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
      const bool APPLY_JPSI_MASS_WINDOW , const bool APPLY_VERTEX_Z_POS_WINDOW, const bool APPLY_PROMPT_JPSI_WINDOW)
    : USE_MC_(USE_MC), APPLY_MUON_MIN_PT_(APPLY_MUON_MIN_PT), APPLY_SOFT_MUONS_(APPLY_SOFT_MUONS), APPLY_DIMUON_VTX_COMPATIBILITY_(APPLY_DIMUON_VTX_COMPATIBILITY),
    APPLY_JPSI_MASS_WINDOW_(APPLY_JPSI_MASS_WINDOW), APPLY_VERTEX_Z_POS_WINDOW_(APPLY_VERTEX_Z_POS_WINDOW), APPLY_PROMPT_JPSI_WINDOW_(APPLY_PROMPT_JPSI_WINDOW) {
      /*
       * Initialize a set of histograms and associate them with a given TDirectory.
       */

      // Set up histograms
      // z_mass_all_
      const std::string z_mass_all_name = "z Mass: All";
      z_mass_all_ = tdir.make<TH1D>(z_mass_all_name.c_str(), z_mass_all_name.c_str(), 300, 0., 300.);
      z_mass_all_->GetXaxis()->SetTitle("m_{ee} [GeV]");
      z_mass_all_->GetYaxis()->SetTitle("Counts / GeV");

      // z_mass_coarse_
      const std::string z_mass_coarse_name = "z Mass: Coarse";
      z_mass_coarse_ = tdir.make<TH1D>(z_mass_coarse_name.c_str(), z_mass_coarse_name.c_str(), 100, 50., 150.);
      z_mass_coarse_->GetXaxis()->SetTitle("m_{ee} [GeV]");
      z_mass_coarse_->GetYaxis()->SetTitle("Counts / GeV");

      // z_mass_fine_
      const std::string z_mass_fine_name = "z Mass: Fine";
      z_mass_fine_ = tdir.make<TH1D>(z_mass_fine_name.c_str(), z_mass_fine_name.c_str(), 80, 80., 100.);
      z_mass_fine_->GetXaxis()->SetTitle("m_{ee} [GeV]");
      z_mass_fine_->GetYaxis()->SetTitle("Counts / 0.25 GeV");

      // z_rapidity_
      const std::string z_rapidity_name = "z Rapidity";
      z_rapidity_ = tdir.make<TH1D>(z_rapidity_name.c_str(), z_rapidity_name.c_str(), 100, -5., 5.);
      z_rapidity_->GetXaxis()->SetTitle("Z_{Y}");
      z_rapidity_->GetYaxis()->SetTitle("Counts");

      // z_pt
      const std::string z_pt_name = "z p_{T}";
      z_pt_ = tdir.make<TH1D>(z_pt_name.c_str(), z_pt_name.c_str(), 200, 0., 200.);
      z_pt_->GetXaxis()->SetTitle("p_{T,Z}");
      z_pt_->GetYaxis()->SetTitle("Counts / GeV");

      // z_vtx_prob_
      const std::string z_vtx_prob_name = "dielectron vertex probability";
      z_vtx_prob_ = tdir.make<TH1D>(z_vtx_prob_name.c_str(), z_vtx_prob_name.c_str(), 200, 0., 1.);
      z_vtx_prob_->GetXaxis()->SetTitle("Probability");
      z_vtx_prob_->GetYaxis()->SetTitle("Probability / 0.005 ");

      // phistar
      const std::string phistar_name = "#phi*";
      phistar_ = tdir.make<TH1D>(phistar_name.c_str(), phistar_name.c_str(), 100, 0., 1.);
      phistar_->GetXaxis()->SetTitle("#phi*");
      phistar_->GetYaxis()->SetTitle("Counts");

      // z_from_muons_mass_all_
      const std::string z_from_muons_mass_all_name = "Z From Muons Mass: All";
      z_from_muons_mass_all_ = tdir.make<TH1D>(z_from_muons_mass_all_name.c_str(), z_from_muons_mass_all_name.c_str(), 300, 0., 300.);
      z_from_muons_mass_all_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      z_from_muons_mass_all_->GetYaxis()->SetTitle("Counts / GeV");

      // z_from_muons_mass_coarse_
      const std::string z_from_muons_mass_coarse_name = "Z From Muons Mass: Coarse";
      z_from_muons_mass_coarse_ = tdir.make<TH1D>(z_from_muons_mass_coarse_name.c_str(), z_from_muons_mass_coarse_name.c_str(), 100, 50., 150.);
      z_from_muons_mass_coarse_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      z_from_muons_mass_coarse_->GetYaxis()->SetTitle("Counts / GeV");

      // z_from_muons_mass_fine_
      const std::string z_from_muons_mass_fine_name = "Z From Muons Mass: Fine";
      z_from_muons_mass_fine_ = tdir.make<TH1D>(z_from_muons_mass_fine_name.c_str(), z_from_muons_mass_fine_name.c_str(), 80, 80., 100.);
      z_from_muons_mass_fine_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      z_from_muons_mass_fine_->GetYaxis()->SetTitle("Counts / 0.25 GeV");

      // z_from_muons_rapidity_
      const std::string z_from_muons_rapidity_name = "Z From Muons Rapidity";
      z_from_muons_rapidity_ = tdir.make<TH1D>(z_from_muons_rapidity_name.c_str(), z_from_muons_rapidity_name.c_str(), 100, -5., 5.);
      z_from_muons_rapidity_->GetXaxis()->SetTitle("Z_{Y}");
      z_from_muons_rapidity_->GetYaxis()->SetTitle("Counts");

      // z_from_muons_pt
      const std::string z_from_muons_pt_name = "Z From Muons p_{T}";
      z_from_muons_pt_ = tdir.make<TH1D>(z_from_muons_pt_name.c_str(), z_from_muons_pt_name.c_str(), 200, 0., 200.);
      z_from_muons_pt_->GetXaxis()->SetTitle("p_{T,Z}");
      z_from_muons_pt_->GetYaxis()->SetTitle("Counts / GeV");

      // z_from_muons_vtx_prob_
      const std::string z_from_muons_vtx_prob_name = "Z From Muons dimuon vertex probability";
      z_from_muons_vtx_prob_ = tdir.make<TH1D>(z_from_muons_vtx_prob_name.c_str(), z_from_muons_vtx_prob_name.c_str(), 200, 0., 1.);
      z_from_muons_vtx_prob_->GetXaxis()->SetTitle("Probability");
      z_from_muons_vtx_prob_->GetYaxis()->SetTitle("Probability / 0.005 ");

      // z_from_muons_phistar
      const std::string z_from_muons_phistar_name = "Z From Muons #phi*";
      z_from_muons_phistar_ = tdir.make<TH1D>(z_from_muons_phistar_name.c_str(), z_from_muons_phistar_name.c_str(), 100, 0., 1.);
      z_from_muons_phistar_->GetXaxis()->SetTitle("#phi*");
      z_from_muons_phistar_->GetYaxis()->SetTitle("Counts");

      // muon0_from_z_pt
      const std::string muon0_from_z_pt_name = "muon0_from_z_pt";
      muon0_from_z_pt_ = tdir.make<TH1D>(muon0_from_z_pt_name.c_str(), muon0_from_z_pt_name.c_str(), 200, 0., 200.);
      muon0_from_z_pt_->GetXaxis()->SetTitle("p_{T,mu_{0}}");
      muon0_from_z_pt_->GetYaxis()->SetTitle("Counts / GeV");

      // muon1_from_z_pt
      const std::string muon1_from_z_pt_name = "muon1_from_z_pt";
      muon1_from_z_pt_ = tdir.make<TH1D>(muon1_from_z_pt_name.c_str(), muon1_from_z_pt_name.c_str(), 200, 0., 200.);
      muon1_from_z_pt_->GetXaxis()->SetTitle("p_{T,mu_{0}}");
      muon1_from_z_pt_->GetYaxis()->SetTitle("Counts / GeV");

      // muon0_from_z_eta_
      const std::string muon0_from_z_eta_name = "muon0_from_z_eta";
      muon0_from_z_eta_ = tdir.make<TH1D>(muon0_from_z_eta_name.c_str(), muon0_from_z_eta_name.c_str(), 50, -5., 5.);
      muon0_from_z_eta_->GetXaxis()->SetTitle("#eta_{mu_{0}}");
      muon0_from_z_eta_->GetYaxis()->SetTitle("Counts");

      // muon1_from_z_eta_
      const std::string muon1_from_z_eta_name = "muon1_from_z_eta";
      muon1_from_z_eta_ = tdir.make<TH1D>(muon1_from_z_eta_name.c_str(), muon1_from_z_eta_name.c_str(), 50, -5., 5.);
      muon1_from_z_eta_->GetXaxis()->SetTitle("#eta_{mu_{1}}");
      muon1_from_z_eta_->GetYaxis()->SetTitle("Counts");

      // muon0_from_z_phi_
      const std::string muon0_from_z_phi_name = "muon0_from_z_phi";
      muon0_from_z_phi_ = tdir.make<TH1D>(muon0_from_z_phi_name.c_str(), muon0_from_z_phi_name.c_str(), 60, -3.15, 3.15);
      muon0_from_z_phi_->GetXaxis()->SetTitle("#phi_{mu_{0}}");
      muon0_from_z_phi_->GetYaxis()->SetTitle("Counts");

      // muon1_from_z_phi_
      const std::string muon1_from_z_phi_name = "muon1_from_z_phi";
      muon1_from_z_phi_ = tdir.make<TH1D>(muon1_from_z_phi_name.c_str(), muon1_from_z_phi_name.c_str(), 50, -3.15, 3.15);
      muon1_from_z_phi_->GetXaxis()->SetTitle("#phi_{mu_{1}}");
      muon1_from_z_phi_->GetYaxis()->SetTitle("counts");

      // muon0_from_z_charge_
      const std::string muon0_from_z_charge_name = "muon0_from_z_charge";
      muon0_from_z_charge_ = tdir.make<TH1D>(muon0_from_z_charge_name.c_str(), muon0_from_z_charge_name.c_str(), 60, -3.15, 3.15);
      muon0_from_z_charge_->GetXaxis()->SetTitle("charge_{mu_{0}}");
      muon0_from_z_charge_->GetYaxis()->SetTitle("Counts");

      // muon1_from_z_charge_
      const std::string muon1_from_z_charge_name = "muon1_from_z_charge";
      muon1_from_z_charge_ = tdir.make<TH1D>(muon1_from_z_charge_name.c_str(), muon1_from_z_charge_name.c_str(), 60, -3.15, 3.15);
      muon1_from_z_charge_->GetXaxis()->SetTitle("charge_{mu_{1}}");
      muon1_from_z_charge_->GetYaxis()->SetTitle("counts");

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
      e0_phi_ = tdir.make<TH1D>(e0_phi_name.c_str(), e0_phi_name.c_str(), 63, -3.15, 3.15);
      e0_phi_->GetXaxis()->SetTitle("#phi_{e_{0}}");
      e0_phi_->GetYaxis()->SetTitle("Counts / 0.1 Rad");

      // e1_phi_
      const std::string e1_phi_name = "#phi_{e_{1}}";
      e1_phi_ = tdir.make<TH1D>(e1_phi_name.c_str(), e1_phi_name.c_str(), 63, -3.15, 3.15);
      e1_phi_->GetXaxis()->SetTitle("#phi_{e_{1}}");
      e1_phi_->GetYaxis()->SetTitle("Counts / 0.1 Rad");

      // e0_charge_
      const std::string e0_charge_name = "charge_{e_{0}}";
      e0_charge_ = tdir.make<TH1D>(e0_charge_name.c_str(), e0_charge_name.c_str(), 3, -1, 2);
      e0_charge_->GetXaxis()->SetTitle("charge_{e_{0}}");
      e0_charge_->GetYaxis()->SetTitle("Counts");

      // e1_charge_
      const std::string e1_charge_name = "charge_{e_{1}}";
      e1_charge_ = tdir.make<TH1D>(e1_charge_name.c_str(), e1_charge_name.c_str(), 3, -1, 2);
      e1_charge_->GetXaxis()->SetTitle("charge_{e_{1}}");
      e1_charge_->GetYaxis()->SetTitle("counts");

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
      const std::string jpsi_mass_fine_name = "jpsi_mass";
      jpsi_mass_fine_ = tdir.make<TH1D>(jpsi_mass_fine_name.c_str(), jpsi_mass_fine_name.c_str(), 1500, 0.0, 15.0);
      jpsi_mass_fine_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      jpsi_mass_fine_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

      // jpsi_four_lepton_mass_
      const std::string jpsi_four_lepton_mass_name = "four_lepton_mass";
      jpsi_four_lepton_mass_ = tdir.make<TH1D>(jpsi_four_lepton_mass_name.c_str(), jpsi_four_lepton_mass_name.c_str(), 300, 50.0, 200.0);
      jpsi_four_lepton_mass_->GetXaxis()->SetTitle("m_{llll} [GeV]");
      jpsi_four_lepton_mass_->GetYaxis()->SetTitle("Counts / 0.5 GeV");

      // jpsi_mass_fine_ptUnder10_
      const std::string jpsi_mass_fine_ptUnder10_name = "jpsi_mass_ptUnder10";
      jpsi_mass_fine_ptUnder10_ = tdir.make<TH1D>(jpsi_mass_fine_ptUnder10_name.c_str(), jpsi_mass_fine_ptUnder10_name.c_str(), 250, 2.5, 5.0);
      jpsi_mass_fine_ptUnder10_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      jpsi_mass_fine_ptUnder10_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

      // jpsi_mass_fine_pt10to15_
      const std::string jpsi_mass_fine_pt10to15_name = "jpsi_mass_pt10to15";
      jpsi_mass_fine_pt10to15_ = tdir.make<TH1D>(jpsi_mass_fine_pt10to15_name.c_str(), jpsi_mass_fine_pt10to15_name.c_str(), 250, 2.5, 5.0);
      jpsi_mass_fine_pt10to15_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      jpsi_mass_fine_pt10to15_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

      // jpsi_mass_fine_pt15to20_
      const std::string jpsi_mass_fine_pt15to20_name = "jpsi_mass_pt15to20";
      jpsi_mass_fine_pt15to20_ = tdir.make<TH1D>(jpsi_mass_fine_pt15to20_name.c_str(), jpsi_mass_fine_pt15to20_name.c_str(), 250, 2.5, 5.0);
      jpsi_mass_fine_pt15to20_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      jpsi_mass_fine_pt15to20_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

      // jpsi_mass_fine_pt20to25_
      const std::string jpsi_mass_fine_pt20to25_name = "jpsi_mass_pt20to25";
      jpsi_mass_fine_pt20to25_ = tdir.make<TH1D>(jpsi_mass_fine_pt20to25_name.c_str(), jpsi_mass_fine_pt20to25_name.c_str(), 250, 2.5, 5.0);
      jpsi_mass_fine_pt20to25_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      jpsi_mass_fine_pt20to25_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

      // jpsi_mass_fine_pt25to30_
      const std::string jpsi_mass_fine_pt25to30_name = "jpsi_mass_pt25to30";
      jpsi_mass_fine_pt25to30_ = tdir.make<TH1D>(jpsi_mass_fine_pt25to30_name.c_str(), jpsi_mass_fine_pt25to30_name.c_str(), 250, 2.5, 5.0);
      jpsi_mass_fine_pt25to30_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      jpsi_mass_fine_pt25to30_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

      // jpsi_mass_fine_ptAbove30_
      const std::string jpsi_mass_fine_ptAbove30_name = "jpsi_mass_ptAbove30";
      jpsi_mass_fine_ptAbove30_ = tdir.make<TH1D>(jpsi_mass_fine_ptAbove30_name.c_str(), jpsi_mass_fine_ptAbove30_name.c_str(), 250, 2.5, 5.0);
      jpsi_mass_fine_ptAbove30_->GetXaxis()->SetTitle("m_{mumu} [GeV]");
      jpsi_mass_fine_ptAbove30_->GetYaxis()->SetTitle("Counts / 0.01 GeV");

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

      // jpsi_efficiency
      const std::string jpsi_efficiency_name = "jpsi efficiency";
      jpsi_efficiency_ = tdir.make<TH1D>(jpsi_efficiency_name.c_str(), jpsi_efficiency_name.c_str(), 100, 0., 1.);
      jpsi_efficiency_->GetXaxis()->SetTitle("Jpsi Efficiency");
      jpsi_efficiency_->GetYaxis()->SetTitle("Counts / 0.01");

      // jpsi_acc_eff
      const std::string jpsi_acc_eff_name = "jpsi_acc_eff";
      jpsi_acc_eff_ = tdir.make<TH1D>(jpsi_acc_eff_name.c_str(), jpsi_acc_eff_name.c_str(), 100, 0., 1.);
      jpsi_acc_eff_->GetXaxis()->SetTitle("Jpsi Acc*Eff");
      jpsi_acc_eff_->GetYaxis()->SetTitle("Counts / 0.01");

      // jpsi_scale_factor
      const std::string jpsi_scale_factor_name = "jpsi scale factor";
      jpsi_scale_factor_ = tdir.make<TH1D>(jpsi_scale_factor_name.c_str(), jpsi_scale_factor_name.c_str(), 200, 0., 2.);
      jpsi_scale_factor_->GetXaxis()->SetTitle("Jpsi Scale factor");
      jpsi_scale_factor_->GetYaxis()->SetTitle("Counts / 0.01");

      // jpsi_reco_pt_vs_jpsi_truth_pt
      const std::string jpsi_reco_pt_vs_jpsi_truth_pt_name = "jpsi_reco_pt_vs_jpsi_truth_pt";
      jpsi_reco_pt_vs_jpsi_truth_pt_ = tdir.make<TH2D>(jpsi_reco_pt_vs_jpsi_truth_pt_name.c_str(), jpsi_reco_pt_vs_jpsi_truth_pt_name.c_str(), 100, 0., 100., 100, 0.0, 100);
      jpsi_reco_pt_vs_jpsi_truth_pt_->GetXaxis()->SetTitle("truth pT [GeV]");
      jpsi_reco_pt_vs_jpsi_truth_pt_->GetYaxis()->SetTitle("reco pT [GeV]");

      // jpsi_truth_pt_minus_jpsi_reco_pt
      const std::string jpsi_truth_pt_minus_jpsi_reco_pt_name = "jpsi_truth_pt_minus_jpsi_reco_pt";
      jpsi_truth_pt_minus_jpsi_reco_pt_ = tdir.make<TH1D>(jpsi_truth_pt_minus_jpsi_reco_pt_name.c_str(), jpsi_truth_pt_minus_jpsi_reco_pt_name.c_str(), 200, -10., 10.);
      jpsi_truth_pt_minus_jpsi_reco_pt_->GetXaxis()->SetTitle("truth pT - reco pT [GeV]");
      jpsi_truth_pt_minus_jpsi_reco_pt_->GetYaxis()->SetTitle("Counts / 0.1 GeV");

      // jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_
      const std::string jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_name = "jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x";
      jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_ = tdir.make<TH1D>(jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_name.c_str(), jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_name.c_str(), 2000, -10., 10.);
      jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_->GetXaxis()->SetTitle("truth vtx_x - reco vtx_y (cm)");
      jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_->GetYaxis()->SetTitle("Counts / 0.01 cm");

      // jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_
      const std::string jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_name = "jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y";
      jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_ = tdir.make<TH1D>(jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_name.c_str(), jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_name.c_str(), 2000, -10., 10.);
      jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_->GetXaxis()->SetTitle("truth vtx_y - reco vtx_y (cm)");
      jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_->GetYaxis()->SetTitle("Counts / 0.01 cm");

      // jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_
      const std::string jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_name = "jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_";
      jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_ = tdir.make<TH1D>(jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_name.c_str(), jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_name.c_str(), 2000, -10., 10.);
      jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_->GetXaxis()->SetTitle("truth vtx_z - reco vtx_z (cm)");
      jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_->GetYaxis()->SetTitle("Counts / 0.01 cm");

      // jpsi_trigger_obj_mu0_pt
      const std::string jpsi_trigger_obj_mu0_pt_name = "jpsi_trigger_obj_mu0_pt";
      jpsi_trigger_obj_mu0_pt_ = tdir.make<TH1D>(jpsi_trigger_obj_mu0_pt_name.c_str(), jpsi_trigger_obj_mu0_pt_name.c_str(), 1000, 0., 100.);
      jpsi_trigger_obj_mu0_pt_->GetXaxis()->SetTitle("jpsi trigger object mu0 pT [GeV]");
      jpsi_trigger_obj_mu0_pt_->GetYaxis()->SetTitle("Counts / 0.1 GeV");

      // jpsi_trigger_obj_mu1_pt
      const std::string jpsi_trigger_obj_mu1_pt_name = "jpsi_trigger_obj_mu1_pt";
      jpsi_trigger_obj_mu1_pt_ = tdir.make<TH1D>(jpsi_trigger_obj_mu1_pt_name.c_str(), jpsi_trigger_obj_mu1_pt_name.c_str(), 1000, 0., 100.);
      jpsi_trigger_obj_mu1_pt_->GetXaxis()->SetTitle("jpsi trigger object mu1 pT [GeV]");
      jpsi_trigger_obj_mu1_pt_->GetYaxis()->SetTitle("Counts / 0.1 GeV");

      // jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt
      const std::string jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt_name = "jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt";
      jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt_ = tdir.make<TH1D>(jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt_name.c_str(), jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt_name.c_str(), 1000, -10., 10.);
      jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt_->GetXaxis()->SetTitle("trigger mu0 difference reco mu0 pT [GeV]");
      jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt_->GetYaxis()->SetTitle("Counts / 0.1 GeV");

      // jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt
      const std::string jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt_name = "jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt";
      jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt_ = tdir.make<TH1D>(jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt_name.c_str(), jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt_name.c_str(), 1000, -10., 10.);
      jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt_->GetXaxis()->SetTitle("trigger mu1 difference reco mu1 pT [GeV]");
      jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt_->GetYaxis()->SetTitle("Counts / 0.02 GeV");

      // jpsi_cos_mu_plus
      const std::string jpsi_cos_mu_plus_name = "jpsi_cos_mu_plus";
      jpsi_cos_mu_plus_ = tdir.make<TH1D>(jpsi_cos_mu_plus_name.c_str(), jpsi_cos_mu_plus_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_->GetXaxis()->SetTitle("jpsi mu_plus cos(theta)");
      jpsi_cos_mu_plus_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_jpsi_pt_8to8p5
      const std::string jpsi_cos_mu_plus_jpsi_pt_8to8p5_name = "jpsi_cos_mu_plus_jpsi_pt_8to8p5";
      jpsi_cos_mu_plus_jpsi_pt_8to8p5_ = tdir.make<TH1D>(jpsi_cos_mu_plus_jpsi_pt_8to8p5_name.c_str(), jpsi_cos_mu_plus_jpsi_pt_8to8p5_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_jpsi_pt_8to8p5_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_8to8p5 cos(theta)");
      jpsi_cos_mu_plus_jpsi_pt_8to8p5_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_jpsi_pt_8p5to9
      const std::string jpsi_cos_mu_plus_jpsi_pt_8p5to9_name = "jpsi_cos_mu_plus_jpsi_pt_8p5to9";
      jpsi_cos_mu_plus_jpsi_pt_8p5to9_ = tdir.make<TH1D>(jpsi_cos_mu_plus_jpsi_pt_8p5to9_name.c_str(), jpsi_cos_mu_plus_jpsi_pt_8p5to9_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_jpsi_pt_8p5to9_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_8p5to9 cos(theta)");
      jpsi_cos_mu_plus_jpsi_pt_8p5to9_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_jpsi_pt_9to10
      const std::string jpsi_cos_mu_plus_jpsi_pt_9to10_name = "jpsi_cos_mu_plus_jpsi_pt_9to10";
      jpsi_cos_mu_plus_jpsi_pt_9to10_ = tdir.make<TH1D>(jpsi_cos_mu_plus_jpsi_pt_9to10_name.c_str(), jpsi_cos_mu_plus_jpsi_pt_9to10_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_jpsi_pt_9to10_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_9to10 cos(theta)");
      jpsi_cos_mu_plus_jpsi_pt_9to10_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_jpsi_pt_10to15
      const std::string jpsi_cos_mu_plus_jpsi_pt_10to15_name = "jpsi_cos_mu_plus_jpsi_pt_10to15";
      jpsi_cos_mu_plus_jpsi_pt_10to15_ = tdir.make<TH1D>(jpsi_cos_mu_plus_jpsi_pt_10to15_name.c_str(), jpsi_cos_mu_plus_jpsi_pt_10to15_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_jpsi_pt_10to15_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_10to15 cos(theta)");
      jpsi_cos_mu_plus_jpsi_pt_10to15_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_jpsi_pt_15to20
      const std::string jpsi_cos_mu_plus_jpsi_pt_15to20_name = "jpsi_cos_mu_plus_jpsi_pt_15to20";
      jpsi_cos_mu_plus_jpsi_pt_15to20_ = tdir.make<TH1D>(jpsi_cos_mu_plus_jpsi_pt_15to20_name.c_str(), jpsi_cos_mu_plus_jpsi_pt_15to20_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_jpsi_pt_15to20_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_15to20 cos(theta)");
      jpsi_cos_mu_plus_jpsi_pt_15to20_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambdaNeg1
      const std::string jpsi_cos_mu_plus_lambdaNeg1_name = "jpsi_cos_mu_plus_lambdaNeg1";
      jpsi_cos_mu_plus_lambdaNeg1_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambdaNeg1_name.c_str(), jpsi_cos_mu_plus_lambdaNeg1_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambdaNeg1_->GetXaxis()->SetTitle("jpsi mu_plus cos(theta)");
      jpsi_cos_mu_plus_lambdaNeg1_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5
      const std::string jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5_name = "jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5";
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5_name.c_str(), jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_8to8p5 cos(theta)");
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9
      const std::string jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9_name = "jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9";
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9_name.c_str(), jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_8p5to9 cos(theta)");
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10
      const std::string jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10_name = "jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10";
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10_name.c_str(), jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_9to10 cos(theta)");
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15
      const std::string jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15_name = "jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15";
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15_name.c_str(), jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_10to15 cos(theta)");
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20
      const std::string jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20_name = "jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20";
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20_name.c_str(), jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_15to20 cos(theta)");
      jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambda1
      const std::string jpsi_cos_mu_plus_lambda1_name = "jpsi_cos_mu_plus_lambda1";
      jpsi_cos_mu_plus_lambda1_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambda1_name.c_str(), jpsi_cos_mu_plus_lambda1_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambda1_->GetXaxis()->SetTitle("jpsi mu_plus cos(theta)");
      jpsi_cos_mu_plus_lambda1_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5
      const std::string jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5_name = "jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5";
      jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5_name.c_str(), jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_8to8p5 cos(theta)");
      jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9
      const std::string jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9_name = "jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9";
      jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9_name.c_str(), jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_8p5to9 cos(theta)");
      jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10
      const std::string jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10_name = "jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10";
      jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10_name.c_str(), jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_9to10 cos(theta)");
      jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15
      const std::string jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15_name = "jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15";
      jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15_name.c_str(), jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_10to15 cos(theta)");
      jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20
      const std::string jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20_name = "jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20";
      jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20_ = tdir.make<TH1D>(jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20_name.c_str(), jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20_->GetXaxis()->SetTitle("jpsi mu_plus_jpsi_pt_15to20 cos(theta)");
      jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20_->GetYaxis()->SetTitle("Counts / 0.01 ");

      // jpsi_cos_mu_minus
      const std::string jpsi_cos_mu_minus_name = "jpsi_cos_mu_minus";
      jpsi_cos_mu_minus_ = tdir.make<TH1D>(jpsi_cos_mu_minus_name.c_str(), jpsi_cos_mu_minus_name.c_str(), 200, -1., 1.);
      jpsi_cos_mu_minus_->GetXaxis()->SetTitle("jpsi mu_minus cos(theta)");
      jpsi_cos_mu_minus_->GetYaxis()->SetTitle("Counts / 0.01 ");

      //TODO variable bin size to use this histogram as an efficiency map
      //double bins[] = {8,8.5,9,9.5,10,12,15,100};
      double bins[] = {8.5,10,14,18,30,100};
      //int bin_number = 7;
      int bin_number = 5;
      // jpsi_pt_vs_rap
      const std::string jpsi_pt_vs_rap_name = "jpsi_pt_vs_rap";
      jpsi_pt_vs_rap_ = tdir.make<TH2D>(jpsi_pt_vs_rap_name.c_str(), jpsi_pt_vs_rap_name.c_str(), 21, -2.1, 2.1, bin_number, bins);
      jpsi_pt_vs_rap_->GetXaxis()->SetTitle("jpsi Rapidity");
      jpsi_pt_vs_rap_->GetYaxis()->SetTitle("jpsi Pt");

      const std::string jpsi_pt_vs_rap_polarization_longname = "jpsi_pt_vs_rap_polarization_long";
      jpsi_pt_vs_rap_polarization_long = tdir.make<TH2D>(jpsi_pt_vs_rap_polarization_longname.c_str(), jpsi_pt_vs_rap_polarization_longname.c_str(), 21, -2.1, 2.1, bin_number, bins);
      jpsi_pt_vs_rap_polarization_long->GetXaxis()->SetTitle("jpsi Rapidity");
      jpsi_pt_vs_rap_polarization_long->GetYaxis()->SetTitle("jpsi Pt");

      const std::string jpsi_pt_vs_rap_polarization_TPlusZeroname = "jpsi_pt_vs_rap_polarization_TPlusZero";
      jpsi_pt_vs_rap_polarization_TPlusZero = tdir.make<TH2D>(jpsi_pt_vs_rap_polarization_TPlusZeroname.c_str(), jpsi_pt_vs_rap_polarization_TPlusZeroname.c_str(), 21, -2.1, 2.1, bin_number, bins);
      jpsi_pt_vs_rap_polarization_TPlusZero->GetXaxis()->SetTitle("jpsi Rapidity");
      jpsi_pt_vs_rap_polarization_TPlusZero->GetYaxis()->SetTitle("jpsi Pt");

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

      // jpsi_zmumupt_difference_
      const std::string jpsi_zmumupt_difference__name = "zmumupt - jpsipT";
      jpsi_zmumupt_difference_ = tdir.make<TH1D>(jpsi_zmumupt_difference__name.c_str(), jpsi_zmumupt_difference__name.c_str(), 200, -100., 100.);
      jpsi_zmumupt_difference_->GetXaxis()->SetTitle("pT [GeV]");
      jpsi_zmumupt_difference_->GetYaxis()->SetTitle("Counts / 1 GeV ");

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
      jpsi_tau_xy_very_fine_ptUnder10_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_ptUnder10_name.c_str(), jpsi_tau_xy_very_fine_ptUnder10_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_ptUnder10_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_ptUnder10_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_pt10to15
      const std::string jpsi_tau_xy_very_fine_pt10to15_name = "jpsi_tau_xy_very_fine_pt_10_to_15";
      jpsi_tau_xy_very_fine_pt10to15_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_pt10to15_name.c_str(), jpsi_tau_xy_very_fine_pt10to15_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_pt10to15_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_pt10to15_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_pt15to20
      const std::string jpsi_tau_xy_very_fine_pt15to20_name = "jpsi_tau_xy_very_fine_pt_15_to_20";
      jpsi_tau_xy_very_fine_pt15to20_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_pt15to20_name.c_str(), jpsi_tau_xy_very_fine_pt15to20_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_pt15to20_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_pt15to20_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_pt20to25
      const std::string jpsi_tau_xy_very_fine_pt20to25_name = "jpsi_tau_xy_very_fine_pt_20_to_25";
      jpsi_tau_xy_very_fine_pt20to25_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_pt20to25_name.c_str(), jpsi_tau_xy_very_fine_pt20to25_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_pt20to25_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_pt20to25_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_pt25to30
      const std::string jpsi_tau_xy_very_fine_pt25to30_name = "jpsi_tau_xy_very_fine_pt_25_to_30";
      jpsi_tau_xy_very_fine_pt25to30_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_pt25to30_name.c_str(), jpsi_tau_xy_very_fine_pt25to30_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_pt25to30_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_pt25to30_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_ptAbove30
      const std::string jpsi_tau_xy_very_fine_ptAbove30_name = "jpsi_tau_xy_very_fine_ptAbove30";
      jpsi_tau_xy_very_fine_ptAbove30_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_ptAbove30_name.c_str(), jpsi_tau_xy_very_fine_ptAbove30_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_ptAbove30_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_ptAbove30_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap0_0to0_3
      const std::string jpsi_tau_xy_very_fine_rap0_0to0_3_name = "jpsi_tau_xy_very_fine_rap_0.0_to_0.3";
      jpsi_tau_xy_very_fine_rap0_0to0_3_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_0to0_3_name.c_str(), jpsi_tau_xy_very_fine_rap0_0to0_3_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap0_0to0_3_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap0_0to0_3_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap0_3to0_6
      const std::string jpsi_tau_xy_very_fine_rap0_3to0_6_name = "jpsi_tau_xy_very_fine_rap_0.3_to_0.6";
      jpsi_tau_xy_very_fine_rap0_3to0_6_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_3to0_6_name.c_str(), jpsi_tau_xy_very_fine_rap0_3to0_6_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap0_3to0_6_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap0_3to0_6_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap0_6to0_9
      const std::string jpsi_tau_xy_very_fine_rap0_6to0_9_name = "jpsi_tau_xy_very_fine_rap_0.6_to_0.9";
      jpsi_tau_xy_very_fine_rap0_6to0_9_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_6to0_9_name.c_str(), jpsi_tau_xy_very_fine_rap0_6to0_9_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap0_6to0_9_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap0_6to0_9_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap0_9to1_2
      const std::string jpsi_tau_xy_very_fine_rap0_9to1_2_name = "jpsi_tau_xy_very_fine_rap_0.9_to_1.2";
      jpsi_tau_xy_very_fine_rap0_9to1_2_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_9to1_2_name.c_str(), jpsi_tau_xy_very_fine_rap0_9to1_2_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap0_9to1_2_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap0_9to1_2_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap1_2to1_5
      const std::string jpsi_tau_xy_very_fine_rap1_2to1_5_name = "jpsi_tau_xy_very_fine_rap_1.2_to_1.5";
      jpsi_tau_xy_very_fine_rap1_2to1_5_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap1_2to1_5_name.c_str(), jpsi_tau_xy_very_fine_rap1_2to1_5_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap1_2to1_5_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap1_2to1_5_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap1_5to1_8
      const std::string jpsi_tau_xy_very_fine_rap1_5to1_8_name = "jpsi_tau_xy_very_fine_rap_1.5_to_1.8";
      jpsi_tau_xy_very_fine_rap1_5to1_8_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap1_5to1_8_name.c_str(), jpsi_tau_xy_very_fine_rap1_5to1_8_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap1_5to1_8_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap1_5to1_8_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap1_8to2_1
      const std::string jpsi_tau_xy_very_fine_rap1_8to2_1_name = "jpsi_tau_xy_very_fine_rap_1.8_to_2.1";
      jpsi_tau_xy_very_fine_rap1_8to2_1_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap1_8to2_1_name.c_str(), jpsi_tau_xy_very_fine_rap1_8to2_1_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap1_8to2_1_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap1_8to2_1_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap2_1to2_4
      const std::string jpsi_tau_xy_very_fine_rap2_1to2_4_name = "jpsi_tau_xy_very_fine_rap_2.1_to_2.4";
      jpsi_tau_xy_very_fine_rap2_1to2_4_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap2_1to2_4_name.c_str(), jpsi_tau_xy_very_fine_rap2_1to2_4_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap2_1to2_4_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap2_1to2_4_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap0_0to0_9_pt10to15
      const std::string jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_name = "jpsi_tau_xy_very_fine_rap_0.0_to_0.9_pt10to15";
      jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_name.c_str(), jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap0_9to1_2_pt10to15
      const std::string jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_name = "jpsi_tau_xy_very_fine_rap_0.9_to_1.2_pt10to15";
      jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_name.c_str(), jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_rap1_2to2_1_pt10to15
      const std::string jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_name = "jpsi_tau_xy_very_fine_rap_1.2_to_2.1_pt10to15";
      jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_name.c_str(), jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_above_12_tracker_layers
      const std::string jpsi_tau_xy_very_fine_above_12_tracker_layers_name = "jpsi_tau_xy_very_fine_above_12_tracker_layers_";
      jpsi_tau_xy_very_fine_above_12_tracker_layers_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_above_12_tracker_layers_name.c_str(), jpsi_tau_xy_very_fine_above_12_tracker_layers_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_above_12_tracker_layers_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_above_12_tracker_layers_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_very_fine_similar_pt_muons_
      const std::string jpsi_tau_xy_very_fine_similar_pt_muons_name = "jpsi_tau_xy_very_fine_similar_pt_muons_";
      jpsi_tau_xy_very_fine_similar_pt_muons_ = tdir.make<TH1D>(jpsi_tau_xy_very_fine_similar_pt_muons_name.c_str(), jpsi_tau_xy_very_fine_similar_pt_muons_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_very_fine_similar_pt_muons_->GetXaxis()->SetTitle("tau_xy_very_fine [ps]");
      jpsi_tau_xy_very_fine_similar_pt_muons_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

      // jpsi_tau_xy_dimuon_continuum_bg_
      const std::string jpsi_tau_xy_dimuon_continuum_bg_name = "jpsi_tau_xy_dimuon_continuum_bg_";
      jpsi_tau_xy_dimuon_continuum_bg_ = tdir.make<TH1D>(jpsi_tau_xy_dimuon_continuum_bg_name.c_str(), jpsi_tau_xy_dimuon_continuum_bg_name.c_str(), 200, -10., 10.);
      jpsi_tau_xy_dimuon_continuum_bg_->GetXaxis()->SetTitle("tau_xy [ps]");
      jpsi_tau_xy_dimuon_continuum_bg_->GetYaxis()->SetTitle("Counts / 0.1 ps ");

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
      jpsi_mass_vs_chi2_ = tdir.make<TH2D>(jpsi_mass_vs_chi2_name.c_str(), jpsi_mass_vs_chi2_name.c_str(), 100, 0., 10., 50, 0., 10.);
      jpsi_mass_vs_chi2_->GetXaxis()->SetTitle("chi2");
      jpsi_mass_vs_chi2_->GetYaxis()->SetTitle("mass");

      // jpsi_tau_xy_vs_tau_z
      const std::string jpsi_tau_xy_vs_tau_z_name = "jpsi tau_xy vs tau_z";
      jpsi_tau_xy_vs_tau_z_ = tdir.make<TH2D>(jpsi_tau_xy_vs_tau_z_name.c_str(), jpsi_tau_xy_vs_tau_z_name.c_str(), 100, -10., 10., 100, -10., 10.);
      jpsi_tau_xy_vs_tau_z_->GetXaxis()->SetTitle("tau_xy [ps]");
      jpsi_tau_xy_vs_tau_z_->GetYaxis()->SetTitle("tau_z [ps]");

      // dimuon_pt_vs_zee_pt
      const std::string dimuon_pt_vs_zee_pt_name = "dimuon_pt_vs_zee_pt";
      dimuon_pt_vs_zee_pt_ = tdir.make<TH2D>(dimuon_pt_vs_zee_pt_name.c_str(), dimuon_pt_vs_zee_pt_name.c_str(), 200, 0.0, 200.0, 200, 0.0, 200.);
      dimuon_pt_vs_zee_pt_->GetXaxis()->SetTitle("J/#psi p_{T} [GeV]");
      dimuon_pt_vs_zee_pt_->GetYaxis()->SetTitle("Z->ee p_{T} [GeV]");

      // dimuon_pt_vs_zmumu_pt
      const std::string dimuon_pt_vs_zmumu_pt_name = "dimuon_pt_vs_zmumu_pt";
      dimuon_pt_vs_zmumu_pt_ = tdir.make<TH2D>(dimuon_pt_vs_zmumu_pt_name.c_str(), dimuon_pt_vs_zmumu_pt_name.c_str(), 200, 0.0, 200.0, 200, 0.0, 200.);
      dimuon_pt_vs_zmumu_pt_->GetXaxis()->SetTitle("J/#psi p_{T} [GeV]");
      dimuon_pt_vs_zmumu_pt_->GetYaxis()->SetTitle("Z->#mu#mu p_{T} [GeV]");

      // dimuon_mass_vs_dimuon_tau_xy
      const std::string dimuon_mass_vs_dimuon_tau_xy_name = "dimuon_mass_vs_dimuon_tau_xy";
      dimuon_mass_vs_dimuon_tau_xy_ = tdir.make<TH2D>(dimuon_mass_vs_dimuon_tau_xy_name.c_str(), dimuon_mass_vs_dimuon_tau_xy_name.c_str(), 200, 2.0, 12.0, 200, -10., 10.);
      dimuon_mass_vs_dimuon_tau_xy_->GetXaxis()->SetTitle("dimuon mass [GeV]");
      dimuon_mass_vs_dimuon_tau_xy_->GetYaxis()->SetTitle("tau_xy [ps]");

      // dimuon_mass_vs_dimuon_tau_xy_fine
      const std::string dimuon_mass_vs_dimuon_tau_xy_fine_name = "dimuon_mass_vs_dimuon_tau_xy_fine";
      dimuon_mass_vs_dimuon_tau_xy_fine_ = tdir.make<TH2D>(dimuon_mass_vs_dimuon_tau_xy_fine_name.c_str(), dimuon_mass_vs_dimuon_tau_xy_fine_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      dimuon_mass_vs_dimuon_tau_xy_fine_->GetXaxis()->SetTitle("dimuon mass [GeV]");
      dimuon_mass_vs_dimuon_tau_xy_fine_->GetYaxis()->SetTitle("tau_xy [ps]");

      // dimuon_mass_vs_dimuon_tau_xy_8p5to10p0
      const std::string dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name = "dimuon_mass_vs_dimuon_tau_xy_8p5to10p0";
      dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_ = tdir.make<TH2D>(dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name.c_str(), dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->GetXaxis()->SetTitle("dimuon mass (pt 8.5 to 10.0) [GeV]");
      dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->GetYaxis()->SetTitle("tau_xy [ps]");

      // dimuon_mass_vs_dimuon_tau_xy_10to14
      const std::string dimuon_mass_vs_dimuon_tau_xy_10to14_name = "dimuon_mass_vs_dimuon_tau_xy_10to14";
      dimuon_mass_vs_dimuon_tau_xy_10to14_ = tdir.make<TH2D>(dimuon_mass_vs_dimuon_tau_xy_10to14_name.c_str(), dimuon_mass_vs_dimuon_tau_xy_10to14_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      dimuon_mass_vs_dimuon_tau_xy_10to14_->GetXaxis()->SetTitle("dimuon mass (pt 10 to 14)[GeV]");
      dimuon_mass_vs_dimuon_tau_xy_10to14_->GetYaxis()->SetTitle("tau_xy [ps]");

      // dimuon_mass_vs_dimuon_tau_xy_14to18
      const std::string dimuon_mass_vs_dimuon_tau_xy_14to18_name = "dimuon_mass_vs_dimuon_tau_xy_14to18";
      dimuon_mass_vs_dimuon_tau_xy_14to18_ = tdir.make<TH2D>(dimuon_mass_vs_dimuon_tau_xy_14to18_name.c_str(), dimuon_mass_vs_dimuon_tau_xy_14to18_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      dimuon_mass_vs_dimuon_tau_xy_14to18_->GetXaxis()->SetTitle("dimuon mass (pt 14 to 18)[GeV]");
      dimuon_mass_vs_dimuon_tau_xy_14to18_->GetYaxis()->SetTitle("tau_xy [ps]");

      // dimuon_mass_vs_dimuon_tau_xy_18to30
      const std::string dimuon_mass_vs_dimuon_tau_xy_18to30_name = "dimuon_mass_vs_dimuon_tau_xy_18to30";
      dimuon_mass_vs_dimuon_tau_xy_18to30_ = tdir.make<TH2D>(dimuon_mass_vs_dimuon_tau_xy_18to30_name.c_str(), dimuon_mass_vs_dimuon_tau_xy_18to30_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      dimuon_mass_vs_dimuon_tau_xy_18to30_->GetXaxis()->SetTitle("dimuon mass (pt 18 to 30) [GeV]");
      dimuon_mass_vs_dimuon_tau_xy_18to30_->GetYaxis()->SetTitle("tau_xy [ps]");

      // dimuon_mass_vs_dimuon_tau_xy_30to100
      const std::string dimuon_mass_vs_dimuon_tau_xy_30to100_name = "dimuon_mass_vs_dimuon_tau_xy_30to100";
      dimuon_mass_vs_dimuon_tau_xy_30to100_ = tdir.make<TH2D>(dimuon_mass_vs_dimuon_tau_xy_30to100_name.c_str(), dimuon_mass_vs_dimuon_tau_xy_30to100_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      dimuon_mass_vs_dimuon_tau_xy_30to100_->GetXaxis()->SetTitle("dimuon mass (pt 30 to 100) [GeV]");
      dimuon_mass_vs_dimuon_tau_xy_30to100_->GetYaxis()->SetTitle("tau_xy [ps]");

      // low_rap_dimuon_mass_vs_dimuon_tau_xy_fine
      const std::string low_rap_dimuon_mass_vs_dimuon_tau_xy_fine_name = "low_rap_dimuon_mass_vs_dimuon_tau_xy_fine";
      low_rap_dimuon_mass_vs_dimuon_tau_xy_fine_ = tdir.make<TH2D>(low_rap_dimuon_mass_vs_dimuon_tau_xy_fine_name.c_str(), low_rap_dimuon_mass_vs_dimuon_tau_xy_fine_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      low_rap_dimuon_mass_vs_dimuon_tau_xy_fine_->GetXaxis()->SetTitle("dimuon mass [GeV]");
      low_rap_dimuon_mass_vs_dimuon_tau_xy_fine_->GetYaxis()->SetTitle("tau_xy [ps]");

      // low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0
      const std::string low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name = "low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0";
      low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_ = tdir.make<TH2D>(low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name.c_str(), low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->GetXaxis()->SetTitle("dimuon mass (pt 8.5 to 10.0) [GeV]");
      low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->GetYaxis()->SetTitle("tau_xy [ps]");

      // low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14
      const std::string low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_name = "low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14";
      low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_ = tdir.make<TH2D>(low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_name.c_str(), low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_->GetXaxis()->SetTitle("dimuon mass (pt 10 to 14)[GeV]");
      low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_->GetYaxis()->SetTitle("tau_xy [ps]");

      // low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18
      const std::string low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_name = "low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18";
      low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_ = tdir.make<TH2D>(low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_name.c_str(), low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_->GetXaxis()->SetTitle("dimuon mass (pt 14 to 18)[GeV]");
      low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_->GetYaxis()->SetTitle("tau_xy [ps]");

      // low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30
      const std::string low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_name = "low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30";
      low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_ = tdir.make<TH2D>(low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_name.c_str(), low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_->GetXaxis()->SetTitle("dimuon mass (pt 18 to 30) [GeV]");
      low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_->GetYaxis()->SetTitle("tau_xy [ps]");

      // low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100
      const std::string low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_name = "low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100";
      low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_ = tdir.make<TH2D>(low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_name.c_str(), low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_->GetXaxis()->SetTitle("dimuon mass (pt 30 to 100) [GeV]");
      low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_->GetYaxis()->SetTitle("tau_xy [ps]");

      // high_rap_dimuon_mass_vs_dimuon_tau_xy_fine
      const std::string high_rap_dimuon_mass_vs_dimuon_tau_xy_fine_name = "high_rap_dimuon_mass_vs_dimuon_tau_xy_fine";
      high_rap_dimuon_mass_vs_dimuon_tau_xy_fine_ = tdir.make<TH2D>(high_rap_dimuon_mass_vs_dimuon_tau_xy_fine_name.c_str(), high_rap_dimuon_mass_vs_dimuon_tau_xy_fine_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      high_rap_dimuon_mass_vs_dimuon_tau_xy_fine_->GetXaxis()->SetTitle("dimuon mass [GeV]");
      high_rap_dimuon_mass_vs_dimuon_tau_xy_fine_->GetYaxis()->SetTitle("tau_xy [ps]");

      // high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0
      const std::string high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name = "high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0";
      high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_ = tdir.make<TH2D>(high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name.c_str(), high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->GetXaxis()->SetTitle("dimuon mass (pt 8.5 to 10.0) [GeV]");
      high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->GetYaxis()->SetTitle("tau_xy [ps]");

      // high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14
      const std::string high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_name = "high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14";
      high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_ = tdir.make<TH2D>(high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_name.c_str(), high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_->GetXaxis()->SetTitle("dimuon mass (pt 10 to 14)[GeV]");
      high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_->GetYaxis()->SetTitle("tau_xy [ps]");

      // high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18
      const std::string high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_name = "high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18";
      high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_ = tdir.make<TH2D>(high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_name.c_str(), high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_->GetXaxis()->SetTitle("dimuon mass (pt 14 to 18)[GeV]");
      high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_->GetYaxis()->SetTitle("tau_xy [ps]");

      // high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30
      const std::string high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_name = "high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30";
      high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_ = tdir.make<TH2D>(high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_name.c_str(), high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_->GetXaxis()->SetTitle("dimuon mass (pt 18 to 30) [GeV]");
      high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_->GetYaxis()->SetTitle("tau_xy [ps]");

      // high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100
      const std::string high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_name = "high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100";
      high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_ = tdir.make<TH2D>(high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_name.c_str(), high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_->GetXaxis()->SetTitle("dimuon mass (pt 30 to 100) [GeV]");
      high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_->GetYaxis()->SetTitle("tau_xy [ps]");

      // dimuon_mass_vs_dimuon_tau_xy_fine_weighted
      const std::string dimuon_mass_vs_dimuon_tau_xy_fine_weighted_name = "dimuon_mass_vs_dimuon_tau_xy_fine_weighted";
      dimuon_mass_vs_dimuon_tau_xy_fine_weighted_ = tdir.make<TH2D>(dimuon_mass_vs_dimuon_tau_xy_fine_weighted_name.c_str(), dimuon_mass_vs_dimuon_tau_xy_fine_weighted_name.c_str(), 50, 2.85, 3.35, 750, -5., 10.);
      dimuon_mass_vs_dimuon_tau_xy_fine_weighted_->GetXaxis()->SetTitle("dimuon mass [GeV]");
      dimuon_mass_vs_dimuon_tau_xy_fine_weighted_->GetYaxis()->SetTitle("tau_xy [ps]");

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
      mu1_pt_->GetXaxis()->SetTitle("p_{T,mu_{1}}");
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

      // mu0_efficiency_
      const std::string mu0_efficiency_name = "Efficiency_{mu_{0}}";
      mu0_efficiency_ = tdir.make<TH1D>(mu0_efficiency_name.c_str(), mu0_efficiency_name.c_str(), 100, 0.0, 1.0);
      mu0_efficiency_->GetXaxis()->SetTitle("Efficiency_{mu_{0}}");
      mu0_efficiency_->GetYaxis()->SetTitle("Counts / 0.2");

      // mu1_efficiency_
      const std::string mu1_efficiency_name = "Efficiency_{mu_{1}}";
      mu1_efficiency_ = tdir.make<TH1D>(mu1_efficiency_name.c_str(), mu1_efficiency_name.c_str(), 100, 0.0, 1.0);
      mu1_efficiency_->GetXaxis()->SetTitle("Efficiency_{mu_{1}}");
      mu1_efficiency_->GetYaxis()->SetTitle("Counts / 0.01");

      // mu0_scale_factor_
      const std::string mu0_scale_factor_name = "scale_factor_{mu_{0}}";
      mu0_scale_factor_ = tdir.make<TH1D>(mu0_scale_factor_name.c_str(), mu0_scale_factor_name.c_str(), 200, 0.0, 2.0);
      mu0_scale_factor_->GetXaxis()->SetTitle("scale_factor_{mu_{0}}");
      mu0_scale_factor_->GetYaxis()->SetTitle("Counts / 0.2");

      // mu1_scale_factor_
      const std::string mu1_scale_factor_name = "scale_factor_{mu_{1}}";
      mu1_scale_factor_ = tdir.make<TH1D>(mu1_scale_factor_name.c_str(), mu1_scale_factor_name.c_str(), 200, 0.0, 2.0);
      mu1_scale_factor_->GetXaxis()->SetTitle("scale_factor_{mu_{1}}");
      mu1_scale_factor_->GetYaxis()->SetTitle("Counts / 0.01");

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

      // mu0_deltaR_z_muon_
      const std::string mu0_deltaR_z_muon_name = "deltaR_z_mu0";
      mu0_deltaR_z_muon_ = tdir.make<TH1D>(mu0_deltaR_z_muon_name.c_str(), mu0_deltaR_z_muon_name.c_str(), 500, 0, 5);
      mu0_deltaR_z_muon_->GetXaxis()->SetTitle("deltaR_z_mu0");
      mu0_deltaR_z_muon_->GetYaxis()->SetTitle("Counts");

      // mu1_deltaR_z_muon_
      const std::string mu1_deltaR_z_muon_name = "deltaR_z_mu1";
      mu1_deltaR_z_muon_ = tdir.make<TH1D>(mu1_deltaR_z_muon_name.c_str(), mu1_deltaR_z_muon_name.c_str(), 500, 0, 5);
      mu1_deltaR_z_muon_->GetXaxis()->SetTitle("deltaR_z_mu1");
      mu1_deltaR_z_muon_->GetYaxis()->SetTitle("Counts");

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
      jet_btag_discriminator_ = tdir.make<TH1D>(jet_btag_discriminator_name.c_str(), jet_btag_discriminator_name.c_str(), 140, -120., 20.);
      jet_btag_discriminator_->GetXaxis()->SetTitle("discriminator");
      jet_btag_discriminator_->GetYaxis()->SetTitle("Counts / 1.0");

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
      muon_jet_btag_discriminator_ = tdir.make<TH1D>(muon_jet_btag_discriminator_name.c_str(), muon_jet_btag_discriminator_name.c_str(), 140, -120., 20.);
      muon_jet_btag_discriminator_->GetXaxis()->SetTitle("discriminator");
      muon_jet_btag_discriminator_->GetYaxis()->SetTitle("Counts / 1.0");

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

      //TODO think about where this should be ordered - organize by folders?
      // z_from_muons_vtx_x_
      const std::string z_from_muons_vtx_x_name = "z_from_muons vertex position x";
      z_from_muons_vtx_x_ = tdir.make<TH1D>(z_from_muons_vtx_x_name.c_str(), z_from_muons_vtx_x_name.c_str(), 2000, -10., 10.);
      z_from_muons_vtx_x_->GetXaxis()->SetTitle("Position [cm]");
      z_from_muons_vtx_x_->GetYaxis()->SetTitle("z_from_muons Vertex / 0.001 cm");

      // z_from_muons_vtx_y_
      const std::string z_from_muons_vtx_y_name = "z_from_muons vertex position y";
      z_from_muons_vtx_y_ = tdir.make<TH1D>(z_from_muons_vtx_y_name.c_str(), z_from_muons_vtx_y_name.c_str(), 2000, -10., 10.);
      z_from_muons_vtx_y_->GetXaxis()->SetTitle("Position [cm]");
      z_from_muons_vtx_y_->GetYaxis()->SetTitle("z_from_muons Vertex / 0.001 cm");

      // z_from_muons_vtx_z_
      const std::string z_from_muons_vtx_z_name = "z_from_muons vertex position z";
      z_from_muons_vtx_z_ = tdir.make<TH1D>(z_from_muons_vtx_z_name.c_str(), z_from_muons_vtx_z_name.c_str(), 200, -10., 10.);
      z_from_muons_vtx_z_->GetXaxis()->SetTitle("Position [cm]");
      z_from_muons_vtx_z_->GetYaxis()->SetTitle("z Vertex / 0.01 cm");

      //primary_vtx_x_minus_z_vtx_x_  
      const std::string primary_vtx_x_minus_z_vtx_x_name = "primary_vtx_x_minus_z_vtx_position_x";
      primary_vtx_x_minus_z_vtx_x_ = tdir.make<TH1D>(primary_vtx_x_minus_z_vtx_x_name.c_str(), primary_vtx_x_minus_z_vtx_x_name.c_str(), 2000, -10., 10.);
      primary_vtx_x_minus_z_vtx_x_->GetXaxis()->SetTitle("#Delta x Position [cm]");
      primary_vtx_x_minus_z_vtx_x_->GetYaxis()->SetTitle("z_from_muons Vertex / 0.001 cm");

      //primary_vtx_y_minus_z_vtx_y_  
      const std::string primary_vtx_y_minus_z_vtx_y_name = "primary_vtx_y_minus_z_vtx_position_y";
      primary_vtx_y_minus_z_vtx_y_ = tdir.make<TH1D>(primary_vtx_y_minus_z_vtx_y_name.c_str(), primary_vtx_y_minus_z_vtx_y_name.c_str(), 2000, -10., 10.);
      primary_vtx_y_minus_z_vtx_y_->GetXaxis()->SetTitle("#Delta y Position [cm]");
      primary_vtx_y_minus_z_vtx_y_->GetYaxis()->SetTitle("z_from_muons Vertex / 0.001 cm");

      //primary_vtx_z_minus_z_vtx_z_ 
      const std::string primary_vtx_z_minus_z_vtx_z_name = "primary_vtx_z_minus_z_vtx_position_z";
      primary_vtx_z_minus_z_vtx_z_ = tdir.make<TH1D>(primary_vtx_z_minus_z_vtx_z_name.c_str(), primary_vtx_z_minus_z_vtx_z_name.c_str(), 2000, -10., 10.);
      primary_vtx_z_minus_z_vtx_z_->GetXaxis()->SetTitle("#Delta z Position [cm]");
      primary_vtx_z_minus_z_vtx_z_->GetYaxis()->SetTitle("z Vertex / 0.01 cm");

      //primary_vtx_x_minus_zmumu_vtx_x_  
      const std::string primary_vtx_x_minus_zmumu_vtx_x_name = "primary_vtx_x_minus_zmumu_vtx_position_x";
      primary_vtx_x_minus_zmumu_vtx_x_ = tdir.make<TH1D>(primary_vtx_x_minus_zmumu_vtx_x_name.c_str(), primary_vtx_x_minus_zmumu_vtx_x_name.c_str(), 2000, -10., 10.);
      primary_vtx_x_minus_zmumu_vtx_x_->GetXaxis()->SetTitle("#Delta x Position [cm]");
      primary_vtx_x_minus_zmumu_vtx_x_->GetYaxis()->SetTitle("z_from_muons Vertex / 0.001 cm");

      //primary_vtx_y_minus_zmumu_vtx_y_  
      const std::string primary_vtx_y_minus_zmumu_vtx_y_name = "primary_vtx_y_minus_zmumu_vtx_position_y";
      primary_vtx_y_minus_zmumu_vtx_y_ = tdir.make<TH1D>(primary_vtx_y_minus_zmumu_vtx_y_name.c_str(), primary_vtx_y_minus_zmumu_vtx_y_name.c_str(), 2000, -10., 10.);
      primary_vtx_y_minus_zmumu_vtx_y_->GetXaxis()->SetTitle("#Delta y Position [cm]");
      primary_vtx_y_minus_zmumu_vtx_y_->GetYaxis()->SetTitle("z_from_muons Vertex / 0.001 cm");

      //primary_vtx_z_minus_zmumu_vtx_z_ 
      const std::string primary_vtx_z_minus_zmumu_vtx_z_name = "primary_vtx_z_minus_zmumu_vtx_position_z";
      primary_vtx_z_minus_zmumu_vtx_z_ = tdir.make<TH1D>(primary_vtx_z_minus_zmumu_vtx_z_name.c_str(), primary_vtx_z_minus_zmumu_vtx_z_name.c_str(), 2000, -10., 10.);
      primary_vtx_z_minus_zmumu_vtx_z_->GetXaxis()->SetTitle("#Delta z Position [cm]");
      primary_vtx_z_minus_zmumu_vtx_z_->GetYaxis()->SetTitle("Events / 0.001 cm");

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

      // primary_vtx_x_vs_zmuons_vtx_x_
      const std::string primary_vtx_x_vs_zmuons_vtx_x_name = "primary vertex x pos vs Zmumu vertex x pos";
      primary_vtx_x_vs_zmuons_vtx_x_ = tdir.make<TH2D>(primary_vtx_x_vs_zmuons_vtx_x_name.c_str(), primary_vtx_x_vs_zmuons_vtx_x_name.c_str(), 200, -1., 1., 200, -1., 1.);
      primary_vtx_x_vs_zmuons_vtx_x_->GetXaxis()->SetTitle("Z Muons Vertex x Pos [cm]");
      primary_vtx_x_vs_zmuons_vtx_x_->GetYaxis()->SetTitle("Primary Vertex x Pos [cm]");

      // primary_vtx_y_vs_zmuons_vtx_y_
      const std::string primary_vtx_y_vs_zmuons_vtx_y_name = "primary vertex y pos vs Zmumu vertex y pos";
      primary_vtx_y_vs_zmuons_vtx_y_ = tdir.make<TH2D>(primary_vtx_y_vs_zmuons_vtx_y_name.c_str(), primary_vtx_y_vs_zmuons_vtx_y_name.c_str(), 200, -1., 1., 200, -1., 1.);
      primary_vtx_y_vs_zmuons_vtx_y_->GetXaxis()->SetTitle("Z Muons Vertex y Pos [cm]");
      primary_vtx_y_vs_zmuons_vtx_y_->GetYaxis()->SetTitle("Primary Vertex y Pos [cm]");

      // primary_vtx_z_vs_zmuons_vtx_z_
      const std::string primary_vtx_z_vs_zmuons_vtx_z_name = "primary vertex z pos vs Zmumu vertex z pos";
      primary_vtx_z_vs_zmuons_vtx_z_ = tdir.make<TH2D>(primary_vtx_z_vs_zmuons_vtx_z_name.c_str(), primary_vtx_z_vs_zmuons_vtx_z_name.c_str(), 200, -10., 10., 200, -10., 10.);
      primary_vtx_z_vs_zmuons_vtx_z_->GetXaxis()->SetTitle("Z Muons Vertex z Pos [cm]");
      primary_vtx_z_vs_zmuons_vtx_z_->GetYaxis()->SetTitle("Primary Vertex z Pos [cm]");

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

      // baseweights
      const std::string baseweights_name = "Base Weight";
      baseweights_ = tdir.make<TH1D>(baseweights_name.c_str(), baseweights_name.c_str(), 500, 0., 5.);
      baseweights_->GetXaxis()->SetTitle("Weight");
      baseweights_->GetYaxis()->SetTitle("Events");

      // pileup
      const std::string pileup_name = "N_{Vertices}";
      pileup_ = tdir.make<TH1D>(pileup_name.c_str(), pileup_name.c_str(), 100, 0., 100.);
      pileup_->GetXaxis()->SetTitle("Number of Vertices");
      pileup_->GetYaxis()->SetTitle("Counts");

    }

  void ZFinderPlotter::Fill(
      const ZFinderEvent& zfe,
      const double EVENT_WEIGHT
      ) {
    /*
     * Given a zfe, fills all the histograms.
     */
    // Z Info
    double event_weight = zfe.event_weight;
    if (!USE_MC_) {
      for (unsigned int i = 0; i < zfe.reco_jpsi.m.size() ; ++i ) {
        if (!zfe.is_real_data && APPLY_SOFT_MUONS_ && zfe.reco_jpsi.has_soft_id_muons.at(i) &&
            zfe.reco_jpsi.has_high_pt_muons.at(i) && zfe.reco_jpsi.has_muons_in_eta_window.at(i) ) {
          event_weight = event_weight * zfe.reco_jpsi.jpsi_scale_factor.at(i);
          //TODO decide how to handle multiple j/psi case!!
          break;
        }
      }

      z_mass_all_->Fill(zfe.reco_z.m, event_weight);
      z_mass_coarse_->Fill(zfe.reco_z.m, event_weight);
      z_mass_fine_->Fill(zfe.reco_z.m, event_weight);
      z_rapidity_->Fill(zfe.reco_z.y, event_weight);
      z_pt_->Fill(zfe.reco_z.pt, event_weight);
      z_vtx_x_->Fill(zfe.reco_z.vtx_x, event_weight);
      z_vtx_y_->Fill(zfe.reco_z.vtx_y, event_weight);
      z_vtx_z_->Fill(zfe.reco_z.vtx_z, event_weight);
      z_vtx_prob_->Fill(zfe.reco_z.vtx_prob, event_weight);
      phistar_->Fill(zfe.reco_z.phistar, event_weight);

      // Fill the histograms with the information from the approriate electron
      if (zfe.e0 != NULL && zfe.e1 != NULL){
        e0_pt_->Fill(zfe.e0->pt, event_weight);
        e0_eta_->Fill(zfe.e0->eta, event_weight);
        e0_phi_->Fill(zfe.e0->phi, event_weight);
        e0_charge_->Fill(zfe.e0->charge, event_weight);
        e1_pt_->Fill(zfe.e1->pt, event_weight);
        e1_eta_->Fill(zfe.e1->eta, event_weight);
        e1_phi_->Fill(zfe.e1->phi, event_weight);
        e1_charge_->Fill(zfe.e1->charge, event_weight);
      }

      //Z To mumu
      z_from_muons_mass_all_->Fill(zfe.reco_z_from_muons.m, event_weight);
      z_from_muons_mass_coarse_->Fill(zfe.reco_z_from_muons.m, event_weight);
      z_from_muons_mass_fine_->Fill(zfe.reco_z_from_muons.m, event_weight);
      z_from_muons_rapidity_->Fill(zfe.reco_z_from_muons.y, event_weight);
      z_from_muons_pt_->Fill(zfe.reco_z_from_muons.pt, event_weight);
      z_from_muons_vtx_x_->Fill(zfe.reco_z_from_muons.vtx_x, event_weight);
      z_from_muons_vtx_y_->Fill(zfe.reco_z_from_muons.vtx_y, event_weight);
      z_from_muons_vtx_z_->Fill(zfe.reco_z_from_muons.vtx_z, event_weight);
      z_from_muons_vtx_prob_->Fill(zfe.reco_z_from_muons.vtx_prob, event_weight);
      z_from_muons_phistar_->Fill(zfe.reco_z_from_muons.phistar, event_weight);

      if (zfe.reco_z_from_muons.m > -1) {
        muon0_from_z_pt_->Fill(zfe.z_muon0.pt(), event_weight);
        muon0_from_z_eta_->Fill(zfe.z_muon0.eta(), event_weight);
        muon0_from_z_phi_->Fill(zfe.z_muon0.phi(), event_weight);
        muon0_from_z_charge_->Fill(zfe.z_muon0.charge(), event_weight);
        muon1_from_z_pt_->Fill(zfe.z_muon1.pt(), event_weight);
        muon1_from_z_eta_->Fill(zfe.z_muon1.eta(), event_weight);
        muon1_from_z_phi_->Fill(zfe.z_muon1.phi(), event_weight);
        muon1_from_z_charge_->Fill(zfe.z_muon1.charge(), event_weight);
      }

      primary_vtx_x_minus_z_vtx_x_->Fill( zfe.reco_vert.primary_x - zfe.reco_z.vtx_x, event_weight );
      primary_vtx_y_minus_z_vtx_y_->Fill( zfe.reco_vert.primary_y - zfe.reco_z.vtx_y, event_weight );
      primary_vtx_z_minus_z_vtx_z_->Fill( zfe.reco_vert.primary_z - zfe.reco_z.vtx_z, event_weight );

      primary_vtx_x_minus_zmumu_vtx_x_->Fill( zfe.reco_vert.primary_x - zfe.reco_z_from_muons.vtx_x, event_weight );
      primary_vtx_y_minus_zmumu_vtx_y_->Fill( zfe.reco_vert.primary_y - zfe.reco_z_from_muons.vtx_y, event_weight );
      primary_vtx_z_minus_zmumu_vtx_z_->Fill( zfe.reco_vert.primary_z - zfe.reco_z_from_muons.vtx_z, event_weight );

      primary_vtx_x_vs_z_vtx_x_->Fill(zfe.reco_z.vtx_x, zfe.reco_vert.primary_x, event_weight );
      primary_vtx_y_vs_z_vtx_y_->Fill(zfe.reco_z.vtx_y, zfe.reco_vert.primary_y, event_weight );
      primary_vtx_z_vs_z_vtx_z_->Fill(zfe.reco_z.vtx_z, zfe.reco_vert.primary_z, event_weight );

      primary_vtx_x_vs_zmuons_vtx_x_->Fill(zfe.reco_z_from_muons.vtx_x, zfe.reco_vert.primary_x, event_weight );
      primary_vtx_y_vs_zmuons_vtx_y_->Fill(zfe.reco_z_from_muons.vtx_y, zfe.reco_vert.primary_y, event_weight );
      primary_vtx_z_vs_zmuons_vtx_z_->Fill(zfe.reco_z_from_muons.vtx_z, zfe.reco_vert.primary_z, event_weight );

      int n_jpsi = 0;
      for (unsigned int i = 0; i < zfe.reco_jpsi.m.size() ; ++i ) {
        //TODO should fold in eta window cut here??
        //TODO should fold in jpsi pT window cut here??
        if (APPLY_MUON_MIN_PT_ && (!zfe.reco_jpsi.has_high_pt_muons.at(i) || !zfe.reco_jpsi.has_muons_in_eta_window.at(i) || 
            !zfe.reco_jpsi.is_high_pt.at(i) ) )  {
          continue;
        }
        if (APPLY_SOFT_MUONS_ && !zfe.reco_jpsi.has_soft_id_muons.at(i) ) {
          continue;
        }
        if (APPLY_JPSI_MASS_WINDOW_ && !zfe.reco_jpsi.is_within_jpsi_mass_window.at(i) ) {
          continue;
        }
        //TODO should this cut be relative to the z position??
        if (APPLY_VERTEX_Z_POS_WINDOW_ && !zfe.reco_jpsi.has_dimuon_vertex_compatible_with_primary_vertex.at(i) ) {
          continue;
        }
        if (APPLY_DIMUON_VTX_COMPATIBILITY_ && !zfe.reco_jpsi.has_muons_with_compatible_vertex.at(i) ) {
          continue;
        }
        n_jpsi++;

        jpsi_mass_all_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        jpsi_mass_coarse_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        jpsi_mass_fine_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        jpsi_four_lepton_mass_->Fill(zfe.reco_jpsi.four_lepton_mass.at(i), event_weight);
        jpsi_rapidity_->Fill(zfe.reco_jpsi.y.at(i), event_weight);
        jpsi_pt_->Fill(zfe.reco_jpsi.pt.at(i), event_weight);
        jpsi_efficiency_->Fill(zfe.reco_jpsi.jpsi_efficiency.at(i), event_weight);
        jpsi_acc_eff_->Fill(zfe.reco_jpsi.jpsi_acc_eff.at(i), event_weight);
        jpsi_scale_factor_->Fill(zfe.reco_jpsi.jpsi_scale_factor.at(i), event_weight);
        jpsi_trigger_obj_mu0_pt_->Fill(zfe.reco_jpsi.trigger_object_mu0_pt.at(i), event_weight);
        jpsi_trigger_obj_mu1_pt_->Fill(zfe.reco_jpsi.trigger_object_mu1_pt.at(i), event_weight);
        jpsi_trigger_obj_mu0_pt_minus_reco_mu0_pt_->Fill(zfe.reco_jpsi.trigger_object_mu0_pt.at(i) - zfe.reco_jpsi.muon0.at(i).pt() , event_weight);
        jpsi_trigger_obj_mu1_pt_minus_reco_mu1_pt_->Fill(zfe.reco_jpsi.trigger_object_mu1_pt.at(i) - zfe.reco_jpsi.muon1.at(i).pt() , event_weight);

        jpsi_cos_mu_plus_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight);
        double lambdaNeg1 = (1.0 - pow(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), 2.0) );
        //divide by 2.0 for lambda positive 1 to keep weight less than 1
        double lambda1 = (1.0 + pow(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), 2.0)) / 2.0;
        jpsi_cos_mu_plus_lambdaNeg1_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambdaNeg1);
        jpsi_cos_mu_plus_lambda1_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambda1);
        if (zfe.reco_jpsi.pt.at(i) >= 8.0 && zfe.reco_jpsi.pt.at(i) < 8.5 ) {
          jpsi_cos_mu_plus_jpsi_pt_8to8p5_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight);
          jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8to8p5_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambdaNeg1);
          jpsi_cos_mu_plus_lambda1_jpsi_pt_8to8p5_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambda1);
        }
        if (zfe.reco_jpsi.pt.at(i) >= 8.5 && zfe.reco_jpsi.pt.at(i) < 9.0 ) {
          jpsi_cos_mu_plus_jpsi_pt_8p5to9_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight);
          jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_8p5to9_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambdaNeg1);
          jpsi_cos_mu_plus_lambda1_jpsi_pt_8p5to9_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambda1);
        }
        if (zfe.reco_jpsi.pt.at(i) >= 9.0 && zfe.reco_jpsi.pt.at(i) < 10.0 ) {
          jpsi_cos_mu_plus_jpsi_pt_9to10_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight);
          jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_9to10_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambdaNeg1);
          jpsi_cos_mu_plus_lambda1_jpsi_pt_9to10_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambda1);
        }
        if (zfe.reco_jpsi.pt.at(i) >= 10.0 && zfe.reco_jpsi.pt.at(i) < 15.0 ) {
          jpsi_cos_mu_plus_jpsi_pt_10to15_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight);
          jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_10to15_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambdaNeg1);
          jpsi_cos_mu_plus_lambda1_jpsi_pt_10to15_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambda1);
        }
        if (zfe.reco_jpsi.pt.at(i) >= 15.0 && zfe.reco_jpsi.pt.at(i) < 20.0 ) {
          jpsi_cos_mu_plus_jpsi_pt_15to20_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight);
          jpsi_cos_mu_plus_lambdaNeg1_jpsi_pt_15to20_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambdaNeg1);
          jpsi_cos_mu_plus_lambda1_jpsi_pt_15to20_->Fill(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), event_weight * lambda1);
        }

        jpsi_cos_mu_minus_->Fill(zfe.reco_jpsi.cos_jpsi_mu_minus.at(i), event_weight);
        if (!zfe.is_real_data) {
          for (unsigned int j = 0; j < zfe.truth_jpsi.m.size() ; ++j ) {
            jpsi_reco_pt_vs_jpsi_truth_pt_->Fill(zfe.truth_jpsi.pt.at(j), zfe.reco_jpsi.pt.at(i), event_weight); 
            jpsi_truth_pt_minus_jpsi_reco_pt_->Fill(zfe.truth_jpsi.pt.at(j) - zfe.reco_jpsi.pt.at(i), event_weight);
            jpsi_truth_vtx_x_minus_jpsi_reco_vtx_x_->Fill(zfe.truth_jpsi.vtx_x.at(j) - zfe.reco_jpsi.vtx_x.at(i), event_weight);
            jpsi_truth_vtx_y_minus_jpsi_reco_vtx_y_->Fill(zfe.truth_jpsi.vtx_y.at(j) - zfe.reco_jpsi.vtx_y.at(i), event_weight);
            jpsi_truth_vtx_z_minus_jpsi_reco_vtx_z_->Fill(zfe.truth_jpsi.vtx_z.at(j) - zfe.reco_jpsi.vtx_z.at(i), event_weight);
          }
        }
        jpsi_pt_vs_rap_->Fill(zfe.reco_jpsi.y.at(i), zfe.reco_jpsi.pt.at(i) , event_weight);
        jpsi_pt_vs_rap_polarization_long->Fill(zfe.reco_jpsi.y.at(i), zfe.reco_jpsi.pt.at(i) , 
            event_weight * (1.0 - pow(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), 2.0) ));
        jpsi_pt_vs_rap_polarization_TPlusZero->Fill(zfe.reco_jpsi.y.at(i), zfe.reco_jpsi.pt.at(i) , 
            event_weight * (1.0 + pow(zfe.reco_jpsi.cos_jpsi_mu_plus.at(i), 2.0) ));
        jpsi_vtx_distance_z_vtx_x_->Fill( zfe.reco_jpsi.distance_x.at(i), event_weight);
        jpsi_vtx_distance_z_vtx_y_->Fill( zfe.reco_jpsi.distance_y.at(i), event_weight);
        jpsi_vtx_distance_z_vtx_z_->Fill( zfe.reco_jpsi.distance_z.at(i), event_weight);
        jpsi_distance_->Fill(zfe.reco_jpsi.distance.at(i), event_weight);
        jpsi_dist_err_->Fill(zfe.reco_jpsi.dist_err.at(i), event_weight);
        jpsi_chi2_->Fill(zfe.reco_jpsi.chi2.at(i), event_weight);
        jpsi_distance_xy_->Fill(zfe.reco_jpsi.distance_xy.at(i), event_weight);
        jpsi_dist_err_xy_->Fill(zfe.reco_jpsi.dist_err_xy.at(i), event_weight);
        jpsi_chi2_xy_->Fill(zfe.reco_jpsi.chi2_xy.at(i), event_weight);
        jpsi_zpt_difference_->Fill( zfe.reco_z.pt - zfe.reco_jpsi.pt.at(i), event_weight);
        jpsi_zmumupt_difference_->Fill( zfe.reco_z_from_muons.pt - zfe.reco_jpsi.pt.at(i), event_weight);
        jpsi_tau_xy_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        jpsi_tau_xy_fine_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        jpsi_tau_xy_very_fine_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps

        //pt
        if ( zfe.reco_jpsi.pt.at(i) < 10.0 ) {
          jpsi_tau_xy_very_fine_ptUnder10_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
          jpsi_mass_fine_ptUnder10_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        }
        if ( zfe.reco_jpsi.pt.at(i) >= 10.0 && zfe.reco_jpsi.pt.at(i) < 15 ) {
          jpsi_tau_xy_very_fine_pt10to15_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
          jpsi_mass_fine_pt10to15_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        }
        if ( zfe.reco_jpsi.pt.at(i) >= 15.0 && zfe.reco_jpsi.pt.at(i) < 20 ) {
          jpsi_tau_xy_very_fine_pt15to20_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
          jpsi_mass_fine_pt15to20_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        }
        if ( zfe.reco_jpsi.pt.at(i) >= 20.0 && zfe.reco_jpsi.pt.at(i) < 25 ) {
          jpsi_tau_xy_very_fine_pt20to25_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
          jpsi_mass_fine_pt20to25_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        }
        if ( zfe.reco_jpsi.pt.at(i) >= 25.0 && zfe.reco_jpsi.pt.at(i) < 30 ) {
          jpsi_tau_xy_very_fine_pt25to30_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
          jpsi_mass_fine_pt25to30_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        }
        if ( zfe.reco_jpsi.pt.at(i) >= 30.0 ) {
          jpsi_tau_xy_very_fine_ptAbove30_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
          jpsi_mass_fine_ptAbove30_->Fill(zfe.reco_jpsi.m.at(i), event_weight);
        }

        //rap
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 0.0 && fabs(zfe.reco_jpsi.y.at(i)) < 0.3 ) {
          jpsi_tau_xy_very_fine_rap0_0to0_3_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 0.3 && fabs(zfe.reco_jpsi.y.at(i)) < 0.6 ) {
          jpsi_tau_xy_very_fine_rap0_3to0_6_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 0.6 && fabs(zfe.reco_jpsi.y.at(i)) < 0.9 ) {
          jpsi_tau_xy_very_fine_rap0_6to0_9_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 0.9 && fabs(zfe.reco_jpsi.y.at(i)) < 1.2 ) {
          jpsi_tau_xy_very_fine_rap0_9to1_2_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 1.2 && fabs(zfe.reco_jpsi.y.at(i)) < 1.5 ) {
          jpsi_tau_xy_very_fine_rap1_2to1_5_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 1.5 && fabs(zfe.reco_jpsi.y.at(i)) < 1.8 ) {
          jpsi_tau_xy_very_fine_rap1_5to1_8_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 1.8 && fabs(zfe.reco_jpsi.y.at(i)) < 2.1 ) {
          jpsi_tau_xy_very_fine_rap1_8to2_1_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 2.1 && fabs(zfe.reco_jpsi.y.at(i)) < 2.4 ) {
          jpsi_tau_xy_very_fine_rap2_1to2_4_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }

        //rap and pt
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 0.0 && fabs(zfe.reco_jpsi.y.at(i)) < 0.9 && (zfe.reco_jpsi.pt.at(i) >= 10.0 && zfe.reco_jpsi.pt.at(i) < 15 ) ) {
          jpsi_tau_xy_very_fine_rap0_0to0_9pt10to15_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 0.9 && fabs(zfe.reco_jpsi.y.at(i)) < 1.2 && (zfe.reco_jpsi.pt.at(i) >= 10.0 && zfe.reco_jpsi.pt.at(i) < 15 ) ) {
          jpsi_tau_xy_very_fine_rap0_9to1_2pt10to15_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        if ( fabs(zfe.reco_jpsi.y.at(i)) >= 1.2 && fabs(zfe.reco_jpsi.y.at(i)) < 2.1 && (zfe.reco_jpsi.pt.at(i) >= 10.0 && zfe.reco_jpsi.pt.at(i) < 15 ) ) {
          jpsi_tau_xy_very_fine_rap1_2to2_1pt10to15_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }

        //similar pt muons
        if ( fabs( zfe.reco_jpsi.muon0.at(i).pt() - zfe.reco_jpsi.muon1.at(i).pt() ) < 2 ) {
          jpsi_tau_xy_very_fine_similar_pt_muons_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }
        //continuum bg
        if ( (zfe.reco_jpsi.m.at(i) < MIN_JPSI_MASS && zfe.reco_jpsi.m.at(i) > 2.5) ||
            (zfe.reco_jpsi.m.at(i) < 5.0 && zfe.reco_jpsi.m.at(i) > MAX_JPSI_MASS) ) {
          jpsi_tau_xy_dimuon_continuum_bg_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }

        //number of muon segments
        bool mu0ID = muon::isGoodMuon(zfe.reco_jpsi.muon0.at(i), muon::TMOneStationTight);
        int mu0_layers = 0;
        if (mu0ID) {
          mu0_layers = zfe.reco_jpsi.muon0.at(i).innerTrack()->hitPattern().trackerLayersWithMeasurement();
        }
        bool mu1ID = muon::isGoodMuon(zfe.reco_jpsi.muon1.at(i), muon::TMOneStationTight);
        int mu1_layers = 0;
        if (mu1ID) {
          mu1_layers = zfe.reco_jpsi.muon1.at(i).innerTrack()->hitPattern().trackerLayersWithMeasurement();
        }
        if ( mu0_layers > 12 && mu1_layers > 12 ) {
          jpsi_tau_xy_very_fine_above_12_tracker_layers_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        }

        jpsi_tau_z_->Fill(zfe.reco_jpsi.tau_z.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        jpsi_tau_z_fine_->Fill(zfe.reco_jpsi.tau_z.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        jpsi_tau_z_very_fine_->Fill(zfe.reco_jpsi.tau_z.at(i) * 1000, event_weight); // multiply by 1000 to go from ns to ps
        jpsi_mass_vs_chi2_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.chi2.at(i), event_weight );

        dimuon_pt_vs_zee_pt_->Fill(zfe.reco_jpsi.pt.at(i) , zfe.reco_z.pt, event_weight );
        dimuon_pt_vs_zmumu_pt_->Fill(zfe.reco_jpsi.pt.at(i) , zfe.reco_z_from_muons.pt, event_weight );

        dimuon_mass_vs_dimuon_tau_xy_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
        dimuon_mass_vs_dimuon_tau_xy_fine_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
        dimuon_mass_vs_dimuon_tau_xy_fine_weighted_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, 
            event_weight / zfe.reco_jpsi.jpsi_acc_eff.at(i) );
        
        if (zfe.reco_jpsi.pt.at(i) >= 8.5 && zfe.reco_jpsi.pt.at(i) < 10.0 ) {
          dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
        }
        if (zfe.reco_jpsi.pt.at(i) >= 10.0 && zfe.reco_jpsi.pt.at(i) < 14.0 ) {
          dimuon_mass_vs_dimuon_tau_xy_10to14_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
        }
        if (zfe.reco_jpsi.pt.at(i) >= 14.0 && zfe.reco_jpsi.pt.at(i) < 18.0 ) {
          dimuon_mass_vs_dimuon_tau_xy_14to18_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
        }
        if (zfe.reco_jpsi.pt.at(i) >= 18.0 && zfe.reco_jpsi.pt.at(i) < 30.0 ) {
          dimuon_mass_vs_dimuon_tau_xy_18to30_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
        }
        if (zfe.reco_jpsi.pt.at(i) >= 30.0 && zfe.reco_jpsi.pt.at(i) < 100.0 ) {
          dimuon_mass_vs_dimuon_tau_xy_30to100_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
        }

        if(zfe.reco_jpsi.y.at(i) < 1.2) {
          low_rap_dimuon_mass_vs_dimuon_tau_xy_fine_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          if (zfe.reco_jpsi.pt.at(i) >= 8.5 && zfe.reco_jpsi.pt.at(i) < 10.0 ) {
            low_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
          if (zfe.reco_jpsi.pt.at(i) >= 10.0 && zfe.reco_jpsi.pt.at(i) < 14.0 ) {
            low_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
          if (zfe.reco_jpsi.pt.at(i) >= 14.0 && zfe.reco_jpsi.pt.at(i) < 18.0 ) {
            low_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
          if (zfe.reco_jpsi.pt.at(i) >= 18.0 && zfe.reco_jpsi.pt.at(i) < 30.0 ) {
            low_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
          if (zfe.reco_jpsi.pt.at(i) >= 30.0 && zfe.reco_jpsi.pt.at(i) < 100.0 ) {
            low_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
        }
        else {
          high_rap_dimuon_mass_vs_dimuon_tau_xy_fine_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          if (zfe.reco_jpsi.pt.at(i) >= 8.5 && zfe.reco_jpsi.pt.at(i) < 10.0 ) {
            high_rap_dimuon_mass_vs_dimuon_tau_xy_8p5to10p0_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
          if (zfe.reco_jpsi.pt.at(i) >= 10.0 && zfe.reco_jpsi.pt.at(i) < 14.0 ) {
            high_rap_dimuon_mass_vs_dimuon_tau_xy_10to14_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
          if (zfe.reco_jpsi.pt.at(i) >= 14.0 && zfe.reco_jpsi.pt.at(i) < 18.0 ) {
            high_rap_dimuon_mass_vs_dimuon_tau_xy_14to18_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
          if (zfe.reco_jpsi.pt.at(i) >= 18.0 && zfe.reco_jpsi.pt.at(i) < 30.0 ) {
            high_rap_dimuon_mass_vs_dimuon_tau_xy_18to30_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
          if (zfe.reco_jpsi.pt.at(i) >= 30.0 && zfe.reco_jpsi.pt.at(i) < 100.0 ) {
            high_rap_dimuon_mass_vs_dimuon_tau_xy_30to100_->Fill(zfe.reco_jpsi.m.at(i) , zfe.reco_jpsi.tau_xy.at(i) * 1000, event_weight );
          }
        }



        jpsi_tau_xy_vs_tau_z_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000 , zfe.reco_jpsi.tau_z.at(i) * 1000, event_weight );
        jpsi_tau_xy_vs_distance_z_->Fill(zfe.reco_jpsi.tau_xy.at(i) * 1000 , (zfe.reco_jpsi.vtx_z.at(i) - zfe.reco_z.vtx_z), event_weight );
        jpsi_tau_z_vs_distance_z_->Fill(zfe.reco_jpsi.tau_z.at(i) * 1000 , (zfe.reco_jpsi.vtx_z.at(i) - zfe.reco_z.vtx_z), event_weight );

        dimuon_vtx_x_->Fill(zfe.reco_jpsi.vtx_x.at(i), event_weight);
        dimuon_vtx_y_->Fill(zfe.reco_jpsi.vtx_y.at(i), event_weight);
        dimuon_vtx_z_->Fill(zfe.reco_jpsi.vtx_z.at(i), event_weight);
        dimuon_vtx_prob_->Fill(zfe.reco_jpsi.vtx_prob.at(i), event_weight);

        z_jpsi_delta_phi_->Fill(zfe.reco_jpsi.z_delta_phi.at(i), event_weight);

        dimuon_delta_phi_->Fill(zfe.reco_jpsi.muons_delta_phi.at(i), event_weight);
        dimuon_deltaR_->Fill(zfe.reco_jpsi.muons_deltaR.at(i), event_weight);
        dimuon_delta_eta_->Fill(zfe.reco_jpsi.muons_delta_eta.at(i), event_weight);

        mu0_pt_->Fill(zfe.reco_jpsi.muon0.at(i).pt(), event_weight);
        mu0_eta_->Fill(zfe.reco_jpsi.muon0.at(i).eta(), event_weight);
        mu0_efficiency_->Fill(zfe.reco_jpsi.muon0_efficiency.at(i), event_weight);
        mu0_scale_factor_->Fill(zfe.reco_jpsi.muon0_scale_factor.at(i), event_weight);
        mu0_phi_->Fill(zfe.reco_jpsi.muon0.at(i).phi(), event_weight);
        mu0_charge_->Fill(zfe.reco_jpsi.muon0.at(i).charge(), event_weight);
        mu0_deltaR_z_muon_->Fill(zfe.reco_jpsi.muon0_deltaR_to_z_muons.at(i), event_weight);
        mu0_deltaR_truth_->Fill(zfe.reco_jpsi.muon0_deltaR_to_truth_muons.at(i), event_weight);
        mu0_tracker_layers_->Fill(mu0_layers, event_weight);

        mu1_pt_->Fill(zfe.reco_jpsi.muon1.at(i).pt(), event_weight);
        mu1_eta_->Fill(zfe.reco_jpsi.muon1.at(i).eta(), event_weight);
        mu1_efficiency_->Fill(zfe.reco_jpsi.muon1_efficiency.at(i), event_weight);
        mu1_scale_factor_->Fill(zfe.reco_jpsi.muon1_scale_factor.at(i), event_weight);
        mu1_phi_->Fill(zfe.reco_jpsi.muon1.at(i).phi(), event_weight);
        mu1_charge_->Fill(zfe.reco_jpsi.muon1.at(i).charge(), event_weight);
        mu1_deltaR_z_muon_->Fill(zfe.reco_jpsi.muon1_deltaR_to_z_muons.at(i), event_weight);
        mu1_deltaR_truth_->Fill(zfe.reco_jpsi.muon1_deltaR_to_truth_muons.at(i), event_weight);
        mu1_tracker_layers_->Fill(mu1_layers, event_weight);

        int n_truth_matched_jpsi_muons = 0;
        if ( zfe.reco_jpsi.muon0_deltaR_to_truth_muons.at(i) <= MAX_DELTAR_TRUTH_MATCHED_JPSI_MUONS ) {
          n_truth_matched_jpsi_muons++;
        }
        if ( zfe.reco_jpsi.muon1_deltaR_to_truth_muons.at(i) <= MAX_DELTAR_TRUTH_MATCHED_JPSI_MUONS ) {
          n_truth_matched_jpsi_muons++;
        }
        n_truth_matched_jpsi_muons_->Fill (n_truth_matched_jpsi_muons, event_weight);

        jpsi_iso_mu0_->Fill(zfe.reco_jpsi.iso_mu0.at(i), event_weight ) ;
        jpsi_iso_sum_charged_hadron_pt_mu0_->Fill(zfe.reco_jpsi.iso_sum_charged_hadron_pt_mu0.at(i), event_weight ) ;
        jpsi_iso_sum_charged_particle_pt_mu0_->Fill(zfe.reco_jpsi.iso_sum_charged_particle_pt_mu0.at(i), event_weight ) ;
        jpsi_iso_sum_neutral_hadron_et_mu0_->Fill(zfe.reco_jpsi.iso_sum_neutral_hadron_et_mu0.at(i), event_weight ) ;
        jpsi_iso_sum_photon_et_mu0_->Fill(zfe.reco_jpsi.iso_sum_photon_et_mu0.at(i), event_weight ) ;
        jpsi_iso_sum_pileup_pt_mu0_->Fill(zfe.reco_jpsi.iso_sum_pileup_pt_mu0.at(i), event_weight ) ;

        jpsi_iso_mu1_->Fill(zfe.reco_jpsi.iso_mu1.at(i), event_weight ) ;
        jpsi_iso_sum_charged_hadron_pt_mu1_->Fill(zfe.reco_jpsi.iso_sum_charged_hadron_pt_mu1.at(i), event_weight ) ;
        jpsi_iso_sum_charged_particle_pt_mu1_->Fill(zfe.reco_jpsi.iso_sum_charged_particle_pt_mu1.at(i), event_weight ) ;
        jpsi_iso_sum_neutral_hadron_et_mu1_->Fill(zfe.reco_jpsi.iso_sum_neutral_hadron_et_mu1.at(i), event_weight ) ;
        jpsi_iso_sum_photon_et_mu1_->Fill(zfe.reco_jpsi.iso_sum_photon_et_mu1.at(i), event_weight ) ;
        jpsi_iso_sum_pileup_pt_mu1_->Fill(zfe.reco_jpsi.iso_sum_pileup_pt_mu1.at(i), event_weight ) ;
      }
      for (unsigned int i = 0; i < zfe.reco_jets.pt.size() ; ++i ) {
        jet_pt_->Fill(zfe.reco_jets.pt.at(i), event_weight);
        jet_eta_->Fill(zfe.reco_jets.eta.at(i), event_weight);
        jet_btag_discriminator_->Fill(zfe.reco_jets.btag_discriminator.at(i), event_weight);
      }
      for (unsigned int i = 0; i < zfe.reco_muon_jets.pt.size() ; ++i ) {
        muon_jet_pt_->Fill(zfe.reco_muon_jets.pt.at(i), event_weight);
        muon_jet_pt_diff_z_pt_->Fill(zfe.reco_z.pt - zfe.reco_muon_jets.pt.at(i), event_weight);
        muon_jet_pt_z_pt_->Fill(zfe.reco_z.pt , zfe.reco_muon_jets.pt.at(i), event_weight);
        muon_jet_phi_z_phi_->Fill(zfe.reco_z.phi , zfe.reco_muon_jets.phi.at(i), event_weight);
        //TODO is this extra for loop needed?
        for (unsigned int j = 0; j < zfe.reco_jpsi.m.size() ; ++j ) {
          //TODO should fold in eta window cut here??
          if (APPLY_MUON_MIN_PT_ && (!zfe.reco_jpsi.has_high_pt_muons.at(j) || !zfe.reco_jpsi.has_muons_in_eta_window.at(j) 
                || !zfe.reco_jpsi.is_high_pt.at(j)) )  {
            continue;
          }
          if (APPLY_SOFT_MUONS_ && !zfe.reco_jpsi.has_soft_id_muons.at(j) ) {
            continue;
          }
          if (APPLY_DIMUON_VTX_COMPATIBILITY_ && !zfe.reco_jpsi.has_muons_with_compatible_vertex.at(j) ) {
            continue;
          }
          if (APPLY_VERTEX_Z_POS_WINDOW_ && !zfe.reco_jpsi.has_dimuon_vertex_compatible_with_primary_vertex.at(j) ) {
            continue;
          }
          if (APPLY_JPSI_MASS_WINDOW_ && !zfe.reco_jpsi.is_within_jpsi_mass_window.at(j) ) {
            continue;
          }
          muon_jet_pt_diff_dimuon_pt_->Fill(zfe.reco_jpsi.pt.at(j) - zfe.reco_muon_jets.pt.at(i), event_weight);
          muon_jet_pt_dimuon_pt_->Fill(zfe.reco_jpsi.pt.at(j) , zfe.reco_muon_jets.pt.at(i), event_weight);
          muon_jet_phi_dimuon_phi_->Fill(zfe.reco_jpsi.phi.at(j) , zfe.reco_muon_jets.phi.at(i), event_weight);
        }
        muon_jet_eta_->Fill(zfe.reco_muon_jets.eta.at(i), event_weight);
        muon_jet_btag_discriminator_->Fill(zfe.reco_muon_jets.btag_discriminator.at(i), event_weight);
      }
      for (unsigned int i = 0; i < zfe.reco_vert.x.size() ; ++i ) {
        vtx_x_->Fill(zfe.reco_vert.x.at(i), event_weight);
        vtx_y_->Fill(zfe.reco_vert.y.at(i), event_weight);
        vtx_z_->Fill(zfe.reco_vert.z.at(i), event_weight);
      }


      // Event Info
      primary_vtx_x_->Fill(zfe.reco_vert.primary_x, event_weight);
      primary_vtx_y_->Fill(zfe.reco_vert.primary_y, event_weight);
      primary_vtx_z_->Fill(zfe.reco_vert.primary_z, event_weight);

      baseweights_->Fill(event_weight); //Don't weight this histogram of event weights

      pileup_->Fill(zfe.reco_vert.num, event_weight);
      nelectrons_->Fill(zfe.n_reco_electrons, event_weight);
      nmuons_->Fill(zfe.n_reco_muons, event_weight);
      njets_->Fill(zfe.n_reco_jets, event_weight);
      n_muonjets_->Fill(zfe.n_reco_muon_jets, event_weight);
      njpsis_->Fill(n_jpsi, event_weight); 
    } 
    else if (USE_MC_ && !zfe.is_real_data) {
      z_mass_all_->Fill(zfe.truth_z.m, event_weight);
      z_mass_coarse_->Fill(zfe.truth_z.m, event_weight);
      z_mass_fine_->Fill(zfe.truth_z.m, event_weight);
      z_rapidity_->Fill(zfe.truth_z.y, event_weight);
      z_pt_->Fill(zfe.truth_z.pt, event_weight);
      phistar_->Fill(zfe.truth_z.phistar, event_weight);

      // Fill the histograms with the information from the approriate electron
      if (zfe.e0_truth != NULL && zfe.e1_truth != NULL) {
        e0_pt_->Fill(zfe.e0_truth->pt, event_weight);
        e0_eta_->Fill(zfe.e0_truth->eta, event_weight);
        e0_phi_->Fill(zfe.e0_truth->phi, event_weight);
        e1_pt_->Fill(zfe.e1_truth->pt, event_weight);
        e1_eta_->Fill(zfe.e1_truth->eta, event_weight);
        e1_phi_->Fill(zfe.e1_truth->phi, event_weight);
      }
      //jpsi plots
      int n_truth_jpsi = 0;
      for (unsigned int i = 0; i < zfe.truth_jpsi.m.size() ; ++i ) {
        if (APPLY_JPSI_MASS_WINDOW_ && !zfe.truth_jpsi.is_within_jpsi_mass_window.at(i) ) {
          continue;
        }
        //TODO fix this
        //if (APPLY_VERTEX_Z_POS_WINDOW_ && fabs(zfe.truth_jpsi.distance_z.at(i)) > MAX_JPSI_VERTEX_Z_DISPLACEMENT ) {
        //  continue;
        //}
        //TODO eta cut needed here?? jpsi pT cut needed here?? || !zfe.truth_jpsi.is_high_pt.at(i) not really needed as included in ZFinderEvent for mc change this??
        if (APPLY_MUON_MIN_PT_ && (!zfe.truth_jpsi.has_high_pt_muons.at(i) || !zfe.truth_jpsi.has_muons_in_eta_window.at(i) ) ) {
          continue;
        }

        n_truth_jpsi++;

        jpsi_mass_all_->Fill(zfe.truth_jpsi.m.at(i), event_weight);
        jpsi_mass_coarse_->Fill(zfe.truth_jpsi.m.at(i), event_weight);
        jpsi_mass_fine_->Fill(zfe.truth_jpsi.m.at(i), event_weight);
        jpsi_rapidity_->Fill(zfe.truth_jpsi.y.at(i), event_weight);
        jpsi_pt_->Fill(zfe.truth_jpsi.pt.at(i), event_weight);
        jpsi_pt_vs_rap_->Fill(zfe.truth_jpsi.y.at(i), zfe.truth_jpsi.pt.at(i) , event_weight);
        jpsi_cos_mu_plus_->Fill(zfe.truth_jpsi.cos_jpsi_mu_plus.at(i), event_weight);
        jpsi_cos_mu_minus_->Fill(zfe.truth_jpsi.cos_jpsi_mu_minus.at(i), event_weight);

        jpsi_pt_vs_rap_polarization_long->Fill(zfe.truth_jpsi.y.at(i), zfe.truth_jpsi.pt.at(i) , 
            event_weight * (1.0 - pow(zfe.truth_jpsi.cos_jpsi_mu_plus.at(i), 2.0) ));
        jpsi_pt_vs_rap_polarization_TPlusZero->Fill(zfe.truth_jpsi.y.at(i), zfe.truth_jpsi.pt.at(i) , 
            //divide by 2.0 as Yuichi suggested to avoid weights greater than 1
            event_weight * (1.0 + pow(zfe.truth_jpsi.cos_jpsi_mu_plus.at(i), 2.0) ) / 2.0);

        mu0_pt_->Fill(zfe.jpsi_muon0.at(i)->pt(), event_weight);
        mu0_eta_->Fill(zfe.jpsi_muon0.at(i)->eta(), event_weight);
        mu0_phi_->Fill(zfe.jpsi_muon0.at(i)->phi(), event_weight);
        mu0_charge_->Fill(zfe.jpsi_muon0.at(i)->charge(), event_weight);
        mu1_pt_->Fill(zfe.jpsi_muon1.at(i)->pt(), event_weight);
        mu1_eta_->Fill(zfe.jpsi_muon1.at(i)->eta(), event_weight);
        mu1_phi_->Fill(zfe.jpsi_muon1.at(i)->phi(), event_weight);
        mu1_charge_->Fill(zfe.jpsi_muon1.at(i)->charge(), event_weight);
      }
      // Event Info
      baseweights_->Fill(event_weight); //Don't weight this histogram of event weights
      pileup_->Fill(zfe.truth_vert.num, event_weight);
      nelectrons_->Fill(2, event_weight);  // We only ever grab the two electrons from the Z
      njpsis_->Fill(n_truth_jpsi, event_weight); 
    }
  }

  //void ZFinderPlotter::CalculateJpsiLifetime(const ZFinderEvent::JPsiData &jpsi_data, const ZFinderEvent::ZFromMuonsData &z_from_muon ) {
  //  for (unsigned int i = 0; i < jpsi_data.m.size() ; ++i ) {
  //  }
  //}
}  // namespace zf
