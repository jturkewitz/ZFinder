#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
void make_images ( )
{
  
  TFile *theFile0 = new TFile("/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/test42.root");

  //TODO put this in rootrc
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();

  TH1D *h_z0_mass_all = (TH1D*) theFile0->Get("Z0 Mass: All");
  TH1D *h_z0_mass_coarse = (TH1D*) theFile0->Get("Z0 Mass: Coarse");
  TH1D *h_z0_mass_fine = (TH1D*) theFile0->Get("Z0 Mass: Fine");
  TH1D *h_z0_rapidity = (TH1D*) theFile0->Get("Z0 Rapidity");
  TH1D *h_z0_pT = (TH1D*) theFile0->Get("Z0 p_{T}");
  TH1D *h_dielectron_vertex_probability = (TH1D*) theFile0->Get("dielectron vertex probability");
  TH1D *h_e0_pT = (TH1D*) theFile0->Get("p_{T,e_{0}}");
  TH1D *h_e1_pT = (TH1D*) theFile0->Get("p_{T,e_{1}}");
  TH1D *h_e0_eta = (TH1D*) theFile0->Get("#eta_{e_{0}}");
  TH1D *h_e1_eta = (TH1D*) theFile0->Get("#eta_{e_{1}}");
  TH1D *h_e0_phi = (TH1D*) theFile0->Get("#phi_{e_{0}}");
  TH1D *h_e1_phi = (TH1D*) theFile0->Get("#phi_{e_{1}}");
  TH1D *h_e0_charge = (TH1D*) theFile0->Get("charge_{e_{0}}");
  TH1D *h_e1_charge = (TH1D*) theFile0->Get("charge_{e_{1}}");
  TH1D *h_z0_phiStar = (TH1D*) theFile0->Get("#phi*");
  TH1D *h_jpsi_mass_all = (TH1D*) theFile0->Get("jpsi0 Mass: All");
  TH1D *h_jpsi_mass_coarse = (TH1D*) theFile0->Get("jpsi0 Mass: Coarse");
  TH1D *h_jpsi_mass_fine = (TH1D*) theFile0->Get("jpsi0 Mass: Fine");
  TH1D *h_jpsi_rapidity = (TH1D*) theFile0->Get("jpsi0 Rapidity");
  TH1D *h_jpsi_pT = (TH1D*) theFile0->Get("jpsi0 p_{T}");
  TH1D *h_jpsi_z_vertex_differece_x = (TH1D*) theFile0->Get("jpsi_vtx_x - z_vtx_x");
  TH1D *h_jpsi_z_vertex_differece_y = (TH1D*) theFile0->Get("jpsi_vtx_y - z_vtx_y");
  TH1D *h_jpsi_z_vertex_differece_z = (TH1D*) theFile0->Get("jpsi_vtx_z - z_vtx_z");
  TH1D *h_jpsi_distance = (TH1D*) theFile0->Get("jpsi0 distance");
  TH1D *h_jpsi_dist_err = (TH1D*) theFile0->Get("jpsi0 dist_err");
  TH1D *h_jpsi_chi2 = (TH1D*) theFile0->Get("jpsi0 chi2");
  TH1D *h_jpsi_distance_xy = (TH1D*) theFile0->Get("jpsi0 distance_xy");
  TH1D *h_jpsi_dist_err_xy = (TH1D*) theFile0->Get("jpsi0 dist_err_xy");
  TH1D *h_jpsi_chi2_xy = (TH1D*) theFile0->Get("jpsi0 chi2_xy");
  TH1D *h_jpsi_z_difference_pT = (TH1D*) theFile0->Get("zPt - jpsipT");
  TH1D *h_jpsi_tau_xy = (TH1D*) theFile0->Get("jpsi0 tau_xy");
  TH1D *h_jpsi_tau_xy_fine = (TH1D*) theFile0->Get("jpsi0 tau_xy_fine");
  TH1D *h_jpsi_tau_z = (TH1D*) theFile0->Get("jpsi0 tau_z");
  TH1D *h_jpsi_tau_z_fine = (TH1D*) theFile0->Get("jpsi0 tau_z_fine");
  TH2D *h_jpsi_mass_vs_chi2 = (TH2D*) theFile0->Get("jpsi0 mass vs chi2");
  TH2D *h_jpsi_tau_xy_vs_tau_z = (TH2D*) theFile0->Get("jpsi0 tau_xy vs tau_z");
  TH2D *h_jpsi_tau_xy_vs_z_vertex_difference = (TH2D*) theFile0->Get("jpsi0 tau_xy vs z vertex difference");
  TH2D *h_jpsi_tau_z_vs_z_vertex_difference = (TH2D*) theFile0->Get("jpsi0 tau_z vs z vertex difference");
  TH1D *h_dimuon_vertex_probability = (TH1D*) theFile0->Get("dimuon vertex probability");
  TH1D *h_dimuon_delta_phi = (TH1D*) theFile0->Get("dimuon delta phi");
  TH1D *h_dimuon_delta_eta = (TH1D*) theFile0->Get("dimuon delta eta");
  TH1D *h_dimuon_deltaR = (TH1D*) theFile0->Get("dimuon deltaR");
  TH1D *h_mu0_pT = (TH1D*) theFile0->Get("p_{T,mu_{0}}");
  TH1D *h_mu1_pT = (TH1D*) theFile0->Get("p_{T,mu_{1}}");
  TH1D *h_mu0_eta = (TH1D*) theFile0->Get("#eta_{mu_{0}}");
  TH1D *h_mu1_eta = (TH1D*) theFile0->Get("#eta_{mu_{1}}");
  TH1D *h_mu0_phi = (TH1D*) theFile0->Get("#phi_{mu_{0}}");
  TH1D *h_mu1_phi = (TH1D*) theFile0->Get("#phi_{mu_{1}}");
  TH1D *h_mu0_charge = (TH1D*) theFile0->Get("charge_{mu_{0}}");
  TH1D *h_mu1_charge = (TH1D*) theFile0->Get("charge_{mu_{1}}");
  TH1D *h_jet_pT = (TH1D*) theFile0->Get("p_{T,jet}");
  TH1D *h_jet_eta = (TH1D*) theFile0->Get("#eta_{jet}");
  TH1D *h_jet_btag_discriminator = (TH1D*) theFile0->Get("jet btag discriminator");
  TH1D *h_muon_jet_pT = (TH1D*) theFile0->Get("p_{T,muon_jet}");
  TH1D *h_z_muon_jet_difference_pT = (TH1D*) theFile0->Get("z p_{T} - muon jet p_{T}");
  TH1D *h_dimuon_muon_jet_difference_pT = (TH1D*) theFile0->Get("dimuon p_{T} - muon jet p_{T}");
  TH2D *h_muon_jet_pT_vs_z_pT = (TH2D*) theFile0->Get("muon jet p_{T} vs z p_{T}");
  TH2D *h_muon_jet_pT_vs_dimuon_pT = (TH2D*) theFile0->Get("muon jet p_{T} vs dimuon p_{T}");
  TH2D *h_muon_jet_phi_vs_z_phi = (TH2D*) theFile0->Get("muon jet phi vs z phi");
  TH2D *h_muon_jet_phi_vs_dimuon_phi = (TH2D*) theFile0->Get("muon jet phi vs dimuon phi");
  TH1D *h_muon_jet_eta = (TH1D*) theFile0->Get("#eta_{muon_jet}");
  TH1D *h_muon_jet_btag_discriminator = (TH1D*) theFile0->Get("muon jet btag discriminator");
  TH1D *h_vertex_position_x = (TH1D*) theFile0->Get("vertex position x");
  TH1D *h_vertex_position_y = (TH1D*) theFile0->Get("vertex position y");
  TH1D *h_vertex_position_z = (TH1D*) theFile0->Get("vertex position z");
  TH1D *h_primary_vertex_position_x = (TH1D*) theFile0->Get("primary vertex position x");
  TH1D *h_primary_vertex_position_y = (TH1D*) theFile0->Get("primary vertex position y");
  TH1D *h_primary_vertex_position_z = (TH1D*) theFile0->Get("primary vertex position z");
  TH1D *h_z_vertex_position_x = (TH1D*) theFile0->Get("z vertex position x");
  TH1D *h_z_vertex_position_y = (TH1D*) theFile0->Get("z vertex position y");
  TH1D *h_z_vertex_position_z = (TH1D*) theFile0->Get("z vertex position z");
  TH2D *h_primary_vertex_x_vs_z_vertex_x_position = (TH2D*) theFile0->Get("primary vertex x pos vs Z vertex x pos");
  TH2D *h_primary_vertex_y_vs_z_vertex_y_position = (TH2D*) theFile0->Get("primary vertex y pos vs Z vertex y pos");
  TH2D *h_primary_vertex_z_vs_z_vertex_z_position = (TH2D*) theFile0->Get("primary vertex z pos vs Z vertex z pos");
  TH1D *h_dimuon_vertex_position_x = (TH1D*) theFile0->Get("dimuon vertex position x");
  TH1D *h_dimuon_vertex_position_y = (TH1D*) theFile0->Get("dimuon vertex position y");
  TH1D *h_dimuon_vertex_position_z = (TH1D*) theFile0->Get("dimuon vertex position z");
  TH1D *h_z_jpsi_delta_phi = (TH1D*) theFile0->Get("z jpsi delta phi");
  TH1D *h_Nmu = (TH1D*) theFile0->Get("N_{mu}");
  TH1D *h_Njet = (TH1D*) theFile0->Get("N_{jet}");
  TH1D *h_Nmuonjet = (TH1D*) theFile0->Get("N_{muon_jet}");
  TH1D *h_Njpsi = (TH1D*) theFile0->Get("N_{jpsi}");
  TH1D *h_Ne = (TH1D*) theFile0->Get("N_{e}");
  TH1D *h_Nvert = (TH1D*) theFile0->Get("N_{Vertices}");


  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();
  h_z0_mass_all->Draw();
  TImage *img = TImage::Create();
  img->FromPad(c1);
  img->WriteImage("TestPlots/z0_mass_all.png");
  delete img;

  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->cd();
  h_z0_mass_all->Draw();
  TImage *img = TImage::Create();
  img->FromPad(c2);
  img->WriteImage("TestPlots/z0_mass_all.png");
  delete img;

}
