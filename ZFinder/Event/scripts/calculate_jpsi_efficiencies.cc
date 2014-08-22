#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
void calculate_jpsi_efficiencies (string file_name )
{
  
  //TODO maybe make the file name an input argument?
  //TFile *theFile0 = new TFile("/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/test9d10.root");
  TFile *theFile0 = new TFile( file_name.c_str() );

  //TODO put this in rootrc
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  TH1D *h_n_truth_matched_muons_all = (TH1D*) theFile0->Get("ZFinder/All/N_truth_matched_jpsi_muons");
  TH1D *h_n_truth_matched_muons_dimuon = (TH1D*) theFile0->Get("ZFinder/Dimuon/N_truth_matched_jpsi_muons");
  TH1D *h_n_truth_matched_muons_dimuon_soft = (TH1D*) theFile0->Get("ZFinder/Dimuon_Soft/N_truth_matched_jpsi_muons");
  TH1D *h_n_truth_matched_muons_dimuon_vtx_comp = (TH1D*) theFile0->Get("ZFinder/Dimuon_Vertex_Compatible/N_truth_matched_jpsi_muons");
  TH1D *h_n_truth_matched_muons_dimuon_primary_vert = (TH1D*) theFile0->Get("ZFinder/Dimuon_Primary_Vertex/N_truth_matched_jpsi_muons");
  TH1D *h_n_truth_matched_muons_jpsi = (TH1D*) theFile0->Get("ZFinder/Jpsi/N_truth_matched_jpsi_muons");

  TH1D *h_jpsi_mass_fine_jpsi_mc = (TH1D*) theFile0->Get("ZFinder/MC_Jpsi/jpsi Mass: Fine");

  double mc_events = h_jpsi_mass_fine_jpsi_mc->Integral();

  double two_truth_matched_muons_dimuon = h_n_truth_matched_muons_dimuon->GetBinContent(3) ;
  double two_truth_matched_muons_dimuon_soft = h_n_truth_matched_muons_dimuon_soft->GetBinContent(3);
  double two_truth_matched_muons_dimuon_vtx_comp = h_n_truth_matched_muons_dimuon_vtx_comp->GetBinContent(3);
  double two_truth_matched_muons_dimuon_primary_vert  = h_n_truth_matched_muons_dimuon_primary_vert->GetBinContent(3);
  double two_truth_matched_muons_jpsi = h_n_truth_matched_muons_jpsi->GetBinContent(3);

  std::cout << "mc_events " << mc_events << std::endl;
  std::cout << "two_truth_matched_muons_dimuon " << two_truth_matched_muons_dimuon << std::endl;
  std::cout << "two_truth_matched_muons_dimuon_soft " << two_truth_matched_muons_dimuon_soft  << std::endl;
  std::cout << "two_truth_matched_muons_dimuon_vtx_comp " << two_truth_matched_muons_dimuon_vtx_comp  << std::endl;
  std::cout << "two_truth_matched_muons_dimuon_primary_vert " << two_truth_matched_muons_dimuon_primary_vert  << std::endl;
  std::cout << "two_truth_matched_muons_jpsi " << two_truth_matched_muons_jpsi  << std::endl;

}
