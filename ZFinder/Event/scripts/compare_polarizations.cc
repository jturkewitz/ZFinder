#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
//void compare_polarizations (string file_name_mc, string file_name_data_zmumu, string file_name_data_zee)
void compare_polarizations (string file_name_mc, string file_name_data_inclusive_jpsi, string out_dir)
{
  
  //TODO give this better documentation/ a more descriptive name
  //TFile *file_mc = new TFile("/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/test9d10.root");
  TFile *file_mc = new TFile( file_name_mc.c_str() );
  TFile *file_data_inclusive_jpsi = new TFile( file_name_data_inclusive_jpsi.c_str() );

  //TODO put this in rootrc
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  TH1D *h_jpsi_cos_mu_plus_pol_0 = (TH1D*) file_mc->Get("ZFinder/Jpsi/jpsi_cos_mu_plus");
  TH1D *h_jpsi_cos_mu_plus_pol_pos = (TH1D*) file_mc->Get("ZFinder/Jpsi/jpsi_cos_mu_plus_lambdaNeg1");
  TH1D *h_jpsi_cos_mu_plus_pol_neg = (TH1D*) file_mc->Get("ZFinder/Jpsi/jpsi_cos_mu_plus_lambda1");
  TH1D *h_jpsi_cos_mu_plus_inc_jpsi = (TH1D*) file_data_inclusive_jpsi->Get("ZFinder/Jpsi/jpsi_cos_mu_plus");

  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->cd();
  //c2->SetLogy();

  //h_jpsi_pt_zjpsi_mpi->Rebin(5);
  //h_jpsi_pt->Rebin(5);
  //h_jpsi_pt_zmumu_jpsi->Rebin(5);
  h_jpsi_cos_mu_plus_pol_0->GetXaxis()->SetTitle("J/#psi cos #mu+");
  h_jpsi_cos_mu_plus_pol_0->GetYaxis()->SetTitle("Norm. Events");
  h_jpsi_cos_mu_plus_pol_0->SetTitle("J/#psi cos #mu+");

  double norm_mc_jpsi_pol_0 = h_jpsi_cos_mu_plus_pol_0->Integral();
  double norm_mc_jpsi_pol_neg = h_jpsi_cos_mu_plus_pol_neg->Integral();
  double norm_mc_jpsi_pol_pos = h_jpsi_cos_mu_plus_pol_pos->Integral();
  double norm_data = h_jpsi_cos_mu_plus_inc_jpsi->Integral();
  //std::cout << norm_mc_jpsi_pol_neg << norm_mc_jpsi_pol_pos << std::endl;

  double norm_data = h_jpsi_cos_mu_plus_inc_jpsi->Integral();

  h_jpsi_cos_mu_plus_pol_0->Scale(1.0/norm_mc_jpsi_pol_0);
  h_jpsi_cos_mu_plus_pol_0->Draw();
  h_jpsi_cos_mu_plus_pol_0->GetYaxis()->SetRangeUser(0.0,0.02);


  h_jpsi_cos_mu_plus_pol_neg->Scale(1.0/norm_mc_jpsi_pol_neg);
  h_jpsi_cos_mu_plus_pol_neg->Draw("same");
  h_jpsi_cos_mu_plus_pol_neg->SetLineColor(kRed-2);

  h_jpsi_cos_mu_plus_pol_pos->Scale(1.0/norm_mc_jpsi_pol_pos);
  h_jpsi_cos_mu_plus_pol_pos->Draw("same");
  h_jpsi_cos_mu_plus_pol_pos->SetLineColor(kGreen-2);

  h_jpsi_cos_mu_plus_inc_jpsi->Scale(1.0/norm_data);
  h_jpsi_cos_mu_plus_inc_jpsi->Draw("same");
  h_jpsi_cos_mu_plus_inc_jpsi->SetLineColor(kBlack);


  //h_jpsi_pt_zmumu_jpsi->Scale(1/norm_data);

  //h_jpsi_pt_zmumu_jpsi->SetLineColor(kRed-2);
  //h_jpsi_pt_zmumu_jpsi->Draw("same");

  //h_jpsi_pt_zjpsi_mpi->Scale(1/norm_mc_zjpsi*norm_data);
  //h_jpsi_pt_zjpsi_mpi->SetLineColor(kGreen-2);
  //h_jpsi_pt_zjpsi_mpi->Draw("same");

  Double_t xl1=.55, yl1=0.70, xl2=xl1+.30, yl2=yl1+.20;
  TLegend *leg2 = new TLegend(xl1,yl1,xl2,yl2);
  leg2->AddEntry(h_jpsi_cos_mu_plus_pol_0,"MC unpolarized","lep");
  leg2->AddEntry(h_jpsi_cos_mu_plus_pol_neg,"MC pol neg 1","lep");
  leg2->AddEntry(h_jpsi_cos_mu_plus_pol_pos,"MC pol pos 1","lep");
  leg2->AddEntry(h_jpsi_cos_mu_plus_inc_jpsi,"Inclusive Jpsi Data","lep");
  //leg2->AddEntry(h_jpsi_pt_zjpsi_mpi,"J/#psi MC (Z+J/#psi MPI)","lep");
  leg2->Draw();

  std::string image_name2 = out_dir;
  image_name2 = image_name2.append("jpsi_cos_mu_plus");
  image_name2.append(".png");
  c2->Print(image_name2.c_str() , "png");
}
