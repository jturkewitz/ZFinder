#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <math.h>
#include <stdio.h>
//void compare_polarizations (string file_name_mc, string file_name_data_zmumu, string file_name_data_zee)
void compare_efficiencies (string file_name_mc, string file_name_data_inclusive_z, string file_name_data_associated_z,  bool use_z_to_ee, string out_dir)
{
   
  for (int i=0;i<6;++i) {
    run_compare_efficiencies(file_name_mc, file_name_data_inclusive_z, file_name_data_associated_z,  use_z_to_ee, i, out_dir);
  }
}
void run_compare_efficiencies(string file_name_mc, string file_name_data_inclusive_z, string file_name_data_associated_z,  bool use_z_to_ee, int pt_slice, string out_dir)
{
  std::cout << "pt_slice " << pt_slice << std::endl;
  
  //TODO give this better documentation/ a more descriptive name
  //TFile *file_mc = new TFile("/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/test9d10.root");
  TFile *file_mc = new TFile( file_name_mc.c_str() );
  TFile *file_data_inclusive_z = new TFile( file_name_data_inclusive_z.c_str() );
  TFile *file_data_associated_z = new TFile( file_name_data_associated_z.c_str() );

  //TODO put this in rootrc
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  TH1D *h_z_pt_mc_truth = (TH1D*) file_mc->Get("ZFinder/MC_All/z p_{T}");

  std::string hist_path_reco;
  std::string hist_path_reco_mc;
  std::string hist_path_truth;
  std::string hist_path_associated;
  if(use_z_to_ee) {
    hist_path_reco = "ZFinder/Z_To_Electrons/z p_{T}";
    hist_path_reco_mc = "ZFinder/Z_To_Electrons/z_pt_mc";
    hist_path_truth = "ZFinder/MC_All/z p_{T}";
    
    if (pt_slice == 0) {
      hist_path_associated = "ZFinder/Z_To_Electrons_And_Jpsi/z p_{T}";
    }
    else if (pt_slice == 1) {
      hist_path_associated = "ZFinder/Z_To_Electrons_And_Jpsi/z_pt_jpsi_8p5to10";
    }
    else if (pt_slice == 2) {
      hist_path_associated = "ZFinder/Z_To_Electrons_And_Jpsi/z_pt_jpsi_10to14";
    }
    else if (pt_slice == 3) {
      hist_path_associated = "ZFinder/Z_To_Electrons_And_Jpsi/z_pt_jpsi_14to18";
    }
    else if (pt_slice == 4) {
      hist_path_associated = "ZFinder/Z_To_Electrons_And_Jpsi/z_pt_jpsi_18to30";
    }
    else if (pt_slice == 5) {
      hist_path_associated = "ZFinder/Z_To_Electrons_And_Jpsi/z_pt_jpsi_30to100";
    }
    else {
      std::cout << "pt_slice out of range" << std::endl;
    }
  }
  else {
    hist_path_reco = "ZFinder/Z_To_Muons/Z From Muons p_{T}";
    hist_path_reco_mc = "ZFinder/Z_To_Muons/Z From Muons p_{T} MC";
    hist_path_truth = "ZFinder/MC_All/Z From Muons p_{T}";
    if (pt_slice == 0) {
      hist_path_associated = "ZFinder/Z_To_Muons_And_Jpsi/Z From Muons p_{T}";
    }
    else if (pt_slice == 1) {
      hist_path_associated = "ZFinder/Z_To_Muons_And_Jpsi/z_from_muons_pt_jpsi_8p5to10";
    }
    else if (pt_slice == 2) {
      hist_path_associated = "ZFinder/Z_To_Muons_And_Jpsi/z_from_muons_pt_jpsi_10to14";
    }
    else if (pt_slice == 3) {
      hist_path_associated = "ZFinder/Z_To_Muons_And_Jpsi/z_from_muons_pt_jpsi_14to18";
    }
    else if (pt_slice == 4) {
      hist_path_associated = "ZFinder/Z_To_Muons_And_Jpsi/z_from_muons_pt_jpsi_18to30";
    }
    else if (pt_slice == 5) {
      hist_path_associated = "ZFinder/Z_To_Muons_And_Jpsi/z_from_muons_pt_jpsi_30to100";
    }
    else {
      std::cout << "pt_slice out of range" << std::endl;
    }
  }

  TH1D *h_z_pt_mc_truth = (TH1D*) file_mc->Get(hist_path_truth.c_str());
  TH1D *h_z_pt_mc_reco = (TH1D*) file_mc->Get(hist_path_reco_mc.c_str());
  TH1D *h_z_pt_data_inclusive = (TH1D*) file_data_inclusive_z->Get(hist_path_reco.c_str());
  TH1D *h_z_pt_data_associated = (TH1D*) file_data_associated_z->Get(hist_path_associated.c_str());

  TH1D *h_z_pt_data_associated_clone = h_z_pt_data_associated->Clone();

  h_z_pt_mc_reco->GetXaxis()->SetRangeUser(0,200);


  double norm_data = h_z_pt_data_inclusive->Integral();
  std::cout << "norm_data " << norm_data << std::endl;
  double norm_data_associated = h_z_pt_data_associated->Integral();
  std::cout << "norm_data_associated " << norm_data_associated << std::endl;


  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->cd();


  //c2->SetLogy();


  //h_jpsi_pt_zjpsi_mpi->Rebin(5);
  //h_jpsi_pt->Rebin(5);
  //h_jpsi_pt_zmumu_jpsi->Rebin(5);

  //h_z_pt_mc_truth->GetXaxis()->SetTitle("J/#psi cos #mu+");
  h_z_pt_mc_reco->GetYaxis()->SetTitle("Efficiency");

  h_z_pt_mc_reco->Divide(h_z_pt_mc_truth);
  if (use_z_to_ee) {
    h_z_pt_mc_reco->SetTitle("Z->ee Madgraph MC");
  }
  else {
    h_z_pt_mc_reco->SetTitle("Z->#mu#mu Madgraph MC");
  }
  h_z_pt_mc_reco->Draw();

  TCanvas *c3 = new TCanvas("c3", "c3");
  c3->cd();

  //h_z_pt_data_inclusive->Scale(1.0/norm_data);
  h_z_pt_data_inclusive->Divide(h_z_pt_mc_reco);
  h_z_pt_data_inclusive->GetYaxis()->SetTitle("Yield");
  if (use_z_to_ee) {
    h_z_pt_data_inclusive->SetTitle("Z->ee Inclusive");
  }
  else {
    h_z_pt_data_inclusive->SetTitle("Z->#mu#mu Inclusive");
  }
  h_z_pt_data_inclusive->Draw();

  TCanvas *c4 = new TCanvas("c4", "c4");
  c4->cd();

  //h_z_pt_data_associated->Scale(1.0/norm_data_associated);
  h_z_pt_data_associated->Divide(h_z_pt_mc_reco);
  h_z_pt_data_associated->GetYaxis()->SetTitle("Yield");
  if (use_z_to_ee) {
    h_z_pt_data_associated->SetTitle("Z->ee + Associated J/#psi");
  }
  else {
    h_z_pt_data_associated->SetTitle("Z->#mu#mu + Associated J/#psi");
  }
  h_z_pt_data_associated->Draw();

  double final_integral = h_z_pt_data_inclusive->Integral();
  double entries = h_z_pt_data_inclusive->GetEntries();
  double efficiency_inclusive = norm_data / final_integral;
  std::cout << "integral inclusive" << final_integral << std::endl;
  //std::cout << "entries " << entries << std::endl;
  std::cout << "efficiency inclusive " << efficiency_inclusive << std::endl;

  double final_integral_associated = h_z_pt_data_associated->Integral();
  double entries_associated = h_z_pt_data_associated->GetEntries();
  double efficiency_associated = norm_data_associated / final_integral_associated;
  std::cout << "integral associated " << final_integral_associated << std::endl;
  //std::cout << "entries_associated " << entries_associated << std::endl;
  std::cout << "efficiency associated " << efficiency_associated << std::endl;

  std::cout << "systematic " << efficiency_associated / efficiency_inclusive << std::endl;
  //std::cout << "eff " << eff << std::endl;

  //Double_t xl1=.55, yl1=0.70, xl2=xl1+.30, yl2=yl1+.20;
  //TLegend *leg2 = new TLegend(xl1,yl1,xl2,yl2);
  //leg2->AddEntry(h_jpsi_cos_mu_plus_pol_0,"MC unpolarized","lep");
  //leg2->AddEntry(h_jpsi_cos_mu_plus_pol_neg,"MC pol neg 1","lep");
  //leg2->AddEntry(h_jpsi_cos_mu_plus_pol_pos,"MC pol pos 1","lep");
  //leg2->AddEntry(h_jpsi_cos_mu_plus_inc_jpsi,"Inclusive Jpsi Data","lep");
  ////leg2->AddEntry(h_jpsi_pt_zjpsi_mpi,"J/#psi MC (Z+J/#psi MPI)","lep");
  //leg2->Draw();

  //int nbinsx = h_z_pt_data_associated->GetXaxis()->GetNbins();
  int nxbins = h_z_pt_data_associated->GetNbinsX();
  //int nybins = calibMapEB_->GetNbinsY();
  double error_sum_sq = 0;
  for (int x=1;x<=nxbins;++x)
  {
    double binentsM = h_z_pt_data_associated_clone->GetBinContent(x);
    if (binentsM == 0) {
      continue;
    }
    double binents_eff = h_z_pt_mc_reco->GetBinContent(x);
    double error = pow( binentsM, 0.5 );


    error_sum_sq = error_sum_sq + pow((final_integral_associated  - 1.0 / binents_eff * norm_data_associated ) / pow(final_integral_associated,2.0) * error, 2.0) ;
    

    //std::cout << binentsM << std::endl; 
  }
  
  std::cout << "error sum sq " << error_sum_sq << std::endl;
  std::cout << "error " << pow(error_sum_sq,0.5) << std::endl;
  std::cout << "error ratio " << pow(error_sum_sq,0.5) / efficiency_inclusive << std::endl;

  

  std::string image_name2 = out_dir;
  if (use_z_to_ee) {
    image_name2 = image_name2.append("zee_efficiency_mc");
  }
  else {
    image_name2 = image_name2.append("zmumu_efficiency_mc");
  }
  image_name2.append(".png");
  c2->Print(image_name2.c_str() , "png");

  std::string image_name3 = out_dir;
  if (use_z_to_ee) {
    image_name3 = image_name3.append("zee_efficiency_data");
  }
  else {
    image_name3 = image_name3.append("zmumu_efficiency_data");
  }
  image_name3.append(".png");
  c3->Print(image_name3.c_str() , "png");

  std::string image_name4 = out_dir;
  if (use_z_to_ee) {
    image_name4 = image_name4.append("zee_associated_efficiency_data");
  }
  else {
    image_name4 = image_name4.append("zmumu_associated_efficiency_data");
  }
  image_name4.append(".png");
  c4->Print(image_name4.c_str() , "png");
}
