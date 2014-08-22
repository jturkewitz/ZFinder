#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
void make_images (string file_name )
{
  
  //TODO maybe make the file name an input argument?
  //TFile *theFile0 = new TFile("/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/test9d10.root");
  TFile *theFile0 = new TFile( file_name.c_str() );

  //TODO put this in rootrc
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);

  TH1D *h_jpsi_mass_fine_all = (TH1D*) theFile0->Get("ZFinder/All/jpsi Mass: Fine");
  TH1D *h_jpsi_mass_fine_dimuon = (TH1D*) theFile0->Get("ZFinder/Dimuon/jpsi Mass: Fine");
  TH1D *h_jpsi_mass_fine_dimuon_soft = (TH1D*) theFile0->Get("ZFinder/Dimuon_Soft/jpsi Mass: Fine");
  TH1D *h_jpsi_mass_fine_dimuon_vtx_comp = (TH1D*) theFile0->Get("ZFinder/Dimuon_Vertex_Compatible/jpsi Mass: Fine");
  TH1D *h_jpsi_mass_fine_jpsi_primary_vert = (TH1D*) theFile0->Get("ZFinder/Dimuon_Primary_Vertex/jpsi Mass: Fine");
  TH1D *h_jpsi_mass_fine_jpsi = (TH1D*) theFile0->Get("ZFinder/Jpsi/jpsi Mass: Fine");

  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();
  c1->SetLogy();
  double h_jpsi_mass_fine_all_max = h_jpsi_mass_fine_all->GetMaximum();
  h_jpsi_mass_fine_all->SetLineColor(kBlue-2);
  h_jpsi_mass_fine_all->GetYaxis()->SetRangeUser(0.5,h_jpsi_mass_fine_all_max*8.0);
  h_jpsi_mass_fine_all->SetLineColor(kBlue-2);
  h_jpsi_mass_fine_all->Draw();
  h_jpsi_mass_fine_dimuon->SetLineColor(kRed-2);
  h_jpsi_mass_fine_dimuon->Draw("same");
  h_jpsi_mass_fine_dimuon_soft->SetLineColor(kCyan-2);
  h_jpsi_mass_fine_dimuon_soft->Draw("same");
  h_jpsi_mass_fine_dimuon_vtx_comp->SetLineColor(kGreen-2);
  h_jpsi_mass_fine_dimuon_vtx_comp->Draw("same");
  h_jpsi_mass_fine_dimuon_primary_vert->SetLineColor(kOrange-2);
  h_jpsi_mass_fine_dimuon_primary_vert->Draw("same");
  h_jpsi_mass_fine_jpsi->SetLineColor(kMagenta-2);
  h_jpsi_mass_fine_jpsi->Draw("same");

  //Double_t xl1=.59, yl1=0.65, xl2=xl1+.30, yl2=yl1+.225;
  Double_t xl1=.65, yl1=0.70, xl2=xl1+.30, yl2=yl1+.225;
  TLegend *leg1 = new TLegend(xl1,yl1,xl2,yl2);
  leg1->AddEntry(h_jpsi_mass_fine_all,"All","lep");
  leg1->AddEntry(h_jpsi_mass_fine_dimuon,"Dimuon","lep");
  leg1->AddEntry(h_jpsi_mass_fine_dimuon_soft,"Dimuon Soft","lep");
  leg1->AddEntry(h_jpsi_mass_fine_dimuon_vtx_comp,"Dimuon Vertex Comp","lep");
  leg1->AddEntry(h_jpsi_mass_fine_dimuon_primary_vert,"Dimuon Primary Vert","lep");
  leg1->AddEntry(h_jpsi_mass_fine_jpsi,"Jpsi","lep");
  leg1->Draw();

  std::string image_name = "";
  image_name.append("/home/user1/turkewitz/public_html/ZPhysics/tmp/TestPlots/");
  image_name.append(file_name);
  image_name.append(".png");
  c1->Print(image_name.c_str() , "png");
  //c1->Close();
  //delete c1;
}
