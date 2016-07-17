#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <RooBinning.h>
#include <RooWorkspace.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <TLatex.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooExtendPdf.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <TCut.h>
#include <TTree.h>
#include <TH1.h>
#include <THStack.h>
#include <TFile.h>
#include <TF1.h>
#include <RooCBShape.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooMCStudy.h>
#include <RooStats/SPlot.h>
#include <RooNumIntConfig.h>
#include <string>
#include "tdrStyle.C"

using namespace RooFit ;
using namespace RooStats ;

//const float PDG_JPSI_MASS = 3.096916;

void make_theory_comparison(std::string image_dir = "~/public_html/ZPhysics/tmp/Test90/", bool use_only_z_to_muons = false, bool use_only_z_to_electrons = false) { 

  //std::cout << "testing" << std::endl;
  //gStyle->SetOptStat(0);
  //gStyle->SetOptFit(0);

  setTDRStyle();
  gStyle->SetLineWidth(2.);
  gStyle->SetHistLineWidth(2.5);
  gROOT->ForceStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.05);
  
  // ************************
  // ************************
  // Start the x-sec plotting
  // ************************
  // ************************

  //RooBinning bins(6, 100) ;
  //bins.addUniform(1, 6, 8.5);
  //bins.addUniform(1, 8.5, 10);
  //bins.addUniform(2, 10, 18);
  //bins.addUniform(1, 18, 30);
  //bins.addUniform(1, 30, 100);

  double inclusive_z_to_muons_events;
  double inclusive_z_to_electrons_events;
  double inclusive_z_to_leptons_events;
  //inclusive_z_to_muons_events     = 7.773024E06;
  //inclusive_z_to_electrons_events = 5.086549E06;
  inclusive_z_to_muons_events     = 6.93504E06;
  inclusive_z_to_electrons_events = 4.63707E06;
  inclusive_z_to_leptons_events = inclusive_z_to_muons_events + inclusive_z_to_electrons_events;

  TLatex mark;
  mark.SetTextSize(0.035);
  mark.SetNDC(true);

  TCanvas* c_jpsi_diff_xsec = new TCanvas("c_jpsi_diff_xsec","J/#psi p_{T}", 1400, 700); c_jpsi_diff_xsec->Divide(2,1);
  // plot the prompt part
  c_jpsi_diff_xsec->cd(1);

  //RooRealVar *onia_pt     = new RooRealVar("onia_pt",  "onia_pt", 0, 100);
  //onia_pt->setBinning(bins);
  //TH1* h_zjpsi_atlas_prompt = zjpsi_prompt_hist.createHistogram("h_zjpsi_atlas_prompt", *onia_pt);

  double bins[] = {8.5,10,14,18,30,100};
  double bins_offset[] = {8.7,10.2,14.2,18.2,30.2,101.5};
  int bin_number = 5;
  TH1F* h_zjpsi_atlas_prompt = new TH1F ("h_zjpsi_atlas_prompt", "h_zjpsi_atlas_prompt", 5, bins_offset);
  h_zjpsi_atlas_prompt->SetBinContent(1, 10.8E-7); h_zjpsi_atlas_prompt->SetBinError(1, TMath::Sqrt(5.6*5.6    +1.9*1.9    )*1E-7);
  h_zjpsi_atlas_prompt->SetBinContent(2,  5.6E-7); h_zjpsi_atlas_prompt->SetBinError(2, TMath::Sqrt(1.9*1.9    +0.8*0.8    )*1E-7);
  h_zjpsi_atlas_prompt->SetBinContent(3,  1.9E-7); h_zjpsi_atlas_prompt->SetBinError(3, TMath::Sqrt(1.1*1.1    +0.1*0.1    )*1E-7);
  h_zjpsi_atlas_prompt->SetBinContent(4, 0.87E-7); h_zjpsi_atlas_prompt->SetBinError(4, TMath::Sqrt(0.37*0.37  +0.12*0.12  )*1E-7);
  h_zjpsi_atlas_prompt->SetBinContent(5, 0.09E-7); h_zjpsi_atlas_prompt->SetBinError(5, TMath::Sqrt(0.037*0.037+0.012*0.012)*1E-7);
  h_zjpsi_atlas_prompt->SetMarkerStyle(24);

  //TH1* h_zjpsi_cms_prompt = zjpsi_prompt_hist.createHistogram("h_zjpsi_cms_prompt", *onia_pt);

  TH1F* h_zjpsi_cms_prompt = new TH1F ("h_zjpsi_cms_prompt", "h_zjpsi_cms_prompt", 5, bins);
  h_zjpsi_cms_prompt->SetBinContent(1, 10.5574277113E-7); h_zjpsi_cms_prompt->SetBinError(1, TMath::Sqrt(6.1043030461**2 + 2.7862418006**2)*1E-7);
  h_zjpsi_cms_prompt->SetBinContent(2, 5.0144490001E-7); h_zjpsi_cms_prompt->SetBinError(2, TMath::Sqrt(1.9208565661**2 + 1.0818902095**2)*1E-7);
  h_zjpsi_cms_prompt->SetBinContent(3, 0.2115452932E-7); h_zjpsi_cms_prompt->SetBinError(3, TMath::Sqrt(0.9112500667**2 + 0.0269464576**2)*1E-7);
  h_zjpsi_cms_prompt->SetBinContent(4, 0.2321810142E-7); h_zjpsi_cms_prompt->SetBinError(4, TMath::Sqrt(0.2788649981**2 + 0.028014271**2)*1E-7);
  h_zjpsi_cms_prompt->SetBinContent(5, 0.0358347174E-7); h_zjpsi_cms_prompt->SetBinError(5, TMath::Sqrt(0.0359961998**2 + 0.005662698**2)*1E-7);
  h_zjpsi_cms_prompt->SetMarkerStyle(25);

  TH1F* h_zjpsi_theory_prompt = new TH1F ("h_zjpsi_theory_prompt", "h_zjpsi_theory_prompt", 5, bins);
  h_zjpsi_theory_prompt->SetBinContent(1, 1.1471544715E-7); h_zjpsi_theory_prompt->SetBinError(1, 0.1147154472E-7);
  h_zjpsi_theory_prompt->SetBinContent(2, 0.7057926829E-7); h_zjpsi_theory_prompt->SetBinError(2, 0.0705792683E-7);
  h_zjpsi_theory_prompt->SetBinContent(3, 0.3445121951E-7); h_zjpsi_theory_prompt->SetBinError(3, 0.0344512195E-7);
  h_zjpsi_theory_prompt->SetBinContent(4, 0.144004065E-7); h_zjpsi_theory_prompt->SetBinError(4, 0.0144004065E-7);
  h_zjpsi_theory_prompt->SetBinContent(5, 0.0126829268E-7); h_zjpsi_theory_prompt->SetBinError(5, 0.0012682927E-7);
  h_zjpsi_theory_prompt->SetMarkerStyle(27);

  TH1F* h_zjpsi_dps_prompt = new TH1F ("h_zjpsi_dps_prompt", "h_zjpsi_dps_prompt", 5, bins);
  h_zjpsi_dps_prompt->SetBinContent(1, 3.8167422744E-7); h_zjpsi_dps_prompt->SetBinError(1, 1.2258395366*1E-7);
  h_zjpsi_dps_prompt->SetBinContent(2, 1.0692782379E-7); h_zjpsi_dps_prompt->SetBinError(2, 0.3434246919*1E-7);
  h_zjpsi_dps_prompt->SetBinContent(3, 0.2095427881E-7); h_zjpsi_dps_prompt->SetBinError(3, 0.0672997588*1E-7);
  h_zjpsi_dps_prompt->SetBinContent(4, 0.0295145771E-7); h_zjpsi_dps_prompt->SetBinError(4, 0.0094793237*1E-7);
  h_zjpsi_dps_prompt->SetBinContent(5, 0.0005232614E-7); h_zjpsi_dps_prompt->SetBinError(5, 0.0001680581*1E-7);
  h_zjpsi_dps_prompt->SetMarkerStyle(26);

  h_zjpsi_cms_prompt->Draw("E1");
  h_zjpsi_atlas_prompt->SetLineColor(kRed - 2);
  h_zjpsi_atlas_prompt->Draw("E1same");
  h_zjpsi_dps_prompt->SetLineColor(kGreen - 2);
  h_zjpsi_dps_prompt->Draw("E1same");
  h_zjpsi_theory_prompt->SetLineColor(kMagenta - 2);
  h_zjpsi_theory_prompt->Draw("E1same");

  Double_t xl1_l5=.60, yl1_l5=0.70, xl2_l5=xl1_l5+.3, yl2_l5=yl1_l5+.2;
  TLegend *leg5 = new TLegend(xl1_l5,yl1_l5,xl2_l5,yl2_l5);
  leg5->SetFillColor(kWhite);
  leg5->AddEntry(h_zjpsi_cms_prompt,"CMS prompt","l");
  leg5->AddEntry(h_zjpsi_atlas_prompt,"ATLAS prompt","l");
  leg5->AddEntry(h_zjpsi_dps_prompt,"DPS prompt","l");
  leg5->AddEntry(h_zjpsi_theory_prompt,"Theory prompt","l");
  leg5->SetShadowColor(0);
  leg5->Draw();
  h_zjpsi_cms_prompt->GetYaxis()->SetRangeUser(1E-11, 1E-5);
  h_zjpsi_cms_prompt->GetXaxis()->SetRangeUser(8.5, 100);
  //splot_tex_pt.DrawLatex(15, 3E-6, "#bf{#it{CMS}} Preliminary, #sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
  //splot_tex_pt.DrawLatex(15, 1E-6 , "#it{pp} #rightarrow prompt #it{J/#psi+Z} : #it{pp} #rightarrow #it{Z}");

  c_jpsi_diff_xsec->cd(1)->SetLogy(1);
  c_jpsi_diff_xsec->cd(1)->SetLogx(1);

  // plot the non-prompt part
  c_jpsi_diff_xsec->cd(2);
  //TH1* h_zjpsi_atlas_nonprompt = zjpsi_nonprompt_hist.createHistogram("h_zjpsi_atlas_nonprompt", *onia_pt);
  TH1F* h_zjpsi_atlas_nonprompt = new TH1F ("h_zjpsi_atlas_nonprompt", "h_zjpsi_atlas_nonprompt", 5, bins_offset);
  h_zjpsi_atlas_nonprompt->SetBinContent(1,   5.1E-7); h_zjpsi_atlas_nonprompt->SetBinError(1, TMath::Sqrt(4.2*4.2    +0.9*0.9    )*1E-7);
  h_zjpsi_atlas_nonprompt->SetBinContent(2,   9.2E-7); h_zjpsi_atlas_nonprompt->SetBinError(2, TMath::Sqrt(2.5*2.5    +1.2*1.2    )*1E-7);
  h_zjpsi_atlas_nonprompt->SetBinContent(3,   3.3E-7); h_zjpsi_atlas_nonprompt->SetBinError(3, TMath::Sqrt(1.2*1.2    +0.4*0.4    )*1E-7);
  h_zjpsi_atlas_nonprompt->SetBinContent(4,  3.04E-7); h_zjpsi_atlas_nonprompt->SetBinError(4, TMath::Sqrt(0.59*0.59  +0.04*0.04  )*1E-7);
  h_zjpsi_atlas_nonprompt->SetBinContent(5, 0.115E-7); h_zjpsi_atlas_nonprompt->SetBinError(5, TMath::Sqrt(0.039*0.039+0.002*0.002)*1E-7);
  h_zjpsi_atlas_nonprompt->SetMarkerStyle(24);

  //TH1* h_zjpsi_cms_nonprompt = zjpsi_nonprompt_hist.createHistogram("h_zjpsi_cms_nonprompt", *onia_pt);

  TH1F* h_zjpsi_cms_nonprompt = new TH1F ("h_zjpsi_cms_nonprompt", "h_zjpsi_cms_nonprompt", 5, bins);
  h_zjpsi_cms_nonprompt->SetBinContent(1, 13.400610612E-7); h_zjpsi_cms_nonprompt->SetBinError(1,TMath::Sqrt( 6.6905206225**2 + 3.5365945628**2)*1E-7);
  h_zjpsi_cms_nonprompt->SetBinContent(2, 5.7336688784E-7); h_zjpsi_cms_nonprompt->SetBinError(2,TMath::Sqrt( 2.1589059133**2 + 1.2370651739**2)*1E-7);
  h_zjpsi_cms_nonprompt->SetBinContent(3, 4.031344658E-7); h_zjpsi_cms_nonprompt->SetBinError(3, TMath::Sqrt( 1.5319663997**2 + 0.5135092165**2)*1E-7);
  h_zjpsi_cms_nonprompt->SetBinContent(4, 1.7323716183E-7); h_zjpsi_cms_nonprompt->SetBinError(4, TMath::Sqrt(0.4914166095**2 + 0.20902281**2)*1E-7);
  h_zjpsi_cms_nonprompt->SetBinContent(5, 0.0157856671E-7); h_zjpsi_cms_nonprompt->SetBinError(5, TMath::Sqrt(0.0337188967**2 + 0.0024944934**2)*1E-7);
  h_zjpsi_cms_nonprompt->SetMarkerStyle(25);

  TH1F* h_zjpsi_dps_nonprompt = new TH1F ("h_zjpsi_dps_nonprompt", "h_zjpsi_dps_nonprompt", 5, bins);
  h_zjpsi_dps_nonprompt->SetBinContent(1, 1.6465232421E-7); h_zjpsi_dps_nonprompt->SetBinError(1, 0.5288209533*1E-7);
  h_zjpsi_dps_nonprompt->SetBinContent(2, 0.5949677939E-7); h_zjpsi_dps_nonprompt->SetBinError(2, 0.191088366*1E-7);
  h_zjpsi_dps_nonprompt->SetBinContent(3, 0.2533485162E-7); h_zjpsi_dps_nonprompt->SetBinError(3, 0.0813690329*1E-7);
  h_zjpsi_dps_nonprompt->SetBinContent(4, 0.0346806073E-7); h_zjpsi_dps_nonprompt->SetBinError(4, 0.0111385199*1E-7);
  h_zjpsi_dps_nonprompt->SetBinContent(5, 0.0009995655E-7); h_zjpsi_dps_nonprompt->SetBinError(5, 0.0003210348*1E-7);
  h_zjpsi_dps_nonprompt->SetMarkerStyle(26);

  //TODO fix this - from DY MC
  TH1F* h_zjpsi_theory_nonprompt = new TH1F ("h_zjpsi_theory_nonprompt", "h_zjpsi_theory_nonprompt", 5, bins);
  h_zjpsi_theory_nonprompt->SetBinContent(1, 6.8176855787E-7); h_zjpsi_theory_nonprompt->SetBinError(1, 6.901411542E-7);
  h_zjpsi_theory_nonprompt->SetBinContent(2, 15.6250588698E-7); h_zjpsi_theory_nonprompt->SetBinError(2, 4.5750258489E-7);
  h_zjpsi_theory_nonprompt->SetBinContent(3, 3.0780504792E-7); h_zjpsi_theory_nonprompt->SetBinError(3, 1.7896424644E-7);
  h_zjpsi_theory_nonprompt->SetBinContent(4, 1.272559886E-7); h_zjpsi_theory_nonprompt->SetBinError(4, 0.611498553E-7);
  h_zjpsi_theory_nonprompt->SetBinContent(5, 0.2078958385E-7); h_zjpsi_theory_nonprompt->SetBinError(5, 0.0848366138E-7);
  h_zjpsi_theory_nonprompt->SetMarkerStyle(27);

  h_zjpsi_cms_nonprompt->SetLineColor(kBlue - 2);
  h_zjpsi_cms_nonprompt->Draw("E1");

  h_zjpsi_atlas_nonprompt->SetLineColor(kRed - 2);
  h_zjpsi_dps_nonprompt->SetLineColor(kGreen - 2);
  h_zjpsi_atlas_nonprompt->Draw("E1same");
  h_zjpsi_dps_nonprompt->Draw("E1same");
  h_zjpsi_theory_nonprompt->SetLineColor(kMagenta - 2);
  h_zjpsi_theory_nonprompt->Draw("E1same");

  Double_t xl1_l6=.60, yl1_l6=0.70, xl2_l6=xl1_l6+.3, yl2_l6=yl1_l6+.2;
  TLegend *leg6 = new TLegend(xl1_l6,yl1_l6,xl2_l6,yl2_l6);
  leg6->SetFillColor(kWhite);
  leg6->AddEntry(h_zjpsi_cms_nonprompt,"CMS nonprompt","l");
  leg6->AddEntry(h_zjpsi_atlas_nonprompt,"ATLAS nonprompt","l");
  leg6->AddEntry(h_zjpsi_dps_nonprompt,"DPS nonprompt","l");
  leg6->AddEntry(h_zjpsi_theory_nonprompt,"DY MC nonprompt","l");
  leg6->SetShadowColor(0);
  leg6->Draw();

  h_zjpsi_cms_nonprompt->GetYaxis()->SetRangeUser(1E-11, 1E-5);
  h_zjpsi_cms_nonprompt->GetXaxis()->SetRangeUser(8.5, 100);
  //splot_tex_pt.DrawLatex(15, 3E-6, "#bf{#it{CMS}} Preliminary, #sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
  //splot_tex_pt.DrawLatex(15, 1E-6 , "#it{pp} #rightarrow non-prompt #it{J/#psi+Z} : #it{pp} #rightarrow #it{Z}");
  c_jpsi_diff_xsec->cd(2)->SetLogy(1);
  c_jpsi_diff_xsec->cd(2)->SetLogx(1);
  //c_jpsi_diff_xsec->SaveAs("c_jpsi_diff_xsec.C");
  //c_jpsi_diff_xsec->Print();

  std::string c_jpsi_diff_xsec_image_name = image_dir;
  c_jpsi_diff_xsec_image_name.append("c_jpsi_diff_xsec");
  c_jpsi_diff_xsec_image_name.append(".png");
  c_jpsi_diff_xsec->Print(c_jpsi_diff_xsec_image_name.c_str() , "png");
  c_jpsi_diff_xsec->Close();

  TCanvas* c_hstack = new TCanvas("c_hstack","J/#psi p_{T}", 1400, 700); c_jpsi_diff_xsec->Divide(2,1);
  THStack* hs = new THStack("hs","Prompt Differential Cross Section");
  THStack* hs_nonprompt = new THStack("hs_nonprompt","NOnPrompt Differential Cross Section");
  //THStack hs_nonprompt("hs_nonprompt","hs_nonprompt");
  c_hstack->cd(1);
  h_zjpsi_dps_prompt->SetFillColor(kGreen - 2);
  h_zjpsi_theory_prompt->SetFillColor(kMagenta - 2);

  h_zjpsi_dps_prompt->GetYaxis()->SetLimits(1E-11, 1E-5);
  h_zjpsi_dps_prompt->GetXaxis()->SetLimits(8.5, 100);

  hs->Add(h_zjpsi_dps_prompt);
  hs->Add(h_zjpsi_theory_prompt);
  //hs.Add(h_zjpsi_cms_prompt);
  //hs.Add(h_zjpsi_atlas_prompt);


  hs->Draw("BAR");
  h_zjpsi_cms_prompt->Draw("E1same");
  //h_zjpsi_cms_prompt->Draw("E1");
  //hs.Draw("BARSAME");
  //hs.Draw("BAR");
  //h_zjpsi_cms_prompt->Draw("E1same");
  h_zjpsi_atlas_prompt->Draw("E1same");

  //hs->GetYaxis()->SetLimits(1E-11, 1E-5);
  hs->SetMaximum(1E-5);
  hs->SetMinimum(1E-11);
  //hs->GetXaxis()->SetLimits(10.5, 100);

  //h_zjpsi_cms_prompt->GetYaxis()->SetRangeUser(1E-11, 1E-5);
  //h_zjpsi_cms_prompt->GetXaxis()->SetRangeUser(8.5, 100);


  Double_t xl1_l7=.60, yl1_l7=0.70, xl2_l7=xl1_l7+.3, yl2_l7=yl1_l5+.2;
  TLegend *leg7 = new TLegend(xl1_l7,yl1_l7,xl2_l7,yl2_l7);
  leg7->SetFillColor(kWhite);
  leg7->AddEntry(h_zjpsi_cms_prompt,"CMS prompt","l");
  leg7->AddEntry(h_zjpsi_atlas_prompt,"ATLAS prompt","l");
  leg7->AddEntry(h_zjpsi_dps_prompt,"DPS prompt","l");
  leg7->AddEntry(h_zjpsi_theory_prompt,"Theory prompt","l");
  leg7->SetShadowColor(0);
  leg7->SetFillStyle(0);
  leg7->SetBorderSize(0);
  leg7->SetLineWidth(1);
  leg7->SetNColumns(1);
  leg7->SetTextFont(42);
  leg7->SetTextSize(0.03);

  leg7->Draw();

  mark.DrawLatex(0.735,0.957,"19.7 fb^{-1} (8 TeV)");
  mark.DrawLatex(0.195,0.89,"CMS");
  mark.DrawLatex(0.195,0.86,"#it{Preliminary}");


  c_hstack->cd(1)->SetLogy(1);
  c_hstack->cd(1)->SetLogx(1);

  // plot the non-prompt part
  c_hstack->cd(2);
  //hs_nonprompt.Add(h_zjpsi_dps_nonprompt);
  //hs_nonprompt.Add(h_zjpsi_cms_nonprompt);
  //hs_nonprompt.Add(h_zjpsi_atlas_nonprompt);
  //hs_nonprompt.Draw();

  h_zjpsi_cms_nonprompt->SetTitle("Nonprompt Differential Cross Section");

  h_zjpsi_dps_nonprompt->SetFillColor(kGreen - 2);
  //h_zjpsi_theory_nonprompt->SetFillColor(kMagenta - 2);
  h_zjpsi_theory_nonprompt->SetFillColorAlpha(kMagenta - 2, 0.35);

  h_zjpsi_dps_nonprompt->GetYaxis()->SetLimits(1E-11, 1E-5);
  h_zjpsi_dps_nonprompt->GetXaxis()->SetLimits(8.5, 100);

  hs_nonprompt->Add(h_zjpsi_dps_nonprompt);
  hs_nonprompt->Add(h_zjpsi_theory_nonprompt);
  hs_nonprompt->Draw("BAR");
  h_zjpsi_cms_nonprompt->Draw("E1same");
  h_zjpsi_atlas_nonprompt->Draw("E1same");

  hs_nonprompt->SetMaximum(1E-5);
  hs_nonprompt->SetMinimum(1E-11);

  //h_zjpsi_cms_nonprompt->SetLineColor(kBlue - 2);
  //h_zjpsi_cms_nonprompt->Draw("E1");

  //h_zjpsi_atlas_nonprompt->SetLineColor(kRed - 2);
  //h_zjpsi_dps_nonprompt->SetLineColor(kGreen - 2);
  //h_zjpsi_atlas_nonprompt->Draw("E1same");
  //h_zjpsi_dps_nonprompt->SetFillColor(kGreen - 2);
  //h_zjpsi_dps_nonprompt->Draw("BARsame");

  Double_t xl1_l8=.60, yl1_l8=0.70, xl2_l8=xl1_l6+.3, yl2_l8=yl1_l6+.2;
  TLegend *leg8 = new TLegend(xl1_l8,yl1_l8,xl2_l8,yl2_l8);
  leg8->SetFillColor(kWhite);
  leg8->AddEntry(h_zjpsi_cms_nonprompt,"CMS nonprompt","l");
  leg8->AddEntry(h_zjpsi_atlas_nonprompt,"ATLAS nonprompt","l");
  leg8->AddEntry(h_zjpsi_dps_nonprompt,"DPS nonprompt","l");
  leg8->AddEntry(h_zjpsi_theory_nonprompt,"DY MC nonprompt","l");
  leg8->SetShadowColor(0);
  leg8->SetFillStyle(0);
  leg8->SetBorderSize(0);
  leg8->SetLineWidth(1);
  leg8->SetNColumns(1);
  leg8->SetTextFont(42);
  leg8->SetTextSize(0.03);
  leg8->Draw();

  mark.DrawLatex(0.735,0.957,"19.7 fb^{-1} (8 TeV)");
  mark.DrawLatex(0.195,0.89,"CMS");
  mark.DrawLatex(0.195,0.86,"#it{Preliminary}");

  //h_zjpsi_cms_nonprompt->GetYaxis()->SetRangeUser(1E-11, 1E-5);
  //h_zjpsi_cms_nonprompt->GetXaxis()->SetRangeUser(8.5, 100);

  c_hstack->cd(2)->SetLogy(1);
  c_hstack->cd(2)->SetLogx(1);

  std::string c_hstack_image_name = image_dir;
  c_hstack_image_name.append("c_hstack");
  c_hstack_image_name.append(".png");
  c_hstack->Print(c_hstack_image_name.c_str() , "png");
  c_hstack->Close();

}
