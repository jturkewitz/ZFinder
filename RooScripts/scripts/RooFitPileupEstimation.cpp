// Standard Library
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>  // std::vector

// ROOT
#include <TCanvas.h>
#include "TH1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"

// RooFit
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitPileupEstimation.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooKeysPdf.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include <TROOT.h>
#include <TStyle.h>
#include "tdrStyle.C"


using namespace RooFit;

int RooFitPileupEstimation(
    const std::string& DATA_FILE_1,
    const std::string& OUT_DIR
    ) {
  // Open the data file
  TFile* f_data_1 = new TFile(DATA_FILE_1.c_str(), "READ");
  if (f_data_1 == NULL) {
    std::cout << "Data file is invalid" << std::endl;
    return 1;
  }
  // Pass the open files to the main RooFitter
  const int RET_CODE = RooFitPileupEstimation(f_data_1, OUT_DIR);

  // Clean up and return the exit code
  delete f_data_1;

  return RET_CODE;
}

int RooFitPileupEstimation(
    TFile* const DATA_FILE_1,
    const std::string& OUT_DIR
    ) {
  // Constants
  //const int N_CPU = 8;
  const int N_CPU = 1;
  setTDRStyle();

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  double distance_z_min = -20.0;
  double distance_z_max = 20.0;
  // Set up the variables we're going to read in from the files
  RooRealVar distance_z("distance_z", "#Delta z" , distance_z_min, distance_z_max, "cm");
  //TODO clean up, use variable names
  distance_z.setRange("negative",distance_z_min,-3.0);
  distance_z.setRange("positive",3.0,distance_z_max);
  // Define a range named "signal" in distance_z from -1,1
  double signal_min = -0.5;
  double signal_max = 0.5;
  distance_z.setRange("signal",signal_min,signal_max) ;
  double range1_min = -0.5;
  double range1_max = 0.5;
  distance_z.setRange("range1",range1_min,range1_max) ;
  double range2_min = -1.0;
  double range2_max = 1.0;
  distance_z.setRange("range2",range2_min,range2_max) ;
  double range3_min = -2.0;
  double range3_max = 2.0;
  distance_z.setRange("range3",range3_min,range3_max) ;
  double range4_min = -3.0;
  double range4_max = 3.0;
  distance_z.setRange("range4",range4_min,range4_max) ;

  //distance_z.setRange("signal",-0.5,0.5) ;
  //distance_z.setRange("signal",-0.5,0.5) ;
  distance_z.setRange("test",distance_z_min,distance_z_max) ;

  //std::string inclusive_jpsi_hist = "ZFinder/Jpsi/";
  std::string inclusive_jpsi_hist = "ZFinder/Dimuon_Jpsi_Vertex_Compatible/";
  //std::string inclusive_jpsi_hist = "ZFinder/Prompt_Jpsi/";
  inclusive_jpsi_hist.append( "jpsi_vtx_z - z_vtx_z" );
  std::string jpsi_hist_name = "";
  jpsi_hist_name.append( "distance_z");

  TH1D *h_jpsi_vertex_distance_z = (TH1D*) DATA_FILE_1->Get( inclusive_jpsi_hist.c_str() );
  h_jpsi_vertex_distance_z->Sumw2();
  h_jpsi_vertex_distance_z->Rebin(4);

  RooDataHist distance_z_data_hist("distance_z_data_hist", jpsi_hist_name.c_str(), distance_z, h_jpsi_vertex_distance_z);

  RooRealVar gauss_pileup_mean("gauss_prompt_mean", "Mean of the Pileup Gaussian", 0., -10., 10.);
  RooRealVar gauss_pileup_sigma("gauss_prompt_sigma", "Width of the Pileup Gaussian", 1., 0.005, 15.);
  RooGaussian pileup_gauss("pileup_gauss", "Gaussian of the Pileup", distance_z, gauss_pileup_mean, gauss_pileup_sigma);


  //RooFitResult *jpsi_fitres = pileup_gauss.fitTo(distance_z_data_hist, Range(distance_z_min, distance_z_max), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), Save());

//  RooFitResult *jpsi_fitres_low = pileup_gauss.fitTo(distance_z_data_hist, Range("negative"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), Save());
//  RooFitResult *jpsi_fitres_high = pileup_gauss.fitTo(distance_z_data_hist, Range("positive"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), Save());

  //RooFitResult *jpsi_fitres = pileup_gauss.fitTo(distance_z_data_hist, Range("negative, positive"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), Save());
  RooFitResult *jpsi_fitres = pileup_gauss.fitTo(distance_z_data_hist, Range("negative,positive"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), Save());
  
  
  //TODO this is just to plot over full range, should be better way to do this
  double jpsi_gaussmean_value = gauss_pileup_mean.getVal();
  double jpsi_gausssigma_value = gauss_pileup_sigma.getVal();

  RooRealVar gauss_pileup_mean_test("gauss_prompt_mean_test", "Mean of the Pileup Gaussian", jpsi_gaussmean_value);
  RooRealVar gauss_pileup_sigma_test("gauss_prompt_sigma_test", "Width of the Pileup Gaussian", jpsi_gausssigma_value);
  RooGaussian pileup_gauss_test("pileup_gauss_test", "Gaussian of the Pileup_test", distance_z, gauss_pileup_mean_test, gauss_pileup_sigma_test );

  // R e t r i e v e   r a w  &   n o r m a l i z e d   v a l u e s   o f   R o o F i t   p . d . f . s
  // --------------------------------------------------------------------------------------------------

  // Return 'raw' unnormalized value of pileup_gauss
  std::cout << "pileup_gauss = " << pileup_gauss.getVal() << std::endl ;
  
  // Return value of pileup_gauss normalized over distance_z in range [-10,10]
  RooArgSet nset(distance_z) ;
  //std::cout << "pileup_gauss_Norm[distance_z] = " << pileup_gauss.getVal(&nset) << std::endl ;

  // Create object representing integral over pileup_gauss
  // which is used to calculate  pileup_gauss_Norm[distance_z] == pileup_gauss / pileup_gauss_Int[distance_z]
  RooAbsReal* ipileup_gauss = pileup_gauss.createIntegral(distance_z) ;
  //std::cout << "pileup_gauss_Int[distance_z] = " << ipileup_gauss->getVal() << std::endl ;


  // I n t e g r a t e   n o r m a l i z e d   p d f   o v e r   s u b r a n g e
  // ----------------------------------------------------------------------------

  
  // Create an integral of pileup_gauss_Norm[distance_z] over distance_z in range "signal"
  // This is the fraction of of p.d.f. pileup_gauss_Norm[distance_z] which is in the
  // range named "signal"
  RooAbsReal* ipileup_gauss_sig = pileup_gauss.createIntegral(distance_z,NormSet(distance_z),Range("signal")) ;
  RooAbsReal* ipileup_gauss_range1 = pileup_gauss.createIntegral(distance_z,NormSet(distance_z),Range("range1")) ;
  RooAbsReal* ipileup_gauss_range2 = pileup_gauss.createIntegral(distance_z,NormSet(distance_z),Range("range2")) ;
  RooAbsReal* ipileup_gauss_range3 = pileup_gauss.createIntegral(distance_z,NormSet(distance_z),Range("range3")) ;
  RooAbsReal* ipileup_gauss_range4 = pileup_gauss.createIntegral(distance_z,NormSet(distance_z),Range("range4")) ;
  RooAbsReal* ipileup_gauss_positive = pileup_gauss.createIntegral(distance_z,NormSet(distance_z),Range("positive")) ;
  RooAbsReal* ipileup_gauss_negative = pileup_gauss.createIntegral(distance_z,NormSet(distance_z),Range("negative")) ;
  RooAbsReal* ipileup_gauss_all = pileup_gauss.createIntegral(distance_z,NormSet(distance_z),Range("test")) ;
  //cout << "Integral signal = " << signal_min << " to " << signal_max << " norm " << ipileup_gauss_sig->getVal() << endl ;
  cout << "Integral range1 = " << range1_min << " to " << range1_max << " norm " << ipileup_gauss_range1->getVal() << endl ;
  cout << "Integral range2 = " << range2_min << " to " << range2_max << " norm " << ipileup_gauss_range2->getVal() << endl ;
  cout << "Integral range3 = " << range3_min << " to " << range3_max << " norm " << ipileup_gauss_range3->getVal() << endl ;
  cout << "Integral range4 = " << range4_min << " to " << range4_max << " norm " << ipileup_gauss_range4->getVal() << endl ;
  cout << "Integral positve 3 to Inf = " << ipileup_gauss_positive->getVal() << endl ;
  cout << "Integral negative 3 to Inf = " << ipileup_gauss_negative->getVal() << endl ;
  cout << "Integral all = " << ipileup_gauss_all->getVal() << endl ;

  std::cout << jpsi_hist_name << std::endl;
  //jpsi_fitres->Print("v");

  jpsi_fitres->Print();
//  jpsi_fitres_low->Print();
//  jpsi_fitres_high->Print();

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1200, 900);

  // Plot the left side
  canvas->cd(1);
  gPad->SetLogy();
  //RooPlot* distance_z_fitframe = distance_z.frame(-0.3,5.0);
  //RooPlot* distance_z_fitframe = distance_z.frame( Title(jpsi_hist_name.c_str()) , Range("negative, positive" ));
  //RooPlot* distance_z_fitframe = distance_z.frame( Title(jpsi_hist_name.c_str()) );
  RooPlot* distance_z_fitframe = distance_z.frame( Title("#Delta z Between Primary Vertex and J/#psi Vertex") );
  //distance_z_fitframe->SetName(0); // Unset title
  distance_z_data_hist.plotOn(distance_z_fitframe);
  pileup_gauss_test.plotOn(distance_z_fitframe, LineColor(kGreen-2), NormRange("negative,positive"));
  pileup_gauss.plotOn(distance_z_fitframe, LineColor(kRed-2));

  pileup_gauss.plotOn(distance_z_fitframe, LineColor(kBlue-2), VLines(), Range("signal"), RooFit::Name("signal"));
  pileup_gauss.plotOn(distance_z_fitframe, LineColor(kRed-2), VLines(), RooFit::Name("background"));
  distance_z_fitframe->SetMinimum(10);
  distance_z_fitframe->SetMaximum(1e7);

  distance_z_fitframe->GetYaxis()->SetTitle("Events / 0.08 cm");
  distance_z_fitframe->GetXaxis()->SetTitleSize(0.055);
  distance_z_fitframe->GetXaxis()->SetLabelSize(0.035);
  //distance_z_fitframe->GetXaxis()->SetNdivisions(5);
  distance_z_fitframe->GetXaxis()->SetTitle("#Delta z Between Primary Vertex and J/#psi Vertex (cm)");
  distance_z_fitframe->Draw();

  Double_t xl1=.7, yl1=0.75, xl2=xl1+.3, yl2=yl1+.15;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillColor(kWhite);
  leg->AddEntry(distance_z_fitframe->findObject("signal"),"signal region","l");
  leg->AddEntry(distance_z_fitframe->findObject("background"),"background region","l");
  
  leg->SetShadowColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineWidth(1);
  leg->SetNColumns(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->Draw();

  TLatex mark;
  mark.SetTextSize(0.035);
  mark.SetNDC(true);
  mark.DrawLatex(0.795,0.957,"19.7 fb^{-1} (8 TeV)");
  mark.DrawLatex(0.195,0.89,"CMS");
  mark.DrawLatex(0.195,0.86,"#it{Preliminary}");


  std::string jpsi_image_name = OUT_DIR;
  jpsi_image_name.append(jpsi_hist_name);
  jpsi_image_name.append(".png");
  canvas->Print(jpsi_image_name.c_str() , "png");
  canvas->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  const int ARGC = 3;
  if (argc < ARGC) {
    std::cout << "Not enough arguments.";
    return 1;
  } else if (argc > ARGC) {
    std::cout << "Too many arguments.";
    return 1;
  } else {
    /* Read in arguments */
    const std::string DATA_FILE_1(argv[1]);
    const std::string OUT_DIR(argv[2]);
    RooFitPileupEstimation(DATA_FILE_1, OUT_DIR);
    return 0;
  }
}
