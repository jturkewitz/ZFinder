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

// RooFit
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitJpsiMass.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooKeysPdf.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"


using namespace RooFit;

int RooFitJpsiMass(
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
  const int RET_CODE = RooFitJpsiMass(f_data_1, OUT_DIR);

  // Clean up and return the exit code
  delete f_data_1;

  return RET_CODE;
}

int RooFitJpsiMass(
    TFile* const DATA_FILE_1,
    const std::string& OUT_DIR
    ) {
  // Constants
  //const int N_CPU = 8;
  const int N_CPU = 1;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  //double dimuon_mass_min = 2.8;
  //double dimuon_mass_max = 3.4;
  double dimuon_mass_min = 3.0;
  double dimuon_mass_max = 3.2;
  //double dimuon_mass_max = 3.2;
  // Set up the variables we're going to read in from the files
  RooRealVar dimuon_mass("dimuon_mass", "dimuon_mass" , dimuon_mass_min, dimuon_mass_max, "GeV");
  //TODO clean up
  //dimuon_mass.setRange("negative",-20.0,-3.0);
  //dimuon_mass.setRange("positive",3.0,20.0);
  //// Define a range named "signal" in dimuon_mass from -1,1
  //dimuon_mass.setRange("signal",-1.0,1.0) ;
  //dimuon_mass.setRange("test",-20.0,20.0) ;

  //std::string inclusive_jpsi_hist = "ZFinder/Jpsi/";
  //std::string inclusive_jpsi_hist = "ZFinder/Dimuon/";
  std::string inclusive_jpsi_hist = "ZFinder/Dimuon_Jpsi_Primary_Vertex/";
  //inclusive_jpsi_hist.append( "jpsi_mass" );
  inclusive_jpsi_hist.append( "jpsi_mass_pt15to20" );
  std::string jpsi_hist_name = "";
  jpsi_hist_name.append( "dimuon_mass");

  TH1D *h_dimuon_mass = (TH1D*) DATA_FILE_1->Get( inclusive_jpsi_hist.c_str() );
  h_dimuon_mass->Sumw2();
  h_dimuon_mass->Rebin(1);

  RooDataHist dimuon_mass_data_hist("dimuon_mass_data_hist", jpsi_hist_name.c_str(), dimuon_mass, h_dimuon_mass);

  //RooRealVar mean("mean", "mean", 3.1, 3.00, 3.18);
  //RooRealVar sigma("sigma", "sigma", 0.05, 0.001, 0.1);
  //RooRealVar alpha("alpha", "alpha", 1.8, 1.0, 2.5);
  //RooRealVar n("n", "n", 2., 1.0, 100.);
  ////RooRealVar n("n", "n", 2.0);
  //RooCBShape crystal_ball ("crystal_ball", "crystal_ball", dimuon_mass, mean, sigma, alpha, n );

  RooRealVar mean("mean", "mean", 3.1, 3.0, 3.2);
  RooRealVar sigma("sigma", "sigma", 0.1, 0.0001, 10.0);
  RooGaussian gauss ("gauss", "gauss", dimuon_mass, mean, sigma);

  RooRealVar slope("slope", "slope", -0.1, -100., 100.);
  RooExponential bg_exponential("bg_exponential", "bg_exponential", dimuon_mass, slope);

  RooRealVar signal_fraction("signal_fraction", "signal_fraction" , 0.9 , 0.0, 1.);
  //RooRealVar signal_fraction("signal_fraction", "signal_fraction" , 1.);
  //RooAddPdf dimuon_mass_fitpdf("dimuon_mass_fitpdf", "dimuon_mass_fitpdf", RooArgList(crystal_ball, bg_exponential), RooArgList(signal_fraction));
  RooAddPdf dimuon_mass_fitpdf("dimuon_mass_fitpdf", "dimuon_mass_fitpdf", RooArgList(gauss, bg_exponential), RooArgList(signal_fraction));

  RooFitResult *jpsi_fitres = dimuon_mass_fitpdf.fitTo(dimuon_mass_data_hist, Range(dimuon_mass_min, dimuon_mass_max), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());

  jpsi_fitres->Print();

  TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 750);

  // Plot the left side
  canvas->cd(1);
  gPad->SetLogy();
  //RooPlot* dimuon_mass_fitframe = dimuon_mass.frame( Title(jpsi_hist_name.c_str()) );
  RooPlot* dimuon_mass_fitframe = dimuon_mass.frame( Title("Inclusive J/Psi Trigger" ));
  //dimuon_mass_fitframe->SetName(0); // Unset title
  dimuon_mass_data_hist.plotOn(dimuon_mass_fitframe);
  //dimuon_mass_fitpdf.plotOn(dimuon_mass_fitframe, Components(voigtian), LineColor(kGreen-2));
  //dimuon_mass_fitpdf.plotOn(dimuon_mass_fitframe, Components(crystal_ball), LineColor(kGreen-2));
  dimuon_mass_fitpdf.plotOn(dimuon_mass_fitframe, Components(gauss), LineColor(kGreen-2));
  dimuon_mass_fitpdf.plotOn(dimuon_mass_fitframe, Components(bg_exponential), LineColor(kBlue-2));
  dimuon_mass_fitpdf.plotOn(dimuon_mass_fitframe, LineColor(kRed-2));
  dimuon_mass_fitframe->SetMinimum(0.5);
  //dimuon_mass_fitframe->SetMaximum(5e4);

  dimuon_mass_fitframe->Draw();

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
    RooFitJpsiMass(DATA_FILE_1, OUT_DIR);
    return 0;
  }
}
