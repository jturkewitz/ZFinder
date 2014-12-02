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
#include "RooFitZMass.h"
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

int RooFitZMass(
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
  const int RET_CODE = RooFitZMass(f_data_1, OUT_DIR);

  // Clean up and return the exit code
  delete f_data_1;

  return RET_CODE;
}

int RooFitZMass(
    TFile* const DATA_FILE_1,
    const std::string& OUT_DIR
    ) {
  // Constants
  //const int N_CPU = 8;
  const int N_CPU = 1;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  double z_mass_min = 50.0;
  //double z_mass_min = 2.85;
  //double z_mass_max = 3.15;
  double z_mass_max = 130.0;
  //double z_mass_max = 3.2;
  // Set up the variables we're going to read in from the files
  RooRealVar z_mass("z_mass", "z_mass" , z_mass_min, z_mass_max, "GeV");
  //TODO clean up
  //z_mass.setRange("negative",-20.0,-3.0);
  //z_mass.setRange("positive",3.0,20.0);
  //// Define a range named "signal" in z_mass from -1,1
  //z_mass.setRange("signal",-1.0,1.0) ;
  //z_mass.setRange("test",-20.0,20.0) ;

  
  //std::string inclusive_z_hist = "ZFinder/Combined Single Reco/6 60 < M_{ee} < 120/";
  std::string inclusive_z_hist = "ZFinder/Dielectron_Z_Good_Compatible_Vertex/";

  //std::string inclusive_z_hist = "ZFinder/Dimuon_Z_Good_Compatible_Vertex/";
  //inclusive_z_hist.append( "z Mass: Fine" );

  inclusive_z_hist.append( "z Mass: Coarse" );
  //inclusive_z_hist.append( "Z0 Mass: Coarse" );

  //inclusive_z_hist.append( "Z From Muons Mass: Coarse" );
  std::string z_hist_name = "";
  z_hist_name.append( "z_mass");

  TH1D *h_z_mass = (TH1D*) DATA_FILE_1->Get( inclusive_z_hist.c_str() );
  h_z_mass->Sumw2();
  h_z_mass->Rebin(1);

  RooDataHist z_mass_data_hist("z_mass_data_hist", z_hist_name.c_str(), z_mass, h_z_mass);

  //RooRealVar mean("mean", "mean", 91.0, 60.0, 120.0);
  ////RooRealVar sigma("sigma", "sigma", 0.1, 0.0001, 10.0);
  //RooRealVar sigma("sigma", "sigma", 10.0, 2, 20);
  //RooRealVar alpha("alpha", "alpha", 1.8, 1.5, 2.5);
  //RooRealVar n("n", "n", 2., 0.0, 8.);
  //RooCBShape crystal_ball ("crystal_ball", "crystal_ball", z_mass, mean, sigma, alpha, n );

  RooRealVar mean("mean","mean", 95.0, 60.0, 120.0);
  RooRealVar width("width","width", 5.0, 0.0, 10.0);
  RooRealVar sigma("sigma","sigma", 5.0, 0.0, 10.0);
  RooVoigtian voigtian("voigtian","voigtian", z_mass, mean, width, sigma);  


  //RooRealVar mean("mean", "mean", 3.1, 3.0, 3.2);
  //RooRealVar sigma("sigma", "sigma", 0.1, 0.0001, 1.0);
  //RooGaussian gauss("gauss", "Gaussian of the Mass", z_mass, mean, sigma);


  //RooRealVar mean("mean", "mean", 3.1, 3.0, 3.2);
  //RooRealVar sigma("sigma", "sigma", 0.1, 0.0001, 10.0);
  //RooGaussian crystal_ball ("crystal_ball", "crystal_ball", z_mass, mean, sigma);

  //RooRealVar sigma2("sigma2", "sigma2", 0.1, 0.0001, 1.0);
  //RooRealVar alpha2("alpha2", "alpha2", 0., 0., 100);
  //RooRealVar n2("n2", "n2", 0., 2., 5.);
  //RooCBShape crystal_ball2 ("crystal_ball2", "crystal_ball2", z_mass, mean, sigma2, alpha2, n2 );

  ////RooRealVar cball_fraction("cball_fraction", "cball_fraction", 0., 0.0, 1.0);
  //RooRealVar cball_fraction("cball_fraction", "cball_fraction", 1.);
  //RooAddPdf crystal_balls("crystal_balls", "crystal_balls", RooArgList( crystal_ball, crystal_ball2), RooArgList(cball_fraction));

  RooRealVar slope("slope", "slope", -0.1, -10., 10.);
  RooExponential bg_exponential("bg_exponential", "bg_exponential", z_mass, slope);

  RooRealVar alpha("alpha","alpha",60.,30,500.);
  RooRealVar gamma("gamma","gamma",0.01,0.001,10.0);
  RooRealVar delta("delta","delta",10.,3.,20.);
  RooFormulaVar var1("var1","(alpha-z_mass)/delta",RooArgSet(alpha,z_mass,delta));
  RooFormulaVar var2("var2","-1.0*gamma*z_mass",RooArgSet(gamma,z_mass));
  RooGenericPdf MyBackgroundPdf("MyBackgroundPdf","ROOT::Math::erfc(var1)*exp(var2)",RooArgSet(var1, var2));

  RooRealVar signal_fraction("signal_fraction", "signal_fraction" , 0.1 , 0.0, 1.);
  //RooRealVar signal_fraction("signal_fraction", "signal_fraction" , 1.);
  //RooAddPdf z_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(crystal_ball, bg_exponential), RooArgList(signal_fraction));
  //RooAddPdf z_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(voigtian, bg_exponential), RooArgList(signal_fraction));
  RooAddPdf z_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(voigtian, MyBackgroundPdf), RooArgList(signal_fraction));
  //RooAddPdf z_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(gauss, bg_exponential), RooArgList(signal_fraction));

  RooFitResult *z_fitres = z_mass_fitpdf.fitTo(z_mass_data_hist, Range(z_mass_min, z_mass_max), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());

  z_fitres->Print();

  TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 750);

  // Plot the left side
  canvas->cd(1);
  gPad->SetLogy();
  RooPlot* z_mass_fitframe = z_mass.frame( Title(z_hist_name.c_str()) );
  //z_mass_fitframe->SetName(0); // Unset title
  z_mass_data_hist.plotOn(z_mass_fitframe);
  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(crystal_ball), LineColor(kGreen-2));
  z_mass_fitpdf.plotOn(z_mass_fitframe, Components(voigtian), LineColor(kGreen-2));
  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(gauss), LineColor(kGreen-2));
  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(bg_exponential), LineColor(kBlue-2));
  z_mass_fitpdf.plotOn(z_mass_fitframe, Components(MyBackgroundPdf), LineColor(kBlue-2));
  z_mass_fitpdf.plotOn(z_mass_fitframe, LineColor(kRed-2));
  z_mass_fitframe->SetMinimum(0.5);
  //z_mass_fitframe->SetMaximum(5e4);

  z_mass_fitframe->Draw();

  std::string z_image_name = OUT_DIR;
  z_image_name.append(z_hist_name);
  z_image_name.append(".png");
  canvas->Print(z_image_name.c_str() , "png");
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
    RooFitZMass(DATA_FILE_1, OUT_DIR);
    return 0;
  }
}
