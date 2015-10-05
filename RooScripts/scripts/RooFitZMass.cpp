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

#include "BW_CB_pdf_class.h"

using namespace RooFit;

int RooFitZMass(
    const std::string& DATA_FILE_1,
    const std::string& DATA_FILE_2,
    const bool USE_Z_TO_EE,
    const std::string& OUT_DIR
    ) {
  // Open the data file
  TFile* f_data_1 = new TFile(DATA_FILE_1.c_str(), "READ");
  if (f_data_1 == NULL) {
    std::cout << "Data file1 is invalid" << std::endl;
    return 1;
  }
  TFile* f_data_2 = new TFile(DATA_FILE_2.c_str(), "READ");
  if (f_data_2 == NULL) {
    std::cout << "Data file2 is invalid" << std::endl;
    return 1;
  }
  // Pass the open files to the main RooFitter
  const int RET_CODE = RooFitZMass(f_data_1, f_data_2, USE_Z_TO_EE, OUT_DIR);

  // Clean up and return the exit code
  delete f_data_1;
  delete f_data_2;

  return RET_CODE;
}

int RooFitZMass(
    TFile* const DATA_FILE_1,
    TFile* const DATA_FILE_2,
    const bool USE_Z_TO_EE,
    const std::string& OUT_DIR
    ) {
  // Constants
  //const int N_CPU = 8;
  const int N_CPU = 1;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  //TODO testing
  double z_mass_min = 50.0;
  double z_mass_max = 150.0;
  //double z_mass_min = 60.0;
  //double z_mass_max = 120.0;
  // Set up the variables we're going to read in from the files
  RooRealVar z_mass("z_mass", "z_mass" , z_mass_min, z_mass_max, "GeV");
  RooRealVar zjpsi_mass("zjpsi_mass", "zjpsi_mass" , z_mass_min, z_mass_max, "GeV");
  //TODO clean up
  //z_mass.setRange("negative",-20.0,-3.0);
  //z_mass.setRange("positive",3.0,20.0);
  //// Define a range named "signal" in z_mass from -1,1
  //z_mass.setRange("signal",-1.0,1.0) ;
  //z_mass.setRange("test",-20.0,20.0) ;

  
  std::string inclusive_z_hist = "";
  std::string z_hist_name = "";
  std::string inclusive_zjpsi_hist = "";
  std::string zjpsi_hist_name = "";
  if (USE_Z_TO_EE) {
    inclusive_z_hist.append( "ZFinder/Z_To_Electrons/z Mass: Coarse" );
    z_hist_name.append("ztoee_mass");
    //inclusive_zjpsi_hist.append( "ZFinder/Z_To_Electrons_And_Jpsi/z Mass: Coarse" );
    //zjpsi_hist_name.append("ztoee_jpsi_mass");

    //TODO TESTING
    //inclusive_zjpsi_hist.append( "ZFinder/Z_To_Electrons_And_Jpsi/z Mass: Coarse" );
    inclusive_zjpsi_hist.append( "ZFinder/Z_To_Electrons_And_Jpsi/z Mass: Coarse" );

    zjpsi_hist_name.append("ztoee_jpsi_mass");
  }
  else {
    inclusive_z_hist.append( "ZFinder/Z_To_Muons/Z From Muons Mass: Coarse" );
    z_hist_name.append("ztomumu_mass");
    //inclusive_zjpsi_hist.append( "ZFinder/Z_To_Muons_And_Jpsi/Z From Muons Mass: Coarse" );
    //zjpsi_hist_name.append("ztomumu_jpsi_mass");
    
    //inclusive_zjpsi_hist.append( "ZFinder/Z_To_Muons_And_Jpsi/Z From Muons Mass: Coarse" );
    //TODO TESTING
    inclusive_zjpsi_hist.append( "ZFinder/Z_To_Muons_And_Jpsi/Z From Muons Mass: Coarse" );


    zjpsi_hist_name.append("ztomumu_jpsi_mass");
  }

  TH1D *h_z_mass = (TH1D*) DATA_FILE_1->Get( inclusive_z_hist.c_str() );
  TH1D *h_zjpsi_mass = (TH1D*) DATA_FILE_2->Get( inclusive_zjpsi_hist.c_str() );
  h_z_mass->Sumw2();
  h_z_mass->Rebin(1);

  RooDataHist z_mass_data_hist("z_mass_data_hist", z_hist_name.c_str(), z_mass, h_z_mass);
  RooDataHist zjpsi_mass_data_hist("zjpsi_mass_data_hist", zjpsi_hist_name.c_str(), zjpsi_mass, h_zjpsi_mass);

  //RooRealVar mean("mean", "mean", 91.0, 60.0, 120.0);
  //RooRealVar sigma("sigma", "sigma", 10.0, 2, 20);
  //RooRealVar alpha("alpha", "alpha", 1.8, 1.5, 2.5);
  //RooRealVar n("n", "n", 2., 0.0, 8.);
  //RooCBShape crystal_ball ("crystal_ball", "crystal_ball", z_mass, mean, sigma, alpha, n );

  RooRealVar mean("mean","mean", 91.0, 60, 120);
  RooRealVar width("width","width", 5.0, 0.0, 10.0);
  RooRealVar sigma("sigma","sigma", 5.0, 0.0, 10.0);
  RooVoigtian voigtian("voigtian","voigtian", z_mass, mean, width, sigma);  

  RooRealVar alpha("alpha","alpha",60.,30,500.);
  RooRealVar gamma("gamma","gamma",0.01,0.001,10.0);
  RooRealVar delta("delta","delta",10.,3.,2000.);
  RooFormulaVar var1("var1","(alpha-z_mass)/delta",RooArgSet(alpha,z_mass,delta));
  RooFormulaVar var2("var2","-1.0*gamma*z_mass",RooArgSet(gamma,z_mass));
  RooGenericPdf MyBackgroundPdf("MyBackgroundPdf","ROOT::Math::erfc(var1)*exp(var2)",RooArgSet(var1, var2));

  RooRealVar signal_fraction("signal_fraction", "signal_fraction" , 0.1 , 0.0, 1.);
  RooAddPdf z_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(voigtian, MyBackgroundPdf), RooArgList(signal_fraction));

  RooFitResult *z_fitres = z_mass_fitpdf.fitTo(z_mass_data_hist, Range(z_mass_min, z_mass_max), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());


  ////RooAbsPdf *model_pdf = new BW_CB_pdf_class(RooRealVar& z_mass, TString name="bw_res", TString title="BW and CB convoluted pdf");
  //RooRealVar z_mass2("z_mass2", "z_mass2" , z_mass_min, z_mass_max, "GeV");
  //RooDataHist z_mass2_data_hist("z_mass2_data_hist", z_hist_name.c_str(), z_mass2, h_z_mass);
  //RooAbsPdf *model_pdf = new BW_CB_pdf_class(RooRealVar& z_mass2, TString name="bw_res", TString title="BW and CB convoluted pdf");
  //void SetFitType(0);

  //RooFitResult *fitres_MC = model_pdf->fitTo(z_mass2_data_hist,RooFit::Save(),
  //    RooFit::NumCPU(1), // numcpu=1
  //    RooFit::Verbose(kFALSE),RooFit::PrintLevel(-1),
  //    RooFit::Warnings(kFALSE),RooFit::PrintEvalErrors(kFALSE),
  //    RooFit::SumW2Error(kFALSE) // _isMCSumW2=true if the fit is unbinned
  //    );
  //fitres_MC->Print("V");

  //TODO testing
  double zjpsi_mean_value = mean.getVal();
  double zjpsi_width_value = width.getVal();
  double zjpsi_sigma_value = sigma.getVal();
  double zjpsi_alpha_value = alpha.getVal();
  double zjpsi_gamma_value = gamma.getVal();
  double zjpsi_delta_value = delta.getVal();
  //double zjpsi_var1_value = var1.getVal();
  //double zjpsi_var2_value = var2.getVal();

  RooRealVar zjpsi_mean("zjpsi_mean","zjpsi_mean", zjpsi_mean_value );
  RooRealVar zjpsi_width("zjpsi_width","zjpsi_width", zjpsi_width_value);
  RooRealVar zjpsi_sigma("zjpsi_sigma","zjpsi_sigma", zjpsi_sigma_value);
  RooVoigtian zjpsi_voigtian("zjpsi_voigtian","zjpsi_voigtian", zjpsi_mass, zjpsi_mean, zjpsi_width, zjpsi_sigma);  
  RooRealVar zjpsi_alpha("zjpsi_alpha","zjpsi_alpha", zjpsi_alpha_value);
  RooRealVar zjpsi_gamma("zjpsi_gamma","zjpsi_gamma", zjpsi_gamma_value);
  RooRealVar zjpsi_delta("zjpsi_delta","zjpsi_delta", zjpsi_delta_value);
  RooFormulaVar zjpsi_var1("zjpsi_var1","(zjpsi_alpha-zjpsi_mass)/zjpsi_delta",RooArgSet(zjpsi_alpha,zjpsi_mass,zjpsi_delta));
  RooFormulaVar zjpsi_var2("zjpsi_var2","-1.0*zjpsi_gamma*zjpsi_mass",RooArgSet(zjpsi_gamma,zjpsi_mass));
  RooGenericPdf zjpsi_MyBackgroundPdf("zjpsi_MyBackgroundPdf","ROOT::Math::erfc(zjpsi_var1)*exp(zjpsi_var2)",RooArgSet(zjpsi_var1, zjpsi_var2));
  RooRealVar zjpsi_signal_fraction("zjpsi_signal_fraction", "zjpsi_signal_fraction" , 0.1 , 0.0, 1.);

  RooAddPdf zjpsi_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(zjpsi_voigtian, zjpsi_MyBackgroundPdf), RooArgList(zjpsi_signal_fraction));

  RooFitResult *zjpsi_fitres = zjpsi_mass_fitpdf.fitTo(zjpsi_mass_data_hist, Range(z_mass_min, z_mass_max), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());




  std::cout << z_hist_name << std::endl;
  if (USE_Z_TO_EE) {
    std::cout << "z to ee" << std::endl;
  }
  else {
    std::cout << "z to mumu" << std::endl;
  }
  z_fitres->Print();

  std::cout << "zjpsi z mass fit: " << std::endl;
  zjpsi_fitres->Print();

  TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 2000, 750);

  canvas1->cd(1);
  //gPad->SetLogy();

  RooPlot* z_mass_fitframe;
  if (USE_Z_TO_EE) { 
    z_mass_fitframe = z_mass.frame( Title("Z->ee") );
  }
  else {
    z_mass_fitframe = z_mass.frame( Title("Z->mumu") );
  }
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
  canvas1->Print(z_image_name.c_str() , "png");
  canvas1->Close();

  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 2000, 750);

  canvas2->cd(1);
  //gPad->SetLogy();

  RooPlot* zjpsi_mass_fitframe;
  if (USE_Z_TO_EE) { 
    //zjpsi_mass_fitframe = zjpsi_mass.frame( Title("Z->ee + Jpsi: Z mass") );
    zjpsi_mass_fitframe = zjpsi_mass.frame( Title("Z->ee + Jpsi: Z mass") );
  }
  else {
    //zjpsi_mass_fitframe = zjpsi_mass.frame( Title("Z->mumu + Jpsi: Z mass") );
    zjpsi_mass_fitframe = zjpsi_mass.frame( Title("Z->mumu + Jpsi: Z mass") );
  }
  zjpsi_mass_data_hist.plotOn(zjpsi_mass_fitframe);
  zjpsi_mass_fitpdf.plotOn(zjpsi_mass_fitframe, Components(zjpsi_voigtian), LineColor(kGreen-2));
  zjpsi_mass_fitpdf.plotOn(zjpsi_mass_fitframe, Components(zjpsi_MyBackgroundPdf), LineColor(kBlue-2));
  zjpsi_mass_fitpdf.plotOn(zjpsi_mass_fitframe, LineColor(kRed-2));
  //zjpsi_mass_fitframe->SetMinimum(0.5);
  //z_mass_fitframe->SetMaximum(5e4);

  zjpsi_mass_fitframe->Draw();

  std::string zjpsi_image_name = OUT_DIR;
  zjpsi_image_name.append(zjpsi_hist_name);
  zjpsi_image_name.append(".png");
  canvas2->Print(zjpsi_image_name.c_str() , "png");
  canvas2->Close();

  return 0;
}

int main(int argc, char* argv[]) {
  const int ARGC = 5;
  if (argc < ARGC) {
    std::cout << "Not enough arguments.";
    return 1;
  } else if (argc > ARGC) {
    std::cout << "Too many arguments.";
    return 1;
  } else {
    /* Read in arguments */
    const std::string DATA_FILE_1(argv[1]);
    const std::string DATA_FILE_2(argv[2]);
    const std::string OUT_DIR(argv[4]);
    bool USE_Z_TO_EE = true;
    std::istringstream ss(argv[3]);
    if (!(ss >> USE_Z_TO_EE ) ) {
      std::cout << "Invalid bool " << argv[3] << std::endl;
      return 1;
    }
    RooFitZMass(DATA_FILE_1, DATA_FILE_2, USE_Z_TO_EE, OUT_DIR);
    return 0;
  }
}
