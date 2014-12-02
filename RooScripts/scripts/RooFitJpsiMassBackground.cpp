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
#include "RooFitJpsiMassBackground.h"
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

int RooFitJpsiMassBackground(
    const std::string& DATA_FILE_1,
    const bool USE_Z_TO_EE,
    const std::string& OUT_DIR
    ) {
  // Open the data file
  TFile* f_data_1 = new TFile(DATA_FILE_1.c_str(), "READ");
  if (f_data_1 == NULL) {
    std::cout << "Data file is invalid" << std::endl;
    return 1;
  }
  // Pass the open files to the main RooFitter
  const int RET_CODE = RooFitJpsiMassBackground(f_data_1, USE_Z_TO_EE, OUT_DIR);

  // Clean up and return the exit code
  delete f_data_1;

  return RET_CODE;
}

int RooFitJpsiMassBackground(
    TFile* const DATA_FILE_1,
    const bool USE_Z_TO_EE,
    const std::string& OUT_DIR
    ) {
  // Constants
  //const int N_CPU = 8;
  const int N_CPU = 1;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  double dimuon_mass_min = 2.5;
  double dimuon_mass_max = 5.0;
  double dimuon_mass_signal_min = 3.0;
  double dimuon_mass_signal_max = 3.2;
  // Set up the variables we're going to read in from the files
  RooRealVar dimuon_mass("dimuon_mass", "dimuon_mass" , dimuon_mass_min, dimuon_mass_max, "GeV");
  dimuon_mass.setRange("low", dimuon_mass_min, dimuon_mass_signal_min);
  dimuon_mass.setRange("high", dimuon_mass_signal_max, dimuon_mass_max);
  dimuon_mass.setRange("signal", dimuon_mass_signal_min, dimuon_mass_signal_max) ;
  dimuon_mass.setRange("all", dimuon_mass_min, dimuon_mass_max) ;

  std::string inclusive_jpsi_hist = "";
  if(USE_Z_TO_EE) {
    inclusive_jpsi_hist.append("ZFinder/Z_To_Electrons_And_Good_Dimuon_Jpsi/");
  }
  else {
    inclusive_jpsi_hist.append("ZFinder/Z_To_Muons_And_Good_Dimuon_Jpsi/");
  }
  inclusive_jpsi_hist.append( "jpsi_mass" );
  std::string jpsi_hist_name = "";
  jpsi_hist_name.append( "dimuon_mass_background");

  TH1D *h_dimuon_mass = (TH1D*) DATA_FILE_1->Get( inclusive_jpsi_hist.c_str() );
  h_dimuon_mass->Sumw2();
  h_dimuon_mass->Rebin(5);

  RooDataHist dimuon_mass_data_hist("dimuon_mass_data_hist", jpsi_hist_name.c_str(), dimuon_mass, h_dimuon_mass);

  RooRealVar slope("slope", "slope", -0.1, -10., 10.);
  RooExponential bg_exponential("bg_exponential", "bg_exponential", dimuon_mass, slope);

  RooFitResult *jpsi_fitres = bg_exponential.fitTo(dimuon_mass_data_hist, Range("low,high"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());
  // R e t r i e v e   r a w  &   n o r m a l i z e d   v a l u e s   o f   R o o F i t   p . d . f . s
  // --------------------------------------------------------------------------------------------------

  // Return 'raw' unnormalized value of bg_exponential
  //std::cout << "bg_exponential = " << bg_exponential.getVal() << std::endl ;
  
  // Return value of bg_exponential normalized over distance_z in range of dimuon_mass
  RooArgSet nset(dimuon_mass) ;
  //std::cout << "bg_exponential_Norm[dimuon_mass] = " << bg_exponential.getVal(&nset) << std::endl ;

  // Create object representing integral over bg_exponential
  // which is used to calculate  bg_exponential_Norm[dimuon_mass] == bg_exponential / bg_exponential_Int[dimuon_mass]
  RooAbsReal* ibg_exponential = bg_exponential.createIntegral(dimuon_mass) ;
  //std::cout << "bg_exponential_Int[dimuon_mass] = " << ibg_exponential->getVal() << std::endl ;


  // I n t e g r a t e   n o r m a l i z e d   p d f   o v e r   s u b r a n g e
  // ----------------------------------------------------------------------------

  
  // Create an integral of bg_exponential_Norm[dimuon_mass] over dimuon_mass in range "signal"
  // This is the fraction of p.d.f. bg_exponential_Norm[dimuon_mass] which is in the
  // range named "signal"
  RooAbsReal* ibg_exponential_sig = bg_exponential.createIntegral(dimuon_mass,NormSet(dimuon_mass),Range("signal")) ;
  RooAbsReal* ibg_exponential_high = bg_exponential.createIntegral(dimuon_mass,NormSet(dimuon_mass),Range("high")) ;
  RooAbsReal* ibg_exponential_low = bg_exponential.createIntegral(dimuon_mass,NormSet(dimuon_mass),Range("low")) ;
  RooAbsReal* ibg_exponential_all = bg_exponential.createIntegral(dimuon_mass,NormSet(dimuon_mass),Range("all")) ;
  std::cout << "bg_exponential_Int[dimuon_mass|signal]_Norm[dimuon_mass] = " << ibg_exponential_sig->getVal() << std::endl ;
  std::cout << "bg_exponential_Int[dimuon_mass|high]_Norm[dimuon_mass] = " << ibg_exponential_high->getVal() << std::endl ;
  std::cout << "bg_exponential_Int[dimuon_mass|low]_Norm[dimuon_mass] = " << ibg_exponential_low->getVal() << std::endl ;
  std::cout << "bg_exponential_Int[dimuon_mass|all]_Norm[dimuon_mass] = " << ibg_exponential_all->getVal() << std::endl ;

  //double bg_events = 0.;
  TAxis *axis = h_dimuon_mass->GetXaxis();

  int bmin_low = axis->FindBin(dimuon_mass_min);
  int bmax_low = axis->FindBin(dimuon_mass_signal_min);
  double integral_low = h_dimuon_mass->Integral(bmin_low,bmax_low);
  integral_low -= h_dimuon_mass->GetBinContent(bmin_low)*(dimuon_mass_min - axis->GetBinLowEdge(bmin_low)) / axis->GetBinWidth(bmin_low);
  integral_low -= h_dimuon_mass->GetBinContent(bmax_low)*(axis->GetBinUpEdge(bmax_low) - dimuon_mass_signal_min)  / axis->GetBinWidth(bmax_low); 
  std::cout << "integral low: " << integral_low << std::endl;

  int bmin_high = axis->FindBin(dimuon_mass_signal_max);
  int bmax_high = axis->FindBin(dimuon_mass_max);
  double integral_high = h_dimuon_mass->Integral(bmin_high,bmax_high);
  integral_high -= h_dimuon_mass->GetBinContent(bmin_high)*(dimuon_mass_signal_max - axis->GetBinLowEdge(bmin_high)) / axis->GetBinWidth(bmin_high);
  integral_high -= h_dimuon_mass->GetBinContent(bmax_high)*(axis->GetBinUpEdge(bmax_high) - dimuon_mass_max)  / axis->GetBinWidth(bmax_high); 
  std::cout << "integral high: " << integral_high << std::endl;
 
  //(i_high + i_low) * bg_in_signal_fraction / (1 - bg_in_signal_fraction)
  double bg_events_in_signal_region = (integral_high + integral_low ) * ibg_exponential_sig->getVal() / (1.0 - ibg_exponential_sig->getVal() ) ;
  std::cout << "bg_events_in_signal_region " << bg_events_in_signal_region << std::endl;
    
  jpsi_fitres->Print();

  TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 750);

  // Plot the left side
  canvas->cd(1);
  //gPad->SetLogy();
  RooPlot* dimuon_mass_fitframe;
  if (USE_Z_TO_EE) {
    dimuon_mass_fitframe = dimuon_mass.frame( Title("Z->ee + Dimuon" ));
  }
  else {
    dimuon_mass_fitframe = dimuon_mass.frame( Title("Z->mumu + Dimuon" ));
  }
  //dimuon_mass_fitframe->SetName(0); // Unset title
  dimuon_mass_data_hist.plotOn(dimuon_mass_fitframe);
  bg_exponential.plotOn(dimuon_mass_fitframe, LineColor(kRed-2));
  //dimuon_mass_fitframe->SetMinimum(0.5);
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
  const int ARGC = 4;
  if (argc < ARGC) {
    std::cout << "Not enough arguments.";
    return 1;
  } else if (argc > ARGC) {
    std::cout << "Too many arguments.";
    return 1;
  } else {
    /* Read in arguments */
    const std::string DATA_FILE_1(argv[1]);
    bool USE_Z_TO_EE = true;
    std::istringstream ss(argv[2]);
    if (!(ss >> USE_Z_TO_EE ) ) {
      std::cout << "Invalid bool " << argv[2] << std::endl;
      return 1;
    }
    const std::string OUT_DIR(argv[3]);
    RooFitJpsiMassBackground(DATA_FILE_1, USE_Z_TO_EE, OUT_DIR);
    return 0;
  }
}
