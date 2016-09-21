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
#include "RooFit4LeptonMass.h"
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

#include "RooEffProd.h"

#include "BW_CB_pdf_class.h"
#include <TROOT.h>
#include <TStyle.h>
#include "tdrStyle.C"

using namespace RooFit;

int RooFit4LeptonMass(
    const std::string& DATA_FILE_1,
    const bool USE_Z_TO_EE,
    const std::string& OUT_DIR
    ) {
  // Open the data file
  TFile* f_data_1 = new TFile(DATA_FILE_1.c_str(), "READ");
  if (f_data_1 == NULL) {
    std::cout << "Data file1 is invalid" << std::endl;
    return 1;
  }
  // Pass the open files to the main RooFitter
  const int RET_CODE = RooFit4LeptonMass(f_data_1, USE_Z_TO_EE, OUT_DIR);

  // Clean up and return the exit code
  delete f_data_1;

  return RET_CODE;
}

int RooFit4LeptonMass(
    TFile* const DATA_FILE_1,
    const bool USE_Z_TO_EE,
    const std::string& OUT_DIR
    ) {
  const int N_CPU = 1;
  setTDRStyle();

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  double z_mass_min = 40.0;
  double z_mass_max = 200.0;
  RooRealVar z_mass("z_mass", "z_mass" , z_mass_min, z_mass_max, "GeV");
  z_mass.setRange("signal",80.0,100.0) ;
  z_mass.setRange("full",z_mass_min,z_mass_max) ;

  
  std::string z_hist = "";
  std::string z_hist_name = "";
  if (USE_Z_TO_EE) {
    z_hist.append( "ZFinder/Z_To_Electrons_And_Jpsi/four_lepton_mass" );
    z_hist_name.append("ztoee_4l_mass");
  }
  else {
    z_hist.append( "ZFinder/Z_To_Muons_And_Jpsi/four_lepton_mass" );
    z_hist_name.append("ztomumu_4l_mass");
  }

  TH1D *h_z_mass = (TH1D*) DATA_FILE_1->Get( z_hist.c_str() );
  h_z_mass->Sumw2();
  h_z_mass->Rebin(10);

  RooDataHist z_mass_data_hist("z_mass_data_hist", z_hist_name.c_str(), z_mass, h_z_mass);

  //RooRealVar mean("mean","mean", 91.0, 70, 150);
  //RooRealVar width("width","width", 5.0, 0.0, 10.0);
  //RooRealVar sigma("sigma","sigma", 5.0, 0.0, 10.0);
  //RooVoigtian voigtian("voigtian","voigtian", z_mass, mean, width, sigma);  

  //RooRealVar mean("mean","mean", 110.0, 70, 130);
  RooRealVar mean_low("mean_low","mean_low", 91.2);
  RooRealVar sigma_low("sigma_low","sigma_low", 3.5, 2.0, 6.0);
  RooGaussian gauss_low("gauss_low","gauss_low", z_mass, mean_low, sigma_low);  

  RooRealVar mean_high("mean_high","mean_high", 110,70,130);
  RooRealVar sigma_high("sigma_high","sigma_high", 5.0, 0.01, 40.0);
  RooGaussian gauss_high("gauss_high","gauss_high", z_mass, mean_high, sigma_high);

  //RooRealVar slope("slope", "slope", -0.1, -50., 50.);
  //RooExponential bg_exponential("bg_exponential", "bg_exponential", z_mass, slope);

  RooRealVar signal_fraction("signal_fraction", "signal_fraction" , 0.5 , 0.0, 1.);
  //RooAddPdf z_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(voigtian, bg_exponential), RooArgList(signal_fraction));
  RooAddPdf z_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(gauss_low, gauss_high), RooArgList(signal_fraction));

  RooFitResult *z_fitres = z_mass_fitpdf.fitTo(z_mass_data_hist, Range(z_mass_min, z_mass_max), NumCPU(N_CPU), Verbose(true), PrintLevel(-1), SumW2Error(kFALSE), Save());
  z_fitres->Print();
  
  TAxis *x_axis = h_z_mass->GetXaxis();
  TAxis *y_axis = h_z_mass->GetYaxis();
  double mass_min = z_mass_min;
  double mass_max = z_mass_max;
  int x_bmin = x_axis->FindBin(mass_min);
  int x_bmax = x_axis->FindBin(mass_max);
  // -1 to fix overflow binning
  //double integral = h_z_dimuon_mass->Integral(x_bmin, x_bmax - 1);
  double integral = h_z_mass->Integral(x_bmin, x_bmax - 1);
  cout << "integral from " << mass_min << " to " << mass_max << " integral " << integral << endl;

  double mass_min_signal = 80;
  double mass_max_signal = 100;
  int x_bmin_signal = x_axis->FindBin(mass_min_signal);
  int x_bmax_signal = x_axis->FindBin(mass_max_signal);
  // -1 to fix overfsignal binning
  //double integral = h_z_dimuon_mass->Integral(x_bmin, x_bmax - 1);
  double integral_signal = h_z_mass->Integral(x_bmin_signal, x_bmax_signal -1 );
  cout << "signal region integral from " << mass_min_signal << " to " << mass_max_signal << " integral_signal " << integral_signal << endl;

  RooAbsReal* ibg_exponential_sig = gauss_low.createIntegral(z_mass,NormSet(z_mass),Range("signal")) ;
  std::cout << "gauss_low_signal_fraction = " << ibg_exponential_sig->getVal() << std::endl ;
  std::cout << "(fraction_gauss_low * gauss_low_80_100_frac * events_to)" << std::endl;
  std::cout << "gauss_low signal events  " << ibg_exponential_sig->getVal() * signal_fraction.getVal() * integral << std::endl;

  std::cout << z_hist_name << std::endl;
  if (USE_Z_TO_EE) {
    std::cout << "z to ee" << std::endl;
  }
  else {
    std::cout << "z to mumu" << std::endl;
  }

  //TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 2000, 750);
  TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 2000, 1000);

  canvas1->cd(1);
  //gPad->SetLogy();

  RooPlot* z_mass_fitframe;
  if (USE_Z_TO_EE) { 
    z_mass_fitframe = z_mass.frame( Title("Z->ee") );
    z_mass_fitframe->GetYaxis()->SetTitle("M_{ee} (GeV)");
  }
  else {
    z_mass_fitframe = z_mass.frame( Title("Z->mumu") );
    z_mass_fitframe->GetYaxis()->SetTitle("M_{#mu#mu} (GeV)");
  }
  z_mass_fitframe->GetYaxis()->SetTitle("Events / 5 GeV");

  z_mass_data_hist.plotOn(z_mass_fitframe);

  z_mass_fitpdf.plotOn(z_mass_fitframe, LineColor(kRed-2), RooFit::Name("total"));
  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(bg_exponential),  LineColor(kGreen-2), RooFit::Name("exponential"));
  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(voigtian), LineColor(kBlue-2), RooFit::Name("voigtian"));
  z_mass_fitpdf.plotOn(z_mass_fitframe, Components(gauss_low),  LineColor(kGreen-2), RooFit::Name("gauss_low"));
  z_mass_fitpdf.plotOn(z_mass_fitframe, Components(gauss_high), LineColor(kBlue-2), RooFit::Name("gauss_high"));

  //Double_t xl1=.58, yl1=0.55, xl2=xl1+.3, yl2=yl1+.325;
  //TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  //leg->SetFillColor(kWhite);
  //leg->AddEntry(z_mass_fitframe->findObject("total"),"total","l");
  //leg->AddEntry(z_mass_fitframe->findObject("gauss_low"),"gauss_low","l");
  //leg->AddEntry(z_mass_fitframe->findObject("gauss_high"),"gauss_high","l");

  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(voigtian), LineColor(kGreen-2));
  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(MyBackgroundPdf), LineColor(kBlue-2));
  //z_mass_fitpdf.plotOn(z_mass_fitframe, LineColor(kRed-2));


  //z_mass_fitframe->SetMinimum(0.5);
  //z_mass_fitframe->SetMaximum(10e2);

  //gPad->SetLogy();
  z_mass_fitframe->Draw();
  if (USE_Z_TO_EE) { 
    z_mass_fitframe->GetXaxis()->SetTitle("M_{ee} (GeV)");
  }
  else {
    z_mass_fitframe->GetXaxis()->SetTitle("M_{#mu#mu} (GeV)");
  }
  z_mass_fitframe->GetYaxis()->SetTitleOffset(0.8);
  //gStyle->SetPadRightMargin(0.35); 
  //gROOT->ForceStyle();

  Double_t xl1=.8, yl1=0.55, xl2=xl1+.3, yl2=yl1+.325;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillColor(kWhite);
  leg->AddEntry(z_mass_fitframe->findObject("total"),"total","l");
  leg->AddEntry(z_mass_fitframe->findObject("gauss_low"),"signal","l");
  leg->AddEntry(z_mass_fitframe->findObject("gauss_high"),"background","l");
  leg->SetFillColor(kWhite);
  
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
  mark.DrawLatex(0.845,0.957,"19.7 fb^{-1} (8 TeV)");
  mark.DrawLatex(0.195,0.89,"CMS");
  mark.DrawLatex(0.195,0.86,"#it{Preliminary}");

  std::string z_image_name = OUT_DIR;
  z_image_name.append(z_hist_name);
  z_image_name.append(".png");
  canvas1->Print(z_image_name.c_str() , "png");
  canvas1->Close();


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
    const std::string OUT_DIR(argv[3]);
    bool USE_Z_TO_EE = true;
    std::istringstream ss(argv[2]);
    if (!(ss >> USE_Z_TO_EE ) ) {
      std::cout << "Invalid bool " << argv[2] << std::endl;
      return 1;
    }
    RooFit4LeptonMass(DATA_FILE_1, USE_Z_TO_EE, OUT_DIR);
    return 0;
  }
}
