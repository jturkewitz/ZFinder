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

// RooFit
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitZMassBackground.h"
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

using namespace RooFit;

int RooFitZMassBackground(
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
  const int RET_CODE = RooFitZMassBackground(f_data_1, USE_Z_TO_EE, OUT_DIR);

  // Clean up and return the exit code
  delete f_data_1;

  return RET_CODE;
}

int RooFitZMassBackground(
    TFile* const DATA_FILE_1,
    const bool USE_Z_TO_EE,
    const std::string& OUT_DIR
    ) {
  // Constants
  //const int N_CPU = 8;
  const int N_CPU = 1;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  //TODO testing
  ////double z_mass_min = 80.0;
  ////double z_mass_max = 100.0;
  //double z_mass_min = 40.0;
  //double z_mass_max = 150.0;
  //// Set up the variables we're going to read in from the files
  //RooRealVar z_mass("z_mass", "z_mass" , z_mass_min, z_mass_max, "GeV");

  double z_mass_min = 40.0;
  double z_mass_min_upper;
  if (USE_Z_TO_EE) {
    z_mass_min_upper = 50;
  }
  else {
    z_mass_min_upper = 45;
  }
  //double z_mass_min_upper = 45.0;
  double z_mass_max = 300.0;
  double z_mass_max_lower = 150.0;
  RooRealVar z_mass("z_mass", "z_mass" , z_mass_min, z_mass_max, "GeV");
  //z_mass.setRange("low",z_mass_min,50.0);
  z_mass.setRange("low",z_mass_min,z_mass_min_upper);
  z_mass.setRange("high",z_mass_max_lower,z_mass_max);
  z_mass.setRange("signal",80.0,100.0) ;
  z_mass.setRange("test",z_mass_min,z_mass_max) ;

  
  std::string inclusive_z_hist = "";
  std::string z_hist_name = "";
  if (USE_Z_TO_EE) {
    //inclusive_z_hist.append( "ZFinder/Z_To_Electrons/z Mass: Coarse" );
    //inclusive_z_hist.append( "ZFinder/Dielectron_Z_Good_Compatible_Vertex/z Mass: Coarse" );
    inclusive_z_hist.append( "ZFinder/Dielectron_Z_Good_Compatible_Vertex/z Mass: All" );
    z_hist_name.append("ztoee_mass");
  }
  else {
    //inclusive_z_hist.append( "ZFinder/Z_To_Muons/Z From Muons Mass: Coarse" );
    //inclusive_z_hist.append( "ZFinder/Dimuon_Z_Good_Compatible_Vertex/Z From Muons Mass: Coarse" );
    inclusive_z_hist.append( "ZFinder/Dimuon_Z_Good_Compatible_Vertex/Z From Muons Mass: All" );
    z_hist_name.append("ztomumu_mass");
  }

  TH1D *h_z_mass = (TH1D*) DATA_FILE_1->Get( inclusive_z_hist.c_str() );
  h_z_mass->Sumw2();
  h_z_mass->Rebin(1);

  RooDataHist z_mass_data_hist("z_mass_data_hist", z_hist_name.c_str(), z_mass, h_z_mass);

  //RooRealVar mean("mean","mean", 91.0, 60, 120);
  //RooRealVar width("width","width", 5.0, 0.0, 10.0);
  //RooRealVar sigma("sigma","sigma", 5.0, 0.0, 10.0);
  //RooVoigtian voigtian("voigtian","voigtian", z_mass, mean, width, sigma);  

  //RooRealVar alpha("alpha","alpha",60.,30,500.);
  //RooRealVar gamma("gamma","gamma",0.01,0.001,10.0);
  //RooRealVar delta("delta","delta",10.,3.,2000.);
  //RooFormulaVar var1("var1","(alpha-z_mass)/delta",RooArgSet(alpha,z_mass,delta));
  //RooFormulaVar var2("var2","-1.0*gamma*z_mass",RooArgSet(gamma,z_mass));
  ////RooGenericPdf MyBackgroundPdf("MyBackgroundPdf","ROOT::Math::erfc(var1)*exp(var2)",RooArgSet(var1, var2));
  //RooGenericPdf bg_exponential("bg_exponential","ROOT::Math::erfc(var1)*exp(var2)",RooArgSet(var1, var2));

  RooRealVar alpha("alpha","alpha",60.,10,500.);
  RooRealVar gamma("gamma","gamma",0.1,0.001,10.0);
  RooRealVar delta("delta","delta",10.,3.,2000.);
  RooFormulaVar var1("var1","(alpha-z_mass)/delta",RooArgSet(alpha,z_mass,delta));
  RooFormulaVar var2("var2","-1.0*gamma*z_mass",RooArgSet(gamma,z_mass));
  //RooGenericPdf MyBackgroundPdf("MyBackgroundPdf","ROOT::Math::erfc(var1)*exp(var2)",RooArgSet(var1, var2));
  RooGenericPdf bg_exponential("bg_exponential","ROOT::Math::erfc(var1)*exp(var2)",RooArgSet(var1, var2));


  // Make pdf
  //RooRealVar slope("slope","slope",-0.1,-10.0,-10.0);
  //RooExponential model("model","model",z_mass,slope);
  //// D e f i n e   e f f i c i e n c y   f u n c t i o n
  //// ---------------------------------------------------
  //// Use error function to simulate turn-on slope
  ////RooFormulaVar eff("eff","0.5*(TMath::Erf((z_mass-1)/0.5)+1)",z_mass);
  //RooFormulaVar eff("eff","0.5*(ROOT::TMath::Erf((z_mass-1)/0.5)+1)",z_mass);
  //// D e f i n e   d e c a y   p d f   w i t h   e f f i c i e n c y 
  //// ---------------------------------------------------------------
  //// Multiply pdf(t) with efficiency in t
  //RooEffProd bg_exponential("bg_exponential","model with efficiency",model,eff);





  //RooRealVar slope("slope", "slope", -0.1, -10., 10.);
  //RooExponential bg_exponential("bg_exponential", "bg_exponential", z_mass, slope);

  //RooRealVar signal_fraction("signal_fraction", "signal_fraction" , 0.1 , 0.0, 1.);
  //RooAddPdf z_mass_fitpdf("z_mass_fitpdf", "z_mass_fitpdf", RooArgList(voigtian, MyBackgroundPdf), RooArgList(signal_fraction));
  

  //RooRealVar slope("slope", "slope", -0.1, -10., 10.);
  //RooRealVar slope("slope", "slope", 0.1, 0.01, 10.);
  //RooGenericPdf bg_exponential("bg_exponential","(1 + slope/z_mass)",RooArgSet(z_mass,slope));

  //RooFitResult *z_fitres = z_mass_fitpdf.fitTo(z_mass_data_hist, Range("low,high"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());
  
  
  RooFitResult *z_fitres = bg_exponential.fitTo(z_mass_data_hist, Range("low,high"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());
  //RooFitResult *z_fitres = MyBackgroundPdf.fitTo(z_mass_data_hist, Range("low,high"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());
  
  
  //RooFitResult *z_fitres = bg_exponential.fitTo(z_mass_data_hist, Range("high"), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());

  ////////////////
  
  //TODO this is just to plot over full range, should be better way to do this
  //double jpsi_gaussmean_value = gauss_pileup_mean.getVal();
  //double jpsi_gausssigma_value = gauss_pileup_sigma.getVal();


  //double slope_value = slope.getVal();
  //RooRealVar slope_test("slope_test", "slope value", slope_value);
  //RooExponential bg_exponential_test("bg_exponential_test", "BG Exponential test", z_mass, slope_test );

  //double slope_value = slope.getVal();
  //RooRealVar slope_test("slope_test", "slope value", slope_value);
  //RooGenericPdf bg_exponential_test("bg_exponential","1 / (z_mass*slope)",RooArgSet(z_mass,slope_test));
  //RooGenericPdf bg_exponential_test("bg_exponential","(1 + slope_test/z_mass)",RooArgSet(z_mass,slope_test));

  //double alpha_value = alpha.getVal();
  //double gamma_value = gamma.getVal();
  //double delta_value = delta.getVal();
  //std::cout << alpha_value << " " << gamma_value << " " << delta_value << std::endl;
  //RooFormulaVar var1_new("var1_new","(alpha_value-z_mass)/delta_value",RooArgSet(z_mass));
  //RooFormulaVar var2_new("var2_new","-1.0*gamma_value*z_mass",RooArgSet(z_mass));

  RooRealVar alpha_value("alpha_value","alpha_value",alpha.getVal());
  RooRealVar gamma_value("gamma_value","gamma_value",gamma.getVal());
  RooRealVar delta_value("delta_value","delta_value",delta.getVal());
  RooFormulaVar var1_new("var1_new","(alpha_value-z_mass)/delta_value",RooArgSet(alpha_value,z_mass,delta_value));
  RooFormulaVar var2_new("var2_new","-1.0*gamma_value*z_mass",RooArgSet(gamma_value,z_mass));
  //TODO fix this hack 43.7494 0.0202972 5.66097
  //RooFormulaVar var1_new("var1_new","(43.7494-z_mass)/5.66097",RooArgSet(z_mass));
  //RooFormulaVar var2_new("var2_new","-1.0*0.0202972*z_mass",RooArgSet(z_mass));
  RooGenericPdf bg_exponential_test("bg_exponential_test","ROOT::Math::erfc(var1_new)*exp(var2_new)",RooArgSet(var1_new, var2_new));



  // R e t r i e v e   r a w  &   n o r m a l i z e d   v a l u e s   o f   R o o F i t   p . d . f . s
  // --------------------------------------------------------------------------------------------------

  // Return 'raw' unnormalized value of pileup_gauss
  std::cout << "exp = " << bg_exponential.getVal() << std::endl ;
  
  // Return value of pileup_gauss normalized over distance_z in range [-10,10]
  RooArgSet nset(z_mass) ;
  std::cout << "bg_exp_Norm[z_mass] = " << bg_exponential.getVal(&nset) << std::endl ;

  // Create object representing integral over bg_exponentail
  // which is used to calculate  pileup_gauss_Norm[distance_z] == pileup_gauss / pileup_gauss_Int[distance_z]
  RooAbsReal* ibg_exponential = bg_exponential.createIntegral(z_mass) ;
  std::cout << "bg_exponential_Int[z_mass] = " << ibg_exponential->getVal() << std::endl ;


  // I n t e g r a t e   n o r m a l i z e d   p d f   o v e r   s u b r a n g e
  // ----------------------------------------------------------------------------

  
  // Create an integral of pileup_gauss_Norm[distance_z] over distance_z in range "signal"
  // This is the fraction of of p.d.f. pileup_gauss_Norm[distance_z] which is in the
  // range named "signal"
  RooAbsReal* ibg_exponential_sig = bg_exponential.createIntegral(z_mass,NormSet(z_mass),Range("signal")) ;
  RooAbsReal* ibg_exponential_high = bg_exponential.createIntegral(z_mass,NormSet(z_mass),Range("high")) ;
  RooAbsReal* ibg_exponential_low = bg_exponential.createIntegral(z_mass,NormSet(z_mass),Range("low")) ;
  RooAbsReal* ibg_exponential_all = bg_exponential.createIntegral(z_mass,NormSet(z_mass),Range("test")) ;
  cout << "bg_exponential_Int[distance_z|signal]_Norm[z_mass] = " << ibg_exponential_sig->getVal() << endl ;
  cout << "bg_exponential_Int[distance_z|high]_Norm[z_mass] = " << ibg_exponential_high->getVal() << endl ;
  cout << "bg_exponential_Int[distance_z|low]_Norm[z_mass] = " << ibg_exponential_low->getVal() << endl ;
  cout << "bg_exponential_Int[distance_z|all]_Norm[z_mass] = " << ibg_exponential_all->getVal() << endl ;
  /////////

  TAxis *x_axis = h_z_mass->GetXaxis();
  TAxis *y_axis = h_z_mass->GetYaxis();
  double mass_min = z_mass_max_lower;
  double mass_max = z_mass_max;
  int x_bmin = x_axis->FindBin(mass_min);
  int x_bmax = x_axis->FindBin(mass_max);
  // -1 to fix overflow binning
  //double integral = h_z_dimuon_mass->Integral(x_bmin, x_bmax - 1);
  double integral = h_z_mass->Integral(x_bmin, x_bmax - 1);
  cout << "high sideband integral from " << mass_min << " to " << mass_max << " integral " << integral << endl;

  double mass_min_low = z_mass_min;
  double mass_max_low = z_mass_min_upper;
  int x_bmin_low = x_axis->FindBin(mass_min_low);
  int x_bmax_low = x_axis->FindBin(mass_max_low);
  // -1 to fix overflow binning
  //double integral = h_z_dimuon_mass->Integral(x_bmin, x_bmax - 1);
  double integral_low = h_z_mass->Integral(x_bmin_low, x_bmax_low - 1);
  cout << "low sideband integral from " << mass_min_low << " to " << mass_max_low << " integral_low " << integral_low << endl;

  double mass_min_signal = 80;
  double mass_max_signal = 100;
  int x_bmin_signal = x_axis->FindBin(mass_min_signal);
  int x_bmax_signal = x_axis->FindBin(mass_max_signal);
  // -1 to fix overfsignal binning
  //double integral = h_z_dimuon_mass->Integral(x_bmin, x_bmax - 1);
  double integral_signal = h_z_mass->Integral(x_bmin_signal, x_bmax_signal -1 );
  cout << "signal sideband integral from " << mass_min_signal << " to " << mass_max_signal << " integral_signal " << integral_signal << endl;

  std::cout << z_hist_name << std::endl;
  if (USE_Z_TO_EE) {
    std::cout << "z to ee" << std::endl;
  }
  else {
    std::cout << "z to mumu" << std::endl;
  }
  z_fitres->Print();

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
  //bg_exponential_test.plotOn(z_mass_fitframe, LineColor(kGreen-2), NormRange("low,high"));
  //bg_exponential_test.plotOn(z_mass_fitframe, LineColor(kGreen-2), NormRange("high"));
  bg_exponential_test.plotOn(z_mass_fitframe, LineColor(kGreen-2), NormRange("low,high"));
  bg_exponential.plotOn(z_mass_fitframe, LineColor(kBlue-2), VLines(), Range("signal"), RooFit::Name("signal"));
  bg_exponential.plotOn(z_mass_fitframe, LineColor(kRed-2), VLines(), RooFit::Name("background"));

  Double_t xl1=.58, yl1=0.55, xl2=xl1+.3, yl2=yl1+.325;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillColor(kWhite);
  leg->AddEntry(z_mass_fitframe->findObject("signal"),"signal region","l");
  leg->AddEntry(z_mass_fitframe->findObject("background"),"background region","l");

  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(voigtian), LineColor(kGreen-2));
  //z_mass_fitpdf.plotOn(z_mass_fitframe, Components(MyBackgroundPdf), LineColor(kBlue-2));
  //z_mass_fitpdf.plotOn(z_mass_fitframe, LineColor(kRed-2));
  z_mass_fitframe->SetMinimum(10);
  z_mass_fitframe->SetMaximum(10e6);

  gPad->SetLogy();
  z_mass_fitframe->Draw();
  leg->Draw();

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
    RooFitZMassBackground(DATA_FILE_1, USE_Z_TO_EE, OUT_DIR);
    return 0;
  }
}
