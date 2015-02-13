// Standard Library
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>  // std::vector

// ROOT
#include <TCanvas.h>
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLegend.h"

// RooFit
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitLifetimeAndMassCrystalBall.h"
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
#include "RooNumIntConfig.h"


using namespace RooFit;

int RooFitLifetimeAndMassCrystalBall(
    const std::string& DATA_FILE_1,
    const std::string& DATA_FILE_2,
    const bool USE_Z_TO_EE,
    const std::string& OUT_DIR
    ) {
  // Open the data file
  TFile* f_data_1 = new TFile(DATA_FILE_1.c_str(), "READ");
  TFile* f_data_2 = new TFile(DATA_FILE_2.c_str(), "READ");
  if (f_data_1 == NULL) {
    std::cout << "Data file_1 is invalid" << std::endl;
    return 1;
  }
  if (f_data_2 == NULL) {
    std::cout << "Data file_2 is invalid" << std::endl;
    return 1;
  }
  // Pass the open files to the main RooFitter
  const int RET_CODE = RooFitLifetimeAndMassCrystalBall(f_data_1, f_data_2, USE_Z_TO_EE, OUT_DIR);

  // Clean up and return the exit code
  delete f_data_1;
  delete f_data_2;

  return RET_CODE;
}

int RooFitLifetimeAndMassCrystalBall(
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
  //Increase default precision of numeric integration
  // as this exercise has high sensitivity to numeric integration precision
  //RooAbsPdf::defaultIntegratorConfig()->setEpsRel(1e-8) ;
  //RooAbsPdf::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
  RooAbsPdf::defaultIntegratorConfig()->setEpsRel(1e-10) ;
  RooAbsPdf::defaultIntegratorConfig()->setEpsAbs(1e-10) ;

  double dimuon_mass_min = 2.85;
  double dimuon_mass_max = 3.35;
  double dimuon_mass_signal_min = 3.0;
  double dimuon_mass_signal_max = 3.2;
  // Set up the variables we're going to read in from the files
  RooRealVar dimuon_mass("dimuon_mass", "dimuon_mass" , dimuon_mass_min, dimuon_mass_max, "GeV");
  dimuon_mass.setRange("low", dimuon_mass_min, dimuon_mass_signal_min);
  dimuon_mass.setRange("high", dimuon_mass_signal_max, dimuon_mass_max);
  dimuon_mass.setRange("signal", dimuon_mass_signal_min, dimuon_mass_signal_max) ;
  dimuon_mass.setRange("all", dimuon_mass_min, dimuon_mass_max) ;

  double tau_xy_min = -0.3;
  double tau_xy_max = 5.0;
  RooRealVar tau_xy("tau_xy", "tau_xy" , tau_xy_min, tau_xy_max, "ps");


  RooRealVar zjpsi_dimuon_mass("zjpsi_dimuon_mass", "zjpsi_dimuon_mass" , dimuon_mass_min, dimuon_mass_max, "GeV");
  RooRealVar zjpsi_tau_xy("zjpsi_tau_xy", "zjpsi_tau_xy" , tau_xy_min, tau_xy_max, "ps");
    
  std::string inclusive_jpsi_hist = "";
  inclusive_jpsi_hist.append("ZFinder/Dimuon_Jpsi_Primary_Vertex/");
  inclusive_jpsi_hist.append("dimuon_mass_vs_dimuon_tau_xy_fine");

  std::string inclusive_jpsi_hist_name = "";
  //inclusive_jpsi_hist_name.append( "inclusive_dimuon_mass_vs_dimuon_tau_xy");
  inclusive_jpsi_hist_name.append( "inclusive_dimuon_mass_vs_dimuon_tau_xy_fine");

  std::string zjpsi_hist = "";
  if(USE_Z_TO_EE) {
    zjpsi_hist.append("ZFinder/Z_To_Electrons_And_Good_Dimuon_Jpsi/");
  }
  else {
    zjpsi_hist.append("ZFinder/Z_To_Muons_And_Good_Dimuon_Jpsi/");
  }
  //zjpsi_hist.append("dimuon_mass_vs_dimuon_tau_xy");
  zjpsi_hist.append("dimuon_mass_vs_dimuon_tau_xy_fine");
  std::string zjpsi_hist_name = "";
  //zjpsi_hist_name.append( "z_and_dimuon_mass_vs_dimuon_tau_xy");
  zjpsi_hist_name.append( "z_and_dimuon_mass_vs_dimuon_tau_xy_fine");

  //2d fit - goal is to fit continuum_bg * prompt, continuum_bg* nonprompt, gauss * prompt (signal), gauss * nonprompt 

  TH2D *h_inclusive_dimuon_mass = (TH2D*) DATA_FILE_1->Get( inclusive_jpsi_hist.c_str() );
  TH2D *h_z_dimuon_mass = (TH2D*) DATA_FILE_2->Get( zjpsi_hist.c_str() );

  h_z_dimuon_mass->RebinY(5);

  //RooDataHist z_dimuon_mass_data_hist("z_dimuon_mass_data_hist", zjpsi_hist_name.c_str(), RooArgSet(dimuon_mass, tau_xy), h_z_dimuon_mass);
  RooDataHist dimuon_mass_data_hist("dimuon_mass_data_hist", inclusive_jpsi_hist_name.c_str(), RooArgSet(dimuon_mass, tau_xy), h_inclusive_dimuon_mass);
  RooDataHist zjpsi_dimuon_mass_data_hist("zjpsi_dimuon_mass_data_hist", zjpsi_hist_name.c_str(), RooArgSet(zjpsi_dimuon_mass, zjpsi_tau_xy), h_z_dimuon_mass);

  RooRealVar gaussmean("gaussmean", "Mean of the smearing Gaussian", 0.);
  RooRealVar gausssigma("gausssigma", "Width of the smearing Gaussian", 0.01, 0.005, 0.3);
  RooGaussModel smear_gauss_model("smear_gauss", "Gaussian used to smear the Exponential Decay", tau_xy, gaussmean, gausssigma);
  RooRealVar decay_lifetime("decay_lifetime", "Tau", 1.3, 0.05, 2.);
  RooDecay decay_exp("decay_exp", "Exponential Decay", tau_xy, decay_lifetime, smear_gauss_model, RooDecay::SingleSided);

  RooRealVar gauss_prompt_mean("gauss_prompt_mean", "Mean of the Prompt Gaussian", 0., -0.1, 0.1);
  RooRealVar gauss_prompt_sigma("gauss_prompt_sigma", "Width of the Prompt Gaussian", 0.007, 0.005, 0.25);
  RooGaussian prompt_gauss("prompt_gauss", "Gaussian of the Prompt Peak", tau_xy, gauss_prompt_mean, gauss_prompt_sigma);

  RooRealVar gauss_prompt_sigma_2("gauss_prompt_sigma_2", "Width of the Prompt Gaussian_2", 0.01, 0.008, 0.4);
  RooGaussian prompt_gauss_2("prompt_gauss_2", "Gaussian_2 of the Prompt Peak", tau_xy, gauss_prompt_mean, gauss_prompt_sigma_2);

  //TODO testing ---------------------------
  RooRealVar frac_prompt_sharp("frac_prompt_sharp", "frac_prompt_sharp" , 0.5 , 0.0, 1.);
//  TODO decide whether or not to use double gaussian fit, if so also figure out how to fix which one is prompt or not
//  RooRealVar prompt_sharp_fraction("prompt_sharp_fraction", "prompt_sharp_fraction" , 1.0);

  RooAddPdf tau_xy_gauss_sum_fitpdf("tau_xy_gauss_sum_fitpdf", "tau_xy_gauss_sum_fitpdf", RooArgList(prompt_gauss, prompt_gauss_2), RooArgList(frac_prompt_sharp));


  //TODO testing c-ball function
  //RooRealVar dimuon_mean("dimuon_mean", "dimuon_mean", 3.1, 3.0, 3.2);
  //RooRealVar dimuon_sigma("dimuon_sigma", "dimuon_sigma", 0.02, 0.001, 0.1);
  //RooGaussian dimuon_gauss ("dimuon_gauss", "dimuon_gauss", dimuon_mass, dimuon_mean, dimuon_sigma);

  //RooRealVar dimuon_slope("dimuon_slope", "dimuon_slope", -0.1, -10., 10.);
  //RooExponential dimuon_bg_exponential("dimuon_bg_exponential", "dimuon_bg_exponential", dimuon_mass, dimuon_slope);

  RooRealVar mean("mean", "mean", 3.1, 3.00, 3.18);
  RooRealVar sigma("sigma", "sigma", 0.05, 0.001, 0.1);
  RooRealVar alpha("alpha", "alpha", 1.8, 1.0, 2.5);
  RooRealVar n("n", "n", 2., 1.0, 100.);
  //RooRealVar n("n", "n", 2.0);
  RooCBShape crystal_ball ("crystal_ball", "crystal_ball", dimuon_mass, mean, sigma, alpha, n );

  RooRealVar dimuon_slope("dimuon_slope", "dimuon_slope", -0.1, -10., 10.);
  RooExponential dimuon_bg_exponential("dimuon_bg_exponential", "dimuon_bg_exponential", dimuon_mass, dimuon_slope);



  //RooRealVar mass_signal_tau_xy_signal_fraction("mass_signal_tau_xy_signal_fraction", "mass_signal_tau_xy_signal_fraction", 0.4, 0.0, 1.0);
  //RooRealVar mass_signal_tau_xy_bg_fraction("mass_signal_tau_xy_bg_fraction", "mass_signal_tau_xy_bg_fraction", 0.4, 0.0, 1.0);
  //RooRealVar mass_bg_tau_xy_bg_fraction("mass_bg_tau_xy_bg_fraction", "mass_bg_tau_xy_bg_fraction", 0.1, 0.0, 1.0);
  //RooRealVar mass_bg_tau_xy_signal_fraction("mass_bg_tau_xy_signal_fraction", "mass_bg_tau_xy_signal_fraction", 0.1, 0.0, 1.0);
  RooRealVar m_sig_tau_sig_frac("m_sig_tau_sig_frac", "m_sig_tau_sig_frac", 0.4, 0.0, 1.0);
  RooRealVar m_sig_tau_bg_frac("m_sig_tau_bg_frac", "m_sig_tau_bg_frac", 0.4, 0.0, 1.0);
  RooRealVar m_bg_tau_bg_frac("m_bg_tau_bg_frac", "m_bg_tau_bg_frac", 0.1, 0.0, 1.0);
  RooRealVar m_bg_tau_sig_frac("m_bg_tau_sig_frac", "m_bg_tau_sig_frac", 0.1, 0.0, 1.0);

  //RooProdPdf m_sig_tau_sig("m_sig_tau_sig", "m_sig_tau_sig", RooArgList(dimuon_gauss, tau_xy_gauss_sum_fitpdf ));
  //RooProdPdf m_sig_tau_bg("m_sig_tau_bg", "m_sig_tau_bg", RooArgList(dimuon_gauss, decay_exp ));
  RooProdPdf m_sig_tau_sig("m_sig_tau_sig", "m_sig_tau_sig", RooArgList(crystal_ball, tau_xy_gauss_sum_fitpdf ));
  RooProdPdf m_sig_tau_bg("m_sig_tau_bg", "m_sig_tau_bg", RooArgList(crystal_ball, decay_exp ));
  RooProdPdf m_bg_tau_bg("m_bg_tau_bg", "m_bg_tau_bg", RooArgList(dimuon_bg_exponential, decay_exp ));
  RooProdPdf m_bg_tau_sig("m_bg_tau_sig", "m_bg_tau_sig", RooArgList(dimuon_bg_exponential, tau_xy_gauss_sum_fitpdf ));


  //RooAddPdf mass_tau_xy_fitpdf("mass_tau_xy_fitpdf","mass_tau_xy_fitpdf",RooArgSet(mass_bg_tau_bg, mass_bg_tau_sig, mass_sig_tau_bg, mass_sig_tau_sig),
  //                                                                       RooArgList(mass_bg_tau_bg_frac, mass_bg_tau_sig_frac, mass_sig_tau_bg_frac, mass_sig_tau_sig_frac));
  //RooAddPdf mass_tau_xy_fitpdf("mass_tau_xy_fitpdf","mass_tau_xy_fitpdf",RooArgSet(mass_bg_tau_bg, mass_bg_tau_sig, mass_sig_tau_bg, mass_sig_tau_sig),
  //                                                                       RooArgList(mass_bg_tau_bg_frac, mass_bg_tau_sig_frac, mass_sig_tau_bg_frac));
  RooAddPdf mass_tau_xy_fitpdf("mass_tau_xy_fitpdf","mass_tau_xy_fitpdf",RooArgSet(m_sig_tau_sig, m_sig_tau_bg, m_bg_tau_sig, m_bg_tau_bg ),
                                                                         RooArgList(m_sig_tau_sig_frac, m_sig_tau_bg_frac, m_bg_tau_sig_frac));

  RooFitResult *jpsi_fitres = mass_tau_xy_fitpdf.fitTo(dimuon_mass_data_hist, NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());
  jpsi_fitres->Print();

  double zjpsi_gaussmean_value = gaussmean.getVal();
  double zjpsi_gausssigma_value = gausssigma.getVal();
  double zjpsi_decay_lifetime_value = decay_lifetime.getVal();
  double zjpsi_gauss_prompt_mean_value = gauss_prompt_mean.getVal();
  double zjpsi_gauss_prompt_sigma_value = gauss_prompt_sigma.getVal();
  double zjpsi_gauss_prompt_sigma_2_value = gauss_prompt_sigma_2.getVal();
  double zjpsi_frac_prompt_sharp_value = frac_prompt_sharp.getVal();
  //double zjpsi_gauss_prompt_sigma_error = gauss_prompt_sigma.getError();

  RooRealVar zjpsi_gaussmean("zjpsi_gaussmean", "Mean of the smearing zjpsi Gaussian", zjpsi_gaussmean_value);
  RooRealVar zjpsi_gausssigma("zjpsi_gausssigma", "Width of the smearing zjpsi Gaussian", zjpsi_gausssigma_value);
  RooGaussModel zjpsi_smear_gauss_model("zjpsi_smear_gauss", "Gaussian used to smear the zjpsi Exponential Decay", zjpsi_tau_xy, zjpsi_gaussmean, zjpsi_gausssigma);
  RooRealVar zjpsi_decay_lifetime("zjpsi_decay_lifetime", "zjpsi_Tau", zjpsi_decay_lifetime_value);
  RooDecay zjpsi_decay_exp("zjpsi_decay_exp", "zjpsi_Exponential Decay", zjpsi_tau_xy, zjpsi_decay_lifetime, zjpsi_smear_gauss_model, RooDecay::SingleSided);

  RooRealVar zjpsi_gauss_prompt_mean("zjpsi_gauss_prompt_mean", "Mean of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_mean_value);
  RooRealVar zjpsi_gauss_prompt_sigma("zjpsi_gauss_prompt_sigma", "Width of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_sigma_value);
  RooGaussian zjpsi_prompt_gauss("zjpsi_prompt_gauss", "Gaussian of the Prompt zjpsi_Peak", zjpsi_tau_xy, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma);

  RooRealVar zjpsi_gauss_prompt_sigma_2("zjpsi_gauss_prompt_sigma_2", "Width of the Prompt zjpsi_Gaussian_2", zjpsi_gauss_prompt_sigma_2_value);
  RooGaussian zjpsi_prompt_gauss_2("zjpsi_prompt_gauss_2", "Gaussian_2 of the Prompt zjpsi_Peak", zjpsi_tau_xy, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma_2);

  RooRealVar zjpsi_frac_prompt_sharp("zjpsi_frac_prompt_sharp", "zjpsi_frac_prompt_sharp" , zjpsi_frac_prompt_sharp_value);
  RooAddPdf zjpsi_tau_xy_gauss_sum_fitpdf("zjpsi_tau_xy_gauss_sum_fitpdf", "zjpsi_tau_xy_gauss_sum_fitpdf", RooArgList(zjpsi_prompt_gauss,zjpsi_prompt_gauss_2), RooArgList(zjpsi_frac_prompt_sharp));

  RooRealVar zjpsi_prompt_fraction("zjpsi_prompt_fraction", "zjpsi_prompt_fraction" , 0.01 , 0.0, 1);
  //RooAddPdf zjpsi_tau_xy_fitpdf("zjpsi_tau_xy_fitpdf", "zjpsi_tau_xy_fitpdf", RooArgList(zjpsi_tau_xy_gauss_sum_fitpdf, zjpsi_decay_exp), RooArgList(zjpsi_prompt_fraction));

  //double zjpsi_dimuon_mean_value = dimuon_mean.getVal();
  //double zjpsi_dimuon_sigma_value = dimuon_sigma.getVal();
  //double zjpsi_dimuon_slope_value = dimuon_slope.getVal();

  double zjpsi_mean_value = mean.getVal();
  double zjpsi_sigma_value = sigma.getVal();
  double zjpsi_alpha_value = alpha.getVal();
  double zjpsi_n_value = n.getVal();
  double zjpsi_dimuon_slope_value = dimuon_slope.getVal();

  //RooRealVar zjpsi_dimuon_mean("zjpsi_dimuon_mean", "zjpsi_dimuon_mean", zjpsi_dimuon_mean_value);
  //RooRealVar zjpsi_dimuon_sigma("zjpsi_dimuon_sigma", "zjpsi_dimuon_sigma", zjpsi_dimuon_sigma_value);
  //RooGaussian zjpsi_dimuon_gauss ("zjpsi_dimuon_gauss", "zjpsi_dimuon_gauss", zjpsi_dimuon_mass, zjpsi_dimuon_mean, zjpsi_dimuon_sigma);

  //RooRealVar zjpsi_dimuon_slope("zjpsi_dimuon_slope", "zjpsi_dimuon_slope", zjpsi_dimuon_slope_value);
  //RooExponential zjpsi_dimuon_bg_exponential("zjpsi_dimuon_bg_exponential", "zjpsi_dimuon_bg_exponential", zjpsi_dimuon_mass, zjpsi_dimuon_slope);

  RooRealVar zjpsi_mean("zjpsi_mean", "zjpsi_mean", zjpsi_mean_value);
  RooRealVar zjpsi_sigma("zjpsi_sigma", "zjpsi_sigma", zjpsi_sigma_value);
  RooRealVar zjpsi_alpha("zjpsi_alpha", "zjpsi_alpha", zjpsi_alpha_value);
  RooRealVar zjpsi_n("zjpsi_n", "zjpsi_n", zjpsi_n_value);
  RooCBShape zjpsi_crystal_ball ("zjpsi_crystal_ball", "zjpsi_crystal_ball", zjpsi_dimuon_mass, zjpsi_mean, zjpsi_sigma, zjpsi_alpha, zjpsi_n );

  RooRealVar zjpsi_dimuon_slope("zjpsi_dimuon_slope", "zjpsi_dimuon_slope", zjpsi_dimuon_slope_value);
  RooExponential zjpsi_dimuon_bg_exponential("zjpsi_dimuon_bg_exponential", "zjpsi_dimuon_bg_exponential", zjpsi_dimuon_mass, zjpsi_dimuon_slope);

  RooRealVar zjpsi_m_sig_tau_sig_frac("zjpsi_m_sig_tau_sig_frac", "zjpsi_m_sig_tau_sig_frac", 0.4, 0.0, 1.0);
  RooRealVar zjpsi_m_sig_tau_bg_frac("zjpsi_m_sig_tau_bg_frac", "zjpsi_m_sig_tau_bg_frac", 0.4, 0.0, 1.0);
  RooRealVar zjpsi_m_bg_tau_bg_frac("zjpsi_m_bg_tau_bg_frac", "zjpsi_m_bg_tau_bg_frac", 0.1, 0.0, 1.0);
  RooRealVar zjpsi_m_bg_tau_sig_frac("zjpsi_m_bg_tau_sig_frac", "zjpsi_m_bg_tau_sig_frac", 0.1, 0.0, 1.0);

  //RooProdPdf zjpsi_m_sig_tau_sig("zjpsi_m_sig_tau_sig", "zjpsi_m_sig_tau_sig", RooArgList(zjpsi_dimuon_gauss, zjpsi_tau_xy_gauss_sum_fitpdf ));
  //RooProdPdf zjpsi_m_sig_tau_bg("zjpsi_m_sig_tau_bg", "zjpsi_m_sig_tau_bg", RooArgList(zjpsi_dimuon_gauss, zjpsi_decay_exp ));
  RooProdPdf zjpsi_m_sig_tau_sig("zjpsi_m_sig_tau_sig", "zjpsi_m_sig_tau_sig", RooArgList(zjpsi_crystal_ball, zjpsi_tau_xy_gauss_sum_fitpdf ));
  RooProdPdf zjpsi_m_sig_tau_bg("zjpsi_m_sig_tau_bg", "zjpsi_m_sig_tau_bg", RooArgList(zjpsi_crystal_ball, zjpsi_decay_exp ));
  RooProdPdf zjpsi_m_bg_tau_bg("zjpsi_m_bg_tau_bg", "zjpsi_m_bg_tau_bg", RooArgList(zjpsi_dimuon_bg_exponential, zjpsi_decay_exp ));
  RooProdPdf zjpsi_m_bg_tau_sig("zjpsi_m_bg_tau_sig", "zjpsi_m_bg_tau_sig", RooArgList(zjpsi_dimuon_bg_exponential, zjpsi_tau_xy_gauss_sum_fitpdf ));

  RooAddPdf zjpsi_mass_tau_xy_fitpdf("zjpsi_mass_tau_xy_fitpdf","zjpsi_mass_tau_xy_fitpdf",RooArgSet(zjpsi_m_sig_tau_sig, zjpsi_m_sig_tau_bg, zjpsi_m_bg_tau_sig, zjpsi_m_bg_tau_bg ),
                                                                         RooArgList(zjpsi_m_sig_tau_sig_frac, zjpsi_m_sig_tau_bg_frac, zjpsi_m_bg_tau_sig_frac));


  RooFitResult *zjpsi_jpsi_fitres = zjpsi_mass_tau_xy_fitpdf.fitTo(zjpsi_dimuon_mass_data_hist, NumCPU(N_CPU), Verbose(false), PrintLevel(-1), SumW2Error(kFALSE), Save());
  zjpsi_jpsi_fitres->Print();

  std::cout << zjpsi_hist_name << std::endl;
  if (USE_Z_TO_EE) {
    std::cout << "z to ee + j/psi" << std::endl;
  }
  else {
    std::cout << "z to mumu + j/psi" << std::endl;
  }

  double zjpsi_prompt_events = 0.;
  double zjpsi_prompt_events_error = 0.;
  double zjpsi_nonprompt_events = 0.;
  double zjpsi_nonprompt_events_error = 0.;
  TAxis *x_axis = h_z_dimuon_mass->GetXaxis();
  TAxis *y_axis = h_z_dimuon_mass->GetYaxis();
  int x_bmin = x_axis->FindBin(dimuon_mass_min);
  int x_bmax = x_axis->FindBin(dimuon_mass_max);
  int y_bmin = y_axis->FindBin(tau_xy_min);
  int y_bmax = y_axis->FindBin(tau_xy_max);
  //TODO testing
  //double integral = h_z_dimuon_mass->Integral(x_bmin, x_bmax, y_bmin, y_bmax);
  double integral = h_z_dimuon_mass->Integral(x_bmin, x_bmax - 1, y_bmin, y_bmax);
  //integral -= h_inclusive_dimuon_mass->GetBinContent(bmin)*(tau_xy_min - axis->GetBinLowEdge(bmin)) / axis->GetBinWidth(bmin);
  //integral -= h_inclusive_dimuon_mass->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax) - tau_xy_max)  / axis->GetBinWidth(bmax); 
  //TODO testing
  std::cout << integral << std::endl;
  

  double prompt_fraction_fit_value = zjpsi_m_sig_tau_sig_frac.getVal();
  double prompt_fraction_fit_err = zjpsi_m_sig_tau_sig_frac.getError();
  double nonprompt_fraction_fit_value = zjpsi_m_sig_tau_bg_frac.getVal();
  double nonprompt_fraction_fit_err = zjpsi_m_sig_tau_bg_frac.getError();
  double integral_error = pow (integral, 0.5);
  zjpsi_prompt_events = integral * prompt_fraction_fit_value;
  // f = A*B, f = zjpsi_prompt events, A = prompt fraction, B = zjpsi_total_events
  //sigma_f = sqrt (A^2 * sigma_b^2 + B^2 * sigma_a^2 + sigma_a^2 * sigma_b^2 )
  zjpsi_prompt_events_error = pow(pow(prompt_fraction_fit_value * integral_error , 2.0) +
                                  pow(integral * prompt_fraction_fit_err, 2.0 ) +
                                  pow(integral_error * prompt_fraction_fit_err, 2.0) , 0.5);

  zjpsi_nonprompt_events = integral * nonprompt_fraction_fit_value;
  //zjpsi_non_prompt_events = integral * (1 - prompt_fraction_fit_value);
  //zjpsi_non_prompt_events_error = pow(pow((1.0 - prompt_fraction_fit_value) * integral_error , 2.0) +
  //                                pow(integral * prompt_fraction_fit_err, 2.0 ) +
  //                                pow(integral_error * prompt_fraction_fit_err, 2.0) , 0.5);
  zjpsi_nonprompt_events_error = pow(pow(nonprompt_fraction_fit_value * integral_error , 2.0) +
                                  pow(integral * nonprompt_fraction_fit_err, 2.0 ) +
                                  pow(integral_error * nonprompt_fraction_fit_err, 2.0) , 0.5);

  std::cout << "zjpsi_prompt_events: " << zjpsi_prompt_events << std::endl;
  std::cout << "zjpsi_prompt_events_error: " << zjpsi_prompt_events_error << std::endl;
  std::cout << "zjpsi_nonprompt_events: " << zjpsi_nonprompt_events << std::endl;
  std::cout << "zjpsi_nonprompt_events_error: " << zjpsi_nonprompt_events_error << std::endl;

  //TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 750);

  TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 500);

  // Plot the left side
  canvas->cd(1);
  //gPad->SetLogy();
  RooPlot* dimuon_mass_fitframe;
  dimuon_mass_fitframe = dimuon_mass.frame( Title("Dimuon" ));
  //dimuon_mass_fitframe->SetName(0); // Unset title
  dimuon_mass_data_hist.plotOn(dimuon_mass_fitframe);
  mass_tau_xy_fitpdf.plotOn(dimuon_mass_fitframe, LineColor(kRed-2), RooFit::Name("total"));
  mass_tau_xy_fitpdf.plotOn(dimuon_mass_fitframe, Components(m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  mass_tau_xy_fitpdf.plotOn(dimuon_mass_fitframe, Components(m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  mass_tau_xy_fitpdf.plotOn(dimuon_mass_fitframe, Components(m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  mass_tau_xy_fitpdf.plotOn(dimuon_mass_fitframe, Components(m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));


  dimuon_mass_fitframe->Draw();
  Double_t xl1=.58, yl1=0.55, xl2=xl1+.3, yl2=yl1+.325;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillColor(kWhite);
  leg->AddEntry(dimuon_mass_fitframe->findObject("total"),"total","l");
  leg->AddEntry(dimuon_mass_fitframe->findObject("prompt j/psi"),"prompt j/psi","l");
  leg->AddEntry(dimuon_mass_fitframe->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg->AddEntry(dimuon_mass_fitframe->findObject("prompt continuum"),"prompt continuum","l");
  leg->AddEntry(dimuon_mass_fitframe->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg->Draw();
  std::string inclusive_jpsi_mass_image_name = OUT_DIR;
  inclusive_jpsi_mass_image_name.append("inclusive_jpsi_mass");
  inclusive_jpsi_mass_image_name.append(".png");
  canvas->Print(inclusive_jpsi_mass_image_name.c_str() , "png");
  canvas->Close();

  TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 1000, 500);
  canvas2->cd(1);
  //gPad->SetLogy();
  RooPlot* tau_xy_fitframe = tau_xy.frame( Title("J/Psi Lifetime") , Range(tau_xy_min, tau_xy_max ));
  dimuon_mass_data_hist.plotOn(tau_xy_fitframe);
  mass_tau_xy_fitpdf.plotOn(tau_xy_fitframe, LineColor(kRed-2), RooFit::Name("total"));
  mass_tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  mass_tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  mass_tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  mass_tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));
  tau_xy_fitframe->Draw();

  TLegend *leg2 = new TLegend(xl1,yl1,xl2,yl2);
  leg2->SetFillColor(kWhite);
  leg2->AddEntry(dimuon_mass_fitframe->findObject("total"),"total","l");
  leg2->AddEntry(dimuon_mass_fitframe->findObject("prompt j/psi"),"prompt j/psi","l");
  leg2->AddEntry(dimuon_mass_fitframe->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg2->AddEntry(dimuon_mass_fitframe->findObject("prompt continuum"),"prompt continuum","l");
  leg2->AddEntry(dimuon_mass_fitframe->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg2->Draw();

  std::string inclusive_jpsi_tau_xy_image_name = OUT_DIR;
  inclusive_jpsi_tau_xy_image_name.append("inclusive_jpsi_tau_xy");
  inclusive_jpsi_tau_xy_image_name.append(".png");
  canvas2->Print(inclusive_jpsi_tau_xy_image_name.c_str() , "png");
  canvas2->Close();

  TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 1000, 500);
  canvas3->cd(1);
  //gPad->SetLogy();
  RooPlot* zjpsi_tau_xy_fitframe;
  if (USE_Z_TO_EE) { 
    zjpsi_tau_xy_fitframe = zjpsi_tau_xy.frame( Title("J/Psi Lifetime (Z->ee + J/Psi->mumu)") );
  }
  else {
    zjpsi_tau_xy_fitframe = zjpsi_tau_xy.frame( Title("J/Psi Lifetime (Z->mumu + J/Psi->mumu)") );
  }
  zjpsi_dimuon_mass_data_hist.plotOn(zjpsi_tau_xy_fitframe);
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, LineColor(kRed-2), RooFit::Name("total"));
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));
  zjpsi_tau_xy_fitframe->Draw();

  TLegend *leg3 = new TLegend(xl1,yl1,xl2,yl2);
  leg3->SetFillColor(kWhite);
  leg3->AddEntry(dimuon_mass_fitframe->findObject("total"),"total","l");
  leg3->AddEntry(dimuon_mass_fitframe->findObject("prompt j/psi"),"prompt j/psi","l");
  leg3->AddEntry(dimuon_mass_fitframe->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg3->AddEntry(dimuon_mass_fitframe->findObject("prompt continuum"),"prompt continuum","l");
  leg3->AddEntry(dimuon_mass_fitframe->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg3->Draw();

  std::string zjpsi_tau_xy_image_name = OUT_DIR;
  if (USE_Z_TO_EE) {
    zjpsi_tau_xy_image_name.append("ztoee_jpsi_tau_xy");
  }
  else {
    zjpsi_tau_xy_image_name.append("ztomumu_jpsi_tau_xy");
  }
  zjpsi_tau_xy_image_name.append(".png");
  canvas3->Print(zjpsi_tau_xy_image_name.c_str() , "png");
  canvas3->Close();

  TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 1000, 500);
  canvas4->cd(1);
  //gPad->SetLogy();
  RooPlot* zjpsi_dimuon_mass_fitframe;
  if (USE_Z_TO_EE) { 
    zjpsi_dimuon_mass_fitframe = zjpsi_dimuon_mass.frame( Title("Dimuon Mass (Z->ee + Dimuon)") );
  }
  else {
    zjpsi_dimuon_mass_fitframe = zjpsi_dimuon_mass.frame( Title("Dimuon Mass (Z->mumu + Dimuon)") );
  }
  zjpsi_dimuon_mass_data_hist.plotOn(zjpsi_dimuon_mass_fitframe);
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_dimuon_mass_fitframe, LineColor(kRed-2), RooFit::Name("total"));
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_dimuon_mass_fitframe, Components(zjpsi_m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_dimuon_mass_fitframe, Components(zjpsi_m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_dimuon_mass_fitframe, Components(zjpsi_m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  zjpsi_mass_tau_xy_fitpdf.plotOn(zjpsi_dimuon_mass_fitframe, Components(zjpsi_m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));
  zjpsi_dimuon_mass_fitframe->Draw();

  TLegend *leg4 = new TLegend(xl1,yl1,xl2,yl2);
  leg4->SetFillColor(kWhite);
  leg4->AddEntry(dimuon_mass_fitframe->findObject("total"),"total","l");
  leg4->AddEntry(dimuon_mass_fitframe->findObject("prompt j/psi"),"prompt j/psi","l");
  leg4->AddEntry(dimuon_mass_fitframe->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg4->AddEntry(dimuon_mass_fitframe->findObject("prompt continuum"),"prompt continuum","l");
  leg4->AddEntry(dimuon_mass_fitframe->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg4->Draw();

  std::string zjpsi_dimuon_mass_image_name = OUT_DIR;
  if (USE_Z_TO_EE) {
    zjpsi_dimuon_mass_image_name.append("ztoee_jpsi_dimuon_mass");
  }
  else {
    zjpsi_dimuon_mass_image_name.append("ztomumu_jpsi_dimuon_mass");
  }
  zjpsi_dimuon_mass_image_name.append(".png");
  canvas4->Print(zjpsi_dimuon_mass_image_name.c_str() , "png");
  canvas4->Close();

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
    bool USE_Z_TO_EE = true;
    std::istringstream ss(argv[3]);
    if (!(ss >> USE_Z_TO_EE ) ) {
      std::cout << "Invalid bool " << argv[3] << std::endl;
      return 1;
    }
    const std::string OUT_DIR(argv[4]);
    RooFitLifetimeAndMassCrystalBall(DATA_FILE_1, DATA_FILE_2, USE_Z_TO_EE, OUT_DIR);
    return 0;
  }
}
