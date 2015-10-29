#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TMath.h"
#include "TGraph.h"
#include "RooPolynomial.h"
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooErrorVar.h"
#include "RooVoigtian.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooNLLVar.h"
#include "TLatex.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TPaveLabel.h"
#include "RooNLLVar.h"
#include "RooProfileLL.h"
#include "RooProdPdf.h"
#include "RooSimultaneous.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooProduct.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "TCut.h"
#include "TTree.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TBranch.h"
#include "TFolder.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "TStyle.h"
#include "RooCBShape.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooHistPdf.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SPlot.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooChebychev.h"
#include "RooDLLSignificanceMCSModule.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/SPlot.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/NumberCountingPdfFactory.h"

using namespace RooFit ;
using namespace RooStats ;

const float PDGMass = 3.096916;

// Enables residuals on fits
//#define DrawResiduals

//void sPlotFit();
//void sPlotFit() { sPlotFit(); }
void sPlotFit() {

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TF1 *f_straighline = new TF1("f_straighline", "0", 0, 1000);
  f_straighline->SetLineColor(kBlack);

  RooRealVar *Oniamass = new RooRealVar("Oniamass", "Oniamass", 2.85, 3.35);
  //RooRealVar *Oniatau  = new RooRealVar("Oniatau", "Oniatau", -5, 10); //this range may be better
  RooRealVar *Oniatau  = new RooRealVar("Oniatau", "Oniatau", -0.3, 5);
  RooRealVar *ZMass    = new RooRealVar("ZMass", "ZMass", 60., 120.);
  RooRealVar *isZee    = new RooRealVar("isZee", "isZee", 0, 1);
  RooRealVar *isZmumu  = new RooRealVar("isZmumu", "isZmumu", 0, 1);

  RooRealVar *inclusive_jpsi_Oniamass = new RooRealVar("inclusive_jpsi_Oniamass", "inclusive_jpsi_Oniamass", 2.85, 3.35);
  //RooRealVar *inclusive_jpsi_Oniatau  = new RooRealVar("inclusive_jpsi_Oniatau", "inclusive_jpsi_Oniatau", -5, 10); //this range may be better
  RooRealVar *inclusive_jpsi_Oniatau  = new RooRealVar("inclusive_jpsi_Oniatau", "inclusive_jpsi_Oniatau", -0.3, 5);
  RooRealVar *inclusive_jpsi_ZMass    = new RooRealVar("inclusive_jpsi_ZMass", "inclusive_jpsi_ZMass", 60., 120.);
  RooRealVar *inclusive_jpsi_isZee    = new RooRealVar("inclusive_jpsi_isZee", "inclusive_jpsi_isZee", 0, 1);
  RooRealVar *inclusive_jpsi_isZmumu  = new RooRealVar("inclusive_jpsi_isZmumu", "inclusive_jpsi_isZmumu", 0, 1);

  RooArgSet Variables(*Oniamass, *Oniatau, *isZee, *isZmumu, *ZMass);
  RooArgSet Inclusive_Jpsi_ArgSet(*inclusive_jpsi_Oniamass, *inclusive_jpsi_Oniatau, *inclusive_jpsi_isZee, *inclusive_jpsi_isZmumu, *inclusive_jpsi_ZMass);

  //Data Set
  //TFile *ntuple = new TFile("../SkimNtuple/ZJpsi.root");

  TFile *ntuple = new TFile("ZJpsi.root");
  TFile *inclusive_jpsi_ntuple = new TFile("Jpsimumu.root");

  TTree* inclusive_jpsi_tree = (TTree*) inclusive_jpsi_ntuple->Get("AUX");
  RooDataSet *inclusive_jpsi_data = new RooDataSet("inclusive_jpsi_data", "inclusive_jpsi_data", inclusive_jpsi_tree, Inclusive_Jpsi_ArgSet);

  TTree* tree = (TTree*) ntuple->Get("AUX");
  RooDataSet *data = new RooDataSet("data", "data", tree, Variables);

  //      TCut SelectionCut = "Oniamass>2.6 && Oniamass<3.6";
  //  RooDataSet *newdata = (RooDataSet*)data->reduce(SelectionCut);

  // *****************************************************************
  //  Jpsi 
  // *****************************************************************

  RooRealVar jpsi_gaussmean("jpsi_gaussmean", "Mean of the smearing jpsi Gaussian", 0.0);
  RooRealVar jpsi_gausssigma("jpsi_gausssigma", "Width of the smearing jpsi Gaussian", 0.01, 0.005, 0.3);
  RooGaussModel jpsi_smear_gauss_model("jpsi_smear_gauss", "Gaussian used to smear the jpsi Exponential Decay", *inclusive_jpsi_Oniatau, jpsi_gaussmean, jpsi_gausssigma);
  RooRealVar jpsi_decay_lifetime("jpsi_decay_lifetime", "jpsi_Tau", 1.3, 0.05, 2.);
  RooDecay jpsi_decay_exp("jpsi_decay_exp", "jpsi_Exponential Decay", *inclusive_jpsi_Oniatau,jpsi_decay_lifetime,jpsi_smear_gauss_model,RooDecay::SingleSided);

  RooRealVar jpsi_gauss_prompt_mean("jpsi_gauss_prompt_mean", "Mean of the Prompt jpsi_Gaussian", 0., -0.1, 0.1);
  RooRealVar jpsi_gauss_prompt_sigma("jpsi_gauss_prompt_sigma", "Width of the Prompt jpsi_Gaussian", 0.007, 0.005, 0.25);
  RooGaussian jpsi_prompt_gauss("jpsi_prompt_gauss", "Gaussian of the Prompt jpsi_Peak", *inclusive_jpsi_Oniatau, jpsi_gauss_prompt_mean, jpsi_gauss_prompt_sigma);

  RooRealVar jpsi_gauss_prompt_sigma_2("jpsi_gauss_prompt_sigma_2", "Width of the Prompt jpsi_Gaussian_2", 0.01, 0.008, 0.4);
  RooGaussian jpsi_prompt_gauss_2("jpsi_prompt_gauss_2", "Gaussian_2 of the Prompt jpsi_Peak", *inclusive_jpsi_Oniatau, jpsi_gauss_prompt_mean, jpsi_gauss_prompt_sigma_2);

  RooRealVar jpsi_frac_prompt_sharp("jpsi_frac_prompt_sharp", "jpsi_frac_prompt_sharp" , 0.5 , 0.0, 1.);
  RooAddPdf jpsi_oniatau_gauss_sum_fitpdf("jpsi_oniatau_gauss_sum_fitpdf", "jpsi_oniatau_gauss_sum_fitpdf", RooArgList(jpsi_prompt_gauss,jpsi_prompt_gauss_2), RooArgList(jpsi_frac_prompt_sharp));

  RooRealVar jpsi_prompt_fraction("jpsi_prompt_fraction", "jpsi_prompt_fraction" , 0.01 , 0.0, 1);

  RooRealVar jpsi_sigma("jpsi_sigma", "jpsi_sigma", 0.05, 0.001, 0.1);
  RooRealVar jpsi_alpha("jpsi_alpha", "jpsi_alpha", 1.8, 1.0, 2.5);
  RooRealVar jpsi_n("jpsi_n", "jpsi_n", 2., 1.0, 80.);
  RooRealVar jpsi_dimuon_mean("jpsi_dimuon_mean", "jpsi_dimuon_mean", 3.1, 3.0, 3.2);
  RooCBShape jpsi_crystal_ball ("jpsi_crystal_ball", "jpsi_crystal_ball", *inclusive_jpsi_Oniamass, jpsi_dimuon_mean, jpsi_sigma, jpsi_alpha, jpsi_n );
  RooRealVar jpsi_dimuon_sigma("jpsi_dimuon_sigma", "jpsi_dimuon_sigma", 0.02, 0.001, 0.1);
  RooGaussian jpsi_dimuon_gauss ("jpsi_dimuon_gauss", "jpsi_dimuon_gauss", *inclusive_jpsi_Oniamass, jpsi_dimuon_mean, jpsi_dimuon_sigma);
  RooRealVar jpsi_frac_cball("jpsi_frac_cball", "jpsi_frac_cball" , 0.5 , 0.0, 1.);
  RooAddPdf jpsi_mass_signal("jpsi_mass_signal", "jpsi_mass_signal", RooArgList(jpsi_crystal_ball, jpsi_dimuon_gauss), RooArgList(jpsi_frac_cball));

  RooRealVar jpsi_dimuon_slope("jpsi_dimuon_slope", "jpsi_dimuon_slope", -0.1, -10., 10.);
  RooExponential jpsi_dimuon_bg_exponential("jpsi_dimuon_bg_exponential", "jpsi_dimuon_bg_exponential", *inclusive_jpsi_Oniamass, jpsi_dimuon_slope);

  RooRealVar jpsi_np_dimuon_slope("jpsi_np_dimuon_slope", "jpsi_np_dimuon_slope", -0.1, -10., 10.);
  RooExponential jpsi_np_dimuon_bg_exponential("jpsi_np_dimuon_bg_exponential", "jpsi_np_dimuon_bg_exponential", *inclusive_jpsi_Oniamass, jpsi_np_dimuon_slope);

  // PDFs
  RooProdPdf jpsi_m_sig_tau_sig("jpsi_m_sig_tau_sig", "jpsi_m_sig_tau_sig", RooArgList(jpsi_mass_signal, jpsi_oniatau_gauss_sum_fitpdf ));
  RooProdPdf jpsi_m_sig_tau_bg ("jpsi_m_sig_tau_bg",  "jpsi_m_sig_tau_bg",  RooArgList(jpsi_mass_signal, jpsi_decay_exp ));
  RooProdPdf jpsi_m_bg_tau_bg  ("jpsi_m_bg_tau_bg",   "jpsi_m_bg_tau_bg",   RooArgList(jpsi_np_dimuon_bg_exponential, jpsi_decay_exp ));
  RooProdPdf jpsi_m_bg_tau_sig ("jpsi_m_bg_tau_sig",  "jpsi_m_bg_tau_sig",  RooArgList(jpsi_dimuon_bg_exponential, jpsi_oniatau_gauss_sum_fitpdf ));

  //yield ranges probably need to be changed?
  // yields
  RooRealVar Njpsi_m_sig_tau_sig("Njpsi_m_sig_tau_sig", "Njpsi_m_sig_tau_sig",10, 0, 100);  
  RooRealVar Njpsi_m_sig_tau_bg ("Njpsi_m_sig_tau_bg",  "Njpsi_m_sig_tau_bg", 10, 0, 100); 
  RooRealVar Njpsi_m_bg_tau_bg  ("Njpsi_m_bg_tau_bg",   "Njpsi_m_bg_tau_bg",  10, 0, 100); 
  RooRealVar Njpsi_m_bg_tau_sig ("Njpsi_m_bg_tau_sig",  "Njpsi_m_bg_tau_sig", 10, 0, 100); 

  // extended PDFs
  RooExtendPdf ejpsi_m_sig_tau_sig("ejpsi_m_sig_tau_sig", "ejpsi_m_sig_tau_sig", jpsi_m_sig_tau_sig, Njpsi_m_sig_tau_sig); 
  RooExtendPdf ejpsi_m_sig_tau_bg ("ejpsi_m_sig_tau_bg",  "ejpsi_m_sig_tau_bg",  jpsi_m_sig_tau_bg , Njpsi_m_sig_tau_bg );
  RooExtendPdf ejpsi_m_bg_tau_bg  ("ejpsi_m_bg_tau_bg",   "ejpsi_m_bg_tau_bg",   jpsi_m_bg_tau_bg  , Njpsi_m_bg_tau_bg  );
  RooExtendPdf ejpsi_m_bg_tau_sig ("ejpsi_m_bg_tau_sig",  "ejpsi_m_bg_tau_sig",  jpsi_m_bg_tau_sig , Njpsi_m_bg_tau_sig );

  //  RooAddPdf model("model", "model", RooArgSet(jpsi_m_sig_tau_sig, jpsi_m_sig_tau_bg, jpsi_m_bg_tau_sig, jpsi_m_bg_tau_bg ), RooArgList(jpsi_m_sig_tau_sig_frac, jpsi_m_sig_tau_bg_frac, jpsi_m_bg_tau_sig_frac), kTRUE);
  RooAddPdf inclusive_jpsi_model("inclusive_jpsi_model", "inclusive_jpsi_model", RooArgList(ejpsi_m_sig_tau_sig, ejpsi_m_sig_tau_bg, ejpsi_m_bg_tau_sig, ejpsi_m_bg_tau_bg));

  RooFitResult *inclusive_jpsi_fr = inclusive_jpsi_model.fitTo(*inclusive_jpsi_data, NumCPU(4, kTRUE), Save());

  //Then would do:

  //double zjpsi_gaussmean_value = jpsi_gaussmean.getVal();
  //double zjpsi_gausssigma_value     = jpsi_gausssigma.getVal();
  //double zjpsi_gausssigma_value_err = jpsi_gausssigma.getError();
  //double zjpsi_decay_lifetime_value     = jpsi_decay_lifetime.getVal();
  //double zjpsi_decay_lifetime_value_err = jpsi_decay_lifetime.getError();
  //double zjpsi_gauss_prompt_mean_value     = jpsi_gauss_prompt_mean.getVal();
  //double zjpsi_gauss_prompt_mean_value_err = jpsi_gauss_prompt_mean.getError();
  //double zjpsi_gauss_prompt_sigma_value     = jpsi_gauss_prompt_sigma.getVal();
  //double zjpsi_gauss_prompt_sigma_value_err = jpsi_gauss_prompt_sigma.getError();
  //double zjpsi_gauss_prompt_sigma_2_value     = jpsi_gauss_prompt_sigma_2.getVal();
  //double zjpsi_gauss_prompt_sigma_2_value_err = jpsi_gauss_prompt_sigma_2.getError();
  //double zjpsi_frac_prompt_sharp_value     = jpsi_frac_prompt_sharp.getVal();
  //double zjpsi_frac_prompt_sharp_value_err = jpsi_frac_prompt_sharp.getError();

  //double zjpsi_sigma_value     = jpsi_sigma.getVal();
  //double zjpsi_sigma_value_err = jpsi_sigma.getError();
  //double zjpsi_alpha_value     = jpsi_alpha.getVal();
  //double zjpsi_alpha_value_err = jpsi_alpha.getError();
  //double zjpsi_n_value     = jpsi_n.getVal();
  //double zjpsi_n_value_err = jpsi_n.getError();
  //double zjpsi_dimuon_slope_value     = jpsi_dimuon_slope.getVal();
  //double zjpsi_dimuon_slope_value_err = jpsi_dimuon_slope.getError();
  //double zjpsi_np_dimuon_slope_value     = jpsi_np_dimuon.getVal();
  //double zjpsi_np_dimuon_slope_value_err = jpsi_np_dimuon_slope.getError();
  //double zjpsi_dimuon_mean_value     = jpsi_dimuon_mean.getVal();
  //double zjpsi_dimuon_mean_value_err = jpsi_dimuon_mean.getError();
  //double zjpsi_dimuon_sigma_value     = jpsi_dimuon_sigma.getVal();
  //double zjpsi_dimuon_sigma_value_err = jpsi_dimuon_sigma.getError();
  //double zjpsi_frac_cball_value     = jpsi_frac_cball.getVal();
  //double zjpsi_frac_cball_value_err = jpsi_frac_cball.getError();

  // ZJPsi Unbinned fit start here

  double zjpsi_gaussmean_value = 0.; //gaussmean.getVal();
  double zjpsi_gausssigma_value     = 1.3496e-01;
  double zjpsi_gausssigma_value_err = 2.41e-03;
  double zjpsi_decay_lifetime_value     = 1.0239;
  double zjpsi_decay_lifetime_value_err = 3.13e-03;
  double zjpsi_gauss_prompt_mean_value     = 1.0000e-01;
  double zjpsi_gauss_prompt_mean_value_err = 1.83e-04;
  double zjpsi_gauss_prompt_sigma_value     = 2.3806e-01;
  double zjpsi_gauss_prompt_sigma_value_err = 1.88e-01;
  double zjpsi_gauss_prompt_sigma_2_value     = 3.8616e-01;
  double zjpsi_gauss_prompt_sigma_2_value_err = 7.18e-03;
  double zjpsi_frac_prompt_sharp_value     = 1.6311e-01;
  double zjpsi_frac_prompt_sharp_value_err = 1.11e-01;

  RooRealVar zjpsi_gaussmean("zjpsi_gaussmean", "Mean of the smearing zjpsi Gaussian", zjpsi_gaussmean_value);
  RooRealVar zjpsi_gausssigma("zjpsi_gausssigma", "Width of the smearing zjpsi Gaussian", zjpsi_gausssigma_value);
  RooGaussModel zjpsi_smear_gauss_model("zjpsi_smear_gauss", "Gaussian used to smear the zjpsi Exponential Decay", *Oniatau, zjpsi_gaussmean, zjpsi_gausssigma);
  RooRealVar zjpsi_decay_lifetime("zjpsi_decay_lifetime", "zjpsi_Tau", zjpsi_decay_lifetime_value, zjpsi_decay_lifetime_value-zjpsi_decay_lifetime_value_err, zjpsi_decay_lifetime_value+zjpsi_decay_lifetime_value_err);
  RooDecay zjpsi_decay_exp("zjpsi_decay_exp", "zjpsi_Exponential Decay", *Oniatau,zjpsi_decay_lifetime,zjpsi_smear_gauss_model,RooDecay::SingleSided);

  RooRealVar zjpsi_gauss_prompt_mean("zjpsi_gauss_prompt_mean", "Mean of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_mean_value, zjpsi_gauss_prompt_mean_value-zjpsi_gauss_prompt_mean_value_err, zjpsi_gauss_prompt_mean_value+zjpsi_gauss_prompt_mean_value_err);
  RooRealVar zjpsi_gauss_prompt_sigma("zjpsi_gauss_prompt_sigma", "Width of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_sigma_value, zjpsi_gauss_prompt_sigma_value-zjpsi_gauss_prompt_sigma_value_err, zjpsi_gauss_prompt_sigma_value+zjpsi_gauss_prompt_sigma_value_err);
  RooGaussian zjpsi_prompt_gauss("zjpsi_prompt_gauss", "Gaussian of the Prompt zjpsi_Peak", *Oniatau, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma);

  RooRealVar zjpsi_gauss_prompt_sigma_2("zjpsi_gauss_prompt_sigma_2", "Width of the Prompt zjpsi_Gaussian_2", zjpsi_gauss_prompt_sigma_2_value, zjpsi_gauss_prompt_sigma_2_value-zjpsi_gauss_prompt_sigma_2_value_err, zjpsi_gauss_prompt_sigma_2_value+zjpsi_gauss_prompt_sigma_2_value_err);
  RooGaussian zjpsi_prompt_gauss_2("zjpsi_prompt_gauss_2", "Gaussian_2 of the Prompt zjpsi_Peak", *Oniatau, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma_2);

  RooRealVar zjpsi_frac_prompt_sharp("zjpsi_frac_prompt_sharp", "zjpsi_frac_prompt_sharp" , zjpsi_frac_prompt_sharp_value, zjpsi_frac_prompt_sharp_value-zjpsi_frac_prompt_sharp_value_err, zjpsi_frac_prompt_sharp_value+zjpsi_frac_prompt_sharp_value_err);
  RooAddPdf Oniatau_gauss_sum_fitpdf("Oniatau_gauss_sum_fitpdf", "Oniatau_gauss_sum_fitpdf", RooArgList(zjpsi_prompt_gauss,zjpsi_prompt_gauss_2), RooArgList(zjpsi_frac_prompt_sharp));

  RooRealVar zjpsi_prompt_fraction("zjpsi_prompt_fraction", "zjpsi_prompt_fraction" , 0.01 , 0.0, 1);

  double zjpsi_sigma_value     = 9.1932e-02;
  double zjpsi_sigma_value_err = 2.27e-03;
  double zjpsi_alpha_value     = 1.9570;
  double zjpsi_alpha_value_err = 1.34e-01 ;
  double zjpsi_n_value     = 3.4482e+00;
  double zjpsi_n_value_err = 1.60e+00;
  double zjpsi_dimuon_slope_value     = -2.7272 ;
  double zjpsi_dimuon_slope_value_err = 4.72;
  double zjpsi_np_dimuon_slope_value     = -1.8674;
  double zjpsi_np_dimuon_slope_value_err = 1.88e-01;
  double zjpsi_dimuon_mean_value     = 3.0903;
  double zjpsi_dimuon_mean_value_err = 1.31e-04;
  double zjpsi_dimuon_sigma_value     = 3.0664e-02;
  double zjpsi_dimuon_sigma_value_err = 4.46e-04;
  double zjpsi_frac_cball_value     = 1.4519e-01;
  double zjpsi_frac_cball_value_err = 2.62e-02;

  RooRealVar zjpsi_sigma("zjpsi_sigma", "zjpsi_sigma", zjpsi_sigma_value, zjpsi_sigma_value-zjpsi_sigma_value_err, zjpsi_sigma_value+zjpsi_sigma_value_err);
  RooRealVar zjpsi_alpha("zjpsi_alpha", "zjpsi_alpha", zjpsi_alpha_value, zjpsi_alpha_value-zjpsi_alpha_value_err, zjpsi_alpha_value+zjpsi_alpha_value_err);
  RooRealVar zjpsi_n("zjpsi_n", "zjpsi_n", zjpsi_n_value, zjpsi_n_value-zjpsi_n_value_err, zjpsi_n_value+zjpsi_n_value_err);
  RooRealVar zjpsi_dimuon_mean("zjpsi_dimuon_mean", "zjpsi_dimuon_mean", zjpsi_dimuon_mean_value, zjpsi_dimuon_mean_value-zjpsi_dimuon_mean_value_err, zjpsi_dimuon_mean_value+zjpsi_dimuon_mean_value_err);
  RooCBShape zjpsi_crystal_ball ("zjpsi_crystal_ball", "zjpsi_crystal_ball", *Oniamass, zjpsi_dimuon_mean, zjpsi_sigma, zjpsi_alpha, zjpsi_n );
  RooRealVar zjpsi_dimuon_sigma("zjpsi_dimuon_sigma", "zjpsi_dimuon_sigma", zjpsi_dimuon_sigma_value, zjpsi_dimuon_sigma_value-zjpsi_dimuon_sigma_value_err, zjpsi_dimuon_sigma_value+zjpsi_dimuon_sigma_value_err);
  RooGaussian zjpsi_dimuon_gauss ("zjpsi_dimuon_gauss", "zjpsi_dimuon_gauss", *Oniamass, zjpsi_dimuon_mean, zjpsi_dimuon_sigma);
  RooRealVar zjpsi_frac_cball("zjpsi_frac_cball", "zjpsi_frac_cball" , zjpsi_frac_cball_value, zjpsi_frac_cball_value-zjpsi_frac_cball_value_err, zjpsi_frac_cball_value+zjpsi_frac_cball_value_err);
  RooAddPdf zjpsi_mass_signal("zjpsi_mass_signal", "zjpsi_mass_signal", RooArgList(zjpsi_crystal_ball, zjpsi_dimuon_gauss), RooArgList(zjpsi_frac_cball));

  RooRealVar zjpsi_dimuon_slope("zjpsi_dimuon_slope", "zjpsi_dimuon_slope", zjpsi_dimuon_slope_value, zjpsi_dimuon_slope_value-zjpsi_dimuon_slope_value_err, zjpsi_dimuon_slope_value+zjpsi_dimuon_slope_value_err);
  RooExponential zjpsi_dimuon_bg_exponential("zjpsi_dimuon_bg_exponential", "zjpsi_dimuon_bg_exponential", *Oniamass, zjpsi_dimuon_slope);

  RooRealVar zjpsi_np_dimuon_slope("zjpsi_np_dimuon_slope", "zjpsi_np_dimuon_slope", zjpsi_np_dimuon_slope_value, zjpsi_np_dimuon_slope_value-zjpsi_np_dimuon_slope_value_err, zjpsi_np_dimuon_slope_value+zjpsi_np_dimuon_slope_value_err);
  RooExponential zjpsi_np_dimuon_bg_exponential("zjpsi_np_dimuon_bg_exponential", "zjpsi_np_dimuon_bg_exponential", *Oniamass, zjpsi_np_dimuon_slope);

  //  RooRealVar zjpsi_m_sig_tau_sig_frac("zjpsi_m_sig_tau_sig_frac", "zjpsi_m_sig_tau_sig_frac", 0.2, 0.0, 1.0);
  //  RooRealVar zjpsi_m_sig_tau_bg_frac("zjpsi_m_sig_tau_bg_frac", "zjpsi_m_sig_tau_bg_frac", 0.6, 0.0, 1.0);
  //  RooRealVar zjpsi_m_bg_tau_sig_frac("zjpsi_m_bg_tau_sig_frac", "zjpsi_m_bg_tau_sig_frac", 0.15, 0.0, 1.0);

  // PDFs
  RooProdPdf zjpsi_m_sig_tau_sig("zjpsi_m_sig_tau_sig", "zjpsi_m_sig_tau_sig", RooArgList(zjpsi_mass_signal, Oniatau_gauss_sum_fitpdf ));
  RooProdPdf zjpsi_m_sig_tau_bg ("zjpsi_m_sig_tau_bg",  "zjpsi_m_sig_tau_bg",  RooArgList(zjpsi_mass_signal, zjpsi_decay_exp ));
  RooProdPdf zjpsi_m_bg_tau_bg  ("zjpsi_m_bg_tau_bg",   "zjpsi_m_bg_tau_bg",   RooArgList(zjpsi_np_dimuon_bg_exponential, zjpsi_decay_exp ));
  RooProdPdf zjpsi_m_bg_tau_sig ("zjpsi_m_bg_tau_sig",  "zjpsi_m_bg_tau_sig",  RooArgList(zjpsi_dimuon_bg_exponential, Oniatau_gauss_sum_fitpdf ));

  // yields
  RooRealVar Nzjpsi_m_sig_tau_sig("Nzjpsi_m_sig_tau_sig", "Nzjpsi_m_sig_tau_sig",10, 0, 100);  
  RooRealVar Nzjpsi_m_sig_tau_bg ("Nzjpsi_m_sig_tau_bg",  "Nzjpsi_m_sig_tau_bg", 10, 0, 100); 
  RooRealVar Nzjpsi_m_bg_tau_bg  ("Nzjpsi_m_bg_tau_bg",   "Nzjpsi_m_bg_tau_bg",  10, 0, 100); 
  RooRealVar Nzjpsi_m_bg_tau_sig ("Nzjpsi_m_bg_tau_sig",  "Nzjpsi_m_bg_tau_sig", 10, 0, 100); 

  // extended PDFs
  RooExtendPdf ezjpsi_m_sig_tau_sig("ezjpsi_m_sig_tau_sig", "ezjpsi_m_sig_tau_sig", zjpsi_m_sig_tau_sig, Nzjpsi_m_sig_tau_sig); 
  RooExtendPdf ezjpsi_m_sig_tau_bg ("ezjpsi_m_sig_tau_bg",  "ezjpsi_m_sig_tau_bg",  zjpsi_m_sig_tau_bg , Nzjpsi_m_sig_tau_bg );
  RooExtendPdf ezjpsi_m_bg_tau_bg  ("ezjpsi_m_bg_tau_bg",   "ezjpsi_m_bg_tau_bg",   zjpsi_m_bg_tau_bg  , Nzjpsi_m_bg_tau_bg  );
  RooExtendPdf ezjpsi_m_bg_tau_sig ("ezjpsi_m_bg_tau_sig",  "ezjpsi_m_bg_tau_sig",  zjpsi_m_bg_tau_sig , Nzjpsi_m_bg_tau_sig );

  //  RooAddPdf model("model", "model", RooArgSet(zjpsi_m_sig_tau_sig, zjpsi_m_sig_tau_bg, zjpsi_m_bg_tau_sig, zjpsi_m_bg_tau_bg ), RooArgList(zjpsi_m_sig_tau_sig_frac, zjpsi_m_sig_tau_bg_frac, zjpsi_m_bg_tau_sig_frac), kTRUE);
  RooAddPdf model("model", "model", RooArgList(ezjpsi_m_sig_tau_sig, ezjpsi_m_sig_tau_bg, ezjpsi_m_bg_tau_sig, ezjpsi_m_bg_tau_bg));

  RooFitResult *fr = model.fitTo(*data, NumCPU(4, kTRUE), Save());

  // *****************************
  // START OF JPSI
  // *****************************

  //----------- Mass Fit ---------------
  //Plots data on to frame
  RooPlot* mFrame = Oniamass->frame(Bins(40));
  data->plotOn(mFrame, Name("PlotData"));

  //Plots full model, prompt and non-prompt models to frame
#ifdef VisualizeError
  model.plotOn(mFrame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("masserror"), VisualizeError(*fr,3));
#endif
  model.plotOn(mFrame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("PlotModel"));
  model.plotOn(mFrame, Components("zjpsi_m_sig_tau_sig,zjpsi_m_bg_tau_sig"), LineWidth(3),LineColor(kAzure),LineStyle(2),Name("MassPrompt"));
  model.plotOn(mFrame, Components("zjpsi_m_bg_tau_sig"), LineWidth(3),LineColor(kAzure),LineStyle(2));
  model.plotOn(mFrame, Components("zjpsi_m_sig_tau_sig"), LineWidth(3),LineColor(kOrange));
  model.plotOn(mFrame, Components("zjpsi_m_sig_tau_bg,zjpsi_m_bg_tau_bg"), LineWidth(3),LineColor(kViolet),LineStyle(2), Name("MassNonPrompt"));
  model.plotOn(mFrame, Name("PlotData"));
  mFrame->GetXaxis()->SetTitleOffset(.9);

  RooCurve * mPJpsi      = mFrame->getCurve("MassPrompt");
  RooCurve * mNPJpsi   = mFrame->getCurve("MassNonPrompt");

  RooPlot* dummy_frame_jpsi = Oniamass->frame(Title("dummy frame to extract residuals"), Bins(40));
  data->plotOn(dummy_frame_jpsi); 
  model.plotOn(dummy_frame_jpsi);

  RooHist* h_residuals_mass_jpsi = dummy_frame_jpsi->pullHist();
  RooPlot* frame_residuals_mass_jpsi = Oniamass->frame(Title("Residual Distribution #mu^{+}#mu^{-} mass"));
  frame_residuals_mass_jpsi->GetYaxis()->SetTitle("(fit - data)/#sigma");
  frame_residuals_mass_jpsi->GetYaxis()->SetTitleSize(.16);
  frame_residuals_mass_jpsi->GetYaxis()->SetTitleOffset(.2);
  frame_residuals_mass_jpsi->addPlotable(h_residuals_mass_jpsi, "P");

  //comment out plots for now

  // Creates and fills canvas with plot and info 
  TCanvas *canvas1 = new TCanvas("MassFit_", "MassFit_", 900, 900); canvas1->cd();
#ifdef DrawResiduals
  TPad *pad1_jpsi = new TPad("pad1_jpsi", "The pad 80% of the height",0.0,0.05,1.0,1.0,21);
  TPad *pad2_jpsi = new TPad("pad2_jpsi", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);
  pad1_jpsi->Draw(); pad2_jpsi->Draw();
  pad1_jpsi->SetFillColor(0); pad2_jpsi->SetFillColor(0);
  pad2_jpsi->cd();
  frame_residuals_mass_jpsi->Draw(); f_straighline->Draw("same");
  pad1_jpsi->cd();
#endif
  mFrame->Draw();
  mFrame->SetTitle("J/#psi Mass [GeV]");

  //----------- Tau Fit ---------------
  //Plots data on to frame
  int timeBins = 20;
  RooPlot* tFrame = Oniatau->frame(Bins(timeBins));
  data->plotOn(tFrame, Name("PlotData"));

  //Plots full model, prompt and non-prompt models to frame


#ifdef VisualizeError
  model.plotOn(tFrame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("Timerror"), VisualizeError(*fr,3));
#endif
  model.plotOn(tFrame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("PlotModelTime"));
  model.plotOn(tFrame, Components("zjpsi_m_sig_tau_sig,zjpsi_m_bg_tau_sig"), LineWidth(3),LineColor(kAzure),LineStyle(2),Name("timePrompt"));
  model.plotOn(tFrame, Components("zjpsi_m_bg_tau_sig"), LineWidth(3),LineColor(kAzure),LineStyle(2));
  model.plotOn(tFrame, Components("zjpsi_m_sig_tau_sig"), LineWidth(3),LineColor(kOrange));
  model.plotOn(tFrame, Components("zjpsi_m_sig_tau_bg,zjpsi_m_bg_tau_bg"), LineWidth(3),LineColor(kViolet),LineStyle(2), Name("timeNonPrompt"));
  model.plotOn(tFrame, Name("PlotData"));
  tFrame->GetXaxis()->SetTitleOffset(.9);

  RooCurve * tPJpsi      = tFrame->getCurve("timePrompt");
  RooCurve * tNPJpsi   = tFrame->getCurve("timeNonPrompt");

  RooPlot* dummy_tframe_jpsi = Oniatau->frame(Title("dummy frame to extract residuals"), Bins(40));
  data->plotOn(dummy_tframe_jpsi); 
  model.plotOn(dummy_tframe_jpsi);

  RooHist* h_residuals_time_jpsi = dummy_tframe_jpsi->pullHist();
  RooPlot* frame_residuals_time_jpsi = Oniatau->frame(Title("Residual Distribution #mu^{+}#mu^{-} time"));
  frame_residuals_time_jpsi->GetYaxis()->SetTitle("(fit - data)/#sigma");
  frame_residuals_time_jpsi->GetYaxis()->SetTitleSize(.16);
  frame_residuals_time_jpsi->GetYaxis()->SetTitleOffset(.2);
  frame_residuals_time_jpsi->addPlotable(h_residuals_time_jpsi, "P");

  // Creates and fills canvas with plot and info 
  TCanvas *canvas2 = new TCanvas("TimeFit_", "TimeFit_", 900, 900); canvas2->cd();
#ifdef DrawResiduals
  TPad *tpad1_jpsi = new TPad("tpad1_jpsi", "The pad 80% of the height",0.0,0.05,1.0,1.0,21);
  TPad *tpad2_jpsi = new TPad("tpad2_jpsi", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);
  tpad1_jpsi->Draw(); tpad2_jpsi->Draw();
  tpad1_jpsi->SetFillColor(0); tpad2_jpsi->SetFillColor(0);
  tpad2_jpsi->cd();
  frame_residuals_time_jpsi->Draw(); f_straighline->Draw("same");
  tpad1_jpsi->cd();
#endif
  tFrame->Draw();
  tFrame->SetTitle("J/#psi Mass [GeV]");
  canvas2->SetLogy(1);




  //// *******************************************************
  //// *******************************************************
  //// Splot
  //// *******************************************************
  //// *******************************************************

  RooWorkspace* ws = new RooWorkspace("myWS");
  ws->import(model);
  ws->import(*data, Rename("sPlotdata"));

  // Jared: are the below correct?

  RooAbsPdf* sPlotmodel = ws->pdf("model");
  RooRealVar* SignalPromptYield = ws->var("Nzjpsi_m_sig_tau_sig");
  RooRealVar* SignalNonPromptYield = ws->var("Nzjpsi_m_sig_tau_bg");
  RooRealVar* BckgPromptYield = ws->var("Nzjpsi_m_bg_tau_sig");
  RooRealVar* BckgNonPromptYield = ws->var("Nzjpsi_m_bg_tau_bg");

  SignalPromptYield->setConstant();
  SignalNonPromptYield->setConstant();
  BckgPromptYield->setConstant();
  BckgNonPromptYield->setConstant();

  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, sPlotmodel, RooArgList(*SignalPromptYield, *SignalNonPromptYield, *BckgPromptYield, *BckgNonPromptYield));

  RooDataSet * data_prompt_weighted     = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "Nzjpsi_m_sig_tau_sig_sw");
  RooDataSet * data_non_prompt_weighted = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "Nzjpsi_m_sig_tau_bg_sw");

  TCut SelecteZee   = "isZee==1";
  TCut SelecteZmumu = "isZmumu==1";

  RooDataSet * data_prompt_Z_ee_weighted       = (RooDataSet*)data_prompt_weighted->reduce(SelecteZee);
  RooDataSet * data_non_prompt_Z_ee_weighted   = (RooDataSet*)data_non_prompt_weighted->reduce(SelecteZee);
  RooDataSet * data_prompt_Z_mumu_weighted     = (RooDataSet*)data_prompt_weighted->reduce(SelecteZmumu);
  RooDataSet * data_non_prompt_Z_mumu_weighted = (RooDataSet*)data_non_prompt_weighted->reduce(SelecteZmumu);

  //comment out drawings to run script faster for now
  int ZBins = 10;
  TCanvas* cdata_Z = new TCanvas("cdata_Z","sPlot Z", 1200, 1200); cdata_Z->Divide(2,2);
  RooPlot* frame_Z_ee_prompt = ZMass->frame(Bins(ZBins)); data_prompt_Z_ee_weighted->plotOn(frame_Z_ee_prompt);
  cdata_Z->cd(1); frame_Z_ee_prompt->Draw();

  RooPlot* frame_Z_ee_non_prompt = ZMass->frame(Bins(ZBins)); data_non_prompt_Z_ee_weighted->plotOn(frame_Z_ee_non_prompt);
  cdata_Z->cd(2); frame_Z_ee_non_prompt->Draw();

  RooPlot* frame_Z_mumu_prompt = ZMass->frame(Bins(ZBins)); data_prompt_Z_mumu_weighted->plotOn(frame_Z_mumu_prompt);
  cdata_Z->cd(3); frame_Z_mumu_prompt->Draw();

  RooPlot* frame_Z_mumu_non_prompt = ZMass->frame(Bins(ZBins)); data_non_prompt_Z_mumu_weighted->plotOn(frame_Z_mumu_non_prompt);
  cdata_Z->cd(4); frame_Z_mumu_non_prompt->Draw();
}
