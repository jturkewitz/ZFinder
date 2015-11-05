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

//In general, better to make these command line options then defines
//same goes for hardcoding file_names, should be a command line option

// Enables residuals on fits
//#define DrawResiduals

// Enables polarisation weights
#define Corrections

// enables the pulls calculation
#define PULLS

void sPlotfit();
void sPlotFit() { sPlotfit(); }
void sPlotfit() {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  TF1 *f_straighline = new TF1("f_straighline", "0", 0, 1000);
  f_straighline->SetLineColor(kBlack);

  RooRealVar *inclusive_jpsi_onia_mass = new RooRealVar("inclusive_jpsi_onia_mass", "inclusive_jpsi_onia_mass", 2.85, 3.35);
  RooRealVar *inclusive_jpsi_onia_tau  = new RooRealVar("inclusive_jpsi_onia_tau", "inclusive_jpsi_onia_tau", -5, 10); // differs from binned fit (-0.3 to 10)
  RooRealVar *inclusive_jpsi_z_mass    = new RooRealVar("inclusive_jpsi_z_mass", "inclusive_jpsi_z_mass", 60., 120.);
  RooRealVar *inclusive_jpsi_is_z_to_electrons    = new RooRealVar("inclusive_jpsi_is_z_to_electrons", "inclusive_jpsi_is_z_to_electrons", 0, 1);
  RooRealVar *inclusive_jpsi_is_z_to_muons  = new RooRealVar("inclusive_jpsi_is_z_to_muons", "inclusive_jpsi_is_z_to_muons", 0, 1);

  RooArgSet zjpsi_argset(*onia_mass, *onia_tau, *is_z_to_electrons, *is_z_to_muons, *z_mass);
  RooArgSet jpsi_argset(*inclusive_jpsi_onia_mass, *inclusive_jpsi_onia_tau, *inclusive_jpsi_is_z_to_electrons, *inclusive_jpsi_is_z_to_muons, *inclusive_jpsi_z_mass);

  RooRealVar *onia_mass = new RooRealVar("onia_mass", "onia_mass", 2.85, 3.35, "GeV");
  RooRealVar *onia_tau  = new RooRealVar("onia_tau", "onia_tau", -5, 10, "ps");
  RooRealVar *z_mass    = new RooRealVar("z_mass", "z_mass", 60., 120., "GeV");
  RooRealVar *is_z_to_electrons    = new RooRealVar("is_z_to_electrons", "is_z_to_electrons", 0, 1);
  RooRealVar *is_z_to_muons  = new RooRealVar("is_z_to_muons", "is_z_to_muons", 0, 1);

  RooRealVar *onia_pt     = new RooRealVar("onia_pt",  "onia_pt", 0, 100);
  RooRealVar *onia_rap    = new RooRealVar("onia_rap",  "onia_rap", -5, 5);
  RooRealVar *onia_mu0_pt  = new RooRealVar("onia_mu0_pt",  "onia_mu0_pt", 0, 100);
  RooRealVar *onia_mu1_pt  = new RooRealVar("onia_mu1_pt",  "onia_mu1_pt", 0, 100);
  RooRealVar *onia_mu0_eta = new RooRealVar("onia_mu0_eta",  "onia_mu0_eta",-3, 3);
  RooRealVar *onia_mu1_eta = new RooRealVar("onia_mu1_eta",  "onia_mu1_eta",-3, 3);

  RooRealVar *unpolarised  = new RooRealVar("unpolarised" , "unpolarised" ,  -1000000000000000000000., 100000000000000000.);
  RooRealVar *longitudinal = new RooRealVar("longitudinal", "longitudinal",  -1000000000000000000000., 100000000000000000.);
  RooRealVar *transverse   = new RooRealVar("transverse"  , "transverse"  ,  -1000000000000000000000., 100000000000000000.);
  RooRealVar *transverse_pos  = new RooRealVar("transverse_pos" , "transverse_pos" ,  -1000000000000000000000., 100000000000000000.);
  RooRealVar *transverse_neg  = new RooRealVar("transverse_neg" , "transverse_neg" ,  -1000000000000000000000., 100000000000000000.);

  RooRealVar *reco_muon1_weight = new RooRealVar("reco_muon1_weight", "reco_muon1_weight", -10, 10);
  RooRealVar *reco_muon2_weight = new RooRealVar("reco_muon2_weight", "reco_muon2_weight", -10, 10);

  RooArgSet zjpsi_argset(*onia_mass, *onia_tau, *is_z_to_electrons, *is_z_to_muons, *z_mass, *onia_pt, *onia_rap);
  zjpsi_argset.add(*onia_mu0_pt); zjpsi_argset.add(*onia_mu1_pt);
  zjpsi_argset.add(*onia_mu0_eta); zjpsi_argset.add(*onia_mu1_eta);
  zjpsi_argset.add(*unpolarised); zjpsi_argset.add(*longitudinal);
  zjpsi_argset.add(*transverse); zjpsi_argset.add(*transverse_pos);
  zjpsi_argset.add(*transverse_neg);
  zjpsi_argset.add(*reco_muon1_weight); zjpsi_argset.add(*reco_muon2_weight);

  //Data Set
  TFile *zjpsi_ntuple = new TFile("../SkimNtuple/ZJpsi.root");
  TTree* zjpsi_tree = (TTree*) zjpsi_ntuple->Get("AUX");
  RooDataSet *zjpsi_data = new RooDataSet("zjpsi_data", "zjpsi_data", zjpsi_tree, zjpsi_argset);

  TFile *jpsi_ntuple = new TFile("ZJpsi.root");
  TFile *inclusive_jpsi_ntuple = new TFile("Jpsimumu.root");

  TTree* inclusive_jpsi_tree = (TTree*) inclusive_jpsi_ntuple->Get("AUX");
  RooDataSet *inclusive_jpsi_data = new RooDataSet("inclusive_jpsi_data", "inclusive_jpsi_data", inclusive_jpsi_tree, jpsi_argset);

  //      TCut SelectionCut = "onia_mass>2.6 && onia_mass<3.6";
  //  RooDataSet *newdata = (RooDataSet*)data->reduce(SelectionCut);

  // *****************************************************************
  //  Jpsi 
  // *****************************************************************

  RooRealVar jpsi_gaussmean("jpsi_gaussmean", "Mean of the smearing jpsi Gaussian", 0.0);
  RooRealVar jpsi_gausssigma("jpsi_gausssigma", "Width of the smearing jpsi Gaussian", 0.01, 0.005, 0.3);
  RooGaussModel jpsi_smear_gauss_model("jpsi_smear_gauss", "Gaussian used to smear the jpsi Exponential Decay", *inclusive_jpsi_onia_tau, jpsi_gaussmean, jpsi_gausssigma);
  RooRealVar jpsi_decay_lifetime("jpsi_decay_lifetime", "jpsi_Tau", 1.3, 0.05, 2.);
  RooDecay jpsi_decay_exp("jpsi_decay_exp", "jpsi_Exponential Decay", *inclusive_jpsi_onia_tau,jpsi_decay_lifetime,jpsi_smear_gauss_model,RooDecay::SingleSided);

  RooRealVar jpsi_gauss_prompt_mean("jpsi_gauss_prompt_mean", "Mean of the Prompt jpsi_Gaussian", 0., -0.1, 0.1);
  RooRealVar jpsi_gauss_prompt_sigma("jpsi_gauss_prompt_sigma", "Width of the Prompt jpsi_Gaussian", 0.007, 0.005, 0.25);
  RooGaussian jpsi_prompt_gauss("jpsi_prompt_gauss", "Gaussian of the Prompt jpsi_Peak", *inclusive_jpsi_onia_tau, jpsi_gauss_prompt_mean, jpsi_gauss_prompt_sigma);

  RooRealVar jpsi_gauss_prompt_sigma_2("jpsi_gauss_prompt_sigma_2", "Width of the Prompt jpsi_Gaussian_2", 0.01, 0.008, 0.4);
  RooGaussian jpsi_prompt_gauss_2("jpsi_prompt_gauss_2", "Gaussian_2 of the Prompt jpsi_Peak", *inclusive_jpsi_onia_tau, jpsi_gauss_prompt_mean, jpsi_gauss_prompt_sigma_2);

  RooRealVar jpsi_frac_prompt_sharp("jpsi_frac_prompt_sharp", "jpsi_frac_prompt_sharp" , 0.5 , 0.0, 1.);
  RooAddPdf jpsi_oniatau_gauss_sum_fitpdf("jpsi_oniatau_gauss_sum_fitpdf", "jpsi_oniatau_gauss_sum_fitpdf", RooArgList(jpsi_prompt_gauss,jpsi_prompt_gauss_2), RooArgList(jpsi_frac_prompt_sharp));

  RooRealVar jpsi_prompt_fraction("jpsi_prompt_fraction", "jpsi_prompt_fraction" , 0.01 , 0.0, 1);

  RooRealVar jpsi_sigma("jpsi_sigma", "jpsi_sigma", 0.05, 0.001, 0.1);
  RooRealVar jpsi_alpha("jpsi_alpha", "jpsi_alpha", 1.8, 1.0, 2.5);
  RooRealVar jpsi_n("jpsi_n", "jpsi_n", 2., 1.0, 80.);
  RooRealVar jpsi_dimuon_mean("jpsi_dimuon_mean", "jpsi_dimuon_mean", 3.1, 3.0, 3.2);
  RooCBShape jpsi_crystal_ball ("jpsi_crystal_ball", "jpsi_crystal_ball", *inclusive_jpsi_onia_mass, jpsi_dimuon_mean, jpsi_sigma, jpsi_alpha, jpsi_n );
  RooRealVar jpsi_dimuon_sigma("jpsi_dimuon_sigma", "jpsi_dimuon_sigma", 0.02, 0.001, 0.1);
  RooGaussian jpsi_dimuon_gauss ("jpsi_dimuon_gauss", "jpsi_dimuon_gauss", *inclusive_jpsi_onia_mass, jpsi_dimuon_mean, jpsi_dimuon_sigma);
  RooRealVar jpsi_frac_cball("jpsi_frac_cball", "jpsi_frac_cball" , 0.5 , 0.0, 1.);
  RooAddPdf jpsi_mass_signal("jpsi_mass_signal", "jpsi_mass_signal", RooArgList(jpsi_crystal_ball, jpsi_dimuon_gauss), RooArgList(jpsi_frac_cball));

  RooRealVar jpsi_dimuon_slope("jpsi_dimuon_slope", "jpsi_dimuon_slope", -0.1, -10., 10.);
  RooExponential jpsi_dimuon_bg_exponential("jpsi_dimuon_bg_exponential", "jpsi_dimuon_bg_exponential", *inclusive_jpsi_onia_mass, jpsi_dimuon_slope);

  RooRealVar jpsi_np_dimuon_slope("jpsi_np_dimuon_slope", "jpsi_np_dimuon_slope", -0.1, -10., 10.);
  RooExponential jpsi_np_dimuon_bg_exponential("jpsi_np_dimuon_bg_exponential", "jpsi_np_dimuon_bg_exponential", *inclusive_jpsi_onia_mass, jpsi_np_dimuon_slope);

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

  //does not work right now
  //RooFitResult *inclusive_jpsi_fr = inclusive_jpsi_model.fitTo(*inclusive_jpsi_data, NumCPU(4, kTRUE), Save());

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

  // *****************************************************************
  //  ZJpsi 
  // *****************************************************************

  // ZJpsi Unbinned fit start here

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
  RooGaussModel zjpsi_smear_gauss_model("zjpsi_smear_gauss", "Gaussian used to smear the zjpsi Exponential Decay", *onia_tau, zjpsi_gaussmean, zjpsi_gausssigma);
  RooRealVar zjpsi_decay_lifetime("zjpsi_decay_lifetime", "zjpsi_Tau", zjpsi_decay_lifetime_value, zjpsi_decay_lifetime_value-zjpsi_decay_lifetime_value_err, zjpsi_decay_lifetime_value+zjpsi_decay_lifetime_value_err);
  RooDecay zjpsi_decay_exp("zjpsi_decay_exp", "zjpsi_Exponential Decay", *onia_tau,zjpsi_decay_lifetime,zjpsi_smear_gauss_model,RooDecay::SingleSided);

  RooRealVar zjpsi_gauss_prompt_mean("zjpsi_gauss_prompt_mean", "Mean of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_mean_value, zjpsi_gauss_prompt_mean_value-zjpsi_gauss_prompt_mean_value_err, zjpsi_gauss_prompt_mean_value+zjpsi_gauss_prompt_mean_value_err);
  RooRealVar zjpsi_gauss_prompt_sigma("zjpsi_gauss_prompt_sigma", "Width of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_sigma_value, zjpsi_gauss_prompt_sigma_value-zjpsi_gauss_prompt_sigma_value_err, zjpsi_gauss_prompt_sigma_value+zjpsi_gauss_prompt_sigma_value_err);
  RooGaussian zjpsi_prompt_gauss("zjpsi_prompt_gauss", "Gaussian of the Prompt zjpsi_Peak", *onia_tau, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma);

  RooRealVar zjpsi_gauss_prompt_sigma_2("zjpsi_gauss_prompt_sigma_2", "Width of the Prompt zjpsi_Gaussian_2", zjpsi_gauss_prompt_sigma_2_value, zjpsi_gauss_prompt_sigma_2_value-zjpsi_gauss_prompt_sigma_2_value_err, zjpsi_gauss_prompt_sigma_2_value+zjpsi_gauss_prompt_sigma_2_value_err);
  RooGaussian zjpsi_prompt_gauss_2("zjpsi_prompt_gauss_2", "Gaussian_2 of the Prompt zjpsi_Peak", *onia_tau, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma_2);

  RooRealVar zjpsi_frac_prompt_sharp("zjpsi_frac_prompt_sharp", "zjpsi_frac_prompt_sharp" , zjpsi_frac_prompt_sharp_value, zjpsi_frac_prompt_sharp_value-zjpsi_frac_prompt_sharp_value_err, zjpsi_frac_prompt_sharp_value+zjpsi_frac_prompt_sharp_value_err);
  RooAddPdf onia_tau_gauss_sum_fitpdf("onia_tau_gauss_sum_fitpdf", "onia_tau_gauss_sum_fitpdf", RooArgList(zjpsi_prompt_gauss,zjpsi_prompt_gauss_2), RooArgList(zjpsi_frac_prompt_sharp));

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
  RooCBShape zjpsi_crystal_ball ("zjpsi_crystal_ball", "zjpsi_crystal_ball", *onia_mass, zjpsi_dimuon_mean, zjpsi_sigma, zjpsi_alpha, zjpsi_n );
  RooRealVar zjpsi_dimuon_sigma("zjpsi_dimuon_sigma", "zjpsi_dimuon_sigma", zjpsi_dimuon_sigma_value, zjpsi_dimuon_sigma_value-zjpsi_dimuon_sigma_value_err, zjpsi_dimuon_sigma_value+zjpsi_dimuon_sigma_value_err);
  RooGaussian zjpsi_dimuon_gauss ("zjpsi_dimuon_gauss", "zjpsi_dimuon_gauss", *onia_mass, zjpsi_dimuon_mean, zjpsi_dimuon_sigma);
  RooRealVar zjpsi_frac_cball("zjpsi_frac_cball", "zjpsi_frac_cball" , zjpsi_frac_cball_value, zjpsi_frac_cball_value-zjpsi_frac_cball_value_err, zjpsi_frac_cball_value+zjpsi_frac_cball_value_err);
  RooAddPdf zjpsi_mass_signal("zjpsi_mass_signal", "zjpsi_mass_signal", RooArgList(zjpsi_crystal_ball, zjpsi_dimuon_gauss), RooArgList(zjpsi_frac_cball));

  RooRealVar zjpsi_dimuon_slope("zjpsi_dimuon_slope", "zjpsi_dimuon_slope", zjpsi_dimuon_slope_value, zjpsi_dimuon_slope_value-zjpsi_dimuon_slope_value_err, zjpsi_dimuon_slope_value+zjpsi_dimuon_slope_value_err);
  RooExponential zjpsi_dimuon_bg_exponential("zjpsi_dimuon_bg_exponential", "zjpsi_dimuon_bg_exponential", *onia_mass, zjpsi_dimuon_slope);

  RooRealVar zjpsi_np_dimuon_slope("zjpsi_np_dimuon_slope", "zjpsi_np_dimuon_slope", zjpsi_np_dimuon_slope_value, zjpsi_np_dimuon_slope_value-zjpsi_np_dimuon_slope_value_err, zjpsi_np_dimuon_slope_value+zjpsi_np_dimuon_slope_value_err);
  RooExponential zjpsi_np_dimuon_bg_exponential("zjpsi_np_dimuon_bg_exponential", "zjpsi_np_dimuon_bg_exponential", *onia_mass, zjpsi_np_dimuon_slope);

//  RooRealVar zjpsi_m_sig_tau_sig_frac("zjpsi_m_sig_tau_sig_frac", "zjpsi_m_sig_tau_sig_frac", 0.2, 0.0, 1.0);
//  RooRealVar zjpsi_m_sig_tau_bg_frac("zjpsi_m_sig_tau_bg_frac", "zjpsi_m_sig_tau_bg_frac", 0.6, 0.0, 1.0);
//  RooRealVar zjpsi_m_bg_tau_sig_frac("zjpsi_m_bg_tau_sig_frac", "zjpsi_m_bg_tau_sig_frac", 0.15, 0.0, 1.0);

  // PDFs
  RooProdPdf zjpsi_m_sig_tau_sig("zjpsi_m_sig_tau_sig", "zjpsi_m_sig_tau_sig", RooArgList(zjpsi_mass_signal, onia_tau_gauss_sum_fitpdf ));
  RooProdPdf zjpsi_m_sig_tau_bg ("zjpsi_m_sig_tau_bg",  "zjpsi_m_sig_tau_bg",  RooArgList(zjpsi_mass_signal, zjpsi_decay_exp ));
  RooProdPdf zjpsi_m_bg_tau_bg  ("zjpsi_m_bg_tau_bg",   "zjpsi_m_bg_tau_bg",   RooArgList(zjpsi_np_dimuon_bg_exponential, zjpsi_decay_exp ));
  RooProdPdf zjpsi_m_bg_tau_sig ("zjpsi_m_bg_tau_sig",  "zjpsi_m_bg_tau_sig",  RooArgList(zjpsi_dimuon_bg_exponential, onia_tau_gauss_sum_fitpdf ));

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

//  RooAddPdf zjpsi_model("zjpsi_model", "zjpsi_model", RooArgSet(zjpsi_m_sig_tau_sig, zjpsi_m_sig_tau_bg, zjpsi_m_bg_tau_sig, zjpsi_m_bg_tau_bg ), RooArgList(zjpsi_m_sig_tau_sig_frac, zjpsi_m_sig_tau_bg_frac, zjpsi_m_bg_tau_sig_frac), kTRUE);
  RooAddPdf zjpsi_model("zjpsi_model", "zjpsi_model", RooArgList(ezjpsi_m_sig_tau_sig, ezjpsi_m_sig_tau_bg, ezjpsi_m_bg_tau_sig, ezjpsi_m_bg_tau_bg));

  RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data, NumCPU(4, kTRUE), Save());

#ifdef PULLS
  int Ntoys = 100;
  RooMCStudy *SimToy = new RooMCStudy(zjpsi_model, RooArgSet(*onia_mass, *onia_tau), Binned(kFALSE), Silence(), Extended(), FitOptions(Save(kTRUE)), PrintEvalErrors(0), Minos(kTRUE));

  SimToy->generateAndFit(Ntoys);

  RooPlot* ngauss_pull_frame        = SimToy->plotPull(Nzjpsi_m_sig_tau_sig, -4, 4, 16, kTRUE);
  RooPlot* nbkgPrompt_pull_frame    = SimToy->plotPull(Nzjpsi_m_bg_tau_sig,    -4, 4, 16, kTRUE);
  RooPlot* nsigNonPrompt_pull_frame = SimToy->plotPull(Nzjpsi_m_sig_tau_bg, -4, 4, 16, kTRUE);
  RooPlot* nbkgNonPrompt_pull_frame = SimToy->plotPull(Nzjpsi_m_bg_tau_bg, -4, 4, 16, kTRUE);

  TCanvas* c2 = new TCanvas("c2","toymc", 800, 800); c2->Divide(2,2);
  c2->cd(1); ngauss_pull_frame->Draw();         ngauss_pull_frame->SetXTitle("Yield Signal Prompt Pull");
  c2->cd(2); nbkgPrompt_pull_frame->Draw();     nbkgPrompt_pull_frame->SetXTitle("Yield Background Prompt Pull");
  c2->cd(3); nsigNonPrompt_pull_frame->Draw();  nsigNonPrompt_pull_frame->SetXTitle("Yield Signal Non Prompt Pull");
  c2->cd(4); nbkgNonPrompt_pull_frame->Draw();  nbkgNonPrompt_pull_frame->SetXTitle("Yield Background Non Prompt Pull");
#endif

	// *****************************
	// START OF JPSI PLOTTING
	// *****************************

  //----------- Mass Fit ---------------
  //Plots zjpsi_data on to frame
  RooPlot* mass_frame = onia_mass->frame(Bins(40));
  zjpsi_data->plotOn(mass_frame, Name("PlotData"));

  //Plots full model, prompt and non-prompt models to frame
#ifdef VisualizeError
  //zjpsi_model.plotOn(mass_frame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("masserror"), VisualizeError(*fr,3));
#endif
  //zjpsi_model.plotOn(mass_frame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("PlotModel"));
  //zjpsi_model.plotOn(mass_frame, Components("zjpsi_m_sig_tau_sig,zjpsi_m_bg_tau_sig"), LineWidth(3),LineColor(kAzure),LineStyle(2),Name("MassPrompt"));
  //zjpsi_model.plotOn(mass_frame, Components("zjpsi_m_bg_tau_sig"), LineWidth(3),LineColor(kAzure),LineStyle(2));
  //zjpsi_model.plotOn(mass_frame, Components("zjpsi_m_sig_tau_sig"), LineWidth(3),LineColor(kOrange));
  //zjpsi_model.plotOn(mass_frame, Components("zjpsi_m_sig_tau_bg,zjpsi_m_bg_tau_bg"), LineWidth(3),LineColor(kViolet),LineStyle(2), Name("MassNonPrompt"));

  zjpsi_model.plotOn(mass_frame, LineColor(kRed-2), RooFit::Name("total"));
  zjpsi_model.plotOn(mass_frame, Components(zjpsi_m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  zjpsi_model.plotOn(mass_frame, Components(zjpsi_m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  zjpsi_model.plotOn(mass_frame, Components(zjpsi_m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  zjpsi_model.plotOn(mass_frame, Components(zjpsi_m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));
  
  
  zjpsi_model.plotOn(mass_frame, Name("PlotData"));
	mass_frame->GetXaxis()->SetTitleOffset(.9);

  //RooCurve * mPJpsi 	 = mass_frame->getCurve("MassPrompt");
  //RooCurve * mNPJpsi   = mass_frame->getCurve("MassNonPrompt");

  RooCurve * mPJpsi 	 = mass_frame->getCurve("prompt j/psi");
  RooCurve * mNPJpsi   = mass_frame->getCurve("non-prompt j/psi");

  RooPlot* dummy_frame_jpsi = onia_mass->frame(Title("dummy frame to extract residuals"), Bins(40));
  zjpsi_data->plotOn(dummy_frame_jpsi); 
  zjpsi_model.plotOn(dummy_frame_jpsi);

  RooHist* h_residuals_mass_jpsi = dummy_frame_jpsi->pullHist();
  RooPlot* frame_residuals_mass_jpsi = onia_mass->frame(Title("Residual Distribution #mu^{+}#mu^{-} mass"));
	frame_residuals_mass_jpsi->GetYaxis()->SetTitle("(fit - data)/#sigma");
	frame_residuals_mass_jpsi->GetYaxis()->SetTitleSize(.16);
	frame_residuals_mass_jpsi->GetYaxis()->SetTitleOffset(.2);
  frame_residuals_mass_jpsi->addPlotable(h_residuals_mass_jpsi, "P");

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
  mass_frame->Draw();
  //mass_frame->SetTitle("J/#psi Mass [GeV]");

  Double_t xl1_l2=.58, yl1_l2=0.55, xl2_l2=xl1_l2+.3, yl2_l2=yl1_l2+.325;
  TLegend *leg2 = new TLegend(xl1_l2,yl1_l2,xl2_l2,yl2_l2);
  leg2->SetFillColor(kWhite);
  leg2->AddEntry(mass_frame->findObject("total"),"total","l");
  leg2->AddEntry(mass_frame->findObject("prompt j/psi"),"prompt j/psi","l");
  leg2->AddEntry(mass_frame->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg2->AddEntry(mass_frame->findObject("prompt continuum"),"prompt continuum","l");
  leg2->AddEntry(mass_frame->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg2->Draw();
  mass_frame->SetTitle("Z + J/#psi");

  //----------- Tau Fit ---------------
  //Plots data on to frame
  int timeBins = 20;
  RooPlot* lifetime_frame = onia_tau->frame(Bins(timeBins));
  zjpsi_data->plotOn(lifetime_frame, Name("Plotzjpsi_data"));

  //Plots full model, prompt and non-prompt models to frame
#ifdef VisualizeError
  //zjpsi_model.plotOn(lifetime_frame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("Timerror"), VisualizeError(*fr,3));
#endif
  //zjpsi_model.plotOn(lifetime_frame, LineColor(kBlue),LineWidth(3),NumCPU(1,kTRUE), Name("Plotzjpsi_modelTime"));
  //zjpsi_model.plotOn(lifetime_frame, Components("zjpsi_m_sig_tau_sig,zjpsi_m_bg_tau_sig"), LineWidth(3),LineColor(kAzure),LineStyle(2),Name("timePrompt"));
  //zjpsi_model.plotOn(lifetime_frame, Components("zjpsi_m_bg_tau_sig"), LineWidth(3),LineColor(kAzure),LineStyle(2));
  //zjpsi_model.plotOn(lifetime_frame, Components("zjpsi_m_sig_tau_sig"), LineWidth(3),LineColor(kOrange));
  //zjpsi_model.plotOn(lifetime_frame, Components("zjpsi_m_sig_tau_bg,zjpsi_m_bg_tau_bg"), LineWidth(3),LineColor(kViolet),LineStyle(2), Name("timeNonPrompt"));

  zjpsi_model.plotOn(lifetime_frame, LineColor(kRed-2), RooFit::Name("total"));
  zjpsi_model.plotOn(lifetime_frame, Components(zjpsi_m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  zjpsi_model.plotOn(lifetime_frame, Components(zjpsi_m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  zjpsi_model.plotOn(lifetime_frame, Components(zjpsi_m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  zjpsi_model.plotOn(lifetime_frame, Components(zjpsi_m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));

  zjpsi_model.plotOn(lifetime_frame, Name("Plotzjpsi_data"));
	lifetime_frame->GetXaxis()->SetTitleOffset(.9);

  //RooCurve * tPJpsi 	 = lifetime_frame->getCurve("timePrompt");
  //RooCurve * tNPJpsi   = lifetime_frame->getCurve("timeNonPrompt");
  RooCurve * tPJpsi 	 = lifetime_frame->getCurve("prompt j/psi");
  RooCurve * tNPJpsi   = lifetime_frame->getCurve("non-prompt j/psi");

  RooPlot* dummy_tframe_jpsi = onia_tau->frame(Title("dummy frame to extract residuals"), Bins(40));
  zjpsi_data->plotOn(dummy_tframe_jpsi); 
  zjpsi_model.plotOn(dummy_tframe_jpsi);

  RooHist* h_residuals_time_jpsi = dummy_tframe_jpsi->pullHist();
  RooPlot* frame_residuals_time_jpsi = onia_tau->frame(Title("Residual Distribution #mu^{+}#mu^{-} time"));
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
  lifetime_frame->Draw();
  Double_t xl1=.58, yl1=0.55, xl2=xl1+.3, yl2=yl1+.325;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillColor(kWhite);
  leg->AddEntry(lifetime_frame->findObject("total"),"total","l");
  leg->AddEntry(lifetime_frame->findObject("prompt j/psi"),"prompt j/psi","l");
  leg->AddEntry(lifetime_frame->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg->AddEntry(lifetime_frame->findObject("prompt continuum"),"prompt continuum","l");
  leg->AddEntry(lifetime_frame->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg->Draw();
  //lifetime_frame->SetTitle("J/#psi Mass [GeV]");
  lifetime_frame->SetTitle("Z + J/#psi");
  canvas2->SetLogy(1);

//// *******************************************************
//// *******************************************************
//// Splot
//// *******************************************************
//// *******************************************************

  RooWorkspace* ws = new RooWorkspace("myWS");
  ws->import(zjpsi_model);
  ws->import(*zjpsi_data, Rename("sPlotzjpsi_data"));
  
  RooAbsPdf* sPlotmodel = ws->pdf("zjpsi_model");
  RooRealVar* SignalPromptYield = ws->var("Nzjpsi_m_sig_tau_sig");
  RooRealVar* SignalNonPromptYield = ws->var("Nzjpsi_m_sig_tau_bg");
  RooRealVar* BckgPromptYield = ws->var("Nzjpsi_m_bg_tau_sig");
  RooRealVar* BckgNonPromptYield = ws->var("Nzjpsi_m_bg_tau_bg");
  
  SignalPromptYield->setConstant();
  SignalNonPromptYield->setConstant();
  BckgPromptYield->setConstant();
  BckgNonPromptYield->setConstant();
  
  RooStats::SPlot* szjpsi_data = new RooStats::SPlot("szjpsi_data","An SPlot", *zjpsi_data, sPlotmodel, RooArgList(*SignalPromptYield, *SignalNonPromptYield, *BckgPromptYield, *BckgNonPromptYield));
  
  RooDataSet * zjpsi_data_prompt_weighted     = new RooDataSet(zjpsi_data->GetName(), zjpsi_data->GetTitle(), zjpsi_data, *zjpsi_data->get(), 0, "Nzjpsi_m_sig_tau_sig_sw");
  RooDataSet * zjpsi_data_non_prompt_weighted = new RooDataSet(zjpsi_data->GetName(), zjpsi_data->GetTitle(), zjpsi_data, *zjpsi_data->get(), 0, "Nzjpsi_m_sig_tau_bg_sw");
  
  TCut SelecteZee   = "is_z_to_electrons==1";
  TCut SelecteZmumu = "is_z_to_muons==1";
  
  RooDataSet * zjpsi_data_prompt_Z_ee_weighted       = (RooDataSet*)zjpsi_data_prompt_weighted->reduce(SelecteZee);
  RooDataSet * zjpsi_data_non_prompt_Z_ee_weighted   = (RooDataSet*)zjpsi_data_non_prompt_weighted->reduce(SelecteZee);
  RooDataSet * zjpsi_data_prompt_Z_mumu_weighted     = (RooDataSet*)zjpsi_data_prompt_weighted->reduce(SelecteZmumu);
  RooDataSet * zjpsi_data_non_prompt_Z_mumu_weighted = (RooDataSet*)zjpsi_data_non_prompt_weighted->reduce(SelecteZmumu);
  
  int ZBins = 10;
  TCanvas* czjpsi_data_Z = new TCanvas("czjpsi_data_Z","sPlot Z", 1200, 1200); czjpsi_data_Z->Divide(2,2);
  RooPlot* frame_Z_ee_prompt = z_mass->frame(Bins(ZBins)); 
  zjpsi_data_prompt_Z_ee_weighted->plotOn(frame_Z_ee_prompt);
  czjpsi_data_Z->cd(1); 
  frame_Z_ee_prompt->SetTitle("Z->ee + Prompt J/#psi");
  frame_Z_ee_prompt->Draw();
  
  RooPlot* frame_Z_ee_non_prompt = z_mass->frame(Bins(ZBins));
  zjpsi_data_non_prompt_Z_ee_weighted->plotOn(frame_Z_ee_non_prompt);
  czjpsi_data_Z->cd(2); 
  frame_Z_ee_non_prompt->SetTitle("Z->ee + Nonprompt J/#psi");
  frame_Z_ee_non_prompt->Draw();
  
  RooPlot* frame_Z_mumu_prompt = z_mass->frame(Bins(ZBins)); 
  zjpsi_data_prompt_Z_mumu_weighted->plotOn(frame_Z_mumu_prompt);
  czjpsi_data_Z->cd(3); 
  frame_Z_mumu_prompt->SetTitle("Z->#mu#mu + Prompt J/#psi");
  frame_Z_mumu_prompt->Draw();
  
  RooPlot* frame_Z_mumu_non_prompt = z_mass->frame(Bins(ZBins)); 
  zjpsi_data_non_prompt_Z_mumu_weighted->plotOn(frame_Z_mumu_non_prompt);
  czjpsi_data_Z->cd(4); 
  frame_Z_mumu_non_prompt->SetTitle("Z->#mu#mu + Nonprompt J/#psi");
  frame_Z_mumu_non_prompt->Draw();
  
  // ************************
  // ************************
  // Start the x-sec plotting
  // ************************
  // ************************
  const RooArgSet* obs = zjpsi_data->get();
  RooRealVar* wtPROMPT = new RooRealVar("wtPROMPT" ,"wtPROMPT"  , -1000, 1000);
  RooRealVar* wtnonPROMPT = new RooRealVar("wtnonPROMPT" ,"wtnonPROMPT"  , -1000, 1000);
  RooDataSet *zjpsi_datatmpPROMPT   = new RooDataSet("zjpsi_datatmpPROMPT", "zjpsi_datatmpPROMPT", RooArgSet(*wtPROMPT));
  RooDataSet *zjpsi_datatmpnonPROMPT   = new RooDataSet("zjpsi_datatmpnonPROMPT", "zjpsi_datatmpnonPROMPT", RooArgSet(*wtnonPROMPT));
  for(Int_t i=0; i < zjpsi_data->numEntries(); i++) {
    obs = zjpsi_data->get(i);
#ifdef Corrections
    RooRealVar* effrow1 = (RooRealVar*) obs->find("reco_muon1_weight");
    RooRealVar* effrow2 = (RooRealVar*) obs->find("reco_muon2_weight");
    RooRealVar* accrow = (RooRealVar*) obs->find("unpolarised");
#else 
    RooRealVar* effrow1 = new RooRealVar("effrow1", "effrow1", 1.);
    RooRealVar* effrow2 = new RooRealVar("effrow2", "effrow2", 1.);
    RooRealVar* accrow = new RooRealVar("accrow", "accrow", 1.);
#endif
    double wtvalPROMPT    = szjpsi_data->GetSWeight(i,"Nzjpsi_m_sig_tau_sig_sw")*accrow->getVal()/(effrow1->getVal()*effrow2->getVal());
    double wtvalnonPROMPT = szjpsi_data->GetSWeight(i,"Nzjpsi_m_sig_tau_bg_sw")*accrow->getVal()/(effrow1->getVal()*effrow2->getVal());

    wtPROMPT->setVal(wtvalPROMPT);
    wtnonPROMPT->setVal(wtvalnonPROMPT);
    zjpsi_datatmpPROMPT->add(RooArgList(*wtPROMPT));
    zjpsi_datatmpnonPROMPT->add(RooArgList(*wtnonPROMPT));
  }
  
  zjpsi_data->merge(zjpsi_datatmpPROMPT); 
  zjpsi_data->merge(zjpsi_datatmpnonPROMPT); 
  RooDataSet * zjpsi_data_weightedPROMPT    = new RooDataSet(zjpsi_data->GetName(), zjpsi_data->GetTitle(), zjpsi_data, *zjpsi_data->get(), 0, wtPROMPT->GetName());  
  RooDataSet * zjpsi_data_weightednonPROMPT = new RooDataSet(zjpsi_data->GetName(), zjpsi_data->GetTitle(), zjpsi_data, *zjpsi_data->get(), 0, wtnonPROMPT->GetName());  

  RooBinning bins(6, 100) ;
  bins.addUniform(1, 6, 8.5);
  bins.addUniform(1, 8.5, 10);
  bins.addUniform(2, 10, 18);
  bins.addUniform(1, 18, 30);
  bins.addUniform(1, 30, 100);

  RooDataHist tmpHistPROMPT   ("tmpHistPROMPT"   , "tmpHistPROMPT", *onia_pt, *zjpsi_data_weightedPROMPT);
  RooDataHist tmpHistnonPROMPT("tmpHistnonPROMPT", "tmpHistnonPROMPT", *onia_pt, *zjpsi_data_weightednonPROMPT);

  TCanvas* c_jpsi_diff_xsec = new TCanvas("c_jpsi_diff_xsec","J/#psi p_{T} sPlot weighted", 1200, 600); c_jpsi_diff_xsec->Divide(2,1);
  // plot the prompt part
  c_jpsi_diff_xsec->cd(1);
  onia_pt->setBinning(bins);
  TH1* hhPROMPT  = tmpHistPROMPT.createHistogram("hhPROMPT", *onia_pt);
  TH1* hh2PROMPT = tmpHistPROMPT.createHistogram("hh2PROMPT", *onia_pt);
  hhPROMPT->SetLineWidth(2);
  hhPROMPT->Sumw2();
  for (int i=1; i<7; i++) {
    hh2PROMPT->SetBinContent(i, hhPROMPT->GetBinContent(i)/hhPROMPT->GetBinWidth(i));
    hh2PROMPT->SetBinError(i, hhPROMPT->GetBinError(i)/hhPROMPT->GetBinWidth(i));
  }
  TLatex splot_tex_pt;
  splot_tex_pt.SetTextSize(0.03);
  hh2PROMPT->Draw("E1");
  hh2PROMPT->GetYaxis()->SetRangeUser(0.01, 30);
  hh2PROMPT->GetXaxis()->SetRangeUser(7, 120);
  splot_tex_pt.DrawLatex(20, 10, "#bf{#it{CMS}} Preliminary" );
  splot_tex_pt.DrawLatex(20, 3 , "#sqrt{s}=8 TeV, #it{L} = 19.7 fb^{-1}" );
  c_jpsi_diff_xsec->cd(1)->SetLogy(1);
  c_jpsi_diff_xsec->cd(1)->SetLogx(1);
  // plot the non-prompt part
  c_jpsi_diff_xsec->cd(2);
  TH1* hhnonPROMPT  = tmpHistnonPROMPT.createHistogram("hhnonPROMPT", *onia_pt);
  TH1* hh2nonPROMPT = tmpHistnonPROMPT.createHistogram("hh2nonPROMPT", *onia_pt);
  hhnonPROMPT->SetLineWidth(2);
  hhnonPROMPT->Sumw2();
  for (int i=1; i<7; i++) {
    hh2nonPROMPT->SetBinContent(i, hhnonPROMPT->GetBinContent(i)/hhnonPROMPT->GetBinWidth(i));
    hh2nonPROMPT->SetBinError(i, hhnonPROMPT->GetBinError(i)/hhnonPROMPT->GetBinWidth(i));
  }
  hh2nonPROMPT->Draw("E1");
  hh2nonPROMPT->GetYaxis()->SetRangeUser(0.01, 30);
  hh2nonPROMPT->GetXaxis()->SetRangeUser(7, 120);
  splot_tex_pt.DrawLatex(20, 10, "#bf{#it{CMS}} Preliminary" );
  splot_tex_pt.DrawLatex(20, 3 , "#sqrt{s}=8 TeV, #it{L} = 19.7 fb^{-1}" );
  c_jpsi_diff_xsec->cd(2)->SetLogy(1);
  c_jpsi_diff_xsec->cd(2)->SetLogx(1);

}
