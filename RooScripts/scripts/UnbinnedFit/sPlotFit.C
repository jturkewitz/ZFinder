#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <RooBinning.h>
#include <RooWorkspace.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <TLatex.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooExtendPdf.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <TCut.h>
#include <TTree.h>
#include <TH1.h>
#include <TFile.h>
#include <TF1.h>
#include <RooCBShape.h>
#include <RooGaussModel.h>
#include <RooDecay.h>
#include <RooMCStudy.h>
#include <RooStats/SPlot.h>
#include <RooNumIntConfig.h>
#include <string>

using namespace RooFit ;
using namespace RooStats ;


//const float PDG_JPSI_MASS = 3.096916;

//In general, better to make these command line options then defines
//same goes for hardcoding file_names, should be a command line option

bool draw_residuals = false;
bool use_polarisation_weights = true;
//bool use_polarisation_weights = false;
bool create_cs_plots = true;
bool do_pull_calculation = false;
//bool use_pileup_correction = false; //subtract events when making ratio if pileup_correction is true
//bool use_pileup_correction = true;


//void sPlotfit(std::string file_name, bool use_pileup_correction = false);
//void sPlotFit(std::string file_name, bool use_pileup_correction = false, std::string image_dir = "~/public_html/ZPhysics/tmp/Test90/", int pt_slice = 0 ) { 
void sPlotFit(std::string file_name, bool use_pileup_correction = false, std::string image_dir = "~/public_html/ZPhysics/tmp/Test90/", bool use_only_z_to_muons = false, bool use_only_z_to_electrons = false) { 

  const double NUM_ZTOMUMU = 6.93504E06;
  const double NUM_ZTOEE = 4.63707E06;
  
  //std::vector<double> zjpsi_0 = do_fit(file_name, use_pileup_correction,image_dir,0); 
  //std::vector<double> zjpsi_1 = do_fit(file_name, use_pileup_correction,image_dir,1); 
  //std::vector<double> zjpsi_2 = do_fit(file_name, use_pileup_correction,image_dir,2); 
  //std::vector<double> zjpsi_3 = do_fit(file_name, use_pileup_correction,image_dir,3); 
  //std::vector<double> zjpsi_4 = do_fit(file_name, use_pileup_correction,image_dir,4); 
  //std::vector<double> zjpsi_5 = do_fit(file_name, use_pileup_correction,image_dir,5); 

  double zjpsi_prompt_events_all = 0.0;
  double zjpsi_prompt_events_all_error = 0.0;
  double zjpsi_prompt_events_pt_summed = 0.0;
  double zjpsi_prompt_events_pt_summed_sq_error = 0.0;

  double zjpsi_non_prompt_events_all = 0.0;
  double zjpsi_non_prompt_events_all_error = 0.0;
  double zjpsi_non_prompt_events_pt_summed = 0.0;
  double zjpsi_non_prompt_events_pt_summed_sq_error = 0.0;

  double zjpsi_prompt_events[6];
  double zjpsi_prompt_events_error[6];
  double zjpsi_non_prompt_events[6];
  double zjpsi_non_prompt_events_error[6];


  for (int i=0 ; i < 6 ; ++i) {
    //for now, if i==0, jpsi_all

      std::vector<double> zjpsi_info = do_fit(file_name, use_pileup_correction,image_dir, use_only_z_to_muons, use_only_z_to_electrons , i); 
      zjpsi_prompt_events[i] = zjpsi_info[0];
      zjpsi_prompt_events_error[i] = zjpsi_info[1];
      zjpsi_non_prompt_events[i] = zjpsi_info[2];
      zjpsi_non_prompt_events_error[i] = zjpsi_info[3];
      if (i == 0 ) {
        zjpsi_prompt_events_all = zjpsi_info[0] ;
        zjpsi_prompt_events_all_error = zjpsi_info[1];
        zjpsi_non_prompt_events_all = zjpsi_info[2] ;
        zjpsi_non_prompt_events_all_error = zjpsi_info[3];
      }
      else {
        zjpsi_prompt_events_pt_summed += zjpsi_info[0] ;
        zjpsi_prompt_events_pt_summed_sq_error += pow(zjpsi_info[1],2.0) ;
        zjpsi_non_prompt_events_pt_summed += zjpsi_info[2] ;
        zjpsi_non_prompt_events_pt_summed_sq_error += pow(zjpsi_info[3], 2.0);
      }
  }

  std::cout << "zjpsi_prompt_events_pt_all: " << zjpsi_prompt_events_all << std::endl;
  std::cout << "zjpsi_prompt_events_pt_all_error: " << zjpsi_prompt_events_all_error << std::endl;
  std::cout << "zjpsi_non_prompt_events_pt_all: " << zjpsi_non_prompt_events_all << std::endl;
  std::cout << "zjpsi_non_prompt_events_pt_all_error: " << zjpsi_non_prompt_events_all_error << std::endl;

  std::cout << "zjpsi_prompt_events_pt_summed: " << zjpsi_prompt_events_pt_summed << std::endl;
  std::cout << "zjpsi_prompt_events_pt_summed_error: " << pow(zjpsi_prompt_events_pt_summed_sq_error , 0.5) << std::endl;
  std::cout << "zjpsi_non_prompt_events_pt_summed: " << zjpsi_non_prompt_events_pt_summed << std::endl;
  std::cout << "zjpsi_non_prompt_events_pt_summed_error: " << pow(zjpsi_non_prompt_events_pt_summed_sq_error , 0.5) << std::endl;

  if (use_only_z_to_electrons) {
    std::cout << "ratio_prompt_pt_summed: " << zjpsi_prompt_events_pt_summed / (NUM_ZTOEE ) << std::endl;
    std::cout << "ratio_prompt_pt_summed_error: " << pow(zjpsi_prompt_events_pt_summed_sq_error , 0.5) / (NUM_ZTOEE ) << std::endl;

    std::cout << "ratio_nonprompt_pt_summed: " << zjpsi_non_prompt_events_pt_summed / (NUM_ZTOEE ) << std::endl;
    std::cout << "ratio_nonprompt_pt_summed_error: " << pow(zjpsi_non_prompt_events_pt_summed_sq_error , 0.5) / (NUM_ZTOEE ) << std::endl;
  }
  else if (use_only_z_to_muons) {

    std::cout << "ratio_prompt_pt_summed: " << zjpsi_prompt_events_pt_summed / (NUM_ZTOMUMU) << std::endl;
    std::cout << "ratio_prompt_pt_summed_error: " << pow(zjpsi_prompt_events_pt_summed_sq_error , 0.5) / (NUM_ZTOMUMU) << std::endl;

    std::cout << "ratio_nonprompt_pt_summed: " << zjpsi_non_prompt_events_pt_summed / (NUM_ZTOMUMU) << std::endl;
    std::cout << "ratio_nonprompt_pt_summed_error: " << pow(zjpsi_non_prompt_events_pt_summed_sq_error , 0.5) / (NUM_ZTOMUMU) << std::endl;
  }
  else {
    std::cout << "ratio_prompt_pt_summed: " << zjpsi_prompt_events_pt_summed / (NUM_ZTOMUMU + NUM_ZTOEE ) << std::endl;
    std::cout << "ratio_prompt_pt_summed_error: " << pow(zjpsi_prompt_events_pt_summed_sq_error , 0.5) / (NUM_ZTOMUMU + NUM_ZTOEE ) << std::endl;

    std::cout << "ratio_nonprompt_pt_summed: " << zjpsi_non_prompt_events_pt_summed / (NUM_ZTOMUMU + NUM_ZTOEE ) << std::endl;
    std::cout << "ratio_nonprompt_pt_summed_error: " << pow(zjpsi_non_prompt_events_pt_summed_sq_error , 0.5) / (NUM_ZTOMUMU + NUM_ZTOEE ) << std::endl;
  }

  std::cout << "i = 0 => all"  << std::endl;
  std::cout << "i = 1 => 8.5-10" << std::endl;
  std::cout << "i = 2 => 10-14" << std::endl;
  std::cout << "i = 3 => 14-18" << std::endl;
  std::cout << "i = 4 => 18-30" << std::endl;
  std::cout << "i = 5 => 30-100" << std::endl;

}
std::vector<double> do_fit(std::string file_name, bool use_pileup_correction = false, std::string image_dir = "~/public_html/ZPhysics/tmp/Test90/", bool use_only_z_to_muons = true, bool use_only_z_to_electrons = true, int pt_slice = 0 ) {
  string PT_SLICE_STRING;          //The string
  ostringstream temp;  //temp as in temporary
  temp<<pt_slice;
  PT_SLICE_STRING=temp.str();      //str is temp as string

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  //RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  //gErrorIgnoreLevel = kError;

  //Increase default precision of numeric integration
  // as this exercise has high sensitivity to numeric integration precision
  //RooAbsPdf::defaultIntegratorConfig()->setEpsRel(1e-8) ;
  //RooAbsPdf::defaultIntegratorConfig()->setEpsAbs(1e-8) ;
  RooAbsPdf::defaultIntegratorConfig()->setEpsRel(1e-12) ;
  RooAbsPdf::defaultIntegratorConfig()->setEpsAbs(1e-12) ;

  TLatex splot_tex_pt;
  splot_tex_pt.SetTextSize(0.04);
  splot_tex_pt.SetTextFont(42);

  TF1 *f_straighline = new TF1("f_straighline", "0", -100, 1000);
  f_straighline->SetLineColor(kBlack);

  RooRealVar *onia_mass = new RooRealVar("onia_mass", "J/#psi#rightarrow#mu#mu [GeV]", 2.85, 3.35);
  //RooRealVar *onia_tau  = new RooRealVar("onia_tau", "J/#psi#rightarrow#mu#mu pseudo-proper time [ps]", -2., 5);
  //RooRealVar *onia_tau  = new RooRealVar("onia_tau", "J/#psi#rightarrow#mu#mu pseudo-proper time [ps]", -0.3, 5);
  RooRealVar *onia_tau  = new RooRealVar("onia_tau", "J/#psi#rightarrow#mu#mu pseudo-proper time [ps]", -1.0, 5);
  //RooRealVar *z_mass    = new RooRealVar("z_mass", "Z Mass [GeV]", 80., 100.);
  //RooRealVar *z_mass    = new RooRealVar("z_mass", "Z Mass [GeV]", 150., 300.);
  //RooRealVar *z_mass    = new RooRealVar("z_mass", "Z Mass [GeV]", 40., 50.);
  //RooRealVar *z_mass    = new RooRealVar("z_mass", "Z Mass [GeV]", 40., 45.);
  RooRealVar *z_mass    = new RooRealVar("z_mass", "Z Mass [GeV]", 80., 100.);
  //RooRealVar *z_mass    = new RooRealVar("z_mass", "Z Mass [GeV]", 40., 100.);
  //RooRealVar *z_mass    = new RooRealVar("z_mass", "Z Mass [GeV]", 50., 150.);
  RooRealVar *is_z_to_electrons    = new RooRealVar("is_z_to_electrons", "is_z_to_electrons", 0, 1);
  RooRealVar *is_z_to_muons  = new RooRealVar("is_z_to_muons", "is_z_to_muons", 0, 1);

  RooRealVar *onia_pt     = new RooRealVar("onia_pt",  "onia_pt", 0, 100);
  RooRealVar *onia_rap    = new RooRealVar("onia_rap",  "onia_rap", -5, 5);
  RooRealVar *onia_mu0_pt  = new RooRealVar("onia_mu0_pt",  "onia_mu0_pt", 0, 100);
  RooRealVar *onia_mu1_pt  = new RooRealVar("onia_mu1_pt",  "onia_mu1_pt", 0, 100);
  RooRealVar *onia_mu0_eta = new RooRealVar("onia_mu0_eta",  "onia_mu0_eta",-3, 3);
  RooRealVar *onia_mu1_eta = new RooRealVar("onia_mu1_eta",  "onia_mu1_eta",-3, 3);

  RooRealVar *unpolarised  = new RooRealVar("unpolarised" , "unpolarised" ,  -10000000000000., 1000000000.);
  RooRealVar *lambda_pos  = new RooRealVar("lambda_pos" , "lambda_pos" ,  -10000000000000., 1000000000.);
  RooRealVar *lambda_neg  = new RooRealVar("lambda_neg" , "lambda_neg" ,  -10000000000000., 1000000000.);
  RooRealVar *longitudinal = new RooRealVar("longitudinal", "longitudinal",  -10000000000000., 1000000000.);
  RooRealVar *transverse   = new RooRealVar("transverse"  , "transverse"  ,  -10000000000000., 1000000000.);
  RooRealVar *transverse_pos  = new RooRealVar("transverse_pos" , "transverse_pos" ,  -10000000000000., 1000000000.);
  RooRealVar *transverse_neg  = new RooRealVar("transverse_neg" , "transverse_neg" ,  -10000000000000., 1000000000.);

  //RooRealVar *reco_muon0_weight = new RooRealVar("reco_muon0_weight", "reco_muon0_weight", -10, 10);
  //RooRealVar *reco_muon1_weight = new RooRealVar("reco_muon1_weight", "reco_muon1_weight", -10, 10);

  RooArgSet jpsi_argset(*onia_mass, *onia_tau, *onia_pt);
  jpsi_argset.add(*unpolarised); 
  jpsi_argset.add(*lambda_pos); 
  jpsi_argset.add(*lambda_neg); 

  RooArgSet zjpsi_argset(*onia_mass, *onia_tau, *is_z_to_electrons, *is_z_to_muons, *z_mass, *onia_pt, *onia_rap);
  zjpsi_argset.add(*onia_mu0_pt); zjpsi_argset.add(*onia_mu1_pt);
  zjpsi_argset.add(*onia_mu0_eta); zjpsi_argset.add(*onia_mu1_eta);
  zjpsi_argset.add(*unpolarised); 
  zjpsi_argset.add(*lambda_pos); 
  zjpsi_argset.add(*lambda_neg); 
  zjpsi_argset.add(*longitudinal);
  zjpsi_argset.add(*transverse); 
  zjpsi_argset.add(*transverse_pos);
  zjpsi_argset.add(*transverse_neg);
  //zjpsi_argset.add(*reco_muon0_weight);
  //zjpsi_argset.add(*reco_muon1_weight);

  //Data Set
  //TFile *zjpsi_ntuple = new TFile("../SkimNtuple/ZJpsi.root");
  //TFile *zjpsi_ntuple = new TFile("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/SkimNtuple/ZJpsi.root");


  TFile *zjpsi_ntuple = new TFile(file_name.c_str());
  //TFile *zjpsi_ntuple = new TFile("ZJpsi_40_300.root");
  //TFile *zjpsi_ntuple = new TFile("ZJpsi.root");
  //TFile *zjpsi_ntuple = new TFile("ZJpsiee.root");
  //TFile *zjpsi_ntuple = new TFile("ZJpsimumu.root");
  
  //TFile *zjpsi_ntuple = new TFile("ZJpsi_pileup.root");
  //TFile *zjpsi_ntuple = new TFile("ZJpsi_sideband.root");
  
  TTree* zjpsi_tree = (TTree*) zjpsi_ntuple->Get("AUX");
  RooDataSet *zjpsi_data     = new RooDataSet("zjpsi_data", "zjpsi_data", zjpsi_tree, zjpsi_argset);

  RooDataSet *zjpsi_data_weighted_testing     = new RooDataSet("zjpsi_data_weighted_testing", "zjpsi_data_weighted_testing", zjpsi_tree, zjpsi_argset,0, "unpolarised");

  RooDataSet *zjpsi_data_fid = new RooDataSet("zjpsi_data_fid", "zjpsi_data_fid", zjpsi_tree, zjpsi_argset);

  //TFile *zjpsi_ntuple = new TFile("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/SkimNtuple/ZJpsi.root");
  //TFile *ntuple = new TFile("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/SkimNtuple/jpsimumu.root");
  TFile *jpsi_ntuple = new TFile("Inclusive_Jpsi.root");

  TTree* jpsi_tree = (TTree*) jpsi_ntuple->Get("AUX");
  RooDataSet *jpsi_data = new RooDataSet("jpsi_data", "jpsi_data", jpsi_tree, jpsi_argset);
  RooDataSet *jpsi_data_weighted_testing     = new RooDataSet("jpsi_data_weighted_testing", "jpsi_data_weighted_testing", jpsi_tree, jpsi_argset, 0, "unpolarised");
  //RooDataSet *jpsi_data_weighted_testing     = new RooDataSet("jpsi_data_weighted_testing", "jpsi_data_weighted_testing", jpsi_tree, jpsi_argset, 0, "lambda_pos");
  //RooDataSet *jpsi_data_weighted_testing     = new RooDataSet("jpsi_data_weighted_testing", "jpsi_data_weighted_testing", jpsi_tree, jpsi_argset, 0, "lambda_neg");

  //      TCut SelectionCut = "onia_mass>2.6 && onia_mass<3.6";
  //  RooDataSet *newdata = (RooDataSet*)data->reduce(SelectionCut);

  TCut SelectionCut_1 = "onia_pt>=8.5 && onia_pt < 10.0";
  TCut SelectionCut_2 = "onia_pt>=10 && onia_pt < 14.0";
  TCut SelectionCut_3 = "onia_pt>=14 && onia_pt < 18.0";
  TCut SelectionCut_4 = "onia_pt>=18 && onia_pt < 30.0";
  TCut SelectionCut_5 = "onia_pt>=30 && onia_pt < 100.0";

  TCut SelectionCut_zmuons = "is_z_to_muons == 1";
  TCut SelectionCut_zelectrons = "is_z_to_electrons == 1";

  if (use_only_z_to_muons) {
    RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*) zjpsi_data_weighted_testing->reduce(SelectionCut_zmuons);
    RooDataSet *zjpsi_data = (RooDataSet*) zjpsi_data->reduce(SelectionCut_zmuons);
  }
  if (use_only_z_to_electrons) {
    RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*) zjpsi_data_weighted_testing->reduce(SelectionCut_zelectrons);
    RooDataSet *zjpsi_data = (RooDataSet*) zjpsi_data->reduce(SelectionCut_zelectrons);
  }

  //TCut SelectionCut_sb_all = "(is_z_to_muons == 1 && z_mass <= 45) || (is_z_to_muons == 1 && z_mass <= 50) || (z_mass >= 150)";
  ////TCut SelectionCut_sb_all = "z_mass >= 150";
  //RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*) zjpsi_data_weighted_testing->reduce(SelectionCut_sb_all);
  //RooDataSet *zjpsi_data = (RooDataSet*) zjpsi_data->reduce(SelectionCut_sb_all);


  if (pt_slice == 0) {
    //continue;
    cout << "don't bin shape parameters in pt" << endl;
  }
  else if (pt_slice == 1) {
    RooDataSet *jpsi_data = (RooDataSet*)jpsi_data->reduce(SelectionCut_1);
    RooDataSet *jpsi_data_weighted_testing = (RooDataSet*)jpsi_data_weighted_testing->reduce(SelectionCut_1);
    RooDataSet *zjpsi_data = (RooDataSet*)zjpsi_data->reduce(SelectionCut_1);
    RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*)zjpsi_data_weighted_testing->reduce(SelectionCut_1);
    cout << "pt bin 8.5-10" << endl;
  }
  else if (pt_slice == 2) {
    RooDataSet *jpsi_data = (RooDataSet*)jpsi_data->reduce(SelectionCut_2);
    RooDataSet *jpsi_data_weighted_testing = (RooDataSet*)jpsi_data_weighted_testing->reduce(SelectionCut_2);
    RooDataSet *zjpsi_data = (RooDataSet*)zjpsi_data->reduce(SelectionCut_2);
    RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*)zjpsi_data_weighted_testing->reduce(SelectionCut_2);
    cout << "pt bin 10-14" << endl;
  }
  else if (pt_slice == 3) {
    RooDataSet *jpsi_data = (RooDataSet*)jpsi_data->reduce(SelectionCut_3);
    RooDataSet *jpsi_data_weighted_testing = (RooDataSet*)jpsi_data_weighted_testing->reduce(SelectionCut_3);
    RooDataSet *zjpsi_data = (RooDataSet*)zjpsi_data->reduce(SelectionCut_3);
    RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*)zjpsi_data_weighted_testing->reduce(SelectionCut_3);
    cout << "pt bin 14-18" << endl;
  }
  else if (pt_slice == 4) {
    RooDataSet *jpsi_data = (RooDataSet*)jpsi_data->reduce(SelectionCut_4);
    RooDataSet *jpsi_data_weighted_testing = (RooDataSet*)jpsi_data_weighted_testing->reduce(SelectionCut_4);
    RooDataSet *zjpsi_data = (RooDataSet*)zjpsi_data->reduce(SelectionCut_4);
    RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*)zjpsi_data_weighted_testing->reduce(SelectionCut_4);
    cout << "pt bin 18-30" << endl;
  }
  else if (pt_slice == 5) {
    RooDataSet *jpsi_data = (RooDataSet*)jpsi_data->reduce(SelectionCut_5);
    RooDataSet *jpsi_data_weighted_testing = (RooDataSet*)jpsi_data_weighted_testing->reduce(SelectionCut_5);
    RooDataSet *zjpsi_data = (RooDataSet*)zjpsi_data->reduce(SelectionCut_5);
    RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*)zjpsi_data_weighted_testing->reduce(SelectionCut_5);
    cout << "pt bin 30-100" << endl;
  }


  //TCut SelectionCut_Z_Mass = "z_mass>=80 && z_mass <= 100.0";
  //RooDataSet *zjpsi_data_weighted_testing = (RooDataSet*)zjpsi_data_weighted_testing->reduce(SelectionCut_Z_Mass);

  // *****************************************************************
  //  Jpsi 
  // *****************************************************************

  RooRealVar jpsi_gaussmean("jpsi_gaussmean", "Mean of the smearing jpsi Gaussian", 0.0);
  RooRealVar jpsi_gausssigma("jpsi_gausssigma", "Width of the smearing jpsi Gaussian", 0.01, 0.00005, 0.3);
  RooGaussModel jpsi_smear_gauss_model("jpsi_smear_gauss", "Gaussian used to smear the jpsi Exponential Decay", *onia_tau, jpsi_gaussmean, jpsi_gausssigma);
  RooRealVar jpsi_decay_time("jpsi_decay_time", "jpsi_Tau", 1.3, 0.5, 2.0);
  RooDecay jpsi_decay_exp("jpsi_decay_exp", "jpsi_Exponential Decay", *onia_tau,jpsi_decay_time,jpsi_smear_gauss_model,RooDecay::SingleSided);

  RooRealVar jpsi_decay_time_continuum("jpsi_decay_time_continuum", "jpsi_Tau_continuum", 1.3, 0.5, 2.0);
  RooDecay jpsi_decay_exp_continuum("jpsi_decay_exp_continuum", "jpsi_Exponential Decay_continuum", *onia_tau,jpsi_decay_time_continuum,jpsi_smear_gauss_model,RooDecay::SingleSided);

  RooRealVar jpsi_gauss_prompt_mean("jpsi_gauss_prompt_mean", "Mean of the Prompt jpsi_Gaussian", 0.0, -0.1, 0.1);
  RooRealVar jpsi_gauss_prompt_sigma("jpsi_gauss_prompt_sigma", "Width of the Prompt jpsi_Gaussian", 0.008, 0.005, 0.25);
  RooGaussian jpsi_prompt_gauss("jpsi_prompt_gauss", "Gaussian of the Prompt jpsi_Peak", *onia_tau, jpsi_gauss_prompt_mean, jpsi_gauss_prompt_sigma);

  RooRealVar jpsi_gauss_prompt_sigma_2("jpsi_gauss_prompt_sigma_2", "Width of the Prompt jpsi_Gaussian_2", 0.25, 0.10, 0.4);
  RooGaussian jpsi_prompt_gauss_2("jpsi_prompt_gauss_2", "Gaussian_2 of the Prompt jpsi_Peak", *onia_tau, jpsi_gauss_prompt_mean, jpsi_gauss_prompt_sigma_2);

  RooRealVar jpsi_frac_prompt_sharp("jpsi_frac_prompt_sharp", "jpsi_frac_prompt_sharp" , 0.5, 0.0, 1.0);
  RooAddPdf jpsi_oniatau_gauss_sum_fitpdf("jpsi_oniatau_gauss_sum_fitpdf", "jpsi_oniatau_gauss_sum_fitpdf", RooArgList(jpsi_prompt_gauss,jpsi_prompt_gauss_2), RooArgList(jpsi_frac_prompt_sharp));

  //RooRealVar jpsi_sigma("jpsi_sigma", "jpsi_sigma", 0.05, 0.0002, 0.1);
  //RooRealVar jpsi_alpha("jpsi_alpha", "jpsi_alpha", 1.8, 1.0, 2.5);
  //RooRealVar jpsi_n("jpsi_n", "jpsi_n", 2., 1.0, 80.);

  RooRealVar jpsi_sigma("jpsi_sigma", "jpsi_sigma", 0.05, 0.008, 0.1);
  RooRealVar jpsi_alpha("jpsi_alpha", "jpsi_alpha", 1.8, 1.0, 5.5);
  RooRealVar jpsi_n("jpsi_n", "jpsi_n", 20., 1.0, 200.);

  RooRealVar jpsi_dimuon_mean("jpsi_dimuon_mean", "jpsi_dimuon_mean", 3.1, 3.0, 3.2);
  RooCBShape jpsi_crystal_ball ("jpsi_crystal_ball", "jpsi_crystal_ball", *onia_mass, jpsi_dimuon_mean, jpsi_sigma, jpsi_alpha, jpsi_n );
  RooRealVar jpsi_dimuon_sigma("jpsi_dimuon_sigma", "jpsi_dimuon_sigma", 0.02, 0.001, 0.1);
  RooGaussian jpsi_dimuon_gauss ("jpsi_dimuon_gauss", "jpsi_dimuon_gauss", *onia_mass, jpsi_dimuon_mean, jpsi_dimuon_sigma);
  RooRealVar jpsi_frac_cball("jpsi_frac_cball", "jpsi_frac_cball" , 0.5, 0.0, 1.0);
  RooAddPdf jpsi_mass_signal("jpsi_mass_signal", "jpsi_mass_signal", RooArgList(jpsi_crystal_ball, jpsi_dimuon_gauss), RooArgList(jpsi_frac_cball));

  RooRealVar jpsi_dimuon_slope("jpsi_dimuon_slope", "jpsi_dimuon_slope", -0.1, -10, 10);
  RooExponential jpsi_dimuon_bg_exponential("jpsi_dimuon_bg_exponential", "jpsi_dimuon_bg_exponential", *onia_mass, jpsi_dimuon_slope);

  RooRealVar jpsi_np_dimuon_slope("jpsi_np_dimuon_slope", "jpsi_np_dimuon_slope", -0.1, -10.0, 10.0);
  RooExponential jpsi_np_dimuon_bg_exponential("jpsi_np_dimuon_bg_exponential", "jpsi_np_dimuon_bg_exponential", *onia_mass, jpsi_np_dimuon_slope);

  // PDFs
  RooProdPdf jpsi_m_sig_tau_sig("jpsi_m_sig_tau_sig", "jpsi_m_sig_tau_sig", RooArgList(jpsi_mass_signal, jpsi_oniatau_gauss_sum_fitpdf ));
  RooProdPdf jpsi_m_sig_tau_bg ("jpsi_m_sig_tau_bg",  "jpsi_m_sig_tau_bg",  RooArgList(jpsi_mass_signal, jpsi_decay_exp ));
  RooProdPdf jpsi_m_bg_tau_bg  ("jpsi_m_bg_tau_bg",   "jpsi_m_bg_tau_bg",   RooArgList(jpsi_np_dimuon_bg_exponential, jpsi_decay_exp_continuum ));
  RooProdPdf jpsi_m_bg_tau_sig ("jpsi_m_bg_tau_sig",  "jpsi_m_bg_tau_sig",  RooArgList(jpsi_dimuon_bg_exponential, jpsi_oniatau_gauss_sum_fitpdf ));

  //yield ranges probably need to be changed?
  // yields
  RooRealVar Njpsi_m_sig_tau_sig("Njpsi_m_sig_tau_sig", "Njpsi_m_sig_tau_sig", 10000, 0, 1000000);  
  RooRealVar Njpsi_m_sig_tau_bg ("Njpsi_m_sig_tau_bg",  "Njpsi_m_sig_tau_bg",  4000, 0, 1000000); 
  RooRealVar Njpsi_m_bg_tau_bg  ("Njpsi_m_bg_tau_bg",   "Njpsi_m_bg_tau_bg",   1000, 0, 1000000); 
  RooRealVar Njpsi_m_bg_tau_sig ("Njpsi_m_bg_tau_sig",  "Njpsi_m_bg_tau_sig",  400, 0, 1000000); 

  // extended PDFs
  RooExtendPdf ejpsi_m_sig_tau_sig("ejpsi_m_sig_tau_sig", "ejpsi_m_sig_tau_sig", jpsi_m_sig_tau_sig, Njpsi_m_sig_tau_sig); 
  RooExtendPdf ejpsi_m_sig_tau_bg ("ejpsi_m_sig_tau_bg",  "ejpsi_m_sig_tau_bg",  jpsi_m_sig_tau_bg , Njpsi_m_sig_tau_bg );
  RooExtendPdf ejpsi_m_bg_tau_bg  ("ejpsi_m_bg_tau_bg",   "ejpsi_m_bg_tau_bg",   jpsi_m_bg_tau_bg  , Njpsi_m_bg_tau_bg  );
  RooExtendPdf ejpsi_m_bg_tau_sig ("ejpsi_m_bg_tau_sig",  "ejpsi_m_bg_tau_sig",  jpsi_m_bg_tau_sig , Njpsi_m_bg_tau_sig );

  //  RooAddPdf model("model", "model", RooArgSet(jpsi_m_sig_tau_sig, jpsi_m_sig_tau_bg, jpsi_m_bg_tau_sig, jpsi_m_bg_tau_bg ), RooArgList(jpsi_m_sig_tau_sig_frac, jpsi_m_sig_tau_bg_frac, jpsi_m_bg_tau_sig_frac), kTRUE);
  RooAddPdf jpsi_model("jpsi_model", "jpsi_model", RooArgList(ejpsi_m_sig_tau_sig, ejpsi_m_sig_tau_bg, ejpsi_m_bg_tau_sig, ejpsi_m_bg_tau_bg));

  //does not work right now
  //RooFitResult *fr_inclusive = jpsi_model.fitTo(*jpsi_data, NumCPU(4, kTRUE), Verbose(false), PrintLevel(-1), Save());
  RooFitResult *fr_inclusive = jpsi_model.fitTo(*jpsi_data_weighted_testing, NumCPU(8, kTRUE), Verbose(false), SumW2Error(kTRUE), PrintLevel(-1), Save());
  fr_inclusive->Print();

  int num_mass_bins = 40;
  int num_time_bins = 40;

  RooPlot* dummy_m = onia_mass->frame(Bins(num_mass_bins));
  jpsi_data->plotOn(dummy_m);
  jpsi_model.plotOn(dummy_m);

  //Plot inclusive j/psi mass frame
  RooPlot* inclusive_jpsi_mass_frame = onia_mass->frame(Bins(num_mass_bins));
  jpsi_data->plotOn(inclusive_jpsi_mass_frame, Name("jpsi_mass_data"));

  //Plots full model, prompt and non-prompt models to frame
  jpsi_model.plotOn(inclusive_jpsi_mass_frame, LineColor(kRed-2), RooFit::Name("total"));
  jpsi_model.plotOn(inclusive_jpsi_mass_frame, Components(jpsi_m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  jpsi_model.plotOn(inclusive_jpsi_mass_frame, Components(jpsi_m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  jpsi_model.plotOn(inclusive_jpsi_mass_frame, Components(jpsi_m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  jpsi_model.plotOn(inclusive_jpsi_mass_frame, Components(jpsi_m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));
  
  inclusive_jpsi_mass_frame->GetXaxis()->SetTitleOffset(.9);

  RooPlot* dummy_frame_jpsi = onia_mass->frame(Title("dummy frame to extract residuals"), Bins(num_mass_bins));
  jpsi_data->plotOn(dummy_frame_jpsi); 
  jpsi_model.plotOn(dummy_frame_jpsi);

  RooHist* h_residuals_mass_jpsi = dummy_frame_jpsi->pullHist();
  RooPlot* frame_residuals_mass_jpsi = onia_mass->frame(Title("Residual Distribution #mu^{+}#mu^{-} mass"));
  frame_residuals_mass_jpsi->GetYaxis()->SetTitle("(fit - data)/#sigma");
  frame_residuals_mass_jpsi->GetYaxis()->SetTitleSize(.16);
  frame_residuals_mass_jpsi->GetYaxis()->SetTitleOffset(.2);
  frame_residuals_mass_jpsi->addPlotable(h_residuals_mass_jpsi, "P");

  // Creates and fills canvas with plot and info 
  TCanvas *canvas1_jpsi = new TCanvas("jpsi_MassFit_", "jpsi_MassFit_", 900, 900);
  canvas1_jpsi->cd();
  if(draw_residuals) {
    TPad *pad1_jpsi = new TPad("pad1_jpsi", "The pad 80% of the height",0.0,0.05,1.0,1.0,21);
    TPad *pad2_jpsi = new TPad("pad2_jpsi", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);
    pad1_jpsi->Draw(); 
    pad2_jpsi->Draw();
    pad1_jpsi->SetFillColor(0); 
    pad2_jpsi->SetFillColor(0);
    pad2_jpsi->cd();
    frame_residuals_mass_jpsi->Draw(); 
    f_straighline->Draw("same");
    pad1_jpsi->cd();
    inclusive_jpsi_mass_frame->Draw();
  }
  if (!draw_residuals) {
    canvas1_jpsi->Close();
  }

  //splot_tex_pt.DrawLatex(2.87, 30, "#bf{#it{CMS}} Preliminary");
  //splot_tex_pt.DrawLatex(2.87, 25, "#sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
 
  //inclusive_jpsi_mass_frame->SetTitle("J/#psi Mass [GeV]");

  Double_t xl1_l3=.60, yl1_l3=0.55, xl2_l3=xl1_l3+.3, yl2_l3=yl1_l3+.325;
  TLegend *leg3 = new TLegend(xl1_l3,yl1_l3,xl2_l3,yl2_l3);
  leg3->SetFillColor(kWhite);
  leg3->AddEntry(inclusive_jpsi_mass_frame->findObject("total"),"total","l");
  leg3->AddEntry(inclusive_jpsi_mass_frame->findObject("prompt j/psi"),"prompt j/psi","l");
  leg3->AddEntry(inclusive_jpsi_mass_frame->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg3->AddEntry(inclusive_jpsi_mass_frame->findObject("prompt continuum"),"prompt continuum","l");
  leg3->AddEntry(inclusive_jpsi_mass_frame->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg3->SetShadowColor(0);
  //leg3->Draw();
  inclusive_jpsi_mass_frame->SetTitle("Inclusive J/#psi Mass");
  //

  //----------- Inclusive Jpsi Tau Fit ---------------
  //Plots data on to frame
  RooPlot* inclusive_jpsi_time_frame = onia_tau->frame(Bins(num_time_bins));
  jpsi_data->plotOn(inclusive_jpsi_time_frame, Name("jpsi_time_data"));

  //Plots full model, prompt and non-prompt models to frame
  jpsi_model.plotOn(inclusive_jpsi_time_frame, LineColor(kRed-2), RooFit::Name("total"));
  jpsi_model.plotOn(inclusive_jpsi_time_frame, Components(jpsi_m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  jpsi_model.plotOn(inclusive_jpsi_time_frame, Components(jpsi_m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  jpsi_model.plotOn(inclusive_jpsi_time_frame, Components(jpsi_m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  jpsi_model.plotOn(inclusive_jpsi_time_frame, Components(jpsi_m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));

  inclusive_jpsi_time_frame->GetXaxis()->SetTitleOffset(.9);

  RooPlot* dummy_tframe_jpsi = onia_tau->frame(Title("dummy frame to extract residuals"), Bins(num_time_bins));
  jpsi_data->plotOn(dummy_tframe_jpsi); 
  jpsi_model.plotOn(dummy_tframe_jpsi);

  RooHist* h_residuals_time_jpsi = dummy_tframe_jpsi->pullHist();
  RooPlot* frame_residuals_time_jpsi = onia_tau->frame(Title("Residual Distribution #mu^{+}#mu^{-} time"));
  frame_residuals_time_jpsi->GetYaxis()->SetTitle("(fit - data)/#sigma");
  frame_residuals_time_jpsi->GetYaxis()->SetTitleSize(.16);
  frame_residuals_time_jpsi->GetYaxis()->SetTitleOffset(.2);
  frame_residuals_time_jpsi->addPlotable(h_residuals_time_jpsi, "P");

  // Creates and fills canvas with plot and info 
  TCanvas *canvas2_jpsi = new TCanvas("jpsi_TimeFit_", "jpsi_TimeFit_", 900, 900);
  canvas2_jpsi->cd();
  if(draw_residuals) {
    TPad *tpad1_jpsi = new TPad("tpad1_jpsi", "The pad 80% of the height",0.0,0.05,1.0,1.0,21); tpad1_jpsi->SetLogy(1);
    TPad *tpad2_jpsi = new TPad("tpad2_jpsi", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);
    tpad1_jpsi->Draw();
    tpad2_jpsi->Draw();
    tpad1_jpsi->SetFillColor(0);
    tpad2_jpsi->SetFillColor(0);
    tpad2_jpsi->cd();
    frame_residuals_time_jpsi->Draw();
    f_straighline->Draw("same");
    tpad1_jpsi->cd();
    inclusive_jpsi_time_frame->Draw();
  }
  if (!draw_residuals) {
    canvas2_jpsi->Close();
  }

  //splot_tex_pt.DrawLatex(2, 115, "#bf{#it{CMS}} Preliminary");
  //splot_tex_pt.DrawLatex(2, 80, "#sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
 
  Double_t xl1_l4=.65, yl1_l4=0.55, xl2_l4=xl1_l4+.2, yl2_l4=yl1_l4+.225;
  TLegend *leg4 = new TLegend(xl1_l4,yl1_l4,xl2_l4,yl2_l4);
  leg4->SetFillColor(kWhite);
  leg4->AddEntry(inclusive_jpsi_time_frame->findObject("total"),"total","l");
  leg4->AddEntry(inclusive_jpsi_time_frame->findObject("prompt j/psi"),"prompt j/psi","l");
  leg4->AddEntry(inclusive_jpsi_time_frame->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg4->AddEntry(inclusive_jpsi_time_frame->findObject("prompt continuum"),"prompt continuum","l");
  leg4->AddEntry(inclusive_jpsi_time_frame->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg4->SetShadowColor(0);
  //leg4->Draw();
  //inclusive_jpsi_time_frame->SetTitle("J/#psi Mass [GeV]");
  inclusive_jpsi_time_frame->SetTitle("Inclusive J/#psi Lifetime");
  canvas2_jpsi->SetLogy(1);

  //
  //RooPlot* dummy_t = onia_tau->frame(Bins(num_time_bins));
  //jpsi_data->plotOn(dummy_t);
  //jpsi_model.plotOn(dummy_t);

  TCanvas *c_dummy = new TCanvas("c_dummy", "c_dummy", 1200, 600); 
  c_dummy->Divide(2,1);
  c_dummy->cd(1);
  //c_dummy->cd(1)->SetLogy(1);
  inclusive_jpsi_mass_frame->Draw();
  leg3->Draw();
  c_dummy->cd(2); 
  inclusive_jpsi_time_frame->Draw(); 
  leg4->Draw();
  c_dummy->cd(2)->SetLogy(1);

  std::string c_dummy_image_name = image_dir;
  c_dummy_image_name.append("inclusive_fit_pt_bin_");
  c_dummy_image_name.append(PT_SLICE_STRING);
  c_dummy_image_name.append(".png");
  //std::cout <<  c_dummy_image_name << std::endl;
  c_dummy->Print(c_dummy_image_name.c_str() , "png");
  c_dummy->Close();


  //TCanvas *c_dummy = new TCanvas("c_dummy", "c_dummy", 1200, 600); c_dummy->Divide(2,1);
  //c_dummy->cd(1); 
  //dummy_m->Draw();
  //c_dummy->cd(2); 
  //dummy_t->Draw(); 
  //c_dummy->cd(2)->SetLogy(1);

  //Then would do:

  double zjpsi_gaussmean_value = jpsi_gaussmean.getVal();
  double zjpsi_gausssigma_value     = jpsi_gausssigma.getVal();
  double zjpsi_gausssigma_value_err = jpsi_gausssigma.getError();
  double zjpsi_decay_time_value     = jpsi_decay_time.getVal();
  double zjpsi_decay_time_value_err = jpsi_decay_time.getError();
  double zjpsi_decay_time_continuum_value     = jpsi_decay_time_continuum.getVal();
  double zjpsi_decay_time_continuum_value_err = jpsi_decay_time_continuum.getError();
  double zjpsi_gauss_prompt_mean_value     = jpsi_gauss_prompt_mean.getVal();
  double zjpsi_gauss_prompt_mean_value_err = jpsi_gauss_prompt_mean.getError();

  double zjpsi_gauss_prompt_sigma_value     = jpsi_gauss_prompt_sigma.getVal();
  double zjpsi_gauss_prompt_sigma_value_err = jpsi_gauss_prompt_sigma.getError();
  double zjpsi_gauss_prompt_sigma_2_value     = jpsi_gauss_prompt_sigma_2.getVal();
  double zjpsi_gauss_prompt_sigma_2_value_err = jpsi_gauss_prompt_sigma_2.getError();

  //double zjpsi_gauss_prompt_sigma_value     = jpsi_gauss_prompt_sigma.getVal();
  //double zjpsi_gauss_prompt_sigma_value_err = jpsi_gauss_prompt_sigma.getError();
  //double zjpsi_gauss_prompt_sigma_2_value     = jpsi_gauss_prompt_sigma_2.getVal();
  //double zjpsi_gauss_prompt_sigma_2_value_err = jpsi_gauss_prompt_sigma_2.getError();
  double zjpsi_frac_prompt_sharp_value     = jpsi_frac_prompt_sharp.getVal();
  double zjpsi_frac_prompt_sharp_value_err = jpsi_frac_prompt_sharp.getError();

  double zjpsi_sigma_value     = jpsi_sigma.getVal();
  double zjpsi_sigma_value_err = jpsi_sigma.getError();
  double zjpsi_alpha_value     = jpsi_alpha.getVal();
  double zjpsi_alpha_value_err = jpsi_alpha.getError();
  double zjpsi_n_value     = jpsi_n.getVal();
  double zjpsi_n_value_err = jpsi_n.getError();
  double zjpsi_dimuon_slope_value     = jpsi_dimuon_slope.getVal();
  double zjpsi_dimuon_slope_value_err = jpsi_dimuon_slope.getError();
  double zjpsi_np_dimuon_slope_value     = jpsi_np_dimuon_slope.getVal();
  double zjpsi_np_dimuon_slope_value_err = jpsi_np_dimuon_slope.getError();
  double zjpsi_dimuon_mean_value     = jpsi_dimuon_mean.getVal();
  double zjpsi_dimuon_mean_value_err = jpsi_dimuon_mean.getError();
  double zjpsi_dimuon_sigma_value     = jpsi_dimuon_sigma.getVal();
  double zjpsi_dimuon_sigma_value_err = jpsi_dimuon_sigma.getError();
  double zjpsi_frac_cball_value     = jpsi_frac_cball.getVal();
  double zjpsi_frac_cball_value_err = jpsi_frac_cball.getError();

  // ZJPsi Unbinned fit start here

  // *****************************************************************
  //  ZJpsi 
  // *****************************************************************

  // ZJpsi Unbinned fit start here

  //RooRealVar zjpsi_gaussmean("zjpsi_gaussmean", "Mean of the smearing zjpsi Gaussian", zjpsi_gaussmean_value);
  //RooRealVar zjpsi_gausssigma("zjpsi_gausssigma", "Width of the smearing zjpsi Gaussian", zjpsi_gausssigma_value);
  //RooGaussModel zjpsi_smear_gauss_model("zjpsi_smear_gauss", "Gaussian used to smear the zjpsi Exponential Decay", *onia_tau, zjpsi_gaussmean, zjpsi_gausssigma);
  //RooRealVar zjpsi_decay_time("zjpsi_decay_time", "zjpsi_Tau", zjpsi_decay_time_value, zjpsi_decay_time_value-zjpsi_decay_time_value_err, zjpsi_decay_time_value+zjpsi_decay_time_value_err);
  //RooDecay zjpsi_decay_exp("zjpsi_decay_exp", "zjpsi_Exponential Decay", *onia_tau,zjpsi_decay_time,zjpsi_smear_gauss_model,RooDecay::SingleSided);

  //RooRealVar zjpsi_gauss_prompt_mean("zjpsi_gauss_prompt_mean", "Mean of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_mean_value, zjpsi_gauss_prompt_mean_value-zjpsi_gauss_prompt_mean_value_err, zjpsi_gauss_prompt_mean_value+zjpsi_gauss_prompt_mean_value_err);
  //RooRealVar zjpsi_gauss_prompt_sigma("zjpsi_gauss_prompt_sigma", "Width of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_sigma_value, zjpsi_gauss_prompt_sigma_value-zjpsi_gauss_prompt_sigma_value_err, zjpsi_gauss_prompt_sigma_value+zjpsi_gauss_prompt_sigma_value_err);
  //RooGaussian zjpsi_prompt_gauss("zjpsi_prompt_gauss", "Gaussian of the Prompt zjpsi_Peak", *onia_tau, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma);

  //RooRealVar zjpsi_gauss_prompt_sigma_2("zjpsi_gauss_prompt_sigma_2", "Width of the Prompt zjpsi_Gaussian_2", zjpsi_gauss_prompt_sigma_2_value, zjpsi_gauss_prompt_sigma_2_value-zjpsi_gauss_prompt_sigma_2_value_err, zjpsi_gauss_prompt_sigma_2_value+zjpsi_gauss_prompt_sigma_2_value_err);
  //RooGaussian zjpsi_prompt_gauss_2("zjpsi_prompt_gauss_2", "Gaussian_2 of the Prompt zjpsi_Peak", *onia_tau, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma_2);

  //RooRealVar zjpsi_frac_prompt_sharp("zjpsi_frac_prompt_sharp", "zjpsi_frac_prompt_sharp" , zjpsi_frac_prompt_sharp_value, zjpsi_frac_prompt_sharp_value-zjpsi_frac_prompt_sharp_value_err, zjpsi_frac_prompt_sharp_value+zjpsi_frac_prompt_sharp_value_err);
  //RooAddPdf onia_tau_gauss_sum_fitpdf("onia_tau_gauss_sum_fitpdf", "onia_tau_gauss_sum_fitpdf", RooArgList(zjpsi_prompt_gauss,zjpsi_prompt_gauss_2), RooArgList(zjpsi_frac_prompt_sharp));

  //RooRealVar zjpsi_sigma("zjpsi_sigma", "zjpsi_sigma", zjpsi_sigma_value, zjpsi_sigma_value-zjpsi_sigma_value_err, zjpsi_sigma_value+zjpsi_sigma_value_err);
  //RooRealVar zjpsi_alpha("zjpsi_alpha", "zjpsi_alpha", zjpsi_alpha_value, zjpsi_alpha_value-zjpsi_alpha_value_err, zjpsi_alpha_value+zjpsi_alpha_value_err);
  //RooRealVar zjpsi_n("zjpsi_n", "zjpsi_n", zjpsi_n_value, zjpsi_n_value-zjpsi_n_value_err, zjpsi_n_value+zjpsi_n_value_err);
  //RooRealVar zjpsi_dimuon_mean("zjpsi_dimuon_mean", "zjpsi_dimuon_mean", zjpsi_dimuon_mean_value, zjpsi_dimuon_mean_value-zjpsi_dimuon_mean_value_err, zjpsi_dimuon_mean_value+zjpsi_dimuon_mean_value_err);
  //RooCBShape zjpsi_crystal_ball ("zjpsi_crystal_ball", "zjpsi_crystal_ball", *onia_mass, zjpsi_dimuon_mean, zjpsi_sigma, zjpsi_alpha, zjpsi_n );
  //RooRealVar zjpsi_dimuon_sigma("zjpsi_dimuon_sigma", "zjpsi_dimuon_sigma", zjpsi_dimuon_sigma_value, zjpsi_dimuon_sigma_value-zjpsi_dimuon_sigma_value_err, zjpsi_dimuon_sigma_value+zjpsi_dimuon_sigma_value_err);
  //RooGaussian zjpsi_dimuon_gauss ("zjpsi_dimuon_gauss", "zjpsi_dimuon_gauss", *onia_mass, zjpsi_dimuon_mean, zjpsi_dimuon_sigma);
  //RooRealVar zjpsi_frac_cball("zjpsi_frac_cball", "zjpsi_frac_cball" , zjpsi_frac_cball_value, zjpsi_frac_cball_value-zjpsi_frac_cball_value_err, zjpsi_frac_cball_value+zjpsi_frac_cball_value_err);
  //RooAddPdf zjpsi_mass_signal("zjpsi_mass_signal", "zjpsi_mass_signal", RooArgList(zjpsi_crystal_ball, zjpsi_dimuon_gauss), RooArgList(zjpsi_frac_cball));

  //RooRealVar zjpsi_dimuon_slope("zjpsi_dimuon_slope", "zjpsi_dimuon_slope", zjpsi_dimuon_slope_value, zjpsi_dimuon_slope_value-zjpsi_dimuon_slope_value_err, zjpsi_dimuon_slope_value+zjpsi_dimuon_slope_value_err);
  //RooExponential zjpsi_dimuon_bg_exponential("zjpsi_dimuon_bg_exponential", "zjpsi_dimuon_bg_exponential", *onia_mass, zjpsi_dimuon_slope);

  //RooRealVar zjpsi_np_dimuon_slope("zjpsi_np_dimuon_slope", "zjpsi_np_dimuon_slope", zjpsi_np_dimuon_slope_value, zjpsi_np_dimuon_slope_value-zjpsi_np_dimuon_slope_value_err, zjpsi_np_dimuon_slope_value+zjpsi_np_dimuon_slope_value_err);
  //RooExponential zjpsi_np_dimuon_bg_exponential("zjpsi_np_dimuon_bg_exponential", "zjpsi_np_dimuon_bg_exponential", *onia_mass, zjpsi_np_dimuon_slope);

  //TODO testing effect of not letting parameters vary ambivalent, want to see if the nan errors disappear
  RooRealVar zjpsi_gaussmean("zjpsi_gaussmean", "Mean of the smearing zjpsi Gaussian", zjpsi_gaussmean_value);
  RooRealVar zjpsi_gausssigma("zjpsi_gausssigma", "Width of the smearing zjpsi Gaussian", zjpsi_gausssigma_value);
  RooGaussModel zjpsi_smear_gauss_model("zjpsi_smear_gauss", "Gaussian used to smear the zjpsi Exponential Decay", *onia_tau, zjpsi_gaussmean, zjpsi_gausssigma);
  RooRealVar zjpsi_decay_time("zjpsi_decay_time", "zjpsi_Tau", zjpsi_decay_time_value);
  RooDecay zjpsi_decay_exp("zjpsi_decay_exp", "zjpsi_Exponential Decay", *onia_tau,zjpsi_decay_time,zjpsi_smear_gauss_model,RooDecay::SingleSided);

  RooRealVar zjpsi_decay_time_continuum("zjpsi_decay_time_continuum", "zjpsi_Tau_continuum", zjpsi_decay_time_continuum_value);
  RooDecay zjpsi_decay_exp_continuum("zjpsi_decay_exp_continuum", "zjpsi_Exponential Decay_continuum", *onia_tau,zjpsi_decay_time_continuum,zjpsi_smear_gauss_model,RooDecay::SingleSided);

  RooRealVar zjpsi_gauss_prompt_mean("zjpsi_gauss_prompt_mean", "Mean of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_mean_value);
  RooRealVar zjpsi_gauss_prompt_sigma("zjpsi_gauss_prompt_sigma", "Width of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_sigma_value);
  RooGaussian zjpsi_prompt_gauss("zjpsi_prompt_gauss", "Gaussian of the Prompt zjpsi_Peak", *onia_tau, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma);

  RooRealVar zjpsi_gauss_prompt_sigma_2("zjpsi_gauss_prompt_sigma_2", "Width of the Prompt zjpsi_Gaussian_2", zjpsi_gauss_prompt_sigma_2_value);
  RooGaussian zjpsi_prompt_gauss_2("zjpsi_prompt_gauss_2", "Gaussian_2 of the Prompt zjpsi_Peak", *onia_tau, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma_2);

  RooRealVar zjpsi_frac_prompt_sharp("zjpsi_frac_prompt_sharp", "zjpsi_frac_prompt_sharp" , zjpsi_frac_prompt_sharp_value);
  RooAddPdf onia_tau_gauss_sum_fitpdf("onia_tau_gauss_sum_fitpdf", "onia_tau_gauss_sum_fitpdf", RooArgList(zjpsi_prompt_gauss,zjpsi_prompt_gauss_2), RooArgList(zjpsi_frac_prompt_sharp));

  RooRealVar zjpsi_sigma("zjpsi_sigma", "zjpsi_sigma", zjpsi_sigma_value);
  RooRealVar zjpsi_alpha("zjpsi_alpha", "zjpsi_alpha", zjpsi_alpha_value);
  RooRealVar zjpsi_n("zjpsi_n", "zjpsi_n", zjpsi_n_value);
  RooRealVar zjpsi_dimuon_mean("zjpsi_dimuon_mean", "zjpsi_dimuon_mean", zjpsi_dimuon_mean_value);
  RooCBShape zjpsi_crystal_ball ("zjpsi_crystal_ball", "zjpsi_crystal_ball", *onia_mass, zjpsi_dimuon_mean, zjpsi_sigma, zjpsi_alpha, zjpsi_n );
  RooRealVar zjpsi_dimuon_sigma("zjpsi_dimuon_sigma", "zjpsi_dimuon_sigma", zjpsi_dimuon_sigma_value);
  RooGaussian zjpsi_dimuon_gauss ("zjpsi_dimuon_gauss", "zjpsi_dimuon_gauss", *onia_mass, zjpsi_dimuon_mean, zjpsi_dimuon_sigma);
  RooRealVar zjpsi_frac_cball("zjpsi_frac_cball", "zjpsi_frac_cball" , zjpsi_frac_cball_value);
  RooAddPdf zjpsi_mass_signal("zjpsi_mass_signal", "zjpsi_mass_signal", RooArgList(zjpsi_crystal_ball, zjpsi_dimuon_gauss), RooArgList(zjpsi_frac_cball));

  RooRealVar zjpsi_dimuon_slope("zjpsi_dimuon_slope", "zjpsi_dimuon_slope", zjpsi_dimuon_slope_value);
  RooExponential zjpsi_dimuon_bg_exponential("zjpsi_dimuon_bg_exponential", "zjpsi_dimuon_bg_exponential", *onia_mass, zjpsi_dimuon_slope);

  RooRealVar zjpsi_np_dimuon_slope("zjpsi_np_dimuon_slope", "zjpsi_np_dimuon_slope", zjpsi_np_dimuon_slope_value);
  RooExponential zjpsi_np_dimuon_bg_exponential("zjpsi_np_dimuon_bg_exponential", "zjpsi_np_dimuon_bg_exponential", *onia_mass, zjpsi_np_dimuon_slope);

//  RooRealVar zjpsi_m_sig_tau_sig_frac("zjpsi_m_sig_tau_sig_frac", "zjpsi_m_sig_tau_sig_frac", 0.2, 0.0, 1.0);
//  RooRealVar zjpsi_m_sig_tau_bg_frac("zjpsi_m_sig_tau_bg_frac", "zjpsi_m_sig_tau_bg_frac", 0.6, 0.0, 1.0);
//  RooRealVar zjpsi_m_bg_tau_sig_frac("zjpsi_m_bg_tau_sig_frac", "zjpsi_m_bg_tau_sig_frac", 0.15, 0.0, 1.0);

  // PDFs
  RooProdPdf zjpsi_m_sig_tau_sig("zjpsi_m_sig_tau_sig", "zjpsi_m_sig_tau_sig", RooArgList(zjpsi_mass_signal, onia_tau_gauss_sum_fitpdf ));
  RooProdPdf zjpsi_m_sig_tau_bg ("zjpsi_m_sig_tau_bg",  "zjpsi_m_sig_tau_bg",  RooArgList(zjpsi_mass_signal, zjpsi_decay_exp ));
  //RooProdPdf zjpsi_m_bg_tau_bg  ("zjpsi_m_bg_tau_bg",   "zjpsi_m_bg_tau_bg",   RooArgList(zjpsi_np_dimuon_bg_exponential, zjpsi_decay_exp ));
  RooProdPdf zjpsi_m_bg_tau_bg  ("zjpsi_m_bg_tau_bg",   "zjpsi_m_bg_tau_bg",   RooArgList(zjpsi_np_dimuon_bg_exponential, zjpsi_decay_exp_continuum ));
  RooProdPdf zjpsi_m_bg_tau_sig ("zjpsi_m_bg_tau_sig",  "zjpsi_m_bg_tau_sig",  RooArgList(zjpsi_dimuon_bg_exponential, onia_tau_gauss_sum_fitpdf ));

  // yields
  RooRealVar Nzjpsi_m_sig_tau_sig("Nzjpsi_m_sig_tau_sig", "Nzjpsi_m_sig_tau_sig", 10, 0, 2000);
  RooRealVar Nzjpsi_m_sig_tau_bg ("Nzjpsi_m_sig_tau_bg",  "Nzjpsi_m_sig_tau_bg",  10, 0, 2000);
  RooRealVar Nzjpsi_m_bg_tau_bg  ("Nzjpsi_m_bg_tau_bg",   "Nzjpsi_m_bg_tau_bg",   10, 0, 2000);
  RooRealVar Nzjpsi_m_bg_tau_sig ("Nzjpsi_m_bg_tau_sig",  "Nzjpsi_m_bg_tau_sig",  10, 0, 2000);

  //RooRealVar Nzjpsi_m_sig_tau_sig("Nzjpsi_m_sig_tau_sig", "Nzjpsi_m_sig_tau_sig", 0.3, 0, 2000);
  //RooRealVar Nzjpsi_m_sig_tau_bg ("Nzjpsi_m_sig_tau_bg",  "Nzjpsi_m_sig_tau_bg",  7, 0, 2000);
  //RooRealVar Nzjpsi_m_bg_tau_bg  ("Nzjpsi_m_bg_tau_bg",   "Nzjpsi_m_bg_tau_bg",   0, 0, 2000);
  //RooRealVar Nzjpsi_m_bg_tau_sig ("Nzjpsi_m_bg_tau_sig",  "Nzjpsi_m_bg_tau_sig",  10, 0, 2000);

  //RooRealVar Nzjpsi_m_sig_tau_sig("Nzjpsi_m_sig_tau_sig", "Nzjpsi_m_sig_tau_sig", 10, -2, 200);
  //RooRealVar Nzjpsi_m_sig_tau_bg ("Nzjpsi_m_sig_tau_bg",  "Nzjpsi_m_sig_tau_bg",  10, -2, 200);
  //RooRealVar Nzjpsi_m_bg_tau_bg  ("Nzjpsi_m_bg_tau_bg",   "Nzjpsi_m_bg_tau_bg",   10, -2, 200);
  //RooRealVar Nzjpsi_m_bg_tau_sig ("Nzjpsi_m_bg_tau_sig",  "Nzjpsi_m_bg_tau_sig",  10, -2, 200);

  // extended PDFs
  RooExtendPdf ezjpsi_m_sig_tau_sig("ezjpsi_m_sig_tau_sig", "ezjpsi_m_sig_tau_sig", zjpsi_m_sig_tau_sig, Nzjpsi_m_sig_tau_sig); 
  RooExtendPdf ezjpsi_m_sig_tau_bg ("ezjpsi_m_sig_tau_bg",  "ezjpsi_m_sig_tau_bg",  zjpsi_m_sig_tau_bg , Nzjpsi_m_sig_tau_bg );
  RooExtendPdf ezjpsi_m_bg_tau_bg  ("ezjpsi_m_bg_tau_bg",   "ezjpsi_m_bg_tau_bg",   zjpsi_m_bg_tau_bg  , Nzjpsi_m_bg_tau_bg  );
  RooExtendPdf ezjpsi_m_bg_tau_sig ("ezjpsi_m_bg_tau_sig",  "ezjpsi_m_bg_tau_sig",  zjpsi_m_bg_tau_sig , Nzjpsi_m_bg_tau_sig );

  RooAddPdf zjpsi_model("zjpsi_model", "zjpsi_model", RooArgList(ezjpsi_m_sig_tau_sig, ezjpsi_m_sig_tau_bg, ezjpsi_m_bg_tau_sig, ezjpsi_m_bg_tau_bg));

  //RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data, NumCPU(4, kTRUE), Verbose(false), PrintLevel(-1), Save());
  //RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data, NumCPU(4, kTRUE), Verbose(false), PrintLevel(-1), SumW2Error(kTRUE), Save());
  //RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data_weighted_testing, NumCPU(4, kTRUE), Extended(kTRUE) ,Verbose(false), PrintLevel(-1), SumW2Error(kTRUE), Save());
  //RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data_weighted_testing, NumCPU(1, kTRUE), Extended(kTRUE) , Verbose(false), PrintLevel(-1), SumW2Error(kTRUE), Save());
  //RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data_weighted_testing, NumCPU(4, kTRUE), Verbose(false), Minimizer(1), PrintLevel(-1), SumW2Error(kTRUE), Save());

  RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data_weighted_testing, NumCPU(1, kTRUE), Verbose(false), PrintLevel(-1), SumW2Error(kTRUE), Save());
  //RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data_weighted_testing, NumCPU(4, kTRUE), Verbose(false), InitialHesse(kTRUE), PrintLevel(-1), SumW2Error(kTRUE), Save());
  //RooFitResult *fr = zjpsi_model.fitTo(*zjpsi_data_weighted_testing, NumCPU(4, kTRUE), Verbose(false), Minos(kTRUE), PrintLevel(-1), SumW2Error(kTRUE), Save());
  fr->Print();

  if (do_pull_calculation) {
    int num_toys = 1000;
    RooMCStudy *toy_simulation = new RooMCStudy(zjpsi_model, RooArgSet(*onia_mass, *onia_tau), Binned(kFALSE), Silence(), Extended(), FitOptions(Save(kTRUE)), PrintEvalErrors(0), Minos(kTRUE));

    toy_simulation->generateAndFit(num_toys);

    RooPlot* sig_prompt_pull_frame        = toy_simulation->plotPull(Nzjpsi_m_sig_tau_sig, -4, 4, 16, kTRUE);
    RooPlot* mass_bkg_prompt_pull_frame    = toy_simulation->plotPull(Nzjpsi_m_bg_tau_sig,    -4, 4, 16, kTRUE);
    RooPlot* sig_nonprompt_pull_frame = toy_simulation->plotPull(Nzjpsi_m_sig_tau_bg, -4, 4, 16, kTRUE);
    RooPlot* mass_bkg_nonprompt_pull_frame = toy_simulation->plotPull(Nzjpsi_m_bg_tau_bg, -4, 4, 16, kTRUE);

    TCanvas* c2 = new TCanvas("c2","toymc", 800, 800); c2->Divide(2,2);
    c2->cd(1); 
    sig_prompt_pull_frame->Draw();
    sig_prompt_pull_frame->SetXTitle("Yield Signal Prompt Pull");
    c2->cd(2); 
    mass_bkg_prompt_pull_frame->Draw();
    mass_bkg_prompt_pull_frame->SetXTitle("Yield Mass Background Prompt Pull");
    c2->cd(3); 
    sig_nonprompt_pull_frame->Draw();
    sig_nonprompt_pull_frame->SetXTitle("Yield Signal Non Prompt Pull");
    c2->cd(4); 
    mass_bkg_nonprompt_pull_frame->Draw();
    mass_bkg_nonprompt_pull_frame->SetXTitle("Yield Mass Background Non Prompt Pull");
  }

  // *****************************
  // START OF JPSI PLOTTING
  // *****************************

  //----------- Mass Fit ---------------
  //Plots zjpsi_data on to frame
  RooPlot* mass_frame = onia_mass->frame(Bins(num_mass_bins));
  zjpsi_data->plotOn(mass_frame, Name("zjpsi_mass_data"));

  //Plots full model, prompt and non-prompt models to frame
  zjpsi_model.plotOn(mass_frame, LineColor(kRed-2), RooFit::Name("total"));
  zjpsi_model.plotOn(mass_frame, Components(zjpsi_m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  zjpsi_model.plotOn(mass_frame, Components(zjpsi_m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  zjpsi_model.plotOn(mass_frame, Components(zjpsi_m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  zjpsi_model.plotOn(mass_frame, Components(zjpsi_m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));
  
  mass_frame->GetXaxis()->SetTitleOffset(.9);

  RooPlot* dummy_frame_zjpsi = onia_mass->frame(Title("dummy frame to extract residuals"), Bins(num_mass_bins));
  zjpsi_data->plotOn(dummy_frame_zjpsi); 
  zjpsi_model.plotOn(dummy_frame_zjpsi);

  RooHist* h_residuals_mass_zjpsi = dummy_frame_zjpsi->pullHist();
  RooPlot* frame_residuals_mass_zjpsi = onia_mass->frame(Title("Residual Distribution #mu^{+}#mu^{-} mass"));
  frame_residuals_mass_zjpsi->GetYaxis()->SetTitle("(fit - data)/#sigma");
  frame_residuals_mass_zjpsi->GetYaxis()->SetTitleSize(.16);
  frame_residuals_mass_zjpsi->GetYaxis()->SetTitleOffset(.2);
  frame_residuals_mass_zjpsi->addPlotable(h_residuals_mass_zjpsi, "P");

  // Creates and fills canvas with plot and info 
  TCanvas *canvas1 = new TCanvas("ZJPsi_MassFit_", "ZJPsi_MassFit_", 900, 900);
  canvas1->cd();
  if(draw_residuals) {
    TPad *pad1_jpsi = new TPad("pad1_jpsi", "The pad 80% of the height",0.0,0.05,1.0,1.0,21);
    TPad *pad2_jpsi = new TPad("pad2_jpsi", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);
    pad1_jpsi->Draw(); 
    pad2_jpsi->Draw();
    pad1_jpsi->SetFillColor(0); 
    pad2_jpsi->SetFillColor(0);
    pad2_jpsi->cd();
    frame_residuals_mass_zjpsi->Draw(); 
    f_straighline->Draw("same");
    pad1_jpsi->cd();
    mass_frame->Draw();
  }
  if (!draw_residuals) {
    canvas1->Close();
  }

  //splot_tex_pt.DrawLatex(2.87, 30, "#bf{#it{CMS}} Preliminary");
  //splot_tex_pt.DrawLatex(2.87, 25, "#sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
 
  //mass_frame->SetTitle("J/#psi Mass [GeV]");

  Double_t xl1_l2=.60, yl1_l2=0.55, xl2_l2=xl1_l2+.3, yl2_l2=yl1_l2+.325;
  TLegend *leg2 = new TLegend(xl1_l2,yl1_l2,xl2_l2,yl2_l2);
  leg2->SetFillColor(kWhite);
  leg2->AddEntry(mass_frame->findObject("total"),"total","l");
  leg2->AddEntry(mass_frame->findObject("prompt j/psi"),"prompt j/psi","l");
  leg2->AddEntry(mass_frame->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg2->AddEntry(mass_frame->findObject("prompt continuum"),"prompt continuum","l");
  leg2->AddEntry(mass_frame->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg2->SetShadowColor(0);
  //leg2->Draw();
  mass_frame->SetTitle("Z + J/#psi Mass");

  //----------- Tau Fit ---------------
  //Plots data on to frame
  RooPlot* time_frame = onia_tau->frame(Bins(num_time_bins));
  zjpsi_data->plotOn(time_frame, Name("zjpsi_time_data"));

  //Plots full model, prompt and non-prompt models to frame
  zjpsi_model.plotOn(time_frame, LineColor(kRed-2), RooFit::Name("total"));
  zjpsi_model.plotOn(time_frame, Components(zjpsi_m_sig_tau_sig), LineColor(kBlue-2), RooFit::Name("prompt j/psi"));
  zjpsi_model.plotOn(time_frame, Components(zjpsi_m_sig_tau_bg), LineColor(kMagenta-2), RooFit::Name("non-prompt j/psi"));
  zjpsi_model.plotOn(time_frame, Components(zjpsi_m_bg_tau_sig), LineColor(kCyan-2), RooFit::Name("prompt continuum"));
  zjpsi_model.plotOn(time_frame, Components(zjpsi_m_bg_tau_bg), LineColor(kGreen-2), RooFit::Name("non-prompt continuum"));

  time_frame->GetXaxis()->SetTitleOffset(.9);

  RooPlot* dummy_tframe_zjpsi = onia_tau->frame(Title("dummy frame to extract residuals"), Bins(num_time_bins));
  zjpsi_data->plotOn(dummy_tframe_zjpsi); 
  zjpsi_model.plotOn(dummy_tframe_zjpsi);

  RooHist* h_residuals_time_zjpsi = dummy_tframe_zjpsi->pullHist();
  RooPlot* frame_residuals_time_zjpsi = onia_tau->frame(Title("Residual Distribution #mu^{+}#mu^{-} time"));
  frame_residuals_time_zjpsi->GetYaxis()->SetTitle("(fit - data)/#sigma");
  frame_residuals_time_zjpsi->GetYaxis()->SetTitleSize(.16);
  frame_residuals_time_zjpsi->GetYaxis()->SetTitleOffset(.2);
  frame_residuals_time_zjpsi->addPlotable(h_residuals_time_zjpsi, "P");

  // Creates and fills canvas with plot and info 
  TCanvas *canvas2 = new TCanvas("ZJpsi_TimeFit_", "ZJpsi_TimeFit_", 900, 900);
  canvas2->cd();
  if(draw_residuals) {
    TPad *tpad1_jpsi = new TPad("tpad1_jpsi", "The pad 80% of the height",0.0,0.05,1.0,1.0,21); tpad1_jpsi->SetLogy(1);
    TPad *tpad2_jpsi = new TPad("tpad2_jpsi", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);
    tpad1_jpsi->Draw();
    tpad2_jpsi->Draw();
    tpad1_jpsi->SetFillColor(0);
    tpad2_jpsi->SetFillColor(0);
    tpad2_jpsi->cd();
    frame_residuals_time_zjpsi->Draw();
    f_straighline->Draw("same");
    tpad1_jpsi->cd();
    time_frame->SetTitle("Z + J/#psi");
    time_frame->Draw();
  }

  //splot_tex_pt.DrawLatex(2, 115, "#bf{#it{CMS}} Preliminary");
  //splot_tex_pt.DrawLatex(2, 80, "#sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
 
  Double_t xl1=.60, yl1=0.70, xl2=xl1+.25, yl2=yl1+.225;
  TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
  leg->SetFillColor(kWhite);
  leg->AddEntry(time_frame->findObject("total"),"total","l");
  leg->AddEntry(time_frame->findObject("prompt j/psi"),"prompt j/psi","l");
  leg->AddEntry(time_frame->findObject("non-prompt j/psi"),"non-prompt j/psi","l");
  leg->AddEntry(time_frame->findObject("prompt continuum"),"prompt continuum","l");
  leg->AddEntry(time_frame->findObject("non-prompt continuum"),"non-prompt continuum","l");
  leg->SetShadowColor(0);
  //leg->Draw();
  //time_frame->SetTitle("J/#psi Mass [GeV]");
  time_frame->SetTitle("Z + J/#psi Lifetime");
  //canvas2->SetLogy(1);
  if (!draw_residuals) {
    canvas2->Close();
  }

  TCanvas *c_dummy_zjpsi = new TCanvas("c_dummy_zjpsi", "c_dummy_zjpsi", 1200, 600);
  c_dummy_zjpsi->Divide(2,1);
  c_dummy_zjpsi->cd(1); 
  //c_dummy_zjpsi->cd(1)->SetLogy(1); 
  mass_frame->Draw();
  leg2->Draw();
  c_dummy_zjpsi->cd(2); 
  time_frame->Draw(); 
  leg->Draw();
  //c_dummy_zjpsi->cd(2)->SetLogy(1);

  std::string c_dummy_zjpsi_image_name = image_dir;
  c_dummy_zjpsi_image_name.append("zjpsi_fit_pt_bin_");
  c_dummy_zjpsi_image_name.append(PT_SLICE_STRING);
  c_dummy_zjpsi_image_name.append(".png");
  //std::cout <<  c_dummy_zjpsi_image_name << std::endl;
  c_dummy_zjpsi->Print(c_dummy_zjpsi_image_name.c_str() , "png");
  c_dummy_zjpsi->Close();

  std::vector<double> zjpsi_info ;
  double zjpsi_prompt_events = Nzjpsi_m_sig_tau_sig.getVal();
  double zjpsi_prompt_events_error = Nzjpsi_m_sig_tau_sig.getError();
  double zjpsi_nonprompt_events = Nzjpsi_m_sig_tau_bg.getVal();
  double zjpsi_nonprompt_events_error = Nzjpsi_m_sig_tau_bg.getError();

  zjpsi_info.push_back(zjpsi_prompt_events);
  zjpsi_info.push_back(zjpsi_prompt_events_error);
  zjpsi_info.push_back(zjpsi_nonprompt_events);
  zjpsi_info.push_back(zjpsi_nonprompt_events_error);
  return zjpsi_info;


//// *******************************************************
//// *******************************************************
//// Splot
//// *******************************************************
//// *******************************************************



//  RooWorkspace* ws = new RooWorkspace("myWS");
//  ws->import(zjpsi_model);
//  ws->import(*zjpsi_data, Rename("sPlot_zjpsi_data"));
//  
//  RooAbsPdf* sPlotmodel = ws->pdf("zjpsi_model");
//  RooRealVar* sig_prompt_yield = ws->var("Nzjpsi_m_sig_tau_sig");
//  RooRealVar* sig_nonprompt_yield = ws->var("Nzjpsi_m_sig_tau_bg");
//  RooRealVar* bkg_prompt_yield = ws->var("Nzjpsi_m_bg_tau_sig");
//  RooRealVar* bkg_nonprompt_yield = ws->var("Nzjpsi_m_bg_tau_bg");
//  
//  sig_prompt_yield->setConstant();
//  sig_nonprompt_yield->setConstant();
//  bkg_prompt_yield->setConstant();
//  bkg_nonprompt_yield->setConstant();
//  
//  RooStats::SPlot* sPlot_zjpsi_data = new RooStats::SPlot("sPlot_zjpsi_data","An SPlot", *zjpsi_data, sPlotmodel, RooArgList(*sig_prompt_yield, *sig_nonprompt_yield, *bkg_prompt_yield, *bkg_nonprompt_yield));
//  
//  RooDataSet * zjpsi_data_prompt_weighted     = new RooDataSet(zjpsi_data->GetName(), zjpsi_data->GetTitle(), zjpsi_data, *zjpsi_data->get(), 0, "Nzjpsi_m_sig_tau_sig_sw");
//  RooDataSet * zjpsi_data_nonprompt_weighted = new RooDataSet(zjpsi_data->GetName(), zjpsi_data->GetTitle(), zjpsi_data, *zjpsi_data->get(), 0, "Nzjpsi_m_sig_tau_bg_sw");
//  
//  TCut select_z_to_electrons  = "is_z_to_electrons==1";
//  TCut select_z_to_muons  = "is_z_to_muons==1";
//  
//  RooDataSet * zjpsi_data_prompt_Z_ee_weighted       = (RooDataSet*)zjpsi_data_prompt_weighted->reduce(select_z_to_electrons);
//  RooDataSet * zjpsi_data_nonprompt_Z_ee_weighted   = (RooDataSet*)zjpsi_data_nonprompt_weighted->reduce(select_z_to_electrons);
//  RooDataSet * zjpsi_data_prompt_Z_mumu_weighted     = (RooDataSet*)zjpsi_data_prompt_weighted->reduce(select_z_to_muons);
//  RooDataSet * zjpsi_data_nonprompt_Z_mumu_weighted = (RooDataSet*)zjpsi_data_nonprompt_weighted->reduce(select_z_to_muons);
//  
//  //int num_z_mass_bins = 10;
//  int num_z_mass_bins = 80;
//  TCanvas* czjpsi_data_Z = new TCanvas("czjpsi_data_Z","sPlot Z", 1200, 1200); 
//  czjpsi_data_Z->Divide(2,2);
//  RooPlot* frame_Z_ee_prompt = z_mass->frame(Bins(num_z_mass_bins)); 
//  zjpsi_data_prompt_Z_ee_weighted->plotOn(frame_Z_ee_prompt);
//
//  int z_mass_min_signal = 80;
//  int z_mass_max_signal = 100;
//  int z_mass_min_low = 40;
//  int z_mass_max_low_ee = 50;
//  int z_mass_max_low_mumu = 45;
//  int z_mass_min_high = 150;
//  int z_mass_max_high = 300;
//  TH1* h_z_ee_prompt_mass = zjpsi_data_prompt_Z_ee_weighted->createHistogram("h_z_ee_prompt", *z_mass);
//  TH1* h_z_ee_nonprompt_mass = zjpsi_data_nonprompt_Z_ee_weighted->createHistogram("h_z_ee_nonprompt", *z_mass);
//  TH1* h_z_mumu_prompt_mass = zjpsi_data_prompt_Z_mumu_weighted->createHistogram("h_z_mumu_prompt", *z_mass);
//  TH1* h_z_mumu_nonprompt_mass = zjpsi_data_nonprompt_Z_mumu_weighted->createHistogram("h_z_mumu_nonprompt", *z_mass);
//
//  TAxis *axis;
//  axis  = h_z_ee_prompt_mass->GetXaxis();
//  int bmin,bmax;
//
//  bmin = axis->FindBin(z_mass_min_signal);
//  bmax = axis->FindBin(z_mass_max_signal);
//  double integral_z_ee_prompt_signal = h_z_ee_prompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_ee prompt signal " << integral_z_ee_prompt_signal << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_signal);
//  bmax = axis->FindBin(z_mass_max_signal);
//  double integral_z_ee_nonprompt_signal = h_z_ee_nonprompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_ee nonprompt signal " << integral_z_ee_nonprompt_signal << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_signal);
//  bmax = axis->FindBin(z_mass_max_signal);
//  double integral_z_mumu_prompt_signal = h_z_mumu_prompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_mumu prompt signal " << integral_z_mumu_prompt_signal << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_signal);
//  bmax = axis->FindBin(z_mass_max_signal);
//  double integral_z_mumu_nonprompt_signal = h_z_mumu_nonprompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_mumu nonprompt signal " << integral_z_mumu_nonprompt_signal << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_low);
//  bmax = axis->FindBin(z_mass_max_low_ee);
//  double integral_z_ee_prompt_low = h_z_ee_prompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_ee prompt low " << integral_z_ee_prompt_low << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_high);
//  bmax = axis->FindBin(z_mass_max_high);
//  double integral_z_ee_prompt_high = h_z_ee_prompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_ee prompt high " << integral_z_ee_prompt_high << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_low);
//  bmax = axis->FindBin(z_mass_max_low_ee);
//  double integral_z_ee_nonprompt_low = h_z_ee_nonprompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_ee nonprompt low " << integral_z_ee_nonprompt_low << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_high);
//  bmax = axis->FindBin(z_mass_max_high);
//  double integral_z_ee_nonprompt_high = h_z_ee_nonprompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_ee nonprompt high " << integral_z_ee_nonprompt_high << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_low);
//  bmax = axis->FindBin(z_mass_max_low_mumu);
//  double integral_z_mumu_prompt_low = h_z_mumu_prompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_mumu prompt low " << integral_z_mumu_prompt_low << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_high);
//  bmax = axis->FindBin(z_mass_max_high);
//  double integral_z_mumu_prompt_high = h_z_mumu_prompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_mumu prompt high " << integral_z_mumu_prompt_high << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_low);
//  bmax = axis->FindBin(z_mass_max_low_mumu);
//  double integral_z_mumu_nonprompt_low = h_z_mumu_nonprompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_mumu nonprompt low " << integral_z_mumu_nonprompt_low << std::endl;
//
//  bmin = axis->FindBin(z_mass_min_high);
//  bmax = axis->FindBin(z_mass_max_high);
//  double integral_z_mumu_nonprompt_high = h_z_mumu_nonprompt_mass->Integral(bmin,bmax - 1);
//  //std::cout << "integral z_mumu nonprompt high " << integral_z_mumu_nonprompt_high << std::endl;
//
//  czjpsi_data_Z->cd(1); 
//  frame_Z_ee_prompt->SetTitle("Z->ee + Prompt J/#psi");
//  frame_Z_ee_prompt->Draw();
//  
//  RooPlot* frame_Z_ee_nonprompt = z_mass->frame(Bins(num_z_mass_bins));
//  zjpsi_data_nonprompt_Z_ee_weighted->plotOn(frame_Z_ee_nonprompt);
//  czjpsi_data_Z->cd(2); 
//  frame_Z_ee_nonprompt->SetTitle("Z->ee + Nonprompt J/#psi");
//  frame_Z_ee_nonprompt->Draw();
//  
//  RooPlot* frame_Z_mumu_prompt = z_mass->frame(Bins(num_z_mass_bins)); 
//  zjpsi_data_prompt_Z_mumu_weighted->plotOn(frame_Z_mumu_prompt);
//  czjpsi_data_Z->cd(3); 
//  frame_Z_mumu_prompt->SetTitle("Z->#mu#mu + Prompt J/#psi");
//  frame_Z_mumu_prompt->Draw();
//  
//  RooPlot* frame_Z_mumu_nonprompt = z_mass->frame(Bins(num_z_mass_bins)); 
//  zjpsi_data_nonprompt_Z_mumu_weighted->plotOn(frame_Z_mumu_nonprompt);
//  czjpsi_data_Z->cd(4); 
//  frame_Z_mumu_nonprompt->SetTitle("Z->#mu#mu + Nonprompt J/#psi");
//  frame_Z_mumu_nonprompt->Draw();
//
//  std::string czjpsi_data_Z_mass_image_name = image_dir;
//  czjpsi_data_Z_mass_image_name.append("z_mass");
//  czjpsi_data_Z_mass_image_name.append(".png");
//  //std::cout <<  czjpsi_data_Z_mass_image_name << std::endl;
//  czjpsi_data_Z->Print(czjpsi_data_Z_mass_image_name.c_str() , "png");
//  czjpsi_data_Z->Close();
//  
//  // ************************
//  // ************************
//  // Start the x-sec plotting
//  // ************************
//  // ************************
//  const RooArgSet* obs_zjpsi = zjpsi_data->get();
//  // inclusive
//  RooRealVar* prompt_weight = new RooRealVar("prompt_weight" ,"prompt_weight"  , -1000, 1000);
//  RooRealVar* nonprompt_weight = new RooRealVar("nonprompt_weight" ,"nonprompt_weight"  , -1000, 1000);
//  //unclear what zjpsi_data_prompt_tmp is ??
//  RooDataSet *zjpsi_data_prompt_tmp   = new RooDataSet("zjpsi_data_prompt_tmp", "zjpsi_data_prompt_tmp", RooArgSet(*prompt_weight));
//  RooDataSet *zjpsi_data_nonprompt_tmp   = new RooDataSet("zjpsi_data_nonprompt_tmp", "zjpsi_data_nonprompt_tmp", RooArgSet(*nonprompt_weight));
//
//  // fiducial
//  RooRealVar* prompt_weight_fid    = new RooRealVar("prompt_weight_fid" ,"prompt_weight_fid"  , -1000, 1000);
//  RooRealVar* nonprompt_weight_fid = new RooRealVar("nonprompt_weight_fid" ,"nonprompt_weight_fid"  , -1000, 1000);
//  //unclear what these datasets are, should have a better name
//  //TODO
//  RooDataSet *zjpsi_data_prompt_tmp_fid    = new RooDataSet("zjpsi_data_prompt_tmp_fid", "zjpsi_data_prompt_tmp_fid", RooArgSet(*prompt_weight_fid));
//  RooDataSet *zjpsi_data_nonprompt_tmp_fid = new RooDataSet("zjpsi_data_nonprompt_tmp_fid", "zjpsi_data_nonprompt_tmp_fid", RooArgSet(*nonprompt_weight_fid));
//
//  for(Int_t i=0; i < zjpsi_data->numEntries(); ++i) {
//    obs_zjpsi = zjpsi_data->get(i);
//    RooRealVar* effrow1 = new RooRealVar("effrow1", "effrow1", 1.);
//    RooRealVar* effrow2 = new RooRealVar("effrow2", "effrow2", 1.);
//    RooRealVar* accrow = new RooRealVar("accrow", "accrow", 1.);
//    if(use_polarisation_weights) {
//      //TODO verify that this is done correctly
//      //unclear what effrow1 and accrow are, or how they are calculated
//      //these varialbes should be renamed to more reflect what they are
//      //RooRealVar* effrow1 = (RooRealVar*) obs_zjpsi->find("reco_muon0_weight");
//      //RooRealVar* effrow2 = (RooRealVar*) obs_zjpsi->find("reco_muon1_weight");
//      //RooRealVar* accrow = (RooRealVar*) obs_zjpsi->find("unpolarised");
//      
//      
//      //effrow1 = (RooRealVar*) obs_zjpsi->find("reco_muon0_weight");
//      //effrow2 = (RooRealVar*) obs_zjpsi->find("reco_muon1_weight");
//      accrow = (RooRealVar*) obs_zjpsi->find("unpolarised");
//    }
//    else {
//      //RooRealVar* effrow1 = new RooRealVar("effrow1", "effrow1", 1.);
//      //RooRealVar* effrow2 = new RooRealVar("effrow2", "effrow2", 1.);
//      //RooRealVar* accrow = new RooRealVar("accrow", "accrow", 1.);
//    }
//    // fiducial calculations
//    double weighted_prompt_fid_value    = sPlot_zjpsi_data->GetSWeight(i,"Nzjpsi_m_sig_tau_sig_sw")/(effrow1->getVal()*effrow2->getVal());
//    double weighted_nonprompt_fid_value = sPlot_zjpsi_data->GetSWeight(i,"Nzjpsi_m_sig_tau_bg_sw")/(effrow1->getVal()*effrow2->getVal());
//    prompt_weight_fid->setVal(weighted_prompt_fid_value);
//    nonprompt_weight_fid->setVal(weighted_nonprompt_fid_value);
//    zjpsi_data_prompt_tmp_fid->add(RooArgList(*prompt_weight_fid));
//    zjpsi_data_nonprompt_tmp_fid->add(RooArgList(*nonprompt_weight_fid));
//
//    // inclusive calculations
//    double weighted_prompt_value    = sPlot_zjpsi_data->GetSWeight(i,"Nzjpsi_m_sig_tau_sig_sw")*accrow->getVal()/(effrow1->getVal()*effrow2->getVal());
//    double weighted_nonprompt_value = sPlot_zjpsi_data->GetSWeight(i,"Nzjpsi_m_sig_tau_bg_sw")*accrow->getVal()/(effrow1->getVal()*effrow2->getVal());
//    prompt_weight->setVal(weighted_prompt_value);
//    nonprompt_weight->setVal(weighted_nonprompt_value);
//    zjpsi_data_prompt_tmp->add(RooArgList(*prompt_weight));
//    zjpsi_data_nonprompt_tmp->add(RooArgList(*nonprompt_weight));
//  }
//  
//  zjpsi_data->merge(zjpsi_data_prompt_tmp); 
//  zjpsi_data->merge(zjpsi_data_nonprompt_tmp); 
//  RooDataSet * zjpsi_data_weighted_prompt    = new RooDataSet(zjpsi_data->GetName(), zjpsi_data->GetTitle(), zjpsi_data, *zjpsi_data->get(), 0, prompt_weight->GetName());  
//  RooDataSet * zjpsi_data_weighted_nonprompt = new RooDataSet(zjpsi_data->GetName(), zjpsi_data->GetTitle(), zjpsi_data, *zjpsi_data->get(), 0, nonprompt_weight->GetName());  
//
//  zjpsi_data_fid->merge(zjpsi_data_prompt_tmp_fid); 
//  zjpsi_data_fid->merge(zjpsi_data_nonprompt_tmp_fid); 
//  RooDataSet * zjpsi_data_weighted_prompt_fid    = new RooDataSet(zjpsi_data_fid->GetName(), zjpsi_data_fid->GetTitle(), zjpsi_data_fid, *zjpsi_data_fid->get(), 0, prompt_weight_fid->GetName());  
//  RooDataSet * zjpsi_data_weighted_nonprompt_fid = new RooDataSet(zjpsi_data_fid->GetName(), zjpsi_data_fid->GetTitle(), zjpsi_data_fid, *zjpsi_data_fid->get(), 0, nonprompt_weight_fid->GetName());  
//
//  RooBinning bins(6, 100) ;
//  bins.addUniform(1, 6, 8.5);
//  bins.addUniform(1, 8.5, 10);
//  bins.addUniform(2, 10, 18);
//  bins.addUniform(1, 18, 30);
//  bins.addUniform(1, 30, 100);
//
//  RooDataHist zjpsi_prompt_hist   ("zjpsi_prompt_hist"   , "zjpsi_prompt_hist", *onia_pt, *zjpsi_data_weighted_prompt);
//  RooDataHist zjpsi_nonprompt_hist("zjpsi_nonprompt_hist", "zjpsi_nonprompt_hist", *onia_pt, *zjpsi_data_weighted_nonprompt);
//
//  RooDataHist zjpsi_prompt_hist_fid   ("zjpsi_prompt_hist_fid"   , "zjpsi_prompt_hist_fid", *onia_pt, *zjpsi_data_weighted_prompt_fid);
//  RooDataHist zjpsi_nonprompt_hist_fid("zjpsi_nonprompt_hist_fid", "zjpsi_nonprompt_hist_fid", *onia_pt, *zjpsi_data_weighted_nonprompt_fid);
//
//  double inclusive_z_to_muons_events;
//  double inclusive_z_to_electrons_events;
//  double inclusive_z_to_leptons_events;
//  if (create_cs_plots) {
//    //inclusive_z_to_muons_events     = 7.773024E06;
//    //inclusive_z_to_electrons_events = 5.086549E06;
//    inclusive_z_to_muons_events     = 6.93504E06;
//    inclusive_z_to_electrons_events = 4.63707E06;
//  }
//  else {
//    inclusive_z_to_muons_events     = 0.5;
//    inclusive_z_to_electrons_events = 0.5;
//  }
//  inclusive_z_to_leptons_events = inclusive_z_to_muons_events + inclusive_z_to_electrons_events;
//
//  //Negative events for pileup make sense?
//
//  //TODO this should be automated if possible
//  //float pileup_inclusive_prompt[5] = {3.56813, 4.67389, 0.794457, 0.0, 0.0};
//  //float pileup_inclusive_nonprompt[5] = {2.59413,5.3454,0.30785,1.51939,0.871015};
//
//
//  //based on inverting at 3 cm for 0.5 cm I believe
//  //float pileup_inclusive_prompt[5] = {1.60329,1.67072,0.323216,0.12009,0};
//  //float pileup_inclusive_nonprompt[5] = {1.15364,2.51752,0.282309,0.189887,0};
//
//  //for pu0.5, inverted at 0.5 cm
//  float pileup_inclusive_prompt[5] = {2.16664,2.14942,0.331587,-0.0822937,-0.035038};
//  float pileup_inclusive_nonprompt[5] = {1.83111,2.52538,0.594624,1.29732,0.272261};
//
//  //for pu1.0, inverted at 1.0 cm
//  //float pileup_inclusive_prompt[5] = {4.30462,4.68395,0.805219,-0.0910758,-0.0902335};
//  //float pileup_inclusive_nonprompt[5] = {2.89074,4.7174,0.631182,1.95426,0.616307};
//
//
//  TCanvas* c_jpsi_diff_xsec = new TCanvas("c_jpsi_diff_xsec","J/#psi p_{T} sPlot weighted", 1400, 700); c_jpsi_diff_xsec->Divide(2,1);
//  // plot the prompt part
//  c_jpsi_diff_xsec->cd(1);
//  onia_pt->setBinning(bins);
//  TH1* h_zjpsi_atlas_prompt = zjpsi_prompt_hist.createHistogram("h_zjpsi_atlas_prompt", *onia_pt);
//  h_zjpsi_atlas_prompt->SetBinContent(2, 10.8E-7); h_zjpsi_atlas_prompt->SetBinError(2, TMath::Sqrt(5.6*5.6    +1.9*1.9    )*1E-7);
//  h_zjpsi_atlas_prompt->SetBinContent(3,  5.6E-7); h_zjpsi_atlas_prompt->SetBinError(3, TMath::Sqrt(1.9*1.9    +0.8*0.8    )*1E-7);
//  h_zjpsi_atlas_prompt->SetBinContent(4,  1.9E-7); h_zjpsi_atlas_prompt->SetBinError(4, TMath::Sqrt(1.1*1.1    +0.1*0.1    )*1E-7);
//  h_zjpsi_atlas_prompt->SetBinContent(5, 0.87E-7); h_zjpsi_atlas_prompt->SetBinError(5, TMath::Sqrt(0.37*0.37  +0.12*0.12  )*1E-7);
//  h_zjpsi_atlas_prompt->SetBinContent(6, 0.09E-7); h_zjpsi_atlas_prompt->SetBinError(6, TMath::Sqrt(0.037*0.037+0.012*0.012)*1E-7);
//  h_zjpsi_atlas_prompt->SetMarkerStyle(24);
//
//  TH1* h_zjpsi_prompt  = zjpsi_prompt_hist.createHistogram("h_zjpsi_prompt", *onia_pt);
//  TH1* h_zjpsi_prompt_fid  = zjpsi_prompt_hist_fid.createHistogram("h_zjpsi_prompt_fid", *onia_pt); 
//  //TODO fix this, make it clearer what it is
//  //h_zjpsi_prompt_cs = hh2PROMPT? no idea what this is or why it is here??
//  TH1* h_zjpsi_prompt_cs = zjpsi_prompt_hist.createHistogram("h_zjpsi_prompt_cs", *onia_pt);
//  h_zjpsi_prompt_cs->SetXTitle("#it{p}_{T}^{J/#psi} [GeV]");
//  h_zjpsi_prompt_cs->SetYTitle("#it{B}(J/#psi#rightarrow#mu#mu) #times #frac{1}{#sigma(Z)} #frac{d#sigma(Z+J/#psi)}{d#it{p}_{T}} [1/GeV]");
//  h_zjpsi_prompt_cs->SetLabelSize(0.03);
//  h_zjpsi_prompt_cs->SetTitleSize(0.04);
//  h_zjpsi_prompt->SetLineWidth(2);
//  h_zjpsi_prompt->SetLineColor(kBlue -2);
//  for (int i=2; i<7; i++) {
//    if (use_pileup_correction) {
//      h_zjpsi_prompt_cs->SetBinContent(i, (h_zjpsi_prompt->GetBinContent(i) - pileup_inclusive_prompt[i-2])/h_zjpsi_prompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//      h_zjpsi_prompt_cs->SetBinError(i, (h_zjpsi_prompt->GetBinError(i) )/h_zjpsi_prompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//    }
//    else {
//      h_zjpsi_prompt_cs->SetBinContent(i, h_zjpsi_prompt->GetBinContent(i)/h_zjpsi_prompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//      h_zjpsi_prompt_cs->SetBinError(i, h_zjpsi_prompt->GetBinError(i)/h_zjpsi_prompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//    }
//  }
//  h_zjpsi_prompt_cs->Draw("E1");
//
//  h_zjpsi_atlas_prompt->SetLineColor(kRed - 2);
//  h_zjpsi_atlas_prompt->Draw("E1same");
//
//  Double_t xl1_l5=.60, yl1_l5=0.70, xl2_l5=xl1_l5+.3, yl2_l5=yl1_l5+.2;
//  TLegend *leg5 = new TLegend(xl1_l5,yl1_l5,xl2_l5,yl2_l5);
//  leg5->SetFillColor(kWhite);
//  leg5->AddEntry(h_zjpsi_prompt_cs,"CMS prompt","l");
//  leg5->AddEntry(h_zjpsi_atlas_prompt,"ATLAS prompt","l");
//  leg5->SetShadowColor(0);
//  leg5->Draw();
//  if (create_cs_plots) {
//    h_zjpsi_prompt_cs->GetYaxis()->SetRangeUser(1E-11, 1E-5);
//  }
//  //TODO should it be 120?
//  h_zjpsi_prompt_cs->GetXaxis()->SetRangeUser(8.5, 100);
//  //splot_tex_pt.DrawLatex(15, 3E-6, "#bf{#it{CMS}} Preliminary, #sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
//  //splot_tex_pt.DrawLatex(15, 1E-6 , "#it{pp} #rightarrow prompt #it{J/#psi+Z} : #it{pp} #rightarrow #it{Z}");
//
//  c_jpsi_diff_xsec->cd(1)->SetLogy(1);
//  c_jpsi_diff_xsec->cd(1)->SetLogx(1);
//  // plot the non-prompt part
//  c_jpsi_diff_xsec->cd(2);
//  TH1* h_zjpsi_atlas_nonprompt = zjpsi_nonprompt_hist.createHistogram("h_zjpsi_atlas_nonprompt", *onia_pt);
//  h_zjpsi_atlas_nonprompt->SetBinContent(2,   5.1E-7); h_zjpsi_atlas_nonprompt->SetBinError(2, TMath::Sqrt(4.2*4.2    +0.9*0.9    )*1E-7);
//  h_zjpsi_atlas_nonprompt->SetBinContent(3,   9.2E-7); h_zjpsi_atlas_nonprompt->SetBinError(3, TMath::Sqrt(2.5*2.5    +1.2*1.2    )*1E-7);
//  h_zjpsi_atlas_nonprompt->SetBinContent(4,   3.3E-7); h_zjpsi_atlas_nonprompt->SetBinError(4, TMath::Sqrt(1.2*1.2    +0.4*0.4    )*1E-7);
//  h_zjpsi_atlas_nonprompt->SetBinContent(5,  3.04E-7); h_zjpsi_atlas_nonprompt->SetBinError(5, TMath::Sqrt(0.59*0.59  +0.04*0.04  )*1E-7);
//  h_zjpsi_atlas_nonprompt->SetBinContent(6, 0.115E-7); h_zjpsi_atlas_nonprompt->SetBinError(6, TMath::Sqrt(0.039*0.039+0.002*0.002)*1E-7);
//  h_zjpsi_atlas_nonprompt->SetMarkerStyle(24);
//
//  TH1* h_zjpsi_nonprompt  = zjpsi_nonprompt_hist.createHistogram("h_zjpsi_nonprompt", *onia_pt);
//  TH1* h_zjpsi_nonprompt_fid  = zjpsi_nonprompt_hist_fid.createHistogram("h_zjpsi_nonprompt_fid", *onia_pt);
//  TH1* h_zjpsi_nonprompt_cs = zjpsi_nonprompt_hist.createHistogram("h_zjpsi_nonprompt_cs", *onia_pt);
//  h_zjpsi_nonprompt_cs->SetXTitle("#it{p}_{T}^{J/#psi} [GeV]");
//  h_zjpsi_nonprompt_cs->SetYTitle("#it{B}(J/#psi#rightarrow#mu#mu) #times #frac{1}{#sigma(Z)} #frac{d#sigma(Z+J/#psi)}{d#it{p}_{T}} [1/GeV]");
//  h_zjpsi_nonprompt->SetLineWidth(2);
//  h_zjpsi_nonprompt->SetLineColor(kBlue -2);
//  for (int i=2; i<7; i++) {
//    if (use_pileup_correction) {
//      h_zjpsi_nonprompt_cs->SetBinContent(i, (h_zjpsi_nonprompt->GetBinContent(i) - pileup_inclusive_nonprompt[i-2])/h_zjpsi_nonprompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//      h_zjpsi_nonprompt_cs->SetBinError(i, (h_zjpsi_nonprompt->GetBinError(i))/h_zjpsi_nonprompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//    }
//    else {
//      h_zjpsi_nonprompt_cs->SetBinContent(i, h_zjpsi_nonprompt->GetBinContent(i)/h_zjpsi_nonprompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//      h_zjpsi_nonprompt_cs->SetBinError(i, h_zjpsi_nonprompt->GetBinError(i)/h_zjpsi_nonprompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//    }
//    //h_zjpsi_nonprompt_cs->SetBinContent(i, h_zjpsi_nonprompt->GetBinContent(i)/h_zjpsi_nonprompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//    //h_zjpsi_nonprompt_cs->SetBinError(i, h_zjpsi_nonprompt->GetBinError(i)/h_zjpsi_nonprompt->GetBinWidth(i)/(inclusive_z_to_leptons_events));
//  }
//  h_zjpsi_nonprompt_cs->Draw("E1");
//
//  h_zjpsi_atlas_nonprompt->SetLineColor(kRed - 2);
//  h_zjpsi_atlas_nonprompt->Draw("E1same");
//
//  Double_t xl1_l6=.60, yl1_l6=0.70, xl2_l6=xl1_l6+.3, yl2_l6=yl1_l6+.2;
//  TLegend *leg6 = new TLegend(xl1_l6,yl1_l6,xl2_l6,yl2_l6);
//  leg6->SetFillColor(kWhite);
//  leg6->AddEntry(h_zjpsi_nonprompt_cs,"CMS nonprompt","l");
//  leg6->AddEntry(h_zjpsi_atlas_nonprompt,"ATLAS nonprompt","l");
//  leg6->SetShadowColor(0);
//  leg6->Draw();
//
//  if (create_cs_plots) {
//    h_zjpsi_nonprompt_cs->GetYaxis()->SetRangeUser(1E-11, 1E-5);
//  }
//  //TODO should it be 120?
//  h_zjpsi_nonprompt_cs->GetXaxis()->SetRangeUser(8.5, 100);
//  //splot_tex_pt.DrawLatex(15, 3E-6, "#bf{#it{CMS}} Preliminary, #sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
//  //splot_tex_pt.DrawLatex(15, 1E-6 , "#it{pp} #rightarrow non-prompt #it{J/#psi+Z} : #it{pp} #rightarrow #it{Z}");
//  c_jpsi_diff_xsec->cd(2)->SetLogy(1);
//  c_jpsi_diff_xsec->cd(2)->SetLogx(1);
//  //c_jpsi_diff_xsec->SaveAs("c_jpsi_diff_xsec.C");
//  //c_jpsi_diff_xsec->Print();
//
//  std::string c_jpsi_diff_xsec_image_name = image_dir;
//  c_jpsi_diff_xsec_image_name.append("c_jpsi_diff_xsec");
//  c_jpsi_diff_xsec_image_name.append(".png");
//  //std::cout <<  c_jpsi_diff_xsec_image_name << std::endl;
//  c_jpsi_diff_xsec->Print(c_jpsi_diff_xsec_image_name.c_str() , "png");
//  c_jpsi_diff_xsec->Close();
//
//  // create total cross-section plots
//  // inclusive cross-section
//  double inclusive_cross_section_prompt=0., inclusive_cross_section_prompt_err=0.;
//  double inclusive_cross_section_nonprompt=0., inclusive_cross_section_nonprompt_err=0.;
//  double fiducial_cross_section_prompt=0., fiducial_cross_section_prompt_err=0.;
//  double fiducial_cross_section_nonprompt=0., fiducial_cross_section_nonprompt_err=0.;
//  
//  //float pileup_factor = 0.1161;
//  //float pileup_correction = pileup_factor / (1.0 - pileup_factor);
//  //for 0.5 pileup cut
//  float pileup_factor = 0.0580493;
//  //for 1 cm pileup
//  //float pileup_factor = 0.115794;
//
//  //for 0.5 pileup cut, but inverted at 3 cm
//  //float pileup_correction = pileup_factor / (0.662);
//  float pileup_correction = pileup_factor / (1 - pileup_factor);
//
//  float bin_size[5] = { 1.5, 4, 4, 12, 70 };
//  
//  float prompt_pileup_values[7];
//  float nonprompt_pileup_values[7];
//  for (int i=2; i<7; ++i) {
//    inclusive_cross_section_prompt         += h_zjpsi_prompt->GetBinContent(i);
//    inclusive_cross_section_prompt_err     += TMath::Power(h_zjpsi_prompt->GetBinError(i),2);
//    inclusive_cross_section_nonprompt     += h_zjpsi_nonprompt->GetBinContent(i);
//    inclusive_cross_section_nonprompt_err += TMath::Power(h_zjpsi_nonprompt->GetBinError(i),2);
//    fiducial_cross_section_prompt          += h_zjpsi_prompt_fid->GetBinContent(i);
//    fiducial_cross_section_prompt_err      += TMath::Power(h_zjpsi_prompt_fid->GetBinError(i),2);
//    fiducial_cross_section_nonprompt      += h_zjpsi_nonprompt_fid->GetBinContent(i);
//    fiducial_cross_section_nonprompt_err  += TMath::Power(h_zjpsi_nonprompt_fid->GetBinError(i),2);
//
//    //cout << "bin " << i << " prompt weighted events " << inclusive_cross_section_prompt / bin_size[i-2] << " +/- " << TMath::Sqrt(inclusive_cross_section_prompt_err) / bin_size[i-2] << endl;
//    //cout << "bin " << i << " nonprompt weighted events " << inclusive_cross_section_nonprompt / bin_size[i-2] << " +/- " << TMath::Sqrt(inclusive_cross_section_nonprompt_err) / bin_size[i-2] << endl;
//
//    //cout << "bin " << i << " prompt weighted " << h_zjpsi_prompt->GetBinContent(i)  << " +/- " << h_zjpsi_prompt->GetBinError(i) << endl;
//    //cout << "bin " << i << " nonprompt weighted " << h_zjpsi_nonprompt->GetBinContent(i) << " +/- " << h_zjpsi_nonprompt->GetBinError(i) << endl;
//
//    //cout << "bin " << i << " prompt weighted " << h_zjpsi_prompt->GetBinContent(i) * pileup_correction << " +/- " << h_zjpsi_prompt->GetBinError(i) * pileup_correction<< endl;
//    //cout << "bin " << i << " nonprompt weighted " << h_zjpsi_nonprompt->GetBinContent(i)* pileup_correction << " +/- " << h_zjpsi_nonprompt->GetBinError(i) * pileup_correction<< endl;
//    
//    prompt_pileup_values[i] = h_zjpsi_prompt->GetBinContent(i) * pileup_correction;
//    nonprompt_pileup_values[i] = h_zjpsi_nonprompt->GetBinContent(i) * pileup_correction;
//  }
//
//  //cout << "prompt corrected pu: " << "{" << prompt_pileup_values[2] << "," << prompt_pileup_values[3] << "," << prompt_pileup_values[4] << "," << prompt_pileup_values[5] << "," << prompt_pileup_values[6] << "}" << endl;
//  //cout << "nonprompt corrected pu: " << "{" << nonprompt_pileup_values[2] << "," << nonprompt_pileup_values[3] << "," << nonprompt_pileup_values[4] << "," << nonprompt_pileup_values[5] << "," << nonprompt_pileup_values[6] << "}" << endl;
//
//  for (int i=2; i<7; ++i) {
//    //cout << "bin " << i << " prompt weighted ratio " << h_zjpsi_prompt->GetBinContent(i) / inclusive_z_to_leptons_events << " +/- " << h_zjpsi_prompt->GetBinError(i) / inclusive_z_to_leptons_events << endl;
//    //cout << "bin " << i << " nonprompt weighted ratio " << h_zjpsi_nonprompt->GetBinContent(i) / inclusive_z_to_leptons_events << " +/- " << h_zjpsi_nonprompt->GetBinError(i) / inclusive_z_to_leptons_events << endl;
//  }
//  cout << "i = 2 => 8.5-10" << endl;
//  cout << "i = 3 => 10-14" << endl;
//  cout << "i = 4 => 14-18" << endl;
//  cout << "i = 5 => 18-30" << endl;
//  cout << "i = 6 => 30-100" << endl;
//  //cout << "total prompt inclusive: " << inclusive_cross_section_prompt << " +/- " << TMath::Sqrt(inclusive_cross_section_prompt_err) << " ratio " << inclusive_cross_section_prompt / inclusive_z_to_leptons_events
//  //  <<  " +/- " << TMath::Sqrt(inclusive_cross_section_prompt_err) / inclusive_z_to_leptons_events << endl;
//  //cout << "total nonprompt inclusive: " << inclusive_cross_section_nonprompt << " +/- " << TMath::Sqrt(inclusive_cross_section_nonprompt_err) << " ratio " << inclusive_cross_section_nonprompt / inclusive_z_to_leptons_events << " +/- " << TMath::Sqrt(inclusive_cross_section_nonprompt_err) / inclusive_z_to_leptons_events << endl;
//
//  inclusive_cross_section_prompt = inclusive_cross_section_prompt/(inclusive_z_to_leptons_events);
//  inclusive_cross_section_nonprompt = inclusive_cross_section_nonprompt/(inclusive_z_to_leptons_events);
//  //TODO why square then take the sqrt of the error here?? At the least varialbe name should be renamed
//  inclusive_cross_section_prompt_err = TMath::Sqrt(inclusive_cross_section_prompt_err)/(inclusive_z_to_leptons_events);
//  inclusive_cross_section_nonprompt_err = TMath::Sqrt(inclusive_cross_section_nonprompt_err)/(inclusive_z_to_leptons_events);
//
//  fiducial_cross_section_prompt         = fiducial_cross_section_prompt/(inclusive_z_to_leptons_events);
//  fiducial_cross_section_nonprompt     = fiducial_cross_section_nonprompt/(inclusive_z_to_leptons_events);
//  //TODO why square then take the sqrt of the error here?? At the least varialbe name should be renamed
//  fiducial_cross_section_prompt_err     = TMath::Sqrt(fiducial_cross_section_prompt_err)/(inclusive_z_to_leptons_events);
//  fiducial_cross_section_nonprompt_err = TMath::Sqrt(fiducial_cross_section_nonprompt_err)/(inclusive_z_to_leptons_events);
//
//  
//
//  TH1F *h_total_cross_section_prompt     = new TH1F("h_total_cross_section_prompt", "h_total_cross_section_prompt", 3, 0, 3.);
//  h_total_cross_section_prompt->GetYaxis()->SetRangeUser(0, 20E-6);
//  h_total_cross_section_prompt->SetYTitle("#it{B}(J/#psi#rightarrow#mu#mu) #times #frac{#sigma(Z+J/#psi)}{#sigma(Z)} [1/GeV]");
//  TH1F *h_total_cross_section_nonprompt = new TH1F("h_total_cross_section_nonprompt", "h_total_cross_section_nonprompt", 3, 0, 3.);
//  h_total_cross_section_nonprompt->GetYaxis()->SetRangeUser(0, 20E-6);
//  h_total_cross_section_nonprompt->SetYTitle("#it{B}(J/#psi#rightarrow#mu#mu) #times #frac{#sigma(Z+J/#psi)}{#sigma(Z)} [1/GeV]");
//
//  TH1F *h_total_cross_section_prompt_ATLAS     = new TH1F("h_total_cross_section_prompt_ATLAS", "h_total_cross_section_prompt_ATLAS", 3, 0, 3.);
//  h_total_cross_section_prompt_ATLAS->SetMarkerStyle(24);
//  h_total_cross_section_prompt_ATLAS->SetBinContent(1, 36.8E-7); 
//  h_total_cross_section_prompt_ATLAS->SetBinError(1, TMath::Sqrt(6.7*6.7+2.5*2.5)*1E-7);
//  h_total_cross_section_prompt_ATLAS->SetBinContent(2, 63E-7); 
//  h_total_cross_section_prompt_ATLAS->SetBinError(2, TMath::Sqrt(13*13+5*5)*1E-7);
//  TH1F *h_total_cross_section_nonprompt_ATLAS = new TH1F("h_total_cross_section_nonprompt_ATLAS", "h_total_cross_section_nonprompt_ATLAS", 3, 0, 3.);
//  h_total_cross_section_nonprompt_ATLAS->SetMarkerStyle(24);
//  h_total_cross_section_nonprompt_ATLAS->SetBinContent(1, 65.8E-7); h_total_cross_section_nonprompt_ATLAS->SetBinError(1, TMath::Sqrt(9.2*9.2+4.2*4.2)*1E-7);
//  h_total_cross_section_nonprompt_ATLAS->SetBinContent(2, 102E-7); h_total_cross_section_nonprompt_ATLAS->SetBinError(2, TMath::Sqrt(15*15+5*5)*1E-7);
//
//  h_total_cross_section_prompt->SetBinContent(1, fiducial_cross_section_prompt);
//  h_total_cross_section_prompt->SetBinError(1, fiducial_cross_section_prompt_err);
//  h_total_cross_section_prompt->SetBinContent(2, inclusive_cross_section_prompt);
//  h_total_cross_section_prompt->SetBinError(2, inclusive_cross_section_prompt_err);
//  h_total_cross_section_prompt->GetXaxis()->SetBinLabel(1, "Fiducial");
//  h_total_cross_section_prompt->GetXaxis()->SetBinLabel(2, "Inclusive");
//  h_total_cross_section_prompt->GetXaxis()->SetBinLabel(3, "DPS-subtracted");
//
//  h_total_cross_section_nonprompt->SetBinContent(1, fiducial_cross_section_nonprompt);
//  h_total_cross_section_nonprompt->SetBinError(1, fiducial_cross_section_nonprompt_err);
//  h_total_cross_section_nonprompt->SetBinContent(2, inclusive_cross_section_nonprompt);
//  h_total_cross_section_nonprompt->SetBinError(2, inclusive_cross_section_nonprompt_err);
//  h_total_cross_section_nonprompt->GetXaxis()->SetBinLabel(1, "Fiducial");
//  h_total_cross_section_nonprompt->GetXaxis()->SetBinLabel(2, "Inclusive");
//  h_total_cross_section_nonprompt->GetXaxis()->SetBinLabel(3, "DPS-subtracted");
//
//  // plotting cross-sections
//
//  TCanvas *c_total_cross_section = new TCanvas("c_total_cross_section", "c_total_cross_section", 1200, 600);
//  c_total_cross_section->Divide(2,1);
//  c_total_cross_section->cd(1);
//  h_total_cross_section_prompt->Draw("e1");
//  h_total_cross_section_prompt_ATLAS->Draw("e1same");
//
//  TLatex *   tex1 = new TLatex(.2,18E-6,"#bf{#it{CMS}}  Preliminary, #sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
//  tex1->SetTextFont(42);
//  tex1->SetTextSize(0.04);
//  tex1->SetLineWidth(2);
//  tex1->Draw();
//  TLatex *   tex2 = new TLatex(.2,16E-6,"#it{pp} #rightarrow prompt #it{J/#psi+Z} : #it{pp} #rightarrow #it{Z}");
//  tex2->SetTextFont(42);
//  tex2->SetTextSize(0.04);
//  tex2->SetLineWidth(2);
//  tex2->Draw();
//  TLatex *   tex3 = new TLatex(.2,14E-6,"|#it{y_{J/#psi}}| < 2.1, 8.5 < #it{p}_{T}^{#it{J/#psi}} < 100 GeV");
//  tex3->SetTextFont(42);
//  tex3->SetTextSize(0.04);
//  tex3->SetLineWidth(2);
//  tex3->Draw();
//
//  c_total_cross_section->cd(2); h_total_cross_section_nonprompt->Draw("e1"); h_total_cross_section_nonprompt_ATLAS->Draw("e1same");
//  TLatex *   tex4 = new TLatex(.2,18E-6,"#bf{#it{CMS}}  Preliminary, #sqrt{#it{s}}=8 TeV, 19.7 fb^{-1}");
//  tex4->SetTextFont(42);               
//  tex4->SetTextSize(0.04);             
//  tex4->SetLineWidth(2);               
//  tex4->Draw();                        
//  TLatex *   tex5 = new TLatex(.2,16E-6,"#it{pp} #rightarrow prompt #it{J/#psi+Z} : #it{pp} #rightarrow #it{Z}");
//  tex5->SetTextFont(42);               
//  tex5->SetTextSize(0.04);             
//  tex5->SetLineWidth(2);               
//  tex5->Draw();                        
//  TLatex *   tex6 = new TLatex(.2,14E-6,"|#it{y_{J/#psi}}| < 2.1, 8.5 < #it{p}_{T}^{#it{J/#psi}} < 100 GeV");
//  tex6->SetTextFont(42);
//  tex6->SetTextSize(0.04);
//  tex6->SetLineWidth(2);
//  tex6->Draw();
//
//  std::string c_total_cross_section_image_name = image_dir;
//  c_total_cross_section_image_name.append("c_total_cross_section");
//  c_total_cross_section_image_name.append(".png");
//  //std::cout <<  c_total_cross_section_image_name << std::endl;
//  c_total_cross_section->Print(c_total_cross_section_image_name.c_str() , "png");
//  c_total_cross_section->Close();
}

