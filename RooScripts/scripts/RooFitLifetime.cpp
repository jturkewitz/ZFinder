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
#include "TGraphErrors.h"

// RooFit
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitLifetime.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooKeysPdf.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"


using namespace RooFit;

TCanvas* get_tcanvas(const int X_DIM, const int Y_DIM) {
  TCanvas* tcan = new TCanvas("canvas", "canvas", X_DIM, Y_DIM);
  tcan->Divide(1); // Split in two side-by-side areas
  return tcan;
}

std::vector<double> RooFitLifetime(
    const std::string& DATA_FILE_1,
    const std::string& DATA_FILE_2,
    const std::string& OUT_DIR,
    const int PT_SLICE,
    const bool USE_PT_SLICES,
    const int RAP_SLICE,
    const bool USE_RAP_SLICES
    ) {
  // Open the data file
  TFile* f_data_1 = new TFile(DATA_FILE_1.c_str(), "READ");
  std::vector<double> zjpsi_info ;
  if (f_data_1 == NULL) {
    std::cout << "Data file is invalid" << std::endl;
    return zjpsi_info;
  }
  TFile* f_data_2 = new TFile(DATA_FILE_2.c_str(), "READ");
  if (f_data_2 == NULL) {
    std::cout << "Data file_2 is invalid" << std::endl;
    return zjpsi_info;
  }
  // Pass the open files to the main RooFitter
  const std::vector<double> RET_CODE = RooFitLifetime(f_data_1, f_data_2, OUT_DIR, PT_SLICE, USE_PT_SLICES, RAP_SLICE, USE_RAP_SLICES);

  // Clean up and return the exit code
  delete f_data_1;
  delete f_data_2;

  return RET_CODE;
}

std::vector<double> RooFitLifetime(
    TFile* const DATA_FILE_1,
    TFile* const DATA_FILE_2,
    const std::string& OUT_DIR,
    const int PT_SLICE,
    const bool USE_PT_SLICES,
    const int RAP_SLICE,
    const bool USE_RAP_SLICES
    ) {
  // Constants
  //const int N_CPU = 8;
  const int N_CPU = 1;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  gErrorIgnoreLevel = kWarning;

  double tau_xy_min = -0.3;
  double tau_xy_max = 5.0;
  // Set up the variables we're going to read in from the files
  RooRealVar tau_xy("tau_xy", "tau_xy" , tau_xy_min, tau_xy_max, "ps");
  //tau_xy.setRange("range_tau_xy", -0.3, 5.0);
  std::vector<double> zjpsi_info ;

  RooRealVar zjpsi_tau_xy("zjpsi_tau_xy", "zjpsi_tau_xy" , tau_xy_min, tau_xy_max, "ps");
  //zjpsi_tau_xy.setRange("zjpsi_range_tau_xy", -0.3, 5.0);

  //TODO name this more script readable - maybe change ZFinder names
  const std::string pt_slice_list[8] = {
    "jpsi_tau_xy_very_fine_all",
    "jpsi_tau_xy_very_fine_ptUnder10",
    "jpsi_tau_xy_very_fine_pt_10_to_15",
    "jpsi_tau_xy_very_fine_pt_15_to_20",
    "jpsi_tau_xy_very_fine_pt_20_to_25",
    "jpsi_tau_xy_very_fine_pt_25_to_30",
    "jpsi_tau_xy_very_fine_ptAbove20",
    "jpsi_tau_xy_very_fine_ptAbove30"};

  const std::string rap_slice_list[8] = {
    "jpsi_tau_xy_very_fine_all",
    "jpsi_tau_xy_very_fine_rap_0.0_to_0.3",
    "jpsi_tau_xy_very_fine_rap_0.3_to_0.6",
    "jpsi_tau_xy_very_fine_rap_0.6_to_0.9",
    "jpsi_tau_xy_very_fine_rap_0.9_to_1.2",
    "jpsi_tau_xy_very_fine_rap_1.2_to_1.5",
    "jpsi_tau_xy_very_fine_rap_1.5_to_1.8",
    "jpsi_tau_xy_very_fine_rap_1.8_to_2.1",
    };

  if(USE_PT_SLICES && ( PT_SLICE >= 8 || PT_SLICE < 0) ) {
    std::cout << "PT_SLICE >=8 exiting" << std::endl;
    return zjpsi_info;
  }
  if(USE_RAP_SLICES && ( RAP_SLICE >= 8 || RAP_SLICE < 0) ) {
    std::cout << "RAP_SLICE >=8 exiting" << std::endl;
    return zjpsi_info;
  }

//  std::string inclusive_jpsi_hist = "ZFinder/Jpsi/";
//  std::string zjpsi_hist = "ZFinder/Jpsi_And_Z/";
  std::string inclusive_jpsi_hist = "ZFinder/Jpsi_Primary_Vertex/";
  std::string zjpsi_hist = "ZFinder/Jpsi_And_Z_Same_Vertex/";
  std::string jpsi_hist_name = "";
  std::string zjpsi_hist_name = "z";

  if (USE_PT_SLICES) {
    inclusive_jpsi_hist.append( pt_slice_list[PT_SLICE] );
    zjpsi_hist.append( pt_slice_list[PT_SLICE] );
    jpsi_hist_name.append(pt_slice_list[PT_SLICE]);
    zjpsi_hist_name.append(pt_slice_list[PT_SLICE]);
  }
  else if (USE_RAP_SLICES) {
    inclusive_jpsi_hist.append( rap_slice_list[RAP_SLICE] );
    zjpsi_hist.append( rap_slice_list[RAP_SLICE] );
    jpsi_hist_name.append(rap_slice_list[RAP_SLICE]);
    zjpsi_hist_name.append(rap_slice_list[RAP_SLICE]);
  }
  else {
    std::cout << "MUST USE RAP OR PT FOR NOW: " << std::endl;
    return zjpsi_info;
  }

  TH1D *h_jpsi_tau_xy_very_fine = (TH1D*) DATA_FILE_1->Get( inclusive_jpsi_hist.c_str() );
  TH1D *h_zjpsi_tau_xy_very_fine = (TH1D*) DATA_FILE_2->Get( zjpsi_hist.c_str() );

  h_jpsi_tau_xy_very_fine->Sumw2();
  h_jpsi_tau_xy_very_fine->Rebin(2);

  h_zjpsi_tau_xy_very_fine->Sumw2();
  h_zjpsi_tau_xy_very_fine->Rebin(10);

  RooDataHist tau_xy_data_hist("tau_xy_data_hist", jpsi_hist_name.c_str(), tau_xy, h_jpsi_tau_xy_very_fine);
  RooDataHist zjpsi_tau_xy_data_hist("zjpsi_tau_xy_data_hist", zjpsi_hist_name.c_str(), zjpsi_tau_xy, h_zjpsi_tau_xy_very_fine);

  RooRealVar gaussmean("gaussmean", "Mean of the smearing Gaussian", 0.);
  RooRealVar gausssigma("gausssigma", "Width of the smearing Gaussian", 0.01, 0.005, 0.5);
  //RooRealVar gausssigma("gausssigma", "Width of the smearing Gaussian", 0.093);
  RooGaussModel smear_gauss_model("smear_gauss", "Gaussian used to smear the Exponential Decay", tau_xy, gaussmean, gausssigma);
  RooRealVar decay_lifetime("decay_lifetime", "Tau", 1.3, 0.05, 2.);
  RooDecay decay_exp("decay_exp", "Exponential Decay", tau_xy, decay_lifetime, smear_gauss_model, RooDecay::SingleSided);

  RooRealVar gauss_prompt_mean("gauss_prompt_mean", "Mean of the Prompt Gaussian", 0., -0.2, 0.2);
  RooRealVar gauss_prompt_sigma("gauss_prompt_sigma", "Width of the Prompt Gaussian", 0.007, 0.005, 0.25);
  RooGaussian prompt_gauss("prompt_gauss", "Gaussian of the Prompt Peak", tau_xy, gauss_prompt_mean, gauss_prompt_sigma);

  RooRealVar gauss_prompt_sigma_2("gauss_prompt_sigma_2", "Width of the Prompt Gaussian_2", 0.01, 0.008, 0.4);
  RooGaussian prompt_gauss_2("prompt_gauss_2", "Gaussian_2 of the Prompt Peak", tau_xy, gauss_prompt_mean, gauss_prompt_sigma_2);

  RooRealVar prompt_sharp_fraction("prompt_sharp_fraction", "prompt_sharp_fraction" , 0.5 , 0.0, 1.);
//  TODO decide whether or not to use double gaussian fit, if so also figure out how to fix which one is prompt or not
//  RooRealVar prompt_sharp_fraction("prompt_sharp_fraction", "prompt_sharp_fraction" , 1.0);
  RooAddPdf tau_xy_gauss_sum_fitpdf("tau_xy_gauss_sum_fitpdf", "tau_xy_gauss_sum_fitpdf", RooArgList(prompt_gauss, prompt_gauss_2), RooArgList(prompt_sharp_fraction));

  RooRealVar prompt_fraction("prompt_fraction", "prompt_fraction" , 0.3 , 0.0, 1.);
  //RooAddPdf tau_xy_fitpdf("tau_xy_fitpdf", "tau_xy_fitpdf", RooArgList(tau_xy_gauss_sum_fitpdf, decay_exp), RooArgList(prompt_fraction));
  //RooAddPdf tau_xy_fitpdf("tau_xy_fitpdf", "tau_xy_fitpdf", RooArgSet(tau_xy_gauss_sum_fitpdf, decay_exp), RooArgList(prompt_fraction));
  RooAddPdf tau_xy_fitpdf("tau_xy_fitpdf", "tau_xy_fitpdf", RooArgList(tau_xy_gauss_sum_fitpdf, decay_exp), RooArgList(prompt_fraction));

  RooFitResult *jpsi_fitres = tau_xy_fitpdf.fitTo(tau_xy_data_hist, Range(tau_xy_min, tau_xy_max), NumCPU(N_CPU), Verbose(false), PrintLevel(-1), Save());
  std::cout << jpsi_hist_name << std::endl;
  //jpsi_fitres->Print("v");
  jpsi_fitres->Print();

  double zjpsi_gaussmean_value = gaussmean.getVal();
  double zjpsi_gausssigma_value = gausssigma.getVal();
  double zjpsi_decay_lifetime_value = decay_lifetime.getVal();
  double zjpsi_gauss_prompt_mean_value = gauss_prompt_mean.getVal();
  double zjpsi_gauss_prompt_sigma_value = gauss_prompt_sigma.getVal();
  double zjpsi_gauss_prompt_sigma_2_value = gauss_prompt_sigma_2.getVal();
  double zjpsi_prompt_sharp_fraction_value = prompt_sharp_fraction.getVal();


  double zjpsi_gauss_prompt_sigma_error = gauss_prompt_sigma.getError();

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

  RooRealVar zjpsi_prompt_sharp_fraction("zjpsi_prompt_sharp_fraction", "zjpsi_prompt_sharp_fraction" , zjpsi_prompt_sharp_fraction_value);
  RooAddPdf zjpsi_tau_xy_gauss_sum_fitpdf("zjpsi_tau_xy_gauss_sum_fitpdf", "zjpsi_tau_xy_gauss_sum_fitpdf", RooArgList(zjpsi_prompt_gauss,zjpsi_prompt_gauss_2), RooArgList(zjpsi_prompt_sharp_fraction));

  RooRealVar zjpsi_prompt_fraction("zjpsi_prompt_fraction", "zjpsi_prompt_fraction" , 0.01 , 0.0, 1);
  RooAddPdf zjpsi_tau_xy_fitpdf("zjpsi_tau_xy_fitpdf", "zjpsi_tau_xy_fitpdf", RooArgList(zjpsi_tau_xy_gauss_sum_fitpdf, zjpsi_decay_exp), RooArgList(zjpsi_prompt_fraction));

  //TODO TESTING allowing zjpsi parameters to float

  //    RooRealVar zjpsi_gaussmean("zjpsi_gaussmean", "Mean of the smearing zjpsi Gaussian", zjpsi_gaussmean_value);
  //    RooRealVar zjpsi_gausssigma("zjpsi_gausssigma", "Width of the smearing zjpsi Gaussian", zjpsi_gausssigma_value, 0, zjpsi_gausssigma_value*3);
  //    RooGaussModel zjpsi_smear_gauss_model("zjpsi_smear_gauss", "Gaussian used to smear the zjpsi Exponential Decay", zjpsi_tau_xy, zjpsi_gaussmean, zjpsi_gausssigma);
  //    RooRealVar zjpsi_decay_lifetime("zjpsi_decay_lifetime", "zjpsi_Tau", zjpsi_decay_lifetime_value, zjpsi_decay_lifetime_value*0.5, zjpsi_decay_lifetime_value*2.0 );
  //    RooDecay zjpsi_decay_exp("zjpsi_decay_exp", "zjpsi_Exponential Decay", zjpsi_tau_xy, zjpsi_decay_lifetime, zjpsi_smear_gauss_model, RooDecay::SingleSided);
  //
  //    RooRealVar zjpsi_gauss_prompt_mean("zjpsi_gauss_prompt_mean", "Mean of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_mean_value, -0.4, 0.4);
  //    RooRealVar zjpsi_gauss_prompt_sigma("zjpsi_gauss_prompt_sigma", "Width of the Prompt zjpsi_Gaussian", zjpsi_gauss_prompt_sigma_value, zjpsi_gauss_prompt_sigma_value*0.5,
  //       zjpsi_gauss_prompt_sigma_value*2.0 );
  //    RooGaussian zjpsi_prompt_gauss("zjpsi_prompt_gauss", "Gaussian of the Prompt zjpsi_Peak", zjpsi_tau_xy, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma);
  //
  //    RooRealVar zjpsi_gauss_prompt_sigma_2("zjpsi_gauss_prompt_sigma_2", "Width of the Prompt zjpsi_Gaussian_2", zjpsi_gauss_prompt_sigma_2_value,zjpsi_gauss_prompt_sigma_2_value*0.5,zjpsi_gauss_prompt_sigma_2_value*2.0 );
  //    RooGaussian zjpsi_prompt_gauss_2("zjpsi_prompt_gauss_2", "Gaussian_2 of the Prompt zjpsi_Peak", zjpsi_tau_xy, zjpsi_gauss_prompt_mean, zjpsi_gauss_prompt_sigma_2);
  //
  //    RooRealVar zjpsi_prompt_sharp_fraction("zjpsi_prompt_sharp_fraction", "zjpsi_prompt_sharp_fraction" , zjpsi_prompt_sharp_fraction_value, 0.0, zjpsi_prompt_sharp_fraction_value*2.0);
  //    RooAddPdf zjpsi_tau_xy_gauss_sum_fitpdf("zjpsi_tau_xy_gauss_sum_fitpdf", "zjpsi_tau_xy_gauss_sum_fitpdf", RooArgList(zjpsi_prompt_gauss,zjpsi_prompt_gauss_2), RooArgList(zjpsi_prompt_sharp_fraction));
  //
  //    RooRealVar zjpsi_prompt_fraction("zjpsi_prompt_fraction", "zjpsi_prompt_fraction" , 0.01 , 0.0, 1);
  //    RooAddPdf zjpsi_tau_xy_fitpdf("zjpsi_tau_xy_fitpdf", "zjpsi_tau_xy_fitpdf", RooArgList(zjpsi_tau_xy_gauss_sum_fitpdf, zjpsi_decay_exp), RooArgList(zjpsi_prompt_fraction));

  //TCanvas* const canvas = get_tcanvas(1500, 750);
  TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 750);

  // Plot the left side
  canvas->cd(1);
  gPad->SetLogy();
  //RooPlot* tau_xy_fitframe = tau_xy.frame(-0.3,5.0);
  RooPlot* tau_xy_fitframe = tau_xy.frame( Title(jpsi_hist_name.c_str()) , Range(tau_xy_min, tau_xy_max ));
  //tau_xy_fitframe->SetName(0); // Unset title
  tau_xy_data_hist.plotOn(tau_xy_fitframe);
  tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(decay_exp), LineColor(kGreen-2));
  tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(prompt_gauss), LineColor(kBlue-2));
  tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(prompt_gauss_2), LineColor(kOrange-2));
  tau_xy_fitpdf.plotOn(tau_xy_fitframe, LineColor(kRed-2));
  //tau_xy_fitpdf.plotOn(tau_xy_fitframe, LineColor(kGreen-2), NumCPU(N_CPU));

  tau_xy_fitframe->Draw();

  std::string jpsi_image_name = OUT_DIR;
  jpsi_image_name.append(jpsi_hist_name);
  jpsi_image_name.append(".png");
  canvas->Print(jpsi_image_name.c_str() , "png");
  canvas->Close();

  RooFitResult * zjpsi_fitres = zjpsi_tau_xy_fitpdf.fitTo(zjpsi_tau_xy_data_hist, Range(tau_xy_min, tau_xy_max), NumCPU(N_CPU), Verbose(false), PrintLevel(-1) , Save());

  std::cout << zjpsi_hist_name << std::endl;
  double zjpsi_prompt_events = 0.;
  double zjpsi_prompt_events_error = 0.;
  TAxis *axis = h_zjpsi_tau_xy_very_fine->GetXaxis();
  int bmin = axis->FindBin(tau_xy_min);
  int bmax = axis->FindBin(tau_xy_max);
  double integral = h_zjpsi_tau_xy_very_fine->Integral(bmin,bmax);
  integral -= h_zjpsi_tau_xy_very_fine->GetBinContent(bmin)*(tau_xy_min - axis->GetBinLowEdge(bmin)) / axis->GetBinWidth(bmin);
  integral -= h_zjpsi_tau_xy_very_fine->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax) - tau_xy_max)  / axis->GetBinWidth(bmax); 

  zjpsi_prompt_events = integral * zjpsi_prompt_fraction.getVal();
  zjpsi_prompt_events_error = integral * zjpsi_prompt_fraction.getError();
  std::cout << "zjpsi_prompt_events: " << zjpsi_prompt_events << std::endl;
  std::cout << "zjpsi_prompt_events_error: " << zjpsi_prompt_events_error << std::endl;
  //zjpsi_fitres->Print("v");
  zjpsi_fitres->Print();

  zjpsi_info.push_back(zjpsi_prompt_events);
  zjpsi_info.push_back(zjpsi_prompt_events_error);
  zjpsi_info.push_back(zjpsi_gauss_prompt_sigma_value);
  zjpsi_info.push_back(zjpsi_gauss_prompt_sigma_error);

  TCanvas *zjpsi_canvas = new TCanvas("zjpsi_canvas", "zjpsi_canvas", 2000, 750);
  zjpsi_canvas->cd(1);
  gPad->SetLogy();
  //RooPlot* tau_xy_fitframe = tau_xy.frame(-0.3,5.0);
  RooPlot* zjpsi_tau_xy_fitframe = zjpsi_tau_xy.frame( Title(zjpsi_hist_name.c_str()) );
  //zjpsi_tau_xy_fitframe->SetName(0); // Unset title
  zjpsi_tau_xy_data_hist.plotOn(zjpsi_tau_xy_fitframe);
  zjpsi_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_decay_exp), LineColor(kGreen-2) );
  zjpsi_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_prompt_gauss), LineColor(kBlue-2) );
  zjpsi_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_prompt_gauss_2), LineColor(kOrange-2));
  zjpsi_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, LineColor(kRed-2));

  zjpsi_tau_xy_fitframe->Draw();

  std::string zjpsi_image_name = OUT_DIR;
  zjpsi_image_name.append(zjpsi_hist_name);
  zjpsi_image_name.append(".png");
  zjpsi_canvas->Print(zjpsi_image_name.c_str(), "png");
  zjpsi_canvas->Close();
  
  //return zjpsi_prompt_events;
  return zjpsi_info;
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
    const std::string DATA_FILE_2(argv[2]);
    const std::string OUT_DIR(argv[3]);

    double zjpsi_prompt_events_pt_all = 0.0;
    double zjpsi_prompt_events_pt_all_error = 0.0;
    double zjpsi_prompt_events_pt_summed = 0.0;
    double zjpsi_prompt_events_pt_summed_sq_error = 0.0;
    double zjpsi_prompt_events_rap_all = 0.0;
    double zjpsi_prompt_events_rap_all_error = 0.0;
    double zjpsi_prompt_events_rap_summed = 0.0;
    double zjpsi_prompt_events_rap_summed_sq_error = 0.0;

    std::vector <double> jpsi_prompt_sigma_pt;
    std::vector <double> jpsi_prompt_sigma_pt_error;
    std::vector <double> jpsi_prompt_sigma_rap;
    std::vector <double> jpsi_prompt_sigma_rap_error;

    for (int i=0 ; i < 8 ; ++i) {
      //for now, if i==0, jpsi_all
      std::vector <double> zjpsi_info = RooFitLifetime( DATA_FILE_1, DATA_FILE_2, OUT_DIR, i, true, 0, false );
      if (zjpsi_info.size() == 4)
      {
        jpsi_prompt_sigma_pt.push_back(zjpsi_info[2]);
        jpsi_prompt_sigma_pt_error.push_back(zjpsi_info[3]);
        if (i == 0 ) {
          zjpsi_prompt_events_pt_all = zjpsi_info[0] ;
          zjpsi_prompt_events_pt_all_error = zjpsi_info[1];
        }
        //above 20 pt skip as use above 30 pt bin instead
        else if (i == 6 ) {
          //do nothing ignore above 20 pt for now
        }
        else {
          zjpsi_prompt_events_pt_summed += zjpsi_info[0] ;
          zjpsi_prompt_events_pt_summed_sq_error += pow(zjpsi_info[1],2) ;
        }
      }
    }
    for (int i=0 ; i < 8 ; ++i) {
      //for now, if i==0, jpsi_all
      std::vector <double> zjpsi_info = RooFitLifetime( DATA_FILE_1, DATA_FILE_2, OUT_DIR, 0, false, i, true );
      if (zjpsi_info.size() == 4)
      {
        jpsi_prompt_sigma_rap.push_back(zjpsi_info[2]);
        jpsi_prompt_sigma_rap_error.push_back(zjpsi_info[3]);
        if (i == 0 ) {
          zjpsi_prompt_events_rap_all = zjpsi_info[0] ;
          zjpsi_prompt_events_rap_all_error = zjpsi_info[1];
        }
        else {
          zjpsi_prompt_events_rap_summed += zjpsi_info[0] ;
          zjpsi_prompt_events_rap_summed_sq_error += pow(zjpsi_info[1],2) ;
        }
      }
    }
    std::cout << "PT bin: 0 --- all " << std::endl;
    std::cout << "1 --- < 10" << std::endl;
    std::cout << "2 --- 10-15" << std::endl;
    std::cout << "3 --- 15-20" << std::endl;
    std::cout << "4 --- 20-25" << std::endl;
    std::cout << "5 --- 25-30" << std::endl;
    std::cout << "6 --- > 20" << std::endl;
    std::cout << "7 --- > 30" << std::endl;
    double pt_bin[6];
    double pt_bin_error[6];
    double pt_prompt_sigma[6];
    double pt_prompt_sigma_error[6];
    for (unsigned int i=0 ; i < jpsi_prompt_sigma_pt.size() ; ++i) {
      std::cout << "pt bin: " << i << " prompt sigma " << jpsi_prompt_sigma_pt.at(i) << std::endl;
      if (i == 0 || i==6) {
        continue;
        //first pt bin is all, 6th bin is above 20 GeV
      }
      if ( i == 1 ) {
        pt_bin[i-1] = 5 ; //first pt pin is under 10
        pt_bin_error[i-1] = 5;
        pt_prompt_sigma[i-1] = jpsi_prompt_sigma_pt.at(i);
        pt_prompt_sigma_error[i-1] = jpsi_prompt_sigma_pt_error.at(i);
      }
      else if (i == 7 ) {
        pt_bin[5] = 45.;
        pt_bin_error[5] = 15.;
        pt_prompt_sigma[5] = jpsi_prompt_sigma_pt.at(i);
        pt_prompt_sigma_error[5] = jpsi_prompt_sigma_pt_error.at(i);
      }

      else {
        pt_bin[i-1] = (i-1.0) * 5 + 7.5;
        pt_bin_error[i-1] = 2.5;
        pt_prompt_sigma[i-1] = jpsi_prompt_sigma_pt.at(i);
        pt_prompt_sigma_error[i-1] = jpsi_prompt_sigma_pt_error.at(i);
      }
    }

    TGraphErrors *pt_graph = new TGraphErrors(6, pt_bin, pt_prompt_sigma, pt_bin_error, pt_prompt_sigma_error);
    TCanvas *pt_graph_canvas = new TCanvas("pt_graph_canvas", "pt_graph_canvas", 2000, 750);
    pt_graph_canvas->cd(1);
    pt_graph->Draw("AP");
    pt_graph->GetXaxis()->SetTitle("Jpsi Pt");
    pt_graph->GetYaxis()->SetTitle("Prompt Sigma [ps]");
    pt_graph->SetTitle("Prompt Gaussian Resolution");
    std::string pt_graph_image_name = OUT_DIR;
    pt_graph_image_name.append("pt_prompt_sigma.png");
    pt_graph_canvas->Print(pt_graph_image_name.c_str(), "png");
    pt_graph_canvas->Close();

    std::cout << "RAP bin: 0 --- all " << std::endl;
    std::cout << "1 --- 0.0-0.3" << std::endl;
    std::cout << "2 --- 0.3-0.6" << std::endl;
    std::cout << "3 --- 0.6-0.9" << std::endl;
    std::cout << "4 --- 0.9-1.2" << std::endl;
    std::cout << "5 --- 1.2-1.5" << std::endl;
    std::cout << "6 --- > 1.5-1.8" << std::endl;
    std::cout << "7 --- > 1.8-2.1" << std::endl;
    double rap_bin[7];
    double rap_bin_error[7];
    double rap_prompt_sigma[7];
    double rap_prompt_sigma_error[7];
    for (unsigned int i=0 ; i < jpsi_prompt_sigma_rap.size() ; ++i) {
      std::cout << "rap bin: " << i << " prompt sigma " << jpsi_prompt_sigma_rap.at(i) << std::endl;
      if (i == 0 ) {
        continue;
        //first rap bin is all
      }
      rap_bin[i-1] = (i - 1.0) * 0.3 + 0.15 ;
      rap_bin_error[i-1] = 0.15;
      rap_prompt_sigma[i-1] = jpsi_prompt_sigma_rap.at(i);
      rap_prompt_sigma_error[i-1] = jpsi_prompt_sigma_rap_error.at(i);
    }
    TGraphErrors *rap_graph = new TGraphErrors(7, rap_bin, rap_prompt_sigma, rap_bin_error, rap_prompt_sigma_error);
    TCanvas *rap_graph_canvas = new TCanvas("rap_graph_canvas", "rap_graph_canvas", 2000, 750);
    rap_graph_canvas->cd(1);
    rap_graph->Draw("AP");
    rap_graph->GetXaxis()->SetTitle("Jpsi Rapidity");
    rap_graph->GetYaxis()->SetTitle("Prompt Sigma [ps]");
    rap_graph->SetTitle("Prompt Gaussian Resolution");
    std::string rap_graph_image_name = OUT_DIR;
    rap_graph_image_name.append("rap_prompt_sigma.png");
    rap_graph_canvas->Print(rap_graph_image_name.c_str(), "png");
    rap_graph_canvas->Close();

    std::cout << "zjpsi_prompt_events_pt_all: " << zjpsi_prompt_events_pt_all << std::endl;
    std::cout << "zjpsi_prompt_events_pt_all_error: " << zjpsi_prompt_events_pt_all_error << std::endl;
    std::cout << "zjpsi_prompt_events_pt_summed: " << zjpsi_prompt_events_pt_summed << std::endl;
    std::cout << "zjpsi_prompt_events_pt_summed_error: " << pow(zjpsi_prompt_events_pt_summed_sq_error , 0.5) << std::endl;
    std::cout << "zjpsi_prompt_events_rap_all: " << zjpsi_prompt_events_rap_all << std::endl;
    std::cout << "zjpsi_prompt_events_rap_all_error: " << zjpsi_prompt_events_rap_all_error << std::endl;
    std::cout << "zjpsi_prompt_events_rap_summed: " << zjpsi_prompt_events_rap_summed << std::endl;
    std::cout << "zjpsi_prompt_events_rap_summed_error: " << pow(zjpsi_prompt_events_rap_summed_sq_error , 0.5) << std::endl;
    return 0;
  }
}
