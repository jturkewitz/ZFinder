// Standard Library
#include <fstream>
#include <iostream>
#include <sstream>

// ROOT
#include <TCanvas.h>
#include "TH1.h"

// RooFit
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFFTConvPdf.h"
#include "RooFitLifetime.h"
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

int RooFitLifetime(
        const std::string& DATA_FILE_1,
        const std::string& DATA_FILE_2,
        const std::string& OUT_DIR
        ) {
    // Open the data file
    TFile* f_data_1 = new TFile(DATA_FILE_1.c_str(), "READ");
    if (f_data_1 == NULL) {
        std::cout << "Data file is invalid" << std::endl;
        return 1;
    }
    TFile* f_data_2 = new TFile(DATA_FILE_2.c_str(), "READ");
    if (f_data_2 == NULL) {
        std::cout << "Data file_2 is invalid" << std::endl;
        return 1;
    }
    // Pass the open files to the main RooFitter
    const int RET_CODE = RooFitLifetime(f_data_1, f_data_2, OUT_DIR);

    // Clean up and return the exit code
    delete f_data_1;
    delete f_data_2;

    return RET_CODE;
}

int RooFitLifetime(
        TFile* const DATA_FILE_1,
        TFile* const DATA_FILE_2,
        const std::string& OUT_DIR
        ) {
    // Constants
    const int N_CPU = 8;

    // Set up the variables we're going to read in from the files
    RooRealVar tau_xy("tau_xy", "tau_xy" , -0.3, 5., "ps");
    tau_xy.setRange("range_tau_xy", -0.3, 5.0);

    RooRealVar zjpsi_tau_xy("zjpsi_tau_xy", "zjpsi_tau_xy" , -0.3, 5., "ps");
    zjpsi_tau_xy.setRange("zjpsi_range_tau_xy", -0.3, 5.0);

    TH1D *h_jpsi_tau_xy_very_fine = (TH1D*) DATA_FILE_1->Get("ZFinder/All/jpsi0 tau_xy_very_fine");
    //TH1D *h_jpsi_tau_xy_very_fine = (TH1D*) DATA_FILE_1->Get("ZFinder/All/jpsi0 tau_xy_very_fine pt 10 to 15 GeV");
    //TH1D *h_jpsi_tau_xy_very_fine = (TH1D*) DATA_FILE_1->Get("ZFinder/All/jpsi0 tau_xy_very_fine pt>30");

    TH1D *h_zjpsi_tau_xy_very_fine = (TH1D*) DATA_FILE_2->Get("ZFinder/Jpsi_And_Z/jpsi0 tau_xy_very_fine");
    //TH1D *h_zjpsi_tau_xy_very_fine = (TH1D*) DATA_FILE_2->Get("ZFinder/Jpsi_And_Z/jpsi0 tau_xy_very_fine");

    h_jpsi_tau_xy_very_fine->Sumw2();
    h_jpsi_tau_xy_very_fine->Rebin(2);

    h_zjpsi_tau_xy_very_fine->Sumw2();
    //h_zjpsi_tau_xy_very_fine->Rebin(2);
    h_zjpsi_tau_xy_very_fine->Rebin(10);

    RooDataHist tau_xy_data_hist("tau_xy_data_hist", "tau_xy_data_hist", tau_xy, h_jpsi_tau_xy_very_fine);
    RooDataHist zjpsi_tau_xy_data_hist("zjpsi_tau_xy_data_hist", "zjpsi_tau_xy_data_hist", zjpsi_tau_xy, h_zjpsi_tau_xy_very_fine);

    RooRealVar gaussmean("gaussmean", "Mean of the smearing Gaussian", 0.);
    RooRealVar gausssigma("gausssigma", "Width of the smearing Gaussian", 0.01, 0.0, 0.5);
    //RooRealVar gausssigma("gausssigma", "Width of the smearing Gaussian", 0.093);
    RooGaussModel smear_gauss_model("smear_gauss", "Gaussian used to smear the Exponential Decay", tau_xy, gaussmean, gausssigma);
    RooRealVar decay_lifetime("decay_lifetime", "Tau", 0.5, 0., 3.);
    //RooDecay decay_exp("decay_exp", "Exponential Decay", tau_xy, decay_lifetime, smear_gauss_model, RooDecay::DoubleSided);
    RooDecay decay_exp("decay_exp", "Exponential Decay", tau_xy, decay_lifetime, smear_gauss_model, RooDecay::SingleSided);

    RooRealVar gauss_prompt_mean("gauss_prompt_mean", "Mean of the Prompt Gaussian", 0., -0.2, 0.2);
    RooRealVar gauss_prompt_sigma("gauss_prompt_sigma", "Width of the Prompt Gaussian", 0.005, 0.005, 0.25);
    RooGaussian prompt_gauss("prompt_gauss", "Gaussian of the Prompt Peak", tau_xy, gauss_prompt_mean, gauss_prompt_sigma);

    RooRealVar gauss_prompt_sigma_2("gauss_prompt_sigma_2", "Width of the Prompt Gaussian_2", 0.01, 0.005, 0.4);
    RooGaussian prompt_gauss_2("prompt_gauss_2", "Gaussian_2 of the Prompt Peak", tau_xy, gauss_prompt_mean, gauss_prompt_sigma_2);

    RooRealVar prompt_sharp_fraction("prompt_sharp_fraction", "prompt_sharp_fraction" , 0.3 , 0.0, 1);
    RooAddPdf tau_xy_gauss_sum_fitpdf("tau_xy_gauss_sum_fitpdf", "tau_xy_gauss_sum_fitpdf", RooArgList(prompt_gauss, prompt_gauss_2), RooArgList(prompt_sharp_fraction));

    RooRealVar prompt_fraction("prompt_fraction", "prompt_fraction" , 0.3 , 0.0, 1);
    RooAddPdf tau_xy_fitpdf("tau_xy_fitpdf", "tau_xy_fitpdf", RooArgList(tau_xy_gauss_sum_fitpdf, decay_exp), RooArgList(prompt_fraction));

    tau_xy_fitpdf.fitTo(tau_xy_data_hist, Range("range_tau_xy"), NumCPU(N_CPU), Verbose(false));
    
    double zjpsi_gaussmean_value = gaussmean.getVal();
    double zjpsi_gausssigma_value = gausssigma.getVal();
    double zjpsi_decay_lifetime_value = decay_lifetime.getVal();
    double zjpsi_gauss_prompt_mean_value = gauss_prompt_mean.getVal();
    double zjpsi_gauss_prompt_sigma_value = gauss_prompt_sigma.getVal();
    double zjpsi_gauss_prompt_sigma_2_value = gauss_prompt_sigma_2.getVal();
    double zjpsi_prompt_sharp_fraction_value = prompt_sharp_fraction.getVal();

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
    //TODO letting fractions vary
    //RooRealVar zjpsi_prompt_sharp_fraction("zjpsi_prompt_sharp_fraction", "zjpsi_prompt_sharp_fraction" , zjpsi_prompt_sharp_fraction_value, 0.0, 1.0);
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
    TCanvas* const canvas = get_tcanvas(2000, 750);

    // Plot the left side
    canvas->cd(1);
    gPad->SetLogy();
    //RooPlot* tau_xy_fitframe = tau_xy.frame(-0.3,5.0);
    RooPlot* tau_xy_fitframe = tau_xy.frame(Range("range_tau_xy"));
    tau_xy_fitframe->SetName(0); // Unset title
    tau_xy_data_hist.plotOn(tau_xy_fitframe);
    tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(decay_exp), LineColor(kGreen-2));
    tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(prompt_gauss), LineColor(kBlue-2));
    tau_xy_fitpdf.plotOn(tau_xy_fitframe, Components(prompt_gauss_2), LineColor(kOrange-2));
    tau_xy_fitpdf.plotOn(tau_xy_fitframe, LineColor(kRed-2));
    //tau_xy_fitpdf.plotOn(tau_xy_fitframe, LineColor(kGreen-2), NumCPU(N_CPU));

    tau_xy_fitframe->Draw();

    canvas->Print("Test.png", "png");

    zjpsi_tau_xy_fitpdf.fitTo(zjpsi_tau_xy_data_hist, Range("zjpsi_range_tau_xy"), NumCPU(N_CPU), Verbose(false));
    TCanvas* const zjpsi_canvas = get_tcanvas(2000, 750);
    zjpsi_canvas->cd(1);
    gPad->SetLogy();
    //RooPlot* tau_xy_fitframe = tau_xy.frame(-0.3,5.0);
    RooPlot* zjpsi_tau_xy_fitframe = zjpsi_tau_xy.frame(Range("zjpsi_range_tau_xy"));
    zjpsi_tau_xy_fitframe->SetName(0); // Unset title
    zjpsi_tau_xy_data_hist.plotOn(zjpsi_tau_xy_fitframe);
    zjpsi_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_decay_exp), LineColor(kGreen-2));
    zjpsi_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_prompt_gauss), LineColor(kBlue-2));
    zjpsi_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, Components(zjpsi_prompt_gauss_2), LineColor(kOrange-2));
    zjpsi_tau_xy_fitpdf.plotOn(zjpsi_tau_xy_fitframe, LineColor(kRed-2));

    zjpsi_tau_xy_fitframe->Draw();
    zjpsi_canvas->Print("Test_2.png", "png");

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
        const std::string DATA_FILE_2(argv[2]);
        const std::string OUT_DIR(argv[3]);

        return RooFitLifetime(
                DATA_FILE_1,
                DATA_FILE_2,
                OUT_DIR
                );
    }
}
