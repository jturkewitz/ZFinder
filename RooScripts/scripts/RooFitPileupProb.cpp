// Standard Library
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>  // std::vector

// ROOT
#include <TCanvas.h>
#include "TH1.h"
#include "TRandom3.h"

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
#include "RooRandom.h"


using namespace RooFit;

int RooFitPileupProb() {

  double sigma = 5.;
  double distance_cut = 1.0;

  TRandom3 random_1(0);
  TRandom3 random_2(0);

  int n_within_distance = 0;
  int n_outside_distance = 0;

  for (int i=0 ; i < 1e9 ; ++i) {
    double x1 = random_1.Gaus(0.0,sigma);
    double x2 = random_2.Gaus(0.0,sigma);
    double dist = fabs(x1 - x2);
    if (dist <= distance_cut ) {
      n_within_distance++;
    }
    else {
      n_outside_distance++;
    }
  }
  float prob = (float) n_within_distance / (n_within_distance + n_outside_distance);
  std::cout << "inside: " << n_within_distance << " outside " << n_outside_distance << std::endl;
  std::cout << "prob: " << prob << std::endl;

  return 0;
}

int main(int argc, char* argv[]) {
  const int ARGC = 1;
  if (argc < ARGC) {
    std::cout << "Not enough arguments.";
    return 1;
  } else if (argc > ARGC) {
    std::cout << "Too many arguments.";
    return 1;
  } else {
    /* Read in arguments */
    //const std::string GAUSSIAN_SIGMA(argv[1]);
    //const std::string DISTANCE(argv[2]);
    RooFitPileupProb();
    return 0;
  }
}
