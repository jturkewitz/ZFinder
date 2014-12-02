#ifndef ROOFITLIFETIME_H_
#define ROOFITLIFETIME_H_

// Standard Library
#include <string>
#include <vector>

// ROOT
#include <TFile.h>

//using namespace RooStats;

// Set up the TCanvas
TCanvas* get_tcanvas(const int X_DIM = 1280, const int Y_DIM = 640);

// Multiple ways to call the function in ROOT
//int RooFitter(
//        const std::string& DATA_FILE,
//        const std::string& OUT_DIR
//        );

std::vector<double> RooFitLifetime(
        TFile* const DATA_FILE_1,
        TFile* const DATA_FILE_2,
        const bool USE_Z_TO_EE,
        const std::string& OUT_DIR,
        const int PT_SLICE,
        const bool USE_PT_SLICES,
        const int RAP_SLICE,
        const bool USE_RAP_SLICES
        );

// Main, used only on the command line
int main(int argc, char* argv[]);
#endif // ROOFITLIFETIME_H_
