#ifndef ROOFITPILEUPESTIMATION_H_
#define ROOFITPILEUPESTIMATION_H_

// Standard Library
#include <string>
#include <vector>

// ROOT
#include <TFile.h>

//using namespace RooStats;

// Multiple ways to call the function in ROOT
//int RooFitter(
//        const std::string& DATA_FILE,
//        const std::string& OUT_DIR
//        );

int RooFitPileupEstimation(
        TFile* const DATA_FILE_1,
        const std::string& OUT_DIR
        );

// Main, used only on the command line
int main(int argc, char* argv[]);
#endif // ROOFITLIFETIME_H_
