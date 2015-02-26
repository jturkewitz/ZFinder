#ifndef ROOFITZMASS_H_
#define ROOFITZMASS_H_

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

int RooFitZMass(
        TFile* const DATA_FILE_1,
        TFile* const DATA_FILE_2,
        const bool USE_Z_TO_EE,
        const std::string& OUT_DIR
        );

// Main, used only on the command line
int main(int argc, char* argv[]);
#endif // ROOFITZMASS_H_
