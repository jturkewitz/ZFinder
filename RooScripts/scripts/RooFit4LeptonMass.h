#ifndef ROOFIT4LEPTONMASS_H_
#define ROOFIT4LEPTONMASS_H_

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

int RooFit4LeptonMass(
        TFile* const DATA_FILE_1,
        const bool USE_Z_TO_EE,
        const std::string& OUT_DIR
        );

// Main, used only on the command line
int main(int argc, char* argv[]);
#endif // ROOFIT4LEPTONMASS_H_
