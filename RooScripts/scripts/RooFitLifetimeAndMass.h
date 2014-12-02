#ifndef ROOFITLIFETIMEANDMASS_H_
#define ROOFITLIFETIMEANDMASS_H_

// Standard Library
#include <string>
#include <vector>

// ROOT
#include <TFile.h>

//using namespace RooStats;

int RooFitLifetimeAndMass(
        TFile* const DATA_FILE_1,
        TFile* const DATA_FILE_2,
        const bool USE_Z_TO_EE,
        const std::string& OUT_DIR
        );

// Main, used only on the command line
int main(int argc, char* argv[]);
#endif // ROOFITLIFETIMEANDMASS_H_
