#ifndef ROOFITLIFETIMEANDMASSSLICES_H_
#define ROOFITLIFETIMEANDMASSSLICES_H_

// Standard Library
#include <string>
#include <vector>

// ROOT
#include <TFile.h>

//using namespace RooStats;

std::vector<double> RooFitLifetimeAndMassSlices(
        TFile* const DATA_FILE_1,
        TFile* const DATA_FILE_2,
        const bool USE_Z_TO_EE,
        const std::string& OUT_DIR,
        const int PT_SLICE
        );

// Main, used only on the command line
int main(int argc, char* argv[]);
#endif // ROOFITLIFETIMEANDMASSSLICES_H_
