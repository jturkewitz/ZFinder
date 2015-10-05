#ifndef ROOFITCOMBINED_H_
#define ROOFITCOMBINED_H_

// Standard Library
#include <string>
#include <vector>

// ROOT
#include <TFile.h>

//using namespace RooStats;

std::vector<double> RooFitCombined(
        TFile* const DATA_FILE_1,
        TFile* const DATA_FILE_2,
        TFile* const DATA_FILE_3,
        const std::string& OUT_DIR,
        const int PT_SLICE
        );

// Main, used only on the command line
int main(int argc, char* argv[]);
#endif // ROOFITCOMBINED_H_
