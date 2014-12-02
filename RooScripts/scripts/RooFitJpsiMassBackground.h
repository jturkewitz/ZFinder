#ifndef ROOFITJPSIMASSBACKGROUND_H_
#define ROOFITJPSIMASSBACKGROUND_H_

// Standard Library
#include <string>
#include <vector>

// ROOT
#include <TFile.h>

//using namespace RooStats;

int RooFitJpsiMassBackground(
        TFile* const DATA_FILE_1,
        const bool USE_Z_TO_EE,
        const std::string& OUT_DIR
        );

// Main, used only on the command line
int main(int argc, char* argv[]);
#endif // ROOFITJPSIMASSBACKGROUND_H_
