#define zfinder_tree_cxx

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "zfinder_tree_muons.h"
//#include "zfinder_tree_electrons.h"
//#include "zfinder_tree_inclusive_jpsi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//TODO refactor this code, for now just make three classes, but better ways to do this exist
//while still avoiding ifdefs, which are difficult to understand and lead to errors
 
//   In a ROOT session, you can do:
//      Root > .L zfinder_tree.C
//      Root > zfinder_tree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

void zfinder_tree::Loop() {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    // if (Cut(ientry) < 0) continue;
  }
}
//void zfinder_tree_electrons::Loop() {
//  if (fChain == 0) return;
//
//  Long64_t nentries = fChain->GetEntriesFast();
//
//  Long64_t nbytes = 0, nb = 0;
//  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//    Long64_t ientry = LoadTree(jentry);
//    if (ientry < 0) break;
//    nb = fChain->GetEntry(jentry);
//    nbytes += nb;
//    // if (Cut(ientry) < 0) continue;
//  }
//}
//void zfinder_tree_jpsi::Loop() {
//  if (fChain == 0) return;
//
//  Long64_t nentries = fChain->GetEntriesFast();
//
//  Long64_t nbytes = 0, nb = 0;
//  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//    Long64_t ientry = LoadTree(jentry);
//    if (ientry < 0) break;
//    nb = fChain->GetEntry(jentry);
//    nbytes += nb;
//    // if (Cut(ientry) < 0) continue;
//  }
//}
