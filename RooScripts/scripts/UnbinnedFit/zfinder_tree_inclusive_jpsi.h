//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 20 11:09:47 2015 by ROOT version 5.34/11
// from TTree zfinder_tree/zfinder_tree
// found on file: ntuples/jpsiTest422_tree_z40DoubleMuon2_5_eta24_2012.root
//////////////////////////////////////////////////////////

#ifndef zfinder_tree_jpsi_h
#define zfinder_tree_jpsi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class zfinder_tree_jpsi {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Double_t        reco_z_z_m;
    Double_t        reco_z_z_pt;
    Double_t        reco_z_z_y;
    Double_t        reco_z_z_phi;
    Double_t        reco_z_z_phistar;
    Double_t        reco_z_z_eta;
    Double_t        reco_z_z_vtx_prob;
    Double_t        reco_z_z_vtx_x;
    Double_t        reco_z_z_vtx_y;
    Double_t        reco_z_z_vtx_z;
    Double_t        reco_z_daughter0_pt;
    Double_t        reco_z_daughter0_eta;
    Double_t        reco_z_daughter0_phi;
    Double_t        reco_z_daughter1_pt;
    Double_t        reco_z_daughter1_eta;
    Double_t        reco_z_daughter1_phi;
    Int_t           reco_z_daughter0_charge;
    Int_t           reco_z_daughter1_charge;
    Double_t        reco_z_from_muons_z_m;
    Double_t        reco_z_from_muons_z_pt;
    Double_t        reco_z_from_muons_z_y;
    Double_t        reco_z_from_muons_z_phi;
    Double_t        reco_z_from_muons_z_phistar;
    Double_t        reco_z_from_muons_z_eta;
    Double_t        reco_z_from_muons_z_vtx_prob;
    Double_t        reco_z_from_muons_z_vtx_x;
    Double_t        reco_z_from_muons_z_vtx_y;
    Double_t        reco_z_from_muons_z_vtx_z;
    Double_t        reco_z_from_muons_daughter0_pt;
    Double_t        reco_z_from_muons_daughter0_eta;
    Double_t        reco_z_from_muons_daughter0_phi;
    Double_t        reco_z_from_muons_daughter1_pt;
    Double_t        reco_z_from_muons_daughter1_eta;
    Double_t        reco_z_from_muons_daughter1_phi;
    Int_t           reco_z_from_muons_daughter0_charge;
    Int_t           reco_z_from_muons_daughter1_charge;
    Double_t        reco_jpsi_jpsi_m;
    Double_t        reco_jpsi_jpsi_pt;
    Double_t        reco_jpsi_jpsi_y;
    Double_t        reco_jpsi_jpsi_phi;
    Double_t        reco_jpsi_jpsi_eta;
    Double_t        reco_jpsi_jpsi_vtx_prob;
    Double_t        reco_jpsi_jpsi_vtx_x;
    Double_t        reco_jpsi_jpsi_vtx_y;
    Double_t        reco_jpsi_jpsi_vtx_z;
    Double_t        reco_jpsi_jpsi_tau_xy;
    Double_t        reco_jpsi_jpsi_tau_z;
    Double_t        reco_jpsi_jpsi_distance_xy;
    Double_t        reco_jpsi_jpsi_distance_z;
    Double_t        reco_jpsi_jpsi_eff;
    Double_t        reco_jpsi_jpsi_acc_eff;
    Double_t        reco_jpsi_jpsi_scale_factor;
    Double_t        reco_jpsi_muon0_pt;
    Double_t        reco_jpsi_muon0_eta;
    Double_t        reco_jpsi_muon0_phi;
    Double_t        reco_jpsi_muon1_pt;
    Double_t        reco_jpsi_muon1_eta;
    Double_t        reco_jpsi_muon1_phi;
    Int_t           reco_jpsi_muon0_charge;
    Int_t           reco_jpsi_muon1_charge;
    Int_t           event_info_event_number;
    Int_t           event_info_run_number;
    Int_t           event_info_n_verts;
    Bool_t          event_info_is_mc;
    Bool_t          event_info_found_high_pt_muons_from_z;
    Bool_t          event_info_found_good_muons_from_z;
    Bool_t          event_info_found_dimuon_z_compatible_vertex;
    Bool_t          event_info_found_z_to_muons_mass;
    Bool_t          event_info_found_high_pt_electrons_from_z;
    Bool_t          event_info_found_good_electrons_from_z;
    Bool_t          event_info_found_dielectron_z_compatible_vertex;
    Bool_t          event_info_found_z_to_electrons_mass;
    Bool_t          event_info_found_dimuon_jpsi_with_muons_in_eta_window;
    Bool_t          event_info_found_dimuon_jpsi_with_high_pt_muons;
    Bool_t          event_info_found_dimuon_jpsi_with_soft_id_and_high_pt_muons;
    Bool_t          event_info_found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex;
    Bool_t          event_info_found_good_dimuon_jpsi_compatible_with_primary_vertex;
    Bool_t          event_info_found_jpsi;

    // List of branches
    TBranch        *b_reco_z;   //!
    TBranch        *b_reco_z_from_muons;   //!
    TBranch        *b_reco_jpsi;   //!
    TBranch        *b_event_info;   //!

    zfinder_tree_jpsi(TTree *tree=0);
    virtual ~zfinder_tree_jpsi();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef zfinder_tree_jpsi_cxx
zfinder_tree_jpsi::zfinder_tree_jpsi(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/ntuples_looser/jpsiTest451_MuOniaPartial2012B_dimuon8_eta24_ptsub25.root");
    if (!f || !f->IsOpen()) {jpsiTest422_tree_z40DoubleMuon2_5_eta24_2012.root
      //f = new TFile("ntuples_looser/jpsiTest425_tree_dimuon8_MuOnia2012B_Partial.root");
      f = new TFile("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/ntuples_looser/jpsiTest451_MuOniaPartial2012B_dimuon8_eta24_ptsub25.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/ntuples_looser/jpsiTest451_MuOniaPartial2012B_dimuon8_eta24_ptsub25.root:/ZFinder");
    dir->GetObject("zfinder_tree",tree);

  }
  Init(tree);
}

zfinder_tree_jpsi::~zfinder_tree_jpsi()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t zfinder_tree_jpsi::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t zfinder_tree_jpsi::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void zfinder_tree_jpsi::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("reco_z", &reco_z_z_m, &b_reco_z);
  fChain->SetBranchAddress("reco_z_from_muons", &reco_z_from_muons_z_m, &b_reco_z_from_muons);
  fChain->SetBranchAddress("reco_jpsi", &reco_jpsi_jpsi_m, &b_reco_jpsi);
  fChain->SetBranchAddress("event_info", &event_info_event_number, &b_event_info);
  Notify();
}

Bool_t zfinder_tree_jpsi::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void zfinder_tree_jpsi::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t zfinder_tree_jpsi::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
//#endif // #ifdef zfinder_tree_jpsi_cxx
