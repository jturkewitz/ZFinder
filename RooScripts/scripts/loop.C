#include <vector>
#include <iostream>
#include <string>
using namespace std;
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include <TStyle.h>
#include "TROOT.h"
#include <map>
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TLatex.h"
#include "TBranch.h"
#include "TTree.h"
#include "TChain.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TProfile.h"
//#include "/Applications/root/macros/AtlasStyle.C"

#define muons
#define SaveNtuple

TFile *root_file;
#include "zfinder_tree.C"

zfinder_tree *jpsi;

void run(const char *  file){
  SetAtlasStyle();

  root_file = new TFile(file,"READ");
  jpsi = new zfinder_tree((TTree*)root_file->Get("zfinder_tree"));

  double JpsiMass, ZMass, JpsiTau, JpsipT, Jpsirap;
  double JpsiMu1pT, JpsiMu2pT, JpsiMu1eta, JpsiMu2eta;
  double unpolarised, longitudinal, transverse, transverse_pos, transverse_neg; 
  double is_z_to_electrons, is_z_to_muons;
  double reco_muon0_weight, reco_muon1_weight;

  int jpsi_entries = (jpsi->fChain)->GetEntries();

  // Load file with acceptance weights
  TFile mainfile_acc( "../AcceptanceMaps/acceptancejpsi_RapMax2.1_PtMax100_RapBins100_PtBins400_Sampling10k.root" );
  TH2* h_weights_acc1 = (TH2*)mainfile_acc.Get("flat");
  TH2* h_weights_acc2 = (TH2*)mainfile_acc.Get("long");
  TH2* h_weights_acc3 = (TH2*)mainfile_acc.Get("trp0");
  TH2* h_weights_acc4 = (TH2*)mainfile_acc.Get("trpp");
  TH2* h_weights_acc5 = (TH2*)mainfile_acc.Get("trpm");

  // Load file with reconstruction weights
  TFile recofile( "../AcceptanceMaps/MuonEfficiencies_Jpsi_13July2014_run2012ABCD_53X.root" );
  TGraphAsymmErrors* h_weights_eta1 = (TGraphAsymmErrors*)recofile.Get("DATA_Loose_pt_abseta<0.9"); 
  TGraphAsymmErrors* h_weights_eta2 = (TGraphAsymmErrors*)recofile.Get("DATA_Loose_pt_abseta0.9-1.2");
  TGraphAsymmErrors* h_weights_eta3 = (TGraphAsymmErrors*)recofile.Get("DATA_Loose_pt_abseta1.2-2.1");

  double *reco_eta1_x = h_weights_eta1->GetX(); double *reco_eta2_x = h_weights_eta2->GetX(); double *reco_eta3_x = h_weights_eta3->GetX(); 
  double *reco_eta1_y = h_weights_eta1->GetY(); double *reco_eta2_y = h_weights_eta2->GetY(); double *reco_eta3_y = h_weights_eta3->GetY(); 

  //  TGraphAsymmErrors* = (TGraphAsymmErrors*)recofile.Get("DATA_Soft_pt_abseta<0.9");    
  //  TGraphAsymmErrors* = (TGraphAsymmErrors*)recofile.Get("DATA_Soft_pt_abseta0.9-1.2"); 
  //  TGraphAsymmErrors* = (TGraphAsymmErrors*)recofile.Get("DATA_Soft_pt_abseta1.2-2.1");
  //  TGraphAsymmErrors* = (TGraphAsymmErrors*)recofile.Get("DATA_Tight_pt_abseta<0.9");
  //  TGraphAsymmErrors* = (TGraphAsymmErrors*)recofile.Get("DATA_Tight_pt_abseta0.9-1.2");
  //  TGraphAsymmErrors* = (TGraphAsymmErrors*)recofile.Get("DATA_Tight_pt_abseta1.2-2.1");

#ifdef muons
  TFile *ntuple = new TFile("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/ntuples_looser/ZJpsimumu.root", "RECREATE");
#else
  TFile *ntuple = new TFile("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/ntuples_looser/ZJpsiee.root", "RECREATE");
#endif
  TTree *aux;
  aux = new TTree("AUX", "AUX");
  aux->Branch("z_mass", &z_mass);
  aux->Branch("onia_mass", &onia_mass);
  aux->Branch("onia_tau", &onia_tau);
  aux->Branch("onia_pt", &onia_pt);
  aux->Branch("onia_rap", &onia_rap);
  aux->Branch("onia_mu0_pt", &onia_mu0_pt); 
  aux->Branch("onia_mu1_pt", &onia_mu1_pt);
  aux->Branch("onia_mu0_eta",&onia_mu0_eta);
  aux->Branch("onia_mu1_eta",&onia_mu1_eta);
  aux->Branch("is_z_to_electrons", &is_z_to_electrons);
  aux->Branch("is_z_to_muons", &is_z_to_muons);
  aux->Branch("unpolarised", &unpolarised);
  aux->Branch("longitudinal",&longitudinal);
  aux->Branch("transverse",  &transverse);
  aux->Branch("transverse_pos", &transverse_pos);
  aux->Branch("transverse_neg", &transverse_neg);
  aux->Branch("reco_muon0_weight", &reco_muon0_weight);
  aux->Branch("reco_muon1_weight", &reco_muon1_weight);

  // find the events with candidates
  for(int iEntry=0; iEntry<(jpsi_entries); iEntry++) {
    if(iEntry%200==0)
      cout << "\r" << (double)iEntry/(double)jpsi_entries*100 << "\% processed" << flush;

    (jpsi->fChain)->GetEntry(iEntry);

    onia_mass = jpsi->reco_jpsi_jpsi_m;
    onia_tau  = jpsi->reco_jpsi_jpsi_tau_xy*1000.;
    onia_pt   = jpsi->reco_jpsi_jpsi_pt;
    onia_rap  = jpsi->reco_jpsi_jpsi_eta;
    onia_mu0_pt  = jpsi->reco_jpsi_muon0_pt;
    onia_mu0_eta = jpsi->reco_jpsi_muon0_eta;
    onia_mu1_pt  = jpsi->reco_jpsi_muon1_pt;
    onia_mu1_eta = jpsi->reco_jpsi_muon1_eta;

#ifdef muons
    z_mass = jpsi->reco_z_from_muons_z_m;
    is_z_to_electrons = 0;
    is_z_to_muons= 1;
    if (jpsi->event_info_found_high_pt_muons_from_z==0)
      continue;
    if (jpsi->event_info_found_good_muons_from_z==0)
      continue;
    if (jpsi->event_info_found_dimuon_z_compatible_vertex==0)
      continue;
    if (jpsi->event_info_found_z_to_muons_mass==0)
      continue;
#else
    z_mass = jpsi->reco_z_z_m;
    is_z_to_electrons = 1;
    is_z_to_muons = 0;
    if (jpsi->event_info_found_high_pt_electrons_from_z==0)
      continue;
    if (jpsi->event_info_found_good_electrons_from_z==0)
      continue;
    if (jpsi->event_info_found_dielectron_z_compatible_vertex==0)
      continue;
    if (jpsi->event_info_found_z_to_electrons_mass==0)
      continue;
#endif
    if (jpsi->event_info_found_dimuon_jpsi_with_muons_in_eta_window==0)
      continue;
    if (jpsi->event_info_found_dimuon_jpsi_with_high_pt_muons==0)
      continue;
    if (jpsi->event_info_found_dimuon_jpsi_with_soft_id_and_high_pt_muons==0)
      continue;
    if (jpsi->event_info_found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex==0)
      continue;
    if (jpsi->event_info_found_good_dimuon_jpsi_compatible_with_primary_vertex==0)
      continue;
    if (jpsi->event_info_found_jpsi==0)
      continue;

    // find acceptance weights
    int binx, biny;
    binx = h_weights_acc1->GetXaxis()->FindBin(onia_rap);
    biny = h_weights_acc1->GetYaxis()->FindBin(onia_pt);
    double weight = 1./h_weights_acc1->GetBinContent(binx, biny);
    if (weight>-100 && weight<100.)
      unpolarised = weight;
    else
      unpolarised = 1.;

    binx = h_weights_acc2->GetXaxis()->FindBin(onia_rap);
    biny = h_weights_acc2->GetYaxis()->FindBin(onia_pt);
    weight = 1./h_weights_acc2->GetBinContent(binx, biny);
    if (weight>-100 && weight<100.)
      longitudinal = weight;
    else
      longitudinal = 1.;

    binx = h_weights_acc3->GetXaxis()->FindBin(onia_rap);
    biny = h_weights_acc3->GetYaxis()->FindBin(onia_pt);
    weight = 1./h_weights_acc3->GetBinContent(binx, biny);
    if (weight>-100 && weight<100.)
      transverse = weight;
    else
      transverse = 1.;


    binx = h_weights_acc4->GetXaxis()->FindBin(onia_rap);
    biny = h_weights_acc4->GetYaxis()->FindBin(onia_pt);
    weight = 1./h_weights_acc4->GetBinContent(binx, biny);
    if (weight>-100 && weight<100.)
      transverse_pos = weight;
    else
      transverse_pos = 1.;


    binx = h_weights_acc5->GetXaxis()->FindBin(onia_rap);
    biny = h_weights_acc5->GetYaxis()->FindBin(onia_pt);
    weight = 1./h_weights_acc5->GetBinContent(binx, biny);
    if (weight>-100 && weight<100.)
      transverse_neg = weight;
    else
      transverse_neg = 1.;

    double xvalue, yvalue;

    // find reco weights
    if (TMath::Abs(onia_mu0_eta)<2.1 && onia_mu0_pt<20.) {
      if (TMath::Abs(onia_mu1_eta)<0.9)
        for (int i = 0; i<h_weights_eta1->GetN(); i++) {
          h_weights_eta1->GetPoint(i, xvalue, yvalue);
          if (JpsiMu1pT >xvalue)
            reco_muon0_weight = yvalue;
        }
      if (TMath::Abs(onia_mu0_eta)>0.9 && TMath::Abs(onia_mu0_eta)<1.2)
        for (int i = 0; i<h_weights_eta2->GetN(); i++) {
          h_weights_eta2->GetPoint(i, xvalue, yvalue);
          if (onia_mu0_pt>xvalue)
            reco_muon0_weight = yvalue;
        }
      if (TMath::Abs(onia_mu0_eta)>1.2)
        for (int i = 0; i<h_weights_eta3->GetN(); i++) {
          h_weights_eta3->GetPoint(i, xvalue, yvalue);
          if (onia_mu0_pt>xvalue)
            reco_muon0_weight = yvalue;
        }
    } else 
      reco_muon0_weight =1.;

    if (TMath::Abs(onia_mu1_eta)<2.1 && onia_mu1_pt<20.) {
      if (TMath::Abs(onia_mu1_eta)<0.9)
        for (int i = 0; i<h_weights_eta1->GetN(); i++) {
          h_weights_eta1->GetPoint(i, xvalue, yvalue);
          if (onia_mu1_pt>xvalue)
            reco_muon1_weight = yvalue;
        }
      if (TMath::Abs(onia_mu1_eta)>0.9 && TMath::Abs(onia_mu1_eta)<1.2)
        for (int i = 0; i<h_weights_eta2->GetN(); i++) {
          h_weights_eta2->GetPoint(i, xvalue, yvalue);
          if (onia_mu1_pt>xvalue)
            reco_muon1_weight = yvalue;
        }
      if (TMath::Abs(onia_mu1_eta)>1.2)
        for (int i = 0; i<h_weights_eta3->GetN(); i++) {
          h_weights_eta3->GetPoint(i, xvalue, yvalue);
          if (onia_mu1_pt>xvalue)
            reco_muon1_weight = yvalue;
        }
    } else 
      reco_muon1_weight =1.;


    aux->Fill();

  }

  // calculate the four lepton invariant mass
  cout << endl;

  ntuple->Write();
  ntuple->Close();

}
