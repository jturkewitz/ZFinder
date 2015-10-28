// -*- C++ -*-
//
// Package:    ZFinder
// Class:      ZFinder
//
/**\class ZFinder ZFinder.cc ZShape/ZFinder/src/ZFinder.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Alexander Gude
//         Created:  Thu Aug  8 15:19:00 CDT 2013
// $Id$
//
//


// system include files
#include <memory>

// standard library files
#include <map>  // std::map
#include <string>  // std::string
#include <utility>  // std::pair
#include <iostream>  // std::cout, std::endl

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CMSSW
#include "FWCore/ServiceRegistry/interface/Service.h" // edm::Service
#include "CommonTools/UtilAlgos/interface/TFileService.h" // TFileService

// ZFinder
#include "ZFinder/Event/interface/ZDefinition.h"  // ZDefinition
#include "ZFinder/Event/interface/ZDefinitionPlotter.h"  // ZDefinitionPlotter
#include "ZFinder/Event/interface/ZDefinitionWorkspace.h"  // ZDefinitionWorkspace
#include "ZFinder/Event/interface/ZFinderEvent.h"  // ZFinderEvent
#include "ZFinder/Event/interface/ZFinderPlotter.h"  // ZFinderPlotter
#include "ZFinder/Event/interface/ZFinderTree.h"  // ZFinderTree
#include "ZFinder/Event/interface/ZFinderCuts.h"  // ZFinderCuts

// HLT information
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//
// class declaration
//

class ZFinder : public edm::EDAnalyzer {
  public:
    explicit ZFinder(const edm::ParameterSet&);
    ~ZFinder();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    //TODO TESTING HLT PRESCALE
    //std::string triggerName ;
    //std::string hltMenuLabel_;
    //HLTConfigProvider hltConfig_;
    //std::pair<int,int> pair2;
    //////////////////////////

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    // ----------member data ---------------------------
    const edm::ParameterSet& iConfig_;
    zf::ZFinderTree *z_tree;

    zf::ZFinderPlotter *zfp_all;
    zf::ZFinderPlotter *zfp_dimuon_jpsi, *zfp_dimuon_jpsi_soft, *zfp_dimuon_jpsi_vtx_compatible, *zfp_dimuon_jpsi_primary_vertex, *zfp_jpsi, *zfp_prompt_jpsi;
    zf::ZFinderPlotter *zfp_dielectron_z, *zfp_dielectron_z_good, *zfp_dielectron_z_good_compatible_vertex, *zfp_z_to_electrons, 
      *zfp_z_to_electrons_and_good_dimuon_jpsi, *zfp_z_to_electrons_and_jpsi, *zfp_z_to_electrons_and_prompt_jpsi;
    zf::ZFinderPlotter *zfp_dimuon_z, *zfp_dimuon_z_good, *zfp_dimuon_z_good_compatible_vertex, *zfp_z_to_muons, 
      *zfp_z_to_muons_and_good_dimuon_jpsi, *zfp_z_to_muons_and_jpsi, *zfp_z_to_muons_and_prompt_jpsi;
    zf::ZFinderPlotter *zfp_all_mc, *zfp_jpsi_mc;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ZFinder::ZFinder(const edm::ParameterSet& iConfig) : iConfig_(iConfig) {
  //now do what ever initialization is needed

  //TODO TESTING HLT PRESCALES ??////////////
  //hltMenuLabel_ = iConfig.getParameter<std::string>("HLTMenuLabel");
  //hltMenuLabel_ = iConfig.getParameter<std::string>("hltTriggerSummaryAOD");
  //hltMenuLabel_ = iConfig.getParameter<std::string>("HLT");
  //std::pair<int,int> prescaleValues(iEvent, iSetup, triggerName);
  /////////////////////////////////////////////////////////

  // Setup plotters
  edm::Service<TFileService> fs;

  TFileDirectory tdir_tree(fs->mkdir("Tree"));

  TFileDirectory tdir_all(fs->mkdir("All"));

  TFileDirectory tdir_dimuon_jpsi(fs->mkdir("Dimuon_Jpsi"));
  TFileDirectory tdir_dimuon_jpsi_soft(fs->mkdir("Dimuon_Jpsi_Soft"));
  TFileDirectory tdir_dimuon_jpsi_vtx_compatible(fs->mkdir("Dimuon_Jpsi_Vertex_Compatible"));
  TFileDirectory tdir_dimuon_jpsi_primary_vertex(fs->mkdir("Dimuon_Jpsi_Primary_Vertex"));
  TFileDirectory tdir_jpsi(fs->mkdir("Jpsi"));
  TFileDirectory tdir_prompt_jpsi(fs->mkdir("Prompt_Jpsi"));

  TFileDirectory tdir_dielectron_z(fs->mkdir("Dielectron_Z"));
  TFileDirectory tdir_dielectron_z_good(fs->mkdir("Dielectron_Z_Good"));
  TFileDirectory tdir_dielectron_z_good_compatible_vertex(fs->mkdir("Dielectron_Z_Good_Compatible_Vertex"));
  TFileDirectory tdir_z_to_electrons(fs->mkdir("Z_To_Electrons"));
  TFileDirectory tdir_z_to_electrons_and_good_dimuon_jpsi(fs->mkdir("Z_To_Electrons_And_Good_Dimuon_Jpsi"));
  TFileDirectory tdir_z_to_electrons_and_jpsi(fs->mkdir("Z_To_Electrons_And_Jpsi"));
  TFileDirectory tdir_z_to_electrons_and_prompt_jpsi(fs->mkdir("Z_To_Electrons_And_Prompt_Jpsi"));

  TFileDirectory tdir_dimuon_z(fs->mkdir("Dimuon_Z"));
  TFileDirectory tdir_dimuon_z_good(fs->mkdir("Dimuon_Z_Good"));
  TFileDirectory tdir_dimuon_z_good_compatible_vertex(fs->mkdir("Dimuon_Z_Good_Compatible_Vertex"));
  TFileDirectory tdir_z_to_muons(fs->mkdir("Z_To_Muons"));
  TFileDirectory tdir_z_to_muons_and_good_dimuon_jpsi(fs->mkdir("Z_To_Muons_And_Good_Dimuon_Jpsi"));
  TFileDirectory tdir_z_to_muons_and_jpsi(fs->mkdir("Z_To_Muons_And_Jpsi"));
  TFileDirectory tdir_z_to_muons_and_prompt_jpsi(fs->mkdir("Z_To_Muons_And_Prompt_Jpsi"));

  TFileDirectory tdir_all_mc(fs->mkdir("MC_All"));
  TFileDirectory tdir_jpsi_mc(fs->mkdir("MC_Jpsi"));

  const bool USE_MC = true;
  const bool APPLY_MUON_MIN_PT = true;
  const bool APPLY_SOFT_MUONS = true;
  const bool APPLY_DIMUON_VTX_COMPATIBILITY = true;
  const bool APPLY_JPSI_MASS_WINDOW = true;
  const bool APPLY_VERTEX_Z_POS_WINDOW = true;
  const bool APPLY_PROMPT_JPSI_WINDOW = true;

  //zf::ZFinderTree* z_tree = new zf::ZFinderTree(*zd_reco, tdir_zd, is_mc_);
  z_tree = new zf::ZFinderTree(tdir_tree, false);
  //z_tuples_.push_back(z_tree);

  //dimuon => muon pt cut for now
  //TODO switch order of cuts? Other improvements?
  zfp_all = new zf::ZFinderPlotter(tdir_all);

  zfp_dimuon_jpsi = new zf::ZFinderPlotter(tdir_dimuon_jpsi, !USE_MC, APPLY_MUON_MIN_PT);
  zfp_dimuon_jpsi_soft = new zf::ZFinderPlotter(tdir_dimuon_jpsi_soft, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS);
  //TODO adding in mass requirement here for acceptance*eff, fix this bit of a kludge
  //zfp_dimuon_jpsi_vtx_compatible = new zf::ZFinderPlotter(tdir_dimuon_jpsi_vtx_compatible, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY);
  zfp_dimuon_jpsi_vtx_compatible = new zf::ZFinderPlotter(tdir_dimuon_jpsi_vtx_compatible, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_dimuon_jpsi_primary_vertex = new zf::ZFinderPlotter(tdir_dimuon_jpsi_primary_vertex, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, !APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);
  //TODO should I apply vertex z pos window here??
  zfp_jpsi = new zf::ZFinderPlotter(tdir_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);
  zfp_prompt_jpsi = new zf::ZFinderPlotter(tdir_prompt_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW, APPLY_PROMPT_JPSI_WINDOW);

  zfp_dielectron_z = new zf::ZFinderPlotter(tdir_dielectron_z);
  zfp_dielectron_z_good = new zf::ZFinderPlotter(tdir_dielectron_z_good);
  zfp_dielectron_z_good_compatible_vertex = new zf::ZFinderPlotter(tdir_dielectron_z_good_compatible_vertex);
  zfp_z_to_electrons = new zf::ZFinderPlotter(tdir_z_to_electrons);
  zfp_z_to_electrons_and_good_dimuon_jpsi = new zf::ZFinderPlotter(tdir_z_to_electrons_and_good_dimuon_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, !APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);
  zfp_z_to_electrons_and_jpsi = new zf::ZFinderPlotter(tdir_z_to_electrons_and_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);
  zfp_z_to_electrons_and_prompt_jpsi = new zf::ZFinderPlotter(tdir_z_to_electrons_and_prompt_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW, APPLY_PROMPT_JPSI_WINDOW);

  zfp_dimuon_z = new zf::ZFinderPlotter(tdir_dimuon_z);
  zfp_dimuon_z_good = new zf::ZFinderPlotter(tdir_dimuon_z_good);
  zfp_dimuon_z_good_compatible_vertex = new zf::ZFinderPlotter(tdir_dimuon_z_good_compatible_vertex);
  zfp_z_to_muons = new zf::ZFinderPlotter(tdir_z_to_muons);
  zfp_z_to_muons_and_good_dimuon_jpsi = new zf::ZFinderPlotter(tdir_z_to_muons_and_good_dimuon_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, !APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);
  zfp_z_to_muons_and_jpsi = new zf::ZFinderPlotter(tdir_z_to_muons_and_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);
  zfp_z_to_muons_and_prompt_jpsi = new zf::ZFinderPlotter(tdir_z_to_muons_and_prompt_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW, APPLY_PROMPT_JPSI_WINDOW);


  //TODO is this needed?

  zfp_all_mc = new zf::ZFinderPlotter(tdir_all_mc, USE_MC);
  zfp_jpsi_mc = new zf::ZFinderPlotter(tdir_jpsi_mc, USE_MC, APPLY_MUON_MIN_PT, !APPLY_SOFT_MUONS, APPLY_JPSI_MASS_WINDOW);
}

ZFinder::~ZFinder() {
  //for (auto& i_zdeft : z_tuples_) {
  //  delete i_zdeft;
  //}
}


//
// member functions
//

// ------------ method called for each event  ------------
void ZFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace zf;


  //TODO TESTING ability to print prescale //
  //// Combined L1T (pair.first) and HLT (pair.second) prescales per HLT path
  ////std::pair<int,int> prescaleValues(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& trigger) const;
  ////hltConfig_.init(iRun,iSetup,processName_,changed)
  ////triggerName = "HLT_Dimuon0_Jpsi_v18";
  //triggerName = "HLT_Dimuon0_Jpsi_v15"; //this one worked! 1 200
  ////pair = hltConfig_.prescaleValues(iEvent, iSetup, triggerName);
  //std::cout << "testA" << std::endl;
  //pair2 = hltConfig_.prescaleValues(iEvent, iSetup, triggerName);
  //std::cout << "testB" << std::endl;
  ////std::cout << "Trigger prescales: " << pair2 << std::endl;
  //std::cout << "Trigger prescales: " << pair2.first << " and second " << pair2.second << std::endl;
  ////pair = HLTConfigProvider::prescaleValues(iEvent, iSetup, triggerName);
  //// any one negative => error in retrieving this (L1T or HLT) prescale
  // ///////////////////////////////////////////////////////////////////

  zf::ZFinderEvent zfe(iEvent, iSetup, iConfig_);

  //Fill Tree
  z_tree->Fill(zfe);


  //Fill histograms
  zfp_all->Fill(zfe);
  zfp_all_mc->Fill(zfe);

  if ( zfe.found_truth_jpsi_with_high_pt_muons ) {
    zfp_jpsi_mc->Fill(zfe);
  }

  //TODO nested if might not be necessary here
  if ( zfe.found_dimuon_jpsi_with_high_pt_muons && zfe.found_dimuon_jpsi_with_muons_in_eta_window ) {
    zfp_dimuon_jpsi->Fill(zfe);
    if ( zfe.found_dimuon_jpsi_with_soft_id_and_high_pt_muons ) {
      zfp_dimuon_jpsi_soft->Fill(zfe);
      if ( zfe.found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex ) {
        zfp_dimuon_jpsi_vtx_compatible->Fill(zfe);
        if (zfe.found_good_dimuon_jpsi_compatible_with_primary_vertex) {
          zfp_dimuon_jpsi_primary_vertex->Fill(zfe);
          if (zfe.found_jpsi) {
            zfp_jpsi->Fill(zfe);
            if (zfe.found_prompt_jpsi) { 
              zfp_prompt_jpsi->Fill(zfe);
            }
          }
        }
      }
    }
  }

  if (zfe.found_high_pt_electrons_from_z) {
    zfp_dielectron_z->Fill(zfe);
    if(zfe.found_good_electrons_from_z ) {
      zfp_dielectron_z_good->Fill(zfe);
      if(zfe.found_dielectron_z_compatible_vertex ) {
        zfp_dielectron_z_good_compatible_vertex->Fill(zfe);
        if(zfe.found_z_to_electrons_mass) {
          zfp_z_to_electrons->Fill(zfe);
          if (zfe.found_good_dimuon_jpsi_compatible_with_primary_vertex) {
            zfp_z_to_electrons_and_good_dimuon_jpsi->Fill(zfe);
            if ( zfe.found_jpsi ) {
              zfp_z_to_electrons_and_jpsi->Fill(zfe);
              if (zfe.found_prompt_jpsi) {
                zfp_z_to_electrons_and_prompt_jpsi->Fill(zfe);
              }
            }
          }
        }
      }
    }
  }

  if (zfe.found_high_pt_muons_from_z) {
    zfp_dimuon_z->Fill(zfe);
    if(zfe.found_good_muons_from_z ) {
      zfp_dimuon_z_good->Fill(zfe);
      if(zfe.found_dimuon_z_compatible_vertex ) {
        zfp_dimuon_z_good_compatible_vertex->Fill(zfe);
        if(zfe.found_z_to_muons_mass) {
          zfp_z_to_muons->Fill(zfe);
          if (zfe.found_good_dimuon_jpsi_compatible_with_primary_vertex) {
            zfp_z_to_muons_and_good_dimuon_jpsi->Fill(zfe);
            if ( zfe.found_jpsi ) {
              zfp_z_to_muons_and_jpsi->Fill(zfe);
              if (zfe.found_prompt_jpsi) {
                zfp_z_to_muons_and_prompt_jpsi->Fill(zfe);
              }
            }
          }
        }
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void ZFinder::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void ZFinder::endJob() {
  //// Since large trees will automatically make new files, we need to get
  //// the current file from each tree and write it, but only if it is new
  //// and hasn't be previously written
  //TFile* file = nullptr;
  //std::vector<TFile*> seen;


  ////for (auto& i_zdeft : z_tuples_) {
  ////  file = i_zdeft->GetCurrentFile();
  ////  // Check if we have seen this file before. If we have not then
  ////  // write the file and add it to the vector
  ////  if (std::find(seen.begin(), seen.end(), file) == seen.end()) {
  ////    file->Write();
  ////    seen.push_back(file);
  ////  }
  ////}
  //file = z_tree->GetCurrentFile();
  //// Check if we have seen this file before. If we have not then
  //// write the file and add it to the vector
  //if (std::find(seen.begin(), seen.end(), file) == seen.end()) {
  //  file->Write();
  //  seen.push_back(file);
  //}
}

// ------------ method called when starting to processes a run  ------------
//void ZFinder::beginRun(edm::Run const&, edm::EventSetup const&) {
void ZFinder::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  //TODO TESTING ??? ///
  //bool changed = true;
  ////string hltMenuLabel2_ = "test";
  ////hltMenuLabel_ = "hltTriggerSummaryAOD";
  //hltMenuLabel_ = "HLT";
  ////if (hltConfig_.init(iRun, iSetup, hltMenuLabel_, changed) ) {
  //if (hltConfig_.init(iRun, iSetup, hltMenuLabel_, changed) ) {
  //
  //}
  /////////////////////////////////////////////////////////////////
}

// ------------ method called when ending the processing of a run  ------------
void ZFinder::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void ZFinder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void ZFinder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ZFinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZFinder);
