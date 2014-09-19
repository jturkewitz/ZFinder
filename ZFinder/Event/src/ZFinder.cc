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
#include "ZFinder/Event/interface/AcceptanceSetter.h"  // AcceptanceSetter
#include "ZFinder/Event/interface/SetterBase.h"  // SetterBase
#include "ZFinder/Event/interface/TruthMatchSetter.h"  // TruthMatchSetter
#include "ZFinder/Event/interface/ZDefinition.h"  // ZDefinition
#include "ZFinder/Event/interface/ZDefinitionPlotter.h"  // ZDefinitionPlotter
#include "ZFinder/Event/interface/ZDefinitionWorkspace.h"  // ZDefinitionWorkspace
#include "ZFinder/Event/interface/ZEfficiencies.h" // ZEfficiencies
#include "ZFinder/Event/interface/ZFinderEvent.h"  // ZFinderEvent
#include "ZFinder/Event/interface/ZFinderPlotter.h"  // ZFinderPlotter
#include "ZFinder/Event/interface/ZFinderCuts.h"  // ZFinderCuts

//
// class declaration
//

class ZFinder : public edm::EDAnalyzer {
  public:
    explicit ZFinder(const edm::ParameterSet&);
    ~ZFinder();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


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
    std::map<std::string, zf::ZFinderPlotter*> z_plotter_map_;
    std::vector<zf::SetterBase*> setters_;
    std::vector<edm::ParameterSet> zdef_psets_;
    std::vector<zf::ZDefinition*> zdefs_;
    std::vector<zf::ZDefinitionPlotter*> zdef_plotters_;
    std::vector<zf::ZDefinitionWorkspace*> zdef_workspaces_;
    zf::ZEfficiencies zeffs_;
    zf::ZFinderPlotter *zfp_all, *zfp_dimuon, *zfp_dimuon_soft, *zfp_dimuon_vtx_compatible, *zfp_dielectron, *zfp_dimuon_and_dielectron, *zfp_jpsi_and_dielectron, *zfp_jpsi_and_z;
    zf::ZFinderPlotter *zfp_jpsi, *zfp_dimuon_primary_vertex;
    zf::ZFinderPlotter *zfp_four_muons, *zfp_jpsi_and_two_muons, *zfp_jpsi_and_two_tight_muons, *zfp_jpsi_and_z_from_muons;
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

  // Setup Cut Setters
  zf::AcceptanceSetter* accset = new zf::AcceptanceSetter();
  setters_.push_back(accset);
  zf::TruthMatchSetter* tmset = new zf::TruthMatchSetter();
  setters_.push_back(tmset);

  // Setup plotters
  edm::Service<TFileService> fs;

  TFileDirectory tdir_all(fs->mkdir("All"));
  TFileDirectory tdir_dimuon(fs->mkdir("Dimuon"));
  TFileDirectory tdir_dimuon_soft(fs->mkdir("Dimuon_Soft"));
  TFileDirectory tdir_dimuon_vtx_compatible(fs->mkdir("Dimuon_Vertex_Compatible"));
  TFileDirectory tdir_dimuon_primary_vertex(fs->mkdir("Dimuon_Primary_Vertex"));
  TFileDirectory tdir_jpsi(fs->mkdir("Jpsi"));
  TFileDirectory tdir_dielectron(fs->mkdir("Dielectron"));
  TFileDirectory tdir_dimuon_and_dielectron(fs->mkdir("Dimuon_And_Dielectron"));
  TFileDirectory tdir_jpsi_and_dielectron(fs->mkdir("Jpsi_And_Dielectron"));
  TFileDirectory tdir_jpsi_and_z(fs->mkdir("Jpsi_And_Z"));

  TFileDirectory tdir_four_muons(fs->mkdir("Four_Muons"));
  TFileDirectory tdir_jpsi_and_two_muons(fs->mkdir("Jpsi_And_Two_Muons"));
  TFileDirectory tdir_jpsi_and_two_tight_muons(fs->mkdir("Jpsi_And_Two_Tight_Muons"));
  TFileDirectory tdir_jpsi_and_z_from_muons(fs->mkdir("Jpsi_And_Z_From_Muons"));

  TFileDirectory tdir_all_mc(fs->mkdir("MC_All"));
  TFileDirectory tdir_jpsi_mc(fs->mkdir("MC_Jpsi"));

  // Make our TFileDirectory for the plotter
  //TFileDirectory t_subdir = tdir_zd.mkdir("JPSI", "JPSI");

  //TFileDirectory t_subdir1 = tdir_zd.mkdir("JPSI_Electron", "JPSI_Electron");
  const bool USE_MC = true;
  const bool APPLY_MUON_MIN_PT = true;
  const bool APPLY_SOFT_MUONS = true;
  const bool APPLY_DIMUON_VTX_COMPATIBILITY = true;
  const bool APPLY_JPSI_MASS_WINDOW = true;
  const bool APPLY_VERTEX_Z_POS_WINDOW = true;

  //dimuon => muon pt cut for now
  //TODO switch order of cuts? Other improvements?
  zfp_all = new zf::ZFinderPlotter(tdir_all);
  zfp_dimuon = new zf::ZFinderPlotter(tdir_dimuon, !USE_MC, APPLY_MUON_MIN_PT);
  zfp_dimuon_soft = new zf::ZFinderPlotter(tdir_dimuon_soft, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS);
  zfp_dimuon_vtx_compatible = new zf::ZFinderPlotter(tdir_dimuon_vtx_compatible, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY);
  zfp_dimuon_primary_vertex = new zf::ZFinderPlotter(tdir_dimuon_primary_vertex, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, !APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);
  zfp_jpsi = new zf::ZFinderPlotter(tdir_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_dielectron = new zf::ZFinderPlotter(tdir_dielectron);
  zfp_dimuon_and_dielectron = new zf::ZFinderPlotter(tdir_dimuon_and_dielectron, !USE_MC, APPLY_MUON_MIN_PT );
  zfp_jpsi_and_dielectron = new zf::ZFinderPlotter(tdir_jpsi_and_dielectron, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_jpsi_and_z = new zf::ZFinderPlotter(tdir_jpsi_and_z, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);

  zfp_four_muons = new zf::ZFinderPlotter(tdir_four_muons, !USE_MC);
  zfp_jpsi_and_two_muons = new zf::ZFinderPlotter(tdir_jpsi_and_two_muons, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_jpsi_and_two_tight_muons = new zf::ZFinderPlotter(tdir_jpsi_and_two_tight_muons, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_jpsi_and_z_from_muons = new zf::ZFinderPlotter(tdir_jpsi_and_z_from_muons, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);

  zfp_all_mc = new zf::ZFinderPlotter(tdir_all_mc, USE_MC);
  zfp_jpsi_mc = new zf::ZFinderPlotter(tdir_jpsi_mc, USE_MC, APPLY_MUON_MIN_PT, !APPLY_SOFT_MUONS, APPLY_JPSI_MASS_WINDOW);

}

ZFinder::~ZFinder() {
}


//
// member functions
//

// ------------ method called for each event  ------------
void ZFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace zf;

  zf::ZFinderEvent zfe(iEvent, iSetup, iConfig_);

  //Fill histograms
  zfp_all->Fill(zfe);
  zfp_all_mc->Fill(zfe);

  if ( zfe.found_truth_jpsi_with_high_pt_muons ) {
    zfp_jpsi_mc->Fill(zfe);
  }

  if (zfe.found_four_muons) {
    zfp_four_muons->Fill(zfe);
  }

  //TODO nested if might not be necessary here
  if ( zfe.found_dimuon_with_high_pt_muons ) {
    zfp_dimuon->Fill(zfe);
    if ( zfe.found_dimuon_with_soft_id_and_high_pt_muons ) {
      zfp_dimuon_soft->Fill(zfe);
      if ( zfe.found_dimuon_with_good_muons_and_compatible_muon_vertex ) {
        zfp_dimuon_vtx_compatible->Fill(zfe);
        if (zfe.found_good_dimuon_compatible_with_primary_vertex) {
          zfp_dimuon_primary_vertex->Fill(zfe);
          if (zfe.found_jpsi) {
            zfp_jpsi->Fill(zfe);
            if (zfe.found_four_muons) {
              zfp_jpsi_and_two_muons->Fill(zfe);
              if (zfe.found_good_muons_from_z) {
                zfp_jpsi_and_two_tight_muons->Fill(zfe);
                if (zfe.found_z_to_muons) {
                  zfp_jpsi_and_z_from_muons->Fill(zfe);
                }
              }
            }
          }
        }
      }
    }
  }

  if (zfe.found_good_electrons_from_z) {
    zfp_dielectron->Fill(zfe);
    if (zfe.found_dimuon_with_high_pt_muons) {
      zfp_dimuon_and_dielectron->Fill(zfe);
      if ( zfe.found_jpsi ) {
        zfp_jpsi_and_dielectron->Fill(zfe);
        if (zfe.found_z_to_electrons) {
          zfp_jpsi_and_z->Fill(zfe);
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
}

// ------------ method called when starting to processes a run  ------------
void ZFinder::beginRun(edm::Run const&, edm::EventSetup const&) {
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
