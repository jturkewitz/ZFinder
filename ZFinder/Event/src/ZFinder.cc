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

// Particle Flow
#include "RecoParticleFlow/Configuration/test/PFCandidateAnalyzer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"

//Muons
#include "DataFormats/MuonReco/interface/Muon.h" // reco::Muon
#include "DataFormats/MuonReco/interface/MuonSelectors.h" //muon::isSoftMuon

// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"  // GsfElectron
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"  // GenParticle

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
    zf::ZFinderPlotter *zfp_jpsi_and_z_same_vertex, *zfp_jpsi, *zfp_dimuon_primary_vertex;
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

  //TODO testing set up Z+jpsi plotter, give better names

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
  TFileDirectory tdir_jpsi_and_z_same_vertex(fs->mkdir("Jpsi_And_Z_Same_Vertex"));

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

  //TODO fix the ugliness of the code
  //go to complicated base ZFinder method?
  //dimuon => muon pt cut for now
  //TODO decide on order of cuts?
  //TODO define an integer that corresponds to a cut level to make this code less ugly
  //TODO order of cuts - JPsi mass cut should probably be the last cut made
  //TODO there must be a better way to do this - very brute forced right now
  zfp_all = new zf::ZFinderPlotter(tdir_all);
  zfp_dimuon = new zf::ZFinderPlotter(tdir_dimuon, !USE_MC, APPLY_MUON_MIN_PT);
  zfp_dimuon_soft = new zf::ZFinderPlotter(tdir_dimuon_soft, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS);
  zfp_dimuon_vtx_compatible = new zf::ZFinderPlotter(tdir_dimuon_vtx_compatible, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY);
  zfp_dimuon_primary_vertex = new zf::ZFinderPlotter(tdir_dimuon_primary_vertex, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, !APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);
  zfp_jpsi = new zf::ZFinderPlotter(tdir_jpsi, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_dielectron = new zf::ZFinderPlotter(tdir_dielectron);
  zfp_dimuon_and_dielectron = new zf::ZFinderPlotter(tdir_dimuon_and_dielectron, !USE_MC, APPLY_MUON_MIN_PT );
  zfp_jpsi_and_dielectron = new zf::ZFinderPlotter(tdir_jpsi_and_dielectron, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_jpsi_and_z = new zf::ZFinderPlotter(tdir_jpsi_and_z, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_jpsi_and_z_same_vertex = new zf::ZFinderPlotter(tdir_jpsi_and_z_same_vertex, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW, APPLY_VERTEX_Z_POS_WINDOW);

  zfp_four_muons = new zf::ZFinderPlotter(tdir_four_muons, !USE_MC);
  zfp_jpsi_and_two_muons = new zf::ZFinderPlotter(tdir_jpsi_and_two_muons, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_jpsi_and_two_tight_muons = new zf::ZFinderPlotter(tdir_jpsi_and_two_tight_muons, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);
  zfp_jpsi_and_z_from_muons = new zf::ZFinderPlotter(tdir_jpsi_and_z_from_muons, !USE_MC, APPLY_MUON_MIN_PT, APPLY_SOFT_MUONS, APPLY_DIMUON_VTX_COMPATIBILITY, APPLY_JPSI_MASS_WINDOW);

  zfp_all_mc = new zf::ZFinderPlotter(tdir_all_mc, USE_MC);
  zfp_jpsi_mc = new zf::ZFinderPlotter(tdir_jpsi_mc, USE_MC, APPLY_MUON_MIN_PT, !APPLY_SOFT_MUONS, APPLY_JPSI_MASS_WINDOW);

  //zf::ZFinderPlotter* zfp_jpsi = new zf::ZFinderPlotter(tdir_z_jpsi_dir, false);  // False = do not plot Truth
  /*
  // Setup ZDefinitions and plotters
  zdef_psets_ = iConfig.getUntrackedParameter<std::vector<edm::ParameterSet> >("ZDefinitions");
  for (auto& i_pset : zdef_psets_) {
  std::string name = i_pset.getUntrackedParameter<std::string>("name");
  std::vector<std::string> cuts0 = i_pset.getUntrackedParameter<std::vector<std::string> >("cuts0");
  std::vector<std::string> cuts1 = i_pset.getUntrackedParameter<std::vector<std::string> >("cuts1");
  double min_mz = i_pset.getUntrackedParameter<double>("min_mz");
  double max_mz = i_pset.getUntrackedParameter<double>("max_mz");
  // Make the ZDef
  zf::ZDefinition* zd = new zf::ZDefinition(name, cuts0, cuts1, min_mz, max_mz);
  zdefs_.push_back(zd);
  // Make the Plotter for the ZDef, and the workstations
  // Reco
  TFileDirectory tdir_zd(fs->mkdir(name + " Reco"));
  zf::ZDefinitionPlotter* zdp_reco = new zf::ZDefinitionPlotter(*zd, tdir_zd, false);  // False = do not plot Truth
  zdef_plotters_.push_back(zdp_reco);
  zf::ZDefinitionWorkspace* zdw_reco = new zf::ZDefinitionWorkspace(*zd, tdir_zd, false, true);  // False = do not use Truth
  zdef_workspaces_.push_back(zdw_reco);
  // MC
  TFileDirectory tdir_zd_truth(fs->mkdir(name + " MC"));
  zf::ZDefinitionPlotter* zdp_truth = new zf::ZDefinitionPlotter(*zd, tdir_zd_truth, true);
  zdef_plotters_.push_back(zdp_truth);
  zf::ZDefinitionWorkspace* zdw_truth = new zf::ZDefinitionWorkspace(*zd, tdir_zd_truth, true, true);
  zdef_workspaces_.push_back(zdw_truth);
  }
  */
}

ZFinder::~ZFinder() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // Delete our heap variables
  // TODO decide if this is needed
  //for (auto& i_set : setters_) {
  //  delete i_set;
  //}
  //for (auto& i_zdef : zdefs_) {
  //  delete i_zdef;
  //}
  //for (auto& i_zdefp : zdef_plotters_) {
  //  delete i_zdefp;
  //}
  //for (auto& i_zdefw : zdef_workspaces_) {
  //  delete i_zdefw;
  //}
}


//
// member functions
//

// ------------ method called for each event  ------------
void ZFinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace zf;

  zf::ZFinderEvent zfe(iEvent, iSetup, iConfig_);

  bool found_dimuon = false;
  bool is_soft_dimuon = false;
  bool is_dimuon_vtx_compatible = false;
  bool is_dimuon_primary_vtx_z_compatible = false;
  bool found_jpsi = false;

  //TODO default for this is true, does this make sense?
  bool found_distinct_z_muons = true;


  //TODO also decide if soft muons is a good criterion, maybe rename the bool
  //TODO decide if want behavior of associating with primary vertex when no z is found
  //TODO move this code to ZFinderPlotter? then initialize Fill with bools that are desired?
  //i.e. zf::Plotter::Fill(zfe,JPSI) would fill at the JPSI cut level
  for (unsigned int i = 0; i < zfe.reco_jpsi.m.size() ; ++i ) {
    if ( zfe.mu0.at(i).pt() >= MIN_MUON_PT && zfe.mu1.at(i).pt() >= MIN_MUON_PT ) {
      found_dimuon = true;
      if (muon::isSoftMuon(zfe.mu0.at(i), zfe.reco_vert.primary_vert ) 
          && muon::isSoftMuon(zfe.mu1.at(i), zfe.reco_vert.primary_vert) ) {
        is_soft_dimuon = true;
        if ( zfe.reco_jpsi.vtx_prob.at(i) >= MIN_VERTEX_PROB ) {
          is_dimuon_vtx_compatible = true;
          if ( fabs(zfe.reco_vert.primary_z - zfe.reco_jpsi.vtx_z.at(i)) <= MAX_JPSI_VERTEX_Z_DISPLACEMENT ) {
            is_dimuon_primary_vtx_z_compatible = true;
            if (zfe.reco_jpsi.m.at(i) <= MAX_JPSI_MASS && zfe.reco_jpsi.m.at(i) >= MIN_JPSI_MASS ) {
              found_jpsi = true;
              if (zfe.reco_jpsi.muon0_deltaR_to_z_muons.at(i) <= MIN_DELTAR_DISTINCT_Z_JPSI_MUONS ||
                  zfe.reco_jpsi.muon1_deltaR_to_z_muons.at(i) <= MIN_DELTAR_DISTINCT_Z_JPSI_MUONS ) {
                found_distinct_z_muons = false;
              }
              //TODO verify this is what is desired -- I think this handles case of multiple jpsi, if z candidate
              //matches one jpsi but not the other than it is not distinct
              //else {
              //  found_distinct_z_muons = false;
              //}
            }
          }
        }
      }
    }
  }

  bool found_truth_dimuon = false;
  bool found_truth_jpsi = false;
  if (zfe.truth_jpsi.m.size() > 0 ) {
    for (unsigned int i = 0; i < zfe.truth_jpsi.m.size() ; ++i ) {
      if ( zfe.jpsi_muon0.at(i)->pt() >= MIN_MUON_PT && zfe.jpsi_muon1.at(i)->pt() >= MIN_MUON_PT ) {
        found_truth_dimuon = true;
        if (zfe.truth_jpsi.m.at(i) <= MAX_JPSI_MASS && zfe.truth_jpsi.m.at(i) >= MIN_JPSI_MASS ) {
          found_truth_jpsi = true;
        }
      }
    }
  }


  //TODO fix case of multiple jpsis handling
  //TODO decide if want behavior of associating with primary vertex when no z is found

  //TODO clean this code up
  bool found_dielectron = false;
  bool found_z = false;
  if (zfe.reco_z.m > -1 && zfe.e0 != NULL && zfe.e1 != NULL) {
    if ( zfe.e0->pt >= MIN_ELECTRON_PT && zfe.e1->pt >= MIN_ELECTRON_PT ) {
      found_dielectron = true;
      if (zfe.reco_z.m >= MIN_Z_MASS && zfe.reco_z.m <= MAX_Z_MASS ) {
        found_z = true;
      }
    }
  }

  bool found_four_muons = false;
  if (zfe.n_reco_muons >= 4) {
    found_four_muons = true;
  }
  bool found_z_from_muons = false;
  //TODO implement good_muons_from_z
  //TODO move requirements into ZFinderCuts.h
  bool found_good_muons_from_z = false;
  if (zfe.reco_z_from_muons.m > -1 ) {
    if (zfe.z_muon0.pt() > 20 && zfe.z_muon1.pt() > 20 && (muon::isTightMuon(zfe.z_muon0, zfe.reco_vert.primary_vert ) &&
          muon::isTightMuon(zfe.z_muon1, zfe.reco_vert.primary_vert ) ) ) {
      found_good_muons_from_z = true;
      if (zfe.reco_z_from_muons.m >= MIN_Z_MASS && zfe.reco_z_from_muons.m <= MAX_Z_MASS ) {
        found_z_from_muons = true;
      }
    }
  }


  //Fill histograms
  zfp_all->Fill(zfe);
  zfp_all_mc->Fill(zfe);

  if ( found_truth_dimuon && found_truth_jpsi ) {
    zfp_jpsi_mc->Fill(zfe);
  }

  if (found_four_muons) {
    zfp_four_muons->Fill(zfe);
  }

  //TODO dimuon has minimum muon pT cut - do we want this?
  //TODO is switch easier or better?
  if ( found_dimuon ) {
    zfp_dimuon->Fill(zfe);
    if ( is_soft_dimuon ) {
      zfp_dimuon_soft->Fill(zfe);
      if ( is_dimuon_vtx_compatible ) {
        zfp_dimuon_vtx_compatible->Fill(zfe);
        if (is_dimuon_primary_vtx_z_compatible) {
          zfp_dimuon_primary_vertex->Fill(zfe);
          if (found_jpsi) { 
            zfp_jpsi->Fill(zfe);
            if (found_four_muons && found_distinct_z_muons) {
              zfp_jpsi_and_two_muons->Fill(zfe);
              if (found_good_muons_from_z) {
                zfp_jpsi_and_two_tight_muons->Fill(zfe);
                if (found_z_from_muons) {
                  zfp_jpsi_and_z_from_muons->Fill(zfe);
                }
              }
            }
          }
        }
      }
    }
  }

  //TODO fix how electron pt is accessed as opposed to how muon pt is accessed () vs no ()
  if (found_dielectron) {
    zfp_dielectron->Fill(zfe);
    if (found_dimuon) {
      zfp_dimuon_and_dielectron->Fill(zfe);
      if (is_dimuon_vtx_compatible && is_soft_dimuon && found_jpsi ) {
        zfp_jpsi_and_dielectron->Fill(zfe);
        if (found_z) {
          zfp_jpsi_and_z->Fill(zfe);
          if ( is_dimuon_primary_vtx_z_compatible ) {
            zfp_jpsi_and_z_same_vertex->Fill(zfe);
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
  // Write all ZDef workspaces
  /*
     for (auto& i_zdefw : zdef_workspaces_) {
     i_zdefw->Write();
     }
     */
  //zfp_jpsi->Print("jpsi");
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
