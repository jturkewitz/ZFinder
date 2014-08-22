#ifndef ZFINDER_ZFINDEREVENT_H_
#define ZFINDER_ZFINDEREVENT_H_

// Standard Library
#include <map>  // std::map
#include <string>  // std::string
#include <utility>  // std::pair
#include <vector>  // std::vector

// CMSSW
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"  // reco::GsfElectron
#include "DataFormats/EgammaCandidates/interface/Photon.h"  // reco::Photon
#include "DataFormats/MuonReco/interface/Muon.h" // reco::Muon
#include "DataFormats/JetReco/interface/PFJet.h" //
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"  // reco::GenParticle
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"  // reco::RecoEcalCandidate
#include "FWCore/Framework/interface/Event.h"  // edm::Event, edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"  // edm::InputTag
#include "DataFormats/HLTReco/interface/TriggerObject.h"  // trigger::TriggerObject
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h" // edm::LumiReWeighting

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// ZFinder
#include "ZFinder/Event/interface/ZFinderElectron.h"  // ZFinderElectron, ZFinderElectron
#include "ZFinder/Event/interface/ZFinderMuon.h"  // ZFinderMuon, ZFinderMuon
#include "ZFinder/Event/interface/ZFinderCuts.h"  // ZFinderMuon, ZFinderMuon

namespace zf {

  // Cut level struct
  struct CutLevel{
    // Constructor sets all values to false
    CutLevel() {
      pass = false;
      t0p1_pass = false;
      t0p1_eff = 1.;
      t1p0_pass = false;
      t1p0_eff = 1.;
      event_weight = 1.;
    }
    bool pass;
    bool t0p1_pass;
    bool t1p0_pass;
    double t0p1_eff;
    double t1p0_eff;
    double event_weight;
  };


  // Used to match cut levels to names
  typedef std::pair<std::string, CutLevel> cutlevel_pair;
  // Used to pass around cut levels
  typedef std::vector<cutlevel_pair> cutlevel_vector;
  // Used to pass around trigger objects for matching
  typedef std::pair<const trigger::TriggerObject*, double> trig_dr_pair;
  typedef std::vector<trig_dr_pair> trig_dr_vec;


  class ZFinderEvent{
    public:
      // Constructor. Although iEvent, iSetup, and iConfig violate our naming
      // convention, they are almost ubiquitous in CMSSW code
      ZFinderEvent() {}
      ZFinderEvent(
          const edm::Event& iEvent,
          const edm::EventSetup& iSetup,
          const edm::ParameterSet& iConfig
          );
      // Destructor
      ~ZFinderEvent();

      // Data or MC
      bool is_real_data;

      // Beam Spot
      struct Beamspot{
        double x;
        double y;
        double z;
      } reco_bs;

      // Primary vertexes
      struct Vertexes{
        unsigned int num;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        double primary_x;
        double primary_y;
        double primary_z;
        reco::Vertex primary_vert;
      } truth_vert, reco_vert;

      // Event ID
      struct EventID{
        unsigned int run_num;
        unsigned int lumi_num;
        unsigned int event_num;
      } id;

      //TODO synchronize how muons and electrons are handled?
      // Z Data
      struct ZData{
        double m;
        double pt;
        double y;
        double phi;
        double phistar;
        double eta;
        double vtx_prob;
        double vtx_x;
        double vtx_y;
        double vtx_z;
        TransientVertex vtx;
      } reco_z, truth_z;
      struct ZFromMuonsData{
        std::vector<double> m;
        std::vector<double> pt;
        std::vector<double> y;
        std::vector<double> phi;
        std::vector<double> phistar;
        std::vector<double> eta;
        std::vector<double> vtx_prob;
        std::vector<double> vtx_x;
        std::vector<double> vtx_y;
        std::vector<double> vtx_z;
        std::vector<TransientVertex> vtx;
      } reco_z_from_muons;

      // JPsi Data
      // TODO isolation - higher threshhold of et for neutral hadrons?
      struct JPsiData{
        std::vector<double> m;
        std::vector<double> pt;
        std::vector<double> y;
        std::vector<double> phistar;
        std::vector<double> eta;
        std::vector<double> phi;
        std::vector<double> tau_xy;
        std::vector<double> tau_z;
        std::vector<double> distance_x;
        std::vector<double> distance_y;
        std::vector<double> distance_z;
        std::vector<double> distance;
        std::vector<double> dist_err;
        std::vector<double> chi2;
        std::vector<double> distance_xy;
        std::vector<double> dist_err_xy;
        std::vector<double> chi2_xy;
        std::vector<double> vtx_x;
        std::vector<double> vtx_y;
        std::vector<double> vtx_z;
        std::vector<double> vtx_prob;
        std::vector<double> muons_delta_phi;
        std::vector<double> muons_delta_eta;
        std::vector<double> muons_deltaR;
        std::vector<double> z_delta_phi;
        std::vector<double> iso_mu0;
        std::vector<double> iso_sum_charged_hadron_pt_mu0;
        std::vector<double> iso_sum_charged_particle_pt_mu0;
        std::vector<double> iso_sum_neutral_hadron_et_mu0;
        std::vector<double> iso_sum_photon_et_mu0;
        std::vector<double> iso_sum_pileup_pt_mu0;
        std::vector<double> iso_mu1;
        std::vector<double> iso_sum_charged_hadron_pt_mu1;
        std::vector<double> iso_sum_charged_particle_pt_mu1;
        std::vector<double> iso_sum_neutral_hadron_et_mu1;
        std::vector<double> iso_sum_photon_et_mu1;
        std::vector<double> iso_sum_pileup_pt_mu1;
      } reco_jpsi, truth_jpsi;

      struct JPsiMuonData{
        std::vector<double> deltaR_truth;
      } reco_jpsi_muon0, reco_jpsi_muon1;

      struct Jets{
        std::vector<double> pt;
        std::vector<double> phi;
        std::vector<double> eta;
        std::vector<double> btag_discriminator;
      } reco_jets, reco_muon_jets;
      // Event weight, used for things like pileup reweighting. Most
      // other weights are cut dependent (efficiencies for example) and
      // so are store with the cuts in each electron, with a combined
      // efficiency calculated in ZDefinition.
      double event_weight;


      // These are the special, selected electrons used to make the Z
      ZFinderElectron* e0;
      ZFinderElectron* e1;
      void set_e0(ZFinderElectron* electron) { e0 = electron; }
      void set_e1(ZFinderElectron* electron) { e1 = electron; }
      void set_both_e(ZFinderElectron* electron0, ZFinderElectron* electron1) { e0 = electron0; e1 = electron1; }
      ZFinderElectron* e0_truth;
      ZFinderElectron* e1_truth;
      void set_e0_truth(ZFinderElectron* electron) { e0_truth = electron; }
      void set_e1_truth(ZFinderElectron* electron) { e1_truth = electron; }
      void set_both_e_truth(ZFinderElectron* electron0, ZFinderElectron* electron1) { e0_truth = electron0; e1_truth = electron1; }
      ZFinderElectron* e0_trig;
      ZFinderElectron* e1_trig;
      void set_e0_trig(ZFinderElectron* electron) { e0_trig = electron; }
      void set_e1_trig(ZFinderElectron* electron) { e1_trig = electron; }
      void set_both_e_trig(ZFinderElectron* electron0, ZFinderElectron* electron1) { e0_trig = electron0; e1_trig = electron1; }

      // These are the special, selected muons used to make the JPsi TODO implement/clean up this
      //ZFinderMuon* mu0;
      //ZFinderMuon* mu1;
      std::vector<const reco::Candidate*> jpsi_muon0;
      std::vector<const reco::Candidate*> jpsi_muon1;
      std::vector<reco::Muon> mu0;
      std::vector<reco::Muon> mu1;
      //void set_mu0(ZFinderMuon* muon) { mu0 = muon; }
      //void set_mu1(ZFinderMuon* muon) { mu1 = muon; }
      //void set_both_mu(ZFinderMuon* muon0, ZFinderMuon* muon1) { mu0 = muon0; mu1 = muon1; }
      //void set_both_mu(const reco::Muon &muon0, const reco::Muon &muon1) { mu0 = muon0; mu1 = muon1; }

      // Access pruned lists of the internal electrons
      std::vector<ZFinderElectron*>* FilteredElectrons();
      std::vector<ZFinderElectron*>* AllElectrons() { return FilteredElectrons(); }
      std::vector<ZFinderElectron*>* FilteredElectrons(const std::string& cut_name);

      //std::vector<reco::PFJet> jets;

      // Number of Electrons
      int n_reco_electrons;

      // Number of Anti Electrons
      int n_reco_anti_electrons;

      // Number of Muons
      int n_reco_muons;

      // Number of Jets
      int n_reco_jets;
      int n_reco_muon_jets;

      //reco::TrackRef GetElectronTrackRef(const reco::GsfElectron & e);
      reco::TrackRef GetMuonTrackRef(const reco::Muon & mu);

      //TODO testing
      //KalmanVertexFitter kalman_fitter;

      // Output
      void PrintElectrons(const int TYPE = 0, const bool PRINT_CUTS = false);  // 0 is reco, 1 is truth, 2 is trig
      void PrintTruthElectrons(const bool PRINT_CUTS = false) { PrintElectrons(1, PRINT_CUTS); }
      void PrintRecoElectrons(const bool PRINT_CUTS = false) { PrintElectrons(0, PRINT_CUTS); }
      void PrintTrigElectrons(const bool PRINT_CUTS = false) { PrintElectrons(2, PRINT_CUTS); }

      // Access ZDefinition information
      void AddZDef(const std::string NAME, cutlevel_vector PASS_OBJ) { zdef_map_[NAME] = PASS_OBJ; }
      const cutlevel_vector* GetZDef(const std::string& NAME) const;
      bool ZDefPassed(const std::string& NAME) const;
      void PrintZDefs(const bool VERBOSE = false) const;

    protected:
      // These variables are defined at the top of ZFinderEvent.cc to
      // avoid compilation issues
      static const double TRIG_DR_;

      void SetEventWeight(const edm::Event& iEvent);

      // Called by the constructor to handle MC and Data separately
      void InitReco(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void InitTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void InitTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      void InitGSFElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void InitHFElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void InitNTElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      //TODO clean this up
      //void InitMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      // Update the Z Info from e0, e1
      void InitZ(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      void InitZFromMuons(const reco::Muon& mu0, const reco::Muon& mu1, const TransientVertex &dimuon_vertex);

      // Update the JPsi Info from two muons
      void InitJPsi(const reco::Muon& mu0, const reco::Muon& mu1, const TransientVertex &dimuon_vertex);

      void InitJets(const edm::Event& iEvent, const edm::EventSetup& iSetup) ;

      // Initialize all variables to safe values
      void InitVariables();

      // Input tags
      struct InputTags{
        edm::InputTag ecal_electron;
        edm::InputTag nt_electron;
        edm::InputTag muon;
        edm::InputTag jet;
        edm::InputTag conversion;
        edm::InputTag beamspot;
        edm::InputTag rho_iso;
        edm::InputTag vertex;
        edm::InputTag pileup;
        edm::InputTag generator;
        std::vector<edm::InputTag> iso_vals;
        edm::InputTag hf_electron;
        edm::InputTag hf_clusters;
      } inputtags_;

      // Find matching trigger objects
      const trig_dr_vec* GetMatchedTriggerObjects(
          const edm::Event& iEvent,
          const std::vector<std::string>& trig_names,
          const double ETA, const double PHI, const double DR_CUT
          );
      const trigger::TriggerObject* GetBestMatchedTriggerObject(
          const edm::Event& iEvent,
          const std::vector<std::string>& trig_names,
          const double ETA, const double PHI
          );
      bool TriggerMatch(
          const edm::Event& iEvent,
          const std::vector<std::string>& trig_names,
          const double ETA, const double PHI, const double DR_CUT
          );

      double JpsiMuonTruthMatch(const reco::Muon&);

      // A list of all electrons, split into reco and gen
      std::vector<ZFinderElectron*> reco_electrons_;
      std::vector<ZFinderElectron*> reco_anti_electrons_;
      ZFinderElectron* AddRecoElectron(reco::GsfElectron electron);
      void AddRecoElectron(zf::ZFinderElectron zf_electron);
      ZFinderElectron* AddRecoElectron(reco::RecoEcalCandidate electron);
      ZFinderElectron* AddRecoElectron(reco::Photon electron);

      std::vector<ZFinderElectron*> truth_electrons_;
      ZFinderElectron* AddTruthElectron(reco::GenParticle electron);

      std::vector<ZFinderElectron*> hlt_electrons_;
      ZFinderElectron* AddHLTElectron(trigger::TriggerObject electron);


      // A list of all muons
      std::vector<ZFinderMuon*> reco_muons_;
      //TODO remove this and below comment
      //ZFinderMuon* AddRecoMuon(reco::Muon muon);
      void AddRecoMuon(reco::Muon muon);

      // Calculate phistar
      static double ReturnPhistar(const double& eta0, const double& phi0, const double& eta1, const double& phi1);

      // Sorting functions
      static bool SortByPTHighLowElectron(const ZFinderElectron* e0, const ZFinderElectron* e1) { return (e0->pt > e1->pt); }
      static bool SortByPTHighLowMuon(const ZFinderMuon* mu0, const ZFinderMuon* mu1) { return (mu0->pt > mu1->pt); }

      // Print cuts
      void PrintCuts(ZFinderElectron* zf_elec);

      // Store ZDefinition Information
      std::map<std::string, cutlevel_vector> zdef_map_;
      // Pileup reweighting
      static edm::LumiReWeighting* lumi_weights_;

  };
}  // namespace zf
#endif  // ZFINDER_ZFINDEREVENT_H_
