#include "ZFinder/Event/interface/ZFinderEvent.h"

// Standard Library
#include <algorithm>  // std::sort, std::swap
#include <iostream>  // std::cout, std::endl

//TODO clean up unneeded includes
// CMSSW
#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"  // reco::PhotonCollection
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"  // reco::HFEMClusterShape
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"  // reco::HFEMClusterShapeAssociationCollection
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeFwd.h"  // reco::HFEMClusterShapeRef,
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"  // reco::SuperClusterCollection, reco::SuperClusterRef
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"  // reco::RecoEcalCandidateCollection
#include "DataFormats/MuonReco/interface/MuonFwd.h" // reco::MuonCollection
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"  // EgammaCutBasedEleId::PassWP, EgammaCutBasedEleId::*
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"  // PileupSummaryInfo
#include "DataFormats/HLTReco/interface/TriggerEvent.h" // trigger::TriggerEvent
#include "DataFormats/TrackReco/interface/Track.h" //reco::Track
#include "DataFormats/JetReco/interface/PFJetCollection.h" //
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

// for vertexing                                                                                                                                                                                        
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

//for b-tagging
#include "DataFormats/BTauReco/interface/JetTag.h"

// ZFinder
#include "ZFinder/Event/interface/PDGID.h"  // PDGID enum (ELECTRON, POSITRON, etc.)
#include "ZFinder/Event/interface/TriggerList.h"  // ET_ET_TIGHT, ET_ET_DZ, ET_ET_LOOSE, ET_NT_ET_TIGHT, ET_HF_ET_TIGHT, ET_HF_ET_LOOSE, ET_HF_HF_TIGHT, ET_HF_HF_LOOSE, SINGLE_ELECTRON_TRIGGER, ALL_TRIGGERS

//Math
#include "Math/GenVector/VectorUtil.h"
#include <math.h>
#include <TMath.h>

#include <TVector.h>
#include <TMatrix.h>

namespace zf {
  /*
   * These variables are hard coded here for easy access, instead of randomly
   * scattering them throughout the code
   */
  // Electrons are considered matched to a trigger object if close than this
  // value
  const double ZFinderEvent::TRIG_DR_ = 0.3;
  const double MIN_MUON_PT = 3.0; //TODO probably should be 5 GeV based on efficiency plot of tight muons, but is eta dependent
  const double MIN_VERTEX_PROB = 0.005; //from double jpsi paper (I think - TODO verify this)

  ZFinderEvent::ZFinderEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::ParameterSet& iConfig) {
    /* Given an event, parses them for the information needed to make the
     * classe.
     *
     * It selects electrons based on a minimum level of hard-coded cuts.
     */
    // Clear Events
    InitVariables();

    // Get event info
    id.run_num = iEvent.run();
    id.lumi_num = iEvent.luminosityBlock();
    id.event_num = iEvent.id().event();

    // Set local is_real_data
    is_real_data = iEvent.isRealData();

    // Get InputTags
    // Reco
    inputtags_.ecal_electron = iConfig.getParameter<edm::InputTag>("ecalElectronsInputTag");
    inputtags_.nt_electron = iConfig.getParameter<edm::InputTag>("ntElectronsInputTag");
    inputtags_.hf_electron = iConfig.getParameter<edm::InputTag>("hfElectronsInputTag");
    inputtags_.muon = iConfig.getParameter<edm::InputTag>("muonsInputTag");
    inputtags_.jet = iConfig.getParameter<edm::InputTag>("ak5PFJetsInputTag");
    inputtags_.hf_clusters = iConfig.getParameter<edm::InputTag>("hfClustersInputTag");
    inputtags_.conversion = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
    inputtags_.beamspot = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
    inputtags_.rho_iso = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
    inputtags_.vertex = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    inputtags_.iso_vals = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");

    // Truth
    inputtags_.pileup = iConfig.getParameter<edm::InputTag>("pileupInputTag");
    inputtags_.generator = iConfig.getParameter<edm::InputTag>("generatorInputTag");

    // Finish initialization of electrons
    InitReco(iEvent, iSetup);  // Data
    if (!is_real_data) {
      InitTruth(iEvent, iSetup);  // MC
    } else {
      InitTrigger(iEvent, iSetup);  // Trigger Matching
    }
  }

  void ZFinderEvent::InitReco(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    /* Count Pile Up and store first vertex location*/
    edm::Handle<reco::VertexCollection> reco_vertices;
    iEvent.getByLabel(inputtags_.vertex, reco_vertices);
    reco_vert.num = 0;
    bool first_vertex = true;
    reco::Vertex primary_vertex;
    for (unsigned int vertex=0; vertex < reco_vertices->size(); ++vertex) {
      if (    // Criteria copied from twiki
          !((*reco_vertices)[vertex].isFake())
          && ((*reco_vertices)[vertex].ndof() > 4)
          && (fabs((*reco_vertices)[vertex].z()) <= 24.0)
          && ((*reco_vertices)[vertex].position().Rho() <= 2.0)
         ) {
        reco_vert.num++;
        reco_vert.x.push_back( (*reco_vertices)[vertex].x() );
        reco_vert.y.push_back( (*reco_vertices)[vertex].y() );
        reco_vert.z.push_back( (*reco_vertices)[vertex].z() );
        // Store first good vertex as "primary"
        // TODO verify that first_vertex is highest pt vertex
        if (first_vertex) {
          first_vertex = false;
          primary_vertex = (*reco_vertices)[vertex];
        }
      }
    }
    reco_vert.primary_x = primary_vertex.position().x();
    reco_vert.primary_y = primary_vertex.position().y();
    reco_vert.primary_z = primary_vertex.position().z();
    reco_vert.primary_vert = primary_vertex;


    /* Beamspot */
    edm::Handle<reco::BeamSpot> beam_spot;
    iEvent.getByLabel(inputtags_.beamspot, beam_spot);
    reco_bs.x = beam_spot->position().X();
    reco_bs.y = beam_spot->position().Y();
    reco_bs.z = beam_spot->position().Z();

    /* Find electrons */
    InitGSFElectrons(iEvent, iSetup);
    //TODO verify this/clean up code
    //for now don't use any non GSF electrons - need tracker to determine vertex
    //InitHFElectrons(iEvent, iSetup);
    //InitNTElectrons(iEvent, iSetup);

    /* Find muons */
    //InitMuons(iEvent, iSetup);

    // Sort our electrons and set e0, e1 as the two with the highest pt
    std::sort(reco_electrons_.begin(), reco_electrons_.end(), SortByPTHighLowElectron);
    std::sort(reco_anti_electrons_.begin(), reco_anti_electrons_.end(), SortByPTHighLowElectron);

    // For Zs
    // TODO clean up vertex code (combine into a function to avoid duplicating code)


    n_reco_electrons = reco_electrons_.size();
    n_reco_anti_electrons = reco_anti_electrons_.size();
    //TODO clean up this code
    //if (n_reco_electrons >= 2) {
    if (n_reco_electrons >= 1 && n_reco_anti_electrons >=1) {
      // Set our internal electrons
      //set_both_e(reco_electrons_[0], reco_electrons_[1]);
      //highest pT electron, and highest pT antielectron
      set_both_e(reco_electrons_[0], reco_anti_electrons_[0]);
      // Set up the Z
      InitZ(iEvent, iSetup);
    }

    //TODO clean this up
    //std::sort(reco_muons_.begin(), reco_muons_.end(), SortByPTHighLowMuon);
    //n_reco_muons = reco_muons_.size();
    //n_reco_muons = reco_muons_.size();
    // For JPsis
    //TODO clean and improve this code

    edm::ESHandle<TransientTrackBuilder> track_builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder);
    std::vector<reco::TransientTrack> transient_tracks;
    //testing - for now don't use ZFinderMuon class

    edm::Handle<reco::MuonCollection> muons_h;
    iEvent.getByLabel(inputtags_.muon, muons_h);
    n_reco_muons = muons_h->size();
    //TODO clean this up
    if (n_reco_muons >= 2 ) {
      //TODO maybe keep highest pt muons - for now keep j/psi muons
      //set_both_mu(reco_muons_[0], reco_muons_[1]);
      //TODO test how list is ordered

      // Set up the JPsi
      // TODO modify this s.t. it takes muons/antimuons
      for ( int i=0 ; i < (n_reco_muons - 1) ; ++i ) {
        const reco::Muon muon0 = muons_h->at(i);
        if ( muon0.pt() < MIN_MUON_PT ) {
          continue;
        }
        if (!muon::isSoftMuon(muon0, primary_vertex)) {
          continue;
        }
        for ( int j=i+1 ; j < n_reco_muons ; ++j ) {
          const reco::Muon muon1 = muons_h->at(j);
          if ( muon1.pt() < MIN_MUON_PT ) {
            continue;
          }
          //TODO testing same sign muons
          if ( muon0.charge() ==  muon1.charge() ) {
            continue;
          }
          if (!muon::isSoftMuon(muon1, primary_vertex)) {
            continue;
          }

          reco::TrackRef muon_track0 = GetMuonTrackRef( muon0 );
          reco::TrackRef muon_track1 = GetMuonTrackRef( muon1 );
          transient_tracks.clear();
          transient_tracks.push_back( (*track_builder).build(muon_track0.get()));
          transient_tracks.push_back( (*track_builder).build(muon_track1.get()));
          TransientVertex dimuon_vertex;
          if (transient_tracks.size() > 1) {
            KalmanVertexFitter kalman_fitter;
            dimuon_vertex = kalman_fitter.vertex(transient_tracks);
            if ( !dimuon_vertex.isValid()) {
              continue;
            }
            double vertex_probability = TMath::Prob(dimuon_vertex.totalChiSquared(), int(dimuon_vertex.degreesOfFreedom()));
            if ( vertex_probability < MIN_VERTEX_PROB ) {
              continue;
            }
          }
          //TODO clean up this code
          //std::cout << "distance " << distance << " dist_err " << dist_err << " chi2 " << chi2 << std::endl;
          //std::cout << "dimuon_vertex.position() " << dimuon_vertex.position() << std::endl;
          InitJPsi( muon0, muon1, dimuon_vertex );
          //TODO - figure out how best to handle events with multiple j/psis
          //set_both_mu(muon0, muon1);
          //TODO how should I handle multiple j/psi in same event?
          mu0.push_back( muon0 );
          mu1.push_back( muon1 );
        }
      }
    }
    InitJets(iEvent, iSetup);
  }

        void ZFinderEvent::InitGSFElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
          // We split this part into a new function because it is very long
          // Most of this code is stolen from the example here:
          // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/src/EGammaCutBasedEleIdAnalyzer.cc?view=markup

          // electrons
          edm::Handle<reco::GsfElectronCollection> els_h;
          iEvent.getByLabel(inputtags_.ecal_electron, els_h);

          // conversions
          edm::Handle<reco::ConversionCollection> conversions_h;
          iEvent.getByLabel(inputtags_.conversion, conversions_h);

          // iso deposits
          typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
          IsoDepositVals isoVals(inputtags_.iso_vals.size());
          for (size_t j = 0; j < inputtags_.iso_vals.size(); ++j) {
            iEvent.getByLabel(inputtags_.iso_vals[j], isoVals[j]);
          }

          // beam spot
          edm::Handle<reco::BeamSpot> beamspot_h;
          iEvent.getByLabel(inputtags_.beamspot, beamspot_h);
          const reco::BeamSpot &beamSpot = *(beamspot_h.product());

          // vertices
          edm::Handle<reco::VertexCollection> vtx_h;
          iEvent.getByLabel(inputtags_.vertex, vtx_h);

          // rho for isolation
          // The python uses:
          // cms.InputTag("kt6PFJetsForIsolation", "rho")
          edm::Handle<double> rho_iso_h;
          iEvent.getByLabel(inputtags_.rho_iso, rho_iso_h);
          const double RHO_ISO = *(rho_iso_h.product());

          // loop on electrons
          for(unsigned int i = 0; i < els_h->size(); ++i) {
            // Get the electron and set put it into the electrons vector
            reco::GsfElectron electron = els_h->at(i);
            //TODO for now hardcoding the requirement for medium electrons nvm

            ZFinderElectron* zf_electron = AddRecoElectron(electron);
            //ZFinderElectron* zf_electron = new ZFinderElectron(electron);

            // get reference to electron and the electron
            reco::GsfElectronRef ele_ref(els_h, i);

            // get particle flow isolation
            const double ISO_CH = (*(isoVals)[0])[ele_ref];
            const double ISO_EM = (*(isoVals)[1])[ele_ref];
            const double ISO_NH = (*(isoVals)[2])[ele_ref];

            // test ID
            // working points
            const bool VETO = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO);
            const bool LOOSE = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO);
            const bool MEDIUM = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO);
            const bool TIGHT = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO);

            // eop/fbrem cuts for extra tight ID
            const bool FBREMEOPIN = EgammaCutBasedEleId::PassEoverPCuts(ele_ref);

            // cuts to match tight trigger requirements
            const bool TRIGTIGHT = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, ele_ref);

            // for 2011 WP70 trigger
            const bool TRIGWP70 = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERWP70, ele_ref);

            // Add the cuts to our electron
            const double WEIGHT = 1.;
            zf_electron->AddCutResult("eg_veto", VETO, WEIGHT);
            zf_electron->AddCutResult("eg_loose", LOOSE, WEIGHT);
            zf_electron->AddCutResult("eg_medium", MEDIUM, WEIGHT);
            zf_electron->AddCutResult("eg_tight", TIGHT, WEIGHT);
            zf_electron->AddCutResult("eg_eop_cut", FBREMEOPIN, WEIGHT);
            zf_electron->AddCutResult("eg_trigtight", TRIGTIGHT, WEIGHT);
            zf_electron->AddCutResult("eg_trigwp70", TRIGWP70, WEIGHT);

            // Check for trigger matching
            const bool EE_TIGHT = TriggerMatch(iEvent, ET_ET_TIGHT, zf_electron->eta, zf_electron->phi, TRIG_DR_);
            const bool EE_LOOSE = TriggerMatch(iEvent, ET_ET_LOOSE, zf_electron->eta, zf_electron->phi, TRIG_DR_);
            const bool EE_DZ = TriggerMatch(iEvent, ET_ET_DZ, zf_electron->eta, zf_electron->phi, TRIG_DR_);
            const bool EENT_TIGHT = TriggerMatch(iEvent, ET_NT_ET_TIGHT, zf_electron->eta, zf_electron->phi, TRIG_DR_);
            const bool EEHF_TIGHT = EENT_TIGHT;
            const bool EEHF_LOOSE = TriggerMatch(iEvent, ET_HF_ET_LOOSE, zf_electron->eta, zf_electron->phi, TRIG_DR_);
            const bool SINGLE_E = TriggerMatch(iEvent, SINGLE_ELECTRON_TRIGGER, zf_electron->eta, zf_electron->phi, TRIG_DR_);

            zf_electron->AddCutResult("trig(et_et_tight)", EE_TIGHT, WEIGHT);
            zf_electron->AddCutResult("trig(et_et_loose)", EE_LOOSE, WEIGHT);
            zf_electron->AddCutResult("trig(et_et_dz)", EE_DZ, WEIGHT);
            zf_electron->AddCutResult("trig(et_nt_etleg)", EENT_TIGHT, WEIGHT);
            zf_electron->AddCutResult("trig(et_hf_tight)", EEHF_TIGHT, WEIGHT);
            zf_electron->AddCutResult("trig(et_hf_loose)", EEHF_LOOSE, WEIGHT);
            zf_electron->AddCutResult("trig(single_ele)", SINGLE_E, WEIGHT);

            //TODO decide what if any cuts to put on electrons
            //if ( zf_electron->CutPassed("eg_medium") ) {
            //    //AddRecoElectron(*zf_electron);
            //    //TODO clean up this code/put into a function
            //    if (zf_electron->charge == 1) {
            //        reco_electrons_.push_back(zf_electron);
            //    }
            //    if (zf_electron->charge == -1) {
            //        reco_anti_electrons_.push_back(zf_electron);
            //    }
            //}
          }
        }

        void ZFinderEvent::InitHFElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
          // HF Electrons
          edm::Handle<reco::RecoEcalCandidateCollection> els_h;
          iEvent.getByLabel(inputtags_.hf_electron, els_h);
          // HF Superclusters
          edm::Handle<reco::SuperClusterCollection> scs_h;
          iEvent.getByLabel(inputtags_.hf_clusters, scs_h);
          edm::Handle<reco::HFEMClusterShapeAssociationCollection> scas_h;
          iEvent.getByLabel(inputtags_.hf_clusters, scas_h);

          // Loop over electrons
          for(unsigned int i = 0; i < els_h->size(); ++i) {
            // Get the electron and set put it into the electrons vector
            reco::RecoEcalCandidate electron = els_h->at(i);
            ZFinderElectron* zf_electron = AddRecoElectron(electron);

            reco::SuperClusterRef cluster_ref = electron.superCluster();
            const reco::HFEMClusterShapeRef CLUSTER_SHAPE_REF = scas_h->find(cluster_ref)->val;
            const reco::HFEMClusterShape& CLUSTER_SHAPE = *CLUSTER_SHAPE_REF;

            const double ECE9 = CLUSTER_SHAPE.eCOREe9();
            const double ESEL = CLUSTER_SHAPE.eSeL();
            const double E9E25 = (CLUSTER_SHAPE.eLong3x3() * 1.0 / CLUSTER_SHAPE.eLong5x5());

            // e9e25 cut
            const bool PASS_E9E25 = (E9E25 > 0.94);

            // HF Tight (as defined in hfRecoEcalCandidate_cfi.py in ZShape)
            const double TIGHT2D = (ECE9 - (ESEL * 0.20));
            const bool HFTIGHT = (TIGHT2D > 0.92);

            // HF Medium
            const double MEDIUM2D = (ECE9 - (ESEL * 0.275));
            const bool HFMEDIUM = (MEDIUM2D > 0.875);

            // HF Loose
            const double LOOSE2D = (ECE9 - (ESEL * 0.475));
            const bool HFLOOSE = (LOOSE2D > 0.815);

            // Add the cuts to our electron
            const double WEIGHT = 1.;
            zf_electron->AddCutResult("hf_e9e25", PASS_E9E25, WEIGHT);
            zf_electron->AddCutResult("hf_2dtight", HFTIGHT, WEIGHT);
            zf_electron->AddCutResult("hf_2dmedium", HFMEDIUM, WEIGHT);
            zf_electron->AddCutResult("hf_2dloose", HFLOOSE, WEIGHT);

            // Check for trigger matching
            const bool HIGHLOW_03 = TriggerMatch(iEvent, ET_HF_HF_LOOSE, zf_electron->eta, zf_electron->phi, TRIG_DR_);
            zf_electron->AddCutResult("trig(hf_loose)", HIGHLOW_03, WEIGHT);

            const bool LOWHIGH_03 = TriggerMatch(iEvent, ET_HF_HF_TIGHT, zf_electron->eta, zf_electron->phi, TRIG_DR_);
            zf_electron->AddCutResult("trig(hf_tight)", LOWHIGH_03, WEIGHT);
          }
        }

        void ZFinderEvent::InitNTElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
          // NT Electrons
          edm::Handle<reco::PhotonCollection> els_h;
          iEvent.getByLabel(inputtags_.nt_electron, els_h);

          // Loop over all electrons
          for(unsigned int i = 0; i < els_h->size(); ++i) {
            reco::Photon electron = els_h->at(i);
            // Because the photon collect is NOT filtered for electrons, we
            // reject all electrons outside of the NT region of ECAL.
            if (2.5 < abs(electron.eta()) && abs(electron.eta()) < 2.850) {
              ZFinderElectron* zf_electron = AddRecoElectron(electron);

              // Apply Alexey's Cuts
              //const double PHOTON_ET = electron.superCluster()->rawEnergy() * sin(electron.superCluster()->theta());
              if (       0.89 < electron.r9() && electron.r9() < 1.02
                  && electron.hadronicOverEm() < 0.05
                  && abs(electron.superCluster()->eta()) > 2.5
                  //&& PHOTON_ET > 20.
                  && electron.sigmaIetaIeta() < 0.029
                  && (electron.ecalRecHitSumEtConeDR03() / electron.pt()) < 0.035
                  && (electron.hcalTowerSumEtConeDR03() / electron.pt()) < 0.11
                 ) {
                const bool PASSED = true;
                const double WEIGHT = 1.;
                zf_electron->AddCutResult("nt_loose", PASSED, WEIGHT);
              }

              // Check for trigger matching
              // HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele15_CaloIdT_CaloIsoVL_trackless_v8
            }
          }
        }
        //TODO implement this if so desired
        //    void ZFinderEvent::InitMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
        //        // Muons
        //        edm::Handle<reco::MuonCollection> muons_h;
        //        iEvent.getByLabel(inputtags_.muon, muons_h);
        //
        //        // Loop over all muons
        //        for(unsigned int i = 0; i < muons_h->size(); ++i) {
        //            reco::Muon muon = muons_h->at(i);
        //            AddRecoMuon(muon);
        //            //TODO remove this and below comment
        //            //ZFinderMuon* zf_muon = AddRecoMuon(muon);
        //        }
        //    }

        void ZFinderEvent::InitZ(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
          if (e0 != NULL && e1 != NULL) {
            //TODO make a function to handle vertex getting
            edm::ESHandle<TransientTrackBuilder> track_builder_e;
            iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_e);
            std::vector<reco::TransientTrack> transient_tracks_e;

            //reco::TrackRef electron_track0 = GetElectronTrackRef( e0->gsf_elec_ );
            //reco::TrackRef electron_track1 = GetElectronTrackRef( e1->gsf_elec_ );

            reco::GsfTrackRef electron_track0 = e0->gsf_elec_.gsfTrack() ;
            reco::GsfTrackRef electron_track1 = e1->gsf_elec_.gsfTrack() ;
            transient_tracks_e.clear();
            transient_tracks_e.push_back( (*track_builder_e).build(electron_track0.get()));
            transient_tracks_e.push_back( (*track_builder_e).build(electron_track1.get()));
            TransientVertex dielectron_vertex;
            if (transient_tracks_e.size() > 1) {
              KalmanVertexFitter kalman_fitter_e;
              dielectron_vertex = kalman_fitter_e.vertex(transient_tracks_e);
              if ( dielectron_vertex.isValid()) {
                reco_z.vtx_prob = TMath::Prob(dielectron_vertex.totalChiSquared(), int(dielectron_vertex.degreesOfFreedom()));
              }
            }
            // Set Z properties
            if (reco_z.vtx_prob >= MIN_VERTEX_PROB) {
              const double ELECTRON_MASS = 5.109989e-4;
              math::PtEtaPhiMLorentzVector e0lv(e0->pt, e0->eta, e0->phi, ELECTRON_MASS);
              math::PtEtaPhiMLorentzVector e1lv(e1->pt, e1->eta, e1->phi, ELECTRON_MASS);
              math::PtEtaPhiMLorentzVector zlv;
              zlv = e0lv + e1lv;

              reco_z.vtx_x = dielectron_vertex.position().x();
              reco_z.vtx_y = dielectron_vertex.position().y();
              reco_z.vtx_z = dielectron_vertex.position().z();

              reco_z.m = zlv.mass();
              reco_z.y = zlv.Rapidity();
              reco_z.phi = zlv.phi();
              reco_z.pt = zlv.pt();
              reco_z.phistar = ReturnPhistar(e0->eta, e0->phi, e1->eta, e1->phi);
              reco_z.eta = zlv.eta();
              reco_z.vtx = dielectron_vertex;
            }
          }
        }

        //void ZFinderEvent::InitJPsi(zf::ZFinderMuon* mu0, zf::ZFinderMuon* mu1) {
        //TODO make this take event, event setup as opposed to two muons
        void ZFinderEvent::InitJPsi(const reco::Muon &mu0, const reco::Muon &mu1, const TransientVertex &dimuon_vertex) {
          const double MUON_MASS = 0.1056583715;
          const double C = 29.979245800; // cm/ns
          math::PtEtaPhiMLorentzVector mu0lv(mu0.pt(), mu0.eta(), mu0.phi(), MUON_MASS);
          math::PtEtaPhiMLorentzVector mu1lv(mu1.pt(), mu1.eta(), mu1.phi(), MUON_MASS);
          math::PtEtaPhiMLorentzVector jpsi_lv;
          jpsi_lv = mu0lv + mu1lv;

          double pos_x = dimuon_vertex.position().x();
          double pos_y = dimuon_vertex.position().y();
          double pos_z = dimuon_vertex.position().z();

          //TODO maybe this is not a good place for this variable but whatever

          double x = -10000;
          double y = -10000;
          double z = -10000;
          double LP_XY = -10000; 
          double tau_xy = -10000 ; // ns
          double LP_Z = -10000; 
          double tau_z = -10000; // ns
          double vertex_probability = -10000;
          double distance = -10000;
          double dist_err = -10000;
          double chi2 = -10000;
          double distance_xy = -10000;
          double dist_err_xy = -10000;
          double chi2_xy = -10000;

          double px = jpsi_lv.px();
          double py = jpsi_lv.py();
          double pz = jpsi_lv.pz();
          double pt = jpsi_lv.pt();

          if ( reco_z.vtx_x != -100 ) {
            x = (pos_x - reco_z.vtx.position().x() );
            y = (pos_y - reco_z.vtx.position().y() );
            z = (pos_z - reco_z.vtx.position().z() );
            LP_XY = ((x * px) + (y * py)) ; // 2d
            tau_xy = jpsi_lv.mass() * LP_XY / (pt * pt) / C; // ns

            LP_Z = z * pz ; // 2d
            tau_z = jpsi_lv.mass() * LP_Z / (pz * pz) / C; // ns

            VertexDistance3D vertTool;
            distance = vertTool.distance(reco_z.vtx, dimuon_vertex).value();
            dist_err = vertTool.distance(reco_z.vtx, dimuon_vertex).error();
            chi2 = vertTool.compatibility(reco_z.vtx, dimuon_vertex);

            VertexDistanceXY vertTool_xy;
            distance_xy = vertTool_xy.distance(reco_z.vtx, dimuon_vertex).value();
            dist_err_xy = vertTool_xy.distance(reco_z.vtx, dimuon_vertex).error();
            chi2_xy = vertTool_xy.compatibility(reco_z.vtx, dimuon_vertex);

            vertex_probability = TMath::Prob(dimuon_vertex.totalChiSquared(), int(dimuon_vertex.degreesOfFreedom()));
          }
          else if ( reco_vert.primary_x != -100 )
          {
            x = (pos_x - reco_vert.primary_vert.position().x() );
            y = (pos_y - reco_vert.primary_vert.position().y() );
            z = (pos_z - reco_vert.primary_vert.position().z() );
            LP_XY = ((x * px) + (y * py)) ; // 2d
            tau_xy = jpsi_lv.mass() * LP_XY / (pt * pt) / C; // ns

            LP_Z = z * pz ; // 2d
            tau_z = jpsi_lv.mass() * LP_Z / (pz * pz) / C; // ns

            VertexDistance3D vertTool;
            distance = vertTool.distance(reco_vert.primary_vert, dimuon_vertex).value();
            dist_err = vertTool.distance(reco_vert.primary_vert, dimuon_vertex).error();
            chi2 = vertTool.compatibility(reco_vert.primary_vert, dimuon_vertex);

            VertexDistanceXY vertTool_xy;
            distance_xy = vertTool_xy.distance(reco_vert.primary_vert, dimuon_vertex).value();
            dist_err_xy = vertTool_xy.distance(reco_vert.primary_vert, dimuon_vertex).error();
            chi2_xy = vertTool_xy.compatibility(reco_vert.primary_vert, dimuon_vertex);

            vertex_probability = TMath::Prob(dimuon_vertex.totalChiSquared(), int(dimuon_vertex.degreesOfFreedom()));
          }

          if (reco_z.phi != -1000) {
            reco_jpsi.z_delta_phi.push_back ( fabs( deltaPhi( reco_z.phi , jpsi_lv.phi() ) ) );
          }
          else {
            reco_jpsi.z_delta_phi.push_back ( -1000 );
          }

          reco_jpsi.tau_xy.push_back (tau_xy);
          reco_jpsi.tau_z.push_back (tau_z);

          reco_jpsi.vtx_prob.push_back (vertex_probability);

          reco_jpsi.distance_x.push_back (x);
          reco_jpsi.distance_y.push_back (y);
          reco_jpsi.distance_z.push_back (z);

          reco_jpsi.distance_xy.push_back (distance_xy);
          reco_jpsi.dist_err_xy.push_back (dist_err_xy);
          reco_jpsi.chi2_xy.push_back (chi2_xy);

          reco_jpsi.distance.push_back (distance);
          reco_jpsi.dist_err.push_back (dist_err);
          reco_jpsi.chi2.push_back (chi2);

          reco_jpsi.vtx_x.push_back (pos_x);
          reco_jpsi.vtx_y.push_back (pos_y);
          reco_jpsi.vtx_z.push_back (pos_z);

          reco_jpsi.m.push_back (jpsi_lv.mass());
          reco_jpsi.y.push_back (jpsi_lv.Rapidity());
          reco_jpsi.pt.push_back (jpsi_lv.pt());
          reco_jpsi.phistar.push_back (ReturnPhistar(mu0.eta(), mu0.phi(), mu1.eta(), mu1.phi()));
          reco_jpsi.eta.push_back (jpsi_lv.eta());
          reco_jpsi.phi.push_back (jpsi_lv.phi());

          reco_jpsi.muons_delta_phi.push_back ( fabs( deltaPhi( mu0.phi() , mu1.phi() ) ) );
          reco_jpsi.muons_delta_eta.push_back ( fabs( mu0.eta() -  mu1.eta() ) );
          reco_jpsi.muons_deltaR.push_back ( fabs( deltaR( mu0.p4() , mu1.p4() ) ) );

          

          //TODO use a better vector class, implement it as the early J/Psi paper did
          //std::cout << "pt " << jpsi_lv.pt() << " px " << jpsi_lv.px() << " py " << jpsi_lv.py() << " pz " << jpsi_lv.pz() << std::endl;
          //TVectorD v1(0 , 1 , jpsi_lv.px() / jpsi_lv.pt() , jpsi_lv.py() / jpsi_lv.pt() , "END");
          //v1.Print();
          //std::cout << "v1.x " << v1(0) << std::endl;
          //TVectorD vx( 0 , 1 , (pos_x - reco_z.vtx.position().x() ) , (pos_y - reco_z.vtx.position().y() ) , "END" );
          //vx.Print();
          //TMatrix m1 = dimuon_vertex.vertexState().fullCovariance(); 
          //m1.Print();
          //TODO do this with actual vectors, which I am pretty sick of right now
          // for now follow Atlas paper as it is simpler



          if(mu0.isPFIsolationValid()){
            reco_jpsi.iso_sum_charged_hadron_pt_mu0.push_back ( mu0.pfIsolationR04().sumChargedHadronPt);
            reco_jpsi.iso_sum_charged_particle_pt_mu0.push_back ( mu0.pfIsolationR04().sumChargedParticlePt);
            reco_jpsi.iso_sum_neutral_hadron_et_mu0.push_back ( mu0.pfIsolationR04().sumNeutralHadronEt);
            reco_jpsi.iso_sum_photon_et_mu0.push_back ( mu0.pfIsolationR04().sumPhotonEtHighThreshold);
            reco_jpsi.iso_sum_pileup_pt_mu0.push_back ( mu0.pfIsolationR04().sumPUPt);
            reco_jpsi.iso_mu0.push_back ( (mu0.pfIsolationR04().sumChargedHadronPt + mu0.pfIsolationR04().sumNeutralHadronEt +
                mu0.pfIsolationR04().sumPhotonEtHighThreshold ) / mu0.pt());
          }
          if(mu1.isPFIsolationValid()){
            reco_jpsi.iso_sum_charged_hadron_pt_mu1.push_back ( mu1.pfIsolationR04().sumChargedHadronPt);
            reco_jpsi.iso_sum_charged_particle_pt_mu1.push_back ( mu1.pfIsolationR04().sumChargedParticlePt);
            reco_jpsi.iso_sum_neutral_hadron_et_mu1.push_back ( mu1.pfIsolationR04().sumNeutralHadronEt);
            reco_jpsi.iso_sum_photon_et_mu1.push_back ( mu1.pfIsolationR04().sumPhotonEtHighThreshold);
            reco_jpsi.iso_sum_pileup_pt_mu1.push_back ( mu1.pfIsolationR04().sumPUPt);
            reco_jpsi.iso_mu1.push_back ( (mu1.pfIsolationR04().sumChargedHadronPt + mu1.pfIsolationR04().sumNeutralHadronEt +
                mu1.pfIsolationR04().sumPhotonEtHighThreshold ) / mu1.pt());
          }
        }
        void ZFinderEvent::InitJets(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
          // Jets
          edm::Handle<reco::PFJetCollection> jets_h;
          iEvent.getByLabel(inputtags_.jet, jets_h);
          edm::Handle<reco::JetTagCollection> bTagHandle;
          iEvent.getByLabel("trackCountingHighPurBJetTags", bTagHandle);
          const reco::JetTagCollection & bTags = *(bTagHandle.product());

          // Loop over jets and study b tag info.
          //for (unsigned int i = 0; i != bTags.size(); ++i) {
          //    std::cout<<" Jet "<< i 
          //    <<" has b tag discriminator = "<<bTags[i].second
          //    << " and jet Pt = "<<bTags[i].first->pt() << std::endl;
          //}
          n_reco_jets = 0;
          n_reco_muon_jets = 0;
          for(unsigned int i = 0; i < jets_h->size(); ++i) {
            reco::PFJet jet = jets_h->at(i);
            

            if ( jet.pt() > 20 && fabs(jet.eta()) <= 2.4 ) {
              //TODO check these - took from Z+jets talk loose jet id
              if ( jet.neutralHadronEnergyFraction() < 0.99 &&
                  jet.neutralEmEnergyFraction()     < 0.99 &&
                  jet.chargedHadronEnergyFraction () > 0   &&
                  jet.chargedMultiplicity()          > 0   &&
                  jet.chargedEmEnergyFraction()      < 0.99)
              {
                //remove electron jets for two highest pt electrons
                if (e0 != NULL && e1 != NULL) {
                  if ( (deltaR( jet.eta() , jet.phi(), e0->eta, e0->phi ) < 0.3) ||
                       (deltaR( jet.eta() , jet.phi(), e1->eta, e1->phi ) < 0.3) )
                  {
                    continue;
                  }
                }


                bool is_jet_matched = false;
                for (unsigned int k = 0; k != bTags.size(); ++k) {
                  if (bTags[k].first->pt() > 10) {
                    float deltaR_tagger = deltaR(jet.p4(), bTags[k].first->p4() );
                    if ( deltaR_tagger < 0.05 && !is_jet_matched ) //matching between default btagged jets and pfjets
                    {
                      //TODO clean this up a bit
                      is_jet_matched = true;
                      reco_jets.btag_discriminator.push_back(bTags[k].second);
                    }
                  }
                }
                if (!is_jet_matched)
                {
                  reco_jets.btag_discriminator.push_back(-1000);
                }

                n_reco_jets++;
                //TODO need to have jet?
                //jets.push_back(jet);
                reco_jets.pt.push_back(jet.pt());
                reco_jets.eta.push_back(jet.eta());
                reco_jets.phi.push_back(jet.phi());
                if ( jet.muonMultiplicity() > 0 ) {
                  for ( unsigned int j = 0; j < reco_jpsi.y.size() ; ++j ) {
                    if ( deltaR ( jet.eta() , jet.phi(), reco_jpsi.y.at(j), reco_jpsi.phi.at(j) ) < 0.3 )
                    {
                      reco_muon_jets.pt.push_back(jet.pt());
                      reco_muon_jets.eta.push_back(jet.eta());
                      reco_muon_jets.phi.push_back(jet.phi());
                      n_reco_muon_jets++;

                      bool is_muon_jet_matched = false;
                      for (unsigned int k = 0; k != bTags.size(); ++k) {
                        if (bTags[k].first->pt() > 10) {
                          float deltaR_tagger = deltaR(jet.p4(), bTags[k].first->p4() );
                          if ( deltaR_tagger < 0.05 && !is_muon_jet_matched ) //matching between default btagged jets and pfjets
                          {
                            //TODO clean this up a bit
                            is_muon_jet_matched = true;
                            reco_muon_jets.btag_discriminator.push_back(bTags[k].second);
                          }
                        }
                      }
                      if (!is_muon_jet_matched)
                      {
                        reco_muon_jets.btag_discriminator.push_back(-1000);
                      }
                    }
                  }
                }
              }
            }
          }
        }

        void ZFinderEvent::InitVariables() {
          // Beamspot
          reco_bs.x = -1000;
          reco_bs.y = -1000;
          reco_bs.z = -1000;

          // Vertexes
          // TODO truth vert not used, clean up this code
          reco_vert.num = -1;
          //reco_vert.x = -1000;
          //reco_vert.y = -1000;
          //reco_vert.z = -1000;
          reco_vert.primary_x = -100;
          reco_vert.primary_y = -100;
          reco_vert.primary_z = -100;
          truth_vert.num = -1;
          //truth_vert.x = -1000;
          //truth_vert.y = -1000;
          //truth_vert.z = -1000;

          // Event ID
          id.run_num = 0;
          id.lumi_num = 0;
          id.event_num = 0;

          // Z Data
          reco_z.m = -1;
          reco_z.y = -1000;
          reco_z.phi = -1000;
          reco_z.pt = -1;
          reco_z.phistar = -1;
          reco_z.eta = -1000;
          reco_z.vtx_prob = -1000;
          reco_z.vtx_x = -100;
          reco_z.vtx_y = -100;
          reco_z.vtx_z = -100;
          truth_z.m = -1;
          truth_z.y = -1000;
          truth_z.pt = -1;
          truth_z.phistar = -1;
          truth_z.eta = -1000;

          // JPsi Data TODO erase this if keep vector method
          //reco_jpsi.m = -1;
          //reco_jpsi.y = -1000;
          //reco_jpsi.pt = -1;
          //reco_jpsi.phistar = -1;
          //reco_jpsi.eta = -1000;
          //reco_jpsi.iso_mu0 = -1.0;
          //reco_jpsi.iso_sum_charged_hadron_pt_mu0 = -1.0;
          //reco_jpsi.iso_sum_charged_particle_pt_mu0 = -1.0;
          //reco_jpsi.iso_sum_neutral_hadron_et_mu0 = -1.0;
          //reco_jpsi.iso_sum_photon_et_mu0 = -1.0;
          //reco_jpsi.iso_sum_pileup_pt_mu0 = -1.0;
          //reco_jpsi.iso_mu1 = -1.0;
          //reco_jpsi.iso_sum_charged_hadron_pt_mu1 = -1.0;
          //reco_jpsi.iso_sum_charged_particle_pt_mu1 = -1.0;
          //reco_jpsi.iso_sum_neutral_hadron_et_mu1 = -1.0;
          //reco_jpsi.iso_sum_photon_et_mu1 = -1.0;
          //reco_jpsi.iso_sum_pileup_pt_mu1 = -1.0;

          // Electrons
          e0 = NULL;
          e1 = NULL;
          n_reco_electrons = -1;
          e0_truth = NULL;
          e1_truth = NULL;
          e0_trig = NULL;
          e1_trig = NULL;

          // Muons
          //mu0 = NULL;
          //mu1 = NULL;
          n_reco_muons = -1;

          // Jets
          n_reco_jets = -1;

          // Is Data
          is_real_data = false;
        }

        void ZFinderEvent::InitTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
          /* Count Pile Up */
          edm::Handle<std::vector<PileupSummaryInfo> > pileup_info;
          iEvent.getByLabel(inputtags_.pileup, pileup_info);
          if (pileup_info.isValid()) {
            truth_vert.num = pileup_info->size();
          } else {
            truth_vert.num = -1;
          }

          /*
           * We don't need to select electrons with cuts, because in Monte Carlo we
           * can just ask for the Z.
           */
          edm::Handle<reco::GenParticleCollection> mc_particles;
          iEvent.getByLabel(inputtags_.generator, mc_particles);

          /* Finding the Z and daughter electrons
           *
           * We loop over all gen particles. If it is a Z, we check its daughters
           * until we find an electron, then we know that it is a Z->ee decay. If
           * this is the first Z we save it. If the particle is an electron, we
           * make sure it came from a Z. This might have problems in ZZ->eeee
           * decays, but we expect those to be impossibly rare.
           */
          const reco::GenParticle* electron_0 = NULL;
          const reco::GenParticle* electron_1 = NULL;
          const reco::GenParticle* z_boson = NULL;

          for(unsigned int i = 0; i < mc_particles->size(); ++i) {
            const reco::GenParticle* gen_particle = &mc_particles->at(i);
            // Is a Z
            if (gen_particle->pdgId() == ZBOSON && z_boson == NULL) {
              for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
                if (gen_particle->daughter(j)->pdgId() == ELECTRON) {
                  z_boson = gen_particle;
                  break;
                }
              }
              // Is an electron
            } else if (   fabs(gen_particle->pdgId()) == ELECTRON  // In pdgId, fabs(POSITRON) == ELECTRON
                && (electron_0 == NULL || electron_1 == NULL)
                ) {
              for (size_t j = 0; j < gen_particle->numberOfMothers(); ++j) {
                if (gen_particle->mother(j)->pdgId() == ZBOSON) {
                  if (electron_0 == NULL) {
                    electron_0 = gen_particle;
                  } else {
                    electron_1 = gen_particle;
                  }
                }
              }
            }
          }

          // Continue only if all particles have been found
          if (z_boson != NULL && electron_0 != NULL && electron_1 != NULL) {
            // We set electron_0 to the higher pt electron
            if (electron_0->pt() < electron_1->pt()) {
              std::swap(electron_0, electron_1);
            }

            // Add electrons
            ZFinderElectron* zf_electron_0 = AddTruthElectron(*electron_0);
            set_e0_truth(zf_electron_0);
            ZFinderElectron* zf_electron_1 = AddTruthElectron(*electron_1);
            set_e1_truth(zf_electron_1);

            // Z Properties
            truth_z.m = z_boson->mass();
            truth_z.pt = z_boson->pt();
            const double ZEPP = z_boson->energy() + z_boson->pz();
            const double ZEMP = z_boson->energy() - z_boson->pz();
            truth_z.y = 0.5 * log(ZEPP / ZEMP);
            truth_z.phistar = ReturnPhistar(electron_0->eta(), electron_0->phi(), electron_1->eta(), electron_1->phi());
            truth_z.eta = z_boson->eta();
          }
        }

        void ZFinderEvent::InitTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
          // Get the trigger objects that are closest in dR to our reco electrons
          if (e0 != NULL && e1 != NULL) {
            const trigger::TriggerObject* trig_obj_0 = GetBestMatchedTriggerObject(iEvent, ALL_TRIGGERS, e0->eta, e0->phi);
            const trigger::TriggerObject* trig_obj_1 = GetBestMatchedTriggerObject(iEvent, ALL_TRIGGERS, e1->eta, e1->phi);

            // If the electrons are good, set them as our trigger electrons
            if (trig_obj_0 != NULL) {
              ZFinderElectron* tmp_e0 = AddHLTElectron(*trig_obj_0);
              set_e0_trig(tmp_e0);
            }
            if (trig_obj_1 != NULL) {
              ZFinderElectron* tmp_e1 = AddHLTElectron(*trig_obj_1);
              set_e1_trig(tmp_e1);
            }
          }
        }

        ZFinderElectron* ZFinderEvent::AddRecoElectron(reco::GsfElectron electron) {
          ZFinderElectron* zf_electron = new ZFinderElectron(electron);
          if (zf_electron->charge == 1) {
            reco_electrons_.push_back(zf_electron);
          }
          if (zf_electron->charge == -1) {
            reco_anti_electrons_.push_back(zf_electron);
          }
          return zf_electron;
        }

        void ZFinderEvent::AddRecoElectron(zf::ZFinderElectron zf_electron) {
          //TODO clean/fix this up
          //ZFinderElectron* zf_electron = new ZFinderElectron(electron);
          if (zf_electron.charge == 1) {
            reco_electrons_.push_back(&zf_electron);
          }
          if (zf_electron.charge == -1) {
            reco_anti_electrons_.push_back(&zf_electron);
          }
        }

        ZFinderElectron* ZFinderEvent::AddRecoElectron(reco::RecoEcalCandidate electron) {
          ZFinderElectron* zf_electron = new ZFinderElectron(electron);
          reco_electrons_.push_back(zf_electron);
          return zf_electron;
        }

        ZFinderElectron* ZFinderEvent::AddRecoElectron(reco::Photon electron) {
          ZFinderElectron* zf_electron = new ZFinderElectron(electron);
          reco_electrons_.push_back(zf_electron);
          return zf_electron;
        }

        ZFinderElectron* ZFinderEvent::AddTruthElectron(reco::GenParticle electron) {
          ZFinderElectron* zf_electron = new ZFinderElectron(electron);
          truth_electrons_.push_back(zf_electron);
          return zf_electron;
        }

        ZFinderElectron* ZFinderEvent::AddHLTElectron(trigger::TriggerObject electron) {
          ZFinderElectron* zf_electron = new ZFinderElectron(electron);
          hlt_electrons_.push_back(zf_electron);
          return zf_electron;
        }

        void ZFinderEvent::AddRecoMuon(reco::Muon muon) {
          ZFinderMuon* zf_muon = new ZFinderMuon(muon);
          reco_muons_.push_back(zf_muon);
          //TODO remove this and below comment
          //return zf_muon;
        }

        //TODO fix this up
        //reco::TrackRef ZFinderEvent::GetMuonTrackRef(const reco::Muon & mu) {
        reco::TrackRef ZFinderEvent::GetMuonTrackRef(const reco::Muon & mu) {
          reco::TrackRef track;
          if(mu.isStandAloneMuon()) {
            track = mu.outerTrack();
          }
          if(mu.isTrackerMuon()) {
            track = mu.innerTrack();
          }
          if(mu.isGlobalMuon()) {
            track = mu.globalTrack();
          }
          return track;
        }
        //TODO improve this code

        //reco::TrackRef ZFinderEvent::GetElectronTrackRef(const zf::ZFinderElectron & e) {
        //reco::TrackRef ZFinderEvent::GetElectronTrackRef(const reco::GsfElectron & e) {
        //    reco::TrackRef track;
        //    track = e.gsfTrack();
        //    return track;
        //}

        double ZFinderEvent::ReturnPhistar(const double& eta0, const double& phi0, const double& eta1, const double& phi1) {
          /* Calculate phi star */
          static const double PI = 3.14159265358979323846;
          double dphi = phi0 - phi1;

          // Properly account for the fact that 2pi == 0.
          if (dphi < 0){
            if (dphi > -PI){
              dphi = fabs(dphi);
            }
            if (dphi < -PI) {
              dphi += 2*PI;
            }
          }
          if (dphi > PI){
            dphi = 2*PI - dphi;
          }

          const double DETA = fabs(eta0 - eta1);

          /* PhiStar */
          return ( 1 / cosh( DETA / 2 ) ) * (1 / tan( dphi / 2 ) );
        }

        void ZFinderEvent::PrintCuts(ZFinderElectron* zf_elec) {
          using std::cout;
          using std::endl;
          // Print all the cuts of the given zf_elec
          for (auto& i_cut : *zf_elec->GetAllCuts()) {
            cout << "\t\t" << i_cut->name << ": pass " << i_cut->passed << " weight " << i_cut->weight << endl;
          }
        }

        void ZFinderEvent::PrintElectrons(const int TYPE, const bool PRINT_CUTS) {
          using std::cout;
          using std::endl;

          enum ETYPE {
            RECO = 0,
            TRUTH = 1,
            TRIG = 2
          };
          /*
           * Loops over the electrons, and prints out the information about them.
           */
          cout << "Run " << id.run_num;
          cout << " event " << id.event_num;
          if (TYPE == RECO) {
            cout << " Reco Z Mass " << reco_z.m << std::endl;
            for (auto& i_elec : reco_electrons_) {
              cout << "\tpt: " << i_elec->pt;
              cout << " eta: " << i_elec->eta;
              cout << " phi: " << i_elec->phi << endl;
              if (PRINT_CUTS) { PrintCuts(i_elec); }
            }
          } else if (TYPE == TRUTH && !is_real_data) {
            if (e0_truth != NULL && e1_truth != NULL) {
              cout << " Truth Z Mass " << truth_z.m << endl;
              cout << "\tpt: " << e0_truth->pt;
              cout << " eta: " << e0_truth->eta;
              cout << " phi: " << e0_truth->phi << endl;
              if (PRINT_CUTS) { PrintCuts(e0_truth); }
              cout << "\tpt: " << e1_truth->pt;
              cout << " eta: " << e1_truth->eta;
              cout << " phi: " << e1_truth->phi << endl;
              if (PRINT_CUTS) { PrintCuts(e1_truth); }
            }
          } else if (TYPE == TRIG) {
            if (e0_trig != NULL || e1_trig != NULL) {
              cout << " Trigger Electrons:" << std::endl;
            }
            if (e0_trig != NULL) {
              cout << "\tpt: " << e0_trig->pt;
              cout << " eta: " << e0_trig->eta;
              cout << " phi: " << e0_trig->phi << endl;
              if (PRINT_CUTS) { PrintCuts(e0_trig); }
            }
            if (e1_trig != NULL) {
              cout << "\tpt: " << e1_trig->pt;
              cout << " eta: " << e1_trig->eta;
              cout << " phi: " << e1_trig->phi << endl;
              if (PRINT_CUTS) { PrintCuts(e1_trig); }
            }
          }
        }

        std::vector<ZFinderElectron*>* ZFinderEvent::FilteredElectrons() {
          /*
           * Return all electrons
           */
          std::vector<ZFinderElectron*>* tmp_vec = new std::vector<ZFinderElectron*>();
          for (auto& i_elec : reco_electrons_) {
            tmp_vec->push_back(i_elec);
          }

          return tmp_vec;
        }

        std::vector<ZFinderElectron*>* ZFinderEvent::FilteredElectrons(const std::string& cut_name) {
          /*
           * Return all electrons that pass a specified cut
           */
          std::vector< ZFinderElectron*>* tmp_vec = new std::vector< ZFinderElectron*>();
          for (auto& i_elec : reco_electrons_) {
            if (i_elec->CutPassed(cut_name)) {
              tmp_vec->push_back(i_elec);
            }
          }

          return tmp_vec;
        }

        bool ZFinderEvent::ZDefPassed(const std::string& NAME) const {
          /*
           * Try to find the ZDef name in the map, if it exists return the pass
           * value, else return false.
           */
          std::map<std::string, cutlevel_vector>::const_iterator it = zdef_map_.find(NAME);
          if (it != zdef_map_.end()) {
            const cutlevel_vector* cuts_vec = &it->second;
            bool has_passed = true;
            for (auto& v_it : *cuts_vec) {
              has_passed = v_it.second.pass && has_passed;
            }
            return has_passed;
          } else {
            return false;
          }
        }

        void ZFinderEvent::PrintZDefs(const bool VERBOSE) const {
          /*
           * Loop over all ZDefs and print the results.
           */
          using std::cout;
          using std::endl;
          cout << "ZDefinitions:" << endl;
          for (auto& i_map : zdef_map_) {
            cout << "\t" << i_map.first << ": ";
            cout << ZDefPassed(i_map.first) << endl;
            // If VERBOSE, print out each cutlevel as well
            if (VERBOSE) {
              const cutlevel_vector* clv = &i_map.second;

              for (auto& i_cutlevel : *clv) {
                cout << "\t\t" << i_cutlevel.first << ": " << i_cutlevel.second.pass;
                cout << "t0p1 " << i_cutlevel.second.t0p1_pass << ' ' << i_cutlevel.second.t0p1_eff;
                cout << "t1p0 " << i_cutlevel.second.t1p0_pass << ' ' << i_cutlevel.second.t1p0_eff << endl;
              }
            }
          }
        }

        const cutlevel_vector* ZFinderEvent::GetZDef(const std::string& NAME) const {
          std::map<std::string, cutlevel_vector>::const_iterator it = zdef_map_.find(NAME);
          if (it != zdef_map_.end()) {
            return &(it->second);
          } else {
            return NULL;
          }
        }

        const trig_dr_vec* ZFinderEvent::GetMatchedTriggerObjects(
            const edm::Event& iEvent,
            const std::vector<std::string>& trig_names,
            const double ETA, const double PHI, const double DR_CUT
            ) {
          /*
           * Find all trigger objects that match a vector of trigger names and
           * are within some minimum dR of a specified eta and phi. Return them
           * as a vector of pairs of the object, and the dr.
           */
          // If our vector is empty or the first item is blank
          if (trig_names.size() == 0 || trig_names[0].size() == 0) {
            return NULL;
          }

          // Load Trigger Objects
          edm::InputTag hltTrigInfoTag("hltTriggerSummaryAOD","","HLT");
          edm::Handle<trigger::TriggerEvent> trig_event;

          iEvent.getByLabel(hltTrigInfoTag, trig_event);
          if (!trig_event.isValid() ){
            std::cout << "No valid hltTriggerSummaryAOD." << std::endl;
            return NULL;
          }

          trig_dr_vec* out_v = new trig_dr_vec();
          // Loop over triggers, filter the objects from these triggers, and then try to match
          for (auto& trig_name : trig_names) {
            // Loop over triggers, filter the objects from these triggers, and then try to match
            // Grab objects that pass our filter
            edm::InputTag filter_tag(trig_name, "", "HLT");
            trigger::size_type filter_index = trig_event->filterIndex(filter_tag);
            if(filter_index < trig_event->sizeFilters()) { // Check that the filter is in triggerEvent
              const trigger::Keys& trig_keys = trig_event->filterKeys(filter_index);
              const trigger::TriggerObjectCollection& trig_obj_collection(trig_event->getObjects());
              // Get the objects from the trigger keys
              for (auto& i_key : trig_keys) {
                const trigger::TriggerObject* trig_obj = &trig_obj_collection[i_key];
                const double DR = deltaR(ETA, PHI, trig_obj->eta(), trig_obj->phi());
                // Do Delta R matching, and assign a new object if we have a
                // better match
                if (DR < DR_CUT) {
                  out_v->push_back(std::make_pair(trig_obj, DR));
                }
              }
            }
          }
          return out_v;
        }

        const trigger::TriggerObject* ZFinderEvent::GetBestMatchedTriggerObject(
            const edm::Event& iEvent,
            const std::vector<std::string>& trig_names,
            const double ETA, const double PHI
            ) {
          /* Given the ETA and PHI of a particle, and a list of trigger paths,
           * returns the trigger object from those paths that is closest to the
           * given coordinates. */
          const double MIN_DR = 0.3;
          const trig_dr_vec* trig_vec = GetMatchedTriggerObjects(iEvent, trig_names, ETA, PHI, MIN_DR);

          double best_dr = 1.;
          const trigger::TriggerObject* trig_obj = NULL;
          for (auto& i_obj : *trig_vec) {
            if (i_obj.second < best_dr) {
              best_dr = i_obj.second;
              trig_obj = i_obj.first;
            }
          }
          return trig_obj;
        }

        bool ZFinderEvent::TriggerMatch(
            const edm::Event& iEvent,
            const std::vector<std::string>& trig_names,
            const double ETA, const double PHI, const double DR_CUT
            ) {
          // Get the vector and see if there are objects
          const trig_dr_vec* zev = GetMatchedTriggerObjects(iEvent, trig_names, ETA, PHI, DR_CUT);
          if (zev != NULL && zev->size() >= 1) {
            return true;
          } else {
            return false;
          }
        }
      }  // namespace zf
