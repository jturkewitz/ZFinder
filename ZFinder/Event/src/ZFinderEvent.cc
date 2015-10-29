#include "ZFinder/Event/interface/ZFinderEvent.h"

// Standard Library
#include <algorithm>  // std::sort, std::swap
#include <iostream>  // std::cout, std::endl

// CMSSW
#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"  // reco::PhotonCollection
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"  // reco::SuperClusterCollection, reco::SuperClusterRef
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"  // reco::RecoEcalCandidateCollection
#include "DataFormats/MuonReco/interface/MuonFwd.h" // reco::MuonCollection
#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"  // EgammaCutBasedEleId::PassWP, EgammaCutBasedEleId::*
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" // GenEventInfoProduct
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

#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

//for b-tagging
#include "DataFormats/BTauReco/interface/JetTag.h"

// ZFinder
#include "ZFinder/Event/interface/PDGID.h"  // PDGID enum (ELECTRON, POSITRON, etc.)
#include "ZFinder/Event/interface/TriggerList.h"  // ET_ET_TIGHT, ET_ET_DZ, ET_ET_LOOSE, ET_NT_ET_TIGHT, ET_HF_ET_TIGHT, ET_HF_ET_LOOSE, ET_HF_HF_TIGHT, ET_HF_HF_LOOSE, SINGLE_ELECTRON_TRIGGER, ALL_TRIGGERS
#include "ZFinder/Event/interface/PileupReweighting.h" // RUN_2012_ABCD_TRUE_PILEUP, SUMMER12_53X_MC_TRUE_PILEUP
#include "ZFinder/Event/interface/MuonEfficiency.h"
#include "ZFinder/Event/interface/JpsiEfficiencyTables.h"
#include "ZFinder/Event/interface/JpsiEfficiencyTablesModified.h"

//Math
#include <math.h>
#include <TMath.h>
#include <TLorentzVector.h>

namespace zf {
  /*
   * These variables are hard coded here for easy access, instead of randomly
   * scattering them throughout the code
   */
  // Electrons are considered matched to a trigger object if close than this
  // value
  const double ZFinderEvent::TRIG_DR_ = 0.3;
  /*
   * The edm::LumiReWeighting constructor spams std::cout like mad, and it
   * can't be turned off. We get around this be constructing one static
   * instance and sharing it with all instances of the class.
   */
  edm::LumiReWeighting* ZFinderEvent::lumi_weights_ = NULL;


  ZFinderEvent::ZFinderEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, const edm::ParameterSet& iConfig) {
    /* Given an event, parses them for the information needed to make the
     * class.
     *
     * It selects Z and/or J/Psi from leptonic decay channels.
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
    inputtags_.muon = iConfig.getParameter<edm::InputTag>("muonsInputTag");
    inputtags_.jet = iConfig.getParameter<edm::InputTag>("ak5PFJetsInputTag");
    inputtags_.conversion = iConfig.getParameter<edm::InputTag>("conversionsInputTag");
    inputtags_.beamspot = iConfig.getParameter<edm::InputTag>("beamSpotInputTag");
    inputtags_.rho_iso = iConfig.getParameter<edm::InputTag>("rhoIsoInputTag");
    inputtags_.vertex = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    inputtags_.iso_vals = iConfig.getParameter<std::vector<edm::InputTag> >("isoValInputTags");

    // Truth
    inputtags_.pileup = iConfig.getParameter<edm::InputTag>("pileupInputTag");
    inputtags_.generator = iConfig.getParameter<edm::InputTag>("generatorInputTag");

    //TODO testing
    //edm::InputTag hltInputTag("TriggerResults","","HLT");
    //edm::Handle<edm::TriggerResults> triggerResults;
    //iEvent.getByLabel(hltInputTag, triggerResults);
    //const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);
    //for (int i = 0; i < (int) triggerResults->size(); ++i) { // you might have to do (int) triggerResults->size() - 1??
    //  if (!(triggerResults->accept(i))) {
    //    std::cout << "Trigger " << names.triggerName(i) << " did not pass" << std::endl;
    //  } else {
    //    std::cout << "Trigger " << names.triggerName(i) << " passed" << std::endl;
    //  }
    //}

    // Set up the lumi reweighting, but only if it is MC.
    if (!is_real_data && lumi_weights_ == NULL) {
      const std::string PILEUP_ERA = iConfig.getParameter<std::string>("pileup_era");
      // We use a flag in the python file to set the pileup reweighting
      // to use. If a blank, or an unrecognized, string is passed, then
      // we use the full ABCD reweighting.
      std::vector<float> pileup_distribution_in_data = RUN_2012_ABCD_TRUE_PILEUP;
      std::cout << "Pileup reweighting using era: " << PILEUP_ERA << std::endl;
      if (PILEUP_ERA == "A") {
        pileup_distribution_in_data = RUN_2012_A_TRUE_PILEUP;
      }
      else if (PILEUP_ERA == "B") {
        pileup_distribution_in_data = RUN_2012_B_TRUE_PILEUP;
      }
      else if (PILEUP_ERA == "C") {
        pileup_distribution_in_data = RUN_2012_C_TRUE_PILEUP;
      }
      else if (PILEUP_ERA == "D") {
        pileup_distribution_in_data = RUN_2012_D_TRUE_PILEUP;
      }
      else {
        std::cout << "Using RUN_2012_ABCD_TRUE_PILEUP" << std::endl;
      }
      lumi_weights_ = new edm::LumiReWeighting(
          SUMMER12_53X_MC_TRUE_PILEUP, // MC distribution
          pileup_distribution_in_data // Data distribution
          );
    }
    // Use the lumi reweighting to set the event weight. It is 1. for data,
    // and dependent on the pileup reweighting for MC.
    event_weight = 1.;
    if (!is_real_data) {
      SetMCEventWeight(iEvent);
      if (lumi_weights_ != NULL) {
        SetLumiEventWeight(iEvent);
      }
    }
    // Finish initialization of electrons
    if (!is_real_data) {
      InitTruth(iEvent, iSetup);  // MC
    }
    InitReco(iEvent, iSetup);  // Data
    InitTrigger(iEvent, iSetup);  // Trigger Matching
  }

  void ZFinderEvent::SetLumiEventWeight(const edm::Event& iEvent) {
    /* Reweight the event to correct for pileup (but only MC). This recipe
     * is give on the Twiki:
     * https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupMCReweightingUtilities
     */
    edm::Handle<std::vector<PileupSummaryInfo> > pileup_info;
    iEvent.getByLabel(inputtags_.pileup, pileup_info);

    // Must be a float because weight() below takes float or int
    float true_number_of_pileup = -1.;
    std::vector<PileupSummaryInfo>::const_iterator PILEUP_ELEMENT;
    for(PILEUP_ELEMENT = pileup_info->begin(); PILEUP_ELEMENT != pileup_info->end(); ++PILEUP_ELEMENT) {
      const int BUNCH_CROSSING = PILEUP_ELEMENT->getBunchCrossing();
      if (BUNCH_CROSSING == 0) {
        true_number_of_pileup = PILEUP_ELEMENT->getTrueNumInteractions();
        truth_vert.num = true_number_of_pileup;
        break;
      }
    }
    event_weight *= lumi_weights_->weight(true_number_of_pileup);
  }

  void ZFinderEvent::SetMCEventWeight(const edm::Event& iEvent) {
    // Some MC is also weighted; multiply by this weight also
    edm::Handle<GenEventInfoProduct> gen_event_info;
    iEvent.getByLabel("generator", gen_event_info);
    event_weight *= gen_event_info->weight();
  }

  void ZFinderEvent::InitReco(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    /* Count Pile Up and store first vertex location*/
    edm::Handle<reco::VertexCollection> reco_vertices;
    iEvent.getByLabel(inputtags_.vertex, reco_vertices);
    reco_vert.num = 0;
    bool first_vertex = true;
    reco::Vertex primary_vertex;
    for (unsigned int vertex=0; vertex < reco_vertices->size(); ++vertex) {
      if (    // Criteria copied from twiki (mythical)
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
        // vertex is ordered by highest pt, in descending order
        //std::cout << iEvent.id().event() << "vertex pT: " << sumPtSquared((*reco_vertices)[vertex]) << std::endl;
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

    // Sort our electrons and set e0, e1 as the two with the highest pt
    std::sort(reco_electrons_.begin(), reco_electrons_.end(), SortByPTHighLowElectron);
    std::sort(reco_anti_electrons_.begin(), reco_anti_electrons_.end(), SortByPTHighLowElectron);

    // For Zs
    n_reco_electrons = reco_electrons_.size();
    n_reco_anti_electrons = reco_anti_electrons_.size();

    edm::Handle<reco::MuonCollection> muons_h;
    iEvent.getByLabel(inputtags_.muon, muons_h);
    n_reco_muons = muons_h->size();

    found_four_muons = false;
    if (n_reco_muons >= 4) {
      found_four_muons = true;
    }

    double leading_muon_pt = -1.0;
    double sub_leading_muon_pt = -1.0;

    //Make a Z candidate from the two highest pT oppositely charged muons
    int leading_z_muon_list_position = -1;
    int sub_leading_z_muon_list_position = -1;
    if ( n_reco_muons >= 2 ) {
      for ( int i=0 ; i < n_reco_muons ; ++i ) {
        const reco::Muon muon_temp = muons_h->at(i);
        if ( muon_temp.pt() >= leading_muon_pt ) {
          z_muon0 = muon_temp;
          leading_muon_pt = z_muon0.pt();
          leading_z_muon_list_position = i;
        }
      }
      for ( int j=0 ; j < n_reco_muons ; ++j ) {
        const reco::Muon muon_temp = muons_h->at(j);
        if (muon_temp.pt() >= sub_leading_muon_pt && muon_temp.pt() < leading_muon_pt && (muon_temp.charge() != z_muon0.charge() )) {
          z_muon1 = muon_temp;
          sub_leading_muon_pt = z_muon1.pt();
          sub_leading_z_muon_list_position = j;
        }
      }
      if (leading_muon_pt > -1 && sub_leading_muon_pt > -1) {
        InitZFromMuons(iEvent, iSetup);
      }
    }

    if (n_reco_electrons >= 1 && n_reco_anti_electrons >=1 ) {
      //TODO coding parity between how electrons and how muons are selected
      // Set our internal electrons
      //highest pT electron, and highest pT antielectron
      if (reco_electrons_[0]->pt >= reco_anti_electrons_[0]->pt) {
        set_both_e(reco_electrons_[0], reco_anti_electrons_[0]);
      }
      else {
        set_both_e(reco_anti_electrons_[0], reco_electrons_[0]);
      }
      // Set up the Z To Electrons
      InitZFromElectrons(iEvent, iSetup);
    }


    //Z to muons cut levels
    found_high_pt_muons_from_z = false;
    found_good_muons_from_z = false;
    found_dimuon_z_compatible_vertex = false;
    found_z_to_muons_mass = false;

    //TODO z rapidity
    if (reco_z_from_muons.m > -1 &&
        z_muon0.pt() >= MIN_Z_MUON_PT && z_muon1.pt() >= MIN_Z_MUON_PT ) {
      found_high_pt_muons_from_z = true;
    }
    if (reco_z_from_muons.m > -1 && 
        muon::isTightMuon(z_muon0, reco_vert.primary_vert ) &&
        muon::isTightMuon(z_muon1, reco_vert.primary_vert ) ) {
        //muon::isTightMuon(z_muon1, reco_vert.primary_vert ) &&
        //TriggerMatch(iEvent, DOUBLE_MUON_TIGHT_LEG_TRIGGER, z_muon0.eta(), z_muon0.phi(), TRIG_DR_) &&
        //TriggerMatch(iEvent, DOUBLE_MUON_LOOSE_LEG_TRIGGER, z_muon1.eta(), z_muon1.phi(), TRIG_DR_) ) {
      found_good_muons_from_z = true;
    }
    if (reco_z_from_muons.m > -1 && 
        reco_z_from_muons.vtx_prob >= MIN_VERTEX_PROB ) {
      found_dimuon_z_compatible_vertex = true;
    }
    if (reco_z_from_muons.m >= MIN_Z_MASS && reco_z_from_muons.m <= MAX_Z_MASS) {
      found_z_to_muons_mass = true;
    }

    //Z to electrons cut levels
    found_high_pt_electrons_from_z = false;
    found_good_electrons_from_z = false;
    found_dielectron_z_compatible_vertex = false;
    found_z_to_electrons_mass = false;

    if (reco_z.m > -1 && e0 != NULL && e1 != NULL) {
      if ( e0->pt >= MIN_ELECTRON_PT && e1->pt >= MIN_ELECTRON_PT ) {
        found_high_pt_electrons_from_z = true;
      }
    }
    if (reco_z.m > -1 && e0 != NULL && e1 != NULL) {
      //somewhat arbitrary decision between eg_medium vs eg_tight, eg_medium should have a greater acceptance
      if (e0->CutPassed("eg_medium") && e1->CutPassed("eg_medium") ) {
        found_good_electrons_from_z = true;
      }
    }
    if (reco_z.m > -1 && e0 != NULL && e1 != NULL) {
      if (reco_z.vtx_prob >= MIN_VERTEX_PROB ) {
        found_dielectron_z_compatible_vertex = true;
      }
    }
    if (reco_z.m > -1 && e0 != NULL && e1 != NULL) {
      if ( reco_z.m >= MIN_Z_MASS && reco_z.m <= MAX_Z_MASS ) {
        found_z_to_electrons_mass = true;
      }
    }

    // For JPsis

    edm::ESHandle<TransientTrackBuilder> track_builder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder);
    std::vector<reco::TransientTrack> transient_tracks;

    if (n_reco_muons >= 2 ) {
      // Set up the JPsi
      for ( int i=0 ; i < (n_reco_muons - 1) ; ++i ) {
        const reco::Muon muon0 = muons_h->at(i);
        //Ensure that muons are not shared between JPsi and Z candidates
        if (found_z_to_muons_mass && (i == leading_z_muon_list_position || i == sub_leading_z_muon_list_position)) {
          continue;
        }
        for ( int j=i+1 ; j < n_reco_muons ; ++j ) {
          const reco::Muon muon1 = muons_h->at(j);
          //Ensure that muons are not shared between JPsi and Z candidates
          if (found_z_to_muons_mass && (j == leading_z_muon_list_position || j == sub_leading_z_muon_list_position)) {
            continue;
          }
          //Ensure oppositely charged muons for jpsi candidates
          if ( muon0.charge() ==  muon1.charge() ) {
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
          }
          //Somewhat redundant, but should extra ensure that muon0 is the higher pT muon
          if (muon0.pt() >= muon1.pt() ) {
            InitJPsi( muon0, muon1, dimuon_vertex );
          }
          else {
            InitJPsi( muon1, muon0, dimuon_vertex );
          }
        }
      }
    }
    for (unsigned int i = 0 ; i < reco_jpsi.muon0.size() ; ++i) {
      reco_jpsi.muon0_deltaR_to_truth_muons.push_back( JpsiMuonTruthMatch(reco_jpsi.muon0.at(i)) ) ;
    }
    for (unsigned int i = 0 ; i < reco_jpsi.muon1.size() ; ++i) {
      reco_jpsi.muon1_deltaR_to_truth_muons.push_back( JpsiMuonTruthMatch(reco_jpsi.muon1.at(i)) ) ;
    }
    InitJets(iEvent, iSetup);
 
    //Set cut level flags for jpsi candidates
    //Note that to pass multiple cut stages the same jpsi candidate should pass all stages
    
    found_dimuon_jpsi_with_muons_in_eta_window = false;
    found_dimuon_jpsi_with_high_pt_muons = false;
    found_dimuon_jpsi_with_soft_id_and_high_pt_muons = false;
    found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex = false;
    found_good_dimuon_jpsi_compatible_with_primary_vertex = false;
    found_jpsi = false;
    found_prompt_jpsi = false;
    for (unsigned int i = 0; i < reco_jpsi.m.size() ; ++i ) {
      if (reco_jpsi.has_muons_in_eta_window.at(i) && reco_jpsi.is_in_rap_window.at(i) ) {
        found_dimuon_jpsi_with_muons_in_eta_window = true;
        //TODO should each cut level have its own flag??
        if ( reco_jpsi.has_high_pt_muons.at(i) && reco_jpsi.is_high_pt.at(i) ) {
          //TODO rename this flag or use an additional flag for high pT j/psi
          found_dimuon_jpsi_with_high_pt_muons = true;
          if (reco_jpsi.has_soft_id_muons.at(i) ) {
            found_dimuon_jpsi_with_soft_id_and_high_pt_muons = true;
            if (reco_jpsi.has_muons_with_compatible_vertex.at(i) ) {
              found_dimuon_jpsi_with_good_muons_and_compatible_muon_vertex = true;
              if (reco_jpsi.has_dimuon_vertex_compatible_with_primary_vertex.at(i)) {
                found_good_dimuon_jpsi_compatible_with_primary_vertex = true;
                if (reco_jpsi.is_within_jpsi_mass_window.at(i) ) {
                  found_jpsi = true;
                  if (reco_jpsi.is_prompt.at(i) ) {
                    found_prompt_jpsi = true;
                  }
                }
              }
            }
          }
        }
      }
    }
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

      //ZFinderElectron* zf_electron = AddRecoElectron(electron);
      ZFinderElectron* zf_electron = new ZFinderElectron(electron);

      // get reference to electron and the electron
      reco::GsfElectronRef ele_ref(els_h, i);

      // get particle flow isolation
      const double ISO_CH = (*(isoVals[0]))[ele_ref];
      const double ISO_EM = (*(isoVals[1]))[ele_ref];
      const double ISO_NH = (*(isoVals[2]))[ele_ref];

      //const double EA_TARGET = GetElectronEffectiveArea(ElectronEffectiveAreaType type, Double_t SCEta, ElectronEffectiveAreaTarget EffectiveAreaTarget = kEleEAData2011)

      // test ID
      // working points
      const bool VETO = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO,ElectronEffectiveArea::kEleEAData2012);
      const bool LOOSE = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO,ElectronEffectiveArea::kEleEAData2012);
      const bool MEDIUM = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO,ElectronEffectiveArea::kEleEAData2012);
      const bool TIGHT = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, ele_ref, conversions_h, beamSpot, vtx_h, ISO_CH, ISO_EM, ISO_NH, RHO_ISO,ElectronEffectiveArea::kEleEAData2012);

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

      //et defined to mean total ecal here
      zf_electron->AddCutResult("trig(et_et_tight)", EE_TIGHT, WEIGHT);
      zf_electron->AddCutResult("trig(et_et_loose)", EE_LOOSE, WEIGHT);
      zf_electron->AddCutResult("trig(et_et_dz)", EE_DZ, WEIGHT);
      zf_electron->AddCutResult("trig(et_nt_etleg)", EENT_TIGHT, WEIGHT);
      zf_electron->AddCutResult("trig(et_hf_tight)", EEHF_TIGHT, WEIGHT);
      zf_electron->AddCutResult("trig(et_hf_loose)", EEHF_LOOSE, WEIGHT);
      zf_electron->AddCutResult("trig(single_ele)", SINGLE_E, WEIGHT);

      if (zf_electron->charge == 1) {
        reco_electrons_.push_back(zf_electron);
      }
      if (zf_electron->charge == -1) {
        reco_anti_electrons_.push_back(zf_electron);
      }
    }
  }

  void ZFinderEvent::InitZFromElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (e0 != NULL && e1 != NULL) {
      edm::ESHandle<TransientTrackBuilder> track_builder_e;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_e);
      std::vector<reco::TransientTrack> transient_tracks_e;

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
      const double ELECTRON_MASS = 5.109989e-4;
      math::PtEtaPhiMLorentzVector e0lv(e0->pt, e0->eta, e0->phi, ELECTRON_MASS);
      math::PtEtaPhiMLorentzVector e1lv(e1->pt, e1->eta, e1->phi, ELECTRON_MASS);
      math::PtEtaPhiMLorentzVector zlv;
      zlv = e0lv + e1lv;
      reco_z.zlv = zlv;
      if ( dielectron_vertex.isValid()) {
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
  //rename reco_z to reco_z_from_electrons or something to avoid conflicts
  void ZFinderEvent::InitZFromMuons (const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::ESHandle<TransientTrackBuilder> track_builder_muon;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_muon);
    std::vector<reco::TransientTrack> transient_tracks_muon;
    reco::TrackRef muon_track0 = GetMuonTrackRef( z_muon0 );
    reco::TrackRef muon_track1 = GetMuonTrackRef( z_muon1 );
    transient_tracks_muon.clear();
    transient_tracks_muon.push_back( (*track_builder_muon).build(muon_track0.get()));
    transient_tracks_muon.push_back( (*track_builder_muon).build(muon_track1.get()));
    TransientVertex dimuon_vertex;
    if (transient_tracks_muon.size() > 1) {
      KalmanVertexFitter kalman_fitter;
      dimuon_vertex = kalman_fitter.vertex(transient_tracks_muon);
      if ( dimuon_vertex.isValid()) {
        reco_z_from_muons.vtx_prob = TMath::Prob(dimuon_vertex.totalChiSquared(), int(dimuon_vertex.degreesOfFreedom()));
      }
    }
    const double MUON_MASS = 0.1056583715;
    math::PtEtaPhiMLorentzVector muon0lv(z_muon0.pt(), z_muon0.eta(), z_muon0.phi(), MUON_MASS);
    math::PtEtaPhiMLorentzVector muon1lv(z_muon1.pt(), z_muon1.eta(), z_muon1.phi(), MUON_MASS);
    math::PtEtaPhiMLorentzVector zlv;
    zlv = muon0lv + muon1lv;
    reco_z_from_muons.zlv = zlv;

    if ( dimuon_vertex.isValid()) {

      reco_z_from_muons.vtx_x = dimuon_vertex.position().x();
      reco_z_from_muons.vtx_y = dimuon_vertex.position().y();
      reco_z_from_muons.vtx_z = dimuon_vertex.position().z();

      reco_z_from_muons.m = zlv.mass();
      reco_z_from_muons.y = zlv.Rapidity();
      reco_z_from_muons.phi = zlv.phi();
      reco_z_from_muons.pt = zlv.pt();
      reco_z_from_muons.phistar = ReturnPhistar(z_muon0.eta(), z_muon0.phi(), z_muon1.eta(), z_muon1.phi());
      reco_z_from_muons.eta = zlv.eta();
      reco_z_from_muons.vtx = dimuon_vertex;
    }
  }

  void ZFinderEvent::InitJPsi(const reco::Muon &mu0, const reco::Muon &mu1, const TransientVertex &dimuon_vertex) {
    const double MUON_MASS = 0.1056583715;
    const double C = 29.979245800; // cm/ns
    math::PtEtaPhiMLorentzVector mu0lv(mu0.pt(), mu0.eta(), mu0.phi(), MUON_MASS);
    math::PtEtaPhiMLorentzVector mu1lv(mu1.pt(), mu1.eta(), mu1.phi(), MUON_MASS);
    math::PtEtaPhiMLorentzVector jpsi_lv;
    math::PtEtaPhiMLorentzVector four_lepton_lv;
    jpsi_lv = mu0lv + mu1lv;

    if(found_z_to_electrons_mass) {
      four_lepton_lv = jpsi_lv + reco_z.zlv;
    }
    else if(found_z_to_muons_mass) {
      four_lepton_lv = jpsi_lv + reco_z_from_muons.zlv;
    }
    else {
      four_lepton_lv = jpsi_lv;
    }

    //math::PtEtaPhiMLorentzVector mu0lv_boosted(mu0.pt(), mu0.eta(), mu0.phi(), MUON_MASS);
    //mu0lv_boosted.boost(jpsi_lv.px() / jpsi_lv.energy(), jpsi_lv.py() / jpsi_lv.energy(), jpsi_lv.pz() / jpsi_lv.energy());
    //math::BetaVector beta_vector;
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector beta_vector;
    beta_vector = jpsi_lv.BoostToCM();
    //std::cout << "px " << jpsi_lv.px() << " bv_x  " << beta_vector.x() << std::endl;
    //std::cout << "py " << jpsi_lv.py() << " bv_y  " << beta_vector.y() << std::endl;
    //std::cout << "pz " << jpsi_lv.pz() << " bv_z  " << beta_vector.z() << std::endl;

    TLorentzVector *mu0_lv2 = new TLorentzVector(mu0lv.px(), mu0lv.py(), mu0lv.pz(), mu0lv.energy());
    mu0_lv2->Boost(beta_vector.x(), beta_vector.y(), beta_vector.z() );

    TLorentzVector *mu1_lv2 = new TLorentzVector(mu1lv.px(), mu1lv.py(), mu1lv.pz(), mu1lv.energy());
    mu1_lv2->Boost(beta_vector.x(), beta_vector.y(), beta_vector.z() );
    //std::cout << "jpsi px " << jpsi_lv.px() << " mu0 px " << mu0lv.px() << " boosted px " << mu0_lv2->Px() << std::endl;

    double mu0_px_boosted =  mu0_lv2->Px();
    double mu0_py_boosted =  mu0_lv2->Py();
    double mu0_pz_boosted =  mu0_lv2->Pz();

    double mu1_px_boosted =  mu1_lv2->Px();
    double mu1_py_boosted =  mu1_lv2->Py();
    double mu1_pz_boosted =  mu1_lv2->Pz();

    double pos_x = -99999;
    double pos_y = -99999;
    double pos_z = -99999;
    if ( dimuon_vertex.isValid()) {
      pos_x = dimuon_vertex.position().x();
      pos_y = dimuon_vertex.position().y();
      pos_z = dimuon_vertex.position().z();
    }

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

    //double mu0_px = mu0.px();
    //double mu0_py = mu0.py();
    //double mu0_pz = mu0.pz();

    //double mu1_px = mu1.px();
    //double mu1_py = mu1.py();
    //double mu1_pz = mu1.pz();

 
    if ( found_z_to_electrons_mass ) {
      x = (pos_x - reco_z.vtx.position().x() );
      y = (pos_y - reco_z.vtx.position().y() );
      z = (pos_z - reco_z.vtx.position().z() );
      LP_XY = ((x * px) + (y * py)) ; // 2d
      tau_xy = jpsi_lv.mass() * LP_XY / (pt * pt) / C; // ns

      LP_Z = z * pz ; // 2d
      tau_z = jpsi_lv.mass() * LP_Z / (pz * pz) / C; // ns

      if ( dimuon_vertex.isValid()) {
        VertexDistance3D vertTool;
        distance = vertTool.distance(reco_z.vtx, dimuon_vertex).value();
        dist_err = vertTool.distance(reco_z.vtx, dimuon_vertex).error();
        //chi2 = vertTool.compatibility(reco_z.vtx, dimuon_vertex);

        VertexDistanceXY vertTool_xy;
        distance_xy = vertTool_xy.distance(reco_z.vtx, dimuon_vertex).value();
        dist_err_xy = vertTool_xy.distance(reco_z.vtx, dimuon_vertex).error();
        chi2_xy = vertTool_xy.compatibility(reco_z.vtx, dimuon_vertex);
      }
    }
    else if ( found_z_to_muons_mass ) {
      x = (pos_x - reco_z_from_muons.vtx.position().x() );
      y = (pos_y - reco_z_from_muons.vtx.position().y() );
      z = (pos_z - reco_z_from_muons.vtx.position().z() );
      LP_XY = ((x * px) + (y * py)) ; // 2d
      tau_xy = jpsi_lv.mass() * LP_XY / (pt * pt) / C; // ns

      LP_Z = z * pz ; // 2d
      tau_z = jpsi_lv.mass() * LP_Z / (pz * pz) / C; // ns

      if ( dimuon_vertex.isValid()) {
        VertexDistance3D vertTool;
        distance = vertTool.distance(reco_z_from_muons.vtx, dimuon_vertex).value();
        dist_err = vertTool.distance(reco_z_from_muons.vtx, dimuon_vertex).error();
        //chi2 = vertTool.compatibility(reco_z_from_muons.vtx, dimuon_vertex);

        VertexDistanceXY vertTool_xy;
        distance_xy = vertTool_xy.distance(reco_z_from_muons.vtx, dimuon_vertex).value();
        dist_err_xy = vertTool_xy.distance(reco_z_from_muons.vtx, dimuon_vertex).error();
        chi2_xy = vertTool_xy.compatibility(reco_z_from_muons.vtx, dimuon_vertex);
      }
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

      if ( dimuon_vertex.isValid()) {
        VertexDistance3D vertTool;
        distance = vertTool.distance(reco_vert.primary_vert, dimuon_vertex).value();
        dist_err = vertTool.distance(reco_vert.primary_vert, dimuon_vertex).error();
        //chi2 = vertTool.compatibility(reco_vert.primary_vert, dimuon_vertex);

        VertexDistanceXY vertTool_xy;
        distance_xy = vertTool_xy.distance(reco_vert.primary_vert, dimuon_vertex).value();
        dist_err_xy = vertTool_xy.distance(reco_vert.primary_vert, dimuon_vertex).error();
        chi2_xy = vertTool_xy.compatibility(reco_vert.primary_vert, dimuon_vertex);
      }
    }

    if (reco_z.m >= MIN_Z_MASS && reco_z.m <= MAX_Z_MASS) {
      reco_jpsi.z_delta_phi.push_back ( fabs( deltaPhi( reco_z.phi , jpsi_lv.phi() ) ) );
    }
    else if (reco_z_from_muons.m >= MIN_Z_MASS && reco_z_from_muons.m <= MAX_Z_MASS ) {
      reco_jpsi.z_delta_phi.push_back ( fabs( deltaPhi( reco_z_from_muons.phi , jpsi_lv.phi() ) ) );
    }
    else {
      reco_jpsi.z_delta_phi.push_back ( -1000 );
    }

    if ( dimuon_vertex.isValid()) {
      vertex_probability = TMath::Prob(dimuon_vertex.totalChiSquared(), int(dimuon_vertex.degreesOfFreedom()));
    }

    //TODO testing
    //double mu0_eff = GetEfficiency (SOFT_MUON_DATA_EFF_TABLE, mu0.eta(), mu0.pt() ) ;
    //double mu1_eff = GetEfficiency (SOFT_MUON_DATA_EFF_TABLE, mu1.eta(), mu1.pt() ) ;
    //double mu0_scale_factor = GetEfficiency (SOFT_MUON_SCALE_FACTOR_TABLE, mu0.eta(), mu0.pt());
    //double mu1_scale_factor = GetEfficiency (SOFT_MUON_SCALE_FACTOR_TABLE, mu1.eta(), mu1.pt());
    //double jpsi_acc_eff = GetAccEff (SOFT_MUON_DATA_ACC_EFF_TABLE, jpsi_lv.Rapidity(), jpsi_lv.pt() );

    double mu0_eff = GetEfficiency (SOFT_MUON_DATA_EFF_TABLE_MODIFIED, mu0.eta(), mu0.pt() ) ;
    double mu1_eff = GetEfficiency (SOFT_MUON_DATA_EFF_TABLE_MODIFIED, mu1.eta(), mu1.pt() ) ;
    double mu0_scale_factor = GetEfficiency (SOFT_MUON_SCALE_FACTOR_TABLE_MODIFIED, mu0.eta(), mu0.pt());
    double mu1_scale_factor = GetEfficiency (SOFT_MUON_SCALE_FACTOR_TABLE_MODIFIED, mu1.eta(), mu1.pt());
    double jpsi_acc_eff = GetAccEff (SOFT_MUON_DATA_ACC_EFF_TABLE_MODIFIED, jpsi_lv.Rapidity(), jpsi_lv.pt() );


    //TODO do this with a function, instead of by hand, cos(theta) = dot_product / (a.len() * b.len() )

    double dot_product_mu0 = px * mu0_px_boosted + py * mu0_py_boosted + pz * mu0_pz_boosted;
    double dot_product_mu1 = px * mu1_px_boosted + py * mu1_py_boosted + pz * mu1_pz_boosted;

    double jpsi_p_mag = pow((px * px + py * py + pz * pz), 0.5);
    double mu0_p_mag = pow((mu0_px_boosted * mu0_px_boosted + mu0_py_boosted * mu0_py_boosted + mu0_pz_boosted * mu0_pz_boosted), 0.5);
    double mu1_p_mag = pow((mu1_px_boosted * mu1_px_boosted + mu1_py_boosted * mu1_py_boosted + mu1_pz_boosted * mu1_pz_boosted), 0.5);
    
    double cos_jpsi_mu_plus = -1000;
    double cos_jpsi_mu_minus = -1000;
    if (mu0.charge() == 1 ) {
      cos_jpsi_mu_plus = dot_product_mu0 / (jpsi_p_mag * mu0_p_mag);
      cos_jpsi_mu_minus = dot_product_mu1 / (jpsi_p_mag * mu1_p_mag);
    }
    else {
      cos_jpsi_mu_plus = dot_product_mu1 / (jpsi_p_mag * mu1_p_mag);
      cos_jpsi_mu_minus = dot_product_mu0 / (jpsi_p_mag * mu0_p_mag);
    }

    reco_jpsi.cos_jpsi_mu_plus.push_back(cos_jpsi_mu_plus);
    reco_jpsi.cos_jpsi_mu_minus.push_back(cos_jpsi_mu_minus);
    
    if (mu0.pt() >= mu1.pt() ) {
      reco_jpsi.muon0.push_back (mu0);
      reco_jpsi.muon1.push_back (mu1);
      reco_jpsi.muon0_efficiency.push_back ( mu0_eff ); 
      reco_jpsi.muon1_efficiency.push_back ( mu1_eff ); 
      reco_jpsi.muon0_scale_factor.push_back ( mu0_scale_factor ); 
      reco_jpsi.muon1_scale_factor.push_back ( mu1_scale_factor ); 
    }
    else {
      reco_jpsi.muon0.push_back (mu1);
      reco_jpsi.muon1.push_back (mu0);
      reco_jpsi.muon0_efficiency.push_back ( mu1_eff ); 
      reco_jpsi.muon1_efficiency.push_back ( mu0_eff ); 
      reco_jpsi.muon0_scale_factor.push_back ( mu1_scale_factor ); 
      reco_jpsi.muon1_scale_factor.push_back ( mu0_scale_factor ); 
    }

    reco_jpsi.four_lepton_mass.push_back(four_lepton_lv.mass());

    reco_jpsi.jpsi_acc_eff.push_back(jpsi_acc_eff);

    reco_jpsi.muon0_deltaR_to_z_muons.push_back(JpsiMuonZMuonMatch(mu0));
    reco_jpsi.muon1_deltaR_to_z_muons.push_back(JpsiMuonZMuonMatch(mu1));

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
    reco_jpsi.jpsi_efficiency.push_back ( mu0_eff * mu1_eff);
    reco_jpsi.jpsi_acc_eff.push_back ( jpsi_acc_eff);
    reco_jpsi.jpsi_scale_factor.push_back ( mu0_scale_factor * mu1_scale_factor);

    reco_jpsi.muons_delta_phi.push_back ( fabs( deltaPhi( mu0.phi() , mu1.phi() ) ) );
    reco_jpsi.muons_delta_eta.push_back ( fabs( mu0.eta() -  mu1.eta() ) );
    reco_jpsi.muons_deltaR.push_back ( fabs( deltaR( mu0.p4() , mu1.p4() ) ) );

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



    //Cut results
    
    if ( jpsi_lv.pt() >= MIN_JPSI_PT ) {
      reco_jpsi.is_high_pt.push_back(true);
    }
    else {
      reco_jpsi.is_high_pt.push_back(false);
    }

    //TODO testing
    if (mu0.pt() > mu1.pt() ) {
      if ( mu0.pt() >= MIN_JPSI_LEADING_MUON_PT && mu1.pt() >= MIN_JPSI_SUBLEADING_MUON_PT ) {
        reco_jpsi.has_high_pt_muons.push_back(true);
      }
      else if ( mu0.pt() >= MIN_JPSI_LEADING_MUON_PT && mu1.pt() >= MIN_JPSI_SUBLEADING_MUON_PT_HIGH_ETA && fabs(mu1.eta()) >= 1.2) {
        reco_jpsi.has_high_pt_muons.push_back(true);
      }
      else {
        reco_jpsi.has_high_pt_muons.push_back(false);
      }
    }
    else {
      if ( mu1.pt() >= MIN_JPSI_LEADING_MUON_PT && mu0.pt() >= MIN_JPSI_SUBLEADING_MUON_PT ) {
        reco_jpsi.has_high_pt_muons.push_back(true);
      }
      else if ( mu1.pt() >= MIN_JPSI_LEADING_MUON_PT && mu0.pt() >= MIN_JPSI_SUBLEADING_MUON_PT_HIGH_ETA && fabs(mu0.eta()) >= 1.2) {
        reco_jpsi.has_high_pt_muons.push_back(true);
      }
      else {
        reco_jpsi.has_high_pt_muons.push_back(false);
      }
    }
    
    if (fabs(mu0.eta()) <= MAX_JPSI_MUON_ETA  && fabs(mu1.eta()) <= MAX_JPSI_MUON_ETA ) {
      reco_jpsi.has_muons_in_eta_window.push_back(true);
    }
    else {
      reco_jpsi.has_muons_in_eta_window.push_back(false);
    }

    if (fabs(jpsi_lv.Rapidity()) < MAX_JPSI_RAP) {
      reco_jpsi.is_in_rap_window.push_back(true);
    }
    else {
      reco_jpsi.is_in_rap_window.push_back(false);
    }

    if ( muon::isSoftMuon(mu0, reco_vert.primary_vert )
        && muon::isSoftMuon(mu1, reco_vert.primary_vert) ) {
      reco_jpsi.has_soft_id_muons.push_back(true);
    }
    else {
      reco_jpsi.has_soft_id_muons.push_back(false);
    }

    if ( vertex_probability >= MIN_VERTEX_PROB ) {
      reco_jpsi.has_muons_with_compatible_vertex.push_back(true);
    }
    else {
      reco_jpsi.has_muons_with_compatible_vertex.push_back(false);
    }

    //TODO primary vertex or vertex from z?
    //if ( fabs(reco_vert.primary_z - pos_z ) <= MAX_JPSI_VERTEX_Z_DISPLACEMENT ) {
    //from Z boson if it exists, otherwise primary vertex
    if ( fabs( z ) <= MAX_JPSI_VERTEX_Z_DISPLACEMENT ) {
      reco_jpsi.has_dimuon_vertex_compatible_with_primary_vertex.push_back(true);
    }
    else {
      reco_jpsi.has_dimuon_vertex_compatible_with_primary_vertex.push_back(false);
    }

    if ( jpsi_lv.mass() <= MAX_JPSI_MASS && jpsi_lv.mass() >= MIN_JPSI_MASS ) {
      reco_jpsi.is_within_jpsi_mass_window.push_back(true);
    }
    else {
      reco_jpsi.is_within_jpsi_mass_window.push_back(false);
    }

    if (tau_xy >= MIN_PROMPT_JPSI_TIME && tau_xy <= MAX_PROMPT_JPSI_TIME) {
      reco_jpsi.is_prompt.push_back(true);
    }
    else {
      reco_jpsi.is_prompt.push_back(false);
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
        //taken from a talk + paper on loose jet Id
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
    // TODO use truth vert?
    reco_vert.num = -1;
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
    math::PtEtaPhiMLorentzVector lv(0,0,0,0);
    reco_z.zlv = lv;

    reco_z_from_muons.m = -1;
    reco_z_from_muons.y = -1000;
    reco_z_from_muons.phi = -1000;
    reco_z_from_muons.pt = -1;
    reco_z_from_muons.phistar = -1;
    reco_z_from_muons.eta = -1000;
    reco_z_from_muons.vtx_prob = -1000;
    reco_z_from_muons.vtx_x = -100;
    reco_z_from_muons.vtx_y = -100;
    reco_z_from_muons.vtx_z = -100;
    reco_z_from_muons.zlv = lv;

    truth_z_electrons.m = -1;
    truth_z_electrons.y = -1000;
    truth_z_electrons.pt = -1;
    truth_z_electrons.phistar = -1;
    truth_z_electrons.eta = -1000;

    truth_z_muons.m = -1;
    truth_z_muons.y = -1000;
    truth_z_muons.pt = -1;
    truth_z_muons.phistar = -1;
    truth_z_muons.eta = -1000;

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
    const reco::GenParticle* muon_0 = NULL;
    const reco::GenParticle* muon_1 = NULL;
    const reco::GenParticle* z_boson_electrons = NULL;
    const reco::GenParticle* z_boson_muons = NULL;

    std::vector<const reco::GenParticle*> jpsi;

    //this code is inherited really convoluted, TODO fix it
    for(unsigned int i = 0; i < mc_particles->size(); ++i) {
      const reco::GenParticle* gen_particle = &mc_particles->at(i);
      // Is a Z
      //if (gen_particle->pdgId() == ZBOSON && z_boson_electrons == NULL) {
      if (gen_particle->pdgId() == ZBOSON && z_boson_electrons == NULL) {
        for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
          if (gen_particle->daughter(j)->pdgId() == ELECTRON) {
            z_boson_electrons = gen_particle;
            break;
          }
        }
        // Is an electron
        // TODO - fix this, should just ask for the daugher particles of the z_boson
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

      //Z->muons 
      if (gen_particle->pdgId() == ZBOSON && z_boson_muons == NULL) {
        for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
          if (gen_particle->daughter(j)->pdgId() == MUON) {
            z_boson_muons = gen_particle;
            break;
          }
        }
        // Is an electron
        // TODO - fix this, should just ask for the daugher particles of the z_boson
      } else if (   fabs(gen_particle->pdgId()) == MUON  // In pdgId, fabs(POSITRON) == ELECTRON
          && (muon_0 == NULL || muon_1 == NULL)
          ) {
        for (size_t j = 0; j < gen_particle->numberOfMothers(); ++j) {
          if (gen_particle->mother(j)->pdgId() == ZBOSON) {
            if (muon_0 == NULL) {
              muon_0 = gen_particle;
            } else {
              muon_1 = gen_particle;
            }
          }
        }
      }

      //selecting jpsi that go to at least one muon
      if ( gen_particle->pdgId() == JPSI ) {
        for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
          if (gen_particle->daughter(j)->pdgId() == MUON) {
            jpsi.push_back(gen_particle);
            break;
          }
        }
      }
    }
    for (unsigned int i = 0; i < jpsi.size() ; ++i ) {
      //TODO it seems that jpsi.size should only be 1, do not use a list here?
      //TODO: decide how to handle FSR!
      //if ( jpsi.at(i)->numberOfDaughters() != 2 ) {
      //  std::cout << "JPSI TO LESS THAN 2 DAUGHTERS: " << jpsi.at(i)->numberOfDaughters() << std::endl;
      //  for (unsigned int j = 0 ; j < jpsi.at(i)->numberOfDaughters() ; ++j) {
      //    std::cout << "daughter id: " << jpsi.at(i)->daughter( j )->pdgId() << std::endl;
      //  }
      //  continue;
      //}
      const reco::Candidate * d_mu0 = jpsi.at(i)->daughter( 0 );
      const reco::Candidate * d_mu1 = jpsi.at(i)->daughter( 1 );
      //TODO decide if need muon/antimuon, or can just make this assumption
      if (abs(d_mu0->pdgId()) != MUON && abs(d_mu1->pdgId()) != MUON) {
        std::cout << "JPSI TO NOT MUONS" << std::endl;
        continue;
      }
      if ( d_mu0->pt() >= d_mu1->pt() ) {
        jpsi_muon0.push_back(d_mu0);
        jpsi_muon1.push_back(d_mu1);
      }
      else {
        jpsi_muon0.push_back(d_mu1);
        jpsi_muon1.push_back(d_mu0);
      }
      truth_jpsi.vtx_x.push_back(jpsi.at(i)->vx());
      truth_jpsi.vtx_y.push_back(jpsi.at(i)->vy());
      truth_jpsi.vtx_z.push_back(jpsi.at(i)->vz());
      truth_jpsi.m.push_back( jpsi.at(i)->mass());
      truth_jpsi.pt.push_back( jpsi.at(i)->pt());
      const double JPSIEPP = jpsi.at(i)->energy() + jpsi.at(i)->pz();
      const double JPSIEMP = jpsi.at(i)->energy() - jpsi.at(i)->pz();
      truth_jpsi.y.push_back( 0.5 * log(JPSIEPP / JPSIEMP));
      truth_jpsi.eta.push_back( jpsi.at(i)->eta());
      truth_jpsi.phi.push_back( jpsi.at(i)->phi());


      ////////////////////////////////////////////////////////////////////////////////////
      //TODO figure out if there is a better way to do this, maybe in a function as right now 
      //I do same process for reco_jpsi and truth_jpsi

      const double MUON_MASS = 0.1056583715;
      math::PtEtaPhiMLorentzVector mu0lv(jpsi_muon0.at(i)->pt(), jpsi_muon0.at(i)->eta(), jpsi_muon0.at(i)->phi() , MUON_MASS);
      math::PtEtaPhiMLorentzVector mu1lv(jpsi_muon1.at(i)->pt(), jpsi_muon1.at(i)->eta(), jpsi_muon1.at(i)->phi() , MUON_MASS);
      math::PtEtaPhiMLorentzVector jpsi_lv;
      jpsi_lv = mu0lv + mu1lv;

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector beta_vector;
      beta_vector = jpsi_lv.BoostToCM();

      TLorentzVector *mu0_lv2 = new TLorentzVector(mu0lv.px(), mu0lv.py(), mu0lv.pz(), mu0lv.energy());
      mu0_lv2->Boost(beta_vector.x(), beta_vector.y(), beta_vector.z() );

      TLorentzVector *mu1_lv2 = new TLorentzVector(mu1lv.px(), mu1lv.py(), mu1lv.pz(), mu1lv.energy());
      mu1_lv2->Boost(beta_vector.x(), beta_vector.y(), beta_vector.z() );

      double px = jpsi_lv.px();
      double py = jpsi_lv.py();
      double pz = jpsi_lv.pz();

      double mu0_px_boosted =  mu0_lv2->Px();
      double mu0_py_boosted =  mu0_lv2->Py();
      double mu0_pz_boosted =  mu0_lv2->Pz();

      double mu1_px_boosted =  mu1_lv2->Px();
      double mu1_py_boosted =  mu1_lv2->Py();
      double mu1_pz_boosted =  mu1_lv2->Pz();

      //TODO do this with a function, instead of by hand, cos(theta) = dot_product / (a.len() * b.len() )

      double dot_product_mu0 = px * mu0_px_boosted + py * mu0_py_boosted + pz * mu0_pz_boosted;
      double dot_product_mu1 = px * mu1_px_boosted + py * mu1_py_boosted + pz * mu1_pz_boosted;

      double jpsi_p_mag = pow((px * px + py * py + pz * pz), 0.5);
      double mu0_p_mag = pow((mu0_px_boosted * mu0_px_boosted + mu0_py_boosted * mu0_py_boosted + mu0_pz_boosted * mu0_pz_boosted), 0.5);
      double mu1_p_mag = pow((mu1_px_boosted * mu1_px_boosted + mu1_py_boosted * mu1_py_boosted + mu1_pz_boosted * mu1_pz_boosted), 0.5);

      double cos_jpsi_mu_plus = -1000;
      double cos_jpsi_mu_minus = -1000;
      if (jpsi_muon0.at(i)->charge() == 1 ) {
        cos_jpsi_mu_plus = dot_product_mu0 / (jpsi_p_mag * mu0_p_mag);
        cos_jpsi_mu_minus = dot_product_mu1 / (jpsi_p_mag * mu1_p_mag);
      }
      else {
        cos_jpsi_mu_plus = dot_product_mu1 / (jpsi_p_mag * mu1_p_mag);
        cos_jpsi_mu_minus = dot_product_mu0 / (jpsi_p_mag * mu0_p_mag);
      }

      truth_jpsi.cos_jpsi_mu_plus.push_back(cos_jpsi_mu_plus);
      truth_jpsi.cos_jpsi_mu_minus.push_back(cos_jpsi_mu_minus);

      ////////////////////////////////////////////////////////////////////////////////////



      //TODO 1.2 should not be hardcoded
      if (jpsi_muon0.at(i)->pt() > jpsi_muon1.at(i)->pt() ) {
        if (jpsi_muon0.at(i)->pt() >= MIN_JPSI_LEADING_MUON_PT && jpsi_muon1.at(i)->pt() >= MIN_JPSI_SUBLEADING_MUON_PT 
            && jpsi.at(i)->pt() >= MIN_JPSI_PT) {
          truth_jpsi.has_high_pt_muons.push_back( true );
        }
        else if (jpsi_muon0.at(i)->pt() >= MIN_JPSI_LEADING_MUON_PT && jpsi_muon1.at(i)->pt() >= MIN_JPSI_SUBLEADING_MUON_PT_HIGH_ETA
            && jpsi.at(i)->pt() >= MIN_JPSI_PT && (fabs(jpsi_muon1.at(i)->eta()) >= 1.2)) {
          truth_jpsi.has_high_pt_muons.push_back( true );
        }
        else {
          truth_jpsi.has_high_pt_muons.push_back( false );
        }
      }
      else {
        if (jpsi_muon1.at(i)->pt() >= MIN_JPSI_LEADING_MUON_PT && jpsi_muon0.at(i)->pt() >= MIN_JPSI_SUBLEADING_MUON_PT 
            && jpsi.at(i)->pt() >= MIN_JPSI_PT) {
          truth_jpsi.has_high_pt_muons.push_back( true );
        }
        else if (jpsi_muon1.at(i)->pt() >= MIN_JPSI_LEADING_MUON_PT && jpsi_muon0.at(i)->pt() >= MIN_JPSI_SUBLEADING_MUON_PT_HIGH_ETA
            && jpsi.at(i)->pt() >= MIN_JPSI_PT && (fabs(jpsi_muon0.at(i)->eta()) >= 1.2)) {
          truth_jpsi.has_high_pt_muons.push_back( true );
        }
        else {
          truth_jpsi.has_high_pt_muons.push_back( false );
        }
      }
      if (fabs(jpsi_muon0.at(i)->eta()) <= MAX_JPSI_MUON_ETA && fabs(jpsi_muon1.at(i)->eta()) <= MAX_JPSI_MUON_ETA ) {
        truth_jpsi.has_muons_in_eta_window.push_back( true );
      }
      else {
        truth_jpsi.has_muons_in_eta_window.push_back( false );
      }

      if (jpsi.at(i)->mass() <= MAX_JPSI_MASS && jpsi.at(i)->mass() >= MIN_JPSI_MASS ) {
        truth_jpsi.is_within_jpsi_mass_window.push_back( true );
      }
      else {
        truth_jpsi.is_within_jpsi_mass_window.push_back( false );
      }
    }
    
    //Set truth jpsi cut flags
    found_truth_jpsi_with_high_pt_muons = false;
    for (unsigned int i = 0; i < truth_jpsi.m.size() ; ++i ) {
      if (truth_jpsi.has_high_pt_muons.at(i) && truth_jpsi.is_within_jpsi_mass_window.at(i)) {
        found_truth_jpsi_with_high_pt_muons = true;
      }
    }

    // Continue only if all particles have been found
    if (z_boson_electrons != NULL && electron_0 != NULL && electron_1 != NULL) {
      // We set electron_0 to the higher pt electron
      if (electron_0->pt() < electron_1->pt()) {
        std::swap(electron_0, electron_1);
      }

      const reco::Candidate * d_e0 = z_boson_electrons->daughter( 0 );
      const reco::Candidate * d_e1 = z_boson_electrons->daughter( 1 );
      if (abs(d_e0->pdgId()) != ELECTRON && abs(d_e1->pdgId()) != ELECTRON) {
        std::cout << "Z TO NOT ELECTRONS" << std::endl;
      }
      if ( d_e0->pt() >= d_e1->pt() ) {
        z_truth_electron0 = d_e0;
        z_truth_electron1 = d_e1;
      }
      else {
        z_truth_electron0 = d_e1;
        z_truth_electron1 = d_e0;
      }

      // Add electrons
      ZFinderElectron* zf_electron_0 = AddTruthElectron(*electron_0);
      set_e0_truth(zf_electron_0);
      ZFinderElectron* zf_electron_1 = AddTruthElectron(*electron_1);
      set_e1_truth(zf_electron_1);

      // Z Properties
      truth_z_electrons.m = z_boson_electrons->mass();
      truth_z_electrons.pt = z_boson_electrons->pt();
      const double ZEPP = z_boson_electrons->energy() + z_boson_electrons->pz();
      const double ZEMP = z_boson_electrons->energy() - z_boson_electrons->pz();
      truth_z_electrons.y = 0.5 * log(ZEPP / ZEMP);
      truth_z_electrons.phistar = ReturnPhistar(electron_0->eta(), electron_0->phi(), electron_1->eta(), electron_1->phi());
      truth_z_electrons.eta = z_boson_electrons->eta();
    }
    if (z_boson_muons != NULL && muon_0 != NULL && muon_1 != NULL) {
      // We set muon_0 to the higher pt muon
      if (muon_0->pt() < muon_1->pt()) {
        std::swap(muon_0, muon_1);
      }

      const reco::Candidate * d_mu0 = z_boson_muons->daughter( 0 );
      const reco::Candidate * d_mu1 = z_boson_muons->daughter( 1 );
      if (abs(d_mu0->pdgId()) != MUON && abs(d_mu1->pdgId()) != MUON) {
        std::cout << "Z TO NOT MUONS" << std::endl;
      }
      if ( d_mu0->pt() >= d_mu1->pt() ) {
        z_truth_muon0 = d_mu0;
        z_truth_muon1 = d_mu1;
      }
      else {
        z_truth_muon0 = d_mu1;
        z_truth_muon1 = d_mu0;
      }

      // Z Properties
      truth_z_muons.m = z_boson_muons->mass();
      truth_z_muons.pt = z_boson_muons->pt();
      const double ZEPP = z_boson_muons->energy() + z_boson_muons->pz();
      const double ZEMP = z_boson_muons->energy() - z_boson_muons->pz();
      truth_z_muons.y = 0.5 * log(ZEPP / ZEMP);
      truth_z_muons.phistar = ReturnPhistar(muon_0->eta(), muon_0->phi(), muon_1->eta(), muon_1->phi());
      truth_z_muons.eta = z_boson_muons->eta();
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
    //TESTING JPSI trigger
    //for (unsigned int i = 0; i < reco_jpsi.m.size() ; ++i ) {
    //  //const trigger::TriggerObject* trig_obj_jpsimuon0 = GetBestMatchedTriggerObject(iEvent, JPSI_TRIGGER, reco_jpsi.muon0.at(i).eta(), reco_jpsi.muon0.at(i).phi());
    //  //const trigger::TriggerObject* trig_obj_jpsimuon1 = GetBestMatchedTriggerObject(iEvent, JPSI_TRIGGER, reco_jpsi.muon1.at(i).eta(), reco_jpsi.muon1.at(i).phi());
    //  const trigger::TriggerObject* trig_obj_jpsimuon0 = GetBestMatchedTriggerObject(iEvent, JPSI_TRIGGER, reco_jpsi.muon0.at(i).eta(), reco_jpsi.muon0.at(i).phi());
    //  const trigger::TriggerObject* trig_obj_jpsimuon1 = GetBestMatchedTriggerObject(iEvent, JPSI_TRIGGER, reco_jpsi.muon1.at(i).eta(), reco_jpsi.muon1.at(i).phi());
    //  if (trig_obj_jpsimuon0 != NULL) {
    //    std::cout << "jpsimuon0 trig_obj pt: " << trig_obj_jpsimuon0->pt() << std::endl;
    //  }
    //  if (trig_obj_jpsimuon1 != NULL) {
    //    std::cout << "jpsimuon1 trig_obj pt: " << trig_obj_jpsimuon1->pt() << std::endl;
    //  }
    //}

    for (unsigned int i = 0; i < reco_jpsi.m.size() ; ++i ) {
      const trigger::TriggerObject* trig_obj_jpsimuon0 = GetBestMatchedTriggerObject(iEvent, JPSI_TRIGGER, reco_jpsi.muon0.at(i).eta(), reco_jpsi.muon0.at(i).phi());
      const trigger::TriggerObject* trig_obj_jpsimuon1 = GetBestMatchedTriggerObject(iEvent, JPSI_TRIGGER, reco_jpsi.muon1.at(i).eta(), reco_jpsi.muon1.at(i).phi());
      if (trig_obj_jpsimuon0 != NULL) {
        reco_jpsi.trigger_object_mu0_pt.push_back(trig_obj_jpsimuon0->pt());
      }
      else {
        reco_jpsi.trigger_object_mu0_pt.push_back( -9000.0);
      }
      if (trig_obj_jpsimuon1 != NULL) {
        reco_jpsi.trigger_object_mu1_pt.push_back(trig_obj_jpsimuon1->pt());
      }
      else {
        reco_jpsi.trigger_object_mu1_pt.push_back( -9000.0);
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
        cout << " Truth Z Mass " << truth_z_electrons.m << endl;
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

  double ZFinderEvent::JpsiMuonTruthMatch(const reco::Muon &muon)
  {
    if (is_real_data) {
      return 9999.0;
    }
    double dr0 = 9999.;
    double dr1 = 9999.;
    //TODO decide what to do about case when have 2 Jpsi MC truth, will this ever happen? 

    if (jpsi_muon0.size() == 1 && jpsi_muon0.at(0)->charge() == muon.charge()) {
      dr0 = deltaR( jpsi_muon0.at(0)->eta(), jpsi_muon0.at(0)->phi(), muon.eta(), muon.phi());
    }
    if (jpsi_muon1.size() == 1 && jpsi_muon1.at(0)->charge() == muon.charge()) {
      dr1 = deltaR( jpsi_muon1.at(0)->eta(), jpsi_muon1.at(0)->phi(), muon.eta(), muon.phi());
    }
    if ( dr0 <= dr1 ) {
      return dr0;
    }
    else {
      return dr1;
    }
  }

  double ZFinderEvent::JpsiMuonZMuonMatch(const reco::Muon &muon)
  {
    if (reco_z_from_muons.m == -1) {
      return 9999.0;
    }
    double dr0 = 9999.;
    double dr1 = 9999.;

    if (muon.charge() == z_muon0.charge() ) {
      dr0 = deltaR( z_muon0.eta(), z_muon0.phi(), muon.eta(), muon.phi());
    }
    if (muon.charge() == z_muon1.charge() ) {
      dr1 = deltaR( z_muon1.eta(), z_muon1.phi(), muon.eta(), muon.phi());
    }
    if ( dr0 <= dr1 ) {
      return dr0;
    }
    else {
      return dr1;
    }
  }


  //From RecoVertex / PrimaryVertexProducer / src / VertexHigherPtSquared.cc
  double ZFinderEvent::sumPtSquared(const reco::Vertex &v)
  {
    double sum = 0.;
    double pT;
    for (reco::Vertex::trackRef_iterator it = v.tracks_begin(); it != v.tracks_end(); it++) {
      pT = (**it).pt();
      double epT=(**it).ptError();
      pT=pT>epT ? pT-epT : 0;
      sum += pT*pT;
    }
    return sum;
  }

  ZFinderEvent::~ZFinderEvent() {
    // Clean up all the heap variables we have declared
    for (auto& i_elec : reco_electrons_) {
      delete i_elec;
    }
    for (auto& i_elec : hlt_electrons_) {
      delete i_elec;
    }
    for (auto& i_elec : truth_electrons_) {
      delete i_elec;
    }
  }
}  // namespace zf
