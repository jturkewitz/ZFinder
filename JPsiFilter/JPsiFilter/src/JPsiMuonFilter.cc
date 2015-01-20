// -*- C++ -*-
//
// Package:    JPsiMuonFilter
// Class:      JPsiMuonFilter
// 
/**\class JPsiMuonFilter JPsiMuonFilter.cc JPsiMuonFilter/JPsiMuonFilter/src/JPsiMuonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jared Turkewitz
//         Created:  Tue May  6 18:57:10 CDT 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "DataFormats/MuonReco/interface/MuonFwd.h" // reco::MuonCollection
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"  // reco::PhotonCollection
#include "DataFormats/EgammaCandidates/interface/Photon.h"  // reco::Photon
#include "DataFormats/MuonReco/interface/Muon.h" // reco::Muon

//
// class declaration
//

class JPsiMuonFilter : public edm::EDFilter {
   public:
      explicit JPsiMuonFilter(const edm::ParameterSet&);
      ~JPsiMuonFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
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
JPsiMuonFilter::JPsiMuonFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


JPsiMuonFilter::~JPsiMuonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
JPsiMuonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   bool found_four_jpsi_muons = false;
   bool found_two_z_muons = false;
   bool found_z_and_jpsi = false;
   //bool found_two_electrons = false;
/*
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/
   edm::Handle<reco::MuonCollection> muons_h;
   iEvent.getByLabel("muons", muons_h);
   if ( !iEvent.getByLabel("muons", muons_h) ) {
     std::cout << "ERROR: Could not find muons" << std::endl ;
     return false ;
   }
   //photon collection is not filtered for electrons
   //edm::Handle<reco::PhotonCollection> photons_h;
   //edm::Handle<reco::GsfElectronCollection> els_h;
   //if ( !iEvent.getByLabel("gsfElectrons", els_h) ) {
   //  std::cout << "ERROR: Could not find photons" << std::endl ;
   //  return false ;
   //}
   //TODO clean up code and decide on good min pt values
   // and any other cuts that are needed

   //int n_electrons = 0;
   //int n_positrons = 0;
   //for(unsigned int i = 0; i < els_h->size(); ++i) {
   //  reco::GsfElectron electron = els_h->at(i);
   //  if ( electron.pt() > 15 ) {
   //    if ( electron.charge() == -1) {
   //        n_electrons++;
   //    }
   //    if ( electron.charge() == 1) {
   //        n_positrons++;
   //    }
   //  }
   //}
   //if ( n_electrons >= 1 && n_positrons >=1 ) {
   //  found_two_electrons = true;
   //}

   int n_jpsi_muons = 0;
   int n_jpsi_antimuons = 0;
   int n_z_muons = 0;
   int n_z_antimuons = 0;
   for(unsigned int i = 0; i < muons_h->size(); ++i) {
     reco::Muon muon = muons_h->at(i);
     //check for muons forming a Z
     if ( muon.pt() > 18 ) {
       if ( muon.charge() == -1) {
           n_z_muons++;
       }
       if ( muon.charge() == 1) {
           n_z_antimuons++;
       }
     }

     //check for muons forming a J/Psi, note that we want a j/psi independent from the z,
     //so crudely we require 4 "j/psi" muons as the j/psi pT requirement is looser
     if ( muon.pt() > 2.5 ) {
       if ( muon.charge() == -1) {
           n_jpsi_muons++;
       }
       if ( muon.charge() == 1) {
           n_jpsi_antimuons++;
       }
     }

   }
   if ( n_jpsi_muons >= 2 && n_jpsi_antimuons >= 2 ) {
     found_four_jpsi_muons = true;
   }
   if (n_z_muons >= 1 && n_z_antimuons >= 1 ) {
     found_two_z_muons = true;
   }
   if (found_four_jpsi_muons && found_two_z_muons ) {
     found_z_and_jpsi = true;
   }

   return (found_z_and_jpsi);
}

// ------------ method called once each job just before starting event loop  ------------
void 
JPsiMuonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiMuonFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
JPsiMuonFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
JPsiMuonFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
JPsiMuonFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
JPsiMuonFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JPsiMuonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(JPsiMuonFilter);
