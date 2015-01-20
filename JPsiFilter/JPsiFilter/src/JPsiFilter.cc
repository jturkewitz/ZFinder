// -*- C++ -*-
//
// Package:    JPsiFilter
// Class:      JPsiFilter
// 
/**\class JPsiFilter JPsiFilter.cc JPsiFilter/JPsiFilter/src/JPsiFilter.cc

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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"  // reco::GsfElectron

//
// class declaration
//

class JPsiFilter : public edm::EDFilter {
   public:
      explicit JPsiFilter(const edm::ParameterSet&);
      ~JPsiFilter();

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
JPsiFilter::JPsiFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


JPsiFilter::~JPsiFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
JPsiFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   bool found_two_muons = false;
   bool found_two_electrons = false;
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
   edm::Handle<reco::GsfElectronCollection> els_h;
   if ( !iEvent.getByLabel("gsfElectrons", els_h) ) {
     std::cout << "ERROR: Could not find photons" << std::endl ;
     return false ;
   }
   //TODO clean up code and decide on good min pt values
   // and any other cuts that are needed

   int n_electrons = 0;
   int n_positrons = 0;
   for(unsigned int i = 0; i < els_h->size(); ++i) {
     reco::GsfElectron electron = els_h->at(i);
     if ( electron.pt() > 18 ) {
       if ( electron.charge() == -1) {
           n_electrons++;
       }
       if ( electron.charge() == 1) {
           n_positrons++;
       }
     }
   }
   if ( n_electrons >= 1 && n_positrons >=1 ) {
     found_two_electrons = true;
   }

   int n_muons = 0;
   int n_antimuons = 0;
   for(unsigned int i = 0; i < muons_h->size(); ++i) {
     reco::Muon muon = muons_h->at(i);
     if ( muon.pt() > 2.5 ) {
       if ( muon.charge() == -1) {
           n_muons++;
       }
       if ( muon.charge() == 1) {
           n_antimuons++;
       }
     }
   }
   if ( n_muons >= 1 && n_antimuons >= 1 ) {
     found_two_muons = true;
   }

   return (found_two_electrons && found_two_muons);
}

// ------------ method called once each job just before starting event loop  ------------
void 
JPsiFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
JPsiFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
JPsiFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
JPsiFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
JPsiFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JPsiFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(JPsiFilter);
