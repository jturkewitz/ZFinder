// -*- C++ -*-
//
// Package:    ZToMuonsSkim
// Class:      ZToMuonsSkim
// 
/**\class ZToMuonsSkim ZToMuonsSkim.cc ZFinder/ZToMuonsSkim/src/ZToMuonsSkim.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jared Turkewitz
//         Created:  Fri Oct 10 15:20:36 CDT 2014
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
#include "DataFormats/MuonReco/interface/Muon.h" // reco::Muon

//
// class declaration
//

class ZToMuonsSkim : public edm::EDFilter {
   public:
      explicit ZToMuonsSkim(const edm::ParameterSet&);
      ~ZToMuonsSkim();

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
ZToMuonsSkim::ZToMuonsSkim(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


ZToMuonsSkim::~ZToMuonsSkim()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ZToMuonsSkim::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   bool found_two_muons = false;
   edm::Handle<reco::MuonCollection> muons_h;
   iEvent.getByLabel("muons", muons_h);
   if ( !iEvent.getByLabel("muons", muons_h) ) {
     std::cout << "ERROR: Could not find muons" << std::endl ;
     return false ;
   }

   int n_muons = 0;
   int n_antimuons = 0;
   for(unsigned int i = 0; i < muons_h->size(); ++i) {
     reco::Muon muon = muons_h->at(i);
     if ( muon.pt() > 18 ) {
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

   return (found_two_muons);
}

// ------------ method called once each job just before starting event loop  ------------
void 
ZToMuonsSkim::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZToMuonsSkim::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
ZToMuonsSkim::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
ZToMuonsSkim::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
ZToMuonsSkim::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
ZToMuonsSkim::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZToMuonsSkim::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ZToMuonsSkim);
