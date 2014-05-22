#include "ZFinder/Event/interface/ZFinderMuon.h"

// Standard Library
#include <iostream>  // std::cout, std::endl;

// ZFinder
#include "ZFinder/Event/interface/PDGID.h"  // PDGID enum (ELECTRON, POSITRON, etc.)


namespace zf {
    ZFinderMuon::ZFinderMuon(reco::Muon input_muon) {
        /* Extract the useful quantities from a Muon */
        pt = input_muon.pt();
        eta = input_muon.eta();
        phi = input_muon.phi();
        charge = input_muon.charge();
        vertex = input_muon.vertex();
    }
}
