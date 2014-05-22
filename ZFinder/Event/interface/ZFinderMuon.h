#ifndef ZFINDER_ZFINDERMUON_H_
#define ZFINDER_ZFINDERMUON_H_

// Standard Library
#include <string>  // std::string
#include <map>  // std::map
#include <vector>  // std::vector

// CMSSW
#include "DataFormats/MuonReco/interface/Muon.h" // reco::Muon

namespace zf {

    class ZFinderMuon {
        public:
            ZFinderMuon() {};
            ZFinderMuon(reco::Muon input_muon);

            // Kinematics variables
            double pt;
            double eta;
            double phi;
            math::XYZPoint vertex;

            // Other physical properties
            int charge;
    };
}  // namespace zfe
#endif  // ZFINDER_ZFINDERMUON_H_
