#ifndef ZFINDER_ZFINDEREVENT_H_
#define ZFINDER_ZFINDEREVENT_H_

// Standard Library
#include <string>  // string
#include <vector>  // vector

// CMSSW
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"  // reco::GsfElectron
#include "FWCore/Framework/interface/Event.h"  // edm::Event, edm::EventSetup
#include "FWCore/ParameterSet/interface/ParameterSet.h"  // edm::ParameterSet
#include "FWCore/Utilities/interface/InputTag.h"  // edm::InputTag
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"  // GenParticle

// ZFinder
#include "ZFinder/ZFinder/interface/ZFinderElectron.h"

struct BasicRequirements{
    // TODO: Fill in defaults
    double ept_min;
    double ept_max;
};

template <class Particle_T>
class ZFinderEvent{
    public:
        // Constructor. Although iEvent, iSetup, and iConfig violate our naming
        // convention, they are almost ubiquitous in CMSSW code
        ZFinderEvent(
                const edm::Event& iEvent, 
                const edm::EventSetup& iSetup, 
                const edm::ParameterSet& iConfig,
                const bool use_truth
                );

        // Data or MC
        bool is_real_data;

        // Beam Spot
        struct Beamspot{
            double x;
            double y;
            double z;
        } bs;

        // Primary vertexes
        struct Vertexes{
            unsigned int num;
            double x;
            double y;
            double z;
        } vert;

        // Event ID
        struct EventID{
            unsigned int run_num;
            unsigned int lumi_num;
            unsigned int event_num;
        } id;

        // Z Data
        struct ZData{
            double m;
            double pt;
            double y;
            double phistar;
        } z;

        // These are the special, selected electrons
        ZFinderElectron<Particle_T>* e0;
        ZFinderElectron<Particle_T>* e1;
        void set_e0(ZFinderElectron<Particle_T>* electron) { e0 = electron; }
        void set_e1(ZFinderElectron<Particle_T>* electron) { e1 = electron; }
        void set_both_e(ZFinderElectron<Particle_T>* electron0, ZFinderElectron<Particle_T>* electron1) { e0 = electron0; e1 = electron1; }

    protected:
        // Called by the constructor to handle MC and Data separately
        void InitReco(const edm::Event& iEvent, const edm::EventSetup& iSetup, const BasicRequirements& cuts);
        void InitTruth(const edm::Event& iEvent, const edm::EventSetup& iSetup, const BasicRequirements& cuts);

        void InitRecoElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup, const BasicRequirements& cuts);
        void InitZ();

        // Input tags
        struct InputTags{
            edm::InputTag electron;
            edm::InputTag conversion;
            edm::InputTag beamspot;
            edm::InputTag rho_iso;
            edm::InputTag vertex;
            edm::InputTag pileup;
            edm::InputTag generator;
            std::vector<edm::InputTag> iso_vals;
        } inputtags_;

        // A list of all electrons
        std::vector<ZFinderElectron<Particle_T> > electrons_;
        ZFinderElectron<Particle_T>* AddElectron(reco::GsfElectron electron);
        ZFinderElectron<Particle_T>* AddElectron(HepMC::GenParticle electron);

        // Calculate phistar
        double ReturnPhistar(const double& eta0, const double& phi0, const double& eta1, const double& phi1);

        // Sorting functions
        static bool SortByPTHighLow(const ZFinderElectron<Particle_T> e0, const ZFinderElectron<Particle_T> e1) { return (e0.pt > e1.pt); }
};
#endif  // ZFINDER_ZFINDEREVENT_H_
