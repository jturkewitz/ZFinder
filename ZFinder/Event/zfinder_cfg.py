import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# Set up message output and logging
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' ##data
#process.GlobalTag.globaltag = 'START53_V29B::All' ##mc

#process.MessageLogger.cerr.FwkReport.reportEvery = 100  # Report status ever 100 events
process.MessageLogger.cerr.FwkReport.reportEvery = 1000  # Report status ever 100 events

# Number of events from each file to process. It should be -1 (all) when
# running for an analysis
N_EVENTS_TO_PROCESS = -1
if N_EVENTS_TO_PROCESS != -1:
    print "NOT RUNNING ON ALL EVENTS IN THE FILE!"
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(N_EVENTS_TO_PROCESS)
        )

process.source = cms.Source("PoolSource",
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/JpsiSkim/jpsiSkimUpdated.root')
    ##fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/ZPhysics/JPsiSkim/DoubleElectron/2012B/jpsiSkimUpdated/jpsiSkimUpdated_999-pool.root')
    fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/ZPhysics/JPsiSkim/MuOnia/Run2012B/jpsiTriggerSkim_198_1_U0I.root')
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/PromptJpsi/JpsiMM_8TeV_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_RECO.root')
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/tmp/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_00037C53-AAD1-E111-B1BE-003048D45F38.root')
    ## fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/ZPhysics/JPsiSkim/DoubleElectron/jpsiSkim/jpsiSkim_000-pool.root')
    ##fileNames = cms.untracked.vstring( 'file:/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/JPsiFilter/jpsiSkim.root')
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/phedex/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/20000/003EC246-5E67-E211-B103-00259059642E.root')
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("test9d10.root")
        )

#
# rho value for isolation
#

from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets  # the 4 references the rParam = 0.4
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

#
# particle flow isolation
#

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.pfiso = cms.Sequence(process.pfParticleSelectionSequence + process.eleIsoSequence)

##for jpsi MuOnia triggerign
##process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
##    ##triggerConditions = cms.vstring(
##    ##  'HLT_Dimuon0_Jpsi_v*',
##    ##  'HLT_Dimuon8_Jpsi_v*',
##    ##  'HLT_Dimuon10_Jpsi_v*'),
##    triggerConditions = cms.vstring(
##      'HLT_Dimuon0_Jpsi_v*'),
##    hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
##    l1tResults = cms.InputTag( "" ),
##    l1tIgnoreMask = cms.bool( False ),
##    l1techIgnorePrescales = cms.bool( False ),
##    daqPartitions = cms.uint32( 1 ),
##    throw = cms.bool( False )
##    ##throw = cms.bool( True )
##)

#
# ZFinder
#

# Import ZDefinitions
from ZFinder.Event.ZDefinitions_cfi import zdefs

process.ZFinder = cms.EDAnalyzer('ZFinder',
        # General tags
        ecalElectronsInputTag  = cms.InputTag("gsfElectrons"),
        muonsInputTag          = cms.InputTag("muons"),
        hfElectronsInputTag    = cms.InputTag("hfRecoEcalCandidate"),
        hfClustersInputTag     = cms.InputTag("hfEMClusters"),
        conversionsInputTag    = cms.InputTag("allConversions"),
        beamSpotInputTag       = cms.InputTag("offlineBeamSpot"),
        rhoIsoInputTag         = cms.InputTag("kt6PFJetsForIsolation", "rho"),
        primaryVertexInputTag  = cms.InputTag("offlinePrimaryVertices"),
        ntElectronsInputTag    = cms.InputTag("photons"),
        ak5PFJetsInputTag      = cms.InputTag("ak5PFJets"),
        isoValInputTags        = cms.VInputTag(
            cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
            cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
            cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')
            ),
        # MC, but still required to be something for data
        pileupInputTag = cms.InputTag("addPileupInfo"),
        generatorInputTag = cms.InputTag("genParticles"),
        # ZDefinitions from ZFinder.ZFinder.ZDefinitions_cfi
        ZDefinitions = zdefs
        )

# RUN
process.p = cms.Path(process.kt6PFJetsForIsolation * process.pfiso * process.ZFinder)
process.schedule = cms.Schedule(process.p)
