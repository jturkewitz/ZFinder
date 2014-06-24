import FWCore.ParameterSet.Config as cms

process = cms.Process("JpsiTriggerSkim")

# Set up message output and logging
#process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' ##data
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 100  # Report status ever 100 events
process.MessageLogger.cerr.FwkReport.reportEvery = 10000  # Report status ever 100 events

# Number of events from each file to process. It should be -1 (all) when
# running for an analysis
N_EVENTS_TO_PROCESS = -1
if N_EVENTS_TO_PROCESS != -1:
    print "NOT RUNNING ON ALL EVENTS IN THE FILE!"
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(N_EVENTS_TO_PROCESS)
        )

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/tmp/Run2012AMuOnia8401CB0B-D634-E211-8247-002618943919.root')
    fileNames = cms.untracked.vstring( 'file:/hdfs/cms/user/turkewitz/ZPhysics/JPsiSkim/DoubleElectron/2012B/jpsiSkimUpdated/jpsiSkimUpdated_999-pool.root')
)
process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(
      'HLT_Dimuon0_Jpsi_v*',
      'HLT_Dimuon8_Jpsi_v*',
      'HLT_Dimuon10_Jpsi_v*'),
    hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( False )
    ##throw = cms.bool( True )
)

##process.TFileService = cms.Service("TFileService",
##        fileName = cms.string("skim.root")
##        )

##process.jpsiFilter = cms.EDFilter("JPsiFilter")

# RUN
process.mypath = cms.Path(process.triggerSelection)
##process.schedule = cms.Schedule(process.p)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('jpsiTriggerSkim.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('mypath')),
                               outputCommands = cms.untracked.vstring(
                                   #"drop *",
                                   "keep *",
                                   #"keep *_*_*_PDF",
                                   #"keep *_genParticles_*_*",
                                   #"keep *_*_*_LHE"
                               )
)
process.output = cms.EndPath(process.out)
