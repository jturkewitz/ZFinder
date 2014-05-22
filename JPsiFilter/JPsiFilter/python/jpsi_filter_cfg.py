import FWCore.ParameterSet.Config as cms

process = cms.Process("JpsiSkim")

# Set up message output and logging
process.load("FWCore.MessageService.MessageLogger_cfi")
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
    fileNames = cms.untracked.vstring( 'file:/local/cms/phedex/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/20000/D21CAFBE-6B67-E211-91E1-002618943970.root')
    #fileNames = cms.untracked.vstring( 'file:/local/cms/phedex/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/20000/003EC246-5E67-E211-B103-00259059642E.root')
)

##process.TFileService = cms.Service("TFileService",
##        fileName = cms.string("skim.root")
##        )

process.jpsiFilter = cms.EDFilter("JPsiFilter")

# RUN
process.mypath = cms.Path(process.jpsiFilter)
##process.schedule = cms.Schedule(process.p)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('jpsiSkimUpdated4.root'),
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
