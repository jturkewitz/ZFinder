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
    fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/tmp/SingleMu_2012B_test_file.root')
    #fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/tmp/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_00037C53-AAD1-E111-B1BE-003048D45F38.root')
    #fileNames = cms.untracked.vstring( 'file:/local/cms/phedex/store/data/Run2012A/DoubleElectron/AOD/22Jan2013-v1/20000/003EC246-5E67-E211-B103-00259059642E.root')
)

##process.TFileService = cms.Service("TFileService",
##        fileName = cms.string("skim.root")
##        )

##TODO make another filter, one for muons, one for electrons
##process.jpsiFilter = cms.EDFilter("JPsiFilter")
process.jpsiFilter = cms.EDFilter("JPsiMuonFilter")

# RUN
process.mypath = cms.Path(process.jpsiFilter)
##process.schedule = cms.Schedule(process.p)

process.out = cms.OutputModule("PoolOutputModule",
                               #fileName = cms.untracked.string('jpsiSkimMCUpdated.root'),
                               fileName = cms.untracked.string('jpsiSkimMuonsUpdated.root'),
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
