import FWCore.ParameterSet.Config as cms

process = cms.Process("JpsiMuMu_Zmumu_Skim")

# Set up message output and logging
process.load("FWCore.MessageService.MessageLogger_cfi")
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
)

process.jpsiFilter = cms.EDFilter("JPsiMuonFilter")

# RUN
process.mypath = cms.Path(process.jpsiFilter)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('jpsiMuMu_Zmumu_Skim.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('mypath')),
                               outputCommands = cms.untracked.vstring(
                                   #"drop *",
                                   "keep *",
                              )
)
process.output = cms.EndPath(process.out)
