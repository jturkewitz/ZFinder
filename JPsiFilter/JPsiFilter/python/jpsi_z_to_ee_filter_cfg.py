import FWCore.ParameterSet.Config as cms

process = cms.Process("JpsiToMuMuZToeeSkim")

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
    fileNames = cms.untracked.vstring( 'file:/hdfs/cms/phedex/store/data/Run2012B/DoubleElectron/AOD/22Jan2013-v1/20000/FE9BF84B-4269-E211-8A41-003048FFD75C.root')
)

# Run only on lumis specified in the lumi file
# Recipe from:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePythonTips#Use_a_JSON_file_of_good_lumi_sec
from FWCore.ParameterSet.Types import untracked, VLuminosityBlockRange
from FWCore.PythonUtilities.LumiList import LumiList
##json_file for electrons
json_file = "/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/Metadata/lumi_json/Run2012ABCD.json" # File location
run_2012abcd_lumis = LumiList(filename = json_file).getCMSSWString().split(',')

process.source.lumisToProcess = untracked(VLuminosityBlockRange(run_2012abcd_lumis))

process.jpsiFilter = cms.EDFilter("JPsiFilter")

# RUN
process.mypath = cms.Path(process.jpsiFilter)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('jpsiMuMu_Zee_Skim.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('mypath')),
                               outputCommands = cms.untracked.vstring(
                                   #"drop *",
                                   "keep *",
                              )
)
process.output = cms.EndPath(process.out)
