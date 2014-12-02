import FWCore.ParameterSet.Config as cms

process = cms.Process("ZToMuMuSkim")

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
    ##fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/ZPhysics/tmp/SingleMu_2012B_test_file.root')
    fileNames = cms.untracked.vstring( 'file:/local/cms/user/turkewitz/tmp/ProblematicZJPsiSkimFiles_try2/0C83A2ED-9FA7-E211-B250-00259073E488_try2.root')
)

# Run only on lumis specified in the lumi file
# Recipe from:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePythonTips#Use_a_JSON_file_of_good_lumi_sec
from FWCore.ParameterSet.Types import untracked, VLuminosityBlockRange
from FWCore.PythonUtilities.LumiList import LumiList
##json_file for muons
json_file = "/home/user1/turkewitz/Work/CMSSW_5_3_13_ZJPsi/src/Metadata/lumi_json/Run2012ABCD_Muons.json" # File location
run_2012abcd_lumis = LumiList(filename = json_file).getCMSSWString().split(',')
process.source.lumisToProcess = untracked(VLuminosityBlockRange(run_2012abcd_lumis))

process.zmumuFilter = cms.EDFilter("ZToMuonsSkim")

# RUN
process.mypath = cms.Path(process.zmumuFilter)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('Zmumu_Skim.root'),
                               # save only events passing the full path
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('mypath')),
                               outputCommands = cms.untracked.vstring(
                                   #"drop *",
                                   "keep *",
                              )
)
process.output = cms.EndPath(process.out)
