# Auto generated configuration file
# using: 
# Revision: 1.381.2.13 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step2 --conditions auto:startup_7E33v4 -n 100000 --eventcontent RAWSIM -s DIGI,L1,DIGI2RAW,HLT:7E33v4 --datatier GEN-SIM-DIGI-RAW --pileup 2012_Summer_50ns_PoissonOOTPU --pileup_input file:/local/cms/user/turkewitz/tmp/minbias_pileup.root --filein file:/local/cms/user/turkewitz/tmp/Jpsi_MM_step1_test15.root --fileout file:/local/cms/user/turkewitz/tmp/Jpsi_MM_step2_test15_pileup.root
import FWCore.ParameterSet.Config as cms
import random
import string

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_2012_Summer_50ns_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_7E33v4_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    ##input = cms.untracked.int32(10000)
    ##input = cms.untracked.int32(1000000)
    ##input = cms.untracked.int32(100000)
    ##input = cms.untracked.int32(500000)
    ##input = cms.untracked.int32(10)
    input = cms.untracked.int32(100)
    ##input = cms.untracked.int32(500)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:/data/whybee0a/user/turkewitz_2/test/turkewitz/temp/MC/DYToMuMu_M_20_Tune4C_8TeV_pythia8_step1.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.13 $'),
    annotation = cms.untracked.string('step2 nevts:100000'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:/data/whybee0a/user/turkewitz_2/test/turkewitz/temp/MC/DYToMuMu_M_20_Tune4C_8TeV_pythia8_step2.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

# Additional output definition

# Other statements
x = random.randrange(1,100000000)
str_x = str(x)
histProbFunction_str = 'histProbFunction'
histProbFunction_str+=str_x
histProbFunction_str+='.root'
process.mix.input.nbPileupEvents.histoFileName = cms.untracked.string(histProbFunction_str)
process.mix.input.fileNames = cms.untracked.vstring(['file:/data/whybee0a/user/turkewitz_2/test/turkewitz/minbias_pileup8.root'])
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup_7E33v4', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWSIMoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
