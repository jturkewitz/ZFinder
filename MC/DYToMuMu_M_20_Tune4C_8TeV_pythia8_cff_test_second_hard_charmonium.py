# Auto generated configuration file
# using: 
# Revision: 1.381.2.28 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: DYToMuMu_M_20_Tune4C_8TeV_pythia8_cff_test.py --conditions auto:startup_7E33v4 -n 10 --eventcontent RAWSIM --step GEN,SIM --datatier GEN-SIM --beamspot Realistic8TeVCollision --datamix NODATAMIXER --pileup NoPileUp --fileout file:/local/cms/user/turkewitz/tmp/Jpsi_MM_step1_test20.root --no_exec
import FWCore.ParameterSet.Config as cms
import random

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

x = random.randrange(1,100000000)
process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(x)


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    ##input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('PYTHIA8 Z/gamma* to mumu, M(mu+mu-) > 20 GeV at sqrt(s) = 8TeV, Tune 4c'),
    name = cms.untracked.string('$Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/GenProduction/python/EightTeV/DYToMuMu_M_20_Tune4C_8TeV_pythia8_cff.py,v $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:/data/whybee0a/user/turkewitz_2/test/turkewitz/temp/MC/DYToMuMu_M_20_Tune4C_8TeV_pythia8_step1.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup_7E33v4', '')

process.generator = cms.EDFilter("Pythia8GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(8000.0),
    crossSection = cms.untracked.double(1300),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring('Main:timesAllowErrors = 10000', 
            'WeakSingleBoson:ffbar2gmZ = on', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10.', 
            'PartonLevel:MPI = on', 
            'SecondHard:generate = on', 
            'SecondHard:TwoJets = off', 
            'SecondHard:Charmonium = on', 
            'PhaseSpace:sameForSecond = off', 
            'PhaseSpace:mHatMin = 20', 
            'PhaseSpace:mHatMinSecond = 2.5', 
            'PhaseSpace:pTHatMinSecond = 8.0', 
            '23:onMode = off', 
            '23:onIfAny = 13', 
            '443:onMode = off', 
            '443:onIfAny = 13', 
            'Tune:pp 5'),
        parameterSets = cms.vstring('processParameters')
    )
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

