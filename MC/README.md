#Generating MC

## Generating Pile-up file
I generated pile-up file with the command:
```bash
cmsDriver.py MinBias_TuneZ2star_8TeV_pythia6_cff.py --conditions auto:mc -n 10000 --eventcontent RAWSIM --step GEN,SIM --datatier GEN-SIM --beamspot Realistic8TeVCollision --fileout file:/data/whybee0a/user/turkewitz_2/test/turkewitz/minbias_pileup8.root
```
This pile-up file is used to add pile-up to the events.

Then I use 3 steps to generate MC reconstructed with CMS - this process was found by trial/error as CMSSW commands to generate MC lack up-to-date documentation. These commands generate DPS (double parton scattering) events with a Z to mumu and a second hard charmonium. Note that this is DPS and not the signal process, but still this MC can be used (with scale factors) to determine the acceptance times efficiency of the J/Psi. 
```bash
condor_run_mc.perl DYToMuMu_M_20_Tune4C_8TeV_pythia8_cff_test_second_hard_charmonium.py count2.txt --batch=1 --jobname=jpsi_mc_pythia8_Test5_second_hard
condor_run_mc.perl step2_DIGI_L1_DIGI2RAW_HLT_PU.py step1_files_second_hard_test5.txt --batch=1 --jobname=jpsi_mc_pythia8_Test5_second_hard_step2_b
condor_run_mc.perl step3_RAW2DIGI_L1Reco_RECO_PU.py step2_files_second_hard_test5_b.txt --batch=1 --jobname=jpsi_mc_pythia8_Test5_step3_b
```
Note that the the condor_run_mc.perl just allows for batch processing on the Minnesota condor cluster. The important physics is contained in DYToMuMu_M_20_Tune4C_8TeV_pythia8_cff_test_second_hard_charmonium.py.

