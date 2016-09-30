# ZFinder

## Setting up a working area

To set up a working area, run the following commands:

```bash
scramv1 project -n CMSSW_5_3_28_ZPhysics CMSSW CMSSW_5_3_28
cd CMSSW_5_3_28_ZPhysics/

git clone -b JPsi --single-branch https://github.com/jturkewitz/ZFinder.git src/
cd src
cmsenv
scram build -j 4
```
Then to run the analysis on for example a test file
```bash
cmsRun ZFinder/Event/zfinder_cfg.py
```
Note that you'll have to change the filename to the file you would like to run over.

The heart of the analysis is in ZFinder/Event/
Specifically ZFinder/Event/src/ZFinderEvent.cc applies the selection criteria to search for Z and JPsi particles decaying leptonically in the event. The selection criteria are specified in ZFinder/Event/interface/ZFinderCuts.h. Also important are ZFinder/Event/src/ZFinder.cc and  ZFinder/Event/src/ZFinderPlotter.cc which together create histograms at various cut levels (e.g. a Z in the event, a Z and a JPsi in the event). The tree which is used for performing the 2d fit is created in this step as well, by ZFinder/Event/src/ZFinderTree.cc.

To run these over condor for the ZJpsi skim I created, I have created a bash script 
```bash
sh submit_jobs2.sh
```
You can edit where the files appear in this script. You will then have to hadd the files together, ahadd is a useful tool for speeding up the hadd process.
```bash
condor_filelist.perl ZFinder/Event/zfinder_dimuon8_cfg.py Metadata/filelists/MuOniaParked2012B.txt --batch=1 --jobname=jpsiTest644_MuOniaPartial2012B_dimuon8_eta21_ptsub21_vertex20p0_prompt_inverted_0p3toInf --xrootd

condor_filelist.perl ZFinder/Event/python/zfinder_mc_cfg.py Metadata/filelists/Summer12_DR53X_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_AODSIM.txt --batch=1 --jobname=jpsiTest580_DYJetsToLL --xrootd

condor_filelist.perl ZFinder/Event/python/zfinder_mc_cfg.py Metadata/filelists/jpsi_mc_generated_privately.txt --batch=1 --jobname=jpsiTest635_Jpsi_MC_Private2_1_eta2_1

```
Above are examples of how to run over the MuOniaParked2012B dataset for inclusive J/psi, the DYJetsToLL for Z MC, and the inclusive J/psi MC I generated. Note that you need a grid certificate and xrootd to run on crab and have the results be transferred to Minnesota.

The plots from the histograms are made by macros in ZFinder/Event/scripts/.

These histograms and/or trees are then used by the scripts in RooScripts/scripts/, and the unbinned fit is performed in  RooScripts/scripts/ and RooScripts/scripts/UnbinnedFit/
```bash
cd Rooscripts/scripts
make
./RooFit4LeptonMass.exe ...
```
To run the unbinned fit, there is a 2 step process. The first step is to create a smaller tree than the one produced by running ZFinder/Event/zfinder_cfg.py. This is done by:
```bash
root -l
 .L loop.C
Run("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/jpsiTest650_Z_40_300_DoubleMuon_2012.root",1,1)
Run("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/jpsiTest651_Z_40_300_DoubleElectron_2012.root",0,1)

Run("/data/whybee0a/user/turkewitz_2/test/turkewitz/TestFiles/jpsiTest615_MuOniaPartial2012B_dimuon8_eta21.root",0,0)
```
Then hadd the output of the doublemuon zjpsi skim and the double electron zjpsi skim together, which will be the input to the next step. You should be using the filename of the zjpsi skim you created, I just showed some example filenames here. The thrid line makes the inclusive J/psi subset, in general you won't need to redo this.

Finally, to do the unbinned fit:
```bash
.x sPlotFit.C("ZJpsi_40_300_PU0p0to0p5_acc_mod.root",0,"image_directory_path_name",0,0) > log_file_path
```
You may need to adjust the last 2 options to this fit, which control whether you run over just Z to muons, just Z to electrons, or both. Sometimes you will have to run the analysis for certain selection criteria, to do this just add a TSelectionCUt (or uncomment one) to the dataset you want to apply the cut to. The relevant outputs are stored in the log file.


If running the analysis more than once, it is useful to skim the double electron and double muon datasets, which can be done with:

JPsiFilter/JPsiFilter/python/jpsi_z_to_ee_filter_cfg.py
JPsiFilter/JPsiFilter/python/jpsi_z_to_mumu_filter_cfg.py
And running via crab (or condor with --xrootd at Minnesota). This will filter out most events, keeping events that are able to have both a Z and a JPsi to mumu. Note that I have already done this step, so unless you want to reskim events it should not be necessary.
