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
