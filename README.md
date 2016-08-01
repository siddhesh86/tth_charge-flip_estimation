# tth_charge-flip_estimation

## Installation
Using combine requires an older version of CMSSW

Instructions for installing combine are taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit#ROOT6_SLC6_release_CMSSW_7_4_X (version 7_4_7 at the time)

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src 
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v6.2.1
```
We also need CombineHarvester:
```bash
cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
cd $CMSSW_BASE/src
scramv1 b clean; scramv1 b -j8
```

Installing the code proceeds as follows:
```bash
git clone git@github.com:HEP-KBFI/tth_charge-flip_estimation.git $CMSSW_BASE/src/tthAnalysis/ChargeFlipEstimation
```

## Running
# Datacard creation
For datacard creation, run:
```bash
ChargeFlipDC
```
# Fit 21 electron pair categories and make plots
```bash
python make_fits.py
```

# Final fit
```bash
FitCF
```
