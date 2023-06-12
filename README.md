# TauFW

Framework for tau analysis using NanoAOD at CMS. Three main packages are
1. [`PicoProducer`](PicoProducer): Tools to process nanoAOD and make custom analysis ntuples.
2. [`Plotter`](Plotter): Tools for further analysis, auxiliary measurements, validation and plotting. [Under development.]
3. [`Fitter`](Fitter): Tools for measurements and fits in combine. [Under development.]

## Installation

### CMSSW environment
First, setup a CMSSW release, for example,
<table>
<tr>
<td> nanoAODv10 </td> <td> older versions </td>
</tr>
<tr>
<td>

```bash
export CMSSW=CMSSW_12_4_8
export SCRAM_ARCH=slc7_amd64_gcc10
cmsrel $CMSSW
cd $CMSSW/src
cmsenv
```
</td>
<td>

```bash
export CMSSW=CMSSW_10_6_13
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel $CMSSW
cd $CMSSW/src
cmsenv
```
</td>
</tr>
</table>

On a Linux 8 node like `lxplus8`, you can use the `el8_amd64_gcc10` architecture instead of `slc7_*`.

Which CMSSW version should matter for post-processing of nanoAOD,
but if you like to use Combine in the same repository,
it is better to use at least the [recommended version](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#setting-up-the-environment-and-installation).

### TauFW
To install `TauFW`:
```
cd $CMSSW_BASE/src/
git clone https://github.com/cms-tau-pog/TauFW TauFW
scram b -j4
```
With each new session, do
```
export SCRAM_ARCH=slc7_amd64_gcc700
cd $CMSSW/src
cmsenv
```

### PicoProducer
If you want to process nanoAOD using `PicoProducer`, install [`NanoAODTools`](https://github.com/cms-nanoAOD/nanoAOD-tools):
```
cd $CMSSW_BASE/src/
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
scram b -j4
```
If you want to use tau ID SF, please install [`TauIDSFs` tool](https://github.com/cms-tau-pog/TauIDSFs):
```
cd $CMSSW_BASE/src
git clone https://github.com/cms-tau-pog/TauIDSFs TauPOG/TauIDSFs
cmsenv
scram b -j4
```

### Fitter and Combine tools
If you want to use the `Combine` tools in `Fitter`, install
[`Combine`](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#setting-up-the-environment-and-installation),
```
cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.1.0 # for CMSSW_10_X
git checkout v9.1.0 # for CMSSW_11_X
```
and then [`CombineHarvester`](https://github.com/cms-analysis/CombineHarvester),
```
cd $CMSSW_BASE/src
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
scramv1 b clean; scramv1 b
git checkout v2.0.0 # for CMSSW_11_X only
```

### TauID Scale Factor (SF) JSON and ROOT files creation
To create JSON files for
[`correctionlib`](https://github.com/cms-nanoAOD/correctionlib),
please follow the instructions
[here](https://gitlab.cern.ch/cms-tau-pog/jsonpog-integration/-/blob/TauPOG_v2/POG/TAU/README4UPDATES.md).
From at least `CMSSW_11_3_X`, `correctionlib` should be pre-installed.


To create ROOT files including the measured SFs please install [`TauIDSFs` tool](https://github.com/cms-tau-pog/TauFW/#picoproducer) as illustrated above.
Modify the `TauIDSFs/utils/createSFFiles.py` script to include your measured SFs into the script. 
Finally, run the `TauFW/scripts/tau_createROOT.sh` to generate your ROOT files. They will be created into `TauFW/scripts/data/`
IMPORTANT: please comment and do not delete older SFs
