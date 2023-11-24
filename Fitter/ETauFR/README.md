# TauFW Fitter ETauFR
## Introduction
The module is intended to evaluate the ETau Fake Rate by using `combine` and `CombineHarvester`. 

## Installation

See [the README.md in the parent directory](../../../../#taufw).

## Usage 
The workflow is divided into three different steps: create inputs, create datacards and perform the fit.
DISCLAIMER: The Fitter code has been tested only for py2 and with `CMSSW_10_6_13`.

### Create inputs for combine.
To do so, please have a look at the `createinputs.py` script. 
In this script you can set up (by hardcoding) which samples set shall be used in the fit, their systematic variations, which observable/variables to be fitted in combine, the event selctions to be applied to the samples, the Pass and Fail regions and the eta regions where perform the fit.
Using the command options, you can choose era and channel.
Example:
```
python createinputs.py -c etau -y UL2017
```
will produce the inputs for UL2017 samples where event selections for etau channel are applied.
Inputs files and plots will be stored into the `input` folder.

### Create datacards
The `combine` tool need specific datacards as input. 
To create such datacards you can use the `writedatacard_DEV.py` script.
Before launching the command, you can set up (by hardcoding) which era, eta range, DeepTau WP and channel (or category) shall be used to create datacards.
No options are implemented for this script. 
Datacards will be stored into the `input/<era>/ETauFR/` folder as `.txt` files.

### Perform the fit (using combine)
Run `text2workspace.py` using the `zttmodels.py`:
```
text2workspace.py -P TauFW.Fitter.ETauFR.zttmodels:ztt_eff --PO "eff={pre-fit FR}" ./input/{era}/ETauFR/{wp}_eta{eta}.txt -o  ./input/{era}/ETauFR/WorkSpace{wp}_eta{eta}.root
```
Then you can run the fitting+plotting or the impacts or everything together, as in the following:
```
./doFitPlots.sh {wp} {eta} {era}
```
```
./doImpacts.sh {wp} {eta} {era}
``` 
```
./runAll.sh {wp} {era}
```
--OLD--
The script `Fitting_DEV.sh` provides several commands to perform the fit and create post-fit histograms for plotting.
Please, have a look at the script to understand where inputs and outputs of each command are stored.
Again, you should set up (by hardcoding) era, eta range and DeepTau WP for which the fit will be performed.
The `--PO` option of the `text2workspace.py` allows you to set up the POI of your fit (for ETau FR, the POI is the pre-fit fake rate).
You can set up other options for your purpose.

### Plotting code (deprecated)
To create pre and post fit plots of the observable you have used into the fit, you can use the `plotpostfit.py` script into the `Plotter` folder.
Please check the documentation in [plotting tool of parent directory](https://github.com/cms-tau-pog/TauFW/blob/master/Fitter/paper/plotpostfit.py) for instructions. To run please use the following command:
```
python plotpostfit.py -c et -y <ERA> -w <Working Point>
```
Plots are created automatically for both PASS and FAIL regions.

#### Plotting code (deprecated)
To create pre and post fit plots of the observable you have used into the fit, you can use the `PlotShapes.C` script into the `Plotter` folder.
To create a single plot, you should set up (by hardcoding) era, eta range, and DeepTau WP, and you can choose Fail or Pass and pre-fit or post-fit plots in `PlotShapes.C`.
```
root -b -q PlotShapes.C
```

You can use the script `FastPlot.sh` to create both pre and post fit plots in both Pass and Fail regions for all the WPs simultaneously.
You should only set up (by hardcoding) era and eta range in `FastPlot.sh`

Plots will be stored as `.pdf` and `.png` files in the `Pre-PostFitPlots` folder.
They will be also stored as `.root` files in the `rootFilesOfPlots` folder.

### Create ROOT and JSON SF files
The scripts `tau_ltf.py` is used to write the measured SFs in ROOT format, while `tau_createJSONs_VSe.py` is used to produce the associated JSON files.
The scripts need to be edited by hardcoding the values of the measured SFs for the appropriate eta region and era.
