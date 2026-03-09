# 2D TES/TauID Fit - Quick Start Guide

## Overview
This guide shows how to run the 2D TES/TauID fits using the automated shell scripts in the Fitter directory.

## Quick Start (3 steps)

### 1. Generate Input Files
```bash
cd ../  # Go to Fitter directory
bash run_inputs.sh
```
This creates input histograms for all working point combinations. Check logs for success.

### 2. Run Full Fit Workflow
```bash
bash run_workflow_zmm.sh
```
This runs the complete analysis pipeline for all jet/electron working point combinations:
- **Step 1:** MultiDimFit (scan over TES/ID parameters)
- **Step 2:** FitDiagnostics + PostFit pulls
- **Step 3:** Generate plots
- **Step 4:** Create correction files (ROOT + JSON)

### 3. Collect Results
Output is organized by working point:
```
output_pt_less_region/againstjet_[JET_WP]/againstelectron_[ELE_WP]/
plots_pt_less_region/againstjet_[JET_WP]/againstelectron_[ELE_WP]/
tau_sf/TauCorrections_2024.json  (merged result)
```

## Configuration Files
- **TES/TauID config:** `TauES_ID/config/config_coarse_TT.yml`
- **Mumu CR config:** `TauES/config/FitSetup_mumu.yml`
- **Year:** 2024 (edit scripts to change)

## Working Points
- **Jet WPs:** VTight, VVLoose, VLoose, Loose, Medium, Tight
- **Electron WPs:** VVLoose, Tight

Edit `run_inputs.sh` and `run_workflow_zmm.sh` to modify which combinations to process.

## Tips
- Logs are created with timestamps: `step*_[JET_WP]_[ELE_WP].log`
- For single working point: modify `J_VALUES` and `E_VALUES` arrays in scripts
- Final JSON corrections are merged in `tau_sf/TauCorrections_2024.json`
