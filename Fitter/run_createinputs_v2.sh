#!/usr/bin/env bash

# ===================================================================
# run_dm_job.sh: Executes createinputs for a single DM on a Condor node
# ===================================================================

# --- Receive arguments from the Condor JDL file ---
# $1: Absolute path to your CMSSW release
# $2: Absolute path to the EOS directory where you want to store the final ROOT file
# start =$(date +%s)
CMSSW_PATH=$1
shift
# EOS_OUTPUT_DIR=$2

echo "========================================="
echo "--- Condor Job Starting ---"
echo "Job running on host: $(hostname)"
echo "OS release: $(cat /etc/redhat-release)"
echo "-----------------------------------------"
echo "CMSSW Release Path: ${CMSSW_PATH}"
# echo "Final Output Directory on EOS: ${EOS_OUTPUT_DIR}"
echo "========================================="

# --- 1. Set up the CMSSW environment ---
# This is a crucial step that must be performed on the Condor worker node
echo ""
echo ">>> Setting up CMSSW environment..."
source /cvmfs/cms.cern.ch/cmsset_default.sh
# Change to your CMSSW release directory
cd ${CMSSW_PATH}
# Set up the CMSSW environment variables (equivalent to cmsenv)
eval `scramv1 runtime -sh`
# Return to the initial working directory of the job
cd -
echo "CMSSW_BASE is now: $CMSSW_BASE"
echo "Environment setup complete."

# --- 2. Run the Python script ---
# Change into the directory containing your script
cd ${CMSSW_PATH}/TauFW/Fitter

# Define a log file specific to this job
LOG_FILE="condor_job_output_mutauMedium.log"
echo ""
echo ">>> Running createinputsTES.py ..."
# echo "    Log file will be: ${LOG_FILE}"

# Execute your command, redirecting all output (stdout and stderr) to the log file
python3 TauES/createinputsTES.py \
    -y 2024 \
    -c TauES_ID/config/config_coarse_PNet_TT.yml \
    --outbase /eos/user/f/fcasalin/CMSSW_14_1_0_pre4/src/TauFW/Fitter/input_pt_less_region \
    "$@"
# python3 TauES/createinputsTES.py \
#     -y 2024 \
#     -c TauES/config/FitSetup_mumu.yml \
#     "$@"

# Record the exit status of the command
CMD_STATUS=$?
if [ ${CMD_STATUS} -ne 0 ]; then
    echo "!!! ERROR: Python script failed with exit code ${CMD_STATUS}."
fi

echo ""
echo "--- Condor Job Finished ---"
# echo time taken $(( $(date +%s) - start )) seconds
echo "========================================="