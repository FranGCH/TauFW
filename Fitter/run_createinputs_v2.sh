#!/bin/bash

# ===================================================================
# run_dm_job.sh: Executes createinputs for a single DM on a Condor node
# ===================================================================

# --- Receive arguments from the Condor JDL file ---
# $1: Absolute path to your CMSSW release
# $2: Absolute path to the EOS directory where you want to store the final ROOT file
# start =$(date +%s)
CMSSW_PATH=$1
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
echo "    Log file will be: ${LOG_FILE}"

# Execute your command, redirecting all output (stdout and stderr) to the log file
# python3 TauES/createinputsTES.py \
#     -y 2024 \
#     -c TauES_ID/config/config_coarse_TT.yml \
#     -j Tight \
#     -e VVLoose > ${LOG_FILE} 2>&1

python3 TauES/createinputsTES.py \
    -y 2024 \
    -c TauES/config/FitSetup_mumu.yml \
    -j Medium \
    -e VVLoose > ${LOG_FILE} 2>&1

# Record the exit status of the command
CMD_STATUS=$?
if [ ${CMD_STATUS} -ne 0 ]; then
    echo "!!! ERROR: Python script failed with exit code ${CMD_STATUS}."
fi

# --- 3. Process and copy the output file ---
# echo ""
# echo ">>> Processing output..."

# # Find the output ROOT filename from the log file
# # This grep command is customized based on your previous log output and is very reliable
# OUTPUT_BASENAME=$(grep "outputfile Name:" ${LOG_FILE} | head -n 1 | awk -F'/' '{print $NF}' | awk '{print $1}')
# LOCAL_OUTPUT_PATH="input_dzytest_upart_v8_inclu_3pt/againstjet_Medium/againstelectron_VVLoose/rebinning/${OUTPUT_BASENAME}"

# # Check if the file was actually created
# if [ -f "${LOCAL_OUTPUT_PATH}" ]; then
#     echo "    Output file found: ${LOCAL_OUTPUT_PATH}"
#     echo "    Copying to EOS directory: ${EOS_OUTPUT_DIR}"
    
#     # Ensure the target directory on EOS exists
#     # mkdir -p ${EOS_OUTPUT_DIR}
    
#     # Copy the file
#     cp ${LOCAL_OUTPUT_PATH} ${EOS_OUTPUT_DIR}/
    
#     # Check if the copy was successful
#     if [ $? -eq 0 ]; then
#         echo "    Successfully copied to ${EOS_OUTPUT_DIR}/${OUTPUT_BASENAME}"
#     else
#         echo "!!! ERROR: Failed to copy output file to EOS!"
#     fi
# else
#     echo "!!! CRITICAL ERROR: Output ROOT file was not found at ${LOCAL_OUTPUT_PATH}!"
#     echo "!!! This likely confirms that running with -d alone produces no output."
#     echo "--- Displaying last 20 lines of log file for debugging ---"
#     tail -n 20 ${LOG_FILE}
# fi

echo ""
echo "--- Condor Job Finished ---"
# echo time taken $(( $(date +%s) - start )) seconds
echo "========================================="