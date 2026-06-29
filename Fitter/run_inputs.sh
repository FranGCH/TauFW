#!/usr/bin/env bash

J_VALUES=("Tight")
E_VALUES=("VVLoose")

YEAR=2024
CONFIG="TauES_ID/config/config_coarse_PNet_TT.yml"
CONFIG_MM="TauES/config/FitSetup_mumu.yml" #CR config

# Create timestamped log file
TIMESTAMP=$(date +"%Y%m%d_%H%M")
LOGFILE="run_TES_${TIMESTAMP}.log"

echo "Logging to ${LOGFILE}"
echo "Run started at $(date)" | tee -a "$LOGFILE"

for j in "${J_VALUES[@]}"; do
    for e in "${E_VALUES[@]}"; do
        if [[ "$j" == "Medium" && "$e" == "Tight" ]]; then
            echo "Skipping combination j=${j}, e=${e} (already done)" | tee -a "$LOGFILE"
            continue
        fi
        
        {
            echo "------------------------------------------------------------"
            echo "Running: -j ${j}, -e ${e}"
            echo "------------------------------------------------------------"
        } | tee -a "$LOGFILE"

        python3 TauES/createinputsTES.py \
            -y ${YEAR} \
            -c ${CONFIG} \
            --outbase input_pt_less_region_pnet \
            -j ${j} \
            -e ${e} 2>&1 | tee -a "$LOGFILE"
        python3 TauES/createinputsTES.py \
            -y ${YEAR} \
            -c ${CONFIG_MM} \
            --outbase input_pt_less_region_pnet \
            -j ${j} \
            -e ${e} 

        STATUS=${PIPESTATUS[0]}   # Correct status when using tee

        if [[ $STATUS -ne 0 ]]; then
            echo "❌ ERROR: iteration failed for j=${j}, e=${e} (exit code $STATUS)" | tee -a "$LOGFILE"
            echo "→ Skipping and continuing..." | tee -a "$LOGFILE"
            continue
        fi

        echo "✔ Completed j=${j}, e=${e}" | tee -a "$LOGFILE"
        echo | tee -a "$LOGFILE"

    done
done

echo "Run finished at $(date)" | tee -a "$LOGFILE"
