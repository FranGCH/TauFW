# #!/bin/bash
# # Save as run_workflow.sh

# # J_VALUES=("VVLoose" "VLoose" "Loose" "Medium" "Tight")
# # E_VALUES=("Tight" "VVLoose")

# # Default values
# JET_WP="Tight"
# ELE_WP="VVLoose"
# YEAR="2024"
# CONFIG_TT="TauES_ID/config/config_coarse_TT.yml"
# CONFIG_MM="TauES/config/FitSetup_mumu.yml"
# # Parse arguments
# while [[ "$#" -gt 0 ]]; do
#     case $1 in
#         -j|--jet) JET_WP="$2"; shift ;;
#         -e|--electron) ELE_WP="$2"; shift ;;
#         *) echo "Unknown parameter passed: $1"; exit 1 ;;
#     esac
#     shift
# done

# echo "Using Jet WP: $JET_WP"
# echo "Using Electron WP: $ELE_WP"

# BASE_INPUT="input_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
# BASE_OUTPUT="output_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
# BASE_PLOTS="plots_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"

# echo "=== Step 1: Running MultiDimFit ==="
# #mumu datacard harvesting
# python3 TauES_ID/harvestDatacards_zmm.py -y $YEAR -c TauES/config/FitSetup_mumu.yml -i ${BASE_INPUT}/ -o ${BASE_OUTPUT}/$YEAR/ 
# #multi-dim fit
# python3 TauES_ID/makecombinedfitTES_SF.py -y $YEAR -c TauES_ID/config/config_coarse_TT.yml -i ${BASE_INPUT}/ --input_file ${BASE_INPUT}/ztt_mt_tes_m_vis.inputs-$YEAR-13TeV_mutau.root -o 3 --mumu_datacard_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt 2>&1 | tee step1_multidimfit.log #-cmm TauES/config/FitSetup_mumu.yml

# echo "=== Step 2: Running FitDiagnostics + PostFit ==="
# #post-fit
# python3 TauES_ID/makecombinedfitTES_SF_postfit.py -y $YEAR -c TauES_ID/config/config_coarse_TT.yml --indir ${BASE_OUTPUT}/ -o 3 --mumu_input_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt -cmm TauES/config/FitSetup_mumu.yml --jet_wp ${JET_WP} --ele_wp ${ELE_WP} 2>&1 | tee step2_postfit.log
# echo "=== Step 3: Running Plots ==="
# #post-fit plots
# python3 python/plot/runpostfit.py -c TauES_ID/config/config_coarse_TT.yml -j ${JET_WP} -e ${ELE_WP} --include-cr 2>&1 | tee step3_plots.log
# # Combine pre/post-fit plots with scans
# python3 pre_post_plot_combiner.py --scan_dir ${BASE_PLOTS}/$YEAR/ --jet_wp ${JET_WP} --ele_wp ${ELE_WP} 
# # Measurement summary plots
# python3 plot_measurements.py --jet_wp ${JET_WP} --ele_wp ${ELE_WP}
# #TODO: add measurement summary plot with ALL WPs together (Hardcore)

# echo "=== Step 4: Correction File Generation ==="
# python3 createroot_TES.py -c TauES_ID/config/config_coarse_TT.yml -j ${JET_WP} -e ${ELE_WP} -f root 
# python3 createroot_TES.py -c TauES_ID/config/config_coarse_TT.yml -j ${JET_WP} -e ${ELE_WP} -f json

# # python3 merge_tau_jsons.py 

# echo "=== Workflow completed ==="


#!/bin/bash
# Save as run_workflow.sh
# Usage: ./run_workflow.sh

# Define working points
J_VALUES=("VVLoose" "VLoose" "Loose" "Medium" "Tight") #add VTight 
E_VALUES=("VVLoose"  "Tight")

YEAR="2024"
CONFIG_TT="TauES_ID/config/config_coarse_TT.yml"
CONFIG_MM="TauES/config/FitSetup_mumu.yml"

# Loop over each combination
for JET_WP in "${J_VALUES[@]}"; do
  for ELE_WP in "${E_VALUES[@]}"; do
    echo "=========================================="
    echo "=== Running Workflow for Jet=$JET_WP, Ele=$ELE_WP ==="
    echo "=========================================="

    BASE_INPUT="input_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
    BASE_OUTPUT="output_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
    BASE_PLOTS="plots_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"

    echo "=== Step 1: Running MultiDimFit ==="
    python3 TauES_ID/harvestDatacards_zmm.py -y $YEAR -c $CONFIG_MM -i ${BASE_INPUT}/ -o ${BASE_OUTPUT}/$YEAR/
    python3 TauES_ID/makecombinedfitTES_SF.py -y $YEAR -c $CONFIG_TT -i ${BASE_INPUT}/ \
      --input_file ${BASE_INPUT}/ztt_mt_tes_m_vis.inputs-$YEAR-13TeV_mutau.root \
      -o 3 --mumu_datacard_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt \
      2>&1 | tee step1_multidimfit_${JET_WP}_${ELE_WP}.log

    echo "=== Step 2: Running FitDiagnostics + PostFit ==="
    python3 TauES_ID/makecombinedfitTES_SF_postfit.py -y $YEAR -c $CONFIG_TT --indir ${BASE_OUTPUT}/ -o 3 \
      --mumu_input_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt \
      -cmm $CONFIG_MM --jet_wp ${JET_WP} --ele_wp ${ELE_WP} 2>&1 | tee step2_postfit_${JET_WP}_${ELE_WP}.log

    echo "=== Step 3: Running Plots ==="
    python3 python/plot/runpostfit.py -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} --include-cr \
      2>&1 | tee step3_plots_${JET_WP}_${ELE_WP}.log
    python3 pre_post_plot_combiner.py --scan_dir ${BASE_PLOTS}/$YEAR/ --jet_wp ${JET_WP} --ele_wp ${ELE_WP}
    python3 plot_measurements.py --jet_wp ${JET_WP} --ele_wp ${ELE_WP}

    echo "=== Step 4: Correction File Generation ==="
    python3 createroot_TES.py -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -f root
    python3 createroot_TES.py -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -f json

    echo "=== Workflow completed for Jet=$JET_WP, Ele=$ELE_WP ==="
    echo
  done
done
# Merge all JSON correction files
echo "=== Merging all JSON correction files ==="
python3 merge_tau_jsons.py --type both -o tau_sf/TauCorrections_2024.json

echo "All workflows completed successfully!"
