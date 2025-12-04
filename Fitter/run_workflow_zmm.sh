#!/bin/bash
# Save as run_workflow.sh

# Default values
JET_WP="Medium"
ELE_WP="VVLoose"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -j|--jet) JET_WP="$2"; shift ;;
        -e|--electron) ELE_WP="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

echo "Using Jet WP: $JET_WP"
echo "Using Electron WP: $ELE_WP"

BASE_INPUT="input_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
BASE_OUTPUT="output_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
BASE_PLOTS="plots_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"

echo "=== Step 1: Running MultiDimFit ==="
python3 TauES_ID/harvestDatacards_zmm.py -y 2024 -c TauES/config/FitSetup_mumu.yml -i ${BASE_INPUT}/ -o ${BASE_OUTPUT}/2024/


python3 TauES_ID/makecombinedfitTES_SF.py -y 2024 -c TauES_ID/config/config_coarse_TT.yml -i ${BASE_INPUT}/ --input_file ${BASE_INPUT}/ztt_mt_tes_m_vis.inputs-2024-13TeV_mutau.root -o 3 --mumu_datacard_file ${BASE_OUTPUT}/2024/ztt_mm_m_vis-baseline_mumu-2024-13TeV.txt 2>&1 | tee step1_multidimfit.log #-cmm TauES/config/FitSetup_mumu.yml

echo "=== Step 2: Running FitDiagnostics + PostFit ==="
python3 TauES_ID/makecombinedfitTES_SF_postfit.py -y 2024 -c TauES_ID/config/config_coarse_TT.yml --indir ${BASE_OUTPUT}/ -o 3 --mumu_input_file ${BASE_OUTPUT}/2024/ztt_mm_m_vis-baseline_mumu-2024-13TeV.txt -cmm TauES/config/FitSetup_mumu.yml --jet_wp ${JET_WP} --ele_wp ${ELE_WP} 2>&1 | tee step2_postfit.log

echo "=== Step 3: Running Plots ==="
python3 python/plot/runpostfit.py -c TauES_ID/config/config_coarse_TT.yml -j ${JET_WP} -e ${ELE_WP} --include-cr 2>&1 | tee step3_plots.log
python3 pre_post_plot_combiner.py --scan_dir ${BASE_PLOTS}/2024/ --jet_wp ${JET_WP} --ele_wp ${ELE_WP} 
python3 plot_measurements.py --jet_wp ${JET_WP} --ele_wp ${ELE_WP}
echo "=== Workflow completed ==="



## MuMu datacard harvesting (if needed)
# python3 TauES_ID/harvestDatacards_zmm.py -y 2024 -c TauES/config/FitSetup_mumu.yml -i input_pt_less_region/againstjet_Medium/againstelectron_VVLoose/ -o output_pt_less_region/againstjet_Medium/againstelectron_VVLoose/2024/



# python3 TauES_ID/makecombinedfitTES_SF.py -y 2024 -c TauES_ID/config/config_coarse_TT.yml -i input_pt_less_region/againstjet_Medium/againstelectron_VVLoose/ --input_file input_pt_less_region/againstjet_Medium/againstelectron_VVLoose/ztt_mt_tes_m_vis.inputs-2024-13TeV_mutau.root -o 1 --mumu_datacard_file output_pt_less_region/againstjet_Medium/againstelectron_VVLoose/2024/ztt_mm_m_vis-baseline_mumu-2024-13TeV.txt && python3 TauES_ID/makecombinedfitTES_SF.py -y 2024 -c TauES_ID/config/config_coarse_TT.yml -i input_pt_less_region/againstjet_Medium/againstelectron_VVLoose/ --input_file input_pt_less_region/againstjet_Medium/againstelectron_VVLoose/ztt_mt_tes_m_vis.inputs-2024-13TeV_mutau.root -o 2 --mumu_datacard_file output_pt_less_region/againstjet_Medium/againstelectron_VVLoose/2024/ztt_mm_m_vis-baseline_mumu-2024-13TeV.txt
