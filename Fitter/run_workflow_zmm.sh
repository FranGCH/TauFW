# Define working points
J_VALUES=("VVLoose" "Loose" "Tight" "VTight" "VLoose" "Medium")
E_VALUES=("VVLoose" "Tight")

YEAR="2024"
CONFIG_TT="TauES_ID/config/config_coarse_PNet_TT.yml"
CONFIG_MM="TauES/config/FitSetup_mumu.yml"


SUFFIX="pnet" #suffix that the input files folder has

case "$SUFFIX" in
  pnet_hps)  TAGGER="PNet_HPS"  ;;
  pnet)      TAGGER="PNet"      ;;
  upart)     TAGGER="UParT"     ;;
  deeptau)   TAGGER="DeepTau"   ;;
  *) echo "unknown SUFFIX=$SUFFIX (use pnet_hps|pnet|upart|deeptau)"; exit 1 ;;
esac


INPUT_ROOT="input_pt_less_region_${SUFFIX}"
OUTPUT_ROOT="output_pt_less_region_${SUFFIX}"
PLOTS_ROOT="plots_pt_less_region_${SUFFIX}"


# Loop over each combination
for JET_WP in "${J_VALUES[@]}"; do
  for ELE_WP in "${E_VALUES[@]}"; do
    echo "=========================================="
    echo "=== Running Workflow for Jet=$JET_WP, Ele=$ELE_WP ==="
    echo "=========================================="

    BASE_INPUT="${INPUT_ROOT}/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
    BASE_OUTPUT="${OUTPUT_ROOT}/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
    BASE_PLOTS="${PLOTS_ROOT}/againstjet_${JET_WP}/againstelectron_${ELE_WP}"

    echo "=== Step 1: Running MultiDimFit ==="
    python3 TauES_ID/harvestDatacards_zmm.py -y $YEAR -c $CONFIG_MM -i ${BASE_INPUT}/ -o ${BASE_OUTPUT}/$YEAR/
    python3 TauES_ID/makecombinedfitTES_SF.py -y $YEAR -c $CONFIG_TT -i ${BASE_INPUT}/ \
      --input_file ${BASE_INPUT}/ztt_mt_tes_m_vis.inputs-$YEAR-13TeV_mutau.root \
      -o 3 --mumu_datacard_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt \
       --extratag ${TAGGER} \
      2>&1 | tee step1_multidimfit_${JET_WP}_${ELE_WP}_nocorrTES.log

    echo "=== Step 2: Running FitDiagnostics + PostFit ==="
    python3 TauES_ID/makecombinedfitTES_SF_postfit.py -y $YEAR -c $CONFIG_TT --indir ${BASE_OUTPUT}/ -o 3 \
      --mumu_input_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt \
      -cmm $CONFIG_MM --jet_wp ${JET_WP} --ele_wp ${ELE_WP} --extratag ${TAGGER} 2>&1 | tee step2_postfit_${JET_WP}_${ELE_WP}_nocorrTES.log

    echo "=== Step 3: Running Plots ==="
    python3 python/plot/runpostfit.py -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -y $YEAR --include-cr --tagger ${SUFFIX} \
      2>&1 | tee step3_plots_${JET_WP}_${ELE_WP}_nocorrTES.log
    python3 pre_post_plot_combiner.py --scan_dir ${BASE_PLOTS}/$YEAR/ --jet_wp ${JET_WP} --ele_wp ${ELE_WP} --tagger ${SUFFIX}
    python3 plot_measurements.py --jet_wp ${JET_WP} --ele_wp ${ELE_WP} --year $YEAR --tagger ${SUFFIX} --extratag ${TAGGER}

    echo "=== Step 4: Correction File Generation ==="
    python3 createroot_TES.py -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -f root -y $YEAR --tagger ${SUFFIX} --extragtag ${TAGGER}
    python3 createroot_TES.py -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -f json -y $YEAR --tagger ${SUFFIX} --extragtag ${TAGGER}

    echo "=== Workflow completed for Jet=$JET_WP, Ele=$ELE_WP ==="
    echo
  done
done

# Build 4x3 scan grid (rows=DM, cols=pt) per WP combo
echo "=== Building scan grids per WP combo ==="
python3 make_scan_grid.py --root ${PLOTS_ROOT} --dms DM0,DM1,DM2,DM10,DM11
python3 make_measurements_grid.py 

# Merge all JSON correction files
echo "=== Merging all JSON correction files ==="
python3 merge_tau_jsons.py --type both -o tau_sf_${SUFFIX}_uncorr/TauCorrections_$YEAR.json

echo "All workflows completed successfully!"

# #!/bin/bash
# Usage: ./run_workflow.sh
