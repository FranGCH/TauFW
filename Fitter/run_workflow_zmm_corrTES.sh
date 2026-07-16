# Workflow for the CORRELATED-TES fit variant.
# One TES POI per DM (correlated across pT bins) + 3 TauID SF POIs per DM (uncorrelated by pT).
# Sibling of run_workflow_zmm.sh (uncorrelated TES). Uses parallel _corrTES scripts and
# a separate output/plot tree so the two fits can be compared side-by-side.

J_VALUES=("VVLoose")
#("VVLoose" "VLoose" "Loose" "Medium" "Tight" "VTight") 
E_VALUES=("Tight") 

YEAR="2024"
CONFIG_TT="TauES_ID/config/config_coarse_PNet_TT_hps.yml"
CONFIG_MM="TauES/config/FitSetup_mumu.yml"

SUFFIX="pnet_hps" #suffix that the input files folder has

case "$SUFFIX" in
  pnet_hps)  TAGGER="PNet_HPS"  ;;
  pnet)      TAGGER="PNet"      ;;
  upart)     TAGGER="UParT"     ;;
  deeptau)   TAGGER="DeepTau"   ;;
  *) echo "unknown SUFFIX=$SUFFIX (use pnet_hps|pnet|upart|deeptau)"; exit 1 ;;
esac

# Separate I/O tree for the corrTES variant — keep uncorrelated outputs intact
INPUT_ROOT="input_pt_less_region_${SUFFIX}"           # inputs are shared with the uncorrelated workflow
OUTPUT_ROOT="output_pt_less_region_corrTES_${SUFFIX}"  # corrTES outputs go here
PLOTS_ROOT="plots_pt_less_region_corrTES_${SUFFIX}"    # corrTES plots go here

for JET_WP in "${J_VALUES[@]}"; do
  for ELE_WP in "${E_VALUES[@]}"; do
    echo "=========================================="
    echo "=== [corrTES] Jet=$JET_WP, Ele=$ELE_WP ==="
    echo "=========================================="

    BASE_INPUT="${INPUT_ROOT}/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
    BASE_OUTPUT="${OUTPUT_ROOT}/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
    BASE_POSTFIT="${OUTPUT_ROOT/output/postfit}/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
    BASE_PLOTS="${PLOTS_ROOT}/againstjet_${JET_WP}/againstelectron_${ELE_WP}"

    echo "=== Step 1: Zmm CR datacard + per-DM MultiDimFit (corrTES) ==="
    python3 TauES_ID/harvestDatacards_zmm.py -y $YEAR -c $CONFIG_MM -i ${BASE_INPUT}/ -o ${BASE_OUTPUT}/$YEAR/
    python3 TauES_ID/makecombinedfitTES_SF_corrTES.py -y $YEAR -c $CONFIG_TT -i ${BASE_INPUT}/ \
      --input_file ${BASE_INPUT}/ztt_mt_tes_m_vis.inputs-$YEAR-13TeV_mutau.root \
      -o 3 --mumu_datacard_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt \
      --extratag ${TAGGER} \
      2>&1 | tee step1_multidimfit_corrTES_${JET_WP}_${ELE_WP}_${SUFFIX}_v2.log

    # echo "=== Step 2: FitDiagnostics + PostFit (corrTES, per-DM) ==="
    # python3 TauES_ID/makecombinedfitTES_SF_postfit_corrTES.py -y $YEAR -c $CONFIG_TT \
    #   --indir ${BASE_OUTPUT}/ -o 3 \
    #   --mumu_input_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt \
    #   -cmm $CONFIG_MM --jet_wp ${JET_WP} --ele_wp ${ELE_WP} --extratag ${TAGGER} \
    #   2>&1 | tee step2_postfit_corrTES_${JET_WP}_${ELE_WP}_${SUFFIX}.log

    # echo "=== Step 2b: Nuisance pull plots (per DM) ==="
    # PULL_TOOL="${CMSSW_BASE}/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py"
    # PULLDIR="${BASE_PLOTS}/${YEAR}/pulls"
    # mkdir -p "${PULLDIR}"

    # for DM in DM0 DM1 DM2 DM10 DM11; do
    #   FD="${BASE_POSTFIT}/${YEAR}/fitDiagnostics.mt_m_vis-${DM}_mutau_${TAGGER}-${YEAR}-13TeV.root"
    #   if [ ! -f "${FD}" ]; then echo "  [pulls] missing ${FD} — skip ${DM}"; continue; fi
    #   PULLTXT="${PULLDIR}/pulls_${JET_WP}_${ELE_WP}_${DM}.txt"
    #   python3 "${PULL_TOOL}" --poi tes_${DM} --vtol=0.1 "${FD}" 2>/dev/null | sed 's/[!,]/ /g' | tail -n +4 > "${PULLTXT}"
    #   python3 scripts/plot_pulls.py -f "${PULLTXT}" -o "${PULLDIR}/pulls_${JET_WP}_${ELE_WP}_${DM}" -t "${JET_WP}/${ELE_WP} ${DM}"
    # done

    # if [[ $SUFFIX ]] ; then
    #   echo "=== Step 3: Plots ==="
    #   python3 python/plot/runpostfit.py --variant corr -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -y $YEAR --include-cr --tagger ${SUFFIX} \
    #     2>&1 | tee step3_plots_corrTES_${JET_WP}_${ELE_WP}_${SUFFIX}.log
    #   python3 pre_post_plot_combiner.py --variant corr --scan_dir ${BASE_PLOTS}/$YEAR/ --jet_wp ${JET_WP} --ele_wp ${ELE_WP} --tagger ${SUFFIX}
    #   python3 plot_measurements.py --variant corr --jet_wp ${JET_WP} --ele_wp ${ELE_WP} --year $YEAR --tagger ${SUFFIX} --extratag ${TAGGER}

    #   echo "=== Step 4: Correction file generation (corrTES) ==="
    #   python3 createroot_TES.py --variant corr -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -f root -y $YEAR --tagger ${SUFFIX} --extratag ${TAGGER}
    #   python3 createroot_TES.py --variant corr -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -f json -y $YEAR --tagger ${SUFFIX} --extratag ${TAGGER}

    # else
    #   echo "=== Step 3: Plots ==="
    #   python3 python/plot/runpostfit.py --variant corr -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -y $YEAR --include-cr \
    #     2>&1 | tee step3_plots_corrTES_${JET_WP}_${ELE_WP}_${SUFFIX}.log
    #   python3 pre_post_plot_combiner.py --variant corr --scan_dir ${BASE_PLOTS}/$YEAR/ --jet_wp ${JET_WP} --ele_wp ${ELE_WP}
    #   python3 plot_measurements.py --variant corr --jet_wp ${JET_WP} --ele_wp ${ELE_WP} --year $YEAR --extratag ${TAGGER}

    #   echo "=== Step 4: Correction file generation (corrTES) ==="
    #   python3 createroot_TES.py --variant corr -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -f root -y $YEAR --extratag ${TAGGER}
    #   python3 createroot_TES.py --variant corr -c $CONFIG_TT -j ${JET_WP} -e ${ELE_WP} -f json -y $YEAR --extratag ${TAGGER}

    # fi

    echo "=== [corrTES] done: Jet=$JET_WP, Ele=$ELE_WP ==="
    echo
  done
done

# echo "=== Scan grids per WP combo ==="
# python3 make_scan_grid.py --variant corr --root ${PLOTS_ROOT} --dms DM0,DM1,DM2,DM10,DM11
# python3 make_measurements_grid.py

# echo "=== Merging JSON correction files (corrTES) ==="
# if [[ $SUFFIX ]] ; then
#   python3 merge_tau_jsons.py --variant corr --type both -c $CONFIG_TT -o tau_sf/TauCorrections_${YEAR}_corrTES_${SUFFIX}.json 2>&1 | tee step_5_merge_json_${SUFFIX}.log
# else
#   python3 merge_tau_jsons.py --variant corr --type both -c $CONFIG_TT -o tau_sf/TauCorrections_${YEAR}_corrTES.json
# fi
# # python3 merge_tau_jsons.py --variant corr --type both -c $CONFIG_TT -o tau_sf/TauCorrections_${YEAR}_corrTES.json

# echo "=== 1D profile NLLs from per-DM joint fits ==="
# python3 plot1D_NLL_profiles.py --variant corr --year $YEAR --indir ${OUTPUT_ROOT} --outdir ${PLOTS_ROOT} --dms DM0,DM1,DM2,DM10,DM11 --extratag ${TAGGER} 2>&1 | tee step_5_plot1D_NLL_profiles_${SUFFIX}.log

# echo "=== Building combined TauEnergy_SF + TauID_SF (corrTES variant) ==="
# python3 make_tid_2025.py --variant corr -y $YEAR -c $CONFIG_TT --extratag ${TAGGER}

echo "All corrTES workflows completed!"