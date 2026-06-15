# WP's to use
JET_VALUES=("VLoose")
#("VVLoose" "VLoose" "Loose" "Medium" "Tight" "VTight")
# JET_VALUES=("Loose")
ELE_VALUES=("VVLoose")
# "Tight")

YEAR="2024"
CONFIG_TT="/eos/user/f/fcasalin/CMSSW_14_1_0_pre4/src/TauFW/Fitter/TauES_ID/config/config_coarse_TT.yml"
CONFIG_MM="/eos/user/f/fcasalin/CMSSW_14_1_0_pre4/src/TauFW/Fitter/TauES/config/FitSetup_mumu.yml"

# loop over each combination

for J_WP in "${JET_VALUES[@]}";do
    for E_WP in "${ELE_VALUES[@]}";do
        # if [ "$J_WP" == "Medium" ] && [ "$E_WP" == "Tight" ]; then
        #     echo "Skipping combination: Jet WP = ${J_WP}, Ele WP = ${E_WP}."
        #     continue
        # elif [ "$J_WP" == "Tight" ] && [ "$E_WP" == "VVLoose" ]; then
        #     echo "Skipping combination: Jet WP = ${J_WP}, Ele WP = ${E_WP}."
        #     continue
        # if [ "$J_WP" == "Medium" ] && [ "$E_WP" == "Tight" ]; then
        #     echo "Skipping combination: Jet WP = ${J_WP}, Ele WP = ${E_WP}."
        #     continue
        # fi
        echo "==========="
        echo "Running for Jet WP: ${J_WP}, Ele WP: ${E_WP}"
        echo "==========="

        BASE_INPUT="/eos/user/f/fcasalin/CMSSW_14_1_0_pre4/src/TauFW/Fitter/input_pt_less_region/againstjet_${J_WP}/againstelectron_${E_WP}"
        BASE_OUTPUT="/eos/user/f/fcasalin/CMSSW_14_1_0_pre4/src/TauFW/Fitter/output_pt_less_region/againstjet_${J_WP}/againstelectron_${E_WP}"
        BASE_PLOTS="plots_pt_less_region/againstjet_${J_WP}/againstelectron_${E_WP}"    

        echo "Step 1: Running MutiDimFit"
        python3 TauES_ID/harvestDatacards_zmm.py -y $YEAR -c $CONFIG_MM -i ${BASE_INPUT} -o ${BASE_OUTPUT}/$YEAR
        python3 TauES_ID/makecombinedfitTES_SF.py -y $YEAR -c $CONFIG_TT -i ${BASE_INPUT}\
            --input_file ${BASE_INPUT}/ztt_mt_tes_m_vis.inputs-$YEAR-13TeV_mutau.root \
            -o 3 --mumu_datacard_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt \
            2>&1 | tee WorkFlow_LOGs/step1_multidimfit_jet${J_WP}_ele${E_WP}_v3.log

        echo "Step 2: Running FitDiagnostics + PostFit"

        python3 TauES_ID/makecombinedfitTES_SF_postfit.py -y $YEAR -c $CONFIG_TT --indir ${BASE_OUTPUT}/ -o 3 \
            --mumu_input_file ${BASE_OUTPUT}/$YEAR/ztt_mm_m_vis-baseline_mumu-$YEAR-13TeV.txt \
            -cmm $CONFIG_MM --jet_wp ${J_WP} --ele_wp ${E_WP} 2>&1 | tee WorkFlow_LOGs/step2_postfit_${J_WP}_${E_WP}_v3.log

        echo "Step 3: Running Plots"
        python3 python/plot/runpostfit.py -c $CONFIG_TT -j ${J_WP} -e ${E_WP} --include-cr \
            2>&1 | tee WorkFlow_LOGs/step3_plots_${J_WP}_${E_WP}_v3.log
        python3 pre_post_plot_combiner.py --scan_dir ${BASE_PLOTS}/$YEAR/ --jet_wp ${J_WP} --ele_wp ${E_WP}
        python3 plot_measurements.py --jet_wp ${J_WP} --ele_wp ${E_WP} 

        # echo "Step 4: Correction File Generation"
        # python3 createroot_TES.py -c $CONFIG_TT -j ${J_WP} -e ${E_WP} -f root
        # python3 createroot_TES.py -c $CONFIG_TT -j ${J_WP} -e ${E_WP} -f json

        # echo "Worflow for Jet WP: ${J_WP}, Ele WP: ${E_WP} completed."
        # echo 
    done
done

# echo "Merge all JSON correction files"
# python3 merge_tau_jsons.py --type both -o tau_sf/TauCorrections_2024.json