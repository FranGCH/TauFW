# WP's to use
JET_VALUES={"VVLoose","VLoose","Loose","Medium","Tight"}
ELE_VALUES={"VVLoose","Tight"}

YEAR="2024"
CONFIT_TT="TauES_ID/config/config_coarse_TT.yml"
CONFIG_MM="TauES/config/FitSetup_mumu.yml"

# loop over each combination

for J_WP IN "${JET_VALUES[@]}";do
    for E_WP IN "${ELE_VALUES[@]}";do
        echo "==========="
        echo "Running for Jet WP: ${J_WP}, Ele WP: ${E_WP}"
        echo "==========="

        BASE_INPUT="input_pt_less_region/againstjet_${J_WP}/againstelectron_${E_WP}"
        BASE_OUTPUT="output_pt_less_region/againstjet_${J_WP}/againstelectron_${E_WP}"
        BASE_PLOTS="plots_pt_less_region/againstjet_${J_WP}/againstelectron_${E_WP}"    

        echo "Step 1: Running MutiDimFit"
        python3 TauES_ID/makecombinedfitTES_SF.py -y $YEAR -c $CONFIT_TT -i ${BASE_INPUT}/\
            --input_file ${BASE_INPUT}/ztt_mt_tes_m_vis.inputs-$YEAR-13TeV_mutau.root \
            -o 3\
            2>&1 | tee step1_multidimfit_jet${J_WP}_ele${E_WP}_nomm.log

        echo "Step 2: Running FitDiagnostics + PostFit"
        python3 TauES_ID/makecombinedfitTES_SF_postfit.py -y $YEAR -c $CONFIG_TT --indir ${BASE_OUTPUT}/ -o 3 \
             --jet_wp ${J_WP} --ele_wp ${E_WP} 2>&1 | tee step2_postfit_${J_WP}_${E_WP}_nomm.log

        echo "Step 3: Running Plots"
        python3 python/plot/runpostfit.py -c $CONFIG_TT -j ${J_WP} -e ${E_WP} \
            2>&1 | tee step3_plots_${J_WP}_${E_WP}_nomm.log

        echo "Worflow for Jet WP: ${J_WP}, Ele WP: ${E_WP} completed."
        echo 
    done
done
