#!/bin/bash

JET_WP=$1
ELE_WP=$2

JET_WP_VALUES=("VVLoose" "VLoose" "Loose" "Medium" "Tight" "VTight")
ELE_WP_VALUES=("VVLoose" "Tight")

for J_WP in "${JET_WP_VALUES[@]}"; do
    if [ "$J_WP" == "$JET_WP" ]; then
    result_jet=1
    break
    fi
done

for E_WP in "${ELE_WP_VALUES[@]}"; do
    if [ "$E_WP" == "$ELE_WP" ]; then
    result_ele=1
    break
    fi
done
if [ "$result_jet" != 1 ] || [ "$result_ele" != 1 ]; then
    echo "Invalid working point(s)."
    echo "Valid JET WP: ${JET_WP_VALUES[*]}"
    echo "Valid ELE WP: ${ELE_WP_VALUES[*]}"
    exit 1
fi

rm -rf "output_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
rm -rf "postfit_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
rm -rf "plots_pt_less_region/againstjet_${JET_WP}/againstelectron_${ELE_WP}"
rm -rf "output_plots/jet_${JET_WP}_ele_${ELE_WP}"
rm -rf "FitDiagnosticsValues/VSjet${JET_WP}_VSele${ELE_WP}"
rm -rf "combined_pre_post/jet_${JET_WP}_ele_${ELE_WP}"
rm -rf "Measurements/VSjet${JET_WP}_VSele${ELE_WP}"