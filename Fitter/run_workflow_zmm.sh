#!/bin/bash
# Save as run_workflow.sh

echo "=== Step 1: Running MultiDimFit ==="
python3 TauES_ID/makecombinedfitTES_SF.py -y 2024 -c TauES_ID/config/config_coarse_TT.yml -i input_pt_less_region/againstjet_Medium/againstelectron_Tight/ --input_file input_pt_less_region/againstjet_Medium/againstelectron_Tight/ztt_mt_tes_m_vis.inputs-2024-13TeV_mutau.root -o 3 --mumu_datacard_file output_pt_less_region/againstjet_Medium/againstelectron_Tight/2024/ztt_mm_m_vis-baseline_mumu-2024-13TeV.txt 2>&1 | tee step1_multidimfit.log #-cmm TauES/config/FitSetup_mumu.yml

echo "=== Step 2: Running FitDiagnostics + PostFit ==="
python3 TauES_ID/makecombinedfitTES_SF_postfit.py -y 2024 -c TauES_ID/config/config_coarse_TT.yml --indir output_pt_less_region/againstjet_Medium/againstelectron_Tight/ -o 3 --mumu_input_file output_pt_less_region/againstjet_Medium/againstelectron_Tight/2024/ztt_mm_m_vis-baseline_mumu-2024-13TeV.txt -cmm TauES/config/FitSetup_mumu.yml 2>&1 | tee step2_postfit.log

echo "=== Step 3: Running Plots ==="
python3 python/plot/runpostfit.py -c TauES_ID/config/config_coarse_TT.yml --include-cr 2>&1 | tee step3_plots.log

echo "=== Workflow completed ==="

# python3 TauES_ID/harvestDatacards_zmm.py -y 2024 -c TauES/config/FitSetup_mumu.yml -i input_pt_less_region/againstjet_Medium/againstelectron_Tight/ -o output_pt_less_region/againstjet_Medium/againstelectron_Tight/2024/



# python3 TauES_ID/makecombinedfitTES_SF.py -y 2024 -c TauES_ID/config/config_coarse_TT.yml -i input_pt_less_region/againstjet_Medium/againstelectron_Tight/ --input_file input_pt_less_region/againstjet_Medium/againstelectron_Tight/ztt_mt_tes_m_vis.inputs-2024-13TeV_mutau.root -o 1 --mumu_datacard_file output_pt_less_region/againstjet_Medium/againstelectron_Tight/2024/ztt_mm_m_vis-baseline_mumu-2024-13TeV.txt && python3 TauES_ID/makecombinedfitTES_SF.py -y 2024 -c TauES_ID/config/config_coarse_TT.yml -i input_pt_less_region/againstjet_Medium/againstelectron_Tight/ --input_file input_pt_less_region/againstjet_Medium/againstelectron_Tight/ztt_mt_tes_m_vis.inputs-2024-13TeV_mutau.root -o 2 --mumu_datacard_file output_pt_less_region/againstjet_Medium/againstelectron_Tight/2024/ztt_mm_m_vis-baseline_mumu-2024-13TeV.txt
