#!/usr/bin/env python3
"""
Script to create combined tau correction JSON files with variation parameters and multiple Working Points.
Author: haawedik, NCBJ. Nov. 2025
Usage: python3 merge_tau_jsons.py --type both -o tau_sf/TauCorrections_2024.json

This script merges individual per-WP JSON files into a single correctionlib file.
Input files should be named like: TauES_SF_dm_DeepTau2018v2p5_2024_VSjetMedium_VSeleVVLoose.json
"""

import json
import glob
import os
import re
from argparse import ArgumentParser
import correctionlib.schemav2 as cs


def extract_genmatch_data(data_node):
    """
    Navigate through the nested category structure to extract the genmatch->DM->syst->pT data.
    The input files have structure: wp_VSmu -> wp_VSe -> wp_VSjet -> genmatch -> ...
    We want to extract just the genmatch part onwards.
    """
    # Navigate: wp_VSmu (Tight) -> wp_VSe (e.g. VVLoose) -> wp_VSjet (e.g. Medium) -> genmatch data
    try:
        # Get the first (and only) wp_VSmu item
        vsmu_content = data_node.content[0].value
        # Get the first (and only) wp_VSe item  
        vse_content = vsmu_content.content[0].value
        # Get the first (and only) wp_VSjet item
        vsjet_content = vse_content.content[0].value
        # This should be the genmatch category
        return vsjet_content
    except (KeyError, IndexError, TypeError, AttributeError) as e:
        print(f"  Warning: Could not extract genmatch data: {e}")
        return None


def create_combined_correction(input_dir, output_filename, correction_type="tes"):
    """Create a single correction with variation parameter from individual JSON files."""
    
    print(f"Creating {correction_type.upper()} correction file...")
    
    # Create the correction
    combined_corr = create_single_correction(input_dir, correction_type)
    
    if not combined_corr:
        print(f"ERROR: Could not create {correction_type.upper()} correction!")
        return
    
    # Create correction set
    cset = cs.CorrectionSet(
        schema_version=2,
        description=combined_corr.description,
        corrections=[combined_corr],
    )
    
    # Write file
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    with open(output_filename, "w") as fout:
        print(f"\n>>> Writing {correction_type.upper()} correction to {output_filename}!")
        fout.write(cset.json(exclude_unset=True))
    
    print(f">>> Successfully created {correction_type.upper()} correction file")
    
    # Test the correction
    test_correction(output_filename, combined_corr)


def test_correction(filename, combined_corr):
    """Test the created correction file."""
    print(f"\n>>> Testing correction...")
    try:
        import correctionlib
        cset = correctionlib.CorrectionSet.from_file(filename)
        evaluator = cset[combined_corr.name]
        
        # Inspect inputs to build test arguments
        input_names = [inp.name for inp in combined_corr.inputs]
        print(f"  Inputs: {input_names}")
        
        # Define test values
        test_values = {
            "genmatch": 5,
            "DM": 0,
            "pT": 50.0,
            "wp_VSjet": "Medium",
            "wp_VSe": "VVLoose",
            "wp_VSmu": "Tight",
            "syst": "nom",
        }
        
        # Build argument list in correct order
        args = []
        syst_idx = -1
        for idx, name in enumerate(input_names):
            if name in test_values:
                args.append(test_values[name])
            else:
                print(f"  Warning: Unknown input '{name}', using default")
                args.append(0)
            
            if name == "syst":
                syst_idx = idx

        # Run test
        if syst_idx >= 0:
            for var in ['nom', 'up', 'down']:
                current_args = list(args)
                current_args[syst_idx] = var
                try:
                    result = evaluator.evaluate(*current_args)
                    print(f"  Test {var}: {result:.6f}")
                except Exception as e:
                    print(f"  Test {var} failed: {e}")
        else:
            result = evaluator.evaluate(*args)
            print(f"  Test nominal: {result:.6f}")
            
        print(f">>> Correction validation successful!")
        
    except Exception as e:
        print(f">>> ERROR testing correction: {e}")
        import traceback
        traceback.print_exc()


def create_combined_both_corrections(input_dir, output_filename):
    """Create a single file with both TES and TauIdSF corrections."""
    
    print("Creating combined file with both TES and TauIdSF corrections...")
    
    # Create TES correction
    print("\n=== Creating TES correction ===")
    tes_corr = create_single_correction(input_dir, "tes")
    
    # Create TauIdSF correction  
    print("\n=== Creating TauIdSF correction ===")
    id_corr = create_single_correction(input_dir, "id")
    
    corrections_list = []
    if tes_corr:
        corrections_list.append(tes_corr)
    if id_corr:
        corrections_list.append(id_corr)
        
    if not corrections_list:
        print("ERROR: Could not create any corrections!")
        return
    
    # Create combined correction set with both corrections
    combined_cset = cs.CorrectionSet(
        schema_version=2,
        description="Tau Energy Scale and ID Scale Factor corrections with systematic variations",
        corrections=corrections_list,
    )
    
    # Write file
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    with open(output_filename, "w") as fout:
        print(f"\n>>> Writing combined corrections to {output_filename}!")
        fout.write(combined_cset.json(exclude_unset=True))
    
    print(f">>> Successfully created combined file with {len(corrections_list)} correction(s)")
    
    # Test each correction
    for corr in corrections_list:
        test_correction(output_filename, corr)


def create_single_correction(input_dir, correction_type):
    """Helper function to create a single correction with all WP combinations merged."""
    
    # Settings based on type
    if correction_type == "tes":
        pattern_base = "TauES"
        final_name = "TauES_2024"
        description = "Tau Energy Scale corrections with all WP combinations"
    else:
        pattern_base = "TauID"
        final_name = "TauIdSF_2024"
        description = "Tau ID Scale Factor corrections with all WP combinations"
    
    # Regex to find files and extract WPs from filename
    # Matches: TauES_SF_dm_DeepTau2018v2p5_2024_VSjetMedium_VSeleVVLoose.json
    wp_regex = re.compile(r".*VSjet(?P<jet>[a-zA-Z]+)_VSele(?P<ele>[a-zA-Z]+)\.json$")
    
    # Find all potential files
    search_pattern = os.path.join(input_dir, f"*{pattern_base}*.json")
    all_files = glob.glob(search_pattern)
    all_files.sort()
    
    if not all_files:
        print(f"No JSON files found for {correction_type} with pattern {search_pattern}")
        return None

    # Structure to hold data: wp_map[jet_wp][ele_wp] = genmatch_data
    wp_map = {}
    reference_inputs = None
    
    print(f"Scanning {len(all_files)} files for {correction_type}...")
    
    for fpath in all_files:
        fname = os.path.basename(fpath)
        
        match = wp_regex.match(fname)
        if match:
            jet_wp = match.group('jet')
            ele_wp = match.group('ele')
            
            print(f"  Found WP: VSjet={jet_wp}, VSele={ele_wp} in {fname}")
            
            try:
                with open(fpath, 'r') as f:
                    cset = cs.CorrectionSet.parse_obj(json.load(f))
                
                # Find the main correction
                found_corr = None
                for corr in cset.corrections:
                    if '_up' not in corr.name and '_down' not in corr.name:
                        found_corr = corr
                        break
                
                if found_corr:
                    # Store inputs from the first valid file (they should all be the same)
                    if reference_inputs is None:
                        reference_inputs = found_corr.inputs
                    
                    # Extract the innermost data (genmatch onwards)
                    genmatch_data = extract_genmatch_data(found_corr.data)
                    
                    if genmatch_data:
                        if jet_wp not in wp_map:
                            wp_map[jet_wp] = {}
                        wp_map[jet_wp][ele_wp] = genmatch_data
                        print(f"    Successfully extracted data for VSjet={jet_wp}, VSele={ele_wp}")
                    else:
                        print(f"    Warning: Could not extract genmatch data from {fname}")
                    
            except Exception as e:
                print(f"  Warning: Failed to read {fpath}: {e}")
                import traceback
                traceback.print_exc()

    if not wp_map:
        print("Error: No valid WP files found/parsed.")
        return None

    if reference_inputs is None:
        print("Error: Could not determine inputs from files.")
        return None

    print(f"\nBuilding combined correction with {sum(len(v) for v in wp_map.values())} WP combinations...")
    
    # Build the Category structure
    # Root: wp_VSmu -> wp_VSe -> wp_VSjet -> genmatch -> DM -> syst -> pT
    # Since we only have Tight for VSmu, we'll keep it simple
    
    # Build wp_VSjet categories for each ele_wp
    jet_cat_items = []
    for jet_wp in sorted(wp_map.keys()):
        ele_map = wp_map[jet_wp]
        
        ele_cat_items = []
        for ele_wp in sorted(ele_map.keys()):
            genmatch_data = ele_map[ele_wp]
            ele_cat_items.append(cs.CategoryItem(key=ele_wp, value=genmatch_data))
        
        # Create wp_VSe category
        ele_cat = cs.Category(
            nodetype="category",
            input="wp_VSe",
            content=ele_cat_items
        )
        
        jet_cat_items.append(cs.CategoryItem(key=jet_wp, value=ele_cat))

    # Create wp_VSjet category
    jet_cat = cs.Category(
        nodetype="category",
        input="wp_VSjet",
        content=jet_cat_items
    )
    
    # Wrap in wp_VSmu (assuming Tight is the only value)
    root_data = cs.Category(
        nodetype="category",
        input="wp_VSmu",
        content=[cs.CategoryItem(key="Tight", value=jet_cat)]
    )

    # Build the correction object
    combined_corr = cs.Correction(
        name=final_name,
        version=1,
        description=description,
        inputs=reference_inputs,
        output=cs.Variable(name="sf", type="real", description=f"{pattern_base} scale factor"),
        data=root_data
    )
    
    return combined_corr


def list_corrections_in_file(filename):
    """List all corrections in a JSON file."""
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
        print(f"\nCorrections in {os.path.basename(filename)}:")
        if 'corrections' in data:
            for corr in data['corrections']:
                print(f"  - {corr.get('name', 'unnamed')}")
                if 'inputs' in corr:
                    print(f"    Inputs: {[inp['name'] for inp in corr['inputs']]}")
    except Exception as e:
        print(f"ERROR reading {filename}: {e}")

def main():
    description = '''Create combined tau correction JSON files with variation parameters.
    
Examples:
  python3 merge_tau_jsons.py --type both                    # Merge all into one file
  python3 merge_tau_jsons.py --type tes -o tau_sf/TES.json  # Only TES
  python3 merge_tau_jsons.py --list-all                     # List all corrections
'''
    parser = ArgumentParser(prog="merge_tau_jsons", description=description)
    parser.add_argument('-i', '--input-dir', dest='input_dir', type=str, default='tau_sf/', 
                        help="Input directory containing JSON files")
    parser.add_argument('-o', '--output', dest='output_file', type=str, 
                        default='tau_sf/TauCorrections_2024.json',
                        help="Output JSON file")
    parser.add_argument('--type', dest='correction_type', type=str, choices=['tes', 'id', 'both'], 
                        default='both', help="Type of corrections to create: tes, id, or both")
    parser.add_argument('--list', dest='list_file', type=str, default=None,
                        help="List corrections in a specific JSON file")
    parser.add_argument('--list-all', dest='list_all', action='store_true',
                        help="List corrections in all found JSON files")
    
    args = parser.parse_args()
    
    if args.list_file:
        list_corrections_in_file(args.list_file)
        return
    
    if args.list_all:
        json_files = glob.glob(os.path.join(args.input_dir, "*.json"))
        for json_file in sorted(json_files):
            list_corrections_in_file(json_file)
        return
    
    # Create corrections based on type
    if args.correction_type == "tes":
        create_combined_correction(args.input_dir, args.output_file, "tes")
    elif args.correction_type == "id":
        create_combined_correction(args.input_dir, args.output_file, "id")
    else:  # both - create single file with both corrections
        create_combined_both_corrections(args.input_dir, args.output_file)


if __name__ == '__main__':
    main()