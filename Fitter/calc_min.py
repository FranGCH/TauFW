import ROOT
import argparse
import ctypes
import numpy as np
import pandas as pd 
import re
import json
import os
import glob

def SampleName(HistName):
    process_name = re.search(r'.+_(.+)',HistName)
    return process_name.group(1) if process_name else HistName

def OpenFile(filepath):
    file = ROOT.TFile.Open(filepath, "READ")
    if not file or file.IsZombie():
        raise FileNotFoundError(f"Could not open file: {filepath}")
    return file


#get the histgram
def Get_histogram(file, hist_name='canvas'):
    hist = file.Get(hist_name)
    if not hist:
        raise ValueError(f"Histogram '{hist_name}' not found in the file.")
    return hist

def RootTH2D_to_PythonHist(filepath, hist_name='hist_2d'):
    """
    Convert a ROOT TH2D histogram to a Python numpy array representation
    """
    file = ROOT.TFile.Open(filepath, "READ")
    if not file or file.IsZombie():
        raise FileNotFoundError(f"Could not open file: {filepath}")
    
    # Get the histogram from canvas if needed
    canvas = file.Get('canvas')
    th2d = None
    
    if canvas:
        primitives = canvas.GetListOfPrimitives()
        for obj in primitives:
            if obj.InheritsFrom("TH2"):
                print("Found TH2D in canvas:", obj.GetName())
                th2d = obj
                break
    else:
        th2d = file.Get(hist_name)
    
    if not th2d:
        raise ValueError(f"TH2D histogram not found")
    
    # Get dimensions
    nx = th2d.GetNbinsX()
    ny = th2d.GetNbinsY()
    
    # Create numpy array to store histogram values
    hist_array = np.zeros((nx, ny))
    
    # Fill the array with bin contents
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            hist_array[i-1, j-1] = th2d.GetBinContent(i, j)
    
    # Also get axis information
    x_axis = th2d.GetXaxis()
    y_axis = th2d.GetYaxis()
    
    x_edges = np.array([x_axis.GetBinLowEdge(i) for i in range(1, nx + 2)])
    y_edges = np.array([y_axis.GetBinLowEdge(i) for i in range(1, ny + 2)])
    
    file.Close()
    
    return {
        'data': hist_array,
        'x_edges': x_edges,
        'y_edges': y_edges,
        'x_label': x_axis.GetTitle(),
        'y_label': y_axis.GetTitle(),
        'title': th2d.GetTitle()
    }

def main():
    j_wp_choices = ['VVLoose', 'VLoose', 'Loose', 'Medium', 'Tight','VTight']
    e_wp_choices = ['VVLoose', 'Tight']
    parser = argparse.ArgumentParser(description="Get minimum value from 2D histogram in a ROOT file.")
    # parser.add_argument("--filepath", type=str, help="Path to the ROOT file.")
    parser.add_argument("--jet_wp", nargs="*", choices=j_wp_choices,default=j_wp_choices, help="Jet working point: [VVLoose, VLoose, Loose, Medium, Tight, VTight].")
    parser.add_argument("--ele_wp", nargs="*", choices=e_wp_choices,default=e_wp_choices, help="Electron working point: [VVLoose, Tight].")
    # parser.add_argument("--hist_name", type=str, help="Name of the 2D histogram.")
    args = parser.parse_args()
    'Fitter/plots_pt_less_region/againstjet_VVLoose/againstelectron_VVLoose/2024/scan_2D_tes_DM1_tid_SF_DM1_mt_DM1_mutaumultidimfit.root'

    jet_wp_list = args.jet_wp
    ele_wp_list = args.ele_wp
    dic_file_info = {}

    for jet_wp in jet_wp_list:
        for ele_wp in ele_wp_list:
            filepath_multy = glob.glob(f'./plots_pt_less_region/againstjet_{jet_wp}/againstelectron_{ele_wp}/2024/scan_2D*mutaumultidimfit.root')
            for filepath in filepath_multy:
                if not os.path.exists(filepath):
                    print(f"File not found: {filepath}")
                    continue

                dic_min_info={}

                py_hist = RootTH2D_to_PythonHist(filepath)
                print(f"Histogram shape: {py_hist['data'].shape}")
                print(f"Min value: {py_hist['data'].min()}")
                for i in range(py_hist['data'].shape[0]):
                    for j in range(py_hist['data'].shape[1]):
                        if py_hist['data'][i,j] == py_hist['data'].min():
                            dic_min_info['bin']=(i+1,j+1)
                            dic_min_info['min']=py_hist['data'][i,j]
                            dic_min_info['TES_edges']=(py_hist['x_edges'][i],py_hist['x_edges'][i+1])
                            dic_min_info['TauID_edges']=(py_hist['y_edges'][j],py_hist['y_edges'][j+1])

                            print(f"Min value located at bin ({i+1}, {j+1}) with edges x:"\
                                f"[{py_hist['x_edges'][i]}, {py_hist['x_edges'][i+1]}],"\
                                    f"y: [{py_hist['y_edges'][j]}, {py_hist['y_edges'][j+1]}]")

                filepath_parts = filepath.split('/')
                for part in filepath_parts:
                    if 'scan' in part:
                        dmX_ptY = re.search(r'scan_2D_tes_(.+)_tid_SF*.', part)
                        dmX_ptY = dmX_ptY.group(1) if dmX_ptY else "unknown"
                
                file_info = f'{dmX_ptY}_jet{jet_wp}_ele{ele_wp}'
                
                dic_file_info[file_info] = dic_min_info
                
    if os.path.exists('multiDimFit_min_info_int2.json'):
        with open('multiDimFit_min_info_int2.json', 'r') as json_file:
            existing_data = json.load(json_file)
        existing_data.update(dic_file_info)
        with open('multiDimFit_min_info_int2.json', 'w') as json_file:
            json.dump(existing_data, json_file, indent=4)
    else:
        with open('multiDimFit_min_info_int2.json', 'w') as json_file:
            json.dump(dic_file_info, json_file, indent=4)


if __name__ == "__main__":
    main()