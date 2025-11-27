#!/usr/bin/env python
# Script to plot correlations and measurements from 2D measurement files
# Creates 3 plots: correlations, TES values, and TauID values

import ROOT
ROOT.gROOT.SetBatch(True)  # Run in batch mode (no GUI)
import os
import glob
import re
from ROOT import TCanvas, TGraph, TGraphAsymmErrors, TLatex, TLegend, TLine, kBlue, kRed, kGreen, kMagenta, kBlack, kOrange, kGray

def load_measurements():
    """Load measurements from 2D measurement files"""
    
    # Pattern to find all measurement files (both inclusive and pt-binned)
    pattern = "output_pt_less_region/againstjet_*/againstelectron_*/2024/FitparameterValues__mutau_DeepTau_2024-13TeV_*.txt"
    all_files = glob.glob(pattern)
    
    print(f"Found {len(all_files)} measurement files (inclusive and pt-binned)")
    
    measurements = []
    
    for filename in all_files:
        if os.path.getsize(filename) == 0:
            print(f"[WARNING] Skipping empty file: {filename}")
            continue
            
        print(f"Processing: {filename}")
        
        try:
            # Extract DM, pt bin, and working points from filename
            basename = os.path.basename(filename)
            tes_dm_match = re.search(r'(DM\d+)(?:_(pt\d+))?', basename)
            jet_wp_match = re.search(r'againstjet_(\w+)', filename)
            ele_wp_match = re.search(r'againstelectron_(\w+)', filename)
            
            if not (tes_dm_match and jet_wp_match and ele_wp_match):
                print(f"Could not parse filename: {basename}")
                continue
                
            dm = tes_dm_match.group(1)
            pt_bin = tes_dm_match.group(2) if tes_dm_match.group(2) else "inclusive"
            jet_wp = jet_wp_match.group(1)
            ele_wp = ele_wp_match.group(1)
            
            if pt_bin == "inclusive":
                region_name = f"{dm} inclusive"
            else:
                # Convert pt bin number to actual pt range (assuming standard binning)
                pt_ranges = {
                    'pt1': '20-40 GeV', 
                    'pt2': '40-60 GeV',
                    'pt3': '60-200 GeV',
                }
                pt_label = pt_ranges.get(pt_bin, pt_bin)
                region_name = f"{dm} {pt_label}"
            
            # Read the measurement file
            with open(filename, 'r') as f:
                lines = f.readlines()
            
            tes_val = None
            tes_err_up = None
            tes_err_down = None
            tid_val = None
            tid_err_up = None
            tid_err_down = None
            correlation = 0.0 # Default to 0.0 as it is missing in the new file format
            
            # Temporary variables for bounds (new format)
            tes_low = None
            tes_high = None
            tid_low = None
            tid_high = None
            
            for line in lines:
                line = line.strip()
                if line.startswith('#') or line == '':
                    continue
                
                # Handle new format: key: value
                if ':' in line:
                    parts = line.split(':')
                    if len(parts) >= 2:
                        key = parts[0].strip()
                        try:
                            val = float(parts[1].strip())
                            
                            if key.startswith('tes_'):
                                if '1sigma_low' in key:
                                    tes_low = val
                                elif '1sigma_high' in key:
                                    tes_high = val
                                else:
                                    tes_val = val
                            elif key.startswith('tid_SF_'):
                                if '1sigma_low' in key:
                                    tid_low = val
                                elif '1sigma_high' in key:
                                    tid_high = val
                                else:
                                    tid_val = val
                        except ValueError:
                            pass

                # Handle old format: name val err_down err_up
                else:
                    parts = line.split()
                    if len(parts) >= 4:
                        param_name = parts[0]
                        try:
                            value = float(parts[1])
                            error_down = float(parts[2])
                            error_up = float(parts[3])
                            
                            if param_name.startswith('tes_'):
                                tes_val = value
                                tes_err_down = error_down
                                tes_err_up = error_up
                            elif param_name.startswith('tid_SF_'):
                                tid_val = value
                                tid_err_down = error_down
                                tid_err_up = error_up
                        except ValueError:
                            pass
                # elif line.startswith('correlation'):
                #     parts = line.split()
                #     if len(parts) >= 2:
                #         correlation = float(parts[1])
            
            # Calculate errors from bounds if using new format
            if tes_val is not None and tes_low is not None and tes_high is not None:
                tes_err_down = abs(tes_val - tes_low)
                tes_err_up = abs(tes_high - tes_val)
                
            if tid_val is not None and tid_low is not None and tid_high is not None:
                tid_err_down = abs(tid_val - tid_low)
                tid_err_up = abs(tid_high - tid_val)
            
            if all(x is not None for x in [tes_val, tid_val, correlation]):
                measurements.append({
                    'region': region_name,
                    'dm': dm,
                    'pt_bin': pt_bin,
                    'jet_wp': jet_wp,
                    'ele_wp': ele_wp,
                    'tes_val': tes_val,
                    'tes_err_down': tes_err_down,
                    'tes_err_up': tes_err_up,
                    'tid_val': tid_val,
                    'tid_err_down': tid_err_down,
                    'tid_err_up': tid_err_up,
                    'correlation': correlation
                })
                print(f"  {region_name}: TES={tes_val:.4f} -{tes_err_down:.4f}/+{tes_err_up:.4f}, TauID={tid_val:.4f} -{tid_err_down:.4f}/+{tid_err_up:.4f}, Corr={correlation:.4f}")
            else:
                print(f"  Incomplete data in {basename}")
                
        except Exception as e:
            print(f"Error processing {filename}: {e}")
    
    print(f"\nLoaded {len(measurements)} complete measurements")
    return measurements

def create_correlation_plot(measurements):
    """Create correlation plot"""
    print("\nCreating correlation plot...")
    
    c1 = TCanvas("c_corr", "TES-TauID Correlations", 800, 600)
    c1.SetMargin(0.45, 0.05, 0.15, 0.08)  # Much larger left margin for longer labels
    c1.SetGrid()
    
    # Create graph with error bars
    n_points = len(measurements)
    gr = TGraphAsymmErrors(n_points)
    
    # Color mapping for decay modes
    dm_colors = {'DM0': kBlack, 'DM1': kRed, 'DM10': kBlue, 'DM11': kGreen}
    
    # Create y-axis labels and set points
    y_labels = []
    for i, meas in enumerate(measurements):
        gr.SetPoint(i, meas['correlation'], i)
        # Set small horizontal error bars (you can adjust these if you have correlation uncertainties)
        gr.SetPointError(i, 0.02, 0.02, 0, 0)  # Small horizontal errors
        
        # Create clean region labels
        if meas['pt_bin'] == 'inclusive':
            label = meas['dm']
        else:
            # Use the actual pt range from the mapping
            pt_ranges = {
                'pt1': '20-40 GeV', 
                'pt2': '40-60 GeV',
                'pt3': '60-200 GeV',
            }
            pt_range = pt_ranges.get(meas['pt_bin'], meas['pt_bin'])
            label = f"{meas['dm']} {pt_range}"
        y_labels.append(label)
    
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1.0)
    gr.SetMarkerColor(kBlack)
    gr.SetLineColor(kBlack)
    gr.SetTitle("")
    gr.GetXaxis().SetTitle("Correlation Coefficient")
    gr.GetXaxis().SetTitleSize(0.045)
    gr.GetXaxis().SetLabelSize(0.04)
    gr.GetXaxis().SetRangeUser(-1.0, 1.0)
    gr.GetYaxis().SetRangeUser(-0.5, len(measurements) - 0.5)
    
    # Remove y-axis labels and ticks
    gr.GetYaxis().SetLabelSize(0)
    gr.GetYaxis().SetTickLength(0)
    
    gr.Draw("AP")
    
    # Add custom y-axis labels positioned within the plot area
    for i, label in enumerate(y_labels):
        text = TLatex()
        text.SetTextSize(0.03)  # Smaller text
        text.SetTextAlign(32)  # Right aligned
        # Position labels at x = -0.95 to ensure they're within the plot bounds
        text.DrawLatex(-0.95, i, label)
    
    # Add correlation values next to each point
    for i, meas in enumerate(measurements):
        corr_text = TLatex()
        corr_text.SetTextSize(0.025)  # Slightly smaller text to fit more digits
        corr_text.SetTextAlign(12)  # Left aligned
        # Position text to the right of the point
        x_pos = meas['correlation'] + 0.05
        if x_pos > 0.9:  # If too close to right edge, put on left
            x_pos = meas['correlation'] - 0.05
            corr_text.SetTextAlign(32)  # Right aligned
        corr_text.DrawLatex(x_pos, i, f"{meas['correlation']:.6f}")
    
    # Add vertical line at correlation = 0
    line_zero = ROOT.TLine(0, -0.5, 0, len(measurements) - 0.5)
    line_zero.SetLineStyle(2)
    line_zero.SetLineColor(ROOT.kGray+1)
    line_zero.Draw()
    
    # Add CMS label
    cms_label = TLatex()
    cms_label.SetNDC()
    cms_label.SetTextFont(61)
    cms_label.SetTextSize(0.05)
    cms_label.DrawLatex(0.46, 0.92, "CMS")
    
    cms_internal = TLatex()
    cms_internal.SetNDC()
    cms_internal.SetTextFont(52)
    cms_internal.SetTextSize(0.04)
    cms_internal.DrawLatex(0.55, 0.92, "Internal")
    
    lumi_text = TLatex()
    lumi_text.SetNDC()
    lumi_text.SetTextFont(42)
    lumi_text.SetTextSize(0.035)
    lumi_text.DrawLatex(0.68, 0.92, "109 fb^{-1} (13.6 TeV)")
    
    c1.SaveAs("correlation_plot.png")
    c1.SaveAs("correlation_plot.pdf")
    c1.SaveAs("correlation_plot.root")
    print("Correlation plot saved as correlation_plot.png/pdf")
    
    return c1

def create_tes_plot(measurements):
    """Create TES measurements plot"""
    print("\nCreating TES plot...")
    
    c2 = TCanvas("c_tes", "TES Measurements", 800, 600)
    c2.SetMargin(0.35, 0.05, 0.15, 0.08)  # Increased left margin for longer labels
    c2.SetGrid()
    
    # Create graph with error bars
    n_points = len(measurements)
    gr = TGraphAsymmErrors(n_points)
    
    # Create y-axis labels and set points
    y_labels = []
    for i, meas in enumerate(measurements):
        gr.SetPoint(i, meas['tes_val'], i)
        gr.SetPointError(i, meas['tes_err_down'], meas['tes_err_up'], 0, 0)
        
        # Create clean region labels
        if meas['pt_bin'] == 'inclusive':
            label = meas['dm']
        else:
            # Use the actual pt range from the mapping
            pt_ranges = {
                'pt1': '20-40 GeV', 
                'pt2': '40-60 GeV',
                'pt3': '60-200 GeV',
            }
            pt_range = pt_ranges.get(meas['pt_bin'], meas['pt_bin'])
            label = f"{meas['dm']} {pt_range}"
        y_labels.append(label)
    
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1.0)
    gr.SetMarkerColor(kBlack)
    gr.SetLineColor(kBlack)
    gr.SetTitle("")
    gr.GetXaxis().SetTitle("TES Scale Factor")
    gr.GetXaxis().SetTitleSize(0.045)
    gr.GetXaxis().SetLabelSize(0.04)
    
    # Find x-range
    tes_values = [m['tes_val'] for m in measurements]
    tes_errors = [max(m['tes_err_up'], m['tes_err_down']) for m in measurements]
    x_min = min(v - e for v, e in zip(tes_values, tes_errors)) * 0.98
    x_max = max(v + e for v, e in zip(tes_values, tes_errors)) * 1.02
    gr.GetXaxis().SetRangeUser(x_min, x_max)
    gr.GetYaxis().SetRangeUser(-0.5, len(measurements) - 0.5)
    
    # Remove y-axis labels and ticks
    gr.GetYaxis().SetLabelSize(0)
    gr.GetYaxis().SetTickLength(0)
    
    gr.Draw("AP")
    
    # Add custom y-axis labels
    for i, label in enumerate(y_labels):
        text = TLatex()
        text.SetTextSize(0.035)
        text.SetTextAlign(32)  # Right aligned
        text.DrawLatex(x_min - 0.05 * (x_max - x_min), i, label)  # More space for labels
    
    # Add vertical line at TES = 1.0
    line_unity = ROOT.TLine(1.0, -0.5, 1.0, len(measurements) - 0.5)
    line_unity.SetLineStyle(2)
    line_unity.SetLineColor(ROOT.kGray+1)
    line_unity.Draw()
    
    # Add CMS label
    cms_label = TLatex()
    cms_label.SetNDC()
    cms_label.SetTextFont(61)
    cms_label.SetTextSize(0.05)
    cms_label.DrawLatex(0.26, 0.92, "CMS")
    
    cms_internal = TLatex()
    cms_internal.SetNDC()
    cms_internal.SetTextFont(52)
    cms_internal.SetTextSize(0.04)
    cms_internal.DrawLatex(0.35, 0.92, "Internal")
    
    lumi_text = TLatex()
    lumi_text.SetNDC()
    lumi_text.SetTextFont(42)
    lumi_text.SetTextSize(0.035)
    lumi_text.DrawLatex(0.65, 0.92, "109 fb^{-1} (13.6 TeV)")
    
    c2.SaveAs("tes_measurements.png")
    c2.SaveAs("tes_measurements.pdf")
    c2.SaveAs("tes_measurements.root")
    print("TES plot saved as tes_measurements.png/pdf")
    
    return c2

def create_tauID_plot(measurements):
    """Create TauID measurements plot"""
    print("\nCreating TauID plot...")
    
    c3 = TCanvas("c_tid", "TauID Measurements", 800, 600)
    c3.SetMargin(0.35, 0.05, 0.15, 0.08)  # Increased left margin for longer labels
    c3.SetGrid()
    
    # Create graph with error bars
    n_points = len(measurements)
    gr = TGraphAsymmErrors(n_points)
    
    # Create y-axis labels and set points
    y_labels = []
    for i, meas in enumerate(measurements):
        gr.SetPoint(i, meas['tid_val'], i)
        gr.SetPointError(i, meas['tid_err_down'], meas['tid_err_up'], 0, 0)
        
        # Create clean region labels
        if meas['pt_bin'] == 'inclusive':
            label = meas['dm']
        else:
            # Use the actual pt range from the mapping
            pt_ranges = {
                'pt1': '20-40 GeV', 
                'pt2': '40-60 GeV',
                'pt3': '60-200 GeV',
            }
            pt_range = pt_ranges.get(meas['pt_bin'], meas['pt_bin'])
            label = f"{meas['dm']} {pt_range}"
        y_labels.append(label)
    
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(1.0)
    gr.SetMarkerColor(kBlack)
    gr.SetLineColor(kBlack)
    gr.SetTitle("")
    gr.GetXaxis().SetTitle("TauID Scale Factor")
    gr.GetXaxis().SetTitleSize(0.045)
    gr.GetXaxis().SetLabelSize(0.04)
    
    # Find x-range
    tid_values = [m['tid_val'] for m in measurements]
    tid_errors = [max(m['tid_err_up'], m['tid_err_down']) for m in measurements]
    x_min = min(v - e for v, e in zip(tid_values, tid_errors)) * 0.95
    x_max = max(v + e for v, e in zip(tid_values, tid_errors)) * 1.05
    gr.GetXaxis().SetRangeUser(x_min, x_max)
    gr.GetYaxis().SetRangeUser(-0.5, len(measurements) - 0.5)
    
    # Remove y-axis labels and ticks
    gr.GetYaxis().SetLabelSize(0)
    gr.GetYaxis().SetTickLength(0)
    
    gr.Draw("AP")
    
    # Add custom y-axis labels
    for i, label in enumerate(y_labels):
        text = TLatex()
        text.SetTextSize(0.035)
        text.SetTextAlign(32)  # Right aligned
        text.DrawLatex(x_min - 0.08 * (x_max - x_min), i, label)  # More space for labels
    
    # Add vertical line at TauID SF = 1.0
    line_unity = ROOT.TLine(1.0, -0.5, 1.0, len(measurements) - 0.5)
    line_unity.SetLineStyle(2)
    line_unity.SetLineColor(ROOT.kGray+1)
    line_unity.Draw()
    
    # Add CMS label
    cms_label = TLatex()
    cms_label.SetNDC()
    cms_label.SetTextFont(61)
    cms_label.SetTextSize(0.05)
    cms_label.DrawLatex(0.26, 0.92, "CMS")
    
    cms_internal = TLatex()
    cms_internal.SetNDC()
    cms_internal.SetTextFont(52)
    cms_internal.SetTextSize(0.04)
    cms_internal.DrawLatex(0.35, 0.92, "Internal")
    
    lumi_text = TLatex()
    lumi_text.SetNDC()
    lumi_text.SetTextFont(42)
    lumi_text.SetTextSize(0.035)
    lumi_text.DrawLatex(0.65, 0.92, "109 fb^{-1} (13.6 TeV)")
    
    c3.SaveAs("tauID_measurements.png")
    c3.SaveAs("tauID_measurements.pdf")
    c3.SaveAs("tauID_measurements.root")
    print("TauID plot saved as tauID_measurements.png/pdf")
    
    return c3

def main():
    """Main function"""
    print("Loading measurements from 2D measurement files...")
    
    # Load measurements
    measurements = load_measurements()
    
    if not measurements:
        print("No measurements found!")
        return
    
    # Sort measurements by decay mode for consistent plotting
    measurements.sort(key=lambda x: (x['dm'], x['jet_wp'], x['ele_wp']))
    
    print(f"\nCreating plots for {len(measurements)} measurements...")
    
    # Create the three plots
    c1 = create_correlation_plot(measurements)
    c2 = create_tes_plot(measurements)
    c3 = create_tauID_plot(measurements)
    
    print("\nAll plots created successfully!")
    print("Files saved:")
    print("  - correlation_plot.png/pdf")
    print("  - tes_measurements.png/pdf") 
    print("  - tauID_measurements.png/pdf")
    
    print("\\nPlots completed successfully!")

if __name__ == "__main__":
    main()