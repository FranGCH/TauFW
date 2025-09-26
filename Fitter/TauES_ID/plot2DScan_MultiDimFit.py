#!/usr/bin/env python
"""
Date : Sept 2025
Author : @haawedik based on plotParabola_POI_region.py by @oponcet
Description :
This script plots 2D parabolas from MultiDimFit output files when you have
scanned over two parameters simultaneously (like option 3 in your workflow).
This is more straightforward than extracting from FitDiagnostics.
"""

import sys
import os
import yaml
import ROOT
import numpy as np
from math import sqrt, pi
from argparse import ArgumentParser
from ROOT import gROOT, gPad, gStyle, TFile, TCanvas, TLegend, TLatex, TF2, TGraph2D, TH2D, TPolyMarker3D, TGraphAsymmErrors, TLine, TEllipse
from ROOT import kBlack, kBlue, kRed, kGreen, kYellow, kOrange, kMagenta, kTeal, kAzure, TMath
from TauFW.Plotter.sample.utils import CMSStyle

# Ensure ROOT runs in batch mode
gROOT.SetBatch(True)
gStyle.SetOptTitle(0)

# CMS style
CMSStyle.setTDRStyle()

def ensureDirectory(dirname):
    """Make directory if it does not exist."""
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def ensureTFile(filename, option='READ'):
    """Open TFile and make sure it exists."""
    if not os.path.isfile(filename):
        print(f"ERROR: File {filename} does not exist!")
        sys.exit(1)
    file = TFile(filename, option)
    if not file or file.IsZombie():
        print(f"ERROR: Could not open file {filename}")
        sys.exit(1)
    return file

def interpolate_scan_data(poi1_vals, poi2_vals, nll_vals, nbins=100):
    """Create a smoother 2D histogram by interpolating the scan data"""
    try:
        from scipy.interpolate import griddata
        from scipy.ndimage import gaussian_filter
    except ImportError:
        raise ImportError("scipy is required for interpolation")
    
    poi1_vals = np.array(poi1_vals)
    poi2_vals = np.array(poi2_vals)
    nll_vals = np.array(nll_vals)
    
    # Create regular grid with higher resolution
    poi1_min, poi1_max = np.min(poi1_vals), np.max(poi1_vals)
    poi2_min, poi2_max = np.min(poi2_vals), np.max(poi2_vals)
    
    # Add smaller padding for better contours
    poi1_range = poi1_max - poi1_min
    poi2_range = poi2_max - poi2_min
    padding = 0.05  # Reduced padding
    poi1_min -= padding * poi1_range
    poi1_max += padding * poi1_range
    poi2_min -= padding * poi2_range
    poi2_max += padding * poi2_range
    
    # Create high-resolution grid
    poi1_grid = np.linspace(poi1_min, poi1_max, nbins)
    poi2_grid = np.linspace(poi2_min, poi2_max, nbins)
    poi1_mesh, poi2_mesh = np.meshgrid(poi1_grid, poi2_grid)
    
    # Interpolate NLL values onto grid
    points = np.column_stack((poi1_vals, poi2_vals))
    
    # First try cubic interpolation
    try:
        nll_interpolated = griddata(points, nll_vals, (poi1_mesh, poi2_mesh), 
                                   method='cubic', fill_value=np.nan)
        # Fill NaN values with linear interpolation
        mask = np.isnan(nll_interpolated)
        if np.any(mask):
            nll_linear = griddata(points, nll_vals, (poi1_mesh, poi2_mesh), 
                                method='linear', fill_value=np.max(nll_vals))
            nll_interpolated[mask] = nll_linear[mask]
    except:
        # Fallback to linear interpolation
        nll_interpolated = griddata(points, nll_vals, (poi1_mesh, poi2_mesh), 
                                   method='linear', fill_value=np.max(nll_vals))
    
    # Apply Gaussian smoothing for even smoother contours
    nll_interpolated = gaussian_filter(nll_interpolated, sigma=1.0)
    
    return poi1_grid, poi2_grid, nll_interpolated

def calculate_correlation_and_uncertainties(poi1_vals, poi2_vals, nll_vals):
    """Calculate correlation and proper uncertainties from scan data - using 3σ filtering like FitDiagnostics"""
    poi1_vals = np.array(poi1_vals)
    poi2_vals = np.array(poi2_vals)
    nll_vals = np.array(nll_vals)
    
    print(f">>> Debug correlation calculation:")
    print(f"    Total scan points: {len(poi1_vals)}")
    print(f"    POI1 range: [{np.min(poi1_vals):.4f}, {np.max(poi1_vals):.4f}]")
    print(f"    POI2 range: [{np.min(poi2_vals):.4f}, {np.max(poi2_vals):.4f}]")
    print(f"    NLL range: [{np.min(nll_vals):.4f}, {np.max(nll_vals):.4f}]")
    
    # Find the best fit point
    min_idx = np.argmin(nll_vals)
    best_poi1 = poi1_vals[min_idx]
    best_poi2 = poi2_vals[min_idx]
    
    print(f"    Best fit: POI1={best_poi1:.4f}, POI2={best_poi2:.4f}")
    
    # **Use same filtering as FitDiagnostics: points within 3σ (deltaNLL < 9)**
    # Try different thresholds, but prefer 3σ
    correlation = 0.0
    poi1_err = 0.0
    poi2_err = 0.0
    
    for threshold in [9.0, 10.0, 20.0]:  # 3σ first, then relax if needed
        mask = nll_vals < threshold
        filtered_poi1 = poi1_vals[mask]
        filtered_poi2 = poi2_vals[mask]
        filtered_nll = nll_vals[mask]
        
        print(f"    Points with deltaNLL < {threshold}: {np.sum(mask)} / {len(poi1_vals)}")
        
        if len(filtered_poi1) < 10:
            print(f"    WARNING: Not enough scan points for reliable correlation calculation at {threshold}σ")
            continue
        
        # **Same method as FitDiagnostics: direct correlation coefficient of filtered points**
        if len(np.unique(filtered_poi1)) > 1 and len(np.unique(filtered_poi2)) > 1:
            correlation = np.corrcoef(filtered_poi1, filtered_poi2)[0, 1]
            print(f"    Raw correlation (deltaNLL < {threshold}): {correlation:.4f}")
            
            if not np.isnan(correlation) and abs(correlation) <= 1.0:
                # Calculate uncertainties using likelihood weighting of filtered points
                weights = np.exp(-0.5 * filtered_nll)
                weights = weights / np.sum(weights)
                
                # Weighted mean and std
                poi1_mean = np.sum(weights * filtered_poi1)
                poi2_mean = np.sum(weights * filtered_poi2)
                
                poi1_var = np.sum(weights * (filtered_poi1 - poi1_mean)**2)
                poi2_var = np.sum(weights * (filtered_poi2 - poi2_mean)**2)
                
                poi1_err = np.sqrt(poi1_var)
                poi2_err = np.sqrt(poi2_var)
                
                # Ensure reasonable minimum uncertainty
                poi1_range = np.max(poi1_vals) - np.min(poi1_vals)
                poi2_range = np.max(poi2_vals) - np.min(poi2_vals)
                poi1_err = max(poi1_err, poi1_range / 100.0)
                poi2_err = max(poi2_err, poi2_range / 100.0)
                
                print(f"    Using threshold deltaNLL < {threshold}")
                break
            else:
                print(f"    Invalid correlation at threshold {threshold}, trying next...")
        else:
            print(f"    Insufficient parameter variation at threshold {threshold}")
    
    # Final fallback if all thresholds failed
    if correlation == 0.0 or np.isnan(correlation) or abs(correlation) > 1.0:
        print("    WARNING: All correlation calculations failed, using fallback method")
        
        # Use all points as last resort
        if len(np.unique(poi1_vals)) > 1 and len(np.unique(poi2_vals)) > 1:
            correlation = np.corrcoef(poi1_vals, poi2_vals)[0, 1]
            if np.isnan(correlation) or abs(correlation) > 1.0:
                correlation = 0.0
        else:
            correlation = 0.0
        
        # Simple uncertainty estimates
        poi1_err = (np.max(poi1_vals) - np.min(poi1_vals)) / 6.0  # Range/6 ≈ 1σ for normal distribution
        poi2_err = (np.max(poi2_vals) - np.min(poi2_vals)) / 6.0
    
    print(f"    Final correlation: {correlation:.4f}")
    print(f"    Uncertainties: {poi1_err:.4f}, {poi2_err:.4f}")
    
    return correlation, poi1_err, poi2_err

def extract_2d_scan_data(multidimfit_file, poi1_name, poi2_name):
    """
    Extract 2D scan data from MultiDimFit output file.
    Returns arrays of parameter values and corresponding NLL values.
    """
    print(f">>> Reading 2D scan data from {multidimfit_file}")
    
    file = ensureTFile(multidimfit_file)
    tree = file.Get('limit')
    if not tree:
        print("ERROR: Could not find 'limit' tree in MultiDimFit file")
        file.Close()
        return None
    
    # First, let's see what branches are available
    print(">>> Available branches:")
    branches = []
    for branch in tree.GetListOfBranches():
        branch_name = branch.GetName()
        branches.append(branch_name)
        print(f"  {branch_name}")
    
    poi1_vals = []
    poi2_vals = []
    nll_vals = []
    
    nentries = tree.GetEntries()
    print(f">>> Processing {nentries} entries from MultiDimFit scan")
    
    # Based on your ROOT output, the branches are directly named
    poi1_branch = poi1_name  # 'tes_DM0'
    poi2_branch = poi2_name  # 'tid_SF_DM0' 
    nll_branch = 'deltaNLL'
    
    # Verify branches exist
    if poi1_branch not in branches:
        print(f"ERROR: Branch '{poi1_branch}' not found")
        file.Close()
        return None
    if poi2_branch not in branches:
        print(f"ERROR: Branch '{poi2_branch}' not found") 
        file.Close()
        return None
    if nll_branch not in branches:
        print(f"ERROR: Branch '{nll_branch}' not found")
        file.Close()
        return None
    
    print(f">>> Using branches: {poi1_branch}, {poi2_branch}, {nll_branch}")
    
    for i in range(nentries):
        tree.GetEntry(i)
        
        # Get parameter values using the found branch names
        poi1_val = getattr(tree, poi1_branch)
        poi2_val = getattr(tree, poi2_branch)
        nll_val = getattr(tree, nll_branch)
        
        poi1_vals.append(poi1_val)
        poi2_vals.append(poi2_val)
        nll_vals.append(nll_val)
        
        # Debug first few entries
        if i < 5:
            print(f"  Entry {i}: {poi1_branch}={poi1_val:.4f}, {poi2_branch}={poi2_val:.4f}, {nll_branch}={nll_val:.4f}")
    
    file.Close()
    
    if len(poi1_vals) == 0:
        print(f"ERROR: No valid data found")
        return None
    
    # Convert to numpy arrays for easier manipulation
    import numpy as np
    poi1_vals = np.array(poi1_vals)
    poi2_vals = np.array(poi2_vals)
    nll_vals = np.array(nll_vals)
    
    # Debug scan pattern
    print(f">>> Scan pattern analysis:")
    print(f"    Unique POI1 values: {len(np.unique(poi1_vals))}")
    print(f"    Unique POI2 values: {len(np.unique(poi2_vals))}")
    unique_poi1 = np.unique(poi1_vals)
    unique_poi2 = np.unique(poi2_vals)
    expected_points = len(unique_poi1) * len(unique_poi2)
    print(f"    Expected grid points: {expected_points}, Actual points: {len(poi1_vals)}")
    print(f"    POI1 range: [{np.min(poi1_vals):.4f}, {np.max(poi1_vals):.4f}]")
    print(f"    POI2 range: [{np.min(poi2_vals):.4f}, {np.max(poi2_vals):.4f}]")
    print(f"    NLL range: [{np.min(nll_vals):.4f}, {np.max(nll_vals):.4f}]")
    
    # Check if this is actually a 2D scan or two 1D scans
    if len(poi1_vals) == len(unique_poi1) + len(unique_poi2):
        print("    WARNING: This looks like two separate 1D scans, not a 2D scan!")
        print("    You need to run: combine -M MultiDimFit --algo=grid --redefineSignalPOIs poi1,poi2")
    elif len(poi1_vals) < expected_points * 0.8:
        print("    WARNING: Scan appears incomplete or not fully gridded")
    else:
        print("    This appears to be a proper 2D grid scan")
    
    # deltaNLL is already the correct quantity (no need to subtract minimum)
    delta_nll_vals = nll_vals
    
    # Find best fit point
    min_idx = np.argmin(delta_nll_vals)
    best_poi1 = poi1_vals[min_idx]
    best_poi2 = poi2_vals[min_idx]
    
    # Calculate correlation and uncertainties
    correlation, poi1_err, poi2_err = calculate_correlation_and_uncertainties(
        poi1_vals, poi2_vals, delta_nll_vals)
    
    print(f">>> Found {len(poi1_vals)} data points")
    print(f">>> Best fit: {poi1_name} = {best_poi1:.4f}, {poi2_name} = {best_poi2:.4f}")
    print(f">>> NLL range: {np.min(delta_nll_vals):.3f} to {np.max(delta_nll_vals):.3f}")
    print(f">>> Calculated correlation: {correlation:.4f}")
    print(f">>> {poi1_name} = {best_poi1:.4f} ± {poi1_err:.4f}")
    print(f">>> {poi2_name} = {best_poi2:.4f} ± {poi2_err:.4f}")
    
    return {
        'poi1_name': poi1_name,
        'poi2_name': poi2_name,
        'poi1_vals': poi1_vals,
        'poi2_vals': poi2_vals,
        'delta_nll_vals': delta_nll_vals,
        'best_poi1': best_poi1,
        'best_poi2': best_poi2,
        'poi1_err': poi1_err,
        'poi2_err': poi2_err,
        'correlation': correlation
    }

def plot_2d_scan(setup, region, year, scan_data, **kwargs):
    """
    Plot 2D parabola from MultiDimFit scan data.
    """
    print(f">>> Plotting 2D scan for {region}")
    
    indir = kwargs.get('indir', f"output_{year}")
    outdir = indir.replace('output', 'plots')
    tag = kwargs.get('tag', "")
    plottag = kwargs.get('plottag', "")
    poi1_name = scan_data['poi1_name']
    poi2_name = scan_data['poi2_name']
    era = f"{year}-13TeV"
    channel = setup["channel"].replace("mu", "m").replace("tau", "t")
    
    ensureDirectory(outdir)
    
    # Canvas name
    canvasname = f"{outdir}/scan_2D_{poi1_name}_{poi2_name}_{channel}_{region}{tag}{plottag}"
    
    # Get data arrays
    poi1_vals = scan_data['poi1_vals']
    poi2_vals = scan_data['poi2_vals']
    delta_nll_vals = scan_data['delta_nll_vals']
    best_poi1 = scan_data['best_poi1']
    best_poi2 = scan_data['best_poi2']
    
    # Determine plot ranges
    poi1_min, poi1_max = np.min(poi1_vals), np.max(poi1_vals)
    poi2_min, poi2_max = np.min(poi2_vals), np.max(poi2_vals)
    
    # Add some margin
    poi1_range = poi1_max - poi1_min
    poi2_range = poi2_max - poi2_min
    margin = 0.1
    
    poi1_min -= margin * poi1_range
    poi1_max += margin * poi1_range
    poi2_min -= margin * poi2_range
    poi2_max += margin * poi2_range
    
    # Create 2D histogram with higher resolution
    nbins_x = 100
    nbins_y = 100
    hist_2d = TH2D("hist_2d", "", nbins_x, poi1_min, poi1_max, nbins_y, poi2_min, poi2_max)
    
    # Always try interpolation for smoother results
    use_interpolation = True
    
    try:
        # Try to use scipy interpolation for smoother results
        poi1_grid, poi2_grid, nll_interpolated = interpolate_scan_data(
            poi1_vals, poi2_vals, delta_nll_vals, nbins_x)
        
        # Fill histogram with interpolated data
        for i in range(nbins_x):
            for j in range(nbins_y):
                hist_2d.SetBinContent(i+1, j+1, nll_interpolated[j, i])
        print(">>> Using interpolated scan data for smoother plot")
        
    except (ImportError, Exception) as e:
        print(f">>> Interpolation failed ({e}), using direct binning with heavy smoothing")
        use_interpolation = False
        
        # Fill histogram using the scan data (direct binning)
        for i in range(len(poi1_vals)):
            bin_x = hist_2d.GetXaxis().FindBin(poi1_vals[i])
            bin_y = hist_2d.GetYaxis().FindBin(poi2_vals[i])
            
            # Set bin content to minimum of current content and new value
            # (in case multiple scan points fall in the same bin)
            current_content = hist_2d.GetBinContent(bin_x, bin_y)
            if current_content == 0 or delta_nll_vals[i] < current_content:
                hist_2d.SetBinContent(bin_x, bin_y, delta_nll_vals[i])
        
        # For empty bins, interpolate from nearby filled bins
        # This is a simple nearest-neighbor interpolation
        for i in range(1, nbins_x + 1):
            for j in range(1, nbins_y + 1):
                if hist_2d.GetBinContent(i, j) == 0:
                    x_center = hist_2d.GetXaxis().GetBinCenter(i)
                    y_center = hist_2d.GetYaxis().GetBinCenter(j)
                    
                    # Find nearest scan point
                    distances = np.sqrt((poi1_vals - x_center)**2 + (poi2_vals - y_center)**2)
                    nearest_idx = np.argmin(distances)
                    hist_2d.SetBinContent(i, j, delta_nll_vals[nearest_idx])
        
        # Apply heavy smoothing to get rid of the grid pattern
        for _ in range(5):  # Multiple smoothing iterations
            hist_2d.Smooth(1)
    
    # Set up canvas
    canvas = TCanvas('canvas', 'canvas', 100, 100, 800, 700)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetTopMargin(0.07)
    canvas.SetBottomMargin(0.12)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.15)  # Space for color scale
    canvas.cd()
    
    # Set axis titles
    if 'tes' in poi1_name.lower():
        x_title = "Tau energy scale"
    elif 'tid_sf' in poi1_name.lower():
        x_title = "Tau ID scale factor"
    else:
        x_title = poi1_name
        
    if 'tes' in poi2_name.lower():
        y_title = "Tau energy scale"
    elif 'tid_sf' in poi2_name.lower():
        y_title = "Tau ID scale factor"
    else:
        y_title = poi2_name
    
    hist_2d.GetXaxis().SetTitle(x_title)
    hist_2d.GetYaxis().SetTitle(y_title)
    hist_2d.GetXaxis().SetTitleSize(0.055)
    hist_2d.GetYaxis().SetTitleSize(0.055)
    hist_2d.GetXaxis().SetLabelSize(0.050)
    hist_2d.GetYaxis().SetLabelSize(0.050)
    hist_2d.GetZaxis().SetTitle("-2#Deltaln(L)")
    hist_2d.GetZaxis().SetTitleSize(0.055)
    hist_2d.GetZaxis().SetLabelSize(0.050)
    
    # Set up color palette for better visualization
    hist_2d.SetMinimum(0)
    max_val = hist_2d.GetMaximum()
    if max_val > 20:
        hist_2d.SetMaximum(20)  # Cap the maximum for better color scale
    
    # Draw the 2D histogram with smoother color transitions
    hist_2d.Draw("COLZ")
    
    # Add contour lines for 1σ, 2σ, 3σ confidence levels for 2D
    # For 2D: 1σ = 2.30, 2σ = 6.18, 3σ = 11.83 (for -2ΔlnL)
    contour_levels = [2.30, 6.18, 11.83]  # 68%, 95%, 99.73% confidence levels for 2D
    
    # Create a separate histogram for contours to avoid interference
    hist_contour = hist_2d.Clone("hist_contour")
    hist_contour.SetContour(len(contour_levels))
    for i, level in enumerate(contour_levels):
        hist_contour.SetContourLevel(i, level)
    
    # Draw smooth contour lines
    hist_contour.SetLineColor(kBlack)
    hist_contour.SetLineWidth(3)
    hist_contour.Draw("CONT3 SAME")
    
    # Overlay the actual scan points
    graph_points = ROOT.TGraph(len(poi1_vals), poi1_vals, poi2_vals)
    graph_points.SetMarkerStyle(20)
    graph_points.SetMarkerSize(0.3)
    graph_points.SetMarkerColor(kBlue)
    graph_points.Draw("P SAME")
    
    # Mark the best fit point
    best_fit_marker = ROOT.TMarker(best_poi1, best_poi2, 29)  # Large star
    best_fit_marker.SetMarkerColor(kRed)
    best_fit_marker.SetMarkerSize(2.0)
    best_fit_marker.Draw("SAME")
    
    # Add text with fit results
    latex = TLatex()
    latex.SetTextSize(0.04)
    latex.SetTextAlign(11)
    latex.SetNDC(True)
    
    text_x = 0.15
    text_y = 0.85
    line_height = 0.05
    
    # Region title
    if region in setup.get("regions", {}):
        region_title = setup["regions"][region]["title"]
    else:
        region_title = region
    latex.DrawLatex(text_x, text_y, f"Region: {region_title}")
    
    # Best fit values with uncertainties
    latex.DrawLatex(text_x, text_y - line_height, 
                   f"{x_title}: {best_poi1:.4f} #pm {scan_data['poi1_err']:.4f}")
    latex.DrawLatex(text_x, text_y - 2*line_height, 
                   f"{y_title}: {best_poi2:.4f} #pm {scan_data['poi2_err']:.4f}")
    latex.DrawLatex(text_x, text_y - 3*line_height, 
                   f"Correlation: {scan_data['correlation']:.3f}")
    latex.DrawLatex(text_x, text_y - 4*line_height, 
                   f"Scan points: {len(poi1_vals)}")
    
    # Add legend for contour lines and points
    legend = TLegend(0.15, 0.50, 0.50, 0.67)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.AddEntry(best_fit_marker, "Best fit", "p")
    legend.AddEntry(graph_points, "Scan points", "p")
    legend.AddEntry(hist_contour, "68%, 95%, 99.7% CL", "l")
    legend.Draw()
    
    # CMS style
    CMSStyle.setCMSLumiStyle(canvas, 0)
    canvas.SetTicks(1, 1)
    canvas.Modified()
    canvas.Update()
    
    # Save the plot
    canvas.SaveAs(canvasname + ".png")
    canvas.SaveAs(canvasname + ".pdf")
    canvas.SaveAs(canvasname + ".root")
    
    print(f">>> Saved 2D scan plot: {canvasname}.png")
    
    canvas.Close()

def main(args):
    """Main function"""
    
    print("Using configuration file: %s" % args.config)
    with open(args.config, 'r') as file:
        setup = yaml.safe_load(file)
    
    channel = setup["channel"].replace("mu", "m").replace("tau", "t")
    tag = setup.get("tag", "")
    era = args.year
    extratag = "_DeepTau"
    
    # Input directory
    if args.indir:
        if not args.indir.rstrip('/').endswith(str(era)):
            indir = os.path.join(args.indir, str(era))
        else:
            indir = args.indir
    else:
        indir = f"output_{era}"
    
    # POI names
    poi1_name = args.poi1  # e.g., "tes_DM0"
    poi2_name = args.poi2  # e.g., "tid_SF_DM0"
    
    # Extract region from POI name (assuming format like "tes_DM0")
    if '_' in poi1_name:
        region = poi1_name.split('_', 1)[1]
    else:
        region = args.region if args.region else "DM0"
    
    # Construct MultiDimFit filename
    # This should match what your makecombinedfitTES_SF.py produces for option 3
    multidimfit_filename = f"{indir}/higgsCombine.{channel}_m_vis-{region}{tag}{extratag}-{era}-13TeV.MultiDimFit.mH90.root"
    
    if not os.path.exists(multidimfit_filename):
        print(f"ERROR: MultiDimFit file not found: {multidimfit_filename}")
        # Try alternative naming
        multidimfit_filename = f"{indir}/higgsCombine.mt_m_vis-{region}{tag}{extratag}-{era}-13TeV.MultiDimFit.mH90.root"
        if not os.path.exists(multidimfit_filename):
            print(f"ERROR: Alternative MultiDimFit file not found: {multidimfit_filename}")
            sys.exit(1)
    
    print(f">>> Using MultiDimFit file: {multidimfit_filename}")
    
    # Extract scan data
    scan_data = extract_2d_scan_data(multidimfit_filename, poi1_name, poi2_name)
    if scan_data is None:
        print("ERROR: Could not extract scan data")
        sys.exit(1)
    
    # Create 2D plot
    plot_2d_scan(setup, region, era, scan_data, 
                 indir=indir, tag=tag, plottag=args.plottag)

if __name__ == '__main__':
    description = '''Plot 2D parabolas from MultiDimFit scan output.'''
    parser = ArgumentParser(prog="plot2DScan_MultiDimFit", 
                          description=description, epilog="Success!")
    
    parser.add_argument('-y', '--year', dest='year', 
                       choices=['2024', '2016', '2017', '2018', 'UL2016_preVFP', 
                               'UL2016_postVFP', 'UL2017', 'UL2018', 'UL2018_v10',
                               '2022_postEE', '2022_preEE', '2023C', '2023D'], 
                       type=str, default='2024', action='store', 
                       help="select year")
    
    parser.add_argument('-c', '--config', dest='config', type=str, 
                       default='TauES_ID/config/config_coarse_TT.yml', 
                       action='store', 
                       help="set config file containing sample & fit setup")
    
    parser.add_argument('--poi1', dest='poi1', type=str, required=True,
                       help="first parameter of interest (e.g., tes_DM0)")
    
    parser.add_argument('--poi2', dest='poi2', type=str, required=True,
                       help="second parameter of interest (e.g., tid_SF_DM0)")
    
    parser.add_argument('-r', '--region', dest='region', type=str,
                       help="region name (if not extractable from POI names)")
    
    parser.add_argument('-i', '--indir', dest='indir', type=str, 
                       help='input directory')
    
    parser.add_argument('-t', '--plottag', dest='plottag', type=str, 
                       default="", help='extra tag for plot filename')
    
    parser.add_argument('-v', '--verbose', dest='verbose', 
                       default=False, action='store_true', help="set verbose")
    
    args = parser.parse_args()
    main(args)