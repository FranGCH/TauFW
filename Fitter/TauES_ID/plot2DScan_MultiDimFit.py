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

def format_region_label(region):
    """Format region name to match CMS style (e.g., DM0, DM1_pt1 -> DM1, pt: 20-40)"""
    if '_pt' in region:
        dm_part, pt_part = region.split('_pt', 1)
        pt_num = pt_part
        # Map pt bins to actual pT ranges
        pt_ranges = {
            '1': '20-40',
            '2': '40-60', 
            '3': '60-200',
            '4': '200+'
        }
        pt_range = pt_ranges.get(pt_num, pt_num)
        
        # Format decay mode 
        if dm_part == 'DM0':
            return f"{dm_part}, pt: {pt_range}"
        elif dm_part == 'DM1':
            return f"{dm_part}, pt: {pt_range}" 
        elif dm_part == 'DM10':
            return f"{dm_part}, pt: {pt_range}"
        elif dm_part == 'DM11':
            return f"{dm_part}, pt: {pt_range}"
        else:
            return f"{dm_part}, pt: {pt_range}"
    else:
        # Just DM labels
        return region

def format_region_for_sorting(region):
    """Create sorting key for regions to match your plot order"""
    # Define the desired order (reverse: DM11, DM10, DM1, DM0)
    order_map = {
        'DM11': 0,
        'DM10': 1, 
        'DM1': 2,
        'DM0': 3
    }
    
    if '_pt' in region:
        dm_part, pt_part = region.split('_pt', 1)
        base_order = order_map.get(dm_part, 999)
        # pt bins in reverse order (pt4, pt3, pt2, pt1)
        pt_order = 10 - int(pt_part) if pt_part.isdigit() else 0
        return (base_order, pt_order)
    else:
        return (order_map.get(region, 999), 0)

def interpolate_scan_data(poi1_vals, poi2_vals, nll_vals, nbins=200):
    """Create a smoother 2D histogram by interpolating the scan data"""
    try:
        from scipy.interpolate import griddata
        from scipy.ndimage import gaussian_filter
        from scipy.interpolate import Rbf
    except ImportError:
        raise ImportError("scipy is required for interpolation (griddata / Rbf / gaussian_filter)")
     
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
    
    max_nll = float(np.nanmax(nll_vals))
    # Try cubic interpolation with a finite fill_value to avoid NaNs
    try:
        nll_interpolated = griddata(points, nll_vals, (poi1_mesh, poi2_mesh),
                                    method='cubic', fill_value=max_nll)
        # If cubic left NaNs, fill those with linear
        if np.any(~np.isfinite(nll_interpolated)):
            nll_linear = griddata(points, nll_vals, (poi1_mesh, poi2_mesh),
                                  method='linear', fill_value=max_nll)
            nll_interpolated[~np.isfinite(nll_interpolated)] = nll_linear[~np.isfinite(nll_interpolated)]
    except Exception:
        # As a more robust fallback try radial-basis interpolation (RBF) which extrapolates better
        try:
            rbf = Rbf(poi1_vals, poi2_vals, nll_vals, function='linear')
            nll_interpolated = rbf(poi1_mesh, poi2_mesh)
        except Exception:
            # Last resort: linear griddata with finite fill
            nll_interpolated = griddata(points, nll_vals, (poi1_mesh, poi2_mesh),
                                        method='linear', fill_value=max_nll)

    # Ensure no NaNs remain (safety) and treat outside-hull as high cost
    nll_interpolated = np.where(np.isfinite(nll_interpolated), nll_interpolated, max_nll)
    # Gentle Gaussian smoothing to keep small-scale features (smaller sigma for higher accuracy)
    nll_interpolated = gaussian_filter(nll_interpolated, sigma=0.8)
     
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
    thresholds = [9.0, 10.0, 20.0, 1.0]
    min_points_required = 5
    for threshold in thresholds:
        mask = nll_vals < threshold
        filtered_poi1 = poi1_vals[mask]
        filtered_poi2 = poi2_vals[mask]
        filtered_nll = nll_vals[mask]
        
        nkept = np.sum(mask)
        print(f"    Points with deltaNLL < {threshold}: {nkept} / {len(poi1_vals)}")
        if nkept < min_points_required:
            print(f"    WARNING: Not enough scan points (<{min_points_required}) at threshold {threshold}, trying next...")
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
                poi1_err = max(poi1_err, poi1_range / 200.0)  # tighter minimum for higher-accuracy scans
                poi2_err = max(poi2_err, poi2_range / 200.0)
                
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
    
    # Based on ROOT output, the branches are directly named
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
    
    # Create 2D histogram with higher resolution for more accurate contours
    nbins_x = 200
    nbins_y = 200
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

def plot_measurement_summary(region_labels, measurements, **kwargs):
    """Create a summary plot showing measurements with error bars"""
    title = kwargs.get('title', "Measurements")
    ylabel = kwargs.get('ylabel', "value")
    outname = kwargs.get('outname', "measurements")
    year = kwargs.get('year', "2024")
    
    n_regions = len(region_labels)
    if n_regions == 0:
        print("No measurements to plot")
        return
    
    # Create canvas
    canvas_height = max(600, 60 + 40*n_regions)
    canvas_width = 800
    canvas = ROOT.TCanvas('canvas_summary', 'canvas_summary', 100, 100, canvas_width, canvas_height)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    
    # Set margins
    top_margin = 0.08
    bottom_margin = 0.12
    left_margin = 0.25
    right_margin = 0.05
    
    canvas.SetTopMargin(top_margin)
    canvas.SetBottomMargin(bottom_margin) 
    canvas.SetLeftMargin(left_margin)
    canvas.SetRightMargin(right_margin)
    canvas.SetGrid(1, 0)
    canvas.cd()
    
    # Determine x-axis range
    values = [m[0] for m in measurements]
    errors_down = [m[1] for m in measurements] 
    errors_up = [m[2] for m in measurements]
    
    x_min = min([v - e for v, e in zip(values, errors_down)])
    x_max = max([v + e for v, e in zip(values, errors_up)])
    x_range = x_max - x_min
    x_margin = 0.15 * x_range
    x_min -= x_margin
    x_max += x_margin
    
    # Create frame
    frame = canvas.DrawFrame(x_min, 0.0, x_max, float(n_regions))
    frame.GetYaxis().SetLabelSize(0.0)
    frame.GetXaxis().SetLabelSize(0.042)
    frame.GetXaxis().SetTitleSize(0.050) 
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetYaxis().SetNdivisions(n_regions, 0, 0, False)
    frame.GetXaxis().SetTitle(ylabel)
    frame.GetXaxis().SetNdivisions(510)
    
    # Create graph with error bars
    graph = ROOT.TGraphAsymmErrors(n_regions)
    
    for i, (region, measurement) in enumerate(zip(region_labels, measurements)):
        y_pos = n_regions - i - 0.5
        val, err_down, err_up = measurement
        graph.SetPoint(i, val, y_pos)
        graph.SetPointError(i, err_down, err_up, 0.1, 0.1)
    
    # Style the graph
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.SetMarkerColor(ROOT.kBlack)
    graph.SetLineColor(ROOT.kBlack)
    graph.SetLineWidth(2)
    
    # Draw the graph
    graph.Draw("PE SAME")
    
    # Add vertical line at 1.0 if appropriate
    if min(values) < 1.0 < max(values):
        line = ROOT.TLine(1.0, 0.0, 1.0, float(n_regions))
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        line.Draw("SAME")
    
    # Add region labels
    latex = ROOT.TLatex()
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.SetTextAlign(32)
    latex.SetNDC(True)
    
    for i, region in enumerate(region_labels):
        y_pos_ndc = 1.0 - top_margin - (i + 0.5) * (1.0 - top_margin - bottom_margin) / n_regions
        latex.DrawLatex(left_margin - 0.02, y_pos_ndc, region)
    
    # Add CMS header
    cms_latex = ROOT.TLatex()
    cms_latex.SetTextSize(0.060)
    cms_latex.SetTextFont(61)
    cms_latex.SetTextAlign(11)
    cms_latex.SetNDC(True)
    cms_latex.DrawLatex(left_margin, 1.0 - top_margin + 0.01, "CMS")
    
    # Add "Internal" label
    internal_latex = ROOT.TLatex()
    internal_latex.SetTextSize(0.045)
    internal_latex.SetTextFont(52)
    internal_latex.SetTextAlign(11)
    internal_latex.SetNDC(True) 
    internal_latex.DrawLatex(left_margin + 0.12, 1.0 - top_margin + 0.01, "Internal")
    
    # Add year and energy
    year_latex = ROOT.TLatex()
    year_latex.SetTextSize(0.045)
    year_latex.SetTextFont(42)
    year_latex.SetTextAlign(31)
    year_latex.SetNDC(True)
    year_latex.DrawLatex(1.0 - right_margin, 1.0 - top_margin + 0.01, f"{year}, 109 fb^{{-1}} (13.6 TeV)")
    
    # Add title
    if title:
        title_latex = ROOT.TLatex()
        title_latex.SetTextSize(0.045)
        title_latex.SetTextFont(42)
        title_latex.SetTextAlign(11)
        title_latex.SetNDC(True)
        title_latex.DrawLatex(left_margin, 1.0 - top_margin - 0.05, title)
    
    canvas.SetTicks(1, 1)
    canvas.Modified()
    canvas.Update()
    
    # Save
    canvas.SaveAs(outname + ".png")
    canvas.SaveAs(outname + ".pdf") 
    canvas.SaveAs(outname + ".root")
    print(f">>> Saved measurement summary: {outname}.png")
    
    canvas.Close()

def plot_scan_correlations(scan_results_all_regions, **kwargs):
    """Create a correlation plot using correlations from scan results"""
    year = kwargs.get('year', '2024')
    indir = kwargs.get('indir', f"output_{year}")
    outdir = indir.replace('output', 'plots')
    tag = kwargs.get('tag', "")
    plottag = kwargs.get('plottag', "")
    outname = f"{outdir}/scan_correlations_multidimfit{tag}{plottag}"
    
    ensureDirectory(outdir)
    
    # Extract correlations and region names directly from scan results
    correlations = []
    region_labels = []
    
    # Sort regions
    sorted_regions = sorted(scan_results_all_regions.items(), key=lambda x: format_region_for_sorting(x[0]))
    
    for region, scan_data in sorted_regions:
        if scan_data is None:
            continue
        correlation = scan_data.get('correlation', 0.0)
        correlations.append(correlation)
        region_labels.append(format_region_label(region))
    
    n_regions = len(region_labels)
    if n_regions == 0:
        print("No correlation data to plot")
        return
    
    # Create canvas
    canvas_height = max(600, 60 + 40*n_regions)
    canvas_width = 800
    canvas = ROOT.TCanvas('canvas_scan_corr', 'canvas_scan_corr', 100, 100, canvas_width, canvas_height)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    
    # Set margins
    top_margin = 0.08
    bottom_margin = 0.12
    left_margin = 0.25
    right_margin = 0.05
    
    canvas.SetTopMargin(top_margin)
    canvas.SetBottomMargin(bottom_margin) 
    canvas.SetLeftMargin(left_margin)
    canvas.SetRightMargin(right_margin)
    canvas.SetGrid(1, 0)
    canvas.cd()
    
    # Determine x-axis range for correlations (-1 to +1)
    x_min = -1.2
    x_max = 1.2
    
    # Create frame
    frame = canvas.DrawFrame(x_min, 0.0, x_max, float(n_regions))
    frame.GetYaxis().SetLabelSize(0.0)
    frame.GetXaxis().SetLabelSize(0.042)
    frame.GetXaxis().SetTitleSize(0.050) 
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetYaxis().SetNdivisions(n_regions, 0, 0, False)
    frame.GetXaxis().SetTitle("TES-TauID Correlation from 2D Scans")
    frame.GetXaxis().SetNdivisions(510)
    
    # Create individual markers for each correlation
    markers = []
    
    for i, (region, correlation) in enumerate(zip(region_labels, correlations)):
        y_pos = n_regions - i - 0.5
        
        # Simple color scheme: blue for all correlations
        color = ROOT.kBlue + 2
        
        # Create marker
        marker = ROOT.TMarker(correlation, y_pos, 20)
        marker.SetMarkerSize(1.2)
        marker.SetMarkerColor(color)
        markers.append(marker)
        marker.Draw("SAME")
    
    # Add vertical line at 0.0 (no correlation)
    line_zero = ROOT.TLine(0.0, 0.0, 0.0, float(n_regions))
    line_zero.SetLineStyle(2)
    line_zero.SetLineColor(ROOT.kGray+2)
    line_zero.SetLineWidth(2)
    line_zero.Draw("SAME")
    
    # Add region labels
    latex = ROOT.TLatex()
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.SetTextAlign(32)
    latex.SetNDC(True)
    
    for i, region in enumerate(region_labels):
        y_pos_ndc = 1.0 - top_margin - (i + 0.5) * (1.0 - top_margin - bottom_margin) / n_regions
        latex.DrawLatex(left_margin - 0.02, y_pos_ndc, region)
    
    # Add correlation values next to points
    corr_latex = ROOT.TLatex()
    corr_latex.SetTextSize(0.030)
    corr_latex.SetTextFont(42)
    corr_latex.SetTextAlign(11)
    
    for i, correlation in enumerate(correlations):
        y_pos = n_regions - i - 0.5
        x_pos = correlation + 0.05 if correlation >= 0 else correlation - 0.05
        corr_latex.DrawLatex(x_pos, y_pos, f"{correlation:.3f}")
    
    # Add CMS header
    cms_latex = ROOT.TLatex()
    cms_latex.SetTextSize(0.060)
    cms_latex.SetTextFont(61)
    cms_latex.SetTextAlign(11)
    cms_latex.SetNDC(True)
    cms_latex.DrawLatex(left_margin, 1.0 - top_margin + 0.01, "CMS")
    
    # Add "Internal" label
    internal_latex = ROOT.TLatex()
    internal_latex.SetTextSize(0.045)
    internal_latex.SetTextFont(52)
    internal_latex.SetTextAlign(11)
    internal_latex.SetNDC(True) 
    internal_latex.DrawLatex(left_margin + 0.12, 1.0 - top_margin + 0.01, "Internal")
    
    # Add year and energy
    year_latex = ROOT.TLatex()
    year_latex.SetTextSize(0.045)
    year_latex.SetTextFont(42)
    year_latex.SetTextAlign(31)
    year_latex.SetNDC(True)
    year_latex.DrawLatex(1.0 - right_margin, 1.0 - top_margin + 0.01, f"{year}, 109 fb^{{-1}} (13.6 TeV)")
    
    # Add title
    title_latex = ROOT.TLatex()
    title_latex.SetTextSize(0.045)
    title_latex.SetTextFont(42)
    title_latex.SetTextAlign(11)
    title_latex.SetNDC(True)
    title_latex.DrawLatex(left_margin, 1.0 - top_margin - 0.05, "TES-TauID Correlations from MultiDimFit Scans")
    
    canvas.SetTicks(1, 1)
    canvas.Modified()
    canvas.Update()
    
    # Save
    canvas.SaveAs(outname + ".png")
    canvas.SaveAs(outname + ".pdf") 
    canvas.SaveAs(outname + ".root")
    print(f">>> Saved scan correlation plot: {outname}.png")
    
    canvas.Close()

def plot_summary_from_multiple_regions(setup, scan_results_all_regions, **kwargs):
    """Create summary plots showing TES and TauID measurements with uncertainties"""
    print(">>> Creating summary plots from MultiDimFit results")
    
    year = kwargs.get('year', '2024')
    indir = kwargs.get('indir', f"output_{year}")
    outdir = indir.replace('output', 'plots')
    tag = kwargs.get('tag', "")
    plottag = kwargs.get('plottag', "")
    
    ensureDirectory(outdir)
    
    # Extract measurements
    tes_measurements = []
    tid_measurements = []
    region_labels = []
    
    # Sort regions
    sorted_regions = sorted(scan_results_all_regions.items(), key=lambda x: format_region_for_sorting(x[0]))
    
    for region, scan_data in sorted_regions:
        if scan_data is None:
            continue
            
        region_label = format_region_label(region)
        region_labels.append(region_label)
        
        poi1_name = scan_data['poi1_name']
        poi2_name = scan_data['poi2_name']
        
        if 'tes' in poi1_name.lower():
            tes_val = scan_data['best_poi1']
            tes_err = scan_data['poi1_err']
            tid_val = scan_data['best_poi2'] 
            tid_err = scan_data['poi2_err']
        else:
            tes_val = scan_data['best_poi2']
            tes_err = scan_data['poi2_err']
            tid_val = scan_data['best_poi1']
            tid_err = scan_data['poi1_err']
            
        tes_measurements.append((tes_val, tes_err, tes_err))
        tid_measurements.append((tid_val, tid_err, tid_err))
    
    # Create plots
    plot_measurement_summary(
        region_labels, tes_measurements,
        title="Tau Energy Scale",
        ylabel="tau energy scale", 
        outname=f"{outdir}/tes_summary_multidimfit{tag}{plottag}",
        year=year
    )
    
    plot_measurement_summary(
        region_labels, tid_measurements,
        title="Tau ID Scale Factor", 
        ylabel="tau ID scale factor",
        outname=f"{outdir}/tid_summary_multidimfit{tag}{plottag}",
        year=year
    )
    
    # Create correlation plot
    plot_scan_correlations(
        scan_results_all_regions,
        year=year, 
        indir=indir, 
        tag=tag, 
        plottag=plottag
    )
    
    print(f">>> Created TES summary plot: {outdir}/tes_summary_multidimfit{tag}{plottag}.png")
    print(f">>> Created TauID summary plot: {outdir}/tid_summary_multidimfit{tag}{plottag}.png")
    print(f">>> Created scan correlation plot: {outdir}/scan_correlations_multidimfit{tag}{plottag}.png")

def write_2d_fit_results(poi1_name, poi2_name, poi1_val, poi1_err_down, poi1_err_up, 
                         poi2_val, poi2_err_down, poi2_err_up, correlation, region, **kwargs):
    """Write 2D fit results to text file, similar to measurepoi() in plotParabola_POI_region.py"""
    year = kwargs.get('year', '2024')
    tag = kwargs.get('tag', '')
    channel = kwargs.get('channel', 'mt')
    outdir = kwargs.get('outdir', 'plots')
    
    ensureDirectory(outdir)
    
    # Create output filename similar to plotParabola_POI_region format
    outfname = f"{outdir}/measurement_2D_{poi1_name}_{poi2_name}_{channel}_{region}{tag}.txt"
    
    print(f">>> Writing 2D fit results to {outfname}")
    
    # Write results to file in same format as plotParabola_POI_region
    with open(outfname, 'w') as file:
        file.write("# 2D Fit Results from MultiDimFit\n")
        file.write(f"# Region: {region}\n")
        file.write(f"# Year: {year}\n")
        file.write(f"# Channel: {channel}\n")
        file.write("# Format: parameter value error_down error_up\n")
        file.write(f"{poi1_name} {poi1_val:.6f} {poi1_err_down:.6f} {poi1_err_up:.6f}\n")
        file.write(f"{poi2_name} {poi2_val:.6f} {poi2_err_down:.6f} {poi2_err_up:.6f}\n")
        file.write(f"correlation {correlation:.6f}\n")
    
    return outfname

def write_summary_results(scan_results_all_regions, **kwargs):
    """Write summary of all 2D fit results to a single file"""
    year = kwargs.get('year', '2024')
    tag = kwargs.get('tag', '')
    channel = kwargs.get('channel', 'mt')
    outdir = kwargs.get('outdir', 'plots')
    
    ensureDirectory(outdir)
    
    # Create summary output filename
    outfname = f"{outdir}/summary_2D_results_{channel}{tag}_{year}.txt"
    
    print(f">>> Writing summary 2D fit results to {outfname}")
    
    with open(outfname, 'w') as file:
        file.write("# Summary of 2D Fit Results from MultiDimFit\n")
        file.write(f"# Year: {year}\n")
        file.write(f"# Channel: {channel}\n")
        file.write("# Format: region poi1_name poi1_val poi1_err poi2_name poi2_val poi2_err correlation\n")
        
        # Sort regions for consistent output
        sorted_regions = sorted(scan_results_all_regions.items(), key=lambda x: format_region_for_sorting(x[0]))
        
        for region, scan_data in sorted_regions:
            if scan_data is None:
                continue
                
            # Get parameter info
            poi1_name = scan_data['poi1_name']
            poi2_name = scan_data['poi2_name']
            poi1_val = scan_data['best_poi1']
            poi2_val = scan_data['best_poi2']
            poi1_err = scan_data['poi1_err']
            poi2_err = scan_data['poi2_err']
            correlation = scan_data.get('correlation', 0.0)
            
            # Write to file
            file.write(f"{region} {poi1_name} {poi1_val:.6f} {poi1_err:.6f} ")
            file.write(f"{poi2_name} {poi2_val:.6f} {poi2_err:.6f} ")
            file.write(f"{correlation:.6f}\n")
    
    return outfname

def main(args):
    """Main function - handle multiple regions and create summary plots"""
    
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
    
    # Process regions
    scan_results_all_regions = {}
    
    if args.poi1 and args.poi2:
        # Single region mode
        poi1_name = args.poi1  # e.g., "tes_DM0"
        poi2_name = args.poi2  # e.g., "tid_SF_DM0"
        
        # Extract region from POI name (assuming format like "tes_DM0")
        if '_' in poi1_name:
            region = poi1_name.split('_', 1)[1]
        else:
            region = args.region if args.region else "DM0"
            
        regions_to_process = [region]
    else:
        # Multiple regions mode - get from config
        try:
            regions_to_process = setup["observables"]["m_vis"]["scanRegions"]
            print(f">>> Processing {len(regions_to_process)} regions from config: {regions_to_process}")
        except KeyError:
            print("ERROR: No regions found in config file. Please specify --poi1 and --poi2 for single region mode.")
            sys.exit(1)
    
    # Process each region
    for region in regions_to_process:
        print(f"\n>>> Processing region: {region}")
        
        # Construct POI names if not provided
        if not args.poi1 or not args.poi2:
            poi1_name = f"tes_{region}"
            poi2_name = f"tid_SF_{region}"
        else:
            poi1_name = args.poi1
            poi2_name = args.poi2
        
        # Construct MultiDimFit filename
        multidimfit_filename = f"{indir}/higgsCombine.{channel}_m_vis-{region}{tag}{extratag}-{era}-13TeV.MultiDimFit.mH90.root"
        
        if not os.path.exists(multidimfit_filename):
            # Try alternative naming
            multidimfit_filename = f"{indir}/higgsCombine.mt_m_vis-{region}{tag}{extratag}-{era}-13TeV.MultiDimFit.mH90.root"
            if not os.path.exists(multidimfit_filename):
                print(f"WARNING: MultiDimFit file not found for {region}")
                scan_results_all_regions[region] = None
                continue
        
        print(f">>> Using: {multidimfit_filename}")
        
        # Extract scan data
        scan_data = extract_2d_scan_data(multidimfit_filename, poi1_name, poi2_name)
        if scan_data is None:
            print(f"ERROR: Could not extract scan data for {region}")
            scan_results_all_regions[region] = None
            continue
            
        scan_results_all_regions[region] = scan_data
        
        # Create individual 2D plot
        plot_2d_scan(setup, region, era, scan_data, 
                     indir=indir, tag=tag, plottag=args.plottag)
        
        # Write individual text results for this region
        outdir = indir.replace('output', 'plots')
        poi1_val = scan_data['best_poi1']
        poi1_err = scan_data['poi1_err']
        poi2_val = scan_data['best_poi2'] 
        poi2_err = scan_data['poi2_err']
        correlation = scan_data.get('correlation', 0.0)
        
        # Write individual region results
        write_2d_fit_results(
            poi1_name, poi2_name, poi1_val, poi1_err, poi1_err,  # symmetric errors
            poi2_val, poi2_err, poi2_err, correlation, region,
            year=era, tag=tag, channel=channel, outdir=outdir
        )
    
    # Create summary plots and text files if we have multiple regions
    if len([r for r in scan_results_all_regions.values() if r is not None]) > 1:
        plot_summary_from_multiple_regions(setup, scan_results_all_regions,
                                         year=era, indir=indir, tag=tag, 
                                         plottag=args.plottag)
        
        # Write summary text results
        outdir = indir.replace('output', 'plots')
        write_summary_results(scan_results_all_regions,
                             year=era, tag=tag, channel=channel, outdir=outdir)
    elif len(scan_results_all_regions) == 1:
        print(">>> Only one region processed, summary plots not created")
    
    print(">>> All plots and text files completed successfully!")

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
    
    parser.add_argument('--poi1', dest='poi1', type=str, required=False,
                       help="first parameter of interest (e.g., tes_DM0). If not provided, will process all regions from config.")
    
    parser.add_argument('--poi2', dest='poi2', type=str, required=False,
                       help="second parameter of interest (e.g., tid_SF_DM0). If not provided, will process all regions from config.")
    
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
# #!/usr/bin/env python
# """
# Date : Sept 2025
# Author : @haawedik based on plotParabola_POI_region.py by @oponcet
# Description :
# This script plots 2D parabolas from MultiDimFit output files when you have
# scanned over two parameters simultaneously (like option 3 in your workflow).
# This is more straightforward than extracting from FitDiagnostics.
# """

# import sys
# import os
# import yaml
# import ROOT
# import numpy as np
# from math import sqrt, pi
# from argparse import ArgumentParser
# from ROOT import gROOT, gPad, gStyle, TFile, TCanvas, TLegend, TLatex, TF2, TGraph2D, TH2D, TPolyMarker3D, TGraphAsymmErrors, TLine, TEllipse
# from ROOT import kBlack, kBlue, kRed, kGreen, kYellow, kOrange, kMagenta, kTeal, kAzure, TMath
# from TauFW.Plotter.sample.utils import CMSStyle

# # Ensure ROOT runs in batch mode
# gROOT.SetBatch(True)
# gStyle.SetOptTitle(0)

# # CMS style
# CMSStyle.setTDRStyle()

# def ensureDirectory(dirname):
#     """Make directory if it does not exist."""
#     if not os.path.exists(dirname):
#         os.makedirs(dirname)

# def ensureDirectory(dirname):
#     """Make directory if it does not exist."""
#     if not os.path.exists(dirname):
#         os.makedirs(dirname)

# def ensureTFile(filename, option='READ'):
#     """Open TFile and make sure it exists."""
#     if not os.path.isfile(filename):
#         print(f"ERROR: File {filename} does not exist!")
#         sys.exit(1)
#     file = TFile(filename, option)
#     if not file or file.IsZombie():
#         print(f"ERROR: Could not open file {filename}")
#         sys.exit(1)
#     return file

# def format_region_label(region):
#     """Format region name to match CMS style (e.g., DM0, DM1_pt1 -> DM1, pt: 20-40)"""
#     if '_pt' in region:
#         dm_part, pt_part = region.split('_pt', 1)
#         pt_num = pt_part
#         # Map pt bins to actual pT ranges
#         pt_ranges = {
#             '1': '20-40',
#             '2': '40-60', 
#             '3': '60-200',
#             '4': '200+'
#         }
#         pt_range = pt_ranges.get(pt_num, pt_num)
        
#         # Format decay mode 
#         if dm_part == 'DM0':
#             return f"{dm_part}, pt: {pt_range}"
#         elif dm_part == 'DM1':
#             return f"{dm_part}, pt: {pt_range}" 
#         elif dm_part == 'DM10':
#             return f"{dm_part}, pt: {pt_range}"
#         elif dm_part == 'DM11':
#             return f"{dm_part}, pt: {pt_range}"
#         else:
#             return f"{dm_part}, pt: {pt_range}"
#     else:
#         # Just DM labels
#         return region

# def format_region_for_sorting(region):
#     """Create sorting key for regions to match your plot order"""
#     # Define the desired order (reverse: DM11, DM10, DM1, DM0)
#     order_map = {
#         'DM11': 0,
#         'DM10': 1, 
#         'DM1': 2,
#         'DM0': 3
#     }
    
#     if '_pt' in region:
#         dm_part, pt_part = region.split('_pt', 1)
#         base_order = order_map.get(dm_part, 999)
#         # pt bins in reverse order (pt4, pt3, pt2, pt1)
#         pt_order = 10 - int(pt_part) if pt_part.isdigit() else 0
#         return (base_order, pt_order)
#     else:
#         return (order_map.get(region, 999), 0)

# def interpolate_scan_data(poi1_vals, poi2_vals, nll_vals, nbins=200):
#     """Create a smoother 2D histogram by interpolating the scan data"""
#     try:
#         from scipy.interpolate import griddata
#         from scipy.ndimage import gaussian_filter
#         from scipy.interpolate import Rbf
#     except ImportError:
#         raise ImportError("scipy is required for interpolation (griddata / Rbf / gaussian_filter)")
     
#     poi1_vals = np.array(poi1_vals)
#     poi2_vals = np.array(poi2_vals)
#     nll_vals = np.array(nll_vals)
    
#     # Create regular grid with higher resolution
#     poi1_min, poi1_max = np.min(poi1_vals), np.max(poi1_vals)
#     poi2_min, poi2_max = np.min(poi2_vals), np.max(poi2_vals)
    
#     # Add smaller padding for better contours
#     poi1_range = poi1_max - poi1_min
#     poi2_range = poi2_max - poi2_min
#     padding = 0.05  # Reduced padding
#     poi1_min -= padding * poi1_range
#     poi1_max += padding * poi1_range
#     poi2_min -= padding * poi2_range
#     poi2_max += padding * poi2_range
    
#     # Create high-resolution grid
#     poi1_grid = np.linspace(poi1_min, poi1_max, nbins)
#     poi2_grid = np.linspace(poi2_min, poi2_max, nbins)
#     poi1_mesh, poi2_mesh = np.meshgrid(poi1_grid, poi2_grid)
    
#     # Interpolate NLL values onto grid
#     points = np.column_stack((poi1_vals, poi2_vals))
    
#     max_nll = float(np.nanmax(nll_vals))
#     # Try cubic interpolation with a finite fill_value to avoid NaNs
#     try:
#         nll_interpolated = griddata(points, nll_vals, (poi1_mesh, poi2_mesh),
#                                     method='cubic', fill_value=max_nll)
#         # If cubic left NaNs, fill those with linear
#         if np.any(~np.isfinite(nll_interpolated)):
#             nll_linear = griddata(points, nll_vals, (poi1_mesh, poi2_mesh),
#                                   method='linear', fill_value=max_nll)
#             nll_interpolated[~np.isfinite(nll_interpolated)] = nll_linear[~np.isfinite(nll_interpolated)]
#     except Exception:
#         # As a more robust fallback try radial-basis interpolation (RBF) which extrapolates better
#         try:
#             rbf = Rbf(poi1_vals, poi2_vals, nll_vals, function='linear')
#             nll_interpolated = rbf(poi1_mesh, poi2_mesh)
#         except Exception:
#             # Last resort: linear griddata with finite fill
#             nll_interpolated = griddata(points, nll_vals, (poi1_mesh, poi2_mesh),
#                                         method='linear', fill_value=max_nll)

#     # Ensure no NaNs remain (safety) and treat outside-hull as high cost
#     nll_interpolated = np.where(np.isfinite(nll_interpolated), nll_interpolated, max_nll)
#     # Gentle Gaussian smoothing to keep small-scale features (smaller sigma for higher accuracy)
#     nll_interpolated = gaussian_filter(nll_interpolated, sigma=0.8)
     
#     return poi1_grid, poi2_grid, nll_interpolated

# def calculate_correlation_and_uncertainties(poi1_vals, poi2_vals, nll_vals):
#     """Calculate correlation and proper uncertainties from scan data - using 3σ filtering like FitDiagnostics"""
#     poi1_vals = np.array(poi1_vals)
#     poi2_vals = np.array(poi2_vals)
#     nll_vals = np.array(nll_vals)
    
#     print(f">>> Debug correlation calculation:")
#     print(f"    Total scan points: {len(poi1_vals)}")
#     print(f"    POI1 range: [{np.min(poi1_vals):.4f}, {np.max(poi1_vals):.4f}]")
#     print(f"    POI2 range: [{np.min(poi2_vals):.4f}, {np.max(poi2_vals):.4f}]")
#     print(f"    NLL range: [{np.min(nll_vals):.4f}, {np.max(nll_vals):.4f}]")
    
#     # Find the best fit point
#     min_idx = np.argmin(nll_vals)
#     best_poi1 = poi1_vals[min_idx]
#     best_poi2 = poi2_vals[min_idx]
    
#     print(f"    Best fit: POI1={best_poi1:.4f}, POI2={best_poi2:.4f}")
    
#     # **Use same filtering as FitDiagnostics: points within 3σ (deltaNLL < 9)**
#     # Try different thresholds, but prefer 3σ
#     thresholds = [9.0, 10.0, 20.0, 1.0]
#     min_points_required = 25
#     for threshold in thresholds:
#         mask = nll_vals < threshold
#         filtered_poi1 = poi1_vals[mask]
#         filtered_poi2 = poi2_vals[mask]
#         filtered_nll = nll_vals[mask]
        
#         nkept = np.sum(mask)
#         print(f"    Points with deltaNLL < {threshold}: {nkept} / {len(poi1_vals)}")
#         if nkept < min_points_required:
#             print(f"    WARNING: Not enough scan points (<{min_points_required}) at threshold {threshold}, trying next...")
#             continue
        
#         # **Same method as FitDiagnostics: direct correlation coefficient of filtered points**
#         if len(np.unique(filtered_poi1)) > 1 and len(np.unique(filtered_poi2)) > 1:
#             correlation = np.corrcoef(filtered_poi1, filtered_poi2)[0, 1]
#             print(f"    Raw correlation (deltaNLL < {threshold}): {correlation:.4f}")
            
#             if not np.isnan(correlation) and abs(correlation) <= 1.0:
#                 # Calculate uncertainties using likelihood weighting of filtered points
#                 weights = np.exp(-0.5 * filtered_nll)
#                 weights = weights / np.sum(weights)
                
#                 # Weighted mean and std
#                 poi1_mean = np.sum(weights * filtered_poi1)
#                 poi2_mean = np.sum(weights * filtered_poi2)
                
#                 poi1_var = np.sum(weights * (filtered_poi1 - poi1_mean)**2)
#                 poi2_var = np.sum(weights * (filtered_poi2 - poi2_mean)**2)
                
#                 poi1_err = np.sqrt(poi1_var)
#                 poi2_err = np.sqrt(poi2_var)
                
#                 # Ensure reasonable minimum uncertainty
#                 poi1_range = np.max(poi1_vals) - np.min(poi1_vals)
#                 poi2_range = np.max(poi2_vals) - np.min(poi2_vals)
#                 poi1_err = max(poi1_err, poi1_range / 200.0)  # tighter minimum for higher-accuracy scans
#                 poi2_err = max(poi2_err, poi2_range / 200.0)
                
#                 print(f"    Using threshold deltaNLL < {threshold}")
#                 break
#             else:
#                 print(f"    Invalid correlation at threshold {threshold}, trying next...")
#         else:
#             print(f"    Insufficient parameter variation at threshold {threshold}")
    
#     # Final fallback if all thresholds failed
#     if correlation == 0.0 or np.isnan(correlation) or abs(correlation) > 1.0:
#         print("    WARNING: All correlation calculations failed, using fallback method")
        
#         # Use all points as last resort
#         if len(np.unique(poi1_vals)) > 1 and len(np.unique(poi2_vals)) > 1:
#             correlation = np.corrcoef(poi1_vals, poi2_vals)[0, 1]
#             if np.isnan(correlation) or abs(correlation) > 1.0:
#                 correlation = 0.0
#         else:
#             correlation = 0.0
        
#         # Simple uncertainty estimates
#         poi1_err = (np.max(poi1_vals) - np.min(poi1_vals)) / 6.0  # Range/6 ≈ 1σ for normal distribution
#         poi2_err = (np.max(poi2_vals) - np.min(poi2_vals)) / 6.0
    
#     print(f"    Final correlation: {correlation:.4f}")
#     print(f"    Uncertainties: {poi1_err:.4f}, {poi2_err:.4f}")
    
#     return correlation, poi1_err, poi2_err

# def extract_2d_scan_data(multidimfit_file, poi1_name, poi2_name):
#     """
#     Extract 2D scan data from MultiDimFit output file.
#     Returns arrays of parameter values and corresponding NLL values.
#     """
#     print(f">>> Reading 2D scan data from {multidimfit_file}")
    
#     file = ensureTFile(multidimfit_file)
#     tree = file.Get('limit')
#     if not tree:
#         print("ERROR: Could not find 'limit' tree in MultiDimFit file")
#         file.Close()
#         return None
    
#     # First, let's see what branches are available
#     print(">>> Available branches:")
#     branches = []
#     for branch in tree.GetListOfBranches():
#         branch_name = branch.GetName()
#         branches.append(branch_name)
#         print(f"  {branch_name}")
    
#     poi1_vals = []
#     poi2_vals = []
#     nll_vals = []
    
#     nentries = tree.GetEntries()
#     print(f">>> Processing {nentries} entries from MultiDimFit scan")
    
#     # Based on ROOT output, the branches are directly named
#     poi1_branch = poi1_name  # 'tes_DM0'
#     poi2_branch = poi2_name  # 'tid_SF_DM0' 
#     nll_branch = 'deltaNLL'
    
#     # Verify branches exist
#     if poi1_branch not in branches:
#         print(f"ERROR: Branch '{poi1_branch}' not found")
#         file.Close()
#         return None
#     if poi2_branch not in branches:
#         print(f"ERROR: Branch '{poi2_branch}' not found") 
#         file.Close()
#         return None
#     if nll_branch not in branches:
#         print(f"ERROR: Branch '{nll_branch}' not found")
#         file.Close()
#         return None
    
#     print(f">>> Using branches: {poi1_branch}, {poi2_branch}, {nll_branch}")
    
#     for i in range(nentries):
#         tree.GetEntry(i)
        
#         # Get parameter values using the found branch names
#         poi1_val = getattr(tree, poi1_branch)
#         poi2_val = getattr(tree, poi2_branch)
#         nll_val = getattr(tree, nll_branch)
        
#         poi1_vals.append(poi1_val)
#         poi2_vals.append(poi2_val)
#         nll_vals.append(nll_val)
        
#         # Debug first few entries
#         if i < 5:
#             print(f"  Entry {i}: {poi1_branch}={poi1_val:.4f}, {poi2_branch}={poi2_val:.4f}, {nll_branch}={nll_val:.4f}")
    
#     file.Close()
    
#     if len(poi1_vals) == 0:
#         print(f"ERROR: No valid data found")
#         return None
    
#     # Convert to numpy arrays for easier manipulation
#     import numpy as np
#     poi1_vals = np.array(poi1_vals)
#     poi2_vals = np.array(poi2_vals)
#     nll_vals = np.array(nll_vals)
    
#     # Debug scan pattern
#     print(f">>> Scan pattern analysis:")
#     print(f"    Unique POI1 values: {len(np.unique(poi1_vals))}")
#     print(f"    Unique POI2 values: {len(np.unique(poi2_vals))}")
#     unique_poi1 = np.unique(poi1_vals)
#     unique_poi2 = np.unique(poi2_vals)
#     expected_points = len(unique_poi1) * len(unique_poi2)
#     print(f"    Expected grid points: {expected_points}, Actual points: {len(poi1_vals)}")
#     print(f"    POI1 range: [{np.min(poi1_vals):.4f}, {np.max(poi1_vals):.4f}]")
#     print(f"    POI2 range: [{np.min(poi2_vals):.4f}, {np.max(poi2_vals):.4f}]")
#     print(f"    NLL range: [{np.min(nll_vals):.4f}, {np.max(nll_vals):.4f}]")
    
#     # Check if this is actually a 2D scan or two 1D scans
#     if len(poi1_vals) == len(unique_poi1) + len(unique_poi2):
#         print("    WARNING: This looks like two separate 1D scans, not a 2D scan!")
#         print("    You need to run: combine -M MultiDimFit --algo=grid --redefineSignalPOIs poi1,poi2")
#     elif len(poi1_vals) < expected_points * 0.8:
#         print("    WARNING: Scan appears incomplete or not fully gridded")
#     else:
#         print("    This appears to be a proper 2D grid scan")
    
#     # deltaNLL is already the correct quantity (no need to subtract minimum)
#     delta_nll_vals = nll_vals
    
#     # Find best fit point
#     min_idx = np.argmin(delta_nll_vals)
#     best_poi1 = poi1_vals[min_idx]
#     best_poi2 = poi2_vals[min_idx]
    
#     # Calculate correlation and uncertainties
#     correlation, poi1_err, poi2_err = calculate_correlation_and_uncertainties(
#         poi1_vals, poi2_vals, delta_nll_vals)
    
#     print(f">>> Found {len(poi1_vals)} data points")
#     print(f">>> Best fit: {poi1_name} = {best_poi1:.4f}, {poi2_name} = {best_poi2:.4f}")
#     print(f">>> NLL range: {np.min(delta_nll_vals):.3f} to {np.max(delta_nll_vals):.3f}")
#     print(f">>> Calculated correlation: {correlation:.4f}")
#     print(f">>> {poi1_name} = {best_poi1:.4f} ± {poi1_err:.4f}")
#     print(f">>> {poi2_name} = {best_poi2:.4f} ± {poi2_err:.4f}")
    
#     return {
#         'poi1_name': poi1_name,
#         'poi2_name': poi2_name,
#         'poi1_vals': poi1_vals,
#         'poi2_vals': poi2_vals,
#         'delta_nll_vals': delta_nll_vals,
#         'best_poi1': best_poi1,
#         'best_poi2': best_poi2,
#         'poi1_err': poi1_err,
#         'poi2_err': poi2_err,
#         'correlation': correlation
#     }

# def plot_2d_scan(setup, region, year, scan_data, **kwargs):
#     """
#     Minimal plot of 2D -2ΔlnL from MultiDimFit scan (no filtering, no smoothing).
#     """
#     print(f">>> Plotting 2D scan for {region} (simple mode)")
    
#     indir = kwargs.get('indir', f"output_{year}")
#     outdir = indir.replace('output', 'plots')
#     tag = kwargs.get('tag', "")
#     plottag = kwargs.get('plottag', "")
#     poi1_name = scan_data['poi1_name']
#     poi2_name = scan_data['poi2_name']
#     channel = setup["channel"].replace("mu", "m").replace("tau", "t")
    
#     ensureDirectory(outdir)
#     canvasname = f"{outdir}/scan_2D_{poi1_name}_{poi2_name}_{channel}_{region}{tag}{plottag}"
    
#     poi1_vals = scan_data['poi1_vals']
#     poi2_vals = scan_data['poi2_vals']
#     delta_nll_vals = scan_data['delta_nll_vals']
#     best_poi1 = scan_data['best_poi1']
#     best_poi2 = scan_data['best_poi2']
    
#     # Determine plot ranges with small margin
#     poi1_min, poi1_max = np.min(poi1_vals), np.max(poi1_vals)
#     poi2_min, poi2_max = np.min(poi2_vals), np.max(poi2_vals)
#     margin = 0.05
#     poi1_range = poi1_max - poi1_min if poi1_max > poi1_min else 1.0
#     poi2_range = poi2_max - poi2_min if poi2_max > poi2_min else 1.0
#     poi1_min -= margin * poi1_range
#     poi1_max += margin * poi1_range
#     poi2_min -= margin * poi2_range
#     poi2_max += margin * poi2_range
    
#     # Histogram resolution (choose grid matching unique values if grid-like)
#     ux = np.unique(np.sort(poi1_vals))
#     uy = np.unique(np.sort(poi2_vals))
#     nbins_x = len(ux) if len(ux) > 1 and len(ux) <= 500 else 200
#     nbins_y = len(uy) if len(uy) > 1 and len(uy) <= 500 else 200
    
#     hist_2d = TH2D("hist_2d", "", nbins_x, poi1_min, poi1_max, nbins_y, poi2_min, poi2_max)
#     hist_2d.GetXaxis().SetTitle(poi1_name)
#     hist_2d.GetYaxis().SetTitle(poi2_name)
#     hist_2d.GetZaxis().SetTitle("-2#Deltaln(L)")
    
#     # Fill histogram with the raw scan values (take minimum if multiple points fall in same bin)
#     for x, y, z in zip(poi1_vals, poi2_vals, delta_nll_vals):
#         bx = hist_2d.GetXaxis().FindBin(x)
#         by = hist_2d.GetYaxis().FindBin(y)
#         if bx < 1 or bx > hist_2d.GetNbinsX() or by < 1 or by > hist_2d.GetNbinsY():
#             continue
#         cur = hist_2d.GetBinContent(bx, by)
#         if cur == 0 or z < cur:
#             hist_2d.SetBinContent(bx, by, z)
    
#     # If there are empty bins, leave them as-is (no smoothing/interpolation)
#     hist_2d.SetMinimum(0)
#     max_val = hist_2d.GetMaximum()
#     if max_val > 20:
#         hist_2d.SetMaximum(20)
    
#     # Canvas and draw
#     c = TCanvas('c_simple', 'c_simple', 800, 700)
#     c.SetTicks(1,1)
#     c.SetBottomMargin(0.12)
#     c.SetLeftMargin(0.12)
#     hist_2d.Draw("COLZ")
    
#     # Draw standard 2D contour levels for -2ΔlnL
#     contour_levels = [2.30, 6.18, 11.83]
#     hist_contour = hist_2d.Clone("hist_contour_simple")
#     hist_contour.SetContour(len(contour_levels))
#     for i, lvl in enumerate(contour_levels):
#         hist_contour.SetContourLevel(i, lvl)
#     hist_contour.SetLineColor(kBlack)
#     hist_contour.SetLineWidth(2)
#     hist_contour.Draw("CONT3 SAME")
    
#     # Overlay raw scan points and best-fit marker
#     gp = ROOT.TGraph(len(poi1_vals), poi1_vals, poi2_vals)
#     gp.SetMarkerStyle(20); gp.SetMarkerSize(0.6); gp.SetMarkerColor(kBlue)
#     gp.Draw("P SAME")
#     best = ROOT.TMarker(best_poi1, best_poi2, 29)
#     best.SetMarkerColor(kRed); best.SetMarkerSize(1.5)
#     best.Draw("SAME")
    
#     # Simple text
#     lat = TLatex()
#     lat.SetNDC(True); lat.SetTextSize(0.035)
#     lat.DrawLatex(0.15, 0.92, f"Region: {region}")
#     lat.DrawLatex(0.15, 0.86, f"Best: {poi1_name}={best_poi1:.4f}, {poi2_name}={best_poi2:.4f}")
    
#     c.Update()
#     c.SaveAs(canvasname + ".png")
#     c.SaveAs(canvasname + ".pdf")
#     c.SaveAs(canvasname + ".root")
#     print(f">>> Saved simple 2D scan plot: {canvasname}.png")
#     c.Close()

# def plot_measurement_summary(region_labels, measurements, **kwargs):
#     """Create a summary plot showing measurements with error bars"""
#     title = kwargs.get('title', "Measurements")
#     ylabel = kwargs.get('ylabel', "value")
#     outname = kwargs.get('outname', "measurements")
#     year = kwargs.get('year', "2024")
    
#     n_regions = len(region_labels)
#     if n_regions == 0:
#         print("No measurements to plot")
#         return
    
#     # Create canvas
#     canvas_height = max(600, 60 + 40*n_regions)
#     canvas_width = 800
#     canvas = ROOT.TCanvas('canvas_summary', 'canvas_summary', 100, 100, canvas_width, canvas_height)
#     canvas.SetFillColor(0)
#     canvas.SetBorderMode(0)
#     canvas.SetFrameFillStyle(0)
#     canvas.SetFrameBorderMode(0)
    
#     # Set margins
#     top_margin = 0.08
#     bottom_margin = 0.12
#     left_margin = 0.25
#     right_margin = 0.05
    
#     canvas.SetTopMargin(top_margin)
#     canvas.SetBottomMargin(bottom_margin) 
#     canvas.SetLeftMargin(left_margin)
#     canvas.SetRightMargin(right_margin)
#     canvas.SetGrid(1, 0)
#     canvas.cd()
    
#     # Determine x-axis range
#     values = [m[0] for m in measurements]
#     errors_down = [m[1] for m in measurements] 
#     errors_up = [m[2] for m in measurements]
    
#     x_min = min([v - e for v, e in zip(values, errors_down)])
#     x_max = max([v + e for v, e in zip(values, errors_up)])
#     x_range = x_max - x_min
#     x_margin = 0.15 * x_range
#     x_min -= x_margin
#     x_max += x_margin
    
#     # Create frame
#     frame = canvas.DrawFrame(x_min, 0.0, x_max, float(n_regions))
#     frame.GetYaxis().SetLabelSize(0.0)
#     frame.GetXaxis().SetLabelSize(0.042)
#     frame.GetXaxis().SetTitleSize(0.050) 
#     frame.GetXaxis().SetTitleOffset(1.1)
#     frame.GetYaxis().SetNdivisions(n_regions, 0, 0, False)
#     frame.GetXaxis().SetTitle(ylabel)
#     frame.GetXaxis().SetNdivisions(510)
    
#     # Create graph with error bars
#     graph = ROOT.TGraphAsymmErrors(n_regions)
    
#     for i, (region, measurement) in enumerate(zip(region_labels, measurements)):
#         y_pos = n_regions - i - 0.5
#         val, err_down, err_up = measurement
#         graph.SetPoint(i, val, y_pos)
#         graph.SetPointError(i, err_down, err_up, 0.1, 0.1)
    
#     # Style the graph
#     graph.SetMarkerStyle(20)
#     graph.SetMarkerSize(1.0)
#     graph.SetMarkerColor(ROOT.kBlack)
#     graph.SetLineColor(ROOT.kBlack)
#     graph.SetLineWidth(2)
    
#     # Draw the graph
#     graph.Draw("PE SAME")
    
#     # Add vertical line at 1.0 if appropriate
#     if min(values) < 1.0 < max(values):
#         line = ROOT.TLine(1.0, 0.0, 1.0, float(n_regions))
#         line.SetLineStyle(2)
#         line.SetLineColor(ROOT.kGray+2)
#         line.Draw("SAME")
    
#     # Add region labels
#     latex = ROOT.TLatex()
#     latex.SetTextSize(0.035)
#     latex.SetTextFont(42)
#     latex.SetTextAlign(32)
#     latex.SetNDC(True)
    
#     for i, region in enumerate(region_labels):
#         y_pos_ndc = 1.0 - top_margin - (i + 0.5) * (1.0 - top_margin - bottom_margin) / n_regions
#         latex.DrawLatex(left_margin - 0.02, y_pos_ndc, region)
    
#     # Add CMS header
#     cms_latex = ROOT.TLatex()
#     cms_latex.SetTextSize(0.060)
#     cms_latex.SetTextFont(61)
#     cms_latex.SetTextAlign(11)
#     cms_latex.SetNDC(True)
#     cms_latex.DrawLatex(left_margin, 1.0 - top_margin + 0.01, "CMS")
    
#     # Add "Internal" label
#     internal_latex = ROOT.TLatex()
#     internal_latex.SetTextSize(0.045)
#     internal_latex.SetTextFont(52)
#     internal_latex.SetTextAlign(11)
#     internal_latex.SetNDC(True) 
#     internal_latex.DrawLatex(left_margin + 0.12, 1.0 - top_margin + 0.01, "Internal")
    
#     # Add year and energy
#     year_latex = ROOT.TLatex()
#     year_latex.SetTextSize(0.045)
#     year_latex.SetTextFont(42)
#     year_latex.SetTextAlign(31)
#     year_latex.SetNDC(True)
#     year_latex.DrawLatex(1.0 - right_margin, 1.0 - top_margin + 0.01, f"{year}, 109 fb^{{-1}} (13.6 TeV)")
    
#     # Add title
#     if title:
#         title_latex = ROOT.TLatex()
#         title_latex.SetTextSize(0.045)
#         title_latex.SetTextFont(42)
#         title_latex.SetTextAlign(11)
#         title_latex.SetNDC(True)
#         title_latex.DrawLatex(left_margin, 1.0 - top_margin - 0.05, title)
    
#     canvas.SetTicks(1, 1)
#     canvas.Modified()
#     canvas.Update()
    
#     # Save
#     canvas.SaveAs(outname + ".png")
#     canvas.SaveAs(outname + ".pdf") 
#     canvas.SaveAs(outname + ".root")
#     print(f">>> Saved measurement summary: {outname}.png")
    
#     canvas.Close()

# def plot_scan_correlations(scan_results_all_regions, **kwargs):
#     """Create a correlation plot using correlations from scan results"""
#     year = kwargs.get('year', '2024')
#     indir = kwargs.get('indir', f"output_{year}")
#     outdir = indir.replace('output', 'plots')
#     tag = kwargs.get('tag', "")
#     plottag = kwargs.get('plottag', "")
#     outname = f"{outdir}/scan_correlations_multidimfit{tag}{plottag}"
    
#     ensureDirectory(outdir)
    
#     # Extract correlations and region names directly from scan results
#     correlations = []
#     region_labels = []
    
#     # Sort regions
#     sorted_regions = sorted(scan_results_all_regions.items(), key=lambda x: format_region_for_sorting(x[0]))
    
#     for region, scan_data in sorted_regions:
#         if scan_data is None:
#             continue
#         correlation = scan_data.get('correlation', 0.0)
#         correlations.append(correlation)
#         region_labels.append(format_region_label(region))
    
#     n_regions = len(region_labels)
#     if n_regions == 0:
#         print("No correlation data to plot")
#         return
    
#     # Create canvas
#     canvas_height = max(600, 60 + 40*n_regions)
#     canvas_width = 800
#     canvas = ROOT.TCanvas('canvas_scan_corr', 'canvas_scan_corr', 100, 100, canvas_width, canvas_height)
#     canvas.SetFillColor(0)
#     canvas.SetBorderMode(0)
#     canvas.SetFrameFillStyle(0)
#     canvas.SetFrameBorderMode(0)
    
#     # Set margins
#     top_margin = 0.08
#     bottom_margin = 0.12
#     left_margin = 0.25
#     right_margin = 0.05
    
#     canvas.SetTopMargin(top_margin)
#     canvas.SetBottomMargin(bottom_margin) 
#     canvas.SetLeftMargin(left_margin)
#     canvas.SetRightMargin(right_margin)
#     canvas.SetGrid(1, 0)
#     canvas.cd()
    
#     # Determine x-axis range for correlations (-1 to +1)
#     x_min = -1.2
#     x_max = 1.2
    
#     # Create frame
#     frame = canvas.DrawFrame(x_min, 0.0, x_max, float(n_regions))
#     frame.GetYaxis().SetLabelSize(0.0)
#     frame.GetXaxis().SetLabelSize(0.042)
#     frame.GetXaxis().SetTitleSize(0.050) 
#     frame.GetXaxis().SetTitleOffset(1.1)
#     frame.GetYaxis().SetNdivisions(n_regions, 0, 0, False)
#     frame.GetXaxis().SetTitle("TES-TauID Correlation from 2D Scans")
#     frame.GetXaxis().SetNdivisions(510)
    
#     # Create individual markers for each correlation
#     markers = []
    
#     for i, (region, correlation) in enumerate(zip(region_labels, correlations)):
#         y_pos = n_regions - i - 0.5
        
#         # Simple color scheme: blue for all correlations
#         color = ROOT.kBlue + 2
        
#         # Create marker
#         marker = ROOT.TMarker(correlation, y_pos, 20)
#         marker.SetMarkerSize(1.2)
#         marker.SetMarkerColor(color)
#         markers.append(marker)
#         marker.Draw("SAME")
    
#     # Add vertical line at 0.0 (no correlation)
#     line_zero = ROOT.TLine(0.0, 0.0, 0.0, float(n_regions))
#     line_zero.SetLineStyle(2)
#     line_zero.SetLineColor(ROOT.kGray+2)
#     line_zero.SetLineWidth(2)
#     line_zero.Draw("SAME")
    
#     # Add region labels
#     latex = ROOT.TLatex()
#     latex.SetTextSize(0.035)
#     latex.SetTextFont(42)
#     latex.SetTextAlign(32)
#     latex.SetNDC(True)
    
#     for i, region in enumerate(region_labels):
#         y_pos_ndc = 1.0 - top_margin - (i + 0.5) * (1.0 - top_margin - bottom_margin) / n_regions
#         latex.DrawLatex(left_margin - 0.02, y_pos_ndc, region)
    
#     # Add correlation values next to points
#     corr_latex = ROOT.TLatex()
#     corr_latex.SetTextSize(0.030)
#     corr_latex.SetTextFont(42)
#     corr_latex.SetTextAlign(11)
    
#     for i, correlation in enumerate(correlations):
#         y_pos = n_regions - i - 0.5
#         x_pos = correlation + 0.05 if correlation >= 0 else correlation - 0.05
#         corr_latex.DrawLatex(x_pos, y_pos, f"{correlation:.3f}")
    
#     # Add CMS header
#     cms_latex = ROOT.TLatex()
#     cms_latex.SetTextSize(0.060)
#     cms_latex.SetTextFont(61)
#     cms_latex.SetTextAlign(11)
#     cms_latex.SetNDC(True)
#     cms_latex.DrawLatex(left_margin, 1.0 - top_margin + 0.01, "CMS")
    
#     # Add "Internal" label
#     internal_latex = ROOT.TLatex()
#     internal_latex.SetTextSize(0.045)
#     internal_latex.SetTextFont(52)
#     internal_latex.SetTextAlign(11)
#     internal_latex.SetNDC(True) 
#     internal_latex.DrawLatex(left_margin + 0.12, 1.0 - top_margin + 0.01, "Internal")
    
#     # Add year and energy
#     year_latex = ROOT.TLatex()
#     year_latex.SetTextSize(0.045)
#     year_latex.SetTextFont(42)
#     year_latex.SetTextAlign(31)
#     year_latex.SetNDC(True)
#     year_latex.DrawLatex(1.0 - right_margin, 1.0 - top_margin + 0.01, f"{year}, 109 fb^{{-1}} (13.6 TeV)")
    
#     # Add title
#     title_latex = ROOT.TLatex()
#     title_latex.SetTextSize(0.045)
#     title_latex.SetTextFont(42)
#     title_latex.SetTextAlign(11)
#     title_latex.SetNDC(True)
#     title_latex.DrawLatex(left_margin, 1.0 - top_margin - 0.05, "TES-TauID Correlations from MultiDimFit Scans")
    
#     canvas.SetTicks(1, 1)
#     canvas.Modified()
#     canvas.Update()
    
#     # Save
#     canvas.SaveAs(outname + ".png")
#     canvas.SaveAs(outname + ".pdf") 
#     canvas.SaveAs(outname + ".root")
#     print(f">>> Saved scan correlation plot: {outname}.png")
    
#     canvas.Close()

# def plot_summary_from_multiple_regions(setup, scan_results_all_regions, **kwargs):
#     """Create summary plots showing TES and TauID measurements with uncertainties"""
#     print(">>> Creating summary plots from MultiDimFit results")
    
#     year = kwargs.get('year', '2024')
#     indir = kwargs.get('indir', f"output_{year}")
#     outdir = indir.replace('output', 'plots')
#     tag = kwargs.get('tag', "")
#     plottag = kwargs.get('plottag', "")
    
#     ensureDirectory(outdir)
    
#     # Extract measurements
#     tes_measurements = []
#     tid_measurements = []
#     region_labels = []
    
#     # Sort regions
#     sorted_regions = sorted(scan_results_all_regions.items(), key=lambda x: format_region_for_sorting(x[0]))
    
#     for region, scan_data in sorted_regions:
#         if scan_data is None:
#             continue
            
#         region_label = format_region_label(region)
#         region_labels.append(region_label)
        
#         poi1_name = scan_data['poi1_name']
#         poi2_name = scan_data['poi2_name']
        
#         if 'tes' in poi1_name.lower():
#             tes_val = scan_data['best_poi1']
#             tes_err = scan_data['poi1_err']
#             tid_val = scan_data['best_poi2'] 
#             tid_err = scan_data['poi2_err']
#         else:
#             tes_val = scan_data['best_poi2']
#             tes_err = scan_data['poi2_err']
#             tid_val = scan_data['best_poi1']
#             tid_err = scan_data['poi1_err']
            
#         tes_measurements.append((tes_val, tes_err, tes_err))
#         tid_measurements.append((tid_val, tid_err, tid_err))
    
#     # Create plots
#     plot_measurement_summary(
#         region_labels, tes_measurements,
#         title="Tau Energy Scale",
#         ylabel="tau energy scale", 
#         outname=f"{outdir}/tes_summary_multidimfit{tag}{plottag}",
#         year=year
#     )
    
#     plot_measurement_summary(
#         region_labels, tid_measurements,
#         title="Tau ID Scale Factor", 
#         ylabel="tau ID scale factor",
#         outname=f"{outdir}/tid_summary_multidimfit{tag}{plottag}",
#         year=year
#     )
    
#     # Create correlation plot
#     plot_scan_correlations(
#         scan_results_all_regions,
#         year=year, 
#         indir=indir, 
#         tag=tag, 
#         plottag=plottag
#     )
    
#     print(f">>> Created TES summary plot: {outdir}/tes_summary_multidimfit{tag}{plottag}.png")
#     print(f">>> Created TauID summary plot: {outdir}/tid_summary_multidimfit{tag}{plottag}.png")
#     print(f">>> Created scan correlation plot: {outdir}/scan_correlations_multidimfit{tag}{plottag}.png")

# def write_2d_fit_results(poi1_name, poi2_name, poi1_val, poi1_err_down, poi1_err_up, 
#                          poi2_val, poi2_err_down, poi2_err_up, correlation, region, **kwargs):
#     """Write 2D fit results to text file, similar to measurepoi() in plotParabola_POI_region.py"""
#     year = kwargs.get('year', '2024')
#     tag = kwargs.get('tag', '')
#     channel = kwargs.get('channel', 'mt')
#     outdir = kwargs.get('outdir', 'plots')
    
#     ensureDirectory(outdir)
    
#     # Create output filename similar to plotParabola_POI_region format
#     outfname = f"{outdir}/measurement_2D_{poi1_name}_{poi2_name}_{channel}_{region}{tag}.txt"
    
#     print(f">>> Writing 2D fit results to {outfname}")
    
#     # Write results to file in same format as plotParabola_POI_region
#     with open(outfname, 'w') as file:
#         file.write("# 2D Fit Results from MultiDimFit\n")
#         file.write(f"# Region: {region}\n")
#         file.write(f"# Year: {year}\n")
#         file.write(f"# Channel: {channel}\n")
#         file.write("# Format: parameter value error_down error_up\n")
#         file.write(f"{poi1_name} {poi1_val:.6f} {poi1_err_down:.6f} {poi1_err_up:.6f}\n")
#         file.write(f"{poi2_name} {poi2_val:.6f} {poi2_err_down:.6f} {poi2_err_up:.6f}\n")
#         file.write(f"correlation {correlation:.6f}\n")
    
#     return outfname

# def write_summary_results(scan_results_all_regions, **kwargs):
#     """Write summary of all 2D fit results to a single file"""
#     year = kwargs.get('year', '2024')
#     tag = kwargs.get('tag', '')
#     channel = kwargs.get('channel', 'mt')
#     outdir = kwargs.get('outdir', 'plots')
    
#     ensureDirectory(outdir)
    
#     # Create summary output filename
#     outfname = f"{outdir}/summary_2D_results_{channel}{tag}_{year}.txt"
    
#     print(f">>> Writing summary 2D fit results to {outfname}")
    
#     with open(outfname, 'w') as file:
#         file.write("# Summary of 2D Fit Results from MultiDimFit\n")
#         file.write(f"# Year: {year}\n")
#         file.write(f"# Channel: {channel}\n")
#         file.write("# Format: region poi1_name poi1_val poi1_err poi2_name poi2_val poi2_err correlation\n")
        
#         # Sort regions for consistent output
#         sorted_regions = sorted(scan_results_all_regions.items(), key=lambda x: format_region_for_sorting(x[0]))
        
#         for region, scan_data in sorted_regions:
#             if scan_data is None:
#                 continue
                
#             # Get parameter info
#             poi1_name = scan_data['poi1_name']
#             poi2_name = scan_data['poi2_name']
#             poi1_val = scan_data['best_poi1']
#             poi2_val = scan_data['best_poi2']
#             poi1_err = scan_data['poi1_err']
#             poi2_err = scan_data['poi2_err']
#             correlation = scan_data.get('correlation', 0.0)
            
#             # Write to file
#             file.write(f"{region} {poi1_name} {poi1_val:.6f} {poi1_err:.6f} ")
#             file.write(f"{poi2_name} {poi2_val:.6f} {poi2_err:.6f} ")
#             file.write(f"{correlation:.6f}\n")
    
#     return outfname

# def main(args):
#     """Main function - handle multiple regions and create summary plots"""
    
#     print("Using configuration file: %s" % args.config)
#     with open(args.config, 'r') as file:
#         setup = yaml.safe_load(file)
    
#     channel = setup["channel"].replace("mu", "m").replace("tau", "t")
#     tag = setup.get("tag", "")
#     era = args.year
#     extratag = "_DeepTau"
    
#     # Input directory
#     if args.indir:
#         if not args.indir.rstrip('/').endswith(str(era)):
#             indir = os.path.join(args.indir, str(era))
#         else:
#             indir = args.indir
#     else:
#         indir = f"output_{era}"
    
#     # Process regions
#     scan_results_all_regions = {}
    
#     if args.poi1 and args.poi2:
#         # Single region mode
#         poi1_name = args.poi1  # e.g., "tes_DM0"
#         poi2_name = args.poi2  # e.g., "tid_SF_DM0"
        
#         # Extract region from POI name (assuming format like "tes_DM0")
#         if '_' in poi1_name:
#             region = poi1_name.split('_', 1)[1]
#         else:
#             region = args.region if args.region else "DM0"
            
#         regions_to_process = [region]
#     else:
#         # Multiple regions mode - get from config
#         try:
#             regions_to_process = setup["observables"]["m_vis"]["scanRegions"]
#             print(f">>> Processing {len(regions_to_process)} regions from config: {regions_to_process}")
#         except KeyError:
#             print("ERROR: No regions found in config file. Please specify --poi1 and --poi2 for single region mode.")
#             sys.exit(1)
    
#     # Process each region
#     for region in regions_to_process:
#         print(f"\n>>> Processing region: {region}")
        
#         # Construct POI names if not provided
#         if not args.poi1 or not args.poi2:
#             poi1_name = f"tes_{region}"
#             poi2_name = f"tid_SF_{region}"
#         else:
#             poi1_name = args.poi1
#             poi2_name = args.poi2
        
#         # Construct MultiDimFit filename
#         multidimfit_filename = f"{indir}/higgsCombine.{channel}_m_vis-{region}{tag}{extratag}-{era}-13TeV.MultiDimFit.mH90.root"
        
#         if not os.path.exists(multidimfit_filename):
#             # Try alternative naming
#             multidimfit_filename = f"{indir}/higgsCombine.mt_m_vis-{region}{tag}{extratag}-{era}-13TeV.MultiDimFit.mH90.root"
#             if not os.path.exists(multidimfit_filename):
#                 print(f"WARNING: MultiDimFit file not found for {region}")
#                 scan_results_all_regions[region] = None
#                 continue
        
#         print(f">>> Using: {multidimfit_filename}")
        
#         # Extract scan data
#         scan_data = extract_2d_scan_data(multidimfit_filename, poi1_name, poi2_name)
#         if scan_data is None:
#             print(f"ERROR: Could not extract scan data for {region}")
#             scan_results_all_regions[region] = None
#             continue
            
#         scan_results_all_regions[region] = scan_data
        
#         # Create individual 2D plot
#         plot_2d_scan(setup, region, era, scan_data, 
#                      indir=indir, tag=tag, plottag=args.plottag)
        
#         # Write individual text results for this region
#         outdir = indir.replace('output', 'plots')
#         poi1_val = scan_data['best_poi1']
#         poi1_err = scan_data['poi1_err']
#         poi2_val = scan_data['best_poi2'] 
#         poi2_err = scan_data['poi2_err']
#         correlation = scan_data.get('correlation', 0.0)
        
#         # Write individual region results
#         write_2d_fit_results(
#             poi1_name, poi2_name, poi1_val, poi1_err, poi1_err,  # symmetric errors
#             poi2_val, poi2_err, poi2_err, correlation, region,
#             year=era, tag=tag, channel=channel, outdir=outdir
#         )
    
#     # Create summary plots and text files if we have multiple regions
#     if len([r for r in scan_results_all_regions.values() if r is not None]) > 1:
#         plot_summary_from_multiple_regions(setup, scan_results_all_regions,
#                                          year=era, indir=indir, tag=tag, 
#                                          plottag=args.plottag)
        
#         # Write summary text results
#         outdir = indir.replace('output', 'plots')
#         write_summary_results(scan_results_all_regions,
#                              year=era, tag=tag, channel=channel, outdir=outdir)
#     elif len(scan_results_all_regions) == 1:
#         print(">>> Only one region processed, summary plots not created")
    
#     print(">>> All plots and text files completed successfully!")

# if __name__ == '__main__':
#     description = '''Plot 2D parabolas from MultiDimFit scan output.'''
#     parser = ArgumentParser(prog="plot2DScan_MultiDimFit", 
#                           description=description, epilog="Success!")
    
#     parser.add_argument('-y', '--year', dest='year', 
#                        choices=['2024', '2016', '2017', '2018', 'UL2016_preVFP', 
#                                'UL2016_postVFP', 'UL2017', 'UL2018', 'UL2018_v10',
#                                '2022_postEE', '2022_preEE', '2023C', '2023D'], 
#                        type=str, default='2024', action='store', 
#                        help="select year")
    
#     parser.add_argument('-c', '--config', dest='config', type=str, 
#                        default='TauES_ID/config/config_coarse_TT.yml', 
#                        action='store', 
#                        help="set config file containing sample & fit setup")
    
#     parser.add_argument('--poi1', dest='poi1', type=str, required=False,
#                        help="first parameter of interest (e.g., tes_DM0). If not provided, will process all regions from config.")
    
#     parser.add_argument('--poi2', dest='poi2', type=str, required=False,
#                        help="second parameter of interest (e.g., tid_SF_DM0). If not provided, will process all regions from config.")
    
#     parser.add_argument('-r', '--region', dest='region', type=str,
#                        help="region name (if not extractable from POI names)")
    
#     parser.add_argument('-i', '--indir', dest='indir', type=str, 
#                        help='input directory')
    
#     parser.add_argument('-t', '--plottag', dest='plottag', type=str, 
#                        default="", help='extra tag for plot filename')
    
#     parser.add_argument('-v', '--verbose', dest='verbose', 
#                        default=False, action='store_true', help="set verbose")
    
#     args = parser.parse_args()
#     main(args)

