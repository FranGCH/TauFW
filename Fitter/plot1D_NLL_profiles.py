#! /usr/bin/env python3
"""Project the 2D MultiDimFit grid onto each POI axis as a profile NLL.

Two output modes:
  * per-fit:   one PNG per 2D fit, two side-by-side panels (TES, TauID SF).
  * overlay:   one PNG per (WP combo, DM), overlaying the 3 pT bins' TES profiles
               to check pT-bin consistency before correlating TES across pT.

The overlay mode also writes a CSV summary with best-fit ± 1σ per pT bin and
pairwise σ-discrepancies between pT bins (useful to decide if correlating TES
across pT is safe)."""
import os, glob, argparse, csv
from array import array
from collections import defaultdict
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def profile_nll(tree, poi):
  """For each unique value of `poi`, take min(deltaNLL) over all other-POI points."""
  prof = defaultdict(lambda: float('inf'))
  n = tree.GetEntries()
  for i in range(n):
    tree.GetEntry(i)
    nll = getattr(tree, 'deltaNLL')
    if nll < 0:
      continue
    key = round(getattr(tree, poi), 6)
    if nll < prof[key]:
      prof[key] = nll
  return sorted(prof.items())


def best_fit_with_errors(xs, ys):
  """From a 1D 2ΔNLL curve `ys` over `xs` (already with ymin=0), return
     (bestfit, err_down, err_up) by linear interpolation of the 2ΔNLL=1 crossings.
     err_down/err_up are positive offsets from bestfit; None if no crossing in range."""
  imin = ys.index(min(ys))
  xbf = xs[imin]
  def cross(rng):
    prev = None
    for i in rng:
      if prev is None or i == imin:
        prev = i
        continue
      if (ys[prev] - 1.0) * (ys[i] - 1.0) <= 0 and ys[i] != ys[prev]:
        t = (1.0 - ys[prev]) / (ys[i] - ys[prev])
        return xs[prev] + t * (xs[i] - xs[prev])
      prev = i
    return None
  x_lo = cross(range(imin, -1, -1))
  x_hi = cross(range(imin, len(xs)))
  err_dn = (xbf - x_lo) if x_lo is not None else None
  err_up = (x_hi - xbf) if x_hi is not None else None
  return xbf, err_dn, err_up


def discrepancy_sigma(bf_a, err_a, bf_b, err_b):
  """Two-sided σ separation: |bf_a-bf_b| / sqrt(σ_a^2 + σ_b^2) using whichever
     side of each error bar points toward the other point."""
  if None in (err_a[0], err_a[1], err_b[0], err_b[1]):
    return None
  err_a_dir = err_a[1] if bf_b > bf_a else err_a[0]
  err_b_dir = err_b[0] if bf_b > bf_a else err_b[1]
  denom = (err_a_dir**2 + err_b_dir**2) ** 0.5
  if denom == 0:
    return None
  return abs(bf_a - bf_b) / denom


# ---------- Per-fit two-panel plot ----------

def draw_profile(pad, xs, ys, poi, title):
  pad.cd()
  n = len(xs)
  g = ROOT.TGraph(n, array('d', xs), array('d', ys))
  g.SetTitle(f"{title};{poi};2#DeltaNLL")
  g.SetLineColor(ROOT.kBlue+1)
  g.SetLineWidth(2)
  g.SetMarkerStyle(20)
  g.SetMarkerSize(0.6)
  g.GetYaxis().SetRangeUser(0, 10.0)
  g.Draw('ALP')
  for lvl, color, lbl in [(1.0, ROOT.kRed, '1#sigma'),
                          (4.0, ROOT.kOrange+1, '2#sigma')]:
    line = ROOT.TLine(xs[0], lvl, xs[-1], lvl)
    line.SetLineColor(color); line.SetLineStyle(2); line.SetLineWidth(2)
    line.DrawClone()
    txt = ROOT.TLatex(xs[-1], lvl, f"  {lbl}")
    txt.SetTextColor(color); txt.SetTextSize(0.035); txt.SetTextAlign(12)
    txt.DrawClone()
  imin = ys.index(min(ys))
  m = ROOT.TMarker(xs[imin], ys[imin], 29)
  m.SetMarkerColor(ROOT.kBlack); m.SetMarkerSize(1.6)
  m.DrawClone()
  pad.SetGridx(); pad.SetGridy()
  return g


def plot_one(filepath, outdir, jet_wp, ele_wp):
  f = ROOT.TFile.Open(filepath)
  t = f.Get('limit')
  if not t:
    f.Close()
    return None
  brs = [b.GetName() for b in t.GetListOfBranches()]
  pois = [b for b in brs if b.startswith('tes_') or b.startswith('tid_SF_')]
  if len(pois) != 2 or 'deltaNLL' not in brs:
    f.Close()
    return None
  pois.sort(key=lambda s: 0 if s.startswith('tes_') else 1)

  base = os.path.basename(filepath)
  label = base.split('mt_m_vis-')[1].split('_mutau')[0]
  title = f"{label}  (VSjet={jet_wp}, VSe={ele_wp})"

  c = ROOT.TCanvas('c', '', 1400, 600)
  c.Divide(2, 1)
  keep = []
  for ipad, poi in enumerate(pois, start=1):
    prof = profile_nll(t, poi)
    if not prof:
      continue
    xs = [p[0] for p in prof]
    ys = [2*p[1] for p in prof]
    ymin = min(ys); ys = [y - ymin for y in ys]
    keep.append(draw_profile(c.cd(ipad), xs, ys, poi, title))

  os.makedirs(outdir, exist_ok=True)
  out = os.path.join(outdir, f"nll1d_{label}.png")
  c.SaveAs(out)
  f.Close()
  return out


# ---------- TES overlay per DM (consistency check) ----------

def collect_profiles(indir, year, jet_wp, ele_wp, dm):
  """For (jet, ele, DM), return {'tes': [...], 'tid': [...]} where each list has
     one profile dict per pt bin: {pt_idx, poi, xs, ys, bf, err_dn, err_up}."""
  out = {'tes': [], 'tid': []}
  for pt_idx in (1, 2, 3):
    fname = (f"higgsCombine.mt_m_vis-{dm}_pt{pt_idx}_mutau_DeepTau-{year}"
             f"-13TeV.MultiDimFit.mH90.root")
    fpath = os.path.join(indir, f'againstjet_{jet_wp}',
                         f'againstelectron_{ele_wp}', year, fname)
    if not os.path.exists(fpath):
      continue
    f = ROOT.TFile.Open(fpath)
    t = f.Get('limit')
    if not t:
      f.Close()
      continue
    brs = [b.GetName() for b in t.GetListOfBranches()]
    pois_by_kind = {
      'tes': [b for b in brs if b.startswith('tes_')],
      'tid': [b for b in brs if b.startswith('tid_SF_')],
    }
    for kind, pois in pois_by_kind.items():
      if not pois:
        continue
      prof = profile_nll(t, pois[0])
      if not prof:
        continue
      xs = [p[0] for p in prof]
      ys = [2*p[1] for p in prof]
      ymin = min(ys); ys = [y - ymin for y in ys]
      bf, dn, up = best_fit_with_errors(xs, ys)
      out[kind].append({'pt_idx': pt_idx, 'poi': pois[0], 'xs': xs, 'ys': ys,
                        'bf': bf, 'err_dn': dn, 'err_up': up})
    f.Close()
  return out


COLORS = [ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kRed+1]


def draw_overlay_panel(pad, profiles, kind_label, axis_label, dm, jet_wp, ele_wp):
  pad.cd()
  pad.SetGridx(); pad.SetGridy()
  leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030)
  keep = []
  xmins, xmaxs = [], []
  for i, p in enumerate(profiles):
    n = len(p['xs'])
    g = ROOT.TGraph(n, array('d', p['xs']), array('d', p['ys']))
    g.SetLineColor(COLORS[p['pt_idx']-1])
    g.SetMarkerColor(COLORS[p['pt_idx']-1])
    g.SetLineWidth(2)
    g.SetMarkerStyle(20); g.SetMarkerSize(0.5)
    if i == 0:
      g.SetTitle(f"{dm}  {kind_label} profile per pT bin"
                 f"  (VSjet={jet_wp}, VSe={ele_wp});{axis_label};2#DeltaNLL")
      g.GetYaxis().SetRangeUser(0, 10.0)
      g.Draw('ALP')
    else:
      g.Draw('LP SAME')
    keep.append(g)
    xmins.append(p['xs'][0]); xmaxs.append(p['xs'][-1])
    dn = p['err_dn'] if p['err_dn'] is not None else float('nan')
    up = p['err_up'] if p['err_up'] is not None else float('nan')
    leg.AddEntry(g, f"pt{p['pt_idx']}: {p['bf']:.4f}_{{-{dn:.4f}}}^{{+{up:.4f}}}", 'LP')
  if not profiles:
    return keep, leg
  xmin = min(xmins); xmax = max(xmaxs)
  for lvl, color, lbl in [(1.0, ROOT.kBlack, '1#sigma'),
                          (4.0, ROOT.kGray+2, '2#sigma')]:
    line = ROOT.TLine(xmin, lvl, xmax, lvl)
    line.SetLineStyle(2); line.SetLineColor(color); line.SetLineWidth(2)
    line.DrawClone()
    txt = ROOT.TLatex(xmax, lvl, f"  {lbl}")
    txt.SetTextColor(color); txt.SetTextSize(0.03); txt.SetTextAlign(12)
    txt.DrawClone()
  leg.Draw()
  keep.append(leg)
  return keep


def plot_overlay(prof_dict, outdir, jet_wp, ele_wp, dm):
  """Two-panel overlay: TES (left), TauID SF (right). 3 pT bins per panel."""
  if not prof_dict['tes'] and not prof_dict['tid']:
    return None
  c = ROOT.TCanvas('c', '', 1600, 700)
  c.Divide(2, 1)
  keep1 = draw_overlay_panel(c.cd(1), prof_dict['tes'], 'TES', f'tes_{dm}',
                             dm, jet_wp, ele_wp)
  keep2 = draw_overlay_panel(c.cd(2), prof_dict['tid'], 'TauID SF',
                             f'tid_SF_{dm}_pt*', dm, jet_wp, ele_wp)
  os.makedirs(outdir, exist_ok=True)
  out = os.path.join(outdir, f'overlay_{dm}.png')
  c.SaveAs(out)
  return out


def summarize_consistency(profiles, jet_wp, ele_wp, dm, csv_writer=None):
  """Print + write pairwise σ-discrepancies between pT bins."""
  if len(profiles) < 2:
    return
  worst = 0.0
  for i in range(len(profiles)):
    for j in range(i+1, len(profiles)):
      a, b = profiles[i], profiles[j]
      sigma = discrepancy_sigma(a['bf'], (a['err_dn'], a['err_up']),
                                b['bf'], (b['err_dn'], b['err_up']))
      if sigma is None:
        continue
      worst = max(worst, sigma)
      flag = 'OK' if sigma < 2.0 else ('WARN' if sigma < 3.0 else 'BAD')
      print(f"  {dm}  pt{a['pt_idx']}({a['bf']:.4f}) "
            f"vs pt{b['pt_idx']}({b['bf']:.4f})  =  {sigma:.2f}σ  [{flag}]")
      if csv_writer is not None:
        csv_writer.writerow([jet_wp, ele_wp, dm, a['pt_idx'], b['pt_idx'],
                             f"{a['bf']:.6f}", f"{b['bf']:.6f}",
                             f"{sigma:.3f}", flag])
  return worst


# ---------- Main ----------

def main():
  ap = argparse.ArgumentParser(description=__doc__)
  ap.add_argument('--indir',  default='output_pt_less_region')
  ap.add_argument('--outdir', default='plots_pt_less_region')
  ap.add_argument('--year',   default='2025')
  ap.add_argument('--jet_wp', default=None)
  ap.add_argument('--ele_wp', default=None)
  ap.add_argument('--mode',   default='both', choices=['per-fit','overlay','both'])
  args = ap.parse_args()

  jet_glob = args.jet_wp if args.jet_wp else '*'
  ele_glob = args.ele_wp if args.ele_wp else '*'

  # Per-fit two-panel plots
  if args.mode in ('per-fit', 'both'):
    pattern = os.path.join(args.indir, f'againstjet_{jet_glob}',
                           f'againstelectron_{ele_glob}', args.year,
                           f'higgsCombine.mt_m_vis-DM*_pt*_mutau_DeepTau-{args.year}-13TeV.MultiDimFit.mH90.root')
    files = sorted(glob.glob(pattern))
    print(f">>> [per-fit] Found {len(files)} MultiDimFit files")
    ok = skip = 0
    for fp in files:
      parts = fp.split(os.sep)
      jet = next((p.replace('againstjet_','') for p in parts if p.startswith('againstjet_')), '?')
      ele = next((p.replace('againstelectron_','') for p in parts if p.startswith('againstelectron_')), '?')
      out_sub = os.path.join(args.outdir, f'againstjet_{jet}',
                             f'againstelectron_{ele}', args.year, 'nll_1d')
      res = plot_one(fp, out_sub, jet, ele)
      if res: ok += 1
      else: skip += 1
    print(f">>> [per-fit] Wrote {ok} plots, skipped {skip}")

  # Overlay TES profiles per (WP combo, DM)
  if args.mode in ('overlay', 'both'):
    jet_dirs = sorted(glob.glob(os.path.join(args.indir, f'againstjet_{jet_glob}')))
    summary_path = os.path.join(args.outdir, f'tes_consistency_{args.year}.csv')
    os.makedirs(args.outdir, exist_ok=True)
    print(f">>> [overlay] writing consistency table to {summary_path}")
    with open(summary_path, 'w', newline='') as fcsv:
      writer = csv.writer(fcsv)
      writer.writerow(['jet_wp','ele_wp','dm','pt_a','pt_b','bf_a','bf_b','sigma','flag'])
      n_plots = 0
      worst_per_combo = []
      for jd in jet_dirs:
        jet_wp = os.path.basename(jd).replace('againstjet_', '')
        ele_dirs = sorted(glob.glob(os.path.join(jd, f'againstelectron_{ele_glob}')))
        for ed in ele_dirs:
          ele_wp = os.path.basename(ed).replace('againstelectron_', '')
          print(f">>> {jet_wp} x {ele_wp}")
          for dm in ('DM0','DM1','DM10','DM11'):
            prof_dict = collect_profiles(args.indir, args.year, jet_wp, ele_wp, dm)
            if not (prof_dict['tes'] or prof_dict['tid']):
              continue
            out_sub = os.path.join(args.outdir, f'againstjet_{jet_wp}',
                                   f'againstelectron_{ele_wp}', args.year, 'nll_1d_overlay')
            if plot_overlay(prof_dict, out_sub, jet_wp, ele_wp, dm):
              n_plots += 1
            worst = summarize_consistency(prof_dict['tes'], jet_wp, ele_wp, dm, csv_writer=writer)
            if worst is not None:
              worst_per_combo.append((jet_wp, ele_wp, dm, worst))
    print(f">>> [overlay] Wrote {n_plots} overlay plots")
    # Top-line: worst offenders
    worst_per_combo.sort(key=lambda x: -x[3])
    print(">>> [overlay] Worst pT-bin discrepancies (top 10):")
    for jet_wp, ele_wp, dm, w in worst_per_combo[:10]:
      flag = 'OK' if w < 2.0 else ('WARN' if w < 3.0 else 'BAD')
      print(f"      {jet_wp:8s}  {ele_wp:8s}  {dm:5s}  {w:.2f}σ  [{flag}]")


if __name__ == '__main__':
  main()
