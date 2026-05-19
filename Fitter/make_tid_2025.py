#! /usr/bin/env python3
"""Load 2025 DM-binned tau SFs from tau_sf/ and produce a combined correctionlib JSON
with TES + TauID, both carrying the per-(DM, pT-bin) uncorrelated systematics.

Variants:
  --variant uncorr (default): per-WP JSONs named ..._DeepTau2018v2p5_<year>_VSjet*_VSele*.json
                              TES + TauID both per (DM, pT bin).
  --variant corr            : per-WP JSONs named ..._<year>_corrTES_VSjet*_VSele*.json,
                              TES is per-DM only (1 inclusive pT bin), TauID still per (DM, pT).
"""
import os, sys, json, argparse
from tau_tid_2024 import makecorr_tid, schema, JSONEncoder

HERE       = os.path.dirname(os.path.abspath(__file__))
SF_DIR     = os.path.join(HERE, 'tau_sf')
VSJET_WPS  = ['VVLoose','VLoose','Loose','Medium','Tight','VTight']
VSE_WPS    = ['VVLoose','Tight']
DMS        = [0,1,10,11]
PT_BINS    = [20.0,40.0,60.0,200.0]  # 3 pT bins for TauID (and TES uncorr)
TES_CORR_BINS = [20.0, 200.0]        # 1 inclusive bin for corrTES TES


def find_dm_node(node):
  """Recursively locate the DM-categorized node in a correctionlib JSON tree."""
  if isinstance(node, dict):
    if node.get('nodetype')=='category' and node.get('input')=='DM':
      return node
    for item in node.get('content',[]):
      r = find_dm_node(item)
      if r is not None: return r
    if 'value' in node:
      r = find_dm_node(node['value'])
      if r is not None: return r
  elif isinstance(node, list):
    for item in node:
      r = find_dm_node(item)
      if r is not None: return r
  return None


def extract_sfs(filepath):
  """Return {dm: [(nom, errup, errdown) per pT bin]} from one 2025 SF JSON."""
  with open(filepath) as f:
    data = json.load(f)
  dm_node = find_dm_node(data['corrections'][0]['data'])
  if dm_node is None:
    raise RuntimeError(f"No DM-categorized node in {filepath}")
  out = {}
  for entry in dm_node['content']:
    syst_vals = {se['key']: se['value']['content'] for se in entry['value']['content']}
    nom  = syst_vals['nom']
    up   = syst_vals['up']
    down = syst_vals['down']
    out[entry['key']] = [(nom[i], up[i]-nom[i], nom[i]-down[i]) for i in range(len(nom))]
  return out


def load_per_wp(prefix, variant='uncorr'):
  """Build dmsfs[wp_VSjet][wp_VSe][dm] = [(nom, errup, errdown) per pT bin] from
     files like {prefix}_DeepTau2018v2p5_2025[_corrTES]_VSjet{X}_VSele{Y}.json.

     Only WPs with files on disk are included — partial WP coverage works."""
  suffix = '_corrTES' if variant == 'corr' else ''
  dmsfs = {}
  for wjet in VSJET_WPS:
    for wse in VSE_WPS:
      fname = os.path.join(SF_DIR,
        f"{prefix}_DeepTau2018v2p5_2025{suffix}_VSjet{wjet}_VSele{wse}.json")
      if not os.path.exists(fname):
        print(f">>> WARNING: missing {fname}")
        continue
      dmsfs.setdefault(wjet, {})[wse] = extract_sfs(fname)
  return dmsfs


def dm_average(dmsfs):
  """Collapse dmsfs[wp][wse][dm] -> ptsfs[wp][wse] by averaging across DMs per pT bin.
     Used for the 'pt' flag (no DM split); 2025 has no separate inclusive pT measurement."""
  ptsfs = {}
  for wjet, vd in dmsfs.items():
    ptsfs[wjet] = {}
    for wse, dd in vd.items():
      n_pt = len(next(iter(dd.values())))
      avg = []
      for i in range(n_pt):
        noms = [dd[d][i][0] for d in dd]
        ups  = [dd[d][i][1] for d in dd]
        dns  = [dd[d][i][2] for d in dd]
        avg.append((sum(noms)/len(noms), sum(ups)/len(ups), sum(dns)/len(dns)))
      ptsfs[wjet][wse] = avg
  return ptsfs


def build_correction(prefix, name, tid_label, info, output_desc, variant='uncorr',
                     pt_bins=None):
  """Read per-WP files matching `prefix`, return a single schema.Correction."""
  dmsfs = load_per_wp(prefix, variant=variant)
  ptsfs = dm_average(dmsfs)
  bins  = pt_bins if pt_bins is not None else PT_BINS
  return makecorr_tid(
    ptsfs       = ptsfs,
    dmsfs       = dmsfs,
    Format      = 'Run3_May24',
    id          = tid_label,
    era         = '2025',
    wps_VSe     = VSE_WPS,
    dms         = DMS,
    bins        = bins,
    dmptbins    = bins,
    name        = name,
    info        = info,
    output_desc = output_desc,
    fname       = '',  # suppress per-correction write
    verb        = 0,
  )


def main():
  ap = argparse.ArgumentParser(description=__doc__)
  ap.add_argument('--variant', choices=['uncorr','corr'], default='uncorr',
                  help="uncorr: TES per (DM, pT); corr: TES per DM only (inclusive pT)")
  ap.add_argument('--out', default=None,
                  help="output JSON (default: data/tau/TauCorrections_2025[_corrTES]_with_uncorrelated_systs.json)")
  args = ap.parse_args()

  variant = args.variant
  # TES axis: 3 pT bins for uncorr, 1 inclusive bin for corr
  tes_pt_bins = TES_CORR_BINS if variant == 'corr' else PT_BINS

  print(f">>> Building TauID correction (variant={variant})...")
  tid_corr = build_correction(
    prefix      = 'TauID_SF_dm',
    name        = 'TauID_SF',
    tid_label   = 'DeepTau2018v2p5VSjet',
    info        = 'Tau ID SFs for DeepTau2018v2p5 in 2025',
    output_desc = 'Tau ID scale factor',
    variant     = variant,
    pt_bins     = PT_BINS,
  )
  print(f">>> Building TES correction (variant={variant})...")
  tes_corr = build_correction(
    prefix      = 'TauES_SF_dm',
    name        = 'TauEnergy_SF',
    tid_label   = 'DeepTau2018v2p5VSjet',
    info        = f"Tau Energy Scale corrections for DeepTau2018v2p5 in 2025 ({'correlated across pT per DM' if variant == 'corr' else 'per (DM, pT)'})",
    output_desc = 'Tau energy scale correction',
    variant     = variant,
    pt_bins     = tes_pt_bins,
  )

  outdir = os.path.join(HERE, 'data', 'tau')
  os.makedirs(outdir, exist_ok=True)
  if args.out is None:
    suffix = '_corrTES' if variant == 'corr' else ''
    args.out = os.path.join(outdir, f'TauCorrections_2025{suffix}_with_uncorrelated_systs.json')
  cset = schema.CorrectionSet(
    schema_version = schema.VERSION,
    description    = f"Tau ES + ID SFs for 2025 with per-(DM, pT-bin) uncorrelated systematic variations ({variant} variant)",
    corrections    = [tes_corr, tid_corr],
  )
  print(f">>> Writing {args.out}...")
  JSONEncoder.write(cset, args.out)


if __name__ == '__main__':
  main()
