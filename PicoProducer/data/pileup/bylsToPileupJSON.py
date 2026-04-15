#!/usr/bin/env python3
# Convert brilcalc --byls CSV into pileupCalc.py input JSON.
#
# Input CSV columns (from: brilcalc lumi --byls --minBiasXsec 69200 ... --output-style csv):
#   run:fill, ls:ls, time, beamstatus, E(GeV), delivered(/ub), recorded(/ub), avgpu, source
#
# Output JSON layout expected by RecoLuminosity.LumiDB.pileupParser:
#   { "<run>": [ [ls, intLumi, RMS/xsec, mean/xsec], ... ], ... }
# where:
#   - intLumi   : recorded lumi per LS in /ub (pileupCalc.py uses it as LSintLumi weight)
#   - RMS/xsec  : 0 (byls does not expose an in-LS RMS)
#   - mean/xsec : avgpu / MINBIAS_XSEC_BRIL, so that pileupCalc.py's
#                 AveNumInt = lumiInfo[2] * user_xsec = avgpu * user_xsec / MINBIAS_XSEC_BRIL
#                 recovers the correct Poisson mean for any user_xsec.

import csv, json, argparse, sys

def main():
  ap = argparse.ArgumentParser()
  ap.add_argument('-i','--input',  required=True, help='byls CSV from brilcalc')
  ap.add_argument('-o','--output', required=True, help='pileup_JSON.txt (pileupCalc.py input)')
  ap.add_argument('-x','--minBiasXsec', type=float, default=69200.0,
                  help='minBiasXsec passed to brilcalc --minBiasXsec (ub, default 69200)')
  args = ap.parse_args()

  byrun = {}  # run -> list of [ls, intLumi, 0, mean/xsec]
  nrows = 0
  with open(args.input, 'r') as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue
      parts = line.split(',')
      if len(parts) < 9:
        continue
      run = int(parts[0].split(':')[0])
      ls  = int(parts[1].split(':')[0])
      try:
        recorded = float(parts[6])
        avgpu    = float(parts[7])
      except ValueError:
        continue
      mean_over_xsec = avgpu / args.minBiasXsec
      byrun.setdefault(run, []).append([ls, recorded, 0.0, mean_over_xsec])
      nrows += 1

  for run in byrun:
    byrun[run].sort(key=lambda r: r[0])

  out = {str(run): rows for run, rows in sorted(byrun.items())}
  with open(args.output, 'w') as f:
    json.dump(out, f, separators=(',', ':'))
    f.write('\n')

  nruns = len(out)
  nls   = sum(len(v) for v in out.values())
  print(f"Wrote {args.output}: {nruns} runs, {nls} LS rows (read {nrows})")

if __name__ == '__main__':
  main()
