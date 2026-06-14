#!/usr/bin/env python3
"""Nuisance pull plotter (matplotlib) for the TES/TauID fit.

Reads the text table produced by combine's diffNuisances.py (after the
`sed 's/[!,]/ /g' | tail -n +4` cleanup used in the workflow), i.e. rows of:

    <name>  <b_val> <b_err>  <sb_val> <sb_err>  <rho>  <impact>

and draws a grouped, readable pull plot:
  * physics nuisances grouped by category with human-readable labels
  * ±1σ (green) / ±2σ (yellow) bands, B-only vs S+B markers
  * labels coloured when |pull|>1 (orange) or |pull|>2 (red)
  * MC-stat (prop_bin*) hidden by default (summarised); --all shows them

Usage:
  plot_pulls.py -f pulls.txt -o out_basename -t "Medium/Tight DM10"
  plot_pulls.py -f pulls.txt -o out_basename --all       # include MC-stat
"""
import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ---- nuisance -> (category, pretty label) -------------------------------
def classify(name):
    """Return (category, pretty_label) for a nuisance branch name."""
    m = re.match(r"prop_bin(.+?)_bin(\d+)$", name)
    if m:
        return "MC stat", f"MC stat {m.group(1)} bin{m.group(2)}"
    table = {
        "lumi":        ("Normalization", "luminosity"),
        "eff_m":       ("Normalization", "muon efficiency"),
        "muonFakerate":("Normalization", "muon fake rate"),
        "norm_qcd":    ("Normalization", "QCD norm"),
        "norm_wj":     ("Normalization", "W norm"),
        "xsec_dy":     ("Cross sections", r"$\sigma$(DY)"),
        "xsec_tt":     ("Cross sections", r"$\sigma$(t$\bar{t}$)"),
        "xsec_st":     ("Cross sections", r"$\sigma$(single t)"),
        "xsec_vv":     ("Cross sections", r"$\sigma$(VV)"),
    }
    if name in table:
        return table[name]
    m = re.match(r"shape_(jTauFake|mTauFake|ttbar)_(DM\d+_pt\d+)$", name)
    if m:
        kind = {"jTauFake": r"j$\to\tau$ shape", "mTauFake": r"$\mu\to\tau$ shape",
                "ttbar": r"t$\bar{t}$ shape"}[m.group(1)]
        return "Shape", f"{kind} {m.group(2)}"
    m = re.match(r"rate_(jTauFake)_(DM\d+_pt\d+)$", name)
    if m:
        return "Rate", rf"j$\to\tau$ rate {m.group(2)}"
    return "Other", name


CATEGORY_ORDER = ["Normalization", "Cross sections", "Shape", "Rate", "Other", "MC stat"]
CATEGORY_COLOR = {
    "Normalization": "#1f77b4", "Cross sections": "#2ca02c", "Shape": "#9467bd",
    "Rate": "#ff7f0e", "Other": "#7f7f7f", "MC stat": "#8c8c8c",
}


def parse(path):
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) < 5:
                continue
            try:
                pb, eb, ps, es = float(t[1]), float(t[2]), float(t[3]), float(t[4])
            except ValueError:
                continue
            cat, label = classify(t[0])
            rows.append(dict(name=t[0], cat=cat, label=label,
                             pb=pb, eb=eb, ps=ps, es=es))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--filename", required=True, help="diffNuisances text table")
    ap.add_argument("-o", "--outname", default="pulls", help="output basename (.png/.pdf)")
    ap.add_argument("-t", "--text", default="", help="title text (e.g. 'Medium/Tight DM10')")
    ap.add_argument("--all", action="store_true", help="also show MC-stat (prop_bin*) nuisances")
    ap.add_argument("--no-bonly", action="store_true", help="show only S+B fit")
    args = ap.parse_args()

    rows = parse(args.filename)
    if not rows:
        print(f">>> plot_pulls: no pulls parsed from {args.filename}")
        return

    n_mcstat = sum(r["cat"] == "MC stat" for r in rows)
    n_mc_1s = sum(r["cat"] == "MC stat" and abs(r["ps"]) > 1 for r in rows)
    n_mc_2s = sum(r["cat"] == "MC stat" and abs(r["ps"]) > 2 for r in rows)
    if not args.all:
        rows = [r for r in rows if r["cat"] != "MC stat"]
    # order by category then by |S+B pull| descending within category
    rows.sort(key=lambda r: (CATEGORY_ORDER.index(r["cat"]), -abs(r["ps"])))

    n = len(rows)
    fig_h = max(2.6, 0.34 * n + 1.5)
    fig, ax = plt.subplots(figsize=(8.4, fig_h))

    # bands
    ax.axvspan(-2, 2, color="#fff3c4", zorder=0)   # ±2σ yellow
    ax.axvspan(-1, 1, color="#c8e6c9", zorder=0)   # ±1σ green
    for x in (-2, -1, 1, 2):
        ax.axvline(x, color="0.7", lw=0.8, ls=":", zorder=1)
    ax.axvline(0, color="0.4", lw=1.0, zorder=1)

    ypos = list(range(n - 1, -1, -1))  # top row first
    off = 0.18
    for y, r in zip(ypos, rows):
        if not args.no_bonly:
            ax.errorbar(r["pb"], y + off, xerr=r["eb"], fmt="s", ms=4,
                        color="0.25", elinewidth=1.3, capsize=2, zorder=4)
        col = "#d62728" if abs(r["ps"]) > 2 else (
            "#ff7f0e" if abs(r["ps"]) > 1 else "#000000")
        ax.errorbar(r["ps"], y - off, xerr=r["es"], fmt="o", ms=5,
                    color=col, elinewidth=1.6, capsize=2, zorder=5)

    # y labels coloured by pull size; category color as a left dash
    ax.set_yticks(ypos)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8)
    for tick, r in zip(ax.get_yticklabels(), rows):
        if abs(r["ps"]) > 2:
            tick.set_color("#d62728"); tick.set_fontweight("bold")
        elif abs(r["ps"]) > 1:
            tick.set_color("#ff7f0e")

    # category separators + labels
    seen = {}
    for y, r in zip(ypos, rows):
        seen.setdefault(r["cat"], []).append(y)
    for cat, ys in seen.items():
        ax.axhline(max(ys) + 0.5, color=CATEGORY_COLOR[cat], lw=2.5, alpha=0.5, zorder=2)
        ax.text(3.02, sum(ys) / len(ys), cat, rotation=90, va="center", ha="left",
                fontsize=7.5, color=CATEGORY_COLOR[cat], fontweight="bold")

    ax.set_xlim(-3, 3)
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_xlabel(r"$(\hat{\theta} - \theta_0)\,/\,\Delta\theta$", fontsize=11)
    title = args.text
    if not args.all and n_mcstat:
        title += f"   [+{n_mcstat} MC-stat hidden: {n_mc_1s}>1σ, {n_mc_2s}>2σ]"
    ax.set_title(title, fontsize=10, loc="left")

    handles = [Line2D([], [], marker="o", color="#000", ls="", label="S+B fit"),
               Line2D([], [], marker="s", color="0.25", ls="", label="B-only fit")]
    ax.legend(handles=handles[:1] if args.no_bonly else handles,
              loc="lower right", fontsize=8, framealpha=0.9)
    ax.grid(axis="x", color="0.9", lw=0.5, zorder=0)
    fig.tight_layout()

    for ext in ("png", "pdf"):
        out = f"{args.outname}.{ext}"
        fig.savefig(out, dpi=140)
    print(f">>> plot_pulls: wrote {args.outname}.png/.pdf  ({n} nuisances shown)")


if __name__ == "__main__":
    main()
