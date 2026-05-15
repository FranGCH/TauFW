#!/usr/bin/env python3
"""Combine per-WP TauID and TES measurement PNGs into 6x2 grids."""
import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

JET_WPS = ["VVLoose", "VLoose", "Loose", "Medium", "Tight", "VTight"]
ELE_WPS = ["VVLoose", "Tight"]
QUANTITIES = ["tauID_measurements", "tes_measurements"]


def build_grid(measurements_dir, quantity, out_path):
    nrows, ncols = len(JET_WPS), len(ELE_WPS)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows))
    if nrows == 1: axes = [axes]
    if ncols == 1: axes = [[ax] for ax in axes]

    missing = []
    for i, jet in enumerate(JET_WPS):
        for j, ele in enumerate(ELE_WPS):
            ax = axes[i][j]
            png = os.path.join(measurements_dir, f"VSjet{jet}_VSele{ele}", f"{quantity}.png")
            if os.path.exists(png):
                ax.imshow(mpimg.imread(png))
            else:
                missing.append(png)
                ax.text(0.5, 0.5, "missing", ha="center", va="center", transform=ax.transAxes, color="red")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values(): spine.set_visible(False)
            if i == 0:
                ax.set_title(f"VSele {ele}", fontsize=20, pad=10)
            if j == 0:
                ax.set_ylabel(f"VSjet {jet}", fontsize=20, labelpad=15, rotation=90)

    title = "TauID SF" if quantity.startswith("tauID") else "TES"
    fig.suptitle(f"{title} measurements per (VSjet, VSele) WP", fontsize=24, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"[grid] wrote {out_path}")
    if missing:
        print(f"[grid] WARNING: {len(missing)} missing PNG(s):")
        for m in missing: print(f"   - {m}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "Measurements"),
                   help="Path to Fitter/Measurements")
    p.add_argument("--out", default=None, help="Output directory (default: same as --dir)")
    args = p.parse_args()

    if not os.path.isdir(args.dir):
        print(f"ERROR: {args.dir} not found"); sys.exit(1)
    out_dir = args.out or args.dir
    os.makedirs(out_dir, exist_ok=True)

    for q in QUANTITIES:
        build_grid(args.dir, q, os.path.join(out_dir, f"grid_{q}.png"))


if __name__ == "__main__":
    main()
