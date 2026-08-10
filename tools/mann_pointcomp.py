#!/usr/bin/env python3
"""Batch driver for the standing point-comparison ladder: the curated
~10-star Teff sequence (mann_lib.STAR_LADDER, 2700-4800 K, GJ names,
Mann+15 M dwarfs + Mann+13 K-dwarf calibrators).

For each star: atlas12 + full-range R=300k synthe (via validate_mann.py,
skip-if-converged), the 3-panel comparison figure vs Mann data + PHOENIX
NewEra (mann_compare_plot.plot_star), and a metrics row.  Ends with a
summary table written to workdir/mann/summary.md.

Usage:
  python3 mann_pointcomp.py                 # whole ladder, 1 job
  python3 mann_pointcomp.py --jobs 4        # 4-way parallel
  python3 mann_pointcomp.py --star GJ887    # subset (comma-separated ok)
  python3 mann_pointcomp.py --force         # rerun atlas12+synthe
"""
import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml

TOOLS = os.path.dirname(os.path.abspath(__file__))

SUMMARY_COLS = [
    ("teff", "Teff", "{:.0f}"), ("logg", "logg", "{:.2f}"),
    ("feh", "[Fe/H]", "{:+.2f}"), ("rms", "RMS%", "{:.2f}"),
    ("integral", "int", "{:.3f}"), ("opt_med", "opt", "{:.3f}"),
    ("nir_med", "NIR", "{:.3f}"), ("tio5", "TiO5", "{:.3f}"),
    ("tio5_data", "TiO5d", "{:.3f}"), ("cah2", "CaH2", "{:.3f}"),
    ("cah2_data", "CaH2d", "{:.3f}"), ("cah3", "CaH3", "{:.3f}"),
    ("cah3_data", "CaH3d", "{:.3f}"),
]


def run_one(name, args):
    """validate_mann.py for one star (subprocess, own log); returns rc."""
    log = os.path.join(ml.RUN_ROOT, f"{name}.pointcomp.log")
    cmd = [sys.executable, os.path.join(TOOLS, "validate_mann.py"),
           "--star", name, "--numit", str(args.numit),
           "--resolu", str(args.resolu)]
    if args.force:
        cmd.append("--force")
    with open(log, "w") as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    status = "ok" if rc == 0 else f"FAILED rc={rc} (see {log})"
    print(f"[{name}] validate: {status}")
    return rc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--star", default=None,
                    help="comma-separated subset (default: whole ladder)")
    ap.add_argument("--jobs", type=int, default=1, help="parallel stars")
    ap.add_argument("--numit", type=int, default=30)
    ap.add_argument("--resolu", type=float, default=300000.0)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-plots", action="store_true",
                    help="skip figures (metrics/summary only)")
    args = ap.parse_args()

    names = ([n.strip() for n in args.star.split(",")] if args.star
             else [n for n, _ in ml.STAR_LADDER])
    os.makedirs(ml.RUN_ROOT, exist_ok=True)

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        rcs = dict(zip(names, ex.map(lambda n: run_one(n, args), names)))

    rows = []
    for name in names:
        s = ml.resolve(name)
        met_path = os.path.join(ml.rundir_for(s), s.name + "_metrics.json")
        met = None
        if rcs[name] == 0 and os.path.isfile(met_path):
            with open(met_path) as fh:
                met = json.load(fh)
            if not args.no_plots:
                try:
                    from mann_compare_plot import plot_star
                    plot_star(s)
                except Exception as e:
                    print(f"[{name}] plot FAILED: {e}")
        rows.append((s, met))

    # ---------------- summary table ----------------
    hdr = ["star"] + [h for _, h, _ in SUMMARY_COLS]
    lines = ["| " + " | ".join(hdr) + " |",
             "|" + "|".join("---" for _ in hdr) + "|"]
    for s, met in rows:
        if met is None:
            lines.append(f"| {s.name} | " +
                         " | ".join("FAIL" if k == "teff" else ""
                                    for k, _, _ in SUMMARY_COLS) + " |")
            continue
        cells = [fmt.format(met[k]) if k in met else "-"
                 for k, _, fmt in SUMMARY_COLS]
        lines.append("| " + " | ".join([s.name] + cells) + " |")
    table = "\n".join(lines)
    print("\n" + table)
    out = os.path.join(ml.RUN_ROOT, "summary.md")
    with open(out, "w") as fh:
        fh.write("# Mann point-comparison ladder\n\n"
                 "model/data metrics; *_d columns are the data-side index "
                 "values.  Full-range R=300k syntheses, measured Mann LSF, "
                 "absolute fluxes (no renormalization).\n\n" + table + "\n")
    print(f"\nsummary: {out}")


if __name__ == "__main__":
    main()
