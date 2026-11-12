#!/usr/bin/env python3
"""Run one ATLAS12 model atmosphere from its stellar parameters.

    run_atlas12_model.py --teff 5750 --logg 4.5 --feh 0.0 --afe 0.0 \\
        --vmic 1.0 --y 0.27 --mlt 2.0 [--start DECK_OR_DIR] [--outdir DIR]

It writes <name>.abund, <name>.atm, <name>.flux, <name>.iter and <name>.log in
--outdir, then prints how well the final model conserves flux.

Composition: Bergemann, Lodders & Palme (2025) solar abundances, scaled by
[Fe/H], with [alpha/Fe] added to O, Ne, Mg, Si, S, Ar, Ca and Ti.  Y is the
helium mass fraction.  Metals are anchored to hydrogen, so [X/H] stays exact as
Y changes (X and Z both fall as Y rises, Z/X fixed).  This is the composition
convention of the CWC v3 atmosphere grid.

Solver settings are those of the CWC v3 grid: the transactional deep-CZ
polish, early stopping (numit is a ceiling), the convection gap gate
(cnv_gap_del=-0.03) and, for log g > 2 only, convection-zone damping
(cz_damp_del=0.1; it degrades low-gravity models).  --atlas-arg passes any
other key=value straight to atlas12c.exe and wins over these.

Start deck: ATLAS rescales the start model to the requested Teff and log g, but
convergence is faster and more reliable from a nearby structure.  --start takes
one deck, or a directory of decks named *_t<Teff>g<logg>.atm, from which the
nearest in (log Teff, log g) is used.  The default is the solar deck in
workdir/.

Needs ATLAS12 built (bin/atlas12c.exe) and its data/ installed, including
data/mol/h2opokazatel.bin for Teff <= 8000 K: without that file ATLAS12
silently runs with no H2O opacity, so this script refuses to start, and checks
the log afterwards.  Python 3.8+, standard library only.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_START = REPO / "workdir" / "ap00t5777g4.44at12.dat"

# Bergemann, Lodders & Palme (2025), A(X) = 12 + log10(N_X/N_H), Z = 1..92;
# LOGEPS_BERG25 in src/atlas12_modules.f90.  -20 marks no solar abundance.
BERG25_AX = [
    12.00, 10.922, 1.04, 1.21, 2.70, 8.51, 7.94, 8.76, 4.40, 8.15,
    6.29, 7.58, 6.43, 7.56, 5.44, 7.16, 5.43, 6.50, 5.09, 6.35,
    3.13, 4.97, 3.89, 5.74, 5.52, 7.51, 4.95, 6.24, 4.24, 4.55,
    3.02, 3.62, 2.34, 3.41, 2.65, 3.31, 2.35, 2.93, 2.30, 2.68,
    1.47, 1.88, -20.0, 1.75, 0.78, 1.57, 0.96, 1.77, 0.80, 2.02,
    1.08, 2.23, 1.76, 2.30, 1.12, 2.27, 1.10, 1.58, 0.75, 1.42,
    -20.0, 0.95, 0.57, 1.08, 0.31, 1.10, 0.48, 0.93, 0.11, 0.85,
    0.10, 0.86, -0.11, 0.79, 0.30, 1.36, 1.42, 1.64, 0.91, 1.14,
    0.95, 1.95, 0.70, -20.0, -20.0, -20.0, -20.0, -20.0, -20.0, 0.09,
    -20.0, -0.50,
]
ALPHA_Z = {8, 10, 12, 14, 16, 18, 20, 22}  # O Ne Mg Si S Ar Ca Ti
# Standard atomic weights for Z = 1..30; beyond Zn ~2.45 Z (< 1e-4 of Z).
ATOMIC_MASS = [
    1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078,
    44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38,
] + [2.45 * z for z in range(31, 93)]

H2O_TEFF_MAX = 8000          # SELECTLINES reads the H2O list at Teff <= this
CZ_DAMP_LOGG_MIN = 2.0       # cz_damp_del=0.1 only above this log g
LINEFORMING_LOGTAU = (-5.0, 1.0)
# Early-stop thresholds: max/p95 |flux error| in percent.
FLUX_MAX, FLUX_P95 = 1.0, 0.5
DECK_RE = re.compile(r"_t(\d+(?:\.\d+)?)g([+-]?\d+(?:\.\d+)?)\.atm$")


def abundances(feh: float, afe: float, y: float) -> tuple[float, dict[int, float], tuple]:
    """N_He/N_total, {Z: log10(N_Z/N_total)} for Z >= 3, and (X, Y, Z)."""
    if not 0.0 < y < 1.0:
        raise SystemExit(f"ERROR: --y {y} must lie in (0, 1)")
    over_h = {}
    for z in range(3, 93):
        ax = BERG25_AX[z - 1]
        if ax > -19.0:
            over_h[z] = ax - 12.0 + feh + (afe if z in ALPHA_Z else 0.0)
    ratio = {z: 10.0 ** v for z, v in over_h.items()}
    metal_mass = sum(r * ATOMIC_MASS[z - 1] for z, r in ratio.items())
    he_over_h = y * (ATOMIC_MASS[0] + metal_mass) / (ATOMIC_MASS[1] * (1.0 - y))
    ntot_over_h = 1.0 + he_over_h + sum(ratio.values())
    total_mass = ATOMIC_MASS[0] + he_over_h * ATOMIC_MASS[1] + metal_mass
    xyz = (ATOMIC_MASS[0] / total_mass, he_over_h * ATOMIC_MASS[1] / total_mass,
           metal_mass / total_mass)
    log_ntot = math.log10(ntot_over_h)
    return he_over_h / ntot_over_h, {z: v - log_ntot for z, v in over_h.items()}, xyz


def write_abund(path: Path, feh: float, afe: float, y: float) -> tuple:
    he, metals, xyz = abundances(feh, afe, y)
    with open(path, "w") as f:
        f.write(f"# ATLAS-12 abund= overrides  (berg25 base, [Fe/H]={feh:+.2f} "
                f"[alpha/Fe]={afe:+.2f})\n")
        f.write(f"# helium: Y = {y:.5f} (mass); metals anchored to H\n")
        f.write(f"# alpha elements (Z, +[alpha/Fe]): {sorted(ALPHA_Z)}\n")
        f.write("# Z   log10(N_Z/N_total)  (Z = 2: N_He/N_total, linear)\n")
        f.write(f"{2:3d}  {he:9.6f}\n")
        for z, v in metals.items():
            f.write(f"{z:3d}  {v:9.4f}\n")
    return xyz


def choose_start(start: Path, teff: float, logg: float) -> Path:
    if start.is_file():
        return start
    if not start.is_dir():
        raise SystemExit(f"ERROR: --start {start} is neither a file nor a directory")
    decks = []
    for path in start.glob("*.atm"):
        m = DECK_RE.search(path.name)
        if m and path.stat().st_size > 0:
            t, g = float(m.group(1)), float(m.group(2))
            d2 = (math.log10(t / teff) / 0.04) ** 2 + ((g - logg) / 0.5) ** 2
            decks.append((d2, path.name, path))
    if not decks:
        raise SystemExit(f"ERROR: no *_t<Teff>g<logg>.atm decks in {start}")
    return min(decks)[2]


def read_blocks(path: Path) -> list[list[list[float]]]:
    """Per-layer rows of every .iter table (J, log10TAU, T, ..., ERROR in col 7)."""
    blocks, rows, in_table = [], [], False
    number = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d{2})?")
    for line in path.read_text(errors="replace").splitlines():
        s = line.lstrip()
        if s.startswith("J") and "log10TAU" in s and "ERROR" in s:
            if rows:
                blocks.append(rows)
            rows, in_table = [], True
            continue
        if not in_table:
            continue
        fields = line.split()
        if len(fields) < 17:
            # The fixed-width writer can run two columns together.
            fields = number.findall(line)
        try:
            values = [float(v) for v in fields[:17]]
        except ValueError:
            values = []
        if len(values) == 17 and float(values[0]).is_integer():
            if rows and values[0] <= rows[-1][0]:
                blocks.append(rows)
                rows = []
            rows.append(values)
        elif rows and not line.strip().startswith(("K", "#")) and line.strip():
            in_table = False
    if rows:
        blocks.append(rows)
    return blocks


def flux_stats(rows: list[list[float]]) -> tuple[float, float]:
    err = sorted(abs(r[7]) for r in rows)
    if not err:
        return math.nan, math.nan
    # p95 by linear interpolation, as numpy.percentile does.
    k = 0.95 * (len(err) - 1)
    lo = math.floor(k)
    p95 = err[lo] + (err[min(lo + 1, len(err) - 1)] - err[lo]) * (k - lo)
    return err[-1], p95


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--teff", type=float, required=True, help="effective temperature, K")
    ap.add_argument("--logg", type=float, required=True, help="log10 surface gravity, cgs")
    ap.add_argument("--feh", type=float, required=True, help="[Fe/H]")
    ap.add_argument("--afe", type=float, required=True, help="[alpha/Fe]")
    ap.add_argument("--vmic", type=float, required=True, help="microturbulence, km/s")
    ap.add_argument("--y", type=float, required=True, help="helium mass fraction Y")
    ap.add_argument("--mlt", type=float, required=True, help="mixing length alpha_MLT")
    ap.add_argument("--start", type=Path, default=DEFAULT_START,
                    help="start deck, or a directory of *_t<Teff>g<logg>.atm decks "
                         "(default: the solar deck in workdir/)")
    ap.add_argument("--outdir", type=Path, default=Path("."), help="output directory")
    ap.add_argument("--name", help="output basename (default: built from the parameters)")
    ap.add_argument("--numit", type=int, default=60, help="iteration ceiling (default 60)")
    ap.add_argument("--atlas12", type=Path, default=Path(os.environ.get("ATLAS12", REPO)),
                    help="ATLAS12 installation with data/ (default: $ATLAS12, else this repo)")
    ap.add_argument("--exe", type=Path, help="atlas12c.exe (default: <atlas12>/bin/atlas12c.exe)")
    ap.add_argument("--atlas-arg", action="append", default=[], metavar="KEY=VALUE",
                    help="extra atlas12c.exe option; repeatable; overrides the defaults")
    args = ap.parse_args()

    if args.mlt <= 0 or args.vmic < 0:
        raise SystemExit("ERROR: --mlt must be > 0 and --vmic >= 0")
    atlas12 = args.atlas12.resolve()
    exe = (args.exe or atlas12 / "bin" / "atlas12c.exe").resolve()
    if not os.access(exe, os.X_OK):
        raise SystemExit(f"ERROR: not an executable: {exe}")
    need_h2o = args.teff <= H2O_TEFF_MAX
    h2o = atlas12 / "data" / "mol" / "h2opokazatel.bin"
    if need_h2o and not (h2o.is_file() and h2o.stat().st_size > 0):
        raise SystemExit(f"ERROR: Teff <= {H2O_TEFF_MAX} K needs the H2O line list, but "
                         f"{h2o} is missing or empty (ATLAS12 would silently omit H2O)")
    start = choose_start(args.start.resolve(), args.teff, args.logg)

    name = args.name or (f"at12_t{args.teff:05.0f}g{args.logg:+.2f}_feh{args.feh:+.2f}"
                         f"_afe{args.afe:+.2f}_v{args.vmic:.2f}_y{args.y:.3f}_mlt{args.mlt:.2f}")
    if "/" in name:
        raise SystemExit("ERROR: --name may not contain '/'")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    abund = outdir / f"{name}.abund"
    x, y, z = write_abund(abund, args.feh, args.afe, args.y)

    options = {
        "teff": f"{args.teff:g}", "logg": f"{args.logg:g}", "vturb": f"{args.vmic:g}",
        "mlt": f"{args.mlt:g}", "solar": "berg25", "abund": abund.name,
        "numit": str(args.numit), "czc_polish": "transactional", "early_stop": "on",
        "minit": "10", "conv_streak": "3", "cnv_gap_del": "-0.03",
        "cz_damp_del": "0.1" if args.logg > CZ_DAMP_LOGG_MIN else "-1",
    }
    for item in args.atlas_arg:
        key, sep, value = item.partition("=")
        if not sep or not key:
            raise SystemExit(f"ERROR: --atlas-arg {item!r} is not KEY=VALUE")
        options[key] = value
    command = [str(exe), str(start), name] + [f"{k}={v}" for k, v in options.items()]

    print(f"{name}\n  Teff={args.teff:g} logg={args.logg:g} [Fe/H]={args.feh:+g} "
          f"[a/Fe]={args.afe:+g} Vmic={args.vmic:g} Y={args.y:g} MLT={args.mlt:g}")
    print(f"  X={x:.5f} Y={y:.5f} Z={z:.6f}")
    print(f"  start deck: {start}")
    print(f"  running: {' '.join(command)}", flush=True)
    log = outdir / f"{name}.log"
    with open(log, "w") as handle:
        status = subprocess.run(command, cwd=outdir, stdout=handle, stderr=subprocess.STDOUT,
                                env={**os.environ, "ATLAS12": str(atlas12) + "/"}).returncode
    atm, iters = outdir / f"{name}.atm", outdir / f"{name}.iter"
    if status != 0 or not atm.is_file() or atm.stat().st_size == 0:
        raise SystemExit(f"ERROR: atlas12c.exe exited {status} without a model; see {log}")
    if need_h2o and "lines from h2opokazatel" not in log.read_text(errors="replace"):
        raise SystemExit(f"ERROR: ATLAS12 did not read the H2O list; {atm.name} has no H2O "
                         f"opacity (see {log})")

    blocks = read_blocks(iters) if iters.is_file() else []
    if not blocks:
        print(f"  WARNING: no iteration tables in {iters}; convergence unknown")
        return
    final = blocks[-1]
    window = [r for r in final if LINEFORMING_LOGTAU[0] <= r[1] <= LINEFORMING_LOGTAU[1]]
    fmax, fp95 = flux_stats(final)
    wmax, wp95 = flux_stats(window)
    print(f"  final flux error: max {fmax:.3f}%  p95 {fp95:.3f}% (all depths); "
          f"max {wmax:.3f}%  p95 {wp95:.3f}% (-5 <= log tau <= 1)")
    if fmax <= FLUX_MAX and fp95 <= FLUX_P95:
        verdict = "converged at all depths"
    elif wmax <= FLUX_MAX and wp95 <= FLUX_P95:
        verdict = "converged in the line-forming layers only; deeper or outer layers are not"
    else:
        verdict = "NOT converged"
    print(f"CONVERGENCE: {verdict} ({len(blocks)} iteration tables; outputs in {outdir})")


if __name__ == "__main__":
    main()
