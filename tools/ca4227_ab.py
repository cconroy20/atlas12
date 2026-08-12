#!/usr/bin/env python3
"""Ca I 4227 resonance-profile A/B driver.

Sets the CA4227_* developer constants in src/atlas12_modules.f90, rebuilds
SYNTHE, and runs a synthesis-only job on a frozen model atmosphere.  Used to
test the Jones et al. (2023, MNRAS 523, 1297) modified-Lorentzian solution
for the 4000-4500 A blue depression of M dwarfs.

  H(x) = exp(-pp x^2) + ps a / x^px,   x = dlam / dlam_Doppler

against (pp, ps, px) = (1, 1, 2) for the conventional far-wing form.  SYNTHE
applies MAX(Voigt, H) and truncates at |dlam| = dlmax, which the x^-px wing
requires (see the CA4227_ block in mod_parameters).

  mode  0  standard Voigt          mode 1  modified       mode -1  suppressed

Usage:
  python3 ca4227_ab.py --tag ca_base  --mode 0
  python3 ca4227_ab.py --tag ca_none  --mode -1
  python3 ca4227_ab.py --tag ca_jones --mode 1 --ps 1e-3 --px 0.5 --dlmax 500
  python3 ca4227_ab.py --reset                 # restore production mode 0

The run directory is workdir/mann/<star>_<tag>/; the .atm is copied in from
the star's baseline run so every variant shares one frozen structure.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml

SRC = os.path.join(ml.REPO, "src", "atlas12_modules.f90")
# pp=1 keeps the true Doppler core; Jones' 1e-4 is a saturated 10-30 A
# flat top that wrecks the solar line and buys nothing (see mod_parameters).
DEFAULTS = dict(mode=0, pp=1.0, ps=8.0e-6, px=0.5, dlmax=750.0)

PATTERNS = {
    "mode":  (r"(INTEGER, PARAMETER :: CA4227_MODE\s*=\s*)(-?\d+)", "{:d}"),
    "pp":    (r"(REAL\(8\), PARAMETER :: CA4227_PP\s*=\s*)(\S+)", "{:.6E}"),
    "ps":    (r"(REAL\(8\), PARAMETER :: CA4227_PS\s*=\s*)(\S+)", "{:.6E}"),
    "px":    (r"(REAL\(8\), PARAMETER :: CA4227_PX\s*=\s*)(\S+)", "{:.6E}"),
    "dlmax": (r"(REAL\(8\), PARAMETER :: CA4227_DLMAX\s*=\s*)(\S+)", "{:.6E}"),
}


def set_constants(vals):
    """Rewrite the CA4227_* PARAMETER values in mod_parameters."""
    with open(SRC) as fh:
        text = fh.read()
    for key, (pat, fmt) in PATTERNS.items():
        new = fmt.format(vals[key])
        if key != "mode":
            new = new.replace("E", "D")          # Fortran REAL(8) literal
        text, n = re.subn(pat, lambda m: m.group(1) + new, text, count=1)
        if n != 1:
            sys.exit(f"error: could not set CA4227_{key.upper()} in {SRC}")
    with open(SRC, "w") as fh:
        fh.write(text)


def build():
    subprocess.run(["make", "synthe"], cwd=os.path.join(ml.REPO, "src"),
                   check=True, stdout=subprocess.DEVNULL)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--star", default="GJ644C")
    ap.add_argument("--from-tag", default=None,
                    help="run dir supplying the frozen .atm (default: baseline)")
    ap.add_argument("--tag", help="run-dir suffix for this variant")
    ap.add_argument("--mode", type=int, default=DEFAULTS["mode"])
    ap.add_argument("--pp", type=float, default=DEFAULTS["pp"])
    ap.add_argument("--ps", type=float, default=DEFAULTS["ps"])
    ap.add_argument("--px", type=float, default=DEFAULTS["px"])
    ap.add_argument("--dlmax", type=float, default=DEFAULTS["dlmax"])
    ap.add_argument("--wlbeg", type=float, default=380.0)
    ap.add_argument("--wlend", type=float, default=520.0)
    ap.add_argument("--resolu", type=float, default=300000.0)
    ap.add_argument("--atm", default=None,
                    help="explicit .atm path, bypassing star resolution")
    ap.add_argument("--outdir", default=None, help="run dir when --atm is used")
    ap.add_argument("--reset", action="store_true",
                    help="restore production constants, rebuild, and exit")
    a = ap.parse_args()

    if a.reset:
        set_constants(DEFAULTS)
        build()
        print("CA4227 constants reset to production (mode 0); synthe rebuilt")
        return
    if not a.tag:
        sys.exit("error: --tag is required")

    vals = dict(mode=a.mode, pp=a.pp, ps=a.ps, px=a.px, dlmax=a.dlmax)
    set_constants(vals)
    build()

    if a.atm:                       # explicit atmosphere (e.g. the Sun)
        # Absolute: synthe.exe is launched with cwd = rundir, so a relative
        # model path would be resolved against the wrong directory.
        base = os.path.splitext(os.path.basename(a.atm))[0]
        src_atm = os.path.abspath(a.atm)
        rundir = os.path.abspath(a.outdir or (a.atm + "_" + a.tag))
    else:
        s = ml.resolve(a.star)
        base = s.name
        src_atm = os.path.join(ml.rundir_for(s, a.from_tag), base + ".atm")
        rundir = ml.rundir_for(s, a.tag)
    os.makedirs(rundir, exist_ok=True)
    atm = os.path.join(rundir, base + ".atm")
    if not os.path.isfile(atm):
        shutil.copy(src_atm, atm)

    print(f"{base}  tag={a.tag}  mode={a.mode}"
          + (f"  pp={a.pp:g} ps={a.ps:g} px={a.px:g} dlmax={a.dlmax:g}"
             if a.mode > 0 else ""))
    print(f"synthe {a.wlbeg:.0f}-{a.wlend:.0f} nm at R={a.resolu:.0f} -> {rundir}")
    t0 = time.time()
    ml.run([os.path.join(ml.REPO, "bin", "synthe.exe"), atm,
            f"wlbeg={a.wlbeg}", f"wlend={a.wlend}", f"resolu={a.resolu:.0f}"],
           rundir, os.path.join(rundir, base + ".synthe.stdout"))
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
