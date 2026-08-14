#!/usr/bin/env python3
"""Build a <model>.dep sidecar for SYNTHE's NLTE_MODE = 2 from a published
departure-coefficient grid.

WHAT THIS DOES, AND WHY IT IS NOT IN SYNTHE
-------------------------------------------
The published grids (Amarsi et al. 2020, as repackaged by Gerber et al. 2023
for Turbospectrum) are multi-gigabyte binaries holding departure coefficients
for thousands of MARCS models, addressed through a separate text index and
carrying their own model-atom level ordering.  A reader for that inside SYNTHE
would bury two lines of physics under several hundred lines of file format, and
would need rewriting every time the grid is reissued.  So the grid I/O, the
stellar-parameter selection, the depth-scale conversion and the level mapping
happen here, and SYNTHE reads one small self-describing text file per model.

THE DEPTH-SCALE PROBLEM, WHICH IS THE REAL WORK
-----------------------------------------------
The grids are tabulated against tau_5000.  SYNTHE has no tau_5000: an ATLAS12
.atm carries RHOX (column mass) and ABROSS, so column mass is exact and
tau_Rosseland is derivable, but the 5000 A optical depth is not stored.

Converting via *our* tau would be wrong anyway.  Optical depth depends on the
continuum opacity, so equating tau between two codes whose opacities differ
silently equates layers at different physical conditions -- the failure that
once faked a 440 K temperature difference between this code and PHOENIX at
2700 K.  Column mass has no such dependence.

So the conversion is done through the atmosphere the grid was ACTUALLY computed
on: the index file names a MARCS model, that model tabulates lgTau5 and RHOX on
the same layers, and those two columns give tau_5000 -> column mass exactly, for
the structure the departure coefficients belong to.  SYNTHE then interpolates
in column mass onto its own layers.  This is why --marcs is required and not a
convenience.

LEVEL SELECTION
---------------
b is per model-atom level, and the atom's level ordering is a property of the
grid, not of us.  --levels states the mapping explicitly and is the reliable
path.  --atom will try to find the levels by matching term energies against the
transition table below, but it prints what it chose and refuses to guess when
the match is ambiguous: silently picking the wrong level would produce a
plausible, entirely wrong spectrum.

Na I D shares a lower level (3s ground) and, in atoms that resolve the 3p fine
structure, has two distinct upper levels.  Atoms that do NOT resolve it give
both D lines the same upper level; that is legitimate and the tool says so.

USAGE
-----
  python3 tools/nlte_make_dep.py OUT.dep \\
      --bin NLTEgrid_Na.bin --aux auxData_Na.txt \\
      --marcs-dir /path/to/marcs/ \\
      --teff 5777 --logg 4.44 --feh 0.0 \\
      --levels 1,3,1,2

  python3 tools/nlte_make_dep.py --selftest

The transition order (--levels, and the sidecar columns) must match
NLTE_TR_* in src/mod_mklinelist.f90: 1 = Na I D2, 2 = Na I D1.

NOTE ON VALIDATION: the binary and index layouts here follow the Turbospectrum
NLTE reader (utilities/read_nlte_grid.py in TSFitPy, and
interpolator/interpol_modeles_nlte_gfort.f in Turbospectrum_NLTE).  --selftest
exercises every code path against fabricated files in that layout, but the tool
has NOT yet been run against a real downloaded grid.  Check the printed level
assignment and the tau/column-mass span the first time you do.
"""
import argparse
import os
import struct
import sys

import numpy as np

# Term energies [cm^-1] of the transitions, in the order SYNTHE expects.
# Must match NLTE_TR_ELO / NLTE_TR_EUP in src/mod_mklinelist.f90.
TRANSITIONS = [
    ("Na I D2  3s 2S1/2 - 3p 2P3/2", 0.000, 16973.366),
    ("Na I D1  3s 2S1/2 - 3p 2P1/2", 0.000, 16956.170),
]
ETOL = 5.0          # cm^-1, matching NLTE_TR_ETOL


# ---------------------------------------------------------------------------
#  Grid index ("auxiliary") file
# ---------------------------------------------------------------------------
def read_aux(path):
    """Index of the binary grid.

    One record per model: id, Teff, log g, [Fe/H], [alpha/Fe], mass, vturb,
    the element abundance the grid was solved at, and the 1-based BYTE offset
    of that model's record in the binary.  The first line is a header and '#'
    lines are comments.
    """
    rows = []
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            t = s.split()
            if len(t) < 9:
                continue                      # header / free text
            try:
                rows.append(dict(id=t[0].strip("'\""), teff=float(t[1]),
                                 logg=float(t[2]),
                                 feh=float(t[3]), alpha=float(t[4]),
                                 mass=float(t[5]), vturb=float(t[6]),
                                 abund=float(t[7]), pos=int(t[8])))
            except ValueError:
                continue                      # header line, not a record
    if not rows:
        sys.exit(f"ERROR: no usable records in {path}")
    return rows


def pick_model(rows, teff, logg, feh, alpha=None, abund=None):
    """Nearest grid node, normalised so no one axis dominates the distance."""
    def d(r):
        x = ((r["teff"] - teff) / 250.0) ** 2 \
            + ((r["logg"] - logg) / 0.5) ** 2 \
            + ((r["feh"] - feh) / 0.5) ** 2
        if alpha is not None:
            x += ((r["alpha"] - alpha) / 0.2) ** 2
        if abund is not None:
            x += ((r["abund"] - abund) / 0.2) ** 2
        return x
    best = min(rows, key=d)
    return best


# ---------------------------------------------------------------------------
#  Binary grid record
# ---------------------------------------------------------------------------
def read_binary_record(path, pos):
    """One model's record, at the 1-based byte offset `pos` from the index.

    Layout (little-endian stream, as written by the grid's Fortran):
        newline-terminated model id string
        int32   ndep
        int32   nlev
        float64 tau(ndep)                 -- tau_5000, linear
        float64 depart(ndep, nlev)        -- Fortran order, depth fastest
    """
    with open(path, "rb") as f:
        f.seek(pos - 1)
        atmos = f.readline(500).decode("utf-8", "ignore").strip()
        ndep = int.from_bytes(f.read(4), "little")
        nlev = int.from_bytes(f.read(4), "little")
        if not (0 < ndep < 10000 and 0 < nlev < 10000):
            sys.exit(f"ERROR: implausible ndep={ndep} nlev={nlev} at pos={pos}; "
                     "index and binary are probably mismatched.")
        tau = np.fromfile(f, dtype="<f8", count=ndep)
        dep = np.fromfile(f, dtype="<f8", count=ndep * nlev).reshape(nlev, ndep)
    return atmos, ndep, nlev, tau, dep


# ---------------------------------------------------------------------------
#  MARCS model: the bridge from tau_5000 to column mass
# ---------------------------------------------------------------------------
def read_marcs(path):
    """Return (lgTau5, rhox) from a MARCS .mod file.

    lgTau5 comes from the first structure table, RHOX from the second; the two
    are on the same layers and are matched by lgTauR, which both carry.  Parsed
    by locating the header lines rather than by fixed line numbers, since the
    preamble length varies between MARCS releases.
    """
    lines = open(path).read().splitlines()
    i1 = i2 = None
    for i, ln in enumerate(lines):
        if "lgTau5" in ln and "lgTauR" in ln:
            i1 = i
        elif "RHOX" in ln and "lgTauR" in ln:
            i2 = i
    if i1 is None or i2 is None:
        sys.exit(f"ERROR: {path} does not look like a MARCS model "
                 "(no lgTau5 / RHOX structure tables found).")

    def block(start):
        out = []
        for ln in lines[start + 1:]:
            t = ln.split()
            if not t:
                break
            try:
                out.append([float(v) for v in t])
            except ValueError:
                break
        return out

    b1, b2 = block(i1), block(i2)
    if len(b1) != len(b2):
        sys.exit(f"ERROR: {path}: structure tables have {len(b1)} and "
                 f"{len(b2)} rows.")
    # columns: k lgTauR lgTau5 ...   and   k lgTauR KappaRoss ... RHOX
    lgtau5 = np.array([r[2] for r in b1])
    rhox = np.array([r[-1] for r in b2])
    tauR1 = np.array([r[1] for r in b1])
    tauR2 = np.array([r[1] for r in b2])
    if not np.allclose(tauR1, tauR2, atol=1e-6):
        sys.exit(f"ERROR: {path}: the two structure tables are not on the "
                 "same lgTauR layers.")
    return lgtau5, rhox


# ---------------------------------------------------------------------------
#  Model atom -> level indices
# ---------------------------------------------------------------------------
def levels_from_atom(path):
    """Best-effort level identification by term energy.

    MULTI-format atoms list levels as  E[cm^-1]  g  'label'  ion  index.
    Returns a list of (lo, up) 1-based level indices, one per transition, or
    exits if any match is missing or ambiguous.
    """
    # The atom declares NK on the line after the "* NK NLIN ..." comment.
    # Counting level records and finding a different number means the record
    # pattern below is matching something it should not, and every index
    # after that point would be silently wrong.
    nk_declared = None
    lines = open(path).read().splitlines()
    for i, ln in enumerate(lines):
        if ln.strip().startswith("*") and "NK" in ln and "NLIN" in ln:
            try:
                nk_declared = int(lines[i + 1].split()[0])
            except (IndexError, ValueError):
                pass
            break

    energies, idx = [], []
    n = 0
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("*"):
            continue
        t = s.split()
        try:
            e = float(t[0]); float(t[1])
        except (ValueError, IndexError):
            continue
        if "'" not in s and '"' not in s:
            continue                     # level records carry a quoted label
        n += 1
        energies.append(e)
        idx.append(n)
    if not energies:
        sys.exit(f"ERROR: no level records recognised in {path}; "
                 "use --levels to state the mapping explicitly.")
    if nk_declared is not None and n != nk_declared:
        sys.exit(f"ERROR: {path} declares NK = {nk_declared} levels but "
                 f"{n} lines matched the level pattern.  The indices would "
                 "be wrong; use --levels.")
    print(f"model atom      : {path} ({n} levels, NK check "
          f"{'ok' if nk_declared == n else 'not declared'})")
    energies = np.array(energies)

    out = []
    for name, elo, eup in TRANSITIONS:
        pair = []
        for what, target in (("lower", elo), ("upper", eup)):
            hit = np.where(np.abs(energies - target) < ETOL)[0]
            if len(hit) == 0:
                sys.exit(f"ERROR: no level within {ETOL} cm^-1 of "
                         f"{target} ({what} of {name}) in {path}. "
                         "Use --levels.")
            if len(hit) > 1:
                sys.exit(f"ERROR: {len(hit)} levels within {ETOL} cm^-1 of "
                         f"{target} ({what} of {name}) in {path}: indices "
                         f"{[idx[h] for h in hit]}. Ambiguous -- use --levels.")
            pair.append(idx[hit[0]])
        out.append(tuple(pair))
    return out


# ---------------------------------------------------------------------------
#  Sidecar
# ---------------------------------------------------------------------------
def write_dep(path, lrhox, b_lo, b_up, provenance):
    n = len(lrhox)
    ntr = len(b_lo)
    order = np.argsort(lrhox)             # SYNTHE requires increasing depth
    with open(path, "w") as f:
        for ln in provenance:
            f.write(f"# {ln}\n")
        f.write(f"NTRANS {ntr}\nNDEPTH {n}\nDEPTHVAR LOGRHOX\n")
        for i in order:
            vals = []
            for k in range(ntr):
                vals += [b_lo[k][i], b_up[k][i]]
            f.write(f"{lrhox[i]:24.17e} "
                    + " ".join(f"{v:24.17e}" for v in vals) + "\n")


def build(args):
    rows = read_aux(args.aux)
    m = pick_model(rows, args.teff, args.logg, args.feh, args.alpha, args.abund)
    print(f"grid index      : {len(rows)} models in {args.aux}")
    print(f"requested       : Teff {args.teff}  logg {args.logg}  "
          f"[Fe/H] {args.feh}"
          + (f"  [a/Fe] {args.alpha}" if args.alpha is not None else "")
          + (f"  A(X) {args.abund}" if args.abund is not None else ""))
    print(f"selected        : Teff {m['teff']}  logg {m['logg']}  "
          f"[Fe/H] {m['feh']}  [a/Fe] {m['alpha']}  A(X) {m['abund']}")
    print(f"                  {m['id']}  (byte {m['pos']})")

    atmos, ndep, nlev, tau5, dep = read_binary_record(args.bin, m["pos"])
    print(f"binary record   : {ndep} depths x {nlev} levels, id '{atmos}'")
    # Whether the record stores tau or log10(tau) has varied between grid
    # releases.  The two cannot be confused: a linear optical depth is
    # strictly positive, so any negative value means the column is already
    # logarithmic.  Decide on that, and say which was found.
    if tau5.min() > 0 and tau5.max() > 10.0:
        lgtau5_grid = np.log10(tau5)
        print("depth column    : linear tau_5000 -> log10")
    elif tau5.min() < 0:
        lgtau5_grid = tau5
        print("depth column    : already log10(tau_5000)")
    elif tau5.min() > 0:
        lgtau5_grid = np.log10(tau5)
        print("depth column    : linear tau_5000 -> log10 (shallow grid)")
    else:
        sys.exit("ERROR: tau column contains zeros; the record layout or the "
                 "byte offset is wrong.")
    if not np.all(np.diff(lgtau5_grid) > 0):
        sys.exit("ERROR: the record's depth column is not monotonic; the byte "
                 "offset is almost certainly wrong.")

    # levels
    if args.levels:
        v = [int(x) for x in args.levels.split(",")]
        if len(v) != 2 * len(TRANSITIONS):
            sys.exit(f"ERROR: --levels needs {2*len(TRANSITIONS)} indices "
                     f"(lo,up per transition), got {len(v)}.")
        pairs = [(v[2 * k], v[2 * k + 1]) for k in range(len(TRANSITIONS))]
    elif args.atom:
        pairs = levels_from_atom(args.atom)
    else:
        sys.exit("ERROR: give --levels or --atom; the level ordering is a "
                 "property of the grid's model atom and must not be guessed.")
    for (name, _, _), (lo, up) in zip(TRANSITIONS, pairs):
        if not (1 <= lo <= nlev and 1 <= up <= nlev):
            sys.exit(f"ERROR: level indices {lo},{up} outside 1..{nlev}.")
        print(f"levels          : {name}  ->  lower {lo}, upper {up}")
    ups = [p[1] for p in pairs]
    if len(set(ups)) < len(ups):
        print("note            : transitions share an upper level -- this "
              "atom does not resolve the 3p fine structure.")

    # tau_5000 -> column mass, through the atmosphere the grid was solved on
    marcs = args.marcs or os.path.join(args.marcs_dir or ".", m["id"])
    if not os.path.exists(marcs) and not marcs.endswith(".mod"):
        marcs += ".mod"
    if not os.path.exists(marcs):
        sys.exit(f"ERROR: MARCS model not found: {marcs}\n"
                 "       It is required: it is the only thing that ties the "
                 "grid's tau_5000 to a column mass.")
    lgtau5_m, rhox_m = read_marcs(marcs)
    print(f"MARCS model     : {marcs}  ({len(rhox_m)} layers)")

    if len(lgtau5_m) == ndep and np.allclose(lgtau5_m, lgtau5_grid, atol=1e-3):
        lrhox = np.log10(rhox_m)          # same layers: direct association
        print("depth mapping   : grid and MARCS layers coincide (direct)")
    else:
        lrhox = np.interp(lgtau5_grid, lgtau5_m, np.log10(rhox_m))
        print(f"depth mapping   : interpolated ({ndep} grid layers onto "
              f"{len(rhox_m)} MARCS layers in lgTau5)")
        if lgtau5_grid.min() < lgtau5_m.min() - 1e-6 or \
           lgtau5_grid.max() > lgtau5_m.max() + 1e-6:
            print("WARNING         : grid tau range exceeds the MARCS model's; "
                  "edge layers are clamped.")

    b_lo = [dep[lo - 1] for lo, _ in pairs]
    b_up = [dep[up - 1] for _, up in pairs]
    for k, (name, _, _) in enumerate(TRANSITIONS):
        print(f"b range         : {name}: b_lo {b_lo[k].min():.4f}"
              f"..{b_lo[k].max():.4f}   b_up {b_up[k].min():.4f}"
              f"..{b_up[k].max():.4f}")
    for k in range(len(TRANSITIONS)):
        if b_lo[k].min() <= 0 or b_up[k].min() <= 0:
            sys.exit(f"ERROR: non-positive departure coefficient for "
                     f"transition {k+1}; SYNTHE interpolates log b and will "
                     "refuse this.")

    prov = [
        "SYNTHE NLTE departure-coefficient sidecar (NLTE_MODE = 2)",
        f"built by tools/nlte_make_dep.py from {os.path.basename(args.bin)}",
        f"grid model  : {m['id']}",
        f"grid params : Teff {m['teff']} logg {m['logg']} [Fe/H] {m['feh']} "
        f"[a/Fe] {m['alpha']} A(X) {m['abund']}",
        f"requested   : Teff {args.teff} logg {args.logg} [Fe/H] {args.feh}",
        f"levels      : " + "  ".join(f"tr{k+1} lo={p[0]} up={p[1]}"
                                      for k, p in enumerate(pairs)),
        f"depth       : log10 column mass [g/cm^2], from lgTau5 via "
        f"{os.path.basename(marcs)}",
        "transition order matches NLTE_TR_* in src/mod_mklinelist.f90",
    ]
    write_dep(args.out, lrhox, b_lo, b_up, prov)
    print(f"wrote           : {args.out}  ({ndep} depths, "
          f"log10 rhox {lrhox.min():.4f} .. {lrhox.max():.4f})")


# ---------------------------------------------------------------------------
#  Self-test: fabricate a grid in the documented layout and round-trip it
# ---------------------------------------------------------------------------
def selftest(tmpdir):
    os.makedirs(tmpdir, exist_ok=True)
    binf = os.path.join(tmpdir, "grid.bin")
    auxf = os.path.join(tmpdir, "aux.txt")
    modf = os.path.join(tmpdir, "testmodel.mod")
    outf = os.path.join(tmpdir, "out.dep")

    ndep, nlev = 56, 5
    lgtau5 = np.linspace(-5.0, 2.0, ndep)
    tau5 = 10.0 ** lgtau5
    lrhox = np.linspace(-3.4, 1.1, ndep)          # the answer we must recover
    rng = np.random.default_rng(12345)
    dep = np.empty((nlev, ndep))
    for k in range(nlev):
        dep[k] = 1.0 + 0.4 * (k + 1) * np.exp(-((lgtau5 + 3.0) / 1.5) ** 2)

    # two models, so selection has something to choose between
    recs = []
    with open(binf, "wb") as f:
        f.write(b"#" + b" " * 998 + b"\n")         # NLTEgrid_header
        for tag, teff, logg, feh in (("testmodel", 5777.0, 4.44, 0.00),
                                     ("othermodel", 4000.0, 4.50, -1.00)):
            pos = f.tell() + 1                     # 1-based, as in the index
            f.write((tag + "\n").encode())
            f.write(struct.pack("<i", ndep))
            f.write(struct.pack("<i", nlev))
            f.write(tau5.astype("<f8").tobytes())
            scale = 1.0 if tag == "testmodel" else 2.0
            f.write((dep * scale).astype("<f8").tobytes(order="C"))
            recs.append((tag, teff, logg, feh, pos))

    with open(auxf, "w") as f:
        f.write("# id teff logg feh alpha mass vturb abund pos\n")
        for tag, teff, logg, feh, pos in recs:
            f.write(f"{tag} {teff} {logg} {feh} 0.0 1.0 1.0 6.17 {pos}\n")

    with open(modf, "w") as f:
        f.write("fabricated MARCS-like model for the self-test\n")
        f.write("Model structure\n")
        f.write(" k lgTauR   lgTau5    Depth     T        Pe  Pg  Prad  Pturb\n")
        for i in range(ndep):
            f.write(f"{i+1:3d} {lgtau5[i]:8.4f} {lgtau5[i]:8.4f} 0.0 5000.0 "
                    f"1.0 1.0 1.0 0.0\n")
        f.write("\n")
        f.write(" k lgTauR    KappaRoss   Density   Mu   Vconv  Fconv/F   RHOX\n")
        for i in range(ndep):
            f.write(f"{i+1:3d} {lgtau5[i]:8.4f} 1.0 1.0 1.0 0.0 0.0 "
                    f"{10.0**lrhox[i]:12.5E}\n")

    args = argparse.Namespace(out=outf, bin=binf, aux=auxf, marcs=modf,
                              marcs_dir=None, teff=5777, logg=4.44, feh=0.0,
                              alpha=None, abund=None, levels="1,3,1,2",
                              atom=None)
    print("=" * 70)
    build(args)
    print("=" * 70)

    # --- verify -----------------------------------------------------------
    rows = [l.split() for l in open(outf)
            if l.strip() and not l.startswith("#")
            and l.split()[0] not in ("NTRANS", "NDEPTH", "DEPTHVAR")]
    got = np.array([[float(v) for v in r] for r in rows])
    ok = True

    # The depth tolerance is set by the MARCS file, not by this tool: RHOX is
    # tabulated with ~6 significant digits there (and in the fabricated model
    # here), so log10(rhox) can only be recovered to ~1e-6 dex.  That is a
    # 4e-6 relative error in column mass -- irrelevant physically, but it is
    # the real floor, so the test states it rather than pretending to 1e-12.
    e = np.abs(got[:, 0] - lrhox).max()
    print(f"depth column recovered   : max |dlog10 rhox| = {e:.3e}"
          f"  (floor ~1e-6, set by MARCS RHOX precision)"
          f"  [{'PASS' if e < 1e-5 else 'FAIL'}]")
    ok &= e < 1e-5

    for k, (lo, up) in enumerate(((1, 3), (1, 2))):
        e1 = np.abs(got[:, 1 + 2 * k] - dep[lo - 1]).max()
        e2 = np.abs(got[:, 2 + 2 * k] - dep[up - 1]).max()
        print(f"transition {k+1} b_lo=level{lo} b_up=level{up} : "
              f"max |db| = {max(e1,e2):.3e}"
              f"  [{'PASS' if max(e1,e2) < 1e-12 else 'FAIL'}]")
        ok &= max(e1, e2) < 1e-12

    inc = np.all(np.diff(got[:, 0]) > 0)
    print(f"depth column increasing  : {inc}  [{'PASS' if inc else 'FAIL'}]")
    ok &= inc

    hdr = open(outf).read()
    hv = "DEPTHVAR LOGRHOX" in hdr and "NTRANS 2" in hdr
    print(f"header written           : {hv}  [{'PASS' if hv else 'FAIL'}]")
    ok &= hv

    # the OTHER model must not have been picked (its b are 2x)
    picked_right = np.abs(got[:, 1] - dep[0]).max() < 1e-12
    print(f"correct model selected   : {picked_right}"
          f"  [{'PASS' if picked_right else 'FAIL'}]")
    ok &= picked_right

    print("=" * 70)
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", nargs="?", help="output .dep path")
    ap.add_argument("--bin", help="NLTE grid binary")
    ap.add_argument("--aux", help="grid index (auxData_*) file")
    ap.add_argument("--marcs", help="MARCS .mod for the selected grid model")
    ap.add_argument("--marcs-dir", help="directory to find it in by name")
    ap.add_argument("--teff", type=float)
    ap.add_argument("--logg", type=float)
    ap.add_argument("--feh", type=float)
    ap.add_argument("--alpha", type=float, default=None)
    ap.add_argument("--abund", type=float, default=None,
                    help="element abundance the grid node was solved at")
    ap.add_argument("--levels", help="lo,up per transition, e.g. 1,3,1,2")
    ap.add_argument("--atom", help="model atom file, to find levels by energy")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--selftest-dir", default="/tmp/nlte_selftest")
    a = ap.parse_args()

    if a.selftest:
        sys.exit(selftest(a.selftest_dir))
    missing = [n for n in ("out", "bin", "aux", "teff", "logg", "feh")
               if getattr(a, n) is None]
    if missing:
        ap.error("missing: " + ", ".join(missing))
    build(a)


if __name__ == "__main__":
    main()
