#!/usr/bin/env python3
"""Convert fresh ExoMol TOTO (2024-update) to the Kurucz MOL2 binary format.

Phases:
  --parse ISO    decompress + parse one isotopologue's trans, save full npz
  --calibrate    strength-threshold study on the parsed 48 npz
  --build THR    apply threshold, merge all 5, sort, write tiototo2024.bin

Conventions (matching data/mol/tiototo.bin and read_molec_bin):
  32-byte records: wl[nm] f32, log10(gf) f32, E f32, E' f32 [cm-1],
  icode=822 i32, iso=46..50 i32, loggr=100*log10(1/tau_upper) i32, lab=0 i32.
  gf is bare per-isotopologue, g_ns divided out (47: 6, 49: 8), abundance
  weighting applied later by molec_dispatch.  lab=0 everywhere (Kurucz's
  lab=1 marked weak high-E variational cross-transitions; broadening choice
  for those is irrelevant).  Header: int32 magic 0x4D4F4C32 ('MOL2'), nrec,
  6 zero pads.
"""
import argparse
import bz2
import os
import subprocess
import sys

import numpy as np

SCR = os.environ.get("TOTO_WORKDIR", os.getcwd())  # states/trans/pf + outputs live here
GNS = {46: 1.0, 47: 6.0, 48: 1.0, 49: 8.0, 50: 1.0}
C2 = 1.4387769
TCAL = 4000.0          # Boltzmann T for the strength cut
GF_CONST = 1.4992e-16

def load_states(iso):
    sid, sE, sg, stau = [], [], [], []
    with bz2.open(f"{SCR}/toto{iso}.states.bz2", "rt") as fh:
        for ln in fh:
            p = ln.split()
            sid.append(int(p[0])); sE.append(float(p[1]))
            sg.append(float(p[2])); stau.append(float(p[5]))
    sid = np.array(sid)
    n = sid.max() + 1
    E = np.zeros(n); G = np.zeros(n); TAU = np.zeros(n)
    E[sid] = sE; G[sid] = sg; TAU[sid] = stau
    return E, G, TAU

def parse_int_cols(b, lo, hi):
    """Vectorized right-justified integer parse of byte columns [lo:hi)."""
    v = np.zeros(len(b), dtype=np.int64)
    for c in range(lo, hi):
        d = b[:, c].astype(np.int64) - 48
        isdig = (d >= 0) & (d <= 9)
        v = v * np.where(isdig, 10, 1) + np.where(isdig, d, 0)
    return v

def parse_trans(iso):
    raw = f"{SCR}/toto{iso}.trans.raw"
    if not os.path.exists(raw):
        with open(raw, "wb") as out:
            subprocess.run(["bzcat", f"{SCR}/toto{iso}.trans.bz2"],
                           stdout=out, check=True)
    sz = os.path.getsize(raw)
    assert sz % 37 == 0, f"unexpected trans size {sz}"
    nrec = sz // 37
    E, G, TAU = load_states(iso)
    out = {k: [] for k in ("wl", "gf", "el", "eu", "loggr")}
    CH = 20_000_000
    mm = np.memmap(raw, dtype=np.uint8, mode="r").reshape(nrec, 37)
    for k0 in range(0, nrec, CH):
        b = np.asarray(mm[k0:k0+CH])
        u = parse_int_cols(b, 0, 12)
        l = parse_int_cols(b, 12, 25)
        # A field: cols 26..35 = m.mmmm'e'sdd
        mant = (b[:, 26] - 48).astype(np.float64)
        scale = 0.1
        for c in range(28, 32):
            mant += (b[:, c] - 48) * scale
            scale *= 0.1
        expo = (b[:, 34] - 48).astype(np.int64) * 10 + (b[:, 35] - 48)
        sign = np.where(b[:, 33] == 45, -1, 1)
        A = mant * 10.0 ** (sign * expo)
        nu = E[u] - E[l]
        good = nu > 1.0e-6
        u, l, A, nu = u[good], l[good], A[good], nu[good]
        wl_nm = 1.0e7 / nu
        gf = GF_CONST * (G[u] / GNS[iso]) * A * (wl_nm * 10.0) ** 2
        tau = TAU[u]
        with np.errstate(divide="ignore"):
            lg = np.where(tau > 0, np.round(100.0 * np.log10(1.0 / np.maximum(tau, 1e-30))), 0.0)
        out["wl"].append(wl_nm.astype(np.float32))
        out["gf"].append(gf.astype(np.float32))
        out["el"].append(E[l].astype(np.float32))
        out["eu"].append(E[u].astype(np.float32))
        out["loggr"].append(np.clip(lg, -999, 999).astype(np.int16))
        print(f"  iso {iso}: {min(k0+CH, nrec)}/{nrec}", flush=True)
    arrs = {k: np.concatenate(v) for k, v in out.items()}
    np.savez(f"{SCR}/toto{iso}_parsed.npz", **arrs)
    os.remove(raw)
    print(f"iso {iso}: saved {len(arrs['wl'])} lines")

def calibrate():
    d = np.load(f"{SCR}/toto48_parsed.npz")
    S = d["gf"].astype(np.float64) * np.exp(-C2 * d["el"] / TCAL)
    tot = S.sum()
    print(f"48TiO: {len(S)} lines, sum S(4000K) = {tot:.4e}")
    for thr in [1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14]:
        keep = S > thr
        print(f"  thr {thr:.0e}: keep {keep.sum():>10d} "
              f"({100*keep.sum()/len(S):5.1f}%)  S kept {S[keep].sum()/tot:.6f}")

def build(thr):
    recs = []
    for iso in [46, 47, 48, 49, 50]:
        d = np.load(f"{SCR}/toto{iso}_parsed.npz")
        S = d["gf"].astype(np.float64) * np.exp(-C2 * d["el"] / TCAL)
        keep = S > thr
        n = keep.sum()
        r = np.zeros(n, dtype=np.dtype([("wl","<f4"),("gf","<f4"),("e","<f4"),
             ("ep","<f4"),("icode","<i4"),("iso","<i4"),("loggr","<i4"),
             ("lab","<i4")]))
        with np.errstate(divide="ignore"):
            r["gf"] = np.log10(d["gf"][keep])
        r["wl"] = d["wl"][keep]; r["e"] = d["el"][keep]; r["ep"] = d["eu"][keep]
        r["icode"] = 822; r["iso"] = iso; r["loggr"] = d["loggr"][keep]
        recs.append(r)
        print(f"iso {iso}: kept {n} of {len(S)}")
    allr = np.concatenate(recs)
    order = np.argsort(allr["wl"], kind="stable")
    allr = allr[order]
    hdr = np.zeros(8, dtype="<i4")
    hdr[0] = 0x4D4F4C32
    hdr[1] = len(allr)
    with open(f"{SCR}/tiototo2024.bin", "wb") as f:
        hdr.tofile(f)
        allr.tofile(f)
    print(f"wrote tiototo2024.bin: {len(allr)} records, "
          f"{(len(allr)*32+32)/1e9:.2f} GB, wl {allr['wl'][0]:.1f}-"
          f"{allr['wl'][-1]:.1f} nm")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parse", type=int)
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--build", type=float)
    a = ap.parse_args()
    if a.parse: parse_trans(a.parse)
    elif a.calibrate: calibrate()
    elif a.build is not None: build(a.build)
