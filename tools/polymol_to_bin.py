#!/usr/bin/env python3
"""Pack Kurucz poly-ascii super-line files into the 8-byte molbin format.

The polyatomic super-line lists (CaOH etc.) are formatted-ascii and their
Fortran formatted read dominates cool-star SYNTHE startup (15-20+ min for
the ~54M-line CaOH pair).  This packs them into the same 8-byte record
layout as h2opokazatel.bin, read by mod_mklinelist.read_molec_poly_bin in
seconds via direct access:

  record 1 (header):  int32 = -1 sentinel, int32 = molecule code (e.g. 10820)
  records 2..N:       int32 iwl  (wl_vac[nm] = exp(iwl * ln(1 + 1/2e6)))
                      packed int32 = (igflog << 16) | ielo
                        ielo   = round(E_low [cm^-1]), int16 > 0
                        igflog = round(1000 * log10 gf) + 16384, int16 > 0

Conventions preserved from read_molec_poly: wl_vac from |E_up| - |E_low|
(NOT the wl column), E_low = min(|e|, |ep|), gf as filed (bare).  Dropped
(negligible for blended super-line haze): per-line loggr (reader uses the
classical radiative default) and the X-X damping branch (reader assumes
electronic-band defaults gammas=3e-5, gammaw=1e-7).

Usage:
  python3 polymol_to_bin.py OUT.bin CODE IN1.dat [IN2.dat ...]
  python3 polymol_to_bin.py data/mol/caoh_oyt6_super.bin 10820 \\
      data/mol/caoh_oyt6_super_blue.dat data/mol/caoh_oyt6_super.dat
"""
import os
import sys

import numpy as np

RATIOLOG = np.log(1.0 + 1.0 / 2.0e6)
LINE_LEN = 77          # 76 chars + newline


def parse_ascii(path):
    print(f"parsing {path} ...", flush=True)
    raw = open(path, "rb").read()
    if len(raw) % LINE_LEN:
        # tolerate a missing trailing newline
        if (len(raw) + 1) % LINE_LEN == 0:
            raw += b"\n"
        else:
            sys.exit(f"{path}: not fixed-width {LINE_LEN - 1}+1 records")
    arr = np.frombuffer(raw, dtype=f"S{LINE_LEN}")
    b = arr.view(np.uint8).reshape(-1, LINE_LEN)

    def fcol(lo, hi):
        s = b[:, lo:hi].tobytes().decode("ascii")
        n = hi - lo
        return np.array(
            np.frombuffer(b[:, lo:hi].tobytes(), dtype=f"S{n}"), dtype=float)

    gflog = fcol(10, 17)
    e = np.abs(fcol(22, 32))
    ep = np.abs(fcol(37, 48))
    nu = np.abs(ep - e)
    elo = np.minimum(e, ep)
    return nu, gflog, elo


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    out, code = sys.argv[1], int(sys.argv[2])
    nus, gfs, els = [], [], []
    for path in sys.argv[3:]:
        nu, gflog, elo = parse_ascii(path)
        nus.append(nu); gfs.append(gflog); els.append(elo)
    nu = np.concatenate(nus)
    gflog = np.concatenate(gfs)
    elo = np.concatenate(els)
    print(f"{len(nu)} lines total")

    wl = 1.0e7 / nu                                   # vacuum nm
    iwl = np.round(np.log(wl) / RATIOLOG).astype(np.int64)
    igf = np.round(1000.0 * gflog).astype(np.int64) + 16384
    ielo = np.clip(np.round(elo).astype(np.int64), 1, 32767)
    ok = (igf >= 1) & (igf <= 32767) & np.isfinite(wl)
    ncut = int((~ok).sum())
    iwl, igf, ielo = iwl[ok], igf[ok], ielo[ok]
    order = np.argsort(iwl, kind="stable")
    iwl, igf, ielo = iwl[order], igf[order], ielo[order]

    rec = np.empty((len(iwl) + 1, 2), dtype="<i4")
    rec[0, 0] = -1
    rec[0, 1] = code
    rec[1:, 0] = iwl
    rec[1:, 1] = (igf.astype(np.int64) << 16) | ielo
    rec.tofile(out)
    print(f"wrote {out}: {len(iwl)} lines (+header), "
          f"{os.path.getsize(out)/1e6:.1f} MB, {ncut} dropped (gf range); "
          f"wl {np.exp(iwl[0]*RATIOLOG):.1f}-{np.exp(iwl[-1]*RATIOLOG):.1f} nm")


if __name__ == "__main__":
    main()
