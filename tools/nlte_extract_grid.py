#!/usr/bin/env python3
"""One-pass extraction of the low-lying levels from a published NLTE grid.

WHY THIS EXISTS
---------------
The Turbospectrum-packaged NLTE grids (Amarsi et al. 2020; Gerber et al. 2023)
ship as a single ZIP64 deflate member -- the Na I 1D grid is 15.9 GB
compressed, 57.1 GB raw, 436255 model records.  Deflate is not seekable, so
reaching a record at uncompressed byte N means decompressing bytes 0..N.
Pulling one record therefore costs a large fraction of the whole file (a
record 44 GB in cost 12.3 GB of transfer), which is fine once and absurd per
model.

This does it ONCE.  It streams the download straight into a decompressor and
keeps only the part of each record that is ever used -- the depth scale and
the first few atomic levels -- discarding the rest without ever storing it.
The result is a small local file with random access, after which the marginal
cost per model is zero.

WHAT MAKES IT CHEAP
-------------------
Two properties of the format, both verified against the index before use:

  * Records are FIXED STRIDE.  All 436254 pointer gaps in auxData are the same
    value, so record i begins at HDR + i*RECSZ and no index lookups are needed
    while streaming; this is a state machine, not a search.

  * The levels we want are at the FRONT of each record.  The departure block is
    Fortran-order (ndep, nlev), so level k is a contiguous run at offset
    k*ndep in the flat stream.  Levels 1..NLEV_KEEP are therefore the first
    NLEV_KEEP runs -- a head read, no seeking inside the record.

Keeping 10 of Na's 290 levels costs 2464 B per record against 130876, a 53x
reduction, and covers every Na I feature likely to be wanted:

    D2 5891.6 / D1 5897.6 A   3s - 3p      levels 1,2,3
    8183 / 8195 A             3p - 3d      levels 2,3,5,6
    1.14 um                   3p - 4s      levels 2,3,4

Choose the level count generously: redoing the pass costs another full
download, whereas storage is linear and small.

OUTPUT
------
  <out>.f32    436255 records x (ndep tau + NLEV_KEEP*ndep b), float32, in
               stream order, which is pointer order
  <out>.npz    per-record Teff, logg, [Fe/H], [alpha/Fe], vturb, A(X), taken
               from the index and ordered to match
  <out>.json   shapes, provenance and the level->transition map

float32 is deliberate: b is a population ratio wanted to ~4 digits, and the
grid's own convergence is far coarser than float32 resolution.

USAGE
-----
  python3 tools/nlte_extract_grid.py OUT_PREFIX --aux auxData_Na.dat \\
      [--url URL | --zip LOCAL.zip] [--levels 10]
"""
import argparse
import json
import os
import sys
import time
import urllib.request
import zlib

import numpy as np


def read_index(path):
    """Index rows sorted by pointer; also the stride check the pass relies on."""
    rows = []
    for ln in open(path):
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        t = s.split()
        if len(t) < 9:
            continue
        try:
            rows.append((int(t[8]), t[0].strip("'\""), float(t[1]), float(t[2]),
                         float(t[3]), float(t[4]), float(t[6]), float(t[7])))
        except ValueError:
            continue
    if not rows:
        sys.exit(f"ERROR: no usable records in {path}")
    rows.sort()
    ptr = np.array([r[0] for r in rows], dtype=np.int64)
    gaps = np.unique(np.diff(ptr))
    if len(gaps) != 1:
        sys.exit(f"ERROR: records are not fixed-stride ({len(gaps)} distinct "
                 "gaps).  This extractor assumes they are; a variable-length "
                 "grid needs the pointer list threaded through the pass.")
    return rows, int(ptr[0]), int(gaps[0])


def open_stream(url=None, zippath=None):
    """Return a file-like positioned at the start of the deflate stream."""
    if zippath:
        fh = open(zippath, "rb")
    else:
        req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
        fh = urllib.request.urlopen(req, timeout=120)
    head = fh.read(30)
    if head[:4] != b"PK\x03\x04":
        sys.exit("ERROR: not a zip local file header")
    method = int.from_bytes(head[8:10], "little")
    nlen = int.from_bytes(head[26:28], "little")
    elen = int.from_bytes(head[28:30], "little")
    name = fh.read(nlen).decode()
    fh.read(elen)
    if method != 8:
        sys.exit(f"ERROR: expected deflate, got compression method {method}")
    return fh, name


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", help="output prefix")
    ap.add_argument("--aux", required=True, help="auxData index file")
    ap.add_argument("--url", help="grid zip URL to stream")
    ap.add_argument("--zip", help="local grid zip instead of a URL")
    ap.add_argument("--levels", type=int, default=10,
                    help="number of low-lying levels to keep (default 10)")
    ap.add_argument("--ndep", type=int, default=56, help="depths per record")
    ap.add_argument("--nlev", type=int, default=290, help="levels per record")
    ap.add_argument("--idlen", type=int, default=500, help="id field bytes")
    a = ap.parse_args()
    if not (a.url or a.zip):
        ap.error("give --url or --zip")

    rows, ptr0, stride = read_index(a.aux)
    nrec = len(rows)
    hdr = ptr0 - 1                       # bytes before the first record
    keep = a.levels

    # Byte layout within one record.
    need = a.idlen + 8 + a.ndep * 8 + keep * a.ndep * 8
    expect = a.idlen + 8 + a.ndep * 8 + a.nlev * a.ndep * 8
    if expect != stride:
        sys.exit(f"ERROR: index stride {stride} != layout {expect} implied by "
                 f"idlen={a.idlen} ndep={a.ndep} nlev={a.nlev}")
    skip_rest = stride - need

    print(f"records   : {nrec}  stride {stride} B, header {hdr} B")
    print(f"keeping   : tau({a.ndep}) + b({keep} of {a.nlev} levels) "
          f"= {need - a.idlen - 8} B/record of {stride}")
    print(f"output    : {nrec * (a.ndep + keep * a.ndep) * 4 / 1e6:.0f} MB "
          f"float32", flush=True)

    fh, name = open_stream(a.url, a.zip)
    print(f"member    : {name}", flush=True)

    d = zlib.decompressobj(-15)
    fout = open(a.out + ".f32", "wb")

    mode, remaining = "skip", hdr
    blk = bytearray()
    n = 0
    bad = 0
    cpos = 0
    t0 = time.time()
    nextrep = 1 << 30

    while n < nrec:
        chunk = fh.read(1 << 23)
        if not chunk:
            break
        cpos += len(chunk)
        out = d.decompress(chunk)
        pos, L = 0, len(out)
        while pos < L and n < nrec:
            take = min(remaining, L - pos)
            if mode == "collect":
                blk += out[pos:pos + take]
            pos += take
            remaining -= take
            if remaining:
                continue
            if mode == "skip":
                mode, remaining = "collect", need
                blk = bytearray()
            else:
                nd = int.from_bytes(blk[a.idlen:a.idlen + 4], "little")
                nl = int.from_bytes(blk[a.idlen + 4:a.idlen + 8], "little")
                if nd != a.ndep or nl != a.nlev:
                    bad += 1
                o = a.idlen + 8
                tau = np.frombuffer(blk, "<f8", a.ndep, o)
                b = np.frombuffer(blk, "<f8", keep * a.ndep, o + a.ndep * 8)
                fout.write(tau.astype("<f4").tobytes())
                fout.write(b.astype("<f4").tobytes())
                n += 1
                mode, remaining = "skip", skip_rest
        if cpos >= nextrep:
            el = time.time() - t0
            print(f"  {cpos/1e9:6.2f} GB in | {n:>7}/{nrec} records | "
                  f"{cpos/el/1e6:5.1f} MB/s | {el/60:5.1f} min", flush=True)
            nextrep += 1 << 30

    fout.close()
    try:
        fh.close()
    except Exception:
        pass

    if bad:
        print(f"WARNING: {bad} records had unexpected ndep/nlev", flush=True)
    if n != nrec:
        print(f"WARNING: extracted {n} of {nrec} records", flush=True)

    np.savez_compressed(
        a.out + ".npz",
        teff=np.array([r[2] for r in rows[:n]], dtype=np.float32),
        logg=np.array([r[3] for r in rows[:n]], dtype=np.float32),
        feh=np.array([r[4] for r in rows[:n]], dtype=np.float32),
        alpha=np.array([r[5] for r in rows[:n]], dtype=np.float32),
        vturb=np.array([r[6] for r in rows[:n]], dtype=np.float32),
        abund=np.array([r[7] for r in rows[:n]], dtype=np.float32),
        model=np.array([r[1] for r in rows[:n]]))
    json.dump({
        "nrecords": n, "ndep": a.ndep, "levels_kept": keep,
        "nlev_in_grid": a.nlev, "dtype": "float32",
        "record_layout": "tau(ndep) then b(levels_kept, ndep), C order",
        "depth_variable": "tau_5000 (linear, as stored in the grid)",
        "source_member": name,
        "index": os.path.basename(a.aux),
        "na_levels": {"1": "3s 2S", "2": "3p 2P1/2", "3": "3p 2P3/2",
                      "4": "4s 2S", "5": "3d 2D", "6": "3d 2D"},
        "na_transitions": {"D2 5891.6": [1, 3], "D1 5897.6": [1, 2],
                           "8183": [2, 5], "8195": [3, 5],
                           "11385": [2, 4], "11407": [3, 4]},
    }, open(a.out + ".json", "w"), indent=2)

    el = time.time() - t0
    print(f"done: {n} records in {el/60:.1f} min, {cpos/1e9:.2f} GB streamed",
          flush=True)
    print(f"  {a.out}.f32  {os.path.getsize(a.out + '.f32')/1e6:.0f} MB")
    print(f"  {a.out}.npz  {os.path.getsize(a.out + '.npz')/1e6:.1f} MB")


if __name__ == "__main__":
    main()
