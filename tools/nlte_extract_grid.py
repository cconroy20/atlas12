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
      [--url URL | --zip LOCAL.zip] [--levels 10] [--element Na] \\
      [--atom atom.na_qmh]

The grid layout is the same for every element in the collection -- only ndep,
nlev and the record stride differ -- so the same pass serves all of them; give
--nlev the element's level count and --levels how many to keep.
"""
import argparse
import json
import os
import re
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
    gaps = np.diff(ptr)
    if len(gaps) and gaps.min() <= 0:
        sys.exit("ERROR: duplicate or non-monotonic pointers in the index")
    vals, cnt = np.unique(gaps, return_counts=True)
    stride = int(vals[cnt.argmax()])          # the record size itself
    # Records are fixed-stride, but the index does not necessarily list every
    # record in the file: Mg's index skips some, so its gaps are multiples of
    # the stride.  Take the modal gap as the record size and let the caller
    # drive the pass off the pointers themselves, which covers both cases.
    extra = int(cnt.sum() - cnt.max())
    if extra:
        bad = int(((vals % stride) != 0).sum())
        if bad:
            sys.exit(f"ERROR: {bad} pointer gaps are not multiples of the "
                     f"record size {stride}; the layout is not what this "
                     "extractor assumes")
        print(f"note      : index skips records -- {extra} of {len(gaps)} "
              f"gaps exceed one record ({int(ptr[-1] - ptr[0]) // stride + 1} "
              f"records span the file, {len(rows)} are indexed)")
    return rows, int(ptr[0]), stride


def read_model_atom(path, keep):
    """The first `keep` levels of a MULTI model atom: energy, g, label, ion.

    The grid stores b by level INDEX, so turning "Ca II H" into a pair of
    indices needs the atom the grid was solved with.  Recording it alongside
    the store keeps that mapping with the data instead of in someone's head.
    The level block starts after the line carrying nlev/nline/ncont and runs
    to the first line that does not parse as `energy g 'label' ion`; comment
    lines (leading *) are skipped, as MULTI allows them anywhere.
    """
    lev, started = [], False
    for ln in open(path, errors="replace"):
        s = ln.strip()
        if not s or s.startswith("*"):
            continue
        t = s.split()
        if not started:
            # the counts line: >=3 integers, the first being nlev
            if len(t) >= 3 and all(x.lstrip("+-").isdigit() for x in t[:3]):
                started = True
            continue
        m = re.match(r"\s*([-\d.eEdD+]+)\s+([-\d.eEdD+]+)\s+'([^']*)'\s+(\d+)", ln)
        if not m:
            break
        lev.append({"index": len(lev) + 1,
                    "energy_cm": float(m.group(1).replace("D", "E")),
                    "g": float(m.group(2).replace("D", "E")),
                    "label": m.group(3).strip(),
                    "ion": int(m.group(4))})
        if len(lev) >= keep:
            break
    return lev


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
    ap.add_argument("--element", default="", help="element symbol, recorded "
                    "in the json (the store itself is element-agnostic)")
    ap.add_argument("--atom", help="MULTI model atom the grid was solved with;"
                    " its level table is copied into the json, which is what "
                    "later maps a transition's energies onto level indices")
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
    # Byte to start each record at, from the index rather than from a running
    # count: the two agree only when every record is indexed.
    ptrs = np.array([r[0] for r in rows], dtype=np.int64) - 1
    skips = ptrs - np.concatenate(([0], ptrs[:-1] + need))
    if skips.min() < 0:
        sys.exit("ERROR: records overlap given the bytes kept per record")

    print(f"records   : {nrec}  stride {stride} B, header {hdr} B, "
          f"{int(skips[1:].sum())} B skipped between records")
    print(f"keeping   : tau({a.ndep}) + b({keep} of {a.nlev} levels) "
          f"= {need - a.idlen - 8} B/record of {stride}")
    print(f"output    : {nrec * (a.ndep + keep * a.ndep) * 4 / 1e6:.0f} MB "
          f"float32", flush=True)

    fh, name = open_stream(a.url, a.zip)
    print(f"member    : {name}", flush=True)

    d = zlib.decompressobj(-15)
    fout = open(a.out + ".f32", "wb")

    mode, remaining = "skip", int(skips[0])
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
                mode = "skip"
                remaining = int(skips[n]) if n < nrec else 0
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
        "element": a.element,
        "model_atom": os.path.basename(a.atom) if a.atom else None,
        "levels": read_model_atom(a.atom, keep) if a.atom else None,
    }, open(a.out + ".json", "w"), indent=2)

    el = time.time() - t0
    print(f"done: {n} records in {el/60:.1f} min, {cpos/1e9:.2f} GB streamed",
          flush=True)
    print(f"  {a.out}.f32  {os.path.getsize(a.out + '.f32')/1e6:.0f} MB")
    print(f"  {a.out}.npz  {os.path.getsize(a.out + '.npz')/1e6:.1f} MB")


if __name__ == "__main__":
    main()
