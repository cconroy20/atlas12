#!/usr/bin/env python3
"""Verification battery for ROSSTAB interpolation modes.

Reads a cache dumped by the Fortran (set ATLAS_ROSSTAB_DUMP=<file>),
re-implements the three lookup modes exactly, and checks:

  R1  Derivative texture: d(log kappa)/d(log T) along a fine T-scan at
      fixed P through the deep-CZ region, per mode.  Reports the RMS of
      the second difference (washboard metric) -- MLS should be far
      smoother than bilinear/Shepard.
  R2  Bubble-probe consistency: the CONVEC probe pattern
      kappa(T +/- dT)/kappa(T) with dT = 0.002*T, swept in T; reports
      the scatter of the implied d(ln kappa)/d(ln T) (the quantity that
      feeds the radiative-loss parameter D).
  R3  Leave-one-out residuals (MLS): predict each stored sample from
      the others; residual scatter should be at the RT-pass noise level
      (~1%) and unbiased.
  R4  Agreement: MLS vs Shepard mean offset over the scan (should be
      ~0 -- same surface, different smoothness).

Usage: python3 rosstab_verify.py [cache_file]
"""
import math
import sys

KFIT = 12
SHEP_EPS = 1e-4


def load(fn):
    with open(fn) as f:
        first = f.readline().split()
        n = int(first[0])
        zt, st, zp, sp = map(float, first[1:5])
        pts = []
        for line in f:
            p = line.split()
            if len(p) == 3:
                pts.append((float(p[0]), float(p[1]), float(p[2])))
    assert len(pts) == n, (len(pts), n)
    return pts, (zt, st, zp, sp)


def knn(pts, qt, qp, k=KFIT, skip=-1):
    d = []
    for i, (t, p, r) in enumerate(pts):
        if i == skip:
            continue
        d.append(((t - qt) ** 2 + (p - qp) ** 2, i))
    d.sort()
    return d[:k]


def shepard(pts, qt, qp, skip=-1):
    nb = knn(pts, qt, qp, skip=skip)
    num = den = 0.0
    for d2, i in nb:
        w = 1.0 / (d2 + SHEP_EPS)
        num += w * pts[i][2]
        den += w
    return num / den


def bilinear(pts, qt, qp):
    # quadrant nearest neighbors
    best = {}
    for i, (t, p, r) in enumerate(pts):
        dt, dp = t - qt, p - qp
        q = (dt >= 0, dp >= 0)
        d2 = dt * dt + dp * dp
        if q not in best or d2 < best[q][0]:
            best[q] = (d2, i)
    if len(best) == 4:
        (tpp, ppp, rpp) = pts[best[(True, True)][1]]
        (tpm, ppm, rpm) = pts[best[(True, False)][1]]
        (tmp_, pmp, rmp) = pts[best[(False, True)][1]]
        (tmm, pmm, rmm) = pts[best[(False, False)][1]]
        rhi = ((qt - tmp_) * rpp + (tpp - qt) * rmp) / (tpp - tmp_)
        phi = ((qt - tmp_) * ppp + (tpp - qt) * pmp) / (tpp - tmp_)
        rlo = ((qt - tmm) * rpm + (tpm - qt) * rmm) / (tpm - tmm)
        plo = ((qt - tmm) * ppm + (tpm - qt) * pmm) / (tpm - tmm)
        return ((qp - plo) * rhi + (phi - qp) * rlo) / (phi - plo)
    num = den = 0.0
    for d2, i in best.values():
        w = 1.0 / (math.sqrt(d2) + 1e-5)
        num += w * pts[i][2]
        den += w
    return num / den


def mls(pts, qt, qp, skip=-1):
    nb = knn(pts, qt, qp, skip=skip)
    d2sup = nb[-1][0] * 1.1025
    A = [[0.0] * 3 for _ in range(3)]
    B = [0.0] * 3
    for d2, i in nb:
        dn = math.sqrt(d2 / d2sup)
        w = (1 - dn) ** 4 * (1 + 4 * dn)
        dt = pts[i][0] - qt
        dp = pts[i][1] - qp
        r = pts[i][2]
        x = (1.0, dt, dp)
        for a in range(3):
            for b in range(3):
                A[a][b] += w * x[a] * x[b]
            B[a] += w * x[a] * r
    ridge = 1e-6 * (A[1][1] + A[2][2]) + 1e-30
    A[1][1] += ridge
    A[2][2] += ridge
    det = (A[0][0] * (A[1][1] * A[2][2] - A[1][2] ** 2)
           - A[0][1] * (A[0][1] * A[2][2] - A[1][2] * A[0][2])
           + A[0][2] * (A[0][1] * A[1][2] - A[1][1] * A[0][2]))
    det1 = (B[0] * (A[1][1] * A[2][2] - A[1][2] ** 2)
            - A[0][1] * (B[1] * A[2][2] - A[1][2] * B[2])
            + A[0][2] * (B[1] * A[1][2] - A[1][1] * B[2]))
    deta = abs(A[0][0] * A[1][1] * A[2][2])
    if abs(det) > 1e-10 * deta and A[0][0] > 0:
        return det1 / det
    return B[0] / A[0][0]


def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else 'rosstab_cache.dat'
    pts, (zt, st, zp, sp) = load(fn)
    print(f'{len(pts)} cache points; logT span {st:.3f} dex, logP span {sp:.3f} dex')

    # pick a deep-CZ anchor: the stored point with the largest T that has
    # neighbors around it (80% up the normalized T range)
    anchor = min(pts, key=lambda x: abs(x[0] - 0.8))
    qp = anchor[1]

    # R1/R2: fine T-scan at fixed P
    ts = [anchor[0] + k * 0.0004 for k in range(-50, 51)]
    modes = {'bilinear': lambda t: bilinear(pts, t, qp),
             'shepard': lambda t: shepard(pts, t, qp),
             'mls': lambda t: mls(pts, t, qp)}
    print('\nR1 washboard metric (RMS second difference of log kappa along T-scan):')
    scans = {}
    for name, f in modes.items():
        v = [f(t) for t in ts]
        scans[name] = v
        d2 = [v[i-1] - 2*v[i] + v[i+1] for i in range(1, len(v)-1)]
        rms2 = math.sqrt(sum(x*x for x in d2)/len(d2))
        print(f'   {name:9s} {rms2:.3e}')

    print('\nR2 bubble-probe d(ln k)/d(ln T) scatter (dT/T=0.002 probes):')
    # convert normalized T step: dlogT = log10(1.002) / st
    dq = math.log10(1.002) / st
    for name, f in modes.items():
        g = [(f(t + dq) - f(t - dq)) / (2 * dq * st) for t in ts[10:-10:4]]
        mean = sum(g)/len(g)
        sd = math.sqrt(sum((x-mean)**2 for x in g)/len(g))
        print(f'   {name:9s} mean={mean:+.3f}  scatter={sd:.4f}')

    # R3: LOO residuals for MLS on a subsample
    sub = range(0, len(pts), max(1, len(pts)//300))
    res = [pts[i][2] - mls(pts, pts[i][0], pts[i][1], skip=i) for i in sub]
    mean = sum(res)/len(res)
    sd = math.sqrt(sum((x-mean)**2 for x in res)/len(res))
    print(f'\nR3 MLS leave-one-out residuals (dex): bias={mean:+.5f} scatter={sd:.5f}'
          f'  ({10**sd - 1:.1%} in kappa)')

    # R4: mode agreement over the scan
    diff = [scans['mls'][i] - scans['shepard'][i] for i in range(len(ts))]
    mean = sum(diff)/len(diff)
    print(f'\nR4 MLS-vs-Shepard mean offset over scan: {mean:+.5f} dex')


if __name__ == '__main__':
    main()
