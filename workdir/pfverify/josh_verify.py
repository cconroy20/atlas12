#!/usr/bin/env python3
"""Verification battery for the JOSH-path numerics (MAP1/DERIV/INTEG).

Reads josh_fortran.csv written by joshtest (which drives the ACTUAL
production routines over perturbation sweeps that are analytic in eps)
and compares against candidate replacement methods implemented here:

  linear   piecewise linear (C0; linear operator in the data)
  steffen  Steffen 1990 monotone C1 cubic (production candidate)
  spline   natural cubic spline (C2; linear operator in the data)
  parab    3-point Lagrange node derivative (Test B only)

Tests
  A   surface flux H = CH_J . remap(S(eps)) with grids FIXED and only
      the source data moving.  Linear/spline are linear-in-data, so
      their H(eps) inherits the analytic smoothness of S(eps) exactly:
      their washboard is the smooth floor.  Excess washboard in the
      MAP1 row is pure blend-weight noise.
  A2  opacity perturbed instead: the tau grid slides under XTAU8
      (INTEG + MAP1 jointly vs candidate quadrature + interpolation).
  B   nested node derivatives (deep-branch Eddington HNU, J-S).
  C   INTEG quadrature accuracy vs analytic integral of exp(x).
  D   JOSH scattering Gauss-Seidel Lambda-iteration (exact replica of
      the production loop, tol 1e-5, max 50 sweeps) vs direct solve of
      (I - diag(alpha) COEFJ) S = (1-alpha) Sbar with COEFJ from
      data/blockj.dat.

Washboard metric: RMS of the second difference of the recorded output
over the eps sweep, normalized by the median |output| (same convention
as rosstab_verify.py R1).

Usage: python3 josh_verify.py [josh_fortran.csv]
"""
import math
import sys
from bisect import bisect_right

HPL = 6.62607015e-27
KB = 1.380649e-16
CL = 2.99792458e10
FREQ = CL / 500.0e-7
NEPS = 401
EPSMAX = 2.0e-3
BLOCKJ = '../../data/blockj.dat'


# ----------------------------------------------------------------------
# interpolation / quadrature primitives
# ----------------------------------------------------------------------
def sgn(v):
    return (v > 0) - (v < 0)


def slopes_steffen(x, y):
    n = len(x)
    h = [x[i + 1] - x[i] for i in range(n - 1)]
    s = [(y[i + 1] - y[i]) / h[i] for i in range(n - 1)]
    d = [0.0] * n
    for i in range(1, n - 1):
        if s[i - 1] * s[i] <= 0.0:
            d[i] = 0.0
        else:
            p = (s[i - 1] * h[i] + s[i] * h[i - 1]) / (h[i - 1] + h[i])
            d[i] = (sgn(s[i - 1]) + sgn(s[i])) * min(
                abs(s[i - 1]), abs(s[i]), 0.5 * abs(p))
    p0 = s[0] * (1.0 + h[0] / (h[0] + h[1])) - s[1] * h[0] / (h[0] + h[1])
    if p0 * s[0] <= 0.0:
        d[0] = 0.0
    elif abs(p0) > 2.0 * abs(s[0]):
        d[0] = 2.0 * s[0]
    else:
        d[0] = p0
    pn = s[-1] * (1.0 + h[-1] / (h[-1] + h[-2])) - s[-2] * h[-1] / (h[-1] + h[-2])
    if pn * s[-1] <= 0.0:
        d[-1] = 0.0
    elif abs(pn) > 2.0 * abs(s[-1]):
        d[-1] = 2.0 * s[-1]
    else:
        d[-1] = pn
    return d


def eval_steffen(x, y, d, xq):
    n = len(x)
    if xq <= x[0]:
        return y[0] + d[0] * (xq - x[0])
    if xq >= x[-1]:
        return y[-1] + d[-1] * (xq - x[-1])
    i = bisect_right(x, xq) - 1
    i = min(i, n - 2)
    h = x[i + 1] - x[i]
    s = (y[i + 1] - y[i]) / h
    a = (d[i] + d[i + 1] - 2.0 * s) / (h * h)
    b = (3.0 * s - 2.0 * d[i] - d[i + 1]) / h
    t = xq - x[i]
    return y[i] + t * (d[i] + t * (b + t * a))


def cumint_steffen(x, y, start):
    d = slopes_steffen(x, y)
    out = [start]
    for i in range(len(x) - 1):
        h = x[i + 1] - x[i]
        s = (y[i + 1] - y[i]) / h
        a = (d[i] + d[i + 1] - 2.0 * s) / (h * h)
        b = (3.0 * s - 2.0 * d[i] - d[i + 1]) / h
        seg = h * (y[i] + h * (d[i] / 2.0 + h * (b / 3.0 + h * a / 4.0)))
        out.append(out[-1] + seg)
    return out


def spline_m(x, y):
    """Natural cubic spline second derivatives (Thomas algorithm)."""
    n = len(x)
    m = [0.0] * n
    if n < 3:
        return m
    cp = [0.0] * n
    dp = [0.0] * n
    for i in range(1, n - 1):
        hl = x[i] - x[i - 1]
        hr = x[i + 1] - x[i]
        rhs = 6.0 * ((y[i + 1] - y[i]) / hr - (y[i] - y[i - 1]) / hl)
        a = hl
        b = 2.0 * (hl + hr)
        c = hr
        if i == 1:
            cp[i] = c / b
            dp[i] = rhs / b
        else:
            den = b - a * cp[i - 1]
            cp[i] = c / den
            dp[i] = (rhs - a * dp[i - 1]) / den
    for i in range(n - 2, 0, -1):
        m[i] = dp[i] - cp[i] * m[i + 1]
    return m


def eval_spline(x, y, m, xq):
    n = len(x)
    if xq <= x[0]:
        h = x[1] - x[0]
        d0 = (y[1] - y[0]) / h - h / 6.0 * (2.0 * m[0] + m[1])
        return y[0] + d0 * (xq - x[0])
    if xq >= x[-1]:
        h = x[-1] - x[-2]
        dn = (y[-1] - y[-2]) / h + h / 6.0 * (2.0 * m[-1] + m[-2])
        return y[-1] + dn * (xq - x[-1])
    i = min(bisect_right(x, xq) - 1, n - 2)
    h = x[i + 1] - x[i]
    A = (x[i + 1] - xq) / h
    B = (xq - x[i]) / h
    return (A * y[i] + B * y[i + 1]
            + ((A ** 3 - A) * m[i] + (B ** 3 - B) * m[i + 1]) * h * h / 6.0)


def slopes_spline(x, y):
    m = spline_m(x, y)
    n = len(x)
    d = [0.0] * n
    for i in range(n - 1):
        h = x[i + 1] - x[i]
        d[i] = (y[i + 1] - y[i]) / h - h / 6.0 * (2.0 * m[i] + m[i + 1])
    h = x[-1] - x[-2]
    d[-1] = (y[-1] - y[-2]) / h + h / 6.0 * (2.0 * m[-1] + m[-2])
    return d


def cumint_spline(x, y, start):
    m = spline_m(x, y)
    out = [start]
    for i in range(len(x) - 1):
        h = x[i + 1] - x[i]
        seg = 0.5 * h * (y[i] + y[i + 1]) - h ** 3 / 24.0 * (m[i] + m[i + 1])
        out.append(out[-1] + seg)
    return out


def eval_linear(x, y, xq):
    n = len(x)
    if xq <= x[0]:
        i = 0
    elif xq >= x[-1]:
        i = n - 2
    else:
        i = min(bisect_right(x, xq) - 1, n - 2)
    w = (xq - x[i]) / (x[i + 1] - x[i])
    return y[i] + w * (y[i + 1] - y[i])


def cumint_trapz(x, y, start):
    out = [start]
    for i in range(len(x) - 1):
        out.append(out[-1] + 0.5 * (x[i + 1] - x[i]) * (y[i] + y[i + 1]))
    return out


def slopes_parab(x, y):
    """Unlimited 3-point Lagrange derivative at nodes (linear in data)."""
    n = len(x)
    d = [0.0] * n
    d[0] = (y[1] - y[0]) / (x[1] - x[0])
    d[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    for i in range(1, n - 1):
        hl = x[i] - x[i - 1]
        hr = x[i + 1] - x[i]
        d[i] = (y[i + 1] * hl * hl - y[i - 1] * hr * hr
                + y[i] * (hr * hr - hl * hl)) / (hl * hr * (hl + hr))
    return d


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def planck(t):
    return 2.0 * HPL * FREQ ** 3 / CL ** 2 / (math.exp(HPL * FREQ / (KB * t)) - 1.0)


def washboard(seq):
    n = len(seq)
    med = sorted(abs(v) for v in seq)[n // 2]
    if med == 0.0:
        med = 1.0
    ss = 0.0
    for i in range(1, n - 1):
        ss += (seq[i + 1] - 2.0 * seq[i] + seq[i - 1]) ** 2
    return math.sqrt(ss / (n - 2)) / med


def solve_dense(a, b):
    """Gaussian elimination with partial pivoting; a, b copied."""
    n = len(b)
    a = [row[:] for row in a]
    b = b[:]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(a[r][i]))
        if p != i:
            a[i], a[p] = a[p], a[i]
            b[i], b[p] = b[p], b[i]
        for r in range(i + 1, n):
            f = a[r][i] / a[i][i]
            for c in range(i, n):
                a[r][c] -= f * a[i][c]
            b[r] -= f * b[i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = b[i] - sum(a[i][c] * x[c] for c in range(i + 1, n))
        x[i] = s / a[i][i]
    return x


# ----------------------------------------------------------------------
# load Fortran results
# ----------------------------------------------------------------------
fn = sys.argv[1] if len(sys.argv) > 1 else 'josh_fortran.csv'
probes = None
xtau, chj = [], []
rhox, tprof, abross, tauros = [], [], [], []
rows_a, rows_a2, rows_b, rows_c = [], [], [], []
rows_c2 = []
aexp_c = 1.0
for line in open(fn):
    p = [q.strip() for q in line.split(',')]
    if p[0] == 'PROBES':
        probes = [int(v) - 1 for v in p[1:4]]
    elif p[0] == 'G':
        xtau.append(float(p[2]))
        chj.append(float(p[3]))
    elif p[0] == 'P':
        rhox.append(float(p[2]))
        tprof.append(float(p[3]))
        abross.append(float(p[4]))
        tauros.append(float(p[5]))
    elif p[0] == 'A':
        rows_a.append([float(v) for v in p[1:]])
    elif p[0] == 'A2':
        rows_a2.append([float(v) for v in p[1:]])
    elif p[0] == 'B':
        rows_b.append([float(v) for v in p[1:]])
    elif p[0] == 'CA':
        aexp_c = float(p[1])
    elif p[0] == 'C':
        rows_c.append([float(v) for v in p[2:]])
    elif p[0] == 'C2':
        rows_c2.append([float(v) for v in p[2:]])

n = len(rhox)
g = [math.exp(-0.5 * ((math.log10(tr) - 0.0) / 0.35) ** 2) for tr in tauros]
eps_grid = [EPSMAX * i / (NEPS - 1) for i in range(NEPS)]
assert len(rows_a) == NEPS and len(rows_a2) == NEPS and len(rows_b) == NEPS

print('JOSH-path numerics verification  (n_layers=%d, probes tau=%s)'
      % (n, ['%.3g' % tauros[j] for j in probes]))
print('=' * 70)

# ----------------------------------------------------------------------
# Test A: fixed grids, S(eps) moves
# ----------------------------------------------------------------------
h_lin, h_stef, h_spl, s05_lin, s05_stef, s05_spl = [], [], [], [], [], []
for eps in eps_grid:
    snu = [planck(tprof[j] * (1.0 + eps * g[j])) for j in range(n)]
    vals = [eval_linear(tauros, snu, xq) for xq in xtau]
    h_lin.append(sum(c * v for c, v in zip(chj, vals)))
    s05_lin.append(vals[25])
    d = slopes_steffen(tauros, snu)
    vals = [eval_steffen(tauros, snu, d, xq) for xq in xtau]
    h_stef.append(sum(c * v for c, v in zip(chj, vals)))
    s05_stef.append(vals[25])
    m = spline_m(tauros, snu)
    vals = [eval_spline(tauros, snu, m, xq) for xq in xtau]
    h_spl.append(sum(c * v for c, v in zip(chj, vals)))
    s05_spl.append(vals[25])

h_map1 = [r[1] for r in rows_a]
s05_map1 = [r[2] for r in rows_a]
print('\nA. Surface flux H(eps), grids fixed (washboard = RMS 2nd diff / med|H|)')
print('   MAP1 (production) : %.3e' % washboard(h_map1))
print('   steffen           : %.3e' % washboard(h_stef))
print('   spline            : %.3e   <- smooth floor (linear operator)' % washboard(h_spl))
print('   linear            : %.3e   <- smooth floor (linear operator)' % washboard(h_lin))
print('   raw interpolant at tau=0.5:')
print('   MAP1 %.3e   steffen %.3e   spline %.3e   linear %.3e'
      % (washboard(s05_map1), washboard(s05_stef),
         washboard(s05_spl), washboard(s05_lin)))
mo = lambda a, b: sum(abs(x - y) / y for x, y in zip(a, b)) / len(a)
print('   accuracy cross-check, mean |rel diff| of H vs MAP1: '
      'steffen %.2e  spline %.2e  linear %.2e'
      % (mo(h_stef, h_map1), mo(h_spl, h_map1), mo(h_lin, h_map1)))

# ----------------------------------------------------------------------
# Test A2: S fixed, tau grid slides (INTEG + MAP1 jointly)
# ----------------------------------------------------------------------
snu0 = [planck(t) for t in tprof]
h2_lin, h2_stef, h2_spl = [], [], []
for eps in eps_grid:
    abp = [abross[j] * (1.0 + eps * g[j]) for j in range(n)]
    start = abp[0] * rhox[0]
    tau = cumint_trapz(rhox, abp, start)
    h2_lin.append(sum(c * eval_linear(tau, snu0, xq) for c, xq in zip(chj, xtau)))
    tau = cumint_steffen(rhox, abp, start)
    d = slopes_steffen(tau, snu0)
    h2_stef.append(sum(c * eval_steffen(tau, snu0, d, xq) for c, xq in zip(chj, xtau)))
    tau = cumint_spline(rhox, abp, start)
    m = spline_m(tau, snu0)
    h2_spl.append(sum(c * eval_spline(tau, snu0, m, xq) for c, xq in zip(chj, xtau)))

h2_f = [r[1] for r in rows_a2]
print('\nA2. Surface flux H(eps), opacity perturbed so tau grid slides')
print('   INTEG+MAP1 (production) : %.3e' % washboard(h2_f))
print('   steffen(int+interp)     : %.3e' % washboard(h2_stef))
print('   spline(int+interp)      : %.3e' % washboard(h2_spl))
print('   trapz+linear            : %.3e' % washboard(h2_lin))
print('   tau(eps) at probe layers, production INTEG washboard: '
      + '  '.join('%.3e' % washboard([r[2 + k] for r in rows_a2]) for k in range(3)))

# ----------------------------------------------------------------------
# Test B: nested node derivatives (deep-branch Eddington)
# ----------------------------------------------------------------------
print('\nB. Nested derivatives at probe layers: HNU=dS/dtau/3, JMS=d(HNU)/dtau')
meth_b = {'parab': slopes_parab, 'steffen': slopes_steffen, 'spline': slopes_spline}
res_b = {k: {'h': [[] for _ in range(3)], 'j': [[] for _ in range(3)]} for k in meth_b}
for eps in eps_grid:
    snu = [planck(tprof[j] * (1.0 + eps * g[j])) for j in range(n)]
    for k, fs in meth_b.items():
        d1 = [v / 3.0 for v in fs(tauros, snu)]
        d2 = fs(tauros, d1)
        for kk in range(3):
            res_b[k]['h'][kk].append(d1[probes[kk]])
            res_b[k]['j'][kk].append(d2[probes[kk]])
hdr = '   %-22s' + '  %11s' * 3
print(hdr % ('washboard', 'tau=0.5', 'tau=2', 'tau=10'))
for tag, cols in (('HNU', (1, 2, 3)), ('JMS', (4, 5, 6))):
    vals = ['%.3e' % washboard([r[c] for r in rows_b]) for c in cols]
    print(('   %-22s' + '  %11s' * 3) % ('DERIV (production) ' + tag, *vals))
    for k in ('parab', 'steffen', 'spline'):
        key = 'h' if tag == 'HNU' else 'j'
        vals = ['%.3e' % washboard(res_b[k][key][kk]) for kk in range(3)]
        print(('   %-22s' + '  %11s' * 3) % ('%s %s' % (k, tag), *vals))

# ----------------------------------------------------------------------
# Test C: quadrature accuracy vs analytic integral of exp(x)
# ----------------------------------------------------------------------
xc = rhox
fexp = [math.exp(aexp_c * v) for v in xc]
exact = [r[1] for r in rows_c]
integ = [r[2] for r in rows_c]
scale = max(exact)


def maxrel(cum):
    return max(abs(c - e) for c, e in zip(cum, exact)) / scale


print('\nC. Cumulative integral of exp(x) on the real RHOX grid, max |err|/max')
print('   INTEG/PARCOE (production) : %.3e' % maxrel(integ))
py_stef = cumint_steffen(xc, fexp, 0.0)
print('   steffen integral          : %.3e' % maxrel(py_stef))
print('   spline integral           : %.3e' % maxrel(cumint_spline(xc, fexp, 0.0)))
print('   trapezoid                 : %.3e' % maxrel(cumint_trapz(xc, fexp, 0.0)))
if rows_c2:
    f_stef = [r[2] for r in rows_c2]
    print('   INTEG quad=1 (production) : %.3e   [Fortran vs Python steffen: %.2e]'
          % (maxrel(f_stef),
             max(abs(a - b) for a, b in zip(f_stef, py_stef)) / scale))

# ----------------------------------------------------------------------
# Test D: scattering Gauss-Seidel Lambda iteration vs direct solve
# ----------------------------------------------------------------------
print('\nD. JOSH scattering iteration vs direct solve (COEFJ from blockj.dat)')
toks = []
with open(BLOCKJ) as f:
    lines = f.readlines()
for line in lines[3:]:
    toks.extend(float(v) for v in line.split())
assert len(toks) >= 2601, len(toks)
coefj = [[toks[i + 51 * j] for j in range(51)] for i in range(51)]  # column-major

tq = [eval_linear(tauros, tprof, xq) for xq in xtau]
sbar = [planck(t) for t in tq]


def gs_replica(alpha):
    xs = sbar[:]
    diag = [1.0 - alpha[l] * coefj[l][l] for l in range(51)]
    diag = [d if abs(d) >= 1e-30 else math.copysign(1e-30, d) for d in diag]
    xsb = [(1.0 - alpha[l]) * sbar[l] for l in range(51)]
    sweeps = 50
    conv = False
    for it in range(50):
        iferr = False
        for k in range(50, -1, -1):
            delxs = sum(coefj[k][mm] * xs[mm] for mm in range(51))
            delxs = (delxs * alpha[k] + xsb[k] - xs[k]) / diag[k]
            if abs(delxs / xs[k]) > 1e-5:
                iferr = True
            xs[k] = max(xs[k] + delxs, 1e-37)
        if not iferr:
            sweeps = it + 1
            conv = True
            break
    return xs, sweeps, conv


cases = [('const 0.30', [0.30] * 51), ('const 0.90', [0.90] * 51),
         ('const 0.99', [0.99] * 51), ('const 0.999', [0.999] * 51),
         ('const 0.9999', [0.9999] * 51),
         ('rayleigh-like', [0.999 * math.exp(-x / 5.0) for x in xtau])]
print('   %-14s %8s %6s %14s' % ('alpha', 'sweeps', 'conv', 'max|GS-dir|/dir'))
for name, alpha in cases:
    xs, sweeps, conv = gs_replica(alpha)
    a = [[(1.0 if i == j else 0.0) - alpha[i] * coefj[i][j] for j in range(51)]
         for i in range(51)]
    b = [(1.0 - alpha[i]) * sbar[i] for i in range(51)]
    xd = solve_dense(a, b)
    dev = max(abs(u - v) / abs(v) for u, v in zip(xs, xd))
    print('   %-14s %8d %6s %14.3e' % (name, sweeps, 'yes' if conv else 'NO', dev))
