#!/usr/bin/env python3
"""Verification battery for the partition-function spline interpolation.

Modes:
    python3 pf_verify.py probes    -> writes pf_probes.txt for the Fortran driver
    python3 pf_verify.py check     -> runs all property tests + Fortran cross-check

Tests (mode 'check'):
  T1  Node exactness: spline passes through every table value.
  T2  Boundary behavior: H2 clamps to end values outside [Tmin,Tmax];
      B&C returns U(T1) below the grid and the -1 sentinel at/above 1e4 K;
      continuity approaching each boundary from inside.
  T3  Overshoot/ringing: max |lnU_spline - lnU_linear| between knots for
      every species/stage; flag anything above LNU_TOL.
  T4  Positivity: spline U > 0 everywhere probed (guaranteed in ln-space,
      verified anyway).
  T5  Derivative smoothness: jump in d(lnU)/d(lnT) across interior knots,
      linear vs spline (the artifact being fixed, quantified).
  T6  EOS-stencil artifact: nabla_ad-style +-0.1% finite difference of
      lnU swept through a knot, linear vs spline kink amplitude.
  T7  Fortran cross-check: production code output equals the independent
      Python implementation at every probe to ~1e-12 relative.
"""
import math
import sys

DATA = '/Users/cconroy/kurucz/atlas12/data'
LNU_TOL = 0.05          # 5% ringing threshold (T3)
SPECIES = [(1, 1), (2, 1), (6, 1), (8, 1), (11, 1), (20, 2),
           (22, 1), (26, 1), (26, 2), (26, 3)]
NEGATIVES = [1, 6, 8]   # H-, C-, O-


# ----------------------------------------------------------------------
# Table parsing (mirrors the Fortran readers)
# ----------------------------------------------------------------------
def read_h2():
    T, Q = [], []
    for line in open(f'{DATA}/partfnh2.dat'):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        p = s.split()
        T.append(float(p[0])); Q.append(float(p[1]))
    return T, Q


ELEM = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S',
        'Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga',
        'Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
        'Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm',
        'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
        'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U']
ROMAN = {1: 'I', 2: 'II', 3: 'III'}


def read_bc():
    lines = [l for l in open(f'{DATA}/partfn_bc2016.dat')
             if l.strip() and not l.lstrip().startswith('#')]
    tgrid = [float(x) for x in lines[2].split()[2:]]   # after 'T [K]'
    tab, neg = {}, {}
    for l in lines[3:]:
        p = l.split()
        name, vals = p[0], [float(x) for x in p[1:]]
        if len(vals) != len(tgrid):
            continue
        if name.endswith('_-I') or name.endswith('-'):
            sym = name.split('_')[0].rstrip('-')
            if sym in ELEM:
                neg[ELEM.index(sym) + 1] = vals
            continue
        if '_' in name:
            sym, stage = name.split('_', 1)
            if sym in ELEM and stage in ('I', 'II', 'III'):
                tab[(ELEM.index(sym) + 1, ('I', 'II', 'III').index(stage) + 1)] = vals
    return tgrid, tab, neg


# ----------------------------------------------------------------------
# Independent implementations (must mirror the Fortran exactly)
# ----------------------------------------------------------------------
def spline_second_derivs(x, y):
    n = len(x)
    y2, c = [0.0] * n, [0.0] * n
    for i in range(1, n - 1):
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
        pden = sig * y2[i-1] + 2.0
        y2[i] = (sig - 1.0) / pden
        c[i] = ((6.0 * ((y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]))
                 / (x[i+1]-x[i-1])) - sig * c[i-1]) / pden
    for i in range(n - 2, 0, -1):
        y2[i] = y2[i] * y2[i+1] + c[i]
    return y2


def uniform_spline_second_derivs(y, h):
    # PARTFNH2's uniform-grid Thomas solve
    n = len(y)
    y2, c = [0.0] * n, [0.0] * n
    for i in range(1, n - 1):
        y2[i] = (6.0 * (y[i+1] - 2.0*y[i] + y[i-1]) / h**2 - y2[i-1]) / (4.0 - c[i-1])
        c[i] = 1.0 / (4.0 - c[i-1])
    for i in range(n - 2, 0, -1):
        y2[i] = y2[i] - c[i] * y2[i+1]
    return y2


def cubic_eval(f, ylo, yhi, y2lo, y2hi, h):
    return ((1-f)*ylo + f*yhi
            + ((1-f)**3 - (1-f)) * y2lo * h*h/6.0
            + (f**3 - f) * y2hi * h*h/6.0)


def h2_eval(T, tt, qq, q2, spline=True):
    tstart, tstep, npf = tt[0], tt[1]-tt[0], len(tt)
    T = min(tstart + (npf-1)*tstep, max(tstart, T))
    n = min(npf - 2, max(0, int((T - tstart)/tstep)))
    f = (T - (tstart + n*tstep)) / tstep
    if not spline:
        f = min(1.0, max(0.0, f))
        return qq[n] + (qq[n+1]-qq[n])*f
    return cubic_eval(f, qq[n], qq[n+1], q2[n], q2[n+1], tstep)


def bc_eval(T, tgrid, vals, y2, spline=True):
    if T <= tgrid[0]:
        return vals[0]
    if T >= tgrid[-1]:
        return -1.0
    lo = 0
    for i in range(len(tgrid)-1):
        if tgrid[i] <= T < tgrid[i+1]:
            lo = i
            break
    lt, ltlo, lthi = math.log(T), math.log(tgrid[lo]), math.log(tgrid[lo+1])
    f = (lt-ltlo)/(lthi-ltlo)
    lu_lo, lu_hi = math.log(vals[lo]), math.log(vals[lo+1])
    if not spline:
        return math.exp(lu_lo + f*(lu_hi-lu_lo))
    return math.exp(cubic_eval(f, lu_lo, lu_hi, y2[lo], y2[lo+1], lthi-ltlo))


# ----------------------------------------------------------------------
def build_probes():
    tt, qq = read_h2()
    tgrid, tab, neg = read_bc()
    probes = []
    # H2: nodes, node+-eps, midpoints, dense sweep, boundary probes
    for t in tt:
        probes += [('H2', 0, 0, t), ('H2', 0, 0, t*(1+1e-9)), ('H2', 0, 0, t-1e-6)]
    t = tt[0] + 0.37
    while t < tt[-1] + 500:
        probes.append(('H2', 0, 0, t)); t += 173.3
    probes += [('H2', 0, 0, 1.0), ('H2', 0, 0, tt[0]), ('H2', 0, 0, tt[-1]),
               ('H2', 0, 0, tt[-1] + 1e4)]
    # BC: per species -- all nodes, perturbed nodes, geometric midpoints,
    #     boundary probes
    for (z, ion) in SPECIES:
        for i, t in enumerate(tgrid):
            probes += [('BC', z, ion, t), ('BC', z, ion, t*(1+1e-9))]
            if i < len(tgrid)-1:
                probes.append(('BC', z, ion, math.sqrt(t*tgrid[i+1])))
        probes += [('BC', z, ion, tgrid[0]*0.5), ('BC', z, ion, tgrid[-1]*(1-1e-12)),
                   ('BC', z, ion, tgrid[-1]), ('BC', z, ion, 1.2e4)]
    for z in NEGATIVES:
        for i, t in enumerate(tgrid):
            probes.append(('NEG', z, 0, t))
            if i < len(tgrid)-1:
                probes.append(('NEG', z, 0, math.sqrt(t*tgrid[i+1])))
    with open('pf_probes.txt', 'w') as f:
        for tag, z, ion, t in probes:
            if tag == 'H2':
                f.write(f'H2 {t:.15e}\n')
            elif tag == 'BC':
                f.write(f'BC {z} {ion} {t:.15e}\n')
            else:
                f.write(f'NEG {z} {t:.15e}\n')
    print(f'{len(probes)} probes written')


def run_checks():
    tt, qq = read_h2()
    tgrid, tab, neg = read_bc()
    ltg = [math.log(t) for t in tgrid]
    q2 = uniform_spline_second_derivs(qq, tt[1]-tt[0])
    y2tab = {}
    for key, vals in {**tab, **{('neg', z): v for z, v in neg.items()}}.items():
        if all(v > 0 for v in vals):
            y2tab[key] = spline_second_derivs(ltg, [math.log(v) for v in vals])
        else:
            y2tab[key] = [0.0]*len(tgrid)
    fails = 0

    # T1 node exactness
    worst = 0.0
    for i, t in enumerate(tt):
        worst = max(worst, abs(h2_eval(t, tt, qq, q2) - qq[i]) / qq[i])
    for key, vals in tab.items():
        y2 = y2tab[key]
        for i, t in enumerate(tgrid[:-1]):   # top node returns sentinel by design
            u = bc_eval(t, tgrid, vals, y2)
            worst = max(worst, abs(u - vals[i]) / vals[i])
    print(f'T1 node exactness: worst rel err = {worst:.2e}',
          'PASS' if worst < 1e-12 else 'FAIL')
    fails += worst >= 1e-12

    # T2 boundaries
    ok = True
    ok &= abs(h2_eval(1.0, tt, qq, q2) - qq[0]) < 1e-12          # below H2 grid
    ok &= abs(h2_eval(tt[-1]+1e4, tt, qq, q2) - qq[-1]) < 1e-12  # above H2 grid
    ok &= abs(h2_eval(tt[-1]-1e-9, tt, qq, q2) - qq[-1]) < 1e-6  # continuity at top
    for key, vals in tab.items():
        y2 = y2tab[key]
        ok &= bc_eval(tgrid[0]*0.5, tgrid, vals, y2) == vals[0]      # below grid
        ok &= bc_eval(tgrid[-1], tgrid, vals, y2) == -1.0            # sentinel
        ok &= bc_eval(1.2e4, tgrid, vals, y2) == -1.0                # above grid
        u_in = bc_eval(tgrid[-1]*(1-1e-12), tgrid, vals, y2)         # just inside
        ok &= abs(u_in - vals[-1]) / vals[-1] < 1e-6                 # continuous to top node
    print('T2 boundary behavior:', 'PASS' if ok else 'FAIL')
    fails += not ok

    # T3 ringing: max |lnU_spline - lnU_linear| between knots
    worst_key, worst_dev = None, 0.0
    for key, vals in list(tab.items()) + [(('neg', z), v) for z, v in neg.items()]:
        y2 = y2tab[key]
        for i in range(len(tgrid)-1):
            for ffrac in (0.1, 0.25, 0.5, 0.75, 0.9):
                t = math.exp(ltg[i] + ffrac*(ltg[i+1]-ltg[i]))
                us = bc_eval(t, tgrid, vals, y2, spline=True)
                ul = bc_eval(t, tgrid, vals, y2, spline=False)
                if us > 0 and ul > 0:
                    d = abs(math.log(us) - math.log(ul))
                    if d > worst_dev:
                        worst_dev, worst_key = d, (key, tgrid[i])
    print(f'T3 max |dlnU| spline-vs-linear = {worst_dev:.4f} at {worst_key}',
          'PASS' if worst_dev < LNU_TOL else 'FAIL')
    fails += worst_dev >= LNU_TOL

    # T4 positivity at dense probes
    ok = True
    for key, vals in tab.items():
        y2 = y2tab[key]
        for i in range(len(tgrid)-1):
            for ffrac in (0.25, 0.5, 0.75):
                t = math.exp(ltg[i] + ffrac*(ltg[i+1]-ltg[i]))
                ok &= bc_eval(t, tgrid, vals, y2) > 0
    print('T4 positivity:', 'PASS' if ok else 'FAIL')
    fails += not ok

    # T5 derivative jump across knots (report only, improvement metric)
    def dln(evalf, t, dt):
        return (math.log(evalf(t+dt)) - math.log(evalf(t-dt))) / (2*dt/t)
    key = (26, 1)
    vals, y2 = tab[key], y2tab[key]
    jl, js = 0.0, 0.0
    for i in range(5, len(tgrid)-2):
        t = tgrid[i]
        eps = t*1e-6
        for sp, acc in ((False, 'lin'), (True, 'spl')):
            lo = (math.log(bc_eval(t-2*eps, tgrid, vals, y2, sp)) -
                  math.log(bc_eval(t-4*eps, tgrid, vals, y2, sp))) / (2*eps/t)
            hi = (math.log(bc_eval(t+4*eps, tgrid, vals, y2, sp)) -
                  math.log(bc_eval(t+2*eps, tgrid, vals, y2, sp))) / (2*eps/t)
            if sp:
                js = max(js, abs(hi-lo))
            else:
                jl = max(jl, abs(hi-lo))
    print(f'T5 Fe I dlnU/dlnT jump across knots: linear {jl:.4f} -> spline {js:.6f}')

    # T6 EOS-stencil kink sweep through the Fe I knot nearest 5040K
    knot = min(tgrid, key=lambda t: abs(t-5040))
    amp = {True: [], False: []}
    for sp in (False, True):
        for k in range(-40, 41):
            t = knot + k*2.0
            up = bc_eval(t*1.001, tgrid, vals, y2, sp)
            dn = bc_eval(t*0.999, tgrid, vals, y2, sp)
            amp[sp].append((math.log(up)-math.log(dn))/0.002)
        a = amp[sp]
        spread = max(a) - min(a)
        print(f'T6 dlnU/dlnT spread over +-80K sweep at {knot:.0f}K knot '
              f'({"spline" if sp else "linear"}): {spread:.4f}')

    # T7 Fortran cross-check
    try:
        rows = [l.split(',') for l in open('pf_fortran.csv')]
    except FileNotFoundError:
        print('T7 SKIPPED (run the Fortran driver first)')
        sys.exit(1 if fails else 0)
    worst, worst_row = 0.0, None
    for tag, z, ion, t, u in rows:
        z, ion, t, u = int(z), int(ion), float(t), float(u)
        if tag == 'H2':
            ref = h2_eval(t, tt, qq, q2)
        elif tag == 'BC':
            key = (z, ion)
            if key not in tab:
                continue
            ref = bc_eval(t, tgrid, tab[key], y2tab[key])
        else:
            if z not in neg:
                continue
            ref = bc_eval(t, tgrid, neg[z], y2tab[('neg', z)])
        denom = max(abs(ref), 1e-30)
        d = abs(u - ref)/denom
        if d > worst:
            worst, worst_row = d, (tag, z, ion, t, u, ref)
    print(f'T7 Fortran-vs-Python worst rel diff = {worst:.2e}',
          'PASS' if worst < 1e-10 else f'FAIL at {worst_row}')
    fails += worst >= 1e-10

    print('=' * 60)
    print('ALL PASS' if fails == 0 else f'{fails} FAILURES')
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'probes':
        build_probes()
    else:
        run_checks()
