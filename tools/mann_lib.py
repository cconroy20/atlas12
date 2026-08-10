#!/usr/bin/env python3
"""Shared library for the Mann point-comparison infrastructure.

Single source of truth for everything the Mann-validation tools used to
duplicate: the M_params table read, GJ<->PM name resolution, spectrum
readers for both libraries (Mann+15 2016 ascii + Mann+13 Calibrators2013
FITS), the measured Mann LSF, flux conventions, band windows, atlas12/
synthe launch helpers, and PHOENIX NewEra LowRes grid access with
trilinear interpolation.

Libraries:
  mann15    ~/sps/SPECTRA/Mann/<NAME>.ascii   (183 M dwarfs, private 2016
            distribution; params from M_params.fits; SNIFS/SpeX splice 9500 A)
  calib2013 ~/sps/SPECTRA/Mann/Calibrators2013/<GJ>.fits  (public Mann+13
            calibrators; params curated in KDWARF_PARAMS, model-free
            theta+Fbol; splice ~8100 A; err==0 points are BT-Settl fill)
"""

import io
import json
import os
import subprocess
import sys

import numpy as np

# ---------------------------------------------------------------- paths
REPO      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANN_DIR  = os.path.expanduser("~/sps/SPECTRA/Mann")
PARAMS    = os.path.join(MANN_DIR, "M_params.fits")
CAL_DIR   = os.path.join(MANN_DIR, "Calibrators2013")
INT_STARS = os.path.join(MANN_DIR, "Interferometric_Stars.txt")
RUN_ROOT  = os.path.join(REPO, "workdir", "mann")

PHOENIX_DIR   = os.path.expanduser(
    "~/kurucz/grids/PHOENIX-NewEraV2-LowRes-SPECTRA")
PHOENIX_CACHE = os.path.join(PHOENIX_DIR, "cache")

C3K_DIR = os.path.expanduser(
    "~/kurucz/grids/c3k_v2.3/at12_feh+0.00_afe+0.0/spec")
C3K_R   = 300000.0

# Starting models for atlas12 (SCALE_MODEL regrids to the target Teff/logg).
# czc ladder below 4000 K; C3K grid atmospheres at 4000-5500 K (validated
# K-dwarf recipe 2026-08-08: t04500g4.50 converged all three K dwarfs).
_C3K_ATM = os.path.expanduser("~/kurucz/grids/test/atm")
START_MODELS = [
    (2500.0, os.path.join(REPO, "workdir", "czc17_2500.atm")),
    (2800.0, os.path.join(REPO, "workdir", "czc19_2800c.atm")),
    (3500.0, os.path.join(REPO, "workdir", "czc19_3500.atm")),
    (4000.0, os.path.join(_C3K_ATM, "t04000g4.50.atm")),
    (4500.0, os.path.join(_C3K_ATM, "t04500g4.50.atm")),
    (5500.0, os.path.join(_C3K_ATM, "t05500g4.50.atm")),
    (5777.0, os.path.join(REPO, "workdir", "czc17_solar.atm")),
]

# ---------------------------------------------------- physical constants (cgs)
GRAV     = 6.67430e-8
MSUN     = 1.98892e33
RSUN     = 6.957e10
PC       = 3.0856776e18
C_ANG    = 2.99792458e18          # speed of light [Angstrom/s]
SIGMA_SB = 5.670374419e-5

# --------------------------------------------------------- measured Mann LSF
# SNIFS effective LSF: constant-FWHM Gaussian, measured 2026-08-09 by
# chunk-fitting R=300k syntheses to GJ887 / GJ380 / GJ105A: dlam ~ 11 A
# FWHM (R rising ~400 blue -> ~850 at the splice).  SpeX side R ~ 2000.
# The SNIFS->SpeX splice wavelength differs by library generation.
SNIFS_DLAM_FWHM = 11.0            # Angstrom
OBS_R_RED       = 2000.0
SPLIT_MANN15    = 9500.0          # Angstrom, 2016 private library
SPLIT_CALIB2013 = 8100.0          # Angstrom, public Calibrators2013

# ------------------------------------------------------------- band windows
# Canonical copies (tmin_fit.py / tmin_rf.py historically carried duplicates).
TIO_BANDS = [(5500, 5700), (5900, 6100), (7053, 7270), (7450, 7650)]
CAH       = {"cah2": (6814., 6846.), "cah3": (6960., 6990.),
             "tio5": (7126., 7135.)}
CAH_DEN   = (7042., 7046.)
CA_CORES  = {"ca8542": (8538., 8551.), "ca8662": (8658., 8671.)}
CA_SIDE   = (8600., 8640.)

# --------------------------------------------------- curated K-dwarf params
# Model-free: Teff from theta_LD + Fbol via Stefan-Boltzmann, dilution
# (theta/2)^2; logg from the literature.  Fbol in erg/s/cm^2 at Earth.
KDWARF_PARAMS = {
    "GJ820A": dict(
        theta_mas=1.775, fbol=37.75e-8, logg=4.63, feh=-0.19, spt="K5V",
        comment="61 Cyg A, original Gaia benchmark. theta Kervella+08 VLTI "
                "(1.775+-0.013); Fbol Boyajian+12. GBSv3 Fbol 39.24e-8 "
                "(-> Teff 4398) is the alternative; B12 value validated "
                "2026-08-08 (integral model/data 0.997)."),
    "GJ892": dict(
        theta_mas=1.106, fbol=21.6e-8, logg=4.55, feh=+0.07, spt="K3V",
        comment="HD 219134. theta CHARA 1.106+-0.007 (34 obs); Fbol GBSv3 "
                "(Soubiran+24) -> Teff ~4800. Boyajian+12 Fbol 19.85e-8 "
                "(-> 4699) produced a uniform 6-7% flux deficit 2026-08-08; "
                "the spectrum integral votes for GBSv3."),
    "GJ105A": dict(
        theta_mas=1.030, fbol=16.70e-8, logg=4.52, feh=-0.08, spt="K3V",
        comment="HD 16160. theta CHARA 1.030 -> Teff 4662 (Fbol back-derived "
                "from that Teff). NOT on the ladder: blue calibration tilt "
                "(0.84-0.86 blue vs ~1.0 red) that no Fbol choice fixes; "
                "nearby M4 companion. Kept resolvable for the LSF work."),
    # Excluded entirely (documented 2026-08-08/09):
    #   GJ570A -- theta disputed: Demory 1.177 vs Rains 1.130 (~150 K).
    #   GJ702B -- 70 Oph B: 2 visibility obs, blended Fbol.
}

# ------------------------------------------------------------- star ladder
# The standing 10-star point-comparison set (user-selected 2026-08-09):
# best-quality spectra spanning the full Teff range, GJ names, Mann+15
# M dwarfs on catalog table params + Mann+13 K dwarfs on interferometric.
STAR_LADDER = [
    ("GJ644C", "vB 8, M7, 2700 K; coolest rung; CaOH super-lines make "
               "synthesis slow"),
    ("GJ905",  "M5.2, 2931 K; NOTE interferometric Teff 2811 disagrees "
               "with table by 120 K (unvetted theta)"),
    ("GJ699",  "Barnard's star, M4, 3224 K"),
    ("GJ411",  "M2, 3548 K; interferometric campaign star"),
    ("GJ15A",  "M3, 3603 K; do not confuse PM_I00184+4401 = GJ15B"),
    ("GJ887",  "M1, 3689 K; anchor star, regression reference"),
    ("GJ205",  "M0, 3825 K; metal-rich [Fe/H]=+0.49"),
    ("GJ380",  "K7.6, 4131 K; 4100 K limit-cycle star, expect ~3% RMS"),
    ("GJ820A", "61 Cyg A, K5, 4361 K; Calibrators2013, model-free params"),
    ("GJ892",  "HD 219134, K3, ~4800 K; Calibrators2013, model-free params"),
]


# ================================================================ table load
_CACHE = {}


def load_mann():
    """Return the Mann+15 table as a dict of 1-D arrays (183 stars)."""
    if "table" in _CACHE:
        return _CACHE["table"]
    from astropy.io import fits
    d = fits.open(PARAMS)[1].data[0]
    name = np.array([n.strip() for n in d["NAME"]])
    cns3 = np.array([n.strip() for n in d["CNS3"]])
    mass, radius = d["MASS"], d["RADIUS"]
    logg = np.log10(GRAV * mass * MSUN / (radius * RSUN) ** 2)
    t = dict(name=name, cns3=cns3,
             spt=np.array([s.strip() for s in d["SPT"]]),
             teff=d["TEFF"], feh=d["FEH"], mh=d["MH"], mass=mass,
             radius=radius, dist=d["DISTANCE"], fbol=d["FBOL"],
             lum=d["LUMINOSITY"], logg=logg)
    _CACHE["table"] = t
    return t


def _gj_key(s):
    """Normalize a Gliese-style designation to compact 'GJ...' form, or None.

    'Gl 644 C' -> 'GJ644C', 'GJ 1111' -> 'GJ1111', 'gj887' -> 'GJ887'.
    """
    c = s.strip().upper().replace(" ", "")
    if c.startswith("GLIESE"):
        c = "GJ" + c[6:]
    elif c.startswith("GL") and len(c) > 2 and not c[2].isalpha():
        c = "GJ" + c[2:]
    return c if c.startswith("GJ") and len(c) > 2 else None


def gj_maps():
    """Return (gj->table index, table index->gj display name).

    Built from the M_params CNS3 column (161/183 populated), augmented by
    Interferometric_Stars.txt (compact GJ names for the LBOI stars).
    """
    if "gjmaps" in _CACHE:
        return _CACHE["gjmaps"]
    t = load_mann()
    gj2i, i2gj = {}, {}
    for i, c in enumerate(t["cns3"]):
        k = _gj_key(c) if c else None
        if k:
            gj2i.setdefault(k, i)
            i2gj.setdefault(i, k)
    if os.path.isfile(INT_STARS):
        name2i = {n: i for i, n in enumerate(t["name"])}
        with open(INT_STARS) as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 2 or parts[0].startswith("PM_Name"):
                    continue
                pm, common = parts[0], parts[1]
                k = _gj_key(common)
                if k and pm in name2i:
                    gj2i.setdefault(k, name2i[pm])
                    i2gj.setdefault(name2i[pm], k)
    _CACHE["gjmaps"] = (gj2i, i2gj)
    return gj2i, i2gj


class Star:
    """Resolved star record: identity + parameters + library dispatch."""

    def __repr__(self):
        return (f"Star({self.name}, lib={self.library}, Teff={self.teff:.0f}, "
                f"logg={self.logg:.2f}, [Fe/H]={self.feh:+.2f})")


def resolve(query):
    """Resolve a star by GJ name ('GJ887', 'Gl 15A') or Mann PM_I name.

    Returns a Star with: name (canonical run-dir name, GJ-form when it
    exists), gj_name, pm_name, library ('mann15'|'calib2013'), spt, teff,
    logg, feh, mh, radius, dist, fbol, theta_mas, dilute, obs_split.
    The 2016 mann15 library wins when a star exists in both (later, better
    calibration); calib2013 is reached for the K-dwarf-only stars.
    """
    t = load_mann()
    gj2i, i2gj = gj_maps()
    q = query.strip()
    idx = None
    if q in t["name"]:
        idx = int(np.where(t["name"] == q)[0][0])
    else:
        k = _gj_key(q)
        if k and k in gj2i:
            idx = gj2i[k]
        elif k and k in KDWARF_PARAMS:
            return _resolve_kdwarf(k)
    if idx is None:
        raise KeyError(
            f"star {query!r} not found: not an M_params NAME, not a known "
            f"GJ alias, not in KDWARF_PARAMS ({', '.join(KDWARF_PARAMS)})")

    s = Star()
    s.library = "mann15"
    s.pm_name = str(t["name"][idx])
    s.gj_name = i2gj.get(idx)
    s.name = s.gj_name or s.pm_name
    s.index = idx
    s.spt = str(t["spt"][idx])
    s.teff = float(t["teff"][idx])
    s.logg = float(t["logg"][idx])
    s.feh = float(t["feh"][idx])
    s.mh = float(t["mh"][idx])
    s.mass = float(t["mass"][idx])
    s.radius = float(t["radius"][idx])
    s.dist = float(t["dist"][idx])
    s.fbol = float(t["fbol"][idx]) * 1.0e-8      # erg/s/cm^2 (Mann convention)
    s.theta_mas = None
    s.dilute = (s.radius * RSUN / (s.dist * PC)) ** 2
    s.obs_split = SPLIT_MANN15
    return s


def _resolve_kdwarf(key):
    p = KDWARF_PARAMS[key]
    s = Star()
    s.library = "calib2013"
    s.pm_name = None
    s.gj_name = key
    s.name = key
    s.index = None
    s.spt = p["spt"]
    s.theta_mas = p["theta_mas"]
    theta_rad = p["theta_mas"] * 1.0e-3 / 206264.806
    s.fbol = p["fbol"]
    s.teff = (4.0 * p["fbol"] / (SIGMA_SB * theta_rad ** 2)) ** 0.25
    s.logg = p["logg"]
    s.feh = p["feh"]
    s.mh = p["feh"]
    s.mass = s.radius = s.dist = None
    s.dilute = (theta_rad / 2.0) ** 2
    s.obs_split = SPLIT_CALIB2013
    return s


def apply_theta(s, theta_mas):
    """Switch a mann15 Star to interferometric params (model-free mode):
    Teff from table Fbol + theta, R = (theta/2) d, logg from table mass."""
    theta_rad = theta_mas * 1.0e-3 / 206264.806
    s.teff = (4.0 * s.fbol / (SIGMA_SB * theta_rad ** 2)) ** 0.25
    s.radius = (theta_rad / 2.0) * s.dist * PC / RSUN
    s.logg = np.log10(GRAV * s.mass * MSUN / (s.radius * RSUN) ** 2)
    s.dilute = (s.radius * RSUN / (s.dist * PC)) ** 2
    s.theta_mas = theta_mas
    return s


# ============================================================ spectrum read
def air_to_vac(wl_A):
    """Air -> vacuum wavelengths [A], Edlen/IAU standard formula."""
    s2 = (1.0e4 / wl_A) ** 2                       # (1/um)^2
    n = (1.0 + 6.4328e-5 + 2.94981e-2 / (146.0 - s2)
         + 2.5540e-4 / (41.0 - s2))
    return wl_A * n


# Per-star registration residual [km/s] AFTER the baseline air->vac
# stretch, measured 2026-08-09 by chunk cross-correlation of the Mann
# spectrum against the star's own R=300k synthesis (SpeX-arm chunks,
# median of chi2-contrast>=1.3 chunks; airvac_check method).  Positive =
# data redward of the model rest frame; read_spectrum divides it out.
# The library is HETEROGENEOUS: most stars are air + RV-removed
# (residuals +-20 km/s = Mann's rest-shift accuracy), but GJ15A is
# clearly already-vacuum (+76 ~ +84 double-shift) and GJ644C is offset
# +31 (weakly constrained -- faint M7).  UVES GJ411 cross-check anchors
# the air+rest baseline.
VEL_RESIDUAL = {
    "GJ887": +4.7, "GJ411": -10.7, "GJ699": -21.8, "GJ905": -8.0,
    "GJ644C": +31.3, "GJ15A": +76.2,
}


def read_spectrum(star, vshift=True):
    """Observed spectrum for a Star (or name) -> (wl_A, flam_earth, err).

    Wavelengths are returned in VACUUM, in the stellar rest frame: the
    Mann distribution is (mostly) AIR wavelengths with the stellar RV
    removed, established 2026-08-09 by model-based velocity fits to
    GJ887/GJ411/GJ699/GJ905 (common -84 km/s offset despite RVs spanning
    120 km/s) plus the GJ411 UVES cross-check.  A baseline air->vac
    stretch is applied, then the measured per-star VEL_RESIDUAL is
    divided out (vshift=False skips the latter, for remeasuring).
    Residual registration after both: ~ +-10 km/s.

    mann15:    <PM name>.ascii (um, f_lambda at Earth, err); f>0 mask.
    calib2013: <GJ name>.fits primary HDU (3, N): um, f_lambda, err;
               err==0 points are BT-Settl model fill (header note) and
               are MASKED OUT along with nonpositive fluxes.  Assumed on
               the same air convention (same SNIFS/SpeX pipeline).
    """
    s = resolve(star) if isinstance(star, str) else star
    if s.library == "mann15":
        path = os.path.join(MANN_DIR, s.pm_name + ".ascii")
        if not os.path.isfile(path):
            path = os.path.join(MANN_DIR, s.pm_name.replace(" ", "") + ".ascii")
        w_um, f, e = np.loadtxt(path, unpack=True)
        ok = np.isfinite(f) & (f > 0)
    else:
        from astropy.io import fits
        d = fits.open(os.path.join(CAL_DIR, s.gj_name + ".fits"))[0].data
        w_um, f, e = d[0], d[1], d[2]
        ok = np.isfinite(f) & (f > 0) & (e > 0)
    wl = air_to_vac(w_um[ok] * 1.0e4)
    if vshift:
        dv = VEL_RESIDUAL.get(s.name, 0.0)
        wl = wl / (1.0 + dv / 2.99792458e5)
    return wl, f[ok], e[ok]


# ================================================================= LSF
def smooth_to_R(flux, model_R, obs_R):
    """Gaussian-smooth a constant-R log-lambda spectrum to resolving power obs_R."""
    from scipy.ndimage import gaussian_filter1d
    sigma_pix = model_R / (obs_R * 2.3548)
    return gaussian_filter1d(flux, sigma_pix, mode="nearest")


def smooth_mann(wave, flux, model_R, split=SPLIT_MANN15):
    """Smooth a constant-R model to the measured Mann LSF: SNIFS constant
    11 A FWHM below the splice, SpeX R=2000 above.  Pass split=
    SPLIT_CALIB2013 (or star.obs_split) for the public 2013 library."""
    from scipy.ndimage import gaussian_filter1d
    step = wave[0] / model_R                    # finest linear spacing
    wl = np.arange(wave[0], min(wave[-1], split + 200.0), step)
    fl = np.interp(wl, wave, flux)
    fsnifs = gaussian_filter1d(fl, SNIFS_DLAM_FWHM / 2.3548 / step,
                               mode="nearest")
    out = smooth_to_R(flux, model_R, OBS_R_RED)
    m = wave < split
    out[m] = np.interp(wave[m], wl, fsnifs)
    return out


# ============================================================ flux + metrics
def hnu_to_flam(wl_A, hnu):
    """Kurucz Eddington flux H_nu -> SURFACE f_lambda [erg/s/cm^2/A]."""
    return 4.0 * np.pi * hnu * C_ANG / wl_A ** 2


def read_spec_file(path, dilute=None):
    """Read a synthe .spec -> (wl_A, flam, flam_cont); surface flux unless
    a dilution factor is given (then f_lambda at Earth)."""
    w, hnu, hcont = np.loadtxt(path, unpack=True)
    f, c = hnu_to_flam(w, hnu), hnu_to_flam(w, hcont)
    if dilute is not None:
        f, c = f * dilute, c * dilute
    return w, f, c


def index_win(w, f, num, den=CAH_DEN):
    """Mean-flux ratio index (numerator window / denominator window)."""
    mn = (w > num[0]) & (w < num[1])
    md = (w > den[0]) & (w < den[1])
    return float(np.mean(f[mn]) / np.mean(f[md]))


def band_metrics(wobs, fobs, fmod):
    """Standard comparison metrics on matched (data, model) arrays."""
    out = {}
    out["integral"] = float(np.trapz(fmod, wobs) / np.trapz(fobs, wobs))
    for lo, hi, lab in [(4000, 9500, "opt_med"), (9500, 24000, "nir_med")]:
        m = (wobs > lo) & (wobs < hi)
        if m.any():
            out[lab] = float(np.median(fmod[m] / fobs[m]))
    for j, (a, b) in enumerate(TIO_BANDS):
        m = (wobs > a) & (wobs < b)
        if m.any():
            out[f"tio{j}"] = float(np.median(fmod[m] / fobs[m]))
    hi_enough = wobs.max() > 7200
    if hi_enough:
        for k, num in CAH.items():
            out[k] = index_win(wobs, fmod, num)
            out[k + "_data"] = index_win(wobs, fobs, num)
    return out


# =========================================================== run helpers
def run(cmd, cwd, logfile):
    env = dict(os.environ, ATLAS12=REPO)
    with open(logfile, "w") as log:
        proc = subprocess.run(cmd, cwd=cwd, env=env, stdout=log,
                              stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        sys.exit(f"error: {' '.join(cmd)} failed (rc={proc.returncode}); "
                 f"see {logfile}")


def iter_convergence(iter_file):
    """Return (rms, max_abs) of the flux ERROR column, final iteration."""
    blocks, cur = [], []
    with open(iter_file) as fh:
        for line in fh:
            if "log10TAU" in line:
                if cur:
                    blocks.append(cur)
                cur = []
                continue
            parts = line.split()
            if not parts or not parts[0].lstrip("-").isdigit():
                continue
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                continue
            if len(vals) >= 13:
                cur.append(vals[7])
    if cur:
        blocks.append(cur)
    if not blocks:
        return np.nan, np.nan
    err = np.array(blocks[-1])
    return float(np.sqrt(np.mean(err ** 2))), float(np.max(np.abs(err)))


def start_model(teff):
    """Nearest starting atmosphere in log Teff."""
    return min(START_MODELS, key=lambda m: abs(np.log(m[0] / teff)))


def rundir_for(star, tag=None):
    s = resolve(star) if isinstance(star, str) else star
    return os.path.join(RUN_ROOT, s.name + (f"_{tag}" if tag else ""))


# ====================================================== PHOENIX NewEra grid
# Per-metallicity 2-line-per-model GAIA-format text files, ~2 GB each.
# We index each file once (byte offset of every model's flux line), then
# extract single models on demand and cache them as npz nodes.
PHX_ZGRID = [-4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5]
PHX_LOGG  = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]


def _phx_zfile(z):
    zstr = f"+{z:.1f}" if z > 0 else f"-{abs(z):.1f}"
    return os.path.join(PHOENIX_DIR,
                        f"PHOENIX-NewEraV2-LowRes-SPECTRA.Z{zstr}.txt")


def _phx_index(z):
    """Byte-offset index {'teff_logg': flux-line offset} for one Z file."""
    key = f"phxidx{z:+.1f}"
    if key in _CACHE:
        return _CACHE[key]
    os.makedirs(PHOENIX_CACHE, exist_ok=True)
    idx_path = os.path.join(PHOENIX_CACHE, f"index.Z{z:+.1f}.json")
    if os.path.isfile(idx_path):
        with open(idx_path) as fh:
            idx = json.load(fh)
    else:
        zfile = _phx_zfile(z)
        if not os.path.isfile(zfile):
            raise FileNotFoundError(zfile)
        print(f"[phoenix] indexing {os.path.basename(zfile)} (one-time) ...")
        idx = {"models": {}, "grid": {}}
        with open(zfile, "rb") as f:
            first = True
            while True:
                hdr = f.readline()
                if not hdr:
                    break
                parts = hdr.split()
                teff, logg = float(parts[12]), float(parts[13])
                if first:
                    idx["grid"] = dict(
                        res=float(parts[7]), nwl=int(parts[8]),
                        wlstart=float(parts[9]), wlend=float(parts[10]))
                    first = False
                pos = f.tell()
                f.readline()                     # skip the flux line
                idx["models"][f"{teff:.0f}_{logg:.2f}"] = pos
        with open(idx_path, "w") as fh:
            json.dump(idx, fh)
    _CACHE[key] = idx
    return idx


def _phx_node(teff, logg, z):
    """One grid node -> (wl_A, flam_surface [erg/s/cm^2/A]), npz-cached."""
    npz = os.path.join(PHOENIX_CACHE,
                       f"t{teff:05.0f}g{logg:.2f}z{z:+.1f}.npz")
    if os.path.isfile(npz):
        d = np.load(npz)
        return d["wl"], d["flam"]
    idx = _phx_index(z)
    key = f"{teff:.0f}_{logg:.2f}"
    if key not in idx["models"]:
        raise KeyError(f"PHOENIX node t{teff:.0f} g{logg:.2f} Z{z:+.1f} "
                       f"not in grid")
    g = idx["grid"]
    with open(_phx_zfile(z), "rb") as f:
        f.seek(idx["models"][key])
        fl = np.loadtxt(io.BytesIO(f.readline()))
    wl = np.linspace(g["wlstart"], g["wlend"], g["nwl"]) * 10.0    # Angstrom
    flam = fl * 1e10 / 1e8                       # erg/s/cm^2/cm -> per A
    os.makedirs(PHOENIX_CACHE, exist_ok=True)
    np.savez_compressed(npz, wl=wl, flam=flam.astype(np.float32))
    return wl, flam


def _bracket(grid, x):
    """(lo, hi, w) with x = (1-w)*lo + w*hi; clamps to grid ends."""
    g = sorted(grid)
    if x <= g[0]:
        return g[0], g[0], 0.0
    if x >= g[-1]:
        return g[-1], g[-1], 0.0
    hi = min(v for v in g if v >= x)
    lo = max(v for v in g if v <= x)
    w = 0.0 if hi == lo else (x - lo) / (hi - lo)
    return lo, hi, w


def phoenix_newera(teff, logg, mh, verbose=True):
    """PHOENIX NewEraV2 LowRes surface spectrum at (teff, logg, [M/H]) by
    trilinear interpolation in log f_lambda over the bracketing grid nodes
    (100 K / 0.5 dex / 0.5 dex spacing; clamped at grid edges).

    Returns (wl_A, flam) with wl 2500-25000 A at 0.1 A; flam is surface
    f_lambda [erg/s/cm^2/A].
    """
    zlo, zhi, wz = _bracket(PHX_ZGRID, mh)
    pieces = []
    for z, wgt_z in [(zlo, 1.0 - wz), (zhi, wz)]:
        if wgt_z == 0.0 and zlo != zhi:
            continue
        idx = _phx_index(z)
        teffs = sorted({float(k.split("_")[0]) for k in idx["models"]})
        loggs = sorted({float(k.split("_")[1]) for k in idx["models"]
                        if float(k.split("_")[1]) in PHX_LOGG})
        tlo, thi, wt = _bracket(teffs, teff)
        glo, ghi, wg = _bracket(loggs, logg)
        for t, wgt_t in [(tlo, 1.0 - wt), (thi, wt)]:
            for g, wgt_g in [(glo, 1.0 - wg), (ghi, wg)]:
                wgt = wgt_z * wgt_t * wgt_g
                if wgt == 0.0:
                    continue
                wl, flam = _phx_node(t, g, z)
                pieces.append((wgt, t, g, z, wl, flam))
    wl0 = pieces[0][4]
    logf = np.zeros_like(pieces[0][5], dtype=float)
    for wgt, t, g, z, wl, flam in pieces:
        assert wl.shape == wl0.shape
        logf += wgt * np.log(np.maximum(flam, 1e-300))
    if verbose:
        used = ", ".join(f"t{t:.0f}/g{g:.1f}/z{z:+.1f}({wgt:.2f})"
                         for wgt, t, g, z, _, _ in pieces)
        print(f"[phoenix] ({teff:.0f}, {logg:.2f}, {mh:+.2f}) <- {used}")
    return wl0, np.exp(logf)


if __name__ == "__main__":
    # smoke: resolve every ladder star and print the resolved params
    for name, note in STAR_LADDER:
        s = resolve(name)
        print(f"{s.name:<8} {s.library:<10} pm={s.pm_name or '-':<18} "
              f"Teff={s.teff:6.0f} logg={s.logg:.2f} [Fe/H]={s.feh:+.2f} "
              f"split={s.obs_split:.0f}  | {note[:60]}")
