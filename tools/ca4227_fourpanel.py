#!/usr/bin/env python3
"""Ca I 4227 profile across the Teff sequence: M7 -> M0.5 -> K5 -> Sun.

Four panels over the blue-depression region, each showing the observed
spectrum, the standard Voigt model, and the modified resonance profile
(Jones et al. 2023 form at p_s = 8e-6, p_x = 0.5, |dlam| < 750 A,
with p_p = 1 so the Doppler core is the true one -- see the CA4227_ block
in mod_parameters for why Jones' p_p = 1e-4 cannot be used).

Everything is put on a common 11 A FWHM (the measured SNIFS LSF) so the
four panels are directly comparable; models keep their native R = 300000
sampling -- they are convolved, never resampled.

Usage:  python3 ca4227_fourpanel.py [OUT.png]
"""
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mann_lib as ml

RSYN = 300000.0
FWHM = ml.SNIFS_DLAM_FWHM                  # 11 A, common display resolution
SOLAR_ATLAS = os.path.expanduser(
    "~/kurucz/SolarAtlas/Kurucz05/irradthuwl.br_05nm")
RSUN_AU = ml.RSUN / 1.495978707e13         # solar dilution^(1/2)

C_DATA, C_VOIGT, C_MOD = "0.15", "#0072B2", "#CC79A7"
XLO, XHI = 0.385, 0.52                     # um


def smooth_const_fwhm(w, f, fwhm=FWHM):
    """Gaussian of constant FWHM in Angstrom, on the input sampling.

    The model grid is logarithmic, so convert to a pixel sigma using the
    local step; over 3850-5200 A that varies by 3 per cent, small enough
    to use the band-centre value.
    """
    step = np.median(np.diff(w))
    return gaussian_filter1d(f, fwhm / 2.3548 / step, mode="nearest")


def model_curve(specpath):
    w, f, _ = ml.read_spec_file(specpath)
    m = (w > 3700) & (w < 5400)
    return w[m], smooth_const_fwhm(w[m], f[m])


def star_panel(name, base_tag, mod_tag):
    s = ml.resolve(name)
    wobs, fobs, _ = ml.read_spectrum(s)
    fsurf = fobs / s.dilute
    d = (wobs, fsurf)
    v = model_curve(os.path.join(ml.rundir_for(s, base_tag), s.name + ".spec"))
    m = model_curve(os.path.join(ml.rundir_for(s, mod_tag), s.name + ".spec"))
    title = (f"{s.name}   {s.spt},  $T_{{\\rm eff}}$ = {s.teff:.0f} K")
    return title, d, v, m


def solar_panel():
    """Kurucz 2005 Thuillier-calibrated irradiance atlas vs the solar model.

    Atlas: W/m^2/nm at 1 AU on 0.001 nm sampling with 0.05 nm FWHM
    broadening; 1 W/m^2/nm = 100 erg/s/cm^2/A.  Converted to SURFACE flux
    with the geometric dilution (Rsun/AU)^2 so the panel matches the stars.
    """
    wl, fl = [], []
    with open(SOLAR_ATLAS) as fh:
        for line in fh:
            p = line.split()
            if len(p) != 2:
                continue
            try:
                a, b = float(p[0]), float(p[1])
            except ValueError:
                continue
            if 370.0 < a < 540.0:
                wl.append(a * 10.0)
                fl.append(b * 100.0 / RSUN_AU ** 2)
    wl, fl = np.array(wl), np.array(fl)
    d = (wl, smooth_const_fwhm(wl, fl))
    root = os.path.join(ml.REPO, "workdir", "solar_check")
    v = model_curve(os.path.join(root, "ca_base", "sun.spec"))
    m = model_curve(os.path.join(root, "ca8e6_pp1", "sun.spec"))
    return "Sun   G2V,  $T_{\\rm eff}$ = 5777 K", d, v, m


def main():
    out = (sys.argv[1] if len(sys.argv) > 1 else
           os.path.join(ml.RUN_ROOT, "ca4227_fourpanel.png"))
    panels = [star_panel("GJ644C", "w_base", "w_pp1"),
              star_panel("GJ887", "ca_base", "ca8p1"),
              star_panel("GJ820A", "ca_base", "ca8p1"),
              solar_panel()]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix",
                         "font.size": 13, "axes.grid": False,
                         "figure.dpi": 600})

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    for ax, (title, d, v, m) in zip(axes.ravel(), panels):
        for lab, (w, f), col, z, lw in [
                ("observed", d, C_DATA, 0, 1.1),
                ("standard Voigt", v, C_VOIGT, 2, 1.3),
                (r"modified Ca I, $p_s$=8e-6", m, C_MOD, 3, 1.3)]:
            k = (w / 1e4 > XLO) & (w / 1e4 < XHI)
            ax.plot(w[k] / 1e4, f[k], color=col, lw=lw, label=lab, zorder=z)
        ax.set_xlim(XLO, XHI)
        ax.set_ylim(bottom=0)
        ax.set_title(title, fontsize=13)
        ax.minorticks_on()
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.tick_params(which="major", length=7)
        ax.tick_params(which="minor", length=3.5)
    for ax in axes[1]:
        ax.set_xlabel(r"wavelength  [$\mu$m]")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$f_\lambda^{\rm surf}$"
                      r"  [erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$]")
    axes[0, 0].legend(frameon=False, fontsize=11, loc="upper left")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    print("wrote", out)

    # --- band table -------------------------------------------------
    print(f"\n{'star':>10} {'window':>12} {'Voigt/obs':>10} {'mod/obs':>9}")
    for title, d, v, m in panels:
        nm = title.split()[0]
        for lo, hi in [(4000, 4300), (4300, 4600), (4600, 5000)]:
            k = (d[0] > lo) & (d[0] < hi)
            do = np.mean(d[1][k])
            kv = (v[0] > lo) & (v[0] < hi)
            km = (m[0] > lo) & (m[0] < hi)
            print(f"{nm:>10} {lo:5d}-{hi:5d} {np.mean(v[1][kv])/do:10.3f} "
                  f"{np.mean(m[1][km])/do:9.3f}")


if __name__ == "__main__":
    main()
