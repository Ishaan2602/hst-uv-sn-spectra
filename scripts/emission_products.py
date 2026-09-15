#!/usr/bin/env python3
# per-source EMISSION products for the public repo (Mg II 2800 + Lya 1216 line flux vs phase).
# measures the continuum-subtracted line flux at every epoch (specutils line_flux, bostroem method),
# fits 4 candidate profile models (gaussian/lorentzian/skew/kwok), picks the best by BIC, and writes a per-SN JSON
# sidecar + a catalog summary + per-epoch scrutiny plots so a follower can eyeball the continuum
# placement / notch interp / shell fit behind every reported flux.
#
# this is the EMISSION thread only (late-time CSM-interaction + photospheric peaks). the absorption /
# ISM foreground work (columns, N(HI), metallicity) lives in absorption_products.py - same 2800A
# wavelength, opposite physics, kept separate on purpose.
#
# usage:  python emission_products.py [SN ...]      (no args -> every SN with spectra under output)

import os, glob, csv, re, json, datetime
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.special import erf
import astropy.units as u
from astropy.io import fits
from dust_extinction.parameter_averages import F19
try:
    from specutils import Spectrum
except ImportError:
    from specutils import Spectrum1D as Spectrum
from specutils import SpectralRegion
from specutils.analysis import line_flux
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import time

def _savefig(fig, path, dpi=110):
    # /mnt/c writes from the WSL Windows-python are intermittently flaky (OSError errno22 on the FS bridge);
    # a bounded retry keeps one transient flake from killing a whole catalog run.
    for i in range(4):
        try:
            fig.savefig(path, dpi=dpi); return
        except OSError:
            if i == 3:
                raise
            time.sleep(0.4)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from paths import OUT, CATALOG
from paths import EMISSION_SUMMARY as SUMMARY
from paths import LYA_NHI_SUMMARY, REFERENCE
import lya_nhi

C_KMS = 2.99792458e5
MG2, LYA = 2799.94, 1215.67
FLAM = u.Unit("erg cm-2 s-1 AA-1")

# --- catalog ---------------------------------------------------------------------------------------
cat = {}
with open(CATALOG) as fh:
    for row in csv.DictReader(fh):
        cat[row["name"].upper()] = row
def _catf(sn, key):
    v = cat[sn.upper()].get(key)
    return float(v) if v not in (None, "") else 0.0     # peripheral SNe can have blank z/ebv in the catalog
zof = lambda sn: _catf(sn, "z")
ebvof = lambda sn: _catf(sn, "ebv")                      # MW foreground only

# host reddening comes from the AUTHORITATIVE reference file, NOT the catalog mirror (which goes stale when
# catalog_clean.py isn't re-run after a host_ebv.csv edit - that bug shipped host=0 for AT2022ACKO/LMC once).
from paths import host_ebv_map
_HOST = host_ebv_map()
for _sn, (_hv, _he, _src) in _HOST.items():             # loud check: warn if the catalog mirror has drifted
    _cv = cat.get(_sn, {}).get("host_ebv")
    if _cv not in (None, "") and abs(float(_cv) - _hv) > 1e-6:
        print(f"  WARN host_ebv drift: {_sn} reference={_hv} catalog={_cv} (run catalog host-sync)")
hostof = lambda sn: _HOST.get(sn.upper(), (0.0, None, None))[0]
hostsrc = lambda sn: (_HOST.get(sn.upper()) or (0.0, None, "none (MW-only)"))[2] or "none (MW-only)"
tnstype = lambda sn: cat[sn.upper()].get("tns_type") or ""


# per-SN logN(HI) for the Lya intrinsic-flux correction (F_ismcorr): curated reference/ism_columns.csv preferred,
# else the automated damped-Lya fit in catalog/lya_nhi_summary.csv. single source of truth shared with the ISM thread.
_ISM_COLUMNS = os.path.join(REFERENCE, "ism_columns.csv")
def _build_nhi():
    m = {}
    if os.path.exists(LYA_NHI_SUMMARY):
        with open(LYA_NHI_SUMMARY) as fh:
            for r in csv.DictReader(fh):
                try:
                    lN = float(r["logN_HI"]); er = float(r.get("logN_HI_syst_vabs") or 0.0) or 0.3
                except (ValueError, KeyError, TypeError):
                    continue
                m[r["sn"].upper()] = (lN, er, "lya_nhi_summary (automated damped-Lya fit)")
    with open(_ISM_COLUMNS) as fh:                              # curated overrides the automated
        for row in csv.reader(fh):
            if len(row) < 4 or row[0].startswith("#") or row[1].strip() != "logN(HI)":
                continue
            try:
                lN = float(row[2]); er = float(row[3]) if row[3] else 0.2
            except ValueError:
                continue
            m[row[0].strip().upper()] = (lN, er, "ism_columns (curated)")
    return m
_NHI = _build_nhi()
nhiof = lambda sn: _NHI.get(sn.upper())                        # (logN, logN_err, source) or None


# --- emission machinery (ported verbatim from emission_investigation2.ipynb) ------------------------
_f19 = F19(Rv=3.1)
def _dered1(wave_A, ebv):
    fac = np.ones_like(wave_A, float)
    if ebv > 0:
        m = (wave_A > 1150) & (wave_A < 33333)
        fac[m] = 1.0 / _f19.extinguish(wave_A[m] * u.AA, Ebv=ebv)
    return fac
def deredden(wrest_A, z, mw_ebv, host_ebv=0.0):
    return _dered1(wrest_A * (1 + z), mw_ebv) * _dered1(wrest_A, host_ebv)   # MW at obs wvl + host at rest


def load_spec(path, z):
    # read a _native.txt 3-column spectrum (obs-frame wvl, flux, err); shift wvl to rest frame
    d = np.loadtxt(path, comments="#")
    w, f = d[:, 0], d[:, 1]
    e = d[:, 2] if d.shape[1] > 2 else np.full_like(f, np.nan)
    return w / (1.0 + z), f, e


AIR_HALF = 3000.0     # km/s, local fit window half-width around geocoronal Lya (v=-cz)
AIR_DEG = 2           # local smooth-SN poly degree in the airglow decomposition (deg2==deg3, deg1 under-fits)

def _airglow_bg(native_path, wrest, z):
    # coadd the sibling per-exposure x1d BACKGROUND onto the native REST grid. geocoronal Lya fills the slit so
    # it lands in the calstis background array regardless of z - a z-independent airglow locator (bostroem note 1).
    # the reduction keeps the x1d next to each native; write_1d just drops the column. None if no sibling x1d.
    xs = sorted(glob.glob(os.path.join(os.path.dirname(native_path), "*_x1d.fits")))
    if not xs:
        return None
    S = []
    for x in xs:
        t = fits.getdata(x, 1)
        wx = np.asarray(t["WAVELENGTH"][0], float) / (1 + z)      # x1d is observed-frame -> rest
        S.append(np.interp(wrest, wx, np.asarray(t["BACKGROUND"][0], float), left=0.0, right=0.0))
    return np.median(S, 0)

def _airglow_subtract(w, f, z, bg, half=AIR_HALF, deg=AIR_DEG):
    # remove geocoronal Lya. two regimes:
    #  (1) STIS low-res MAMA (G140L/G230L): the airglow shows as a spike in the sibling x1d BACKGROUND -> additive
    #      decomposition, model = smooth SN (local poly) + alpha*airglow(bg-shaped), subtract ONLY alpha*bg so the
    #      real SN emission UNDER the airglow survives (a straight bridge inflates the 2023ixf shell / eats 2005ip).
    #  (2) echelle / COS / flat-bg: the airglow is NOT in the background -> fall back to a narrow velocity notch at
    #      v=-cz (the pre-bg behavior), only when the airglow is separated from the SN Lya core (|vgeo|>800).
    v = (w - LYA) / LYA * C_KMS
    vgeo = -C_KMS * z / (1 + z)
    if bg is not None:
        reg = (v > vgeo - half) & (v < vgeo + half) & np.isfinite(bg) & np.isfinite(f)
        if reg.sum() >= 8 and np.nanmax(bg[reg]) > 0:
            base = np.nanmedian(bg[reg]); mad = np.nanmedian(np.abs(bg[reg] - base)) + 1e-30
            if (np.nanmax(bg[reg]) - base) > 5 * mad:            # a real airglow spike is present in the bg
                bn = np.clip(bg, 0.0, None) / (np.nanmax(bg[reg]) + 1e-30)
                x = (v - vgeo) / 1000.0
                A = np.column_stack([np.ones(reg.sum())] + [x[reg] ** p for p in range(1, deg + 1)] + [bn[reg]])
                coef, *_ = np.linalg.lstsq(A, f[reg], rcond=None)
                cleaned = f[reg] - coef[-1] * bn[reg]
                # the bg template is slightly broader/shifted vs the real net airglow residual, so a free alpha
                # mis-subtracts the core: a positive airglow leaves a negative over-sub spike, a negative airglow
                # (over-subtracted dip) gets OVER-filled into a spurious positive peak at v=-cz (fake detection).
                # in the airglow core (bn>0.1) hold the cleaned flux to the smooth SN model +-2*noise: real Lya is
                # BROAD and lives in `smooth`, only the narrow at-vcz artifact is capped.
                smooth = (A @ coef) - coef[-1] * bn[reg]
                sig = float(np.nanstd(f[reg] - (A @ coef)))
                cc = bn[reg] > 0.1
                cleaned[cc] = np.clip(cleaned[cc], smooth[cc] - 2.0 * sig, smooth[cc] + 2.0 * sig)
                out = f.copy(); out[reg] = cleaned
                return out
    if abs(vgeo) > 800:                                          # fallback velocity notch (echelle / COS / flat bg)
        emisw = (v > -8000) & (v < 5000)
        notch = emisw & (np.abs(v - vgeo) < 350)
        if notch.any() and (emisw & ~notch).sum() > 2:
            out = f.copy()
            out[notch] = np.interp(w[notch], w[emisw & ~notch], f[emisw & ~notch])
            return out
    return f


def phase_of(p):
    m = re.search(r"day([0-9.]+)", p)
    return float(m.group(1)) if m else np.nan


def _central_notch(v, fc, emis, lam0, sigma=None, z=0.0):
    # conditional narrow-absorption notch: mask a genuine central Mg II absorption (2796/2803 doublet, near
    # systemic) ONLY where it actually cuts into emission, over its real extent. replaces the old unconditional
    # +-650 floor, which over-masked narrow-line IIn and invented flux over a photospheric P-Cygni trough.
    # gate: both flanks must be in emission (>2 sigma) so we never repair a photospheric absorption; a real
    # dip = the central min sits below the outer-flank emission envelope by >2 sigma. detection is local-min
    # based (a straight-floor over-masked). geocoronal-Lya airglow is handled upstream now (see _airglow_subtract).
    notch = np.zeros_like(v, bool)
    if sigma is None:
        sigma = float(np.nanstd(fc[emis])) or (0.05 * np.nanmax(fc[emis]))
    search = 1300.0 if lam0 == MG2 else 1600.0
    inner = emis & (np.abs(v) < search); outer = emis & (np.abs(v) >= search)
    lsh = fc[emis & (v > -1700) & (v < -search)]; rsh = fc[emis & (v > search) & (v < 1700)]
    if inner.sum() >= 4 and outer.sum() >= 4 and lsh.size and rsh.size \
            and np.nanmedian(lsh) > 2 * sigma and np.nanmedian(rsh) > 2 * sigma:
        env = np.interp(v, v[outer], fc[outer])         # emission shape with the center bridged (follows asymmetry)
        core = emis & (np.abs(v) < 600)                 # a real absorption min, if any, sits near systemic
        if core.any():
            ci = np.where(core)[0][np.argmin(fc[core])]
            lo = emis & (v > v[ci] - 800) & (v < v[ci] - 300); hi = emis & (v > v[ci] + 300) & (v < v[ci] + 800)
            valley = lo.any() and hi.any() and fc[ci] < np.nanmedian(fc[lo]) - sigma and fc[ci] < np.nanmedian(fc[hi]) - sigma
            if valley and fc[ci] < env[ci] - 2 * sigma:   # a genuine absorption VALLEY (below both sides + the envelope)
                notch[ci] = True                          # grow over the contiguous below-envelope trough (capped)
                j = ci - 1
                while j >= 0 and emis[j] and fc[j] < env[j] - 0.5 * sigma and abs(v[j] - v[ci]) < 1400:
                    notch[j] = True; j -= 1
                j = ci + 1
                while j < len(v) and emis[j] and fc[j] < env[j] - 0.5 * sigma and abs(v[j] - v[ci]) < 1400:
                    notch[j] = True; j += 1
    return notch


def bostroem_flux(w, f, e=None, lam0=MG2, vline=(-10000, 6000), vcont=(-16000, 16000), vabs=350, deg=1, smooth=0, z=0.0):
    # continuum-subtract and integrate the line flux (bostroem+2026 method). fit a deg-1 continuum on the vcont
    # flanks, then integrate the residual over vline. returns a dict with F, sigF, continuum, and diagnostic arrays.
    v = (w - lam0) / lam0 * C_KMS
    ff = savgol_filter(f, smooth, 2) if smooth > 2 else f.copy()
    emis = (v > vline[0]) & (v < vline[1]); contfit = (v > vcont[0]) & (v < vcont[1]) & ~emis
    cont = np.polyval(np.polyfit(w[contfit], ff[contfit], deg), w); fc = ff - cont
    sigma = float(np.nanstd(fc[contfit])) or (0.05 * np.nanmax(fc[emis]))    # continuum scatter = notch noise scale
    notch = _central_notch(v, fc, emis, lam0, sigma=sigma, z=z); fc[notch] = np.interp(w[notch], w[emis & ~notch], fc[emis & ~notch])
    lo, hi = lam0 * (1 + vline[0] / C_KMS), lam0 * (1 + vline[1] / C_KMS)
    reg = SpectralRegion(lo * u.AA, hi * u.AA)
    F = line_flux(Spectrum(spectral_axis=w * u.AA, flux=fc * FLAM), reg).to("erg cm-2 s-1").value
    sigF = np.nan
    if e is not None:                          # pixel-noise flux error over the same window (cont-fit err not included)
        win = (w >= lo) & (w <= hi); dw = np.gradient(w)
        var = (e[win] * dw[win]) ** 2
        if np.isfinite(var).any():
            sigF = float(np.sqrt(np.nansum(var)))
    return dict(F=F, sigF=sigF, cont=cont, fc=fc, v=v, ff=ff, emis=emis, contfit=contfit, notch=notch, deg=deg)


def kwok_shell(v, A, mu, fwhm, vc, vin):
    # optically thin expanding shell (kwok 1994): gaussian emissivity sphere with an off-center spherical hole.
    # produces the flat-topped boxy profile seen in broad CSM shells; mu = centroid offset (blueshift).
    sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
    g = np.exp(-0.5 * ((v - mu) / sig) ** 2)
    vh = mu + vc
    supp = np.ones_like(v, float)
    inside = np.abs(v - vh) < vin
    supp[inside] = np.exp(-0.5 * (vin ** 2 - (v[inside] - vh) ** 2) / sig ** 2)
    return A * g * supp


# kwok fit bounds on [A, mu, fwhm, vc, vin] (A on the O(1)-normalized flux). shared so we can flag pegging.
KWOK_LO = [0, -8000, 3000, -3000, 1000]
KWOK_HI = [3, 2000, 15000, 6000, 9000]


def fit_kwok(vv, ff):
    # raw flux ~1e-15 makes curve_fit's squared resid sit under gtol -> it "converges" at p0 and never moves
    # the shape. fit the flux normalized to O(1), then scale the amplitude back.
    A = np.nanmax(ff)
    p = curve_fit(kwok_shell, vv, ff / A, p0=[1, -2995, 7580, 1750, 5000], bounds=(KWOK_LO, KWOK_HI), maxfev=30000)[0]
    p = np.asarray(p, float); p[0] *= A
    return p


# candidate profile models (the flux stays the model-free direct integral; these only pick the SHAPE label).
def gaussv(v, A, mu, sig):                        # symmetric peak (k=3 free params)
    return A * np.exp(-0.5 * ((v - mu) / sig) ** 2)


def lorentzian(v, A, mu, gam):                    # narrow-line IIn scattering wings (k=3)
    return A * gam ** 2 / ((v - mu) ** 2 + gam ** 2)


def skewg(v, A, mu, sig, al):                     # mild asymmetry (k=4)
    t = (v - mu) / sig
    return A * np.exp(-0.5 * t ** 2) * (1 + erf(al * t / np.sqrt(2)))


# shape-model parameter names + bounds. width floors are 100 km/s (were 500/300): a fixed 500 km/s
# floor pegged genuinely narrow lines on COS/echelle (SN2010jl Lya d595 on G130M) and a floor is not what keeps a
# low-res fit honest anyway -- the generic peg guard in _fit_profiles is. kwok keeps its 3000 km/s FWHM floor
# because the shell model is only meaningful for a broad shell (below that it degenerates into a narrow gaussian).
MODEL_SPECS = {
    "gaussian":   (gaussv,     ["A", "mu", "sig"],              [1, -3000, 4000],          [0, -9000, 100],       [5, 4000, 12000]),
    "lorentzian": (lorentzian, ["A", "mu", "gam"],              [1, -2000, 3000],          [0, -9000, 100],       [5, 4000, 12000]),
    "skew":       (skewg,      ["A", "mu", "sig", "al"],        [1, -3000, 4000, -2],      [0, -9000, 100, -20],  [5, 4000, 12000, 20]),
    "kwok":       (kwok_shell, ["A", "mu", "fwhm", "vc", "vin"], [1, -2995, 7580, 1750, 5000], KWOK_LO,          KWOK_HI),
}
PEG_TOL = 0.02          # a parameter within 2% of its bound span of either bound counts as pegged
DBIC_STRONG = 6.0       # kass & raftery 'strong' evidence; a k>3 model must beat the best k=3 model by this to win
COH_MIN = 6.0           # a k>3 model may only win BY DEFAULT (no unpegged k=3 competitor) if the profile is this coherent

# instrumental LSF sigma sets the width floor for the narrow models. high-resolution gratings (COS medium, STIS
# medium) resolve genuinely narrow lines that a fixed 100 km/s floor would peg and then bar (SN2010jl Lya d595 on
# G130M was a real narrow line pegging an arbitrary floor); low-resolution (G140L, G230L, CCD) never produce a line
# narrower than ~100 km/s. kwok keeps its 3000 km/s FWHM floor (the shell model is only meaningful for a broad shell).
WIDTH_FLOOR = {"G130M": 15, "G160M": 15, "G140M": 15, "G185M": 15, "G225M": 15, "G285M": 15, "G230M": 15, "G230MB": 15}
def _width_floor(instr):
    return WIDTH_FLOOR.get(instr, 100)

# extra p0 seeds for the non-convex models (skew, kwok) -- a single start lands in a different local minimum on some
# epochs (GGI d232 kwok<->skew flipped when only a bound moved). multi-start + keep-lowest-RSS makes the fit stable.
P0_SEEDS = {
    "skew": [[1, -3000, 4000, -2], [1, -5000, 3000, 2], [1, -1500, 5000, 0]],
    "kwok": [[1, -2995, 7580, 1750, 5000], [1, -4500, 5000, 1000, 3000], [1, -1000, 9000, 3000, 6000]],
}


def _lo_with_floor(name, lo, wfloor):
    # replace the width lower bound (sig/gam, index 2) with the instrument-aware floor for the narrow models.
    lo = list(lo)
    if name in ("gaussian", "lorentzian", "skew"):
        lo[2] = wfloor
    return lo


def _pegged(name, p, wfloor=100):
    # which SHAPE parameters (everything but the amplitude, which is rescaled after the fit) sit on a bound.
    # a pegged parameter means the optimizer wanted to leave the allowed region: the shape is not constrained
    # by the data inside the model's domain, so that model must not be *selected* (it stays in models{} for the record).
    _, names, _, lo, hi = MODEL_SPECS[name]; lo = _lo_with_floor(name, lo, wfloor); out = []
    for nm, x, l, h in zip(names[1:], p[1:], lo[1:], hi[1:]):
        span = h - l
        if x <= l + PEG_TOL * span:
            out.append(f"{nm}@lo")
        elif x >= h - PEG_TOL * span:
            out.append(f"{nm}@hi")
    return out


def _fit_models(vv, ff, wfloor=100):
    # every fit on O(1)-normalized flux then amplitude scaled back. fitting raw ~1e-15 froze the bounded
    # lorentzian at p0 (curve_fit hits its gradient tol at iter 0); normalizing is the same fix kwok used.
    # skew/kwok are multi-started (P0_SEEDS) and the lowest-RSS fit is kept to remove local-minimum flips.
    A = np.nanmax(ff); out = {}
    for name, (fn, _names, p0, lo, hi) in MODEL_SPECS.items():
        lo = _lo_with_floor(name, lo, wfloor)
        best = None
        for seed in P0_SEEDS.get(name, [p0]):
            try:
                p = curve_fit(fn, vv, ff / A, p0=seed, bounds=(lo, hi), maxfev=40000)[0]
            except Exception:
                continue
            rss = float(np.nansum((ff / A - fn(vv, *np.asarray(p, float))) ** 2))
            if best is None or rss < best[1]:
                best = (np.asarray(p, float), rss)
        if best is not None:
            p = best[0]; p[0] *= A; out[name] = (fn, p)
    return out


def _vphot_normalized(w, f, lam0):
    # APPROXIMATE photospheric velocity for a P-Cygni epoch: deg-2 pseudo-continuum through the flanks
    # (feature masked), normalize, take the blueward absorption minimum. the UV is a forest of overlapping
    # Fe/Mg absorption, so this is a rough diagnostic, not a precision measurement.
    v = (w - lam0) / lam0 * C_KMS
    reg = (v > -22000) & (v < 12000); feat = (v > -16000) & (v < 7000)
    cf = reg & ~feat
    if cf.sum() < 5:
        return None
    cont = np.polyval(np.polyfit(w[cf], f[cf], 2), w)
    with np.errstate(invalid="ignore", divide="ignore"):
        fn = f / cont
    blue = (v > -18000) & (v < -1000) & np.isfinite(fn)
    if not blue.any():
        return None
    return int(round(float(v[blue][np.argmin(fn[blue])])))


def _model_ic(vv, ff, fn, p, sigma):
    # BIC/AICc use RSS directly (no sigma); redchi2 needs the noise (continuum scatter, rough for bright profiles).
    r = ff - fn(vv, *p); rss = float(np.nansum(r ** 2)); N = len(vv); k = len(p)
    return {"params": [round(float(x), 2) for x in p], "k": k,
            "bic": round(N * np.log(rss / N) + k * np.log(N), 1),
            "aicc": round(N * np.log(rss / N) + 2 * k + 2 * k * (k + 1) / max(N - k - 1, 1), 1),
            "redchi2": round(rss / (sigma ** 2 * max(N - k, 1)), 2)}


def _profile_coherence(v, ff):
    # smoothed-profile amplitude / pixel residual scatter. HIGH = a coherent absorption structure (a genuine
    # photospheric P-Cygni), LOW = noise / no clean line. amplitude/sigma was backwards (noise scored high),
    # coherence is the honest confidence for a deep-central-absorption epoch. splits: >=12 photospheric,
    # 6-12 marginal (low-SNR gray zone), <6 low_snr.
    if len(ff) < 7:
        return 0.0
    dv = np.median(np.abs(np.diff(v))) or 1.0
    win = max(5, min(int(1200 / dv) // 2 * 2 + 1, (len(ff) // 2) * 2 - 1))
    if win < 5 or len(ff) < win:
        return 0.0
    sm = savgol_filter(ff, win, 2); resid = float(np.nanstd(ff - sm))
    return float((np.nanmax(sm) - np.nanmin(sm)) / resid) if resid else 0.0



# Mg II emission window (km/s). the ASYMMETRIC default (more blue, since the CSM shell is blueshifted) is the
# right call; it reproduces bostroem+2026 table 4 to 1.00 +- 0.03 on every 2023ixf epoch AND on SN2024ggi once the
# comparison is done on OBSERVED flux; her table 4 is not dereddened.
# the old per-SN override {"SN2024GGI": (-6000, 4000)} was fitted to a phantom 2.6x offset that was really the
# dereddening convention; it captured only 41-64% of the GGI flux and is REMOVED. the dict stays as the mechanism.
# NOTE: a data-driven adaptive-tightening window was tried and REJECTED - it is noise-fragile (it cut real emission on
# 2023ixf d66 429->312); bostroem's own 99%-enclosed-flux edges are likewise set by her manual outer limits (A3), and
# the window-spill census over every clean epoch (A4) shows <= 10% outside the default. so: fixed default, no override.
MG2_WINDOW = {}
_MG2_DEFAULT = (-10000, 6000)

# gratings the emission thread measures on, and the echelle modes it skips. echelle spectra
# (E230M/E140M/E230H/E140H) sample only the Mg II / Lya CORE at R~30000 in narrow orders whose continuum is
# poorly defined over a +-10000 km/s window; the 1998S E230M "emission" epochs were really the narrow ISM Mg II
# absorption seen against the broad SN line. those spectra belong to the ISM thread.
EMIS_GRATINGS = ("G230LB", "G230L", "G140L", "G230M", "G130M", "G160M", "G185M", "G225M", "G285M", "G230MB")
ECHELLE_GRATINGS = ("E230M", "E140M", "E230H", "E140H")

# Lya blue edge -8000 -> -10000, shared with Mg II (bostroem fig 5, cool dense shell radiates both). the
# C III 1176 scale-check (emission_investigation3 A8/A9) found no distinct C III feature in -11000..-9000 for any
# science epoch; the wider edge adds real Lya blue-wing toe (spill 1.00-1.13) and de-inflates the window-perturbation
# systematic (the -8000 edge sat on the steep toe, which had flagged the strong 2023ixf d183 line unreliable). red
# edge stays +5000 (N V 1238/1242 walk in beyond that).
_LYA_DEFAULT = (-10000, 5000)


def _edges99(v, fc, wide, frac=0.99, win_px=5):
    # bostroem table-5 method: smooth the cont-sub profile, integrate over a generous window, step inward from each
    # side to the velocities enclosing `frac` of that flux. DIAGNOSTIC ONLY - A3 showed these edges are set by the
    # outer manual window (the last 0.5%/side sits in noise), not stable; reported for comparison, not used for flux.
    m = (v > wide[0]) & (v < wide[1]) & np.isfinite(fc)
    vv, ff = v[m], fc[m]
    if len(ff) <= win_px:
        return None, None
    ff = savgol_filter(ff, win_px, 2)
    dv = np.gradient(vv); tot = float(np.sum(ff * dv))
    if tot <= 0:
        return None, None
    cb = np.cumsum(ff * dv); cr = np.cumsum((ff * dv)[::-1])[::-1]
    ib = int(np.argmax(cb >= (1 - frac) / 2 * tot))
    ir = len(vv) - 1 - int(np.argmax(cr[::-1] >= (1 - frac) / 2 * tot))
    return float(vv[ib]), float(vv[ir])


def _diag_keys(bf, w, f, lam0, win, deg, z, wide):
    # per-epoch window diagnostics, DIAGNOSTIC ONLY: the 99%-enclosed-flux edges (bostroem table 5 method)
    # and the spill = F(default widened 3000 km/s each side) / F(default). spill is the cleaner one: ~1 means the
    # default window captures the line, >> 1 means real flux sits outside it.
    vb, vr = _edges99(bf["v"], bf["fc"], wide)
    outer = (win[0] - 3000, win[1] + 3000)
    vc = (min(-16000, outer[0] - 3000), max(16000, outer[1] + 3000))
    Fo = bostroem_flux(w, f, lam0=lam0, vline=outer, vcont=vc, deg=deg, z=z)["F"]
    Fd = bf["F"]
    return {"v_blue99": round(vb) if vb is not None else None,
            "v_red99": round(vr) if vr is not None else None,
            "spill_frac": round(Fo / Fd, 3) if (Fd and Fd > 0 and Fo > 0) else None}   # None when the wider window is absorption-dominated (uninterpretable as a spill)


def _lya_ismcorr_flux(w, bf, logN, b=25.0, Tmin=0.3, vline=_LYA_DEFAULT):
    # divide the airglow-subtracted, continuum-subtracted Lya profile by the foreground ISM H I transmission
    # exp(-tau) from the measured N(HI); bridge only the black core (T<=Tmin). returns (F, core_halfwidth) or None.
    T = np.exp(-lya_nhi.lya_tau(w, logN, b, 0.0))
    fc = bf["ff"] - bf["cont"]
    good = T > Tmin; hole = ~good
    if hole.all():
        return None
    fcc = np.where(good, fc / np.clip(T, Tmin, 1.0), np.nan)
    if hole.any():
        fcc[hole] = np.interp(w[hole], w[~hole], fcc[~hole])
    v = bf["v"]
    hw = 0.5 * float(v[hole].max() - v[hole].min()) if hole.any() else 0.0
    m = (v > vline[0]) & (v < vline[1])
    return float(np.trapezoid(fcc[m], w[m])), hw


def _ismcorr_keys(w, bf, sn, vline=_LYA_DEFAULT):
    # secondary INTRINSIC Lya flux (foreground H I divided out) with a logN-propagated error. keeps `flux` primary;
    # this is a different, model-dependent quantity (A7) and swings ~0.9-1.6 over logN +-0.3. only when N(HI) exists.
    hi = nhiof(sn)
    if hi is None:
        return {}
    logN, logN_err, src = hi
    base = _lya_ismcorr_flux(w, bf, logN, vline=vline)
    if base is None:
        return {}
    F, hw = base
    fp = _lya_ismcorr_flux(w, bf, logN + logN_err, vline=vline)
    fm = _lya_ismcorr_flux(w, bf, logN - logN_err, vline=vline)
    err = abs(fp[0] - fm[0]) / 2.0 * 1e15 if (fp and fm) else None
    return {"F_ismcorr": round(F * 1e15, 1),
            "F_ismcorr_err": round(err, 1) if err is not None else None,
            "lya_core_halfwidth_kms": round(hw),
            "logN_HI_used": logN, "logN_HI_source": src}


def _fit_profiles(bf, instr=None):
    # multi-model fit + P-Cygni detection on the cont-sub emission profile (exclude the central notch).
    # returns {models:{name:{params,bic,aicc,redchi2}}, best_model, pcygni, _fits} or None if no emission.
    fitm = bf["emis"] & ~bf["notch"]
    vv, ff = bf["v"][fitm], bf["fc"][fitm]
    if len(vv) < 20 or np.nanmax(ff) <= 0:      # absorption / no emission at this epoch
        return None
    if np.nanmin(ff) / np.nanmax(ff) < -0.3:    # deep central absorption vs peak: photospheric P-Cygni OR low-SNR
        coh = _profile_coherence(vv, ff)        # honest confidence: coherent structure vs noise (amp/sigma is backwards)
        reason = "photospheric" if coh >= 12 else ("marginal" if coh >= 6 else "low_snr")
        return {"models": None, "best_model": None, "pcygni": True, "pcygni_reason": reason,
                "coherence": round(coh, 1), "_fits": {}}
    sigma = float(np.nanstd(bf["fc"][bf["contfit"]])) or (0.05 * np.nanmax(ff))
    wfloor = _width_floor(instr)
    fits = _fit_models(vv, ff, wfloor)
    coh = round(_profile_coherence(vv, ff), 1)      # honest confidence for EVERY epoch (amp/scatter), for the detection gate
    if not fits:
        return {"models": None, "best_model": None, "pcygni": False, "coherence": coh, "_fits": {}}
    models = {}
    for name, (fn, p) in fits.items():
        ic = _model_ic(vv, ff, fn, p, sigma)
        ic["pegged"] = _pegged(name, p, wfloor)                  # shape params on a fit bound (empty = well inside)
        if name == "kwok":
            ic["inner_v"] = round(float(p[1] + p[3] - p[4]))    # mu+vc-vin = shell blue edge
        models[name] = ic
    # generic peg guard. a model whose shape parameters peg a bound is
    # barred from SELECTION (kept in models{} with its `pegged` list for the record). census before the change:
    # 8/55 shipped best models were pegged (1993J Lya mu@-9000, 2010jl Lya sig@floor on G130M, kwok vc@bound ...).
    # then the parsimony rule: a k>3 model (skew k=4, kwok k=5) only wins if it beats the best unpegged k=3 model
    # by >= DBIC_STRONG (kass & raftery 'strong'); a thinner margin is window-flippable (SN2005ip d3065 kwok won by
    # dBIC 7.5 and flipped at a wider window). if every model pegs, best_model=None and shape_note says so --
    # the flux is model-independent and unaffected.
    unpeg = {n: m for n, m in models.items() if not m["pegged"]}
    note = None
    if not unpeg:
        best = None; note = "all shape models peg a fit bound: shape unconstrained (flux unaffected)"
    else:
        best = min(unpeg, key=lambda n: unpeg[n]["bic"])
        simple = [n for n in unpeg if unpeg[n]["k"] <= 3]
        if unpeg[best]["k"] > 3 and simple:
            best3 = min(simple, key=lambda n: unpeg[n]["bic"])
            if unpeg[best3]["bic"] - unpeg[best]["bic"] < DBIC_STRONG:
                note = f"{best} beat {best3} by dBIC<{DBIC_STRONG:.0f}: kept the simpler model"; best = best3
        elif unpeg[best]["k"] > 3 and not simple and coh < COH_MIN:
            # no unpegged k<=3 competitor (the narrow models all pegged) AND the profile is not coherent: a k>3
            # shape would be selected BY DEFAULT on a noise feature (AT2022acko Lya d18). do not label it.
            note = f"only a k>3 model ({best}) is unpegged and coherence {coh}<{COH_MIN:.0f}: shape unconstrained (flux unaffected)"
            best = None
        if best is not None:
            barred = [n for n in models if n not in unpeg]
            if barred and models[min(models, key=lambda n: models[n]["bic"])]["pegged"]:
                note = (note + "; " if note else "") + f"lowest-BIC model {min(models, key=lambda n: models[n]['bic'])} barred (pegged)"
    return {"models": models, "best_model": best, "shape_note": note,
            "pcygni": False, "coherence": coh, "_fits": fits}


def _edge_flux_frac(bf, ew=2000.0):
    # emission quality: a real line tapers to ~0 at the window edges. edge/peak ~0.2-0.6 = genuinely broad
    # CSM shell; >0.7 = continuum-subtraction artifact or photospheric P-Cygni (the deg-1 continuum slope
    # leaks into the window). computed from the same bf, no re-fit.
    v, fc = bf["v"], bf["fc"]
    lo, hi = v[bf["emis"]].min(), v[bf["emis"]].max()
    core = (v > lo + ew) & (v < hi - ew)
    if not core.any():
        return np.nan
    peak = np.nanmax(fc[core])
    if not np.isfinite(peak) or peak <= 0:
        return np.nan
    blue = (v >= lo) & (v < lo + ew); red = (v > hi - ew) & (v <= hi)
    eB = np.nanmedian(fc[blue]) if blue.any() else 0.0
    eR = np.nanmedian(fc[red]) if red.any() else 0.0
    return float(max(abs(eB), abs(eR)) / peak)


def _is_spike(v, fc, emis):
    # unresolved single/double-pixel cosmic-ray / hot-pixel spike, NOT a resolved emission line. a real UV
    # line spans many pixels above half-max; a CR is 1-2 px with the neighbors near zero. (MAMA is photon-
    # counting with no CR rejection, so these survive to the 1D.) reject as a false detection.
    ff = fc[emis]
    if len(ff) < 5:
        return False
    ip = int(np.argmax(ff)); pk = ff[ip]
    if pk <= 0:
        return False
    n = 1; j = ip - 1
    while j >= 0 and ff[j] > 0.5 * pk:
        n += 1; j -= 1
    j = ip + 1
    while j < len(ff) and ff[j] > 0.5 * pk:
        n += 1; j += 1
    far = np.concatenate([ff[max(0, ip - 6):max(0, ip - 2)], ff[ip + 3:ip + 7]])   # flux a few px off the peak
    isolated = far.size == 0 or np.nanmax(far) < 0.3 * pk
    return n <= 2 and isolated


def _flux_syst(w, f, lam0, vline, vcont, deg, base_F, z=0.0):
    # bostroem+2026 / sembach&savage systematic flux error: vary the integration limits, the continuum
    # window, and smoothing; take the LARGEST deviation from nominal. this is the DOMINANT uncertainty
    # (continuum/window placement); a photon-noise MC understates it ~10-25x. in 1e-15 units.
    dv = 1500.0
    devs = []
    for vl in ((vline[0] + dv, vline[1]), (vline[0] - dv, vline[1]), (vline[0], vline[1] + dv), (vline[0], vline[1] - dv)):
        Fv = bostroem_flux(w, f, lam0=lam0, vline=vl, vcont=vcont, deg=deg, z=z)["F"]
        if np.isfinite(Fv):
            devs.append(abs(Fv - base_F))
    for dc in (dv, -dv):
        Fv = bostroem_flux(w, f, lam0=lam0, vline=vline, vcont=(vcont[0] - dc, vcont[1] + dc), deg=deg, z=z)["F"]
        if np.isfinite(Fv):
            devs.append(abs(Fv - base_F))
    Fs = bostroem_flux(w, f, lam0=lam0, vline=vline, vcont=vcont, deg=deg, smooth=5, z=z)["F"]
    if np.isfinite(Fs):
        devs.append(abs(Fs - base_F))
    return max(devs) * 1e15 if devs else np.nan


def _emis_rec(bf, ph, instr, prof=None, vphot=None, syst_err=None):
    # two-layer detection gate. the marginal real-vs-noise boundary is genuinely fuzzy
    # (coherence / syst-SNR overlap between faint real lines and noise), so we do NOT force a perfect binary.
    # layer 1: photon-significance floor F > 3*sigma_photon -- rejects low-count junk (few counts -> large error).
    # layer 2: drop only when BOTH metrics agree it is noise: shape coherence low_snr (<6) AND syst-SNR F/flux_err
    # is <3. everything surviving keeps its coherence grade + flux_reliable flag.
    F, sigF = bf["F"] * 1e15, bf["sigF"] * 1e15
    has_err = np.isfinite(sigF) and sigF > 0
    coh = prof.get("coherence") if prof else None
    ferr = syst_err if (syst_err is not None and np.isfinite(syst_err)) else (sigF if has_err else None)
    ssnr = (F / ferr) if (ferr and ferr > 0) else np.inf
    if not (F > 3 * sigF if has_err else F > 0):
        return None                                 # photon-significance floor: rejects low-count junk
    if (coh is None or coh < 6) and ssnr < 3:
        return None                                 # + drop unambiguous noise: both the shape AND the systematic fail
    if _is_spike(bf["v"], bf["fc"], bf["emis"]):
        return None                                 # unresolved cosmic-ray/hot-pixel spike, not a real line
    edge = _edge_flux_frac(bf)
    pcyg = bool(prof["pcygni"]) if prof else False
    reason = prof.get("pcygni_reason") if prof else None
    vph = vphot if reason == "photospheric" else None   # only carry v_phot for a CONFIDENT photospheric feature
    # flux_reliable also requires the SYSTEMATIC significance F/flux_err >= 3. the layer-1 gate is photon
    # SNR (which understates the real error 10-25x), so an epoch could pass detection with a flux smaller than its own
    # continuum-placement error and still be flagged reliable. nothing is dropped; the flag becomes honest.
    unreliable = []
    if pcyg:
        unreliable.append("pcygni: deep central absorption, continuum-placement dominated")
    if ssnr < 3:
        unreliable.append(f"low syst-SNR: F/flux_err={ssnr:.1f}<3")
    rec = {"phase": ph, "flux": round(F, 1), "flux_err": round(ferr, 1) if ferr is not None else None,
           "flux_err_photon": round(sigF, 1) if has_err else None, "instr": instr,
           "flux_reliable": not unreliable,
           "flux_reliable_reason": "; ".join(unreliable) or None,
           "edge_flux_frac": round(edge, 2) if np.isfinite(edge) else None,
           "pcygni": pcyg, "pcygni_reason": reason,   # photospheric / marginal / low_snr (was an overclaiming bool)
           "coherence": coh,                          # smoothed-amp / residual scatter: the honest confidence number
           "v_phot_kms": vph,                         # only for a confident photospheric epoch (approx, normalized abs min)
           "best_model": prof["best_model"] if prof else None,
           "shape_note": prof.get("shape_note") if prof else None,   # peg-guard / parsimony decisions, if any
           "models": prof["models"] if prof else None}
    return rec


def _emis_diag(bf, prof, sn, instr, ph, line, plotdir):
    # per-epoch scrutiny panel so a follower can eyeball the steps behind every reported flux:
    # (top) dereddened flux + the fitted continuum + the flank pts it used + line window + masked core;
    # (bottom) continuum-subtracted profile + the kwok shell overlay where the fit is well-constrained.
    v, ff, cont, fc = bf["v"], bf["ff"], bf["cont"], bf["fc"]
    emis, contfit, notch = bf["emis"], bf["contfit"], bf["notch"]
    vlo, vhi = v[emis].min(), v[emis].max()
    band = contfit | emis                        # top panel = continuum+line region only (skip the far-blue CCD spike)
    vsmin, vsmax = v[band].min(), v[band].max()
    show = (v >= vsmin) & (v <= vsmax)
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(7.2, 6.4))
    a0.plot(v[show], ff[show] * 1e15, color="0.45", lw=0.8, label="dereddened flux")
    a0.plot(v[show], cont[show] * 1e15, color="crimson", lw=1.3, label=f"deg-{bf['deg']} continuum")
    a0.plot(v[contfit], ff[contfit] * 1e15, ".", color="crimson", ms=2.5, label="continuum-fit flanks")
    a0.axvspan(vlo, vhi, color="gold", alpha=0.10, label="line window")
    if notch.any():
        a0.axvspan(v[notch].min(), v[notch].max(), color="dodgerblue", alpha=0.14, label="masked core (interp)")
    a0.set_xlim(vsmin, vsmax)
    a0.set_ylabel("flux [1e-15]"); a0.legend(fontsize=7, loc="upper right")
    a0.set_title(f"{sn}  {line}  {instr}  day{ph:.0f}", fontsize=9)
    a1.axhline(0, color="0.7", lw=0.7)
    a1.plot(v[emis], fc[emis] * 1e15, color="navy", lw=1.0, label="continuum-subtracted")
    _mcol = {"gaussian": "tab:green", "lorentzian": "tab:purple", "skew": "tab:orange", "kwok": "crimson"}
    if prof and prof.get("_fits"):
        vv = np.linspace(vlo, vhi, 500)
        for name, (fn, pp) in prof["_fits"].items():
            best = name == prof["best_model"]
            a1.plot(vv, fn(vv, *pp) * 1e15, color=_mcol.get(name, "gray"),
                    lw=2.2 if best else 0.9, label=name + (" *" if best else ""))
        km = prof["models"].get("kwok", {}) if prof["models"] else {}
        note = f"best={prof['best_model']}" + (f"\nkwok inner_v={km['inner_v']}" if "inner_v" in km else "")
        a1.text(0.02, 0.96, note, transform=a1.transAxes, fontsize=7, va="top")
    elif prof and prof.get("pcygni"):
        rn = prof.get("pcygni_reason", "photospheric"); ch = prof.get("coherence")
        a1.text(0.02, 0.96, f"deep central absorption\n({rn}, coherence={ch})\nnot a reliable emission line",
                transform=a1.transAxes, fontsize=7, va="top", color="crimson")
    if notch.any():
        a1.axvspan(v[notch].min(), v[notch].max(), color="dodgerblue", alpha=0.14)
    F, sigF = bf["F"] * 1e15, bf["sigF"] * 1e15
    lab = f"F={F:.0f}" + (f" +-{sigF:.0f}" if np.isfinite(sigF) else "") + " e-15"
    a1.text(0.02, 0.80, lab, transform=a1.transAxes, fontsize=7, va="top", color="navy")
    a1.set_xlabel("velocity [km/s]"); a1.set_ylabel("flux [1e-15]"); a1.legend(fontsize=7, loc="upper right")
    a1.set_xlim(vlo, vhi)
    os.makedirs(plotdir, exist_ok=True)
    tag = line.replace(" ", "").lower()
    instr_safe = instr.replace("?", "unk")      # '?' is invalid in Windows filenames
    fig.tight_layout(); _savefig(fig, os.path.join(plotdir, f"{sn}_{instr_safe}_day{ph:.0f}_{tag}_diag.png"))
    plt.close(fig)


def _emis_summary(sn, mg2, lya, mg2_prof, lya_prof, plotdir):
    # per-SN roll-up: (left) line flux vs phase for Mg II + Lya; (right) peak-normalized profiles (B10)
    # so the shape evolution is visible independent of the fading flux.
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(11, 4.2))
    for recs, col, lab in ((mg2, "darkorange", "Mg II"), (lya, "steelblue", "Lya")):
        if not recs:
            continue
        rel = [r for r in recs if r.get("flux_reliable", True)]     # clean emission line
        unr = [r for r in recs if not r.get("flux_reliable", True)]  # photospheric P-Cygni, flux not trustworthy
        if rel:
            a0.errorbar([r["phase"] for r in rel], [r["flux"] for r in rel], yerr=[r["flux_err"] or 0 for r in rel],
                        fmt="o-", color=col, ms=4, lw=1, capsize=2, label=lab)
        if unr:
            a0.plot([r["phase"] for r in unr], [r["flux"] for r in unr], "x", color=col, ms=6, alpha=0.55,
                    label=f"{lab} flux unreliable (pcygni / low syst-SNR)")
    a0.set_xlabel("phase [day]"); a0.set_ylabel("line flux [1e-15]"); a0.set_yscale("log")
    a0.set_title(f"{sn}  emission-line flux vs phase", fontsize=9); a0.legend(fontsize=8)
    allp = mg2_prof + lya_prof
    if allp:
        phs = [pp[0] for pp in allp]; pmin, pmax = min(phs), max(phs)
        norm = matplotlib.colors.Normalize(vmin=pmin, vmax=pmax)
        for prof, ls in ((mg2_prof, "-"), (lya_prof, "--")):
            for ph, vv, fcv in prof:
                pk = np.nanmax(fcv)
                if pk <= 0:
                    continue
                a1.plot(vv, fcv / pk, color=plt.cm.viridis(norm(ph)), lw=0.9, ls=ls)
        a1.axvline(0, color="0.7", lw=0.7); a1.set_ylim(-0.3, 1.15)
        a1.set_xlabel("velocity [km/s]"); a1.set_ylabel("peak-normalized flux")
        a1.set_title("profiles peak=1 (Mg II solid, Lya dashed)", fontsize=9)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm); sm.set_array([])
        fig.colorbar(sm, ax=a1, label="phase [day]")   # was 'color=phase' with no key; give the reader the mapping
    os.makedirs(plotdir, exist_ok=True)
    fig.tight_layout(); _savefig(fig, os.path.join(plotdir, f"{sn}_emission_summary.png"))
    plt.close(fig)


def compute_emission(sn):
    z, mw, host = zof(sn), ebvof(sn), hostof(sn)
    mg2, lya, seen = [], [], set()
    mg2_prof, lya_prof = [], []                 # (phase, v, fc) over the line window, for the peak-norm montage
    plotdir = os.path.join(OUT, sn, "emission")
    for old in glob.glob(os.path.join(plotdir, "*_diag.png")):   # drop stale per-epoch diags so non-detections leave no orphan PNG
        os.remove(old)
    for pf in sorted(glob.glob(f"{OUT}/{sn}/**/*_native.txt", recursive=True), key=phase_of):
        pth = pf.replace("\\", "/")
        ph = phase_of(pth)
        if "epochcoadd" in pth or "/epochs/" in pth or not np.isfinite(ph):
            continue                            # skip stitched epoch/coadd files -> use per-grating native
        instr = next((g for g in EMIS_GRATINGS + ECHELLE_GRATINGS if f"/{g}/" in pth), "?")
        if instr in ECHELLE_GRATINGS:
            continue                            # echelle is an ISM/absorption product, not an emission one (see provenance)
        w, f, e = load_spec(pf, z)
        fac = deredden(w, z, mw, host); f = f * fac; e = e * fac       # deredden the errors too
        ok = np.isfinite(f); w, f, e = w[ok], f[ok], e[ok]
        if len(w) < 50:
            continue
        if w.min() < 2680 and w.max() > 2900 and ("mg", round(ph)) not in seen:
            seen.add(("mg", round(ph)))
            win = MG2_WINDOW.get(sn.upper(), _MG2_DEFAULT)
            bf = bostroem_flux(w, f, e=e, lam0=MG2, vline=win, deg=1, z=z)
            # item-1: fit the SHAPE over the full default window; a narrow per-SN flux override (GGI (-6000,4000))
            # truncates a broad shell into a lorentzian. the flux stays on `win`; only the shape label uses the wide.
            bf_shape = bf if win == _MG2_DEFAULT else bostroem_flux(w, f, e=e, lam0=MG2, vline=_MG2_DEFAULT, deg=1, z=z)
            prof = _fit_profiles(bf_shape, instr)
            vphot = _vphot_normalized(w, f, MG2) if (prof and prof.get("pcygni")) else None
            syst = _flux_syst(w, f, MG2, win, (-16000, 16000), 1, bf["F"], z=z)
            rec = _emis_rec(bf, ph, instr, prof, vphot, syst)
            if rec:
                rec.update(_diag_keys(bf, w, f, MG2, win, 1, z, (-14000, 9000)))
                mg2.append(rec)
                _emis_diag(bf_shape, prof, sn, instr, ph, "Mg II", plotdir)
                mg2_prof.append((ph, bf["v"][bf["emis"]], bf["fc"][bf["emis"]]))
        if w.min() < 1185 and w.max() > 1270 and ("ly", round(ph)) not in seen:
            seen.add(("ly", round(ph)))
            fly = _airglow_subtract(w, f, z, _airglow_bg(pth, w, z))    # remove geocoronal Lya via the x1d background
            bf = bostroem_flux(w, fly, e=e, lam0=LYA, vline=_LYA_DEFAULT, deg=0, vabs=800, vcont=(-13000, 13000), z=z)
            prof = _fit_profiles(bf, instr)       # same 4-model selection on Lya (kwok validated there too)
            vphot = _vphot_normalized(w, fly, LYA) if (prof and prof.get("pcygni")) else None
            syst = _flux_syst(w, fly, LYA, _LYA_DEFAULT, (-13000, 13000), 0, bf["F"])
            rec = _emis_rec(bf, ph, instr, prof, vphot, syst)
            if rec:
                rec.update(_diag_keys(bf, w, fly, LYA, _LYA_DEFAULT, 0, z, (-13000, 8000)))
                if not rec["pcygni"]:                          # F_ismcorr is meaningless on a core-peaked P-Cygni profile (the bridge cuts the peak)
                    rec.update(_ismcorr_keys(w, bf, sn))
                lya.append(rec)
                _emis_diag(bf, prof, sn, instr, ph, "Lya", plotdir)
                lya_prof.append((ph, bf["v"][bf["emis"]], bf["fc"][bf["emis"]]))
    if mg2 or lya:
        _emis_summary(sn, mg2, lya, mg2_prof, lya_prof, plotdir)
    return mg2, lya


# --- assemble + write ------------------------------------------------------------------------------
def build_emission(sn):
    mg2, lya = compute_emission(sn)
    peak = max((r["flux"] for r in mg2 + lya), default=0.0)
    flags = []
    if peak > 1e4:      # >1e-11; a normal extragalactic UV SN line is ~1-1e3 e-15. flags nearby/resolved objects
        flags.append("flux_scale_outlier: nearby/resolved object (aperture-dependent flux, not comparable to extragalactic SNe)")
    prod = {
        "sn": sn.upper(),
        "sn_type": tnstype(sn),
        "generated": datetime.date.today().isoformat(),
        "flags": flags,
        "provenance": {
            "z": zof(sn), "ebv_mw": ebvof(sn), "host_ebv": hostof(sn), "host_ebv_src": hostsrc(sn),
            "dered": "F19 Rv=3.1, MW at observed wvl + host at rest",
            "flux_units": "1e-15 erg s-1 cm-2",
            "emission_method": "continuum-subtracted emission-line flux (specutils line_flux over the line window)",
            "emission_content": "Mg II 2800 and Ly-a 1216 fluxes are reported for detected epochs. Broad late-time shell emission and near-maximum photospheric P-Cygni peaks can both appear here, so interpret each epoch with sn_type, phase, best_model, and pcygni.",
            "emission_epochs": "Epochs pass a two-step gate: F > 3 sigma_photon first, then only clearly noisy cases are dropped when coherence < 6 and F/flux_err < 3. flux_err is the continuum/window systematic from varying the measurement setup; flux_err_photon is the pixel-noise term kept for reference.",
            "emission_quality": "flux_reliable=false marks epochs where the quoted line flux is not cleanly separated from continuum placement, either because a deep P-Cygni trough drives the measurement or because F/flux_err < 3. flux_reliable_reason records which case. pcygni_reason grades the profile as photospheric, marginal, or low_snr from the reported coherence value. v_phot_kms is reported only for confident photospheric cases, and edge_flux_frac is kept as a window-edge diagnostic.",
            "csm_benchmark": "For SN2023ixf and SN2024ggi, the Mg II comparison to Bostroem+2026 Table 4 is an observed-flux comparison. The stored `flux` values here are MW+host dereddened, so they sit above Table 4 by the extinction factor. For Ly-a, the remaining offset is a core-treatment difference: that comparison bridges the airglow/ISM core, while this product subtracts airglow and keeps the damped foreground trough.",
            "lya_airglow": "geocoronal Lya airglow (observed 1215.67 A = v approx -cz in the SN frame) is removed at product-build time: for STIS low-resolution the sibling x1d BACKGROUND array localizes it (redshift-independent) and a local poly2 + alpha*background decomposition subtracts only the airglow while preserving the real SN line; COS/flat-background epochs fall back to a narrow velocity notch. the damped foreground ISM Lya absorption trough (0..+1000 km/s for 2023ixf, N(HI)~1e21) is NOT filled: `flux` is the Lya that reaches us after H I absorption, minus airglow.",
            "lya_window": "the Ly-a integration window is (-10000, +5000) km/s. The wider blue edge follows the same shell wing seen in Mg II and does not pick up a separate C III 1176 feature; the red edge stays at +5000 to avoid N V 1238/1242. Narrow host-ISM absorption can still carve structure into the broad line, so the direct-integral flux is measured as seen while the shape models remain approximate over those notches.",
            "lya_ismcorr": "F_ismcorr is a SECONDARY, model-dependent Lya flux: the INTRINSIC SN Lya before the foreground H I absorbed it, obtained by dividing the airglow-subtracted continuum-subtracted profile by the ISM transmission exp(-tau) from the measured N(HI) (b=25 km/s, host frame) where transmission > 0.3, bridging the opaque core. `flux` stays primary (what reaches us). F_ismcorr swings ~0.9-1.6x over logN(HI) +-0.3 because the damping wings multiply the whole blue peak, so F_ismcorr_err is the logN-propagated half-range and logN_HI_source records whether N(HI) is curated (ism_columns.csv) or the automated damped-Lya fit (lya_nhi_summary.csv). lya_core_halfwidth_kms is the bridged (T<=0.3) core half-width. present only when an N(HI) exists for the SN, and only for non-P-Cygni epochs (the bridge cuts the peak of a core-peaked profile).",
            "window_diagnostics": "v_blue99/v_red99 are the velocities enclosing 99% of the smoothed continuum-subtracted flux, reported only as a comparison metric. spill_frac compares the default window to one widened by 3000 km/s on each side; values near 1 mean the default window already captures the line.",
            "gratings": "measured on low/medium-resolution first-order gratings only (STIS G230L/G230LB/G140L/G230M/G230MB, COS G130M/G160M/G185M/G225M/G285M). Echelle modes (E230M/E140M/E230H/E140H) are excluded because narrow orders do not define a stable continuum over a +-10000 km/s window, and in cases like SN1998S they isolate narrow ISM absorption rather than the broad emission component.",
            "shell_model": "clean-emission epochs are fit with four models (gaussian, lorentzian, skew-normal, kwok off-center-hole shell), each on O(1)-normalized flux. models{} stores each fit's params, BIC, AICc, redchi2, and any `pegged` shape parameters. best_model is chosen from the unpegged fits, with simpler k=3 models preferred unless a higher-parameter model wins by dBIC>=6; if every model pegs, best_model=null and shape_note explains why. The flux itself is always the direct integral, so these models describe shape only.",
        },
        "emission": {"mg2": mg2, "lya": lya},
    }
    dst = os.path.join(OUT, sn, f"{sn}_emission.json")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as fh:
        json.dump(prod, fh, indent=2)
    return prod, dst


def summary_row(prod):
    mg2 = prod["emission"]["mg2"]; lya = prod["emission"]["lya"]
    peak = max((r["flux"] for r in mg2), default=None)
    return {
        "sn": prod["sn"], "sn_type": prod["sn_type"], "z": prod["provenance"]["z"],
        "host_ebv": prod["provenance"]["host_ebv"],
        "n_mg2_epochs": len(mg2), "n_mg2_kwok": sum(1 for r in mg2 if r.get("best_model") == "kwok"),
        "n_mg2_pcygni": sum(1 for r in mg2 if r.get("pcygni")),
        "n_lya_epochs": len(lya), "mg2_peak_e15": peak, "flag": ";".join(prod["flags"]),
    }


def main(names):
    if not names:
        names = sorted(d for d in os.listdir(OUT)
                       if os.path.isdir(os.path.join(OUT, d))
                       and glob.glob(f"{OUT}/{d}/**/*_native.txt", recursive=True))
    summ = []
    for sn in names:
        if sn.upper() not in cat:
            print(f"  skip {sn}: not in catalog"); continue
        prod, dst = build_emission(sn)
        summ.append(summary_row(prod))
        nmg, nly = len(prod["emission"]["mg2"]), len(prod["emission"]["lya"])
        print(f"  {sn:14s} mg2={nmg:2d} lya={nly:2d}  -> {os.path.relpath(dst, ROOT)}")
    if summ:
        with open(SUMMARY, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(summ[0].keys()))
            wr.writeheader(); wr.writerows(summ)
        print(f"\nsummary -> {os.path.relpath(SUMMARY, ROOT)} ({len(summ)} sne)")


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
