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
# usage:  python emission_products.py [SN ...]      (no args -> every SN with spectra under output5)

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
SUMMARY = os.path.join(OUT, "emission_summary.csv")

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
    # kwok+2023/24 asymmetric shell: gaussian-emissivity sphere with an off-center spherical hole cut out.
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
def gaussv(v, A, mu, sig):                        # symmetric peak (k=3)
    return A * np.exp(-0.5 * ((v - mu) / sig) ** 2)


def lorentzian(v, A, mu, gam):                    # narrow-line IIn scattering wings (k=3)
    return A * gam ** 2 / ((v - mu) ** 2 + gam ** 2)


def skewg(v, A, mu, sig, al):                     # mild asymmetry (k=4)
    t = (v - mu) / sig
    return A * np.exp(-0.5 * t ** 2) * (1 + erf(al * t / np.sqrt(2)))


def _fit_models(vv, ff):
    # every fit on O(1)-normalized flux then amplitude scaled back. fitting raw ~1e-15 froze the bounded
    # lorentzian at p0 (curve_fit hits its gradient tol at iter 0); normalizing is the same fix kwok used.
    A = np.nanmax(ff); out = {}
    specs = [("gaussian", gaussv, [1, -3000, 4000], ([0, -9000, 500], [5, 4000, 12000])),
             ("lorentzian", lorentzian, [1, -2000, 3000], ([0, -9000, 300], [5, 4000, 12000])),
             ("skew", skewg, [1, -3000, 4000, -2], ([0, -9000, 500, -20], [5, 4000, 12000, 20])),
             ("kwok", kwok_shell, [1, -2995, 7580, 1750, 5000], (KWOK_LO, KWOK_HI))]
    for name, fn, p0, bnds in specs:
        try:
            p = curve_fit(fn, vv, ff / A, p0=p0, bounds=bnds, maxfev=40000)[0]
            p = np.asarray(p, float); p[0] *= A; out[name] = (fn, p)
        except Exception:
            pass
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
# right call and reproduces bostroem+2026 for 2023ixf; phase-2's "symmetric +-10000" note is superseded.
# MG2_WINDOW is a curated per-SN override for the few broad-blue-wing SNe the default over-captures (GGI ~2x).
# NOTE (aug 17 phase-3.5): a data-driven adaptive-tightening window was tried and REJECTED - it is noise-fragile
# (it cut real emission on 2023ixf d66 429->312 and did not fix GGI), so the fixed default + curated override
# is the robust choice. a fit-based window (integrate over the fitted profile extent) is the future direction.
MG2_WINDOW = {"SN2024GGI": (-6000, 4000)}
_MG2_DEFAULT = (-10000, 6000)


def _fit_profiles(bf):
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
    fits = _fit_models(vv, ff)
    coh = round(_profile_coherence(vv, ff), 1)      # honest confidence for EVERY epoch (amp/scatter), for the detection gate
    if not fits:
        return {"models": None, "best_model": None, "pcygni": False, "coherence": coh, "_fits": {}}
    models = {}
    for name, (fn, p) in fits.items():
        ic = _model_ic(vv, ff, fn, p, sigma)
        if name == "kwok":
            ic["inner_v"] = round(float(p[1] + p[3] - p[4]))    # mu+vc-vin = shell blue edge
        models[name] = ic
    # item-2 guard: don't let kwok WIN on a pegged FWHM floor (KWOK_LO[2]=3000) or a razor-thin, window-flippable
    # BIC margin over lorentzian (SN2005ip d3065 pegged at 3000 won by dBIC 7.5, flips at a wider window). keep
    # kwok in models{} for the record; only bar it from the SELECTION. real broad shells (2023ixf, GGI, 1998S,
    # 2026ayt) have FWHM >> the floor and beat lorentzian by a wide margin, so they keep kwok.
    sel = dict(models)
    if "kwok" in sel and "lorentzian" in sel:
        pegged = fits["kwok"][1][2] <= KWOK_LO[2] * 1.02
        thin = (sel["lorentzian"]["bic"] - sel["kwok"]["bic"]) < 6.0
        if pegged or thin:
            sel = {k: v for k, v in sel.items() if k != "kwok"}
    return {"models": models, "best_model": min(sel, key=lambda n: sel[n]["bic"]),
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
    # (continuum/window placement); a photon-noise MC understates it ~10-25x (phase-2 B1/B7). in 1e-15 units.
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
    # HONEST detection gate (phase-4 item 3). the marginal real-vs-noise boundary is genuinely fuzzy
    # (coherence / syst-SNR overlap between faint real lines and noise), so we do NOT force a perfect binary.
    # two layers: (1) KEEP the phase-3.5 photon-significance floor F > 3*sigma_photon -- it legitimately rejects
    # low-count junk (few counts -> large photon error); (2) additionally drop an epoch only when BOTH honest
    # metrics agree it is noise: the shape coherence is low_snr (<6) AND the systematic significance F/flux_err
    # is <3. everything surviving keeps its coherence grade + flux_reliable flag. see docs/analysis_phase4.md sec 2.
    F, sigF = bf["F"] * 1e15, bf["sigF"] * 1e15
    has_err = np.isfinite(sigF) and sigF > 0
    coh = prof.get("coherence") if prof else None
    ferr = syst_err if (syst_err is not None and np.isfinite(syst_err)) else (sigF if has_err else None)
    ssnr = (F / ferr) if (ferr and ferr > 0) else np.inf
    if not (F > 3 * sigF if has_err else F > 0):
        return None                                 # photon-significance floor (phase 3.5): rejects low-count junk
    if (coh is None or coh < 6) and ssnr < 3:
        return None                                 # + drop unambiguous noise: both the shape AND the systematic fail
    if _is_spike(bf["v"], bf["fc"], bf["emis"]):
        return None                                 # unresolved cosmic-ray/hot-pixel spike, not a real line
    edge = _edge_flux_frac(bf)
    pcyg = bool(prof["pcygni"]) if prof else False
    reason = prof.get("pcygni_reason") if prof else None
    vph = vphot if reason == "photospheric" else None   # only carry v_phot for a CONFIDENT photospheric feature
    rec = {"phase": ph, "flux": round(F, 1), "flux_err": round(ferr, 1) if ferr is not None else None,
           "flux_err_photon": round(sigF, 1) if has_err else None, "instr": instr,
           "flux_reliable": not pcyg,                 # deep-central-absorption flux is continuum-placement dominated
           "edge_flux_frac": round(edge, 2) if np.isfinite(edge) else None,
           "pcygni": pcyg, "pcygni_reason": reason,   # photospheric / marginal / low_snr (was an overclaiming bool)
           "coherence": coh,                          # smoothed-amp / residual scatter: the honest confidence number
           "v_phot_kms": vph,                         # only for a confident photospheric epoch (approx, normalized abs min)
           "best_model": prof["best_model"] if prof else None,
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
                    label=f"{lab} pcygni (unreliable)")
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
        instr = next((g for g in ("G230LB", "G230L", "G140L", "E230M", "G230M") if f"/{g}/" in pth), "?")
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
            prof = _fit_profiles(bf_shape)
            vphot = _vphot_normalized(w, f, MG2) if (prof and prof.get("pcygni")) else None
            syst = _flux_syst(w, f, MG2, win, (-16000, 16000), 1, bf["F"], z=z)
            rec = _emis_rec(bf, ph, instr, prof, vphot, syst)
            if rec:
                mg2.append(rec)
                _emis_diag(bf_shape, prof, sn, instr, ph, "Mg II", plotdir)
                mg2_prof.append((ph, bf["v"][bf["emis"]], bf["fc"][bf["emis"]]))
        if w.min() < 1185 and w.max() > 1270 and ("ly", round(ph)) not in seen:
            seen.add(("ly", round(ph)))
            fly = _airglow_subtract(w, f, z, _airglow_bg(pth, w, z))    # remove geocoronal Lya via the x1d background
            bf = bostroem_flux(w, fly, e=e, lam0=LYA, vline=(-8000, 5000), deg=0, vabs=800, vcont=(-13000, 13000), z=z)
            prof = _fit_profiles(bf)       # same 4-model selection on Lya (kwok validated there too)
            vphot = _vphot_normalized(w, fly, LYA) if (prof and prof.get("pcygni")) else None
            syst = _flux_syst(w, fly, LYA, (-8000, 5000), (-13000, 13000), 0, bf["F"])
            rec = _emis_rec(bf, ph, instr, prof, vphot, syst)
            if rec:
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
            "emission_content": "Mg II 2800 + Lya 1216 emission-line flux for any epoch where the line is detected; NOT restricted to CSM - includes CSM-interaction shells (broad, strengthen at late phase, best_model kwok/skew) AND photospheric P-Cygni peaks (near max, pcygni=true). use sn_type + phase + best_model + pcygni to interpret.",
            "emission_epochs": "epochs kept by a two-layer detection gate: (1) the photon-significance floor F > 3 sigma_photon; (2) additionally dropped ONLY when BOTH the shape coherence is low_snr (<6) AND the systematic significance F/flux_err is <3 -- i.e. only unambiguous noise is removed. marginal gray-zone epochs are KEPT carrying their coherence grade + flux_reliable flag for downstream filtering, because the real-vs-noise boundary is genuinely fuzzy (we grade rather than force-classify). flux_err is the SYSTEMATIC error (bostroem+2026 / sembach&savage 1992 method: vary the integration limits, the continuum window, and smoothing, take the LARGEST deviation) - the dominant continuum-placement uncertainty. flux_err_photon is the pixel-noise MC, which understates the real error ~10-25x and is kept only for reference.",
            "emission_quality": "flux_reliable=false marks an epoch with a deep central absorption vs the peak (pcygni=true, min/max<-0.3): its continuum-subtracted flux is dominated by continuum placement over the SN photosphere, so it is kept but is NOT a reliable emission-line flux. those epochs get no clean shape model. pcygni_reason grades what it actually is, from the profile COHERENCE (smoothed-amplitude / pixel residual scatter, reported as `coherence`): 'photospheric' (coherence>=12, a genuine coherent P-Cygni), 'marginal' (6-12, a low-SNR gray-zone feature), 'low_snr' (<6, noise / no clean line). this replaces an overclaiming pcygni bool that labeled noise as photospheric. v_phot_kms carries an APPROXIMATE photospheric velocity (absorption minimum of the continuum-normalized spectrum, UV is a forest so rough) ONLY for a confident 'photospheric' epoch. edge_flux_frac (|cont-sub flux at window edges|/peak) is a reported diagnostic. judge each epoch with pcygni_reason + coherence + best_model + sn_type + phase.",
            "csm_benchmark": "our MW+host-dereddened Mg II 2800 and Lya 1216 fluxes sit ~15-25% above bostroem+2026 table 4 for SN2023IXF. this is NOT an extinction-convention difference: bostroem applies the SAME MW+host extinction (her MW 0.0076, host 0.031, Fitzpatrick) as we do (MW 0.0089, host 0.031, F19), so the offset is a FLUX-MEASUREMENT difference (continuum placement + fixed integration window vs her exact per-epoch recipe), within our own systematic error (flux_err). note: applying MW-only extinction appears to match by coincidence -- our flux runs ~15% high before the host term, not because the conventions agree. our custom reduction reproduces the MAST default x1d to 1-3%, so it is not a flux-cal error. SN2024GGI reads ~3x high from a genuine adopted-extinction disagreement (our total E(B-V)=0.154 from Jacobson-Galan+2024a/Chen+2024 vs bostroem's ~0.046), not a method error.",
            "lya_airglow": "geocoronal Lya airglow (observed 1215.67 A = v approx -cz in the SN frame) is removed at product-build time: for STIS low-resolution the sibling x1d BACKGROUND array localizes it (redshift-independent) and a local poly2 + alpha*background decomposition subtracts only the airglow while preserving the real SN line; echelle/COS/flat-background epochs fall back to a narrow velocity notch.",
            "shell_model": "clean-emission epochs are fit with 4 models (gaussian, lorentzian, skew-normal, kwok off-center-hole shell), each fit on O(1)-normalized flux (fixes a prior freeze that pinned the lorentzian at its guess); models{} lists every fit's params + BIC/AICc/redchi2 and best_model is the BIC winner. kwok carries inner_v = mu+vc-vin (shell blue edge). the flux is model-independent (direct integral); the models only describe the shape. best model varies by regime: broad CSM shell -> kwok/skew, narrow-line IIn -> gaussian or lorentzian (electron-scattering wings). photospheric P-Cygni epochs get no shape model (see emission_quality).",
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
