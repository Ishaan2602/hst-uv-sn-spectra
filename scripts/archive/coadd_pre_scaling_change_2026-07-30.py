import warnings
import numpy as np

# shared coadd + scaling for the stis/cos product building.
# common wavelength axis: finer in the uv, coarser to the red. extended blue for the fuv modes.
COMMON_AXIS = np.concatenate([
    np.arange(1100, 1650, 1.0),
    np.arange(1650, 3050, 1.4),
    np.arange(3050, 5600, 2.7),
    np.arange(5600, 10260, 4.9)])
# drop the occasional float-precision duplicate at a segment join so the axis is STRICTLY increasing
# (the FluxConservingResampler requires it - this is the splice-edge issue that bit us in phase 1).
COMMON_AXIS = COMMON_AXIS[np.concatenate([[True], np.diff(COMMON_AXIS) > 1e-6])]

# inter-grating overlap windows (A) used to align each grating to the g230lb anchor.
OVERLAP = {'G430L': (2900, 3150), 'G430M': (2900, 3150),
           'G750L': (5300, 5750), 'G750M': (5300, 5750)}


def clean_wf(w, f, e=None, dq=None):
    # dq: drop ONLY genuinely-bad px (512 bad-ref, 256 satur, 4 bad-det). KEEP bit 16 "high dark":
    # it flags ~half the ccd in broad coherent swaths and on a bright source the value is fine -
    # dropping it gutted real structure (the Mg II 2800 two-teeth). finite + nonzero; sort + dedup.
    w = np.asarray(w, float); f = np.asarray(f, float)
    m = np.isfinite(w) & np.isfinite(f) & (f != 0)
    if dq is not None:
        dq = np.asarray(dq)
        m &= (dq & 512 != 512) & (dq & 256 != 256) & (dq & 4 != 4)
    w, f = w[m], f[m]
    e = np.asarray(e, float)[m] if e is not None else None
    if len(w) == 0:
        return w, f, e
    o = np.argsort(w); w, f = w[o], f[o]
    if e is not None: e = e[o]
    keep = np.concatenate([[True], np.diff(w) > 0])
    return w[keep], f[keep], (e[keep] if e is not None else None)


def resample(w, f, e=None, axis=COMMON_AXIS):
    fo = np.interp(axis, w, f, left=np.nan, right=np.nan)
    eo = np.interp(axis, w, e, left=np.nan, right=np.nan) if e is not None else None
    return fo, eo


def ivar_combine(fluxes, errors):
    # inverse-variance weighted mean on a shared axis. no deviation reject: in the low-snr uv the
    # scatter dwarfs the flux (and the running median sits near zero, so a relative cut divides by
    # ~0 and rejects everything), which is what was gutting the cos/mama uv. dq is already applied
    # upstream in clean_wf; keep every finite point and let the 1/err^2 weighting handle the noise.
    F = np.vstack(fluxes).astype(float)
    E = np.vstack(errors).astype(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        wgt = 1.0 / np.where(E > 0, E ** 2, np.nan)
        num = np.nansum(np.where(np.isfinite(F), F * wgt, 0.0), axis=0)
        den = np.nansum(np.where(np.isfinite(F), wgt, 0.0), axis=0)
        comb = np.where(den > 0, num / den, np.nan)
        cerr = np.where(den > 0, np.sqrt(1.0 / den), np.nan)
    return comb, cerr


def _scale_leg(ref, f):
    # scale leg f onto the running reference by a deg-1 flux-ratio fit in their overlap. the linear
    # trend is applied only across the overlap (held flat outside via the wavelength clip) so it
    # removes the level+slope kink at the join without runaway extrapolation down the leg. falls
    # back to a constant-% median ratio, then to 1.0 (flagged). KEY POINT (phase 7): the deg-1 fit
    # + fallback is the confirmed default; the noted alternatives are constant-% only (Bostroem) and
    # a hard splice at a chosen wavelength per pair.
    ov = np.isfinite(ref) & np.isfinite(f) & (ref > 0) & (f > 0)
    n = int(ov.sum())
    if n >= 20:
        x = COMMON_AXIS[ov]; r = ref[ov] / f[ov]
        med = np.nanmedian(r); sd = np.nanstd(r) or 1.0
        good = np.abs(r - med) < 3 * sd
        if int(good.sum()) >= 10:
            c = np.polyfit(x[good], r[good], 1)
            xcl = np.clip(COMMON_AXIS, x[good].min(), x[good].max())   # hold flat beyond the overlap
            scale = np.clip(np.polyval(c, xcl), 0.3, 3.0)
            return f * scale, {'scale': 'linear', 'slope': float(c[0]), 'intercept': float(c[1])}
    if n >= 5:
        r = np.nanmedian(ref[ov] / f[ov])
        if np.isfinite(r) and 0.2 < r < 5.0:
            return f * float(r), {'scale': 'const', 'factor': float(r)}
    return f, {'scale': 'none', 'factor': 1.0}


def align_and_merge(legs, prefer=()):
    # legs: list of (label, flux_on_COMMON_AXIS). align each leg onto a running reference with
    # _scale_leg (blue/anchor-fit, linear-in-overlap), then nanmedian. anchor = first leg whose
    # label contains a `prefer` substring (e.g. G230LB for stis, STIS CCD for the cross-instrument
    # merge), else the widest-coverage leg. greedy order (most overlap-with-reference first) makes
    # the chain work, e.g. g750l joins through the already-scaled g430l. KEY POINT (phase 7):
    # nanmedian-after-scaling is the confirmed combine; noted alternative is a linear ramp crossfade.
    arrs = [(str(lab), np.asarray(f, float)) for lab, f in legs]
    if not arrs:
        return np.full_like(COMMON_AXIS, np.nan), []
    ai = next((i for i, (lab, _) in enumerate(arrs) if any(p in lab for p in prefer)), None)
    if ai is None:
        ai = int(np.argmax([np.isfinite(f).sum() for _, f in arrs]))
    ref = arrs[ai][1].copy()
    stack = [ref.copy()]
    used = {ai}
    recs = [{'leg': arrs[ai][0], 'scale': 'anchor'}]
    while len(used) < len(arrs):
        cand = [(i, int((np.isfinite(ref) & np.isfinite(f) & (ref > 0) & (f > 0)).sum()))
                for i, (lab, f) in enumerate(arrs) if i not in used]
        i = max(cand, key=lambda t: t[1])[0]
        lab, f = arrs[i]
        sc, rec = _scale_leg(ref, f)
        stack.append(sc)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            ref = np.nanmedian(np.vstack(stack), axis=0)
        used.add(i); recs.append({'leg': lab, **rec})
    return ref, recs


def scale_factor(anchor_w, anchor_f, w, f, overlap):
    # constant-% factor to bring (w,f) onto the anchor: median flux ratio in the overlap window.
    # guard to 1.0 (flagged) if the overlap is faint/edge or the ratio is implausible.
    lo, hi = overlap
    ma = (anchor_w >= lo) & (anchor_w <= hi) & np.isfinite(anchor_f) & (anchor_f > 0)
    mb = (w >= lo) & (w <= hi) & np.isfinite(f) & (f > 0)
    if ma.sum() < 5 or mb.sum() < 5:
        return 1.0, True
    fa = np.interp(w[mb], anchor_w[ma], anchor_f[ma])
    r = np.nanmedian(fa / f[mb])
    if not np.isfinite(r) or not (0.2 < r < 5.0):
        return 1.0, True
    return float(r), False


# per-grating resolution element (resel) = the flux-conserving bin that preserves the native resolving
# power. the COS medium-res gratings (G130M/G160M, R~15k) MUST be coadded at their resel, not the coarse
# display axis, or their lines wash into noise (COS handbook / CLASSY James 2022). FUV resel = 6 px,
# NUV = 3 px, x the grating dispersion.
COS_RESEL = {'G130M': 0.060, 'G160M': 0.073, 'G140L': 0.48, 'G230L': 1.17,
             'G185M': 0.11, 'G225M': 0.10, 'G285M': 0.12}

# stis low-res modes already sit at the native px scale on COMMON_AXIS (tier-1 validated), so they
# coadd there. everything else (stis medium-res G*M / G230MB + echelle E*) is 9-68x downsampled on
# COMMON_AXIS and gets a resel coadd instead, same as the cos medium-res gratings.
LOWRES_STIS = {'G230LB', 'G430L', 'G750L', 'G140L', 'G230L'}


def resel_step(grating):
    return COS_RESEL.get(str(grating).upper())


def resample_fcr(w, f, grid, e=None):
    # the PI's flux-conserving resample (specutils FluxConservingResampler, nan_fill), per exposure.
    # w must be monotonic increasing (clean_wf guarantees it) - that avoids the splice-edge failure
    # that forced the np.interp fallback in phase 1. lazy import so non-coadd scripts don't pay for it.
    import astropy.units as u
    from astropy.nddata import StdDevUncertainty
    from specutils.manipulation import FluxConservingResampler
    try:
        from specutils import Spectrum
    except ImportError:
        from specutils import Spectrum1D as Spectrum
    fu = u.Unit('erg cm-2 s-1 AA-1')
    unc = StdDevUncertainty(np.asarray(e, float) * fu) if e is not None else None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sp = Spectrum(spectral_axis=np.asarray(w, float) * u.angstrom,
                      flux=np.asarray(f, float) * fu, uncertainty=unc)
        ns = FluxConservingResampler(extrapolation_treatment='nan_fill')(sp, np.asarray(grid, float) * u.angstrom)
    fo = np.asarray(ns.flux.value, float)
    # FCR hands the uncertainty back as inverse-variance; convert to std-dev before returning.
    eo = None
    if unc is not None and ns.uncertainty is not None:
        eo = np.asarray(ns.uncertainty.represent_as(StdDevUncertainty).array, float)
    return fo, eo


def coadd_resel(specs, grating):
    # COS per-grating coadd at the resel (the science product). specs = [(w,f,e,dq), ...] observed-frame.
    # build the resel grid over the data range, FCR each exposure onto it, inverse-variance combine.
    cleaned = []
    for w, f, e, dq in specs:
        cw, cf, ce = clean_wf(w, f, e, dq)
        if len(cw) >= 2:
            if ce is None:
                ce = np.full_like(cf, np.nanmedian(np.abs(cf)) or 1.0)
            cleaned.append((cw, cf, ce))
    if not cleaned:
        return np.array([]), np.array([]), np.array([])
    step = resel_step(grating)
    if step is None:   # stis medium-res/echelle: 2 native px (the resel) measured off the data
        d = np.nanmedian([np.nanmedian(np.abs(np.diff(c[0]))) for c in cleaned if len(c[0]) > 1])
        step = float(2 * d) if np.isfinite(d) and d > 0 else 0.1
    wmin = min(c[0].min() for c in cleaned); wmax = max(c[0].max() for c in cleaned)
    grid = np.arange(np.floor(wmin), np.ceil(wmax) + step, step)
    fl, er = [], []
    for cw, cf, ce in cleaned:
        rf, re = resample_fcr(cw, cf, grid, ce)
        fl.append(rf); er.append(re if re is not None else np.full_like(rf, np.nan))
    if len(fl) == 1:
        return grid, fl[0], er[0]
    comb, cerr = ivar_combine(fl, er)
    return grid, comb, cerr


def coadd_exposures(specs, axis=COMMON_AXIS):
    # specs: list of (w, f, e, dq). per-exposure flux-conserving resample (the PI's FCR, like the cos
    # path) onto the axis, then inverse-variance combine. FCR conserves flux + propagates error right
    # (np.interp interpolated the error array, wrong); no-op on native low-res, cleans the 1.7x G140L.
    fl, er = [], []
    for w, f, e, dq in specs:
        cw, cf, ce = clean_wf(w, f, e, dq)
        if len(cw) < 2:
            continue
        if ce is None:
            ce = np.full_like(cf, np.nanmedian(np.abs(cf)) or 1.0)
        rf, re = resample_fcr(cw, cf, axis, ce)
        fl.append(rf); er.append(re if re is not None else np.full_like(rf, np.nan))
    if not fl:
        return axis, np.full_like(axis, np.nan), np.full_like(axis, np.nan)
    if len(fl) == 1:
        return axis, fl[0], er[0]
    comb, cerr = ivar_combine(fl, er)
    return axis, comb, cerr
