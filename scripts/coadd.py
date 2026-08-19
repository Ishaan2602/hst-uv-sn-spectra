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

# a grid px is "in a gap" if the input pair bracketing it is wider than GAP_K native px. the FCR only
# nan-fills the OUTER range, so without this it bridges internal gaps (cos stripe tiling, dead stripe
# edges, echelle order gaps) with bogus flux. 5 px keeps small dq dropouts bridged, nan's the real gaps.
GAP_K = 5.0


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


# per-grating trim boundaries (option A, OFF by default): cut each grating just before its throughput
# dies so the drooped edge can't bleed into the join (the reference cuts). the ivar combine gets the same
# effect for free (the dead edge carries a huge error and self-down-weights), so this stays here only as
# a selectable alternative.
STIS_TRIM = {'G230LB': (None, 3050), 'G430L': (3040, None), 'G750L': (5400, None)}


def align_and_merge(legs, prefer=(), method='ivar', trim=None):
    # legs: (label, w, f) or (label, w, f, e), each on its OWN grid (native or resel). build a union
    # grid from the legs' own points (native-fine in every region, no coarsening), FCR each leg onto
    # it, combine. NO inter-grating scaling - the flux cal agrees in the overlaps (checked broadly on
    # every grating pair; the reference doesn't scale, it trims + median-combines).
    #   'ivar'      - inverse-variance: a drooped low-throughput edge self-down-weights at the join
    #   'nanmedian' - plain nanmedian (the reference combine)
    # `prefer` only tags the record. `trim`={label-substr:(lo,hi)} hard-cuts a leg first.
    # returns (grid, merged, merged_err, recs).
    prepared = []
    for leg in legs:
        lab = str(leg[0]); w = np.asarray(leg[1], float); f = np.asarray(leg[2], float)
        e = np.asarray(leg[3], float) if len(leg) > 3 and leg[3] is not None else None
        m = np.isfinite(w) & np.isfinite(f)
        if trim:
            for g, (lo, hi) in trim.items():
                if g in lab:
                    if lo is not None:
                        m &= w >= lo
                    if hi is not None:
                        m &= w <= hi
        w, f = w[m], f[m]; e = e[m] if e is not None else None
        if len(w) < 2:
            continue
        o = np.argsort(w); w, f = w[o], f[o]; e = e[o] if e is not None else None
        keep = np.concatenate([[True], np.diff(w) > 0])
        w, f = w[keep], f[keep]; e = e[keep] if e is not None else None
        prepared.append((lab, w, f, e))
    if not prepared:
        return np.array([]), np.array([]), np.array([]), []
    recs = [{'leg': lab, 'scale': 'none'} for lab, _, _, _ in prepared]
    grid = np.unique(np.concatenate([w for _, w, _, _ in prepared]))
    grid = grid[np.concatenate([[True], np.diff(grid) > 1e-6])]   # strictly increasing for the FCR
    fl, er = [], []
    for lab, w, f, e in prepared:
        rf, re = resample_fcr(w, f, grid, e)
        fl.append(rf); er.append(re if re is not None else np.full_like(rf, np.nan))
    F = np.vstack(fl)
    if method == 'ivar' and all(np.any(np.isfinite(x)) for x in er):
        E = np.vstack(er)
        with np.errstate(invalid='ignore', divide='ignore'):
            wgt = 1.0 / np.where(E > 0, E ** 2, np.nan)
            ok = np.isfinite(F) & np.isfinite(wgt)
            num = np.nansum(np.where(ok, F * wgt, 0.0), axis=0)
            den = np.nansum(np.where(ok, wgt, 0.0), axis=0)
            merged = np.where(den > 0, num / den, np.nan)
            merr = np.where(den > 0, np.sqrt(1.0 / den), np.nan)
    else:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            merged = np.nanmedian(F, axis=0)
            merr = np.full(merged.shape, np.nan)
    return grid, merged, merr, recs


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
    # nan the internal gaps the FCR just bridged (see GAP_K) so the ivar combine skips them, not averages
    # a bogus interpolated value across a stripe/order gap.
    ws = np.asarray(w, float)
    if len(ws) >= 2:
        gg = np.asarray(grid, float)
        dw = np.median(np.diff(ws))
        idx = np.clip(np.searchsorted(ws, gg), 1, len(ws) - 1)
        gap = (gg >= ws[0]) & (gg <= ws[-1]) & ((ws[idx] - ws[idx - 1]) > GAP_K * dw)
        if gap.any():
            fo = fo.copy(); fo[gap] = np.nan
            if eo is not None:
                eo = eo.copy(); eo[gap] = np.nan
    return fo, eo


def _clean_specs(specs):
    # clean + prep each exposure (dq/finite/sort/dedup) and give every leg a std-dev error (median
    # |flux| placeholder where the x1d had none) so the ivar combine always has weights.
    cleaned = []
    for w, f, e, dq in specs:
        cw, cf, ce = clean_wf(w, f, e, dq)
        if len(cw) >= 2:
            if ce is None:
                ce = np.full_like(cf, np.nanmedian(np.abs(cf)) or 1.0)
            cleaned.append((cw, cf, ce))
    return cleaned


def _native_disp(cleaned):
    # median native px dispersion across the cleaned exposures (A/px).
    d = np.nanmedian([np.nanmedian(np.abs(np.diff(c[0]))) for c in cleaned if len(c[0]) > 1])
    return float(d) if np.isfinite(d) and d > 0 else None


def _coadd_on_step(cleaned, step):
    # FCR each cleaned exposure onto a uniform grid of the given step over the data range, ivar combine.
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


def coadd_native(specs, grating=None):
    # native (raw) per-grating coadd: FCR each exposure onto a 1-native-px grid, ivar combine. the
    # priority product (raw/native over rebinning). specs = [(w,f,e,dq), ...] observed-frame.
    cleaned = _clean_specs(specs)
    if not cleaned:
        return np.array([]), np.array([]), np.array([])
    d = _native_disp(cleaned)
    if d is None:
        return np.array([]), np.array([]), np.array([])
    return _coadd_on_step(cleaned, d)


def coadd_resel(specs, grating, cos=False):
    # per-grating coadd at the resel. COS (cos=True) uses the handbook resel (6px FUV / 3px NUV) from
    # COS_RESEL; STIS uses 2 native px (near-Nyquist). the grating names G230L/G140L exist on BOTH
    # instruments, so only trust COS_RESEL when the caller says it's cos - else always 2x native.
    cleaned = _clean_specs(specs)
    if not cleaned:
        return np.array([]), np.array([]), np.array([])
    step = resel_step(grating) if cos else None
    if step is None:   # stis (any grating), or a cos grating missing from the dict: 2 native px
        d = _native_disp(cleaned)
        step = float(2 * d) if d is not None else 0.1
    return _coadd_on_step(cleaned, step)


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
