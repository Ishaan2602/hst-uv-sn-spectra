import os
import numpy as np
import matplotlib
matplotlib.use('Agg')   # headless, we only save pngs
import matplotlib.pyplot as plt

# per-grating extrsize. stis ccd psf broadens toward the red (charge diffusion + optical psf),
# so a fixed box doesn't fit every grating. measured off 2024iss enclosed-flux curves:
# g230lb es7 ~94%, g430l es9 ~96%, g750l es7 clips to 80% at the red end -> need es15 (91-98%).
# mama source core is fwhm~4px on finer pixels -> es11.
EXTRSIZE = {'G230LB': 7, 'G430L': 9, 'G750L': 15, 'G750M': 15, 'G430M': 9,
            'G140L': 11, 'G230L': 11, 'G140M': 11, 'G230M': 11, 'G230MB': 7}

# adaptive-widen clamps for the ccd gratings (default anchor -> up to the ceiling). from the
# enclosed-flux curves: g750l needs the room to reach the red-end plateau, blue gratings less so.
EXTRSIZE_CLAMP = {'G230LB': (7, 15), 'G430L': (9, 17), 'G750L': (15, 25), 'G750M': (15, 25),
                  'G430M': (9, 17), 'G230MB': (7, 15)}

ECHELLE = ('E140M', 'E230M', 'E140H', 'E230H')


def is_echelle(grating):
    return str(grating).upper() in ECHELLE


def extrsize_for(grating, detector='CCD'):
    return EXTRSIZE.get(str(grating).upper(), 7 if detector == 'CCD' else 11)


def trace_window(detector, grating):
    # spatial row window to search for the peak. ccd trace sits at E1 (~890-912) or center (~510),
    # so keep it broad. mama splits by band: fuv ~360-435, nuv ~469-552 (tighten so host/glow at
    # other rows can't pull the center off the source -- this is the SN2005IP G230L fix).
    g = str(grating or '').upper()
    if detector == 'CCD':
        return (100, 1000)
    if g in ('G140L', 'G140M'):
        return (320, 470)
    return (430, 590)


def adaptive_extrsize(sci, center, grating, detector, target=0.98, flags=None):
    # hybrid: per-grating default is the floor; for a cleanly-detected ccd trace, widen to enclose
    # `target` of the cross-dispersion flux, clamped per grating. mama/other gratings stay fixed
    # (airglow pedestal makes the enclosed curve meaningless). falls back to default on any finder flag.
    g = str(grating).upper()
    default = extrsize_for(g, detector)
    bad = {'finder_fail', 'extended_host', 'low_detect'}
    if detector != 'CCD' or g not in EXTRSIZE_CLAMP or (flags and bad & set(flags)):
        return default, 'fixed'
    prof = np.nanmedian(sci[:, 300:800], axis=1)
    sub = np.clip(prof - np.nanmedian(prof), 0, None)
    c = int(round(center)); win = 30
    tot = sub[c - win:c + win + 1].sum()
    if tot <= 0:
        return default, 'fixed(no_flux)'
    hw = next((h for h in range(1, win) if sub[c - h:c + h + 1].sum() / tot >= target), win)
    lo, hi = EXTRSIZE_CLAMP[g]
    es = int(np.clip(max(2 * hw + 1, default), lo, hi))
    return es, ('adaptive_clamphi' if es == hi and 2 * hw + 1 > hi else f'adaptive_es{es}@{target:.2f}')


def find_trace(sci, detector='CCD', extrsize=None, refcol=512, step=40, grating=None):
    # per-band track the trace across dispersion, fit a deg-2 tilt, report a2center at col 512.
    # ccd: median a 300-800 col band (drops warm cols). mama: sum all cols (fuv flux is line-only).
    ny, nx = sci.shape
    extrsize = extrsize or (7 if detector == 'CCD' else 11)
    coarse = np.nanmedian(sci[:, 300:800], axis=1) if detector == 'CCD' else np.nansum(sci, axis=1)
    slo, shi = trace_window(detector, grating)
    shi = min(shi, ny)
    c0 = int(np.argmax(coarse[slo:shi]) + slo)
    bg = np.nanmedian(coarse); sd = np.nanstd(coarse)
    detect = (coarse[c0] - bg) / (sd + 1e-9)
    bx, bc = [], []
    for x0 in range(step, nx - step, step):
        sub = (np.nanmedian(sci[:, x0 - step // 2:x0 + step // 2], axis=1) if detector == 'CCD'
               else np.nansum(sci[:, x0 - step // 2:x0 + step // 2], axis=1))
        lo, hi = max(0, c0 - 25), min(ny, c0 + 25)
        p = int(np.argmax(sub[lo:hi]) + lo)
        if (sub[p] - np.nanmedian(sub)) > 3 * np.nanstd(sub):
            bx.append(x0); bc.append(p)
    if len(bx) >= 4:
        coef = np.polyfit(bx, bc, 2)
        curve = np.polyval(coef, np.arange(nx))
    else:
        curve = np.full(nx, float(c0))
    # read a2center off the trace curve at the reference column. 512 (detector centre) for a full
    # 1024 frame, which is where x1d anchors. only if that is out of range (a narrow subarray) fall
    # back to the subarray's OWN centre, not its last column (the edge, where the trace narrows).
    if refcol >= nx:
        refcol = (nx - 1) // 2
    a2center = float(curve[refcol])
    sub = np.clip(coarse - bg, 0, None); win = 25
    tot = sub[c0 - win:c0 + win + 1].sum()
    encl = sub[c0 - extrsize // 2:c0 + extrsize // 2 + 1].sum() / tot if tot > 0 else 0.0
    cc = np.cumsum(sub[c0 - win:c0 + win + 1]) / tot if tot > 0 else np.zeros(2 * win + 1)
    rows = np.arange(c0 - win, c0 + win + 1)
    hw90 = (rows[np.searchsorted(cc, 0.95)] - rows[np.searchsorted(cc, 0.05)]) / 2 if tot > 0 else np.nan
    flags = []
    if detector == 'MAMA' and detect < 4.5: flags.append('low_detect')
    if detector == 'CCD' and encl < 0.5:    flags.append('finder_fail')
    if hw90 > 6 and detector == 'CCD':      flags.append('extended_host')
    return {'a2center': round(a2center, 1), 'curve': curve, 'nbands': len(bx),
            'detect': round(detect, 1), 'encl': round(encl, 3),
            'hw90': round(hw90, 1) if not np.isnan(hw90) else np.nan, 'flags': flags}


def adaptive_bg(sci, center, extrsize, size=10, gap=8, reach=40):
    # close straddling, just outside the aperture (matches the +-14 we liked on 2024iss).
    # nudges out only if the close window is contaminated (host edge, nearby source).
    prof = np.nanmedian(sci[:, 300:800], axis=1)
    floor = np.nanpercentile(prof, 20)
    scat = np.nanstd(prof[prof < np.nanpercentile(prof, 80)])
    base = extrsize // 2 + gap
    def side(sign):
        for off in range(base, reach):
            c = center + sign * off
            if c - size // 2 < 0 or c + size // 2 > len(prof) - 1:
                continue
            lvl = np.nanmedian(prof[c - size // 2:c + size // 2 + 1])
            if lvl < floor + 3 * scat:
                return c
        return center + sign * base
    return side(-1), side(+1)


def draw_extraction(sci, res, detector='CCD', extrsize=None, title='', save=None):
    # curved trace + aperture + bg on the 2d, plus the ref-col spatial profile. saves as png.
    extrsize = extrsize or (7 if detector == 'CCD' else 11)
    curve = res['curve']; a2 = res['a2center']; nx = sci.shape[1]; xs = np.arange(nx)
    b1, b2 = adaptive_bg(sci, int(round(a2)), extrsize)
    o1, o2 = b1 - int(round(a2)), b2 - int(round(a2))
    fig, axes = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={'width_ratios': [3, 1]})
    vlo, vhi = np.nanpercentile(sci, [5, 98])
    axes[0].imshow(sci, origin='lower', aspect='auto', vmin=vlo, vmax=vhi, cmap='gray')
    axes[0].plot(xs, curve, color='crimson', lw=1.2, label=f'trace (a2c {a2:.1f})')
    axes[0].plot(xs, curve - extrsize / 2, color='crimson', lw=0.6, ls='--')
    axes[0].plot(xs, curve + extrsize / 2, color='crimson', lw=0.6, ls='--')
    for off, c, lab in [(o1, 'orange', f'bg1 ({o1:+d})'), (o2, 'gold', f'bg2 ({o2:+d})')]:
        axes[0].plot(xs, curve + off - 5, color=c, lw=0.6)
        axes[0].plot(xs, curve + off + 5, color=c, lw=0.6, label=lab)
    axes[0].set_ylim(a2 - 45, a2 + 45); axes[0].set_xlim(0, nx)
    axes[0].set_xlabel('column'); axes[0].set_ylabel('row')
    axes[0].legend(fontsize=8, loc='upper right')
    flagtxt = (' FLAGS: ' + ','.join(res['flags'])) if res['flags'] else ''
    axes[0].set_title(f'{title}  detect {res["detect"]}{flagtxt}', fontsize=9)
    # rectify each column by the trace tilt before collapsing, otherwise a tilted trace smears the
    # profile and its peak drifts off the a2center-drawn aperture (the collapse is over many columns).
    shifts = np.round(a2 - curve).astype(int)
    rect = np.empty_like(sci)
    for j in range(nx):
        rect[:, j] = np.roll(sci[:, j], shifts[j])
    prof = np.nanmedian(rect[:, 300:800], axis=1) if detector == 'CCD' else np.nansum(rect, axis=1)
    axes[1].plot(prof, np.arange(len(prof)), 'k', lw=0.7)
    axes[1].axhspan(a2 - extrsize / 2, a2 + extrsize / 2, color='crimson', alpha=0.2)
    axes[1].axhspan(a2 + o1 - 5, a2 + o1 + 5, color='orange', alpha=0.2)
    axes[1].axhspan(a2 + o2 - 5, a2 + o2 + 5, color='gold', alpha=0.2)
    axes[1].set_ylim(a2 - 45, a2 + 45); axes[1].set_xlabel('counts'); axes[1].set_title('profile', fontsize=9)
    plt.tight_layout()
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        plt.savefig(save, dpi=90, bbox_inches='tight')
    plt.close(fig)
    return (o1, o2)
