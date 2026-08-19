import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits

# g750l vs g230lb cross-dispersion width: does es7 clip g750l flux?
# measure enclosed-flux vs aperture at several columns for all three gratings on 2024iss.
g750 = fits.getdata('data/2024iss/mastDownload/HST/of8b02040/of8b02040_crj.fits', 1)
g230 = fits.getdata('data/2024iss/mastDownload/HST/of8b02010/of8b02010_crj.fits', 1)
g430 = fits.getdata('data/2024iss/mastDownload/HST/of8b02030/of8b02030_flt.fits', 1)

def prof_at(sci, col, half=15):
    band = np.nanmedian(sci[:, col-half:col+half], axis=1)
    bg = np.nanmedian(band)
    return band - bg

def enclosed(band, peak, hws):
    # fraction within +/-hw of peak, normalized to +/-20
    tot = band[peak-20:peak+21].clip(0).sum()
    return [band[peak-hw:peak+hw+1].clip(0).sum()/tot if tot>0 else 0 for hw in hws]

cols = [100, 300, 500, 700, 900]
hws = np.arange(0, 16)
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, sci, name in [(axes[0], g230, 'G230LB'), (axes[1], g430, 'G430L'), (axes[2], g750, 'G750L')]:
    for col in cols:
        band = prof_at(sci, col)
        peak = int(np.argmax(band[800:1000]) + 800)   # E1 region ~894
        enc = enclosed(band, peak, hws)
        ax.plot(2*hws+1, enc, marker='.', label=f'col {col} (pk {peak})')
    for es, lab in [(7,'es7'),(11,'es11'),(15,'es15')]:
        ax.axvline(es, color='gray', ls=':', lw=0.8)
        ax.text(es, 0.55, lab, fontsize=8, rotation=90)
    ax.axhline(0.95, color='red', ls='--', lw=0.6)
    ax.set_xlabel('extrsize (2*hw+1)'); ax.set_ylabel('enclosed fraction')
    ax.set_title(f'2024iss {name} enclosed flux vs aperture'); ax.legend(fontsize=8); ax.set_ylim(0.4, 1.02)
plt.tight_layout()
plt.savefig('output/g750l_width_check.png', dpi=95, bbox_inches='tight')

# print the es7/es11/es15 enclosed at each column
for sci, name in [(g230,'G230LB'), (g430,'G430L'), (g750,'G750L')]:
    print(f'\n{name}: enclosed at es7 / es11 / es15 by column')
    for col in cols:
        band = prof_at(sci, col); peak = int(np.argmax(band[800:1000])+800)
        e = enclosed(band, peak, [3,5,7])
        print(f'  col {col:4d} peak {peak}: {e[0]:.3f} / {e[1]:.3f} / {e[2]:.3f}')
print('\nsaved output/g750l_width_check.png')
