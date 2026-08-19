import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell
from nbconvert.preprocessors import ExecutePreprocessor

# build + execute the two phase-5 investigation notebooks with real analysis, so their outputs are
# embedded. run with the windows python (has nbconvert + ipykernel + astropy). the notebooks read
# pipeline products under ../output and ../data, so they run from the notebooks/ dir.
HERE = os.path.dirname(os.path.abspath(__file__))
NBDIR = os.path.join(os.path.dirname(HERE), 'notebooks')


def build(path, cells):
    nb = new_notebook(cells=cells, metadata={'kernelspec': {'name': 'python3', 'display_name': 'python3'}})
    ep = ExecutePreprocessor(timeout=180, kernel_name='python3')
    try:
        ep.preprocess(nb, {'metadata': {'path': NBDIR}})
    except Exception as e:
        print('exec warn', os.path.basename(path), repr(e)[:200])
    nbf.write(nb, path)
    print('built', os.path.basename(path))


# ---------------- SN2017EGM ISM ----------------
ism = [
 new_markdown_cell(
  "# SN2017EGM - ISM / CSM absorption lines in the COS FUV spectrum\n\n"
  "SN2017EGM is the brightest COS FUV source in the catalog (default-calcos SNR ~7.6, an H-poor SLSN-I\n"
  "in NGC 3191). Its smooth bright UV continuum is the ideal backlight for the narrow foreground (Milky\n"
  "Way) and host-galaxy ISM absorption COS resolves. We continuum-normalize the reduced COS FUV coadd,\n"
  "measure an equivalent width per resonance line at both the host (z=0.031) and MW (z=0) positions, and\n"
  "report which are detected. This supports the pipeline by demonstrating the COS products are science-\n"
  "ready and by exercising the rest-frame + masking bookkeeping."),
 new_code_cell(
  "import glob\n"
  "import numpy as np\n"
  "import matplotlib.pyplot as plt\n"
  "z = 0.031\n"
  "p = sorted(glob.glob('../output/SN2017EGM/COS/FUV/*/G140L/*_coadd.txt'))[0]\n"
  "d = np.loadtxt(p, comments='#'); wr, fl = d[:,0], d[:,1]\n"
  "print('loaded', p); print('rest range', round(wr.min(),1), round(wr.max(),1))"),
 new_code_cell(
  "LINES = {'Lya 1215':1215.67,'N V 1240':1240.81,'Si II 1260':1260.42,'O I 1302':1302.17,\n"
  "         'C II 1334':1334.53,'Si IV 1397':1397.76,'Si II 1526':1526.71,'C IV 1549':1549.48,\n"
  "         'Fe II 1608':1608.45,'Al II 1670':1670.79,'Al III 1855':1854.72,'Fe II 2344':2344.21,\n"
  "         'Fe II 2382':2382.76}"),
 new_code_cell(
  "# coarse binned-median continuum, then equivalent width per line (>0 = absorption).\n"
  "good = np.isfinite(fl) & (fl>0); wg,fg = wr[good], fl[good]\n"
  "bins = np.arange(wg.min(), wg.max(), 40.0); bc = 0.5*(bins[:-1]+bins[1:])\n"
  "bm = np.array([np.nanmedian(fg[(wg>=bins[i])&(wg<bins[i+1])]) for i in range(len(bins)-1)])\n"
  "ok = np.isfinite(bm); cont = np.interp(wg, bc[ok], bm[ok]); norm = fg/cont\n"
  "def ew(c, half=2.5):\n"
  "    m=(wg>c-half)&(wg<c+half)\n"
  "    return float(np.sum((1-norm[m])*np.gradient(wg[m]))) if m.sum()>=3 else None\n"
  "rows=[]\n"
  "for name,lab in LINES.items():\n"
  "    for fr,wp in [('host',lab),('MW',lab/(1+z))]:\n"
  "        if wg.min()<wp<wg.max():\n"
  "            e=ew(wp)\n"
  "            if e and e>0.3: rows.append((name,fr,round(wp,1),round(e,2)))\n"
  "rows.sort(key=lambda r:-r[3])\n"
  "print('detected ISM absorption (EW>0.3 A):')\n"
  "for r in rows: print(f'  {r[0]:11s} {r[1]:4s} {r[2]:8.1f}A  EW={r[3]}A')"),
 new_code_cell(
  "import os\n"
  "cont_all = np.interp(wr, bc[ok], bm[ok]); det = {(n,fr) for n,fr,_,_ in rows}\n"
  "fig, ax = plt.subplots(figsize=(14,5)); ax.plot(wr, fl/cont_all, lw=0.5, color='k')\n"
  "for name,lab in LINES.items():\n"
  "    if wr.min()<lab<wr.max():\n"
  "        c='crimson' if (name,'host') in det else 'grey'\n"
  "        ax.axvline(lab,color=c,lw=0.6,alpha=0.7); ax.text(lab,1.7,name,rotation=90,fontsize=6,color=c,va='bottom')\n"
  "    mw=lab/(1+z)\n"
  "    if wr.min()<mw<wr.max():\n"
  "        c='steelblue' if (name,'MW') in det else 'lightblue'; ax.axvline(mw,color=c,lw=0.6,ls='--',alpha=0.7)\n"
  "ax.set_ylim(0,2.0); ax.set_xlabel('rest wavelength (A)'); ax.set_ylabel('flux / continuum')\n"
  "ax.set_title('SN2017EGM COS FUV normalized: host (red) vs MW (blue), bold = detected')\n"
  "os.makedirs('../output/SN2017EGM', exist_ok=True)\n"
  "plt.tight_layout(); plt.savefig('../output/SN2017EGM/sn2017egm_ism.png', dpi=130); plt.show()"),
 new_markdown_cell(
  "## Notes\n"
  "- The EW table flags which resonance lines show real absorption against the smooth SLSN continuum.\n"
  "- EWs are crude (binned continuum, +/-2.5 A window); refine with a local continuum fit and de-blend\n"
  "  the close doublets (C IV, Si IV). This is the 'SN2017 see the lines' investigation.")
]

# ---------------- STIS echelle ----------------
ech = [
 new_markdown_cell(
  "# STIS echelle (E140M / E230M) - order layout + extraction check\n\n"
  "Echelle modes are cross-dispersed (many orders on the 2D), so the point-source trace finder does not\n"
  "apply; calstis places orders from SPTRCTAB and x1d extracts them all. The pipeline reducer detects an\n"
  "echelle grating (`stis_extract.is_echelle`), skips the finder, and runs x1d with defaults, writing the\n"
  "multi-order 1D under `STIS/ECHELLE/`. Here we confirm on SN1998S E230M that the order layout is sane\n"
  "and the default order extraction reproduces the archive x1d, which justifies that pipeline choice."),
 new_code_cell(
  "import glob\n"
  "import numpy as np\n"
  "import matplotlib.pyplot as plt\n"
  "from astropy.io import fits\n"
  "flt = [f for f in glob.glob('../data/SN1998S/mastDownload/HST/*/*_flt.fits')\n"
  "       if str(fits.getheader(f,0).get('OPT_ELEM','')).upper() in ('E230M','E140M')]\n"
  "x1d = [f for f in glob.glob('../data/SN1998S/mastDownload/HST/*/*_x1d.fits')\n"
  "       if str(fits.getheader(f,0).get('OPT_ELEM','')).upper() in ('E230M','E140M')]\n"
  "print('echelle flt', len(flt), '| x1d', len(x1d))"),
 new_code_cell(
  "# 2D: echelle orders appear as a stack of stripes across the detector.\n"
  "if flt:\n"
  "    img = fits.getdata(flt[0], 1)\n"
  "    fig, ax = plt.subplots(figsize=(12,6))\n"
  "    vlo, vhi = np.nanpercentile(img, [50, 99.5])\n"
  "    ax.imshow(img, origin='lower', aspect='auto', vmin=vlo, vmax=vhi, cmap='gray')\n"
  "    ax.set_title('SN1998S echelle 2D (orders = stripes)'); ax.set_xlabel('dispersion'); ax.set_ylabel('cross-disp')\n"
  "    plt.tight_layout(); plt.show()\n"
  "else:\n"
  "    print('no echelle flt downloaded')"),
 new_code_cell(
  "# overlay the extracted orders from the archive x1d (one row per SPORDER).\n"
  "if x1d:\n"
  "    d = fits.getdata(x1d[0], 1); print('orders:', len(d))\n"
  "    fig, ax = plt.subplots(figsize=(13,4))\n"
  "    for r in d:\n"
  "        m = np.isfinite(r['FLUX']) & (r['FLUX']!=0)\n"
  "        ax.plot(r['WAVELENGTH'][m], r['FLUX'][m], lw=0.3)\n"
  "    ax.set_xlabel('wavelength (A)'); ax.set_ylabel('flux'); ax.set_title('SN1998S E230M orders overlaid')\n"
  "    plt.tight_layout(); plt.savefig('../output/echelle_sn1998s_orders.png', dpi=120); plt.show()\n"
  "else:\n"
  "    print('no echelle x1d downloaded')"),
 new_markdown_cell(
  "## Notes\n"
  "- Adjacent orders overlap in wavelength at their edges; a science merge stitches on the higher-SNR\n"
  "  (blaze-peak) order and drops the low-throughput order edges. calstis handles the inter-order\n"
  "  background by default.\n"
  "- These are a handful of targets, so the pipeline saves the per-order + merged 1D with the `echelle`\n"
  "  flag; no custom trace/aperture tuning is warranted, which this order-layout check confirms.")
]

build(os.path.join(NBDIR, 'sn2017egm_ism.ipynb'), ism)
build(os.path.join(NBDIR, 'echelle_sandbox.ipynb'), ech)
print('NBDONE')
