# HST UV Supernova Spectra

Uniformly re-reduced UV spectra of supernovae observed by the Hubble Space Telescope (STIS and COS). Unofficial data release from a Caltech SURF 2026 project; pipeline and products are updated as the survey expands.

---

## Layout

| Path | Contents |
|------|----------|
| `scripts/` | Reduction pipeline and analysis code |
| `catalog/` | SN target catalog and ISM curve-of-growth summary |
| `reference/` | Hand-curated inputs: host reddening, ISM column densities |
| `linelists/` | ISM and CSM line wavelength tables |
| `output5/` | Reduced spectra and per-SN data products (canonical) |

Raw HST FITS exposures are not included; they are publicly available on MAST and can be fetched with `scripts/download_all.py`.

---

## Spectra Format

Spectrum files (`*.txt`) are 3-column whitespace-delimited ASCII:

```
# obs_wvl  flux  error
1663.0  9.88e-14  2.12e-14
```

- **wavelength**: observer-frame Angstroms
- **flux / error**: erg s⁻¹ cm⁻² Å⁻¹, 1-sigma statistical

Filename suffix `_native` = instrument native pixel sampling.  
Filename suffix `_resel` = 2-pixel (resolution element) binning.

---

## Output Tree

```
output5/
  {SN}/
    {SN}_manifest.json                         # epoch list, instruments, n_epochs
    {SN}_{stis,cos}_manifest.json              # per-instrument manifests
    {SN}_scaling.csv                           # inter-grating flux scale factors
    {SN}_emission.json                         # emission line measurements
    {SN}_absorption.json                       # absorption line measurements
    STIS/
      CCD/{date_dayN}/
        {GRATING}/
          {root}_1d.txt                        # single-exposure extracted spectrum
          {SN}_{date}_{GRATING}_native.txt     # grating coadd, native sampling
          {SN}_{date}_{GRATING}_resel.txt      # grating coadd, resel sampling
        {SN}_{date}_epochcoadd_native.txt      # all-grating epoch coadd, native
        {SN}_{date}_epochcoadd_resel.txt       # all-grating epoch coadd, resel
      MAMA/...                                 # same layout, MAMA detector
    COS/...                                    # same layout
  absorption_summary.csv                       # catalog-level ISM fit results
```

Each epoch coadd also ships as a `.fits` alongside the `.txt`, and per-SN
diagnostic plots (`.png`: extraction traces, coadds, time series, line fits) sit
in the same tree.

---

## Catalog

`catalog/uv_sn_catalog_clean.csv` is the main target table. Key columns:

| Column | Description |
|--------|-------------|
| `name` | SN name (canonical) |
| `ra`, `dec` | J2000 coordinates (degrees) |
| `z` | Heliocentric redshift |
| `host_ebv` | Adopted host reddening E(B-V) |
| `instr` | Instruments with UV coverage |
| `gratings` | Gratings observed (semicolon-separated) |
| `n_spec` | Total number of spectra |
| `has_uv` | True if UV data is available |
| `flags` | Data quality notes |

`catalog/ism_cog_summary.csv` contains ISM curve-of-growth fit results per target.

---

## Reproducing the Reduction

Requires a working CRDS cache configured for HST:

```bash
export CRDS_SERVER_URL=https://hst-crds.stsci.edu
export CRDS_PATH=/path/to/crds_cache

pip install stistools calcos costools astropy astroquery specutils matplotlib

python scripts/download_all.py --sn SN2023IXF   # fetch raw data for one target
python scripts/run_full_catalog.py               # reduce the full catalog
```

Output lands in a new versioned directory. To promote a new reduction, update `CANONICAL` in `scripts/paths.py` — all analysis scripts read from there automatically.

---

## Authors

Ishaan Gurazada, Wynn Jacobson-Galán, Mansi Kasliwal

## License

BSD 3-Clause — see [LICENSE](LICENSE).
