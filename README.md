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
| `output/` | Reduced spectra and per-SN data products (canonical) |

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
output/
  {SN}/
    {SN}_products.json                         # index of every product available for this SN (read this first)
    {SN}_manifest.json                         # epoch list, instruments, n_epochs
    {SN}_{stis,cos}_manifest.json              # per-instrument manifests
    {SN}_scaling.csv                           # inter-grating flux scale factors
    {SN}_emission.json                         # Mg II 2800 + Ly-a 1216 line fluxes + shape models
    {SN}_absorption.json                       # curated ISM absorption columns / metallicity
    {SN}_fuv_aod.json                          # resolved COS FUV apparent-optical-depth ISM columns
    {SN}_lya_nhi.json                          # N(HI) from the damped Ly-a fit
    {SN}_naid_ebv.json                         # Na I D host reddening E(B-V) (survey SNe)
    {SN}_timeseries.png                        # UV time-series waterfall
    epochs/
      {SN}_{date}_dayN_native.txt              # per-epoch coadd, native sampling
      {SN}_{date}_dayN_resel.txt               # per-epoch coadd, resel sampling
    emission/
      {SN}_{grating}_dayN_{mgii,lya}_diag.png  # per-epoch emission-line fit diagnostics
      {SN}_emission_summary.png
    absorption/
      {SN}_{grating}_dayN_cog.csv              # per-ion ISM columns; the Fe II row carries
      {SN}_{grating}_dayN_lines.csv            #   adopted_eligible + exclude_reason so a photospheric
      {SN}_{grating}_dayN_{cog,cont}.png       #   fit is never mistaken for a foreground ISM column
      {SN}_{grating}_dayN_aod.{csv,png}        # FUV apparent-optical-depth (COS G130M/G160M)
      {SN}_lya_nhi_fit.png                     # damped Ly-a N(HI) fit diagnostic
    reddening/
      {SN}_naid_ew.png                         # Na I D equivalent-width fit diagnostic
    STIS/  COS/                                # raw per-grating reduction tree (1d, native, resel, epochcoadd)
```

`{SN}_products.json` lists which analysis threads produced a product for the SN and the path to each
file - read it to discover what is available instead of learning the naming grammar. All product files
key off the SN directory name.

Analysis summary tables live in `catalog/` alongside the main target table:
`emission_summary.csv`, `absorption_summary.csv`, `ism_cog_summary.csv`, `lya_nhi_summary.csv`,
`fuv_aod_summary.csv`, `naid_ebv_summary.csv`.

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
`catalog/lya_nhi_summary.csv` contains automated N(HI) measurements from damped Lya.

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
