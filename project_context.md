# Project Context: Comprehensive UV Spectroscopic Repository of Supernovae

## 1. Executive Summary
This project aims to construct a unified, uniformly reduced, and scientifically validated repository of all supernova (SN) ultraviolet (UV) spectra observed by the Hubble Space Telescope (HST). The repository spans every HST UV spectrograph, principally the Space Telescope Imaging Spectrograph (STIS) and the Cosmic Origins Spectrograph (COS), and is explicitly NOT restricted to STIS. Early-time UV spectroscopy captures critical high-ionization diagnostics (e.g., N IV, C IV, He II) that reveal the mass-loss histories of Red Supergiant (RSG) progenitor stars in the years immediately preceding core collapse. This legacy dataset will provide empirical constraints to resolve "missing links" in stellar evolution and serve as groundwork for upcoming missions like the Ultraviolet Explorer (UVEX).

---

## 2. Scientific & Modeling Framework
* **Target Domain:** Type II Supernovae (SNe II) and associated Core-Collapse events showing early-time interaction.
* **Key Observational Method:** **Flash Spectroscopy**. The initial supernova shockwave ionizes the immediate Circumstellar Material (CSM), acting as a flash lamp that illuminates the dense, ambient stellar wind.
* **Physical Parameters to Extract:**
    * Mass-loss rates ($\dot{M}$)
    * Circumstellar Material radii ($R_{CSM}$)
    * Chemical abundance variations in the progenitor's local environment.
* **Modeling Framework:** Custom data-reduced spectra will be systematically compared to synthetic spectral models generated via **CMFGEN**, a non-LTE radiative transfer code, to break degeneracies between velocity, density, and ionization states.

---

## 3. Core Software & Pipeline Stack
* **Data Access Environment:** Virtual Observatory queries via `astroquery.mast`, plus direct MAST archive retrieval and the higher-level products from HASP and the HSLA (see `data_access.md` for how these relate).
* **Reduction Engines (instrument dependent):**
    * STIS: STScI `stistools` core package, utilizing `stistools.ocrreject` (combining split exposures, masking cosmic ray signatures) and `stistools.x1d` (1D spatial extraction, background modeling, and calibration from 2D images).
    * COS: STScI `calcos` pipeline and `costools` for custom extraction; COS is photon-counting so its reduction path differs from STIS.
* **Analysis Libraries:** `astropy.table`, `astropy.io.fits`, `specutils`, and standard scientific Python environments (`numpy`, `scipy`, `matplotlib`).

---

## 4. Key Programmatic Constraints & References
* **Primary Contact/Mentor:** Dr. Wynn Jacobson-Galán
* **Faculty Sponsor:** Prof. Mansi Kasliwal
* **Milestones & Deliverables:**
    * **Interim Report 1:** July 3, 2026 (End of Week 3)
    * **Interim Report 2 / Abstract:** July 31, 2026 (End of Week 7)
    * Theres also a **SURF Oral Presentation:** and **Final Report:** towards end of 10-week SURF.