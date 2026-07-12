# Supernova Spectral Physics & P-Cygni Formation

## Core Concept: P-Cygni Profile Formation
*(Note: These notes discuss how a single isolated atomic transition forms a P-Cygni profile)*

* **Definition:** The characteristic spectral signature of an expanding shell of gas surrounding a central continuum source.
* **Supernova Application:** Traces the boundary interaction between the optically thick inner ejecta (photosphere) and the optically thin outer ejecta (line-forming envelope).
* **Profile Morphology:** Composed of a blueshifted absorption trough paired with an emission peak centered near the rest wavelength.

## Physical Setup & Geometry
1. **Central Photosphere:** The dense, innermost core where optical depth $\tau \approx 2/3$. It emits a continuous, blackbody-like spectrum that provides the background light source.
2. **Expanding Envelope (Line Opacity Region):** Optically thin outer gas expanding radially outward at high velocities (10,000 to 30,000 km/s). Contains specific ions undergoing bound-bound atomic transitions.

**Coordinate Framework:**
* **Impact Parameter ($p$):** Perpendicular distance from the line of sight to the center.
* **Line-of-Sight Velocity ($z$):** Velocity component directed toward or away from the observer.

## Formation Mechanisms

### 1. The Blueshifted Absorption Trough
* **Origin:** Gas positioned directly between the observer and the photosphere (the "Absorption Region" where $p < R_{\text{photosphere}}$ and $z < 0$).
* **Mechanism:** This gas is moving rapidly toward the observer. It absorbs continuum photons from the photosphere at the rest wavelength of the atomic transition, but due to the Doppler effect, this absorption is blueshifted to shorter wavelengths.

### 2. The Emission Peak
* **Origin:** Gas in the surrounding "Emission Region" (where $p > R_{\text{photosphere}}$) that is not directly in front of or behind the photosphere.
* **Mechanism:** Photons absorbed by the envelope are re-emitted isotropically. The observer sees emission from gas moving both toward (blueshifted) and away (redshifted) from them. This creates a broad emission profile roughly centered on the rest wavelength.

### 3. The Occluded Region
* **Origin:** Gas behind the photosphere (where $p < R_{\text{photosphere}}$ and $z > 0$).
* **Mechanism:** Any emission from this region is physically blocked by the optically thick photosphere and does not reach the observer. This truncates the redshifted "red tail" of the emission peak.

---

## Important Spectral Lines & Rest Wavelengths ($\AA$)

### Calcium (Ca II)
* **Ca II H and K:**
  * 3933.66
  * 3968.47
* **Ca II IR Triplet:**
  * 8498.02
  * 8542.09
  * 8662.14

### Sulfur (S II)
* 4815.55
* 4924.11
* 5032.43
* 5428.65
* 5432.80
* 5453.85
* 5473.61
* 5606.15
* 5640.35

### Iron (Fe II)
*(Note: Fe II features heavily in SN spectra across the optical/UV, exact transitions to be appended as needed based on spectral range).*

---

## Instruments: COS vs STIS resolution

* **COS generally has slightly higher spectral resolution than STIS** (in the UV, comparing like-for-like modes). COS was built for faint UV point sources at high throughput/resolution.
* This higher resolution is why COS is the go-to for resolving the narrow lines of **hot stars** (e.g. hot subdwarfs / OB stars in the Milky Way) and narrow ISM / intergalactic absorption. STIS trades some UV resolution for its much wider wavelength grasp and spatial (long-slit) coverage.
* Practical takeaway for us: for a compact/point-like SN where we care about narrow features, COS resolves them better; for broad SN features + full NUV-NIR SED coverage, STIS is the workhorse.
