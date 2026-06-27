*Source: Circumstellar Interaction in the Ultraviolet Spectra of SN 2023ixf 14-66 Days After Explosion (Bostroem et al., 2024)*
# Review of Section 1: Introduction

This section establishes the theoretical context of the study, highlighting the persistent uncertainties in stellar evolution and why SN 2023ixf is a benchmark event for understanding pre-supernova mass loss.

* **The Progenitor Problem:** Type II supernovae originate from red supergiants with initial masses between 8 and 25 solar masses. However, stellar evolution theory lacks a consensus on quantitative mass-loss rates for these stars prior to core collapse.
* **Signatures of Mass Loss:** When a red supergiant explodes into a dense circumstellar material (CSM), it produces "flash features"—narrow emission lines caused by the shock breakout ionizing the surrounding gas. These lines typically disappear within days to weeks.
* **Long-Term Interaction Effects:** Even after flash features fade, the collision between the supernova ejecta and the CSM converts kinetic energy into thermal energy. This interaction leads to:
    * A boosted overall luminosity.
    * A bluer spectral energy distribution.
    * High-velocity absorption features in the optical (e.g., H-alpha) and distinct emission lines in the ultraviolet (UV), such as Mg II.
* **The Significance of SN 2023ixf:** Located in the nearby galaxy M101 and caught within a day of explosion, this is the best-observed supernova with flash features to date. Early multi-wavelength data (X-ray, radio, UV, optical) suggest the progenitor was surrounded by a confined, dense CSM (Radius $\le 10^{15}$ cm) with a high mass-loss rate (approx. $10^{-2}$ to $10^{-4}$ solar masses per year).

# Review of Section 2: Observations and Data Reduction

### 1. Observation Strategy and Instrument Setup
* **Instrument:** Hubble Space Telescope (HST) using the Space Telescope Imaging Spectrograph (STIS) CCD detector.
* **Aperture & Positioning:** Utilized a $52" \times 0.2"$ slit at the E1 position (specifically chosen to mitigate flux loss from charge transfer efficiency degradation).
* **Gratings Used:** G230LB (UV), G430L (Optical), and G750L (Red/NIR).
* **Timeline:** Primary observations were taken on days 14, 19, 24, and 66 post-explosion. (A planned day 50 observation failed due to guide star acquisition issues and was successfully replaced by the day 66 visit).
* **Cosmic Ray Mitigation:** A minimum of four exposures were taken with each grating to allow for automatic cosmic-ray rejection.

### 2. Standard Pipeline Processing & Extraction
* **Automated Reduction:** Bias subtraction, flat-fielding, and cosmic-ray rejection were handled automatically by the Mikulski Archive for Space Telescopes (MAST) pipeline prior to download.
* **1D Extraction:** Most G230LB and G430L data were automatically extracted, wavelength-calibrated, and flux-calibrated.
* **Manual Intervention:** The day 66 G230LB observation had extremely low flux at the blue end, causing the automated pipeline to fail. This was manually extracted using the `stistools.calstis` pipeline by tracing the red end of the spectrum.

### 3. Grating-Specific Corrections (G750L Fringing)
* The G750L grating experiences known fringing at its reddest wavelengths.
* This was corrected by applying a contemporaneous fringe flat (using the `stistools.defringe` module) taken through a smaller $0.3" \times 0.09"$ aperture (standard for E1).
* Following the defringing process, 1D spectra for this grating were extracted using the `stistools.xld` routine.

### 4. Data Quality Control and Artifact Removal
* **Pixel Masking:** Bad pixels were strictly eliminated using data-quality flags 16 (high dark rate) and 512 (bad reference pixel).
* **False Emission Lines:** The authors carefully noted that several apparent "narrow emission features" present across all spectra were actually localized artifacts with high dark rates. By verifying against the 2D dark images, these artifacts were confidently flagged and excluded from scientific analysis to prevent false constraints on circumstellar material (CSM).

### 5. Integration of Supplementary Early-Time Data (Days 3-11)
* The study incorporated earlier HST observations (from GO-17205) to provide baseline context. However, these data suffered from pointing, acquisition, and CCD saturation issues (causing the supernova to be misaligned in the slit).
* **Flux Calibration Workarounds:**
    * To correct for variable flux loss across exposures in a single visit, the authors scaled the affected spectra to the spectrum with the highest flux (assuming it was the best-centered) using a first-degree polynomial fit.
    * G430L and G750L observations were aligned to the G230LB spectra using a constant percentage multiple.
    * The aligned fluxes were then median-combined (excluding bad pixels).
* **Saturation Limitations:** The authors acknowledge that major wavelength ranges in these early epochs (e.g., 3200-5000 Å in visit 1) are saturated and unreliable for absolute flux analysis. Consequently, these early spectra are utilized strictly for identifying the presence and timing of spectral features, not for precision flux measurements.

# Review of Section 3: Spectral Evolution

This section details the unprecedented, high-cadence timeline of SN 2023ixf's near-ultraviolet to near-infrared spectral changes between days 3 and 66.

### 1. Early Phase (Days 3-11): The Flash Features
* Initial spectra show narrow emission lines (e.g., N IV 1718 Å, C III 2297 Å) indicating the presence of dense, unshocked CSM.
* These narrow lines fade completely by day 8.

### 2. Transition & Middle Phase (Days 14-25): UV Persistence
* **Day 14:** The optical spectrum evolves into that of a typical young Type II supernova, dominated by broad P Cygni profiles (H I, He I). Unusually, significant UV continuum flux remains strong, which is rare as UV flux typically plummets quickly in classical Type II supernovae.
* **Day 20-25:** The UV continuum heavily fades by day 20, but absorption features (Mg II, Fe II/III) persist. A prominent apparent "emission" line near 1925 Å emerges, peaking around day 25.

### 3. Late Phase (Day 66): The Mg II Emission
* By day 66, ejecta cooling has shifted the spectral energy distribution to longer wavelengths, leaving almost no UV flux. H-alpha has developed a fully broad P Cygni profile.
* The most critical feature is the Mg II 2796, 2802 Å doublet, which appears as a broad, boxy emission line with a slanted top.

### 4. Kinematics and Modeling of the Mg II Line
* **Origin:** Comparing the velocity spreads of Mg II and H-alpha confirms that the Mg II emission originates exclusively from the "cool dense shell" forming at the very outer edge of the ejecta. Since it lacks underlying UV continuum, Mg II appears purely in emission.
* **Profile Shape:** The boxy, slanted-top shape of the Mg II line is analyzed using various models. It is best described as an asymmetric shell with an inner velocity of -6500 km/s.
* **Optical Depth Effects:** The asymmetry (where the blue side is brighter than the red side) is likely not due to an inherently lopsided explosion, but rather an optical depth effect—more photons escape from the approaching (blue) side of the shell than the receding (red) side, which is attenuated by the intervening ejecta.

# Review of Section 4: Comparison with Literature Models

This section compares the observational data of SN 2023ixf to preexisting 1D spherically symmetric radiative-transfer models to estimate the energy generated by the ejecta crashing into the CSM.

### 1. The Model Setup
* **Baseline:** The models start with a 15 solar mass progenitor at solar metallicity, evolved and exploded in simulations, then mapped to code at day 10.
* **The Cool Dense Shell:** To simulate interaction, 0.1 solar masses of material is placed at $11,700 \text{ km/s}$ to act as the cool dense shell.
* **Power Deposition:** Shock power is continuously deposited into this shell at various constant rates, ranging from 0 (no interaction) up to $1 \times 10^{43} \text{ erg/s}$.

### 2. Caveats & Inherent Discrepancies
* The authors stress these models were not custom-built for SN 2023ixf.
* **Velocity Mismatch:** Because SN 2023ixf had a dense, confined CSM, its outer ejecta was slowed down to approx. $7200 \text{ km/s}$. The models lack this dense early CSM, resulting in model velocities that are too fast. Consequently, features in the models appear broader and more blueshifted than in the actual observations.
* **Temperature Mismatch:** In the early epochs (days 14 and 19), the models are too faint in the near-UV and too bright in the optical. The authors suggest the model temperatures in the emitting region are too low.

### 3. Key Finding: Decreasing Interaction Power
* Despite the caveats, the models successfully reproduce many observed features without requiring ad-hoc scaling.
* By matching the observed flux to the models across the timeline, the authors find a progressive drop in interaction power over time.
* **Physical Implication:** This decreasing power indicates that the ejecta is sweeping through a CSM that becomes less dense at larger radii, translating to a drop in the progenitor's mass-loss rate the further you get from the star.

### 4. Decomposing the UV Spectrum (The Iron Forest)
* Because the UV spectrum is a blended "forest" of absorption lines, the authors ran modified models that intentionally omitted specific atomic transitions to isolate which elements were doing the absorbing.
* **Heavy Metal Dominance:** They found the UV absorption is overwhelmingly driven by iron-group elements—specifically Fe and Ni on day 14, with Cr joining strongly by day 19.
* **The 1924 Å Illusion:** This decomposition proved that the prominent feature at 1924 Å, which looks like an emission line, is actually a "window" of minimal absorption squeezed between two massive iron absorption troughs.

### 5. Day 66: Confirming Late-Stage Interaction
* The models dictate that the broad Mg II emission feature only appears after day 50 if there is ongoing CSM interaction. Its presence in SN 2023ixf definitively proves interaction is still happening at day 66.
* **Energy Estimate:** The day 66 flux sits cleanly between the $1 \times 10^{40}$ and $1 \times 10^{41} \text{ erg/s}$ shock power models.
* **Profile Shape:** While the model Mg II profile is a bit boxier and more blueshifted than the actual observation, the underlying physics match: the shape is caused by a shell of material where the emission from the receding (red) side is obscured by the inner ejecta.

# Review of Section 5: Comparison with Other Supernovae

This section contextualizes SN 2023ixf by comparing its UV spectra during the plateau phase to a small historical sample of Type II supernovae with available UV data.

### 1. Day ~10: UV Brightness as a CSM Indicator
* **High UV Flux:** Excess UV flux is a primary signature of CSM interaction. On day 10, SN 2023ixf is among the brightest in the sample at near-UV wavelengths.
* **Absorption Features:** Several supernovae in the sample share a distinct Fe absorption complex between 1750 Å and 2050 Å.
* **Missing Features as Clues:** The models dictate that an Fe absorption complex redward of 3000 Å only appears when CSM interaction is low. The fact that this feature is missing in SN 2023ixf supports the conclusion that it is experiencing higher levels of interaction.

### 2. Day ~20: Flux Evolution and the 1924 Å Window
* **Relative Fading:** By day 20, UV flux decreases across all supernovae in the sample. However, SN 2023ixf remains exceptionally UV-bright for its age compared to its peers.
* **The "Emission" Window:** The false emission line at 1924 Å clearly emerges in multiple supernovae in the sample at this phase.
* **Model Limitations:** The generic models struggle to perfectly reproduce the specific slope and continuum of SN 2023ixf at this epoch, highlighting the need for customized, time-dependent models.

### 3. The Ubiquity of Late-Time Mg II Emission
* The broad, boxy, blueshifted Mg II emission feature is a known hallmark of older, strongly interacting supernovae. It is confirmed in SN 2023ixf at day 66. It is notably absent in supernovae with blue supergiant progenitors (like SN 1987A) due to a lack of UV flux.

### 4. Key Conclusion: The Cause of Mg II Asymmetry
* By comparing the Mg II profiles across the historical sample, a universal pattern emerges:
    * If the profile evolves over time, it starts out symmetric and boxy, and then transitions to asymmetric (where the blue side is brighter than the red side).
    * For supernovae that always show an asymmetric profile, the asymmetry is always higher on the blue side.
* **Physical Implication:** Because the asymmetry universally favors the blue (approaching) side, it cannot be the result of a lopsided CSM viewed from random angles. Instead, it proves that the asymmetry is a universal optical depth effect—the red (receding) side of the shell is consistently blocked by the intervening ejecta, while photons from the blue side escape freely.

# Review of Section 6: Discussion

This section synthesizes the observations to test theoretical predictions, proposes a custom model to resolve discrepancies, and reconstructs the mass-loss history of the progenitor star.

### 1. Testing the Model Predictions
The time-series observations allowed the authors to evaluate key predictions for supernovae interacting with a low-density CSM:
* **Prediction 1 (Higher Luminosity):** Confirmed. The UV luminosity is consistently higher and matches interaction models at all epochs.
* **Prediction 2 (Late-Time Mg II Emission):** Confirmed. The models predict broad, boxy Mg II emission appearing ~50 days post-explosion. SN 2023ixf clearly shows this at day 66.
* **Prediction 3 (Weakened Absorption):** Confirmed. Photospheric absorption components (like in H-alpha) are very shallow and slow to develop.
* **Prediction 4 (Sharp Red Shelf in H-alpha):** Inconclusive. The signal-to-noise ratio in the data is too low to confirm.
* **Prediction 5 (Strong Narrow Absorption at Shell Velocity):** Failed/Weaker than expected. The observed absorption is significantly weaker than predicted, possibly due to a lower CSM density or asymmetrical geometry.

### 2. A Custom Model for SN 2023ixf
Because the generic models had velocities that were too fast and temperatures that were too low, the authors generated a custom model for the day 14 observation:
* **Velocity Adjustment:** They lowered the shell velocity to $8500 \text{ km/s}$.
* **X-Ray Reprocessing:** Instead of directly injecting power, they applied a uniform X-ray field (total power $5 \times 10^{42} \text{ erg/s}$) and added higher-energy ions to simulate shock physics.
* **Clumping:** They increased the volume filling factor (clumping) outside of $8000 \text{ km/s}$ from 1% to 10%.
* **Result:** This custom model successfully matched the velocities and the overall spectral energy distribution much better, proving the underlying interaction physics are correct.

### 3. Progenitor Mass-Loss History & Timeline
Using the models, the authors traced the density of the CSM backward in time (assuming a shock velocity of $10,000 \text{ km/s}$ and a wind velocity of $55 \text{ km/s}$):
* **Day 14 (Lookback: ~7 years prior):** Density was $5.2 \times 10^{-16} \text{ g/cm}^3$, equating to a high mass-loss rate of approx. $8.7 \times 10^{-4}$ solar masses per year.
* **Day 66 (Lookback: ~33-35 years prior):** Density dropped drastically to $4.8 \times 10^{-20} \text{ g/cm}^3$, which aligns with a standard, quiescent red supergiant mass-loss rate of $1.7 \times 10^{-6}$ solar masses per year.
* **Conclusion:** The star had a normal wind until about 35 years before it died. Over its final decades, the mass-loss rate ramped up dramatically before core collapse.

# Review of Section 7: Summary

This section briefly reiterates the core findings of the paper as a consolidated conclusion.

* **UV as a Critical Diagnostic:** SN 2023ixf remained exceptionally UV bright through day 24, with its spectrum dominated by iron, nickel, magnesium, and chromium absorption. UV observations are highly sensitive to weak CSM interaction that optical spectra might miss.
* **Late-Time Interaction Validation:** The broad, boxy, asymmetric Mg II emission on day 66 definitively proves that interaction with the CSM was still ongoing months after the explosion, long after the early "flash features" faded.
* **The Big Picture of SN 2023ixf's Progenitor:** While early flash spectroscopy probes the immediate, highly confined CSM (the final few years of life), these UV spectra probe the lower-density CSM further out (up to a radius of $\sim 5.7 \times 10^{15} \text{ cm}$). Combined, the data paint a clear picture of a red supergiant that transitioned from a quiescent state ($\sim 10^{-6}$ solar masses per year) to extreme, eruptive mass loss ($\sim 10^{-2}$ solar masses per year) over the last ~40 years of its life.