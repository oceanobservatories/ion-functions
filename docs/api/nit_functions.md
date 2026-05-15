# Nitrate Functions

## Background

The Ocean Observatories Initiative deploys the Sea-Bird Scientific SUNA V2
(Submersible Ultraviolet Nitrate Analyzer) as its nitrate sensor. The SUNA V2
is an optical, chemical-free sensor derived from the ISUS (In Situ Ultraviolet
Spectroscopy) technology developed at MBARI. OOI instrument class NUTNR
(Nutrient Sensor) covers both the original ISUS instruments and the SUNA V2.

`nit_functions.py` computes a single L2 data product, NITRTSC_L2, from raw
UV absorption spectra produced by NUTNR instruments.

### Primary Sources

No OOI Data Product Specification exists for the NUTNR data products. The algorithm 
documented here is derived from the source code, its embedded comments, and the two 
journal articles cited in the code.

---

### NITRTSC_L2 — Temperature and Salinity Corrected Dissolved Nitrate

NITRTSC_L2 is the dissolved nitrate concentration ($\mu$M) corrected for
the temperature and salinity dependence of the seawater UV absorption
spectrum. It is produced by the NUTNR instrument class using the
Sakamoto et al. (2009) algorithm.

#### L0 Data and Frame Types

The SUNA V2 produces two frame types in each measurement cycle:

- **Light frames** (frame type `SLB`) — UV absorption spectra measured
  through the seawater sample. These are used to compute NITRTSC_L2.
- **Dark frames** (frame type `SDB`, `SDF`, or `NDF`) — baseline
  measurements with the UV lamp off. These are used to correct light
  frames for electronic noise and are not converted to nitrate
  concentrations; they are filled with `NaN` on output.

Each frame contains a 256-element UV absorption spectrum (NITROPT_L0)
covering approximately 190–400 nm, sampled at roughly 0.8 nm intervals.
The instrument also provides a dark current scalar (`dark_value`) averaged
from dark frame measurements.

#### Calibration Coefficients

Four sets of 256-element wavelength-dependent arrays are supplied from
factory calibration sheets:

- `wl` — wavelength bins (nm) for each spectral channel
- `eno3` — nitrate molar absorptivity (extinction coefficients) at each
  wavelength, determined from laboratory standards
- `eswa` — seawater extinction coefficients at each wavelength, at a
  reference salinity of 35 and the calibration temperature `cal_temp`
- `di` — deionized water reference spectrum, measured in the laboratory
  at the time of calibration

A scalar calibration temperature `cal_temp` ($^\circ$C) specifies the
temperature at which `eswa` was determined. All calibration coefficients
are time-vectorized by the OOI data management system: each coefficient
array is tiled to shape $(N, 256)$ where $N$ is the number of data packets,
allowing coefficients to vary across deployments within a single input array.

Two wavelength limit parameters bound the spectral window used in the
regression:

- `wllower` — lower wavelength limit (nm); default 217 nm for the
  1-cm pathlength probe tip or 220 nm for the 4-cm pathlength probe tip
- `wlupper` — upper wavelength limit (nm); default 240 nm for the
  1-cm pathlength probe tip or 245 nm for the 4-cm pathlength probe tip

#### NITRTSC_L2 Algorithm

The algorithm follows Sakamoto et al. (2009). For each light frame $i$,
the computation proceeds as follows.

**Step 1 — Select the spectral fitting window.** Wavelength bins and their
corresponding calibration arrays are subset to channels where
$\text{wllower} \leq \lambda \leq \text{wlupper}$:

$$\lambda,\ ENO3,\ ESWA,\ DI \gets \text{channels where } \text{wllower} \leq wl \leq \text{wlupper}$$

**Step 2 — Dark-correct the measured spectrum** by subtracting the dark
current scalar from the raw spectral counts:

$$SW_{corr} = SW - dark\_value$$

**Step 3 — Compute absorbance** from the ratio of the deionized water
reference to the dark-corrected seawater spectrum:

$$A = \log_{10}\!\left(\frac{DI}{SW_{corr}}\right)$$

**Step 4 — Temperature-correct the seawater extinction coefficients.**
The in situ temperature ($T$) and calibration temperature ($T_{cal}$) are
used to scale `eswa` following Sakamoto et al. (2009), Equation 4. The
coefficients $A_{sak}$, $B_{sak}$, $C_{sak}$, $D_{sak}$ are fixed
constants from Sakamoto et al. (2009):

$$\begin{align}
A_{sak} &= 1.1500276 \\
B_{sak} &= 0.02840 \\
C_{sak} &= -0.3101349 \\
D_{sak} &= 0.001222
\end{align}$$

$$ESWA_T = ESWA \times \frac{A_{sak} + B_{sak} \times T}{A_{sak} + B_{sak} \times T_{cal}}
\times \exp\!\left(D_{sak} \times (T - T_{cal}) \times (\lambda - 210)\right)$$

**Step 5 — Subtract the seawater bromide contribution** from the measured
absorbance. The seawater absorbance is proportional to practical salinity
$S_P$:

$$A_{swa} = S_P \times ESWA_T$$

$$A_{comp} = A - A_{swa}$$

**Step 6 — Solve for nitrate concentration** by linear least squares. The
model matrix $M$ includes the nitrate extinction spectrum, a constant
baseline offset, and a linear wavelength baseline:

$$M = \begin{bmatrix} ENO3 & \mathbf{1}/100 & \lambda/1000 \end{bmatrix}$$

The solution vector $C$ is obtained via the Moore-Penrose pseudoinverse:

$$C = M^{+} \times A_{comp}$$

The first element $C[0]$ is NITRTSC_L2, the nitrate concentration in
$\mu$M. The second and third elements capture a linear baseline that
accounts for residual spectral interference.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.nit_functions.ts_corrected_nitrate

#### Additional Notes
Dark frames (frame types `SDB`, `SDF`, and `NDF`) produce `NaN` on output.
Light frames (frame type `SLB`) are processed through the full Sakamoto
et al. (2009) TS-corrected algorithm.

The `wllower` and `wlupper` arguments default to 217 nm and 240 nm
respectively, matching the 1-cm pathlength probe tip specification. For
4-cm pathlength probe tips, the OOI system supplies 220 nm and 245 nm via
the time-vectorized calibration coefficient mechanism.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-22 | Craig Risien | Initial code |
| 2014-05-27 | Craig Risien | Added light vs dark frame filtering |
| 2015-04-09 | Russell Desiderio | Revised for time-vectorized calibration coefficient arrays; changed fill value to NaN |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

Johnson, K. S. and Coletti, L. J. (2002). In situ ultraviolet
spectrophotometry for high resolution and long-term monitoring of nitrate,
bromide and bisulfide in the ocean. Deep-Sea Research I, 49, 1291-1305.

Sakamoto, C. M., Johnson, K. S., and Coletti, L. J. (2009). Improved
algorithm for the computation of nitrate concentrations in seawater using
an in situ ultraviolet spectrophotometer. Limnology and Oceanography:
Methods, 7, 132-143.

OOI (2014). Data Product Specification for NUTNR Data Products. Document
Control Number 1341-00620. (Not released.)
