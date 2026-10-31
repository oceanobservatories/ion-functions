# Nitrate Functions

## Background

The Ocean Observatories Initiative deploys the Sea-Bird Scientific SUNA V2
(Submersible Ultraviolet Nitrate Analyzer) as its nitrate sensor. The SUNA V2
is an optical, chemical-free sensor derived from the ISUS (In Situ Ultraviolet
Spectroscopy) technology developed at MBARI. Prior to 2017, OOI deployed the ISUS
sensor on fixed depth platforms in the uncabled arrays. These were all converted
to the SUNA V2. The OOI instrument class NUTNR (Nutrient Sensor) covers both the 
original ISUS instruments and the SUNA V2.

`nit_functions.py` computes a single L2 data product, NITRTSC_L2, from raw
UV absorption spectra produced by NUTNR instruments. Within the OOI data
system, NIT functions fall under the Water Column science regime and the
Nitrate category.

### Primary Sources

No OOI Data Product Specification exists for the NUTNR data products. The algorithm 
documented here is derived from the source code, its embedded comments, and the two 
journal articles cited below.

---

### NITRTSC_L2 -- Temperature and Salinity Corrected Dissolved Nitrate

NITRTSC_L2 is the dissolved nitrate concentration ($\mu$M) corrected for
the temperature and salinity dependence of the seawater UV absorption
spectrum. It is produced by the NUTNR instrument class using the
Sakamoto et al. (2009) algorithm, with the temperature correction of the
seawater extinction coefficients updated to Plant et al. (2023).

#### L0 Data and Frame Types

The SUNA V2 produces two frame types in each measurement cycle:

- **Light frames** (frame type `SLB`) -- UV absorption spectra measured
  through the seawater sample. These are used to compute NITRTSC_L2.
- **Dark frames** (frame type `SDB`, `SDF`, or `NDF`) -- baseline
  measurements with the UV lamp off. These are not used by OOI; they are 
  filled with `NaN` on output.

Each frame contains a 256-element UV absorption spectrum (NITROPT_L0)
covering approximately 190-400 nm, sampled at roughly 0.8 nm intervals.
The instrument also provides a dark current scalar (`dark_value`) averaged
from the previous dark frame measurements.

#### Calibration Coefficients

Four sets of 256-element wavelength-dependent arrays are supplied from
factory calibration sheets:

- `wl` -- wavelength bins (nm) for each spectral channel
- `eno3` -- nitrate molar absorptivity (extinction coefficients) at each
  wavelength, determined from laboratory standards
- `eswa` -- seawater extinction coefficients at each wavelength, at a
  reference salinity of 35 and the calibration temperature `cal_temp`
- `di` -- deionized water reference spectrum, measured in the laboratory
  at the time of calibration

A scalar calibration temperature `cal_temp` ($^\circ$C) specifies the
temperature at which `eswa` was determined. All calibration coefficients
are time-vectorized by the OOI data management system: each coefficient
array is tiled to shape $(N, 256)$ where $N$ is the number of data packets,
allowing coefficients to vary across deployments within a single input array.

Two wavelength limit parameters bound the spectral window used in the
regression:

- `wllower` -- lower wavelength limit (nm); default 217 nm for the
  1-cm pathlength probe tip or 220 nm for the 4-cm pathlength probe tip
- `wlupper` -- upper wavelength limit (nm); default 240 nm for the
  1-cm pathlength probe tip or 245 nm for the 4-cm pathlength probe tip

#### NITRTSC_L2 Algorithm

The algorithm follows Sakamoto et al. (2009), with the temperature correction in Step 4 from Plant et al. (2023).
For each light frame $i$, the computation proceeds as follows.

**Step 1 -- Select the spectral fitting window.** Wavelength bins and their
corresponding calibration arrays are subset to channels where
$\text{wllower} \leq \lambda \leq \text{wlupper}$:

$$\lambda,\ ENO3,\ ESWA,\ DI \gets \text{channels where } \text{wllower} \leq wl \leq \text{wlupper}$$

**Step 2 -- Dark-correct the measured spectrum** by subtracting the dark
current scalar from the raw spectral counts:

$$SW_{corr} = SW - dark_value$$

**Step 3 -- Compute absorbance** from the ratio of the deionized water
reference to the dark-corrected seawater spectrum:

$$A = \log_{10}\!\left(\frac{DI}{SW_{corr}}\right)$$

**Step 4 -- Temperature-correct the seawater extinction coefficients.**
The seawater extinction coefficients, determined at the calibration temperature $T_{cal}$, are corrected to the in
situ temperature $T$ from the co-located CTD following Plant et al. (2023), Equations 6 and 8. The temperature
response is a fourth-order polynomial in the scaled wavelength $\lambda' = \lambda - 210$:

$$\begin{align}
f(\lambda') &= c_0 + c_1 \times \lambda' + c_2 \times {\lambda'}^2 \\
&\quad + c_3 \times {\lambda'}^3 + c_4 \times {\lambda'}^4
\end{align}$$

with the fixed coefficients from Plant et al. (2023):

$$\begin{align}
c_0 &= 1.46380 \times 10^{-2} \\
c_1 &= 1.67660 \times 10^{-3} \\
c_2 &= 2.91898 \times 10^{-5} \\
c_3 &= -7.56395 \times 10^{-6} \\
c_4 &= 1.27353 \times 10^{-7}
\end{align}$$

$$ESWA_T = ESWA \times \exp\!\left(f(\lambda') \times (T - T_{cal})\right)$$

Plant et al. (2023) developed this model to replace the temperature correction of Sakamoto et al. (2009), Equation
4, which can bias nitrate high by 1-2 $\mu$mol kg$^{-1}$ in waters warmer than the calibration temperature; the
updated model reduces this bias by at least 50%. The pressure correction in Plant et al. (2023), Equation 9, is not
applied.

**Step 5 -- Subtract the seawater bromide contribution** from the measured
absorbance. The seawater absorbance is proportional to practical salinity
$S_P$:

$$A_{swa} = S_P \times ESWA_T$$

$$A_{comp} = A - A_{swa}$$

**Step 6 -- Solve for nitrate concentration** by linear least squares. The
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
Light frames (frame type `SLB`) are processed through the algorithm described in the Background.

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
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |
| 2026-10-27 | Christopher Wingard | Replaced the Sakamoto et al. (2009) temperature correction of the seawater extinction coefficients with Plant et al. (2023) |

---

## References

Johnson, K. S. and Coletti, L. J. (2002). In situ ultraviolet
spectrophotometry for high resolution and long-term monitoring of nitrate,
bromide and bisulfide in the ocean. Deep-Sea Research I, 49, 1291-1305.

OOI (2014). Data Product Specification for NUTNR Data Products. Document
Control Number 1341-00620. (Not released.)

Plant, J. N., Sakamoto, C. M., Johnson, K. S., Maurer, T. L., and Bif, M. B. (2023). Updated temperature correction
for computing seawater nitrate with in situ ultraviolet spectrophotometer and submersible ultraviolet nitrate analyzer
nitrate sensors. Limnology and Oceanography: Methods, 21, 581-593.

Sakamoto, C. M., Johnson, K. S., and Coletti, L. J. (2009). Improved
algorithm for the computation of nitrate concentrations in seawater using
an in situ ultraviolet spectrophotometer. Limnology and Oceanography:
Methods, 7, 132-143.
