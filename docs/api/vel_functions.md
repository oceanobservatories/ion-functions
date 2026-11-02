# VEL Functions

## Background

The Ocean Observatories Initiative deploys two families of point velocity
instruments across its moored and profiling platforms: **VELPTMN** (Mean
Point Water Velocity) and **VELPTTU** (Turbulent Point Water Velocity). Both
families produce three-component velocity data products -- eastward (VLE),
northward (VLN), and upward (VLU) -- referenced to true Earth coordinates
after correction for magnetic declination using the International Geomagnetic
Reference Field (IGRF-14). The table below lists the OOI instrument classes
covered by this module.

| Class | Hardware | Platform | Designator Meaning |
|---|---|---|---|
| VELPT | Nortek Aquadopp | Moored (fixed depth) | Velocity Point |
| VEL3D-B | Nobska MAVS-4 | Moored (fixed depth) | Velocity 3D |
| VEL3D-C/D | Nortek Vector | Moored (fixed depth) | Velocity 3D |
| VEL3D-A | FSI ACM-3D-MP (RSN) | Cabled deep profiler | Velocity 3D |
| VEL3D-K | Nortek Aquadopp II (MMP) | Wire following profiler | Velocity 3D |
| VEL3D-L | FSI ACM-Plus (Scripps) | Global wire following profiler | Velocity 3D |

`vel_functions.py` converts L0 raw velocity data from each instrument class
into calibrated L1 velocity data products in m/s, referenced to true North
and true East. All magnetic declination corrections are computed from the IGRF
using the deployment location (latitude, longitude, depth) and NTP timestamp.

Within the OOI data system, VEL data products fall under the Water Column [sampling
regime](https://oceanobservatories.org/data-products/) (see the [complete list of OOI data
products](https://oceanobservatories.org/data-products/data-products-complete-list/)).

Apart from the VEL3D-A compass calibration, no calibration coefficients are applied by these functions; the deployment
latitude and longitude are used for the magnetic declination correction. VELPT and VEL3D-B/C/D/L instruments report
velocities in instrument earth coordinates (magnetic); VEL3D-K beam velocities are transformed to earth coordinates by
ion-functions using the fixed instrument and profiler geometry.

For where to find the calibration coefficients and raw inputs for an OOI dataset, and how the variables in an OOI data
file map to these functions' arguments, see [Calibration Coefficients](../calibrations.md) and [Using These Functions
with OOI Data](../using-ooi-data.md).

### Primary Sources

| DCN | Document                                                                                                                                                                                             |
|---|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1341-00780 | [OOI (2012). Data Product Specification for Turbulent Point Water Velocity (Nortek Vector).](https://oceanobservatories.org/wp-content/uploads/2023/09/Data_Product_SPEC_VELPTTU_combined.pdf)       |
| 1341-00781 | [OOI (2012). Data Product Specification for Turbulent Point Water Velocity (Nobska MAVS-4).](https://oceanobservatories.org/wp-content/uploads/2023/09/Data_Product_SPEC_VELPTTU_combined.pdf)       |
| 1341-00782 | [OOI (2012). Data Product Specification for Turbulent Point Water Velocity (ACM and Aquadopp II).](https://oceanobservatories.org/wp-content/uploads/2023/09/Data_Product_SPEC_VELPTTU_combined.pdf) |
| 1341-00790 | [OOI (2013). Data Product Specification for Mean Point Water Velocity.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00790_Data_Product_SPEC_VELPTMN_OOI.pdf)                      |
| 1341-00792 | OOI (2014). Data Product Specification for Mean Point Water Velocity Data from FSI Acoustic Current Meters. |

---

### VELPTMN_L1 -- Mean Point Water Velocity

VELPTMN_L1 is the mean point water velocity (m/s) measured by moored and
profiling current meters that average over an ensemble interval. The product
comprises three components: eastward (VLE), northward (VLN), and upward
(VLU), all referenced to true Earth coordinates.

**VELPT-A/B/D/J (Nortek Aquadopp):** The Aquadopp outputs velocity components
already converted to earth coordinates (East, North, Up) onboard using its
internal heading, pitch, and roll sensors. The L0 data are in mm/s. The L1
product is computed by `velpt_mag_corr_east` and `velpt_mag_corr_north`, which
call `vel_mag_correction` to apply the IGRF magnetic declination rotation and
convert units to m/s. The upward component requires only a unit conversion via
`velpt_up_vel`.

**VEL3D-A (FSI ACM-3D-MP, RSN cabled deep profiler):** The FSI ACM stinger
configuration has four acoustic transducers at the tips of four prongs
oriented 45$^\circ$ from a central post. Four raw beam velocities (vp1-vp4)
are produced; vp1 and vp3 lie in the horizontal plane and vp2 and vp4 lie in
the vertical plane. The core function `fsi_acm_horz_vel` converts the
horizontal beam pair to eastward and northward velocity in m/s using the
beam-to-instrument transform:

$$V_x = \frac{-\text{vp1} - \text{vp3}}{\sqrt{2}}$$

$$V_y = \frac{\text{vp1} - \text{vp3}}{\sqrt{2}}$$

The speed and cartesian heading are then combined with the nautical heading
and the IGRF magnetic declination to produce true-Earth-referenced components.
Because pitch and roll are negligible for profiler-mounted instruments, no
tilt correction is applied.

For RSN deployments (VEL3D-A), the wrapper functions `fsi_acm_rsn_east` and
`fsi_acm_rsn_north` first derive the instrument heading from field direction
cosines (hx, hy) corrected by an 8-point spin test calibration (via
`fsi_acm_nautical_heading` and `fsi_acm_compass_cal`), then call
`fsi_acm_horz_vel`. The calibration computes offsets, scaling factors, and a
compass bias following the procedure in the McLane Moored Profiler (MMP)
manual Appendix D.

#### Calibration Coefficients

| Function Argument | Calibration Coefficient | Source |
|---|---|---|
| `hdg_cal` | `CC_hdg_cal` | Compass calibration spin test (MMP manual, Appendix D) |
| `hx_cal` | `CC_hx_cal` | Compass calibration spin test (MMP manual, Appendix D) |
| `hy_cal` | `CC_hy_cal` | Compass calibration spin test (MMP manual, Appendix D) |

OOI calibration coefficient names are explained in [Using These Functions with OOI Data](../using-ooi-data.md).

For Scripps deployments (VEL3D-L), the wrapper functions `fsi_acm_sio_east`
and `fsi_acm_sio_north` pass the instrument-reported heading directly to
`fsi_acm_horz_vel`; no calibration correction is needed.

**Vertical velocity (VEL3D-A and VEL3D-L):** Two vertical velocity products
are computed depending on profiler travel direction, to avoid contamination
from the sheet-flow wake of the stinger's central post. During ascent (vp2
contaminated), `fsi_acm_up_profiler_ascending` uses vp1, vp3, and vp4:

$$w_{asc} = \frac{1}{100}\left(-\frac{\text{vp1}+\text{vp3}}{\sqrt{2}}
            - \sqrt{2}\,\text{vp4}\right)$$

During descent (vp4 contaminated), `fsi_acm_up_profiler_descending` uses
vp1, vp2, and vp3:

$$w_{dsc} = \frac{1}{100}\left(\frac{\text{vp1}+\text{vp3}}{\sqrt{2}}
            + \sqrt{2}\,\text{vp2}\right)$$

The factor $1/100$ converts from cm/s to m/s. The MMP manual (Rev E,
p. 8-30) contains errors in both upward velocity formulas; the expressions
above reflect the corrected implementations in the code.

---

### VELPTTU_L1 -- Turbulent Point Water Velocity

VELPTTU_L1 is the turbulent point water velocity (m/s) comprising eastward
(VLE), northward (VLN), and upward (VLU) components in true Earth coordinates.

**VEL3D-B (Nobska MAVS-4):** The MAVS-4 outputs velocity in earth coordinates
in cm/s. The wrapper functions `nobska_mag_corr_east` and
`nobska_mag_corr_north` call `vel_mag_correction` to apply the IGRF rotation
and convert to m/s. The upward component requires only a unit conversion via
`nobska_scale_up_vel`.

**VEL3D-C/D (Nortek Vector):** The Vector outputs velocity in mm/s. The
wrapper functions `nortek_mag_corr_east` and `nortek_mag_corr_north` convert
units and call `vel_mag_correction`. An additional scaling factor of 0.1 is
applied when bit 1 of the instrument status code is set, via the helper
`nortek_scale_velocity`. The upward component is handled by `nortek_up_vel`.

**VEL3D-K (Nortek Aquadopp II on McLane WFP):** The Aquadopp II provides
individual beam velocities that require a full beam-to-ENU coordinate
transformation via the core function `vel3dk_transform`. The transformation
proceeds in two stages:

*Stage 1 -- Beam to profiler XYZ:* Beam velocities are converted to Cartesian
coordinates in the profiler frame using a pre-computed transformation matrix
$T_{beam2XYZ}$ selected by `get_XYZ_transform`. Three beam configurations are
supported: upward-traveling profiles (beams 1, 2, 4), downward-traveling
profiles (beams 2, 3, 4), and stationary measurements (all 4 beams). The
matrices are computed once at module import by `generate_beam_transforms`.

The azimuthal angles of the four beams increase in the counter-clockwise
direction (1 $\to$ 4 $\to$ 3 $\to$ 2), with elevation angles of 47.5$^\circ$
(beams 1 and 3) and 25$^\circ$ (beams 2 and 4) from vertical. The stationary
(4-beam) transform has not been exercised in OOI deployments and is
unverified.

*Stage 2 -- Profiler XYZ to Earth ENU:* The rotation matrix $T_{XYZ2ENU}$ is
constructed by `generate_ENU_transform` from the instrument attitude at each
sample:

$$T_{XYZ2ENU} = R_z \cdot R_y \cdot R_x$$

where $R_x$, $R_y$, $R_z$ are the un-roll, un-pitch, and un-heading rotation
matrices. Because the Aquadopp II measures heading as the direction of the
instrument X-axis relative to magnetic north, 90$^\circ$ is subtracted from
heading before constructing $R_z$. Attitude variables in the raw binary data
file are in deci-degrees (0.1$^\circ$) and are converted to degrees before
the rotation is applied.

The wrapper functions `vel3dk_east`, `vel3dk_north`, and `vel3dk_up` each call
`vel3dk_transform` and extract the appropriate component. Magnetic declination
correction via `vel_mag_correction` is applied to the horizontal components.

In 2017 the VEL3D-K velocities were found to be aliased by an ambiguity-velocity problem. In diagnosing it, data were
sent to Nortek, who returned corrected data and code that differed significantly from the 2014 code used for the
original implementation; the first set returned had a sign error that produced downward 'current' velocities during
downward profiles. The 2018 revision matches Nortek's second set, and was independently derived from the transducer
geometry, the radial beam-velocity sign convention, and the X-axis heading reference.

The 2018 revision to `generate_beam_transforms` and `generate_ENU_transform`
corrected three errors in the original 2014 Nortek-supplied code: (1) the
transducer azimuthal angles were defined clockwise, producing a left-handed
coordinate system; (2) pitch values were incorrectly negated; (3) heading was
not offset by 90$^\circ$ to account for the Aquadopp II heading reference
being the instrument X-axis rather than the Y-axis.

---

### Magnetic Declination Correction

The core function `vel_mag_correction` applies the IGRF rotation to any
horizontal velocity pair:

$$\begin{align}
u_{cor} &= u \cos\theta + v \sin\theta \\
v_{cor} &= v \cos\theta - u \sin\theta
\end{align}$$

where $\theta$ is the magnetic declination (degrees, positive east) computed
from `generic_functions.magnetic_declination` using the deployment latitude,
longitude, depth, and NTP timestamp. All instrument-specific wrapper functions
that produce horizontal L1 products call `vel_mag_correction`.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.vel_functions.vel_mag_correction

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2013-04-24 | Stuart Pearce | Generalized for all velocity instruments. |
| 2014-02-05 | Christopher Wingard | Refactored to use generic_functions module. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.fsi_acm_horz_vel

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-18 | Russell Desiderio | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.fsi_acm_up_profiler_ascending

#### Additional Notes
The MMP manual (Rev E, p. 8-30) contains an error in the formula for upward
velocity during ascent. The formula implemented here follows the corrected
algorithm derived independently of the manual.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-13 | Russell Desiderio | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.fsi_acm_up_profiler_descending

#### Additional Notes
The MMP manual (Rev E, p. 8-30) contains an error in the formula for upward
velocity during descent. The formula implemented here follows the corrected
algorithm derived independently of the manual.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-13 | Russell Desiderio | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.nobska_scale_up_vel

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.nortek_up_vel

#### Additional Notes
As of 2015-06-08, the DPS (1341-00780) incorrectly states the Nortek Vector
outputs velocities in m/s. The instrument outputs velocities in mm/s; the
code implements the correct unit conversion.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2015-06-08 | Russell Desiderio | Corrected input units to mm/s. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.velpt_up_vel

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.vel3dk_transform

#### Additional Notes
The stationary (4-beam) transform case has not been exercised in OOI McLane
profiler deployments and should be treated as unverified per developer
warning in the code (R. Desiderio, 2018-06-20).

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-13 | Stuart Pearce | Reimplementation of Nortek (Rusello) code. |
| 2015-06-02 | Russell Desiderio | Trap fill values in beam config; return NaNs for those records. |
| 2018-06-19 | Russell Desiderio | Documented that outputs are ENU uncorrected for magnetic declination. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

## Helper Functions

::: ion_functions.data.vel_functions.fsi_acm_nautical_heading

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-18 | Russell Desiderio | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.fsi_acm_compass_cal

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-16 | Russell Desiderio | Initial code. |
| 2015-06-01 | Russell Desiderio | Adjusted to process vertically stacked input variables. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.generate_beam_transforms

#### Additional Notes
The stationary transform (all 4 beams) has not been used in OOI McLane
profiler deployments and is unverified per developer warning (R. Desiderio,
2018-06-20). The 2018 revision corrected three errors in the original 2014
Nortek-supplied code; see the Background section for details.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-13 | Stuart Pearce | Reimplementation of Nortek (Rusello) code. |
| 2018-06-18 | Russell Desiderio | Corrected transducer azimuth direction, pitch sign, and heading offset. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.get_XYZ_transform

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-13 | Stuart Pearce | Reimplementation of Nortek (Rusello) code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.generate_ENU_transform

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-13 | Stuart Pearce | Reimplementation of Nortek (Rusello) code. |
| 2018-06-18 | Russell Desiderio | Corrected heading offset and pitch sign. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.nortek_scale_velocity

#### History
| Date | Author | Change |
|---|---|---|
| 2015-06-08 | Russell Desiderio | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.valid_lat

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.valid_lon

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

## Wrapper Functions

::: ion_functions.data.vel_functions.fsi_acm_rsn_east

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-13 | Russell Desiderio | Initial code. |
| 2015-05-29 | Russell Desiderio | Time-vectorized calcoeffs. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.fsi_acm_rsn_north

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-13 | Russell Desiderio | Initial code. |
| 2015-05-29 | Russell Desiderio | Time-vectorized calcoeffs. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.fsi_acm_sio_east

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-18 | Russell Desiderio | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.fsi_acm_sio_north

#### History
| Date | Author | Change |
|---|---|---|
| 2015-02-18 | Russell Desiderio | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.nobska_mag_corr_east

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.nobska_mag_corr_north

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.nortek_mag_corr_east

#### Additional Notes
As of 2015-06-08, the DPS (1341-00780) incorrectly states the Nortek Vector
outputs velocities in m/s. The instrument outputs velocities in mm/s; the
code implements the correct unit conversion.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2015-06-08 | Russell Desiderio | Corrected input units to mm/s. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.nortek_mag_corr_north

#### Additional Notes
As of 2015-06-08, the DPS (1341-00780) incorrectly states the Nortek Vector
outputs velocities in m/s. The instrument outputs velocities in mm/s; the
code implements the correct unit conversion.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2015-06-08 | Russell Desiderio | Corrected input units to mm/s. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.velpt_mag_corr_east

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.velpt_mag_corr_north

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-17 | Stuart Pearce | Initial code. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.vel3dk_east

#### Additional Notes
The stationary (4-beam) transform case has not been exercised in OOI McLane
profiler deployments and should be treated as unverified per developer
warning in the code (R. Desiderio, 2018-06-20).

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-13 | Stuart Pearce | Initial code (Rusello implementation). |
| 2015-06-03 | Russell Desiderio | Fill value handling for beam config. |
| 2018-06-18 | Russell Desiderio | Corrected transform errors; added documentation. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.vel3dk_north

#### Additional Notes
See Additional Notes for `vel3dk_east`; the same transform corrections and
caveats apply.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-13 | Stuart Pearce | Initial code (Rusello implementation). |
| 2015-06-03 | Russell Desiderio | Fill value handling for beam config. |
| 2018-06-18 | Russell Desiderio | Corrected transform errors; added documentation. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.vel_functions.vel3dk_up

#### Additional Notes
See Additional Notes for `vel3dk_east`; the same transform corrections and
caveats about the stationary case apply.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-13 | Stuart Pearce | Initial code (Rusello implementation). |
| 2018-06-18 | Russell Desiderio | Corrected transform errors; added documentation. |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

## References

Beggan, C.D., et al. (2026). International Geomagnetic Reference Field: the fourteenth
generation. *Earth Planets Space*, 78, 127. <https://doi.org/10.1186/s40623-025-02360-0>

International Association of Geomagnetism and Aeronomy (2024). [IGRF-14.
Zenodo.](https://doi.org/10.5281/zenodo.14012302)

[OOI (2012). Data Product Specification for Turbulent Point Water Velocity
(Nobska MAVS-4). Document Control Number 1341-00781.](https://oceanobservatories.org/wp-content/uploads/2023/09/Data_Product_SPEC_VELPTTU_combined.pdf)

[OOI (2012). Data Product Specification for Turbulent Point Water Velocity
(Nortek Vector). Document Control Number 1341-00780.](https://oceanobservatories.org/wp-content/uploads/2023/09/Data_Product_SPEC_VELPTTU_combined.pdf)

[OOI (2013). Data Product Specification for Mean Point Water Velocity. Document
Control Number 1341-00790.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00790_Data_Product_SPEC_VELPTMN_OOI.pdf)

OOI (2014). Data Product Specification for Mean Point Water Velocity Data from
FSI Acoustic Current Meters. Document Control Number 1341-00792.

Strom, K.M., and Reistad, H. (2024). [ppigrf: Python package for computing
the International Geomagnetic Reference Field
(IGRF).](https://github.com/IAGA-VMOD/ppigrf)
