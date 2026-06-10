# ADCP Functions

!!! note "Module Note"
    This page covers two instrument families that share a single module:
    **ADCPS/ADCPT/ADCPA** (Teledyne RDI Workhorse, producing VELPROF and
    ECHOINT) and **VADCP** (Teledyne RDI 5-Beam Workhorse, producing VELTURB
    and ECHOINT). The deprecated original VADCP functions are listed under
    [Deprecated Functions](#deprecated-functions).

## Background

The ADCP instrument family covers Teledyne RDI Workhorse acoustic Doppler
current profilers (ADCPs) deployed across the OOI network. Three instrument
classes produce the VELPROF and ECHOINT data products, and one class produces
the VELTURB and ECHOINT data products.

| Class | Hardware | Platform | Designator Meaning |
|---|---|---|---|
| ADCPA | TRDI Explorer PA DVL | Coastal glider | ADCP, glider-mounted |
| ADCPS | TRDI Workhorse Long Ranger (75 kHz) | Subsurface mooring | ADCP, subsurface |
| ADCPT | TRDI Workhorse Quartermaster (150 kHz) | Surface mooring | ADCP, surface |
| VADCP | TRDI Workhorse Sentinel, 5-beam (600 kHz) | RSN 200 m platform | Vertical ADCP |
| VADCP-B | Nortek Signature 55 | RSN 200 m platform | Vertical ADCP, 2nd generation |

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00750 | [OOI (2020). Data Product Specification for Velocity Profile and Echo Intensity. Document Control Number 1341-00750.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00750_Data_Product_SPEC_VELPROF_ECHOINT_OOI.pdf) |
| 1341-00760 | [OOI (2013). Data Product Specification for Turbulent Velocity Profile and Echo Intensity. Document Control Number 1341-00760.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00760_Data_Product_SPEC_VELTURB_ECHOINT_OOI.pdf) |

### VELPROF_L1 and ECHOINT_L1 -- Velocity Profile and Echo Intensity

VELPROF_L1 is the velocity profile (m/s) in Earth coordinates with magnetic
variation correction applied, produced by ADCPA, ADCPS, and ADCPT instruments.
ECHOINT_L1 is the echo intensity (dB) for each of the four beams.

Each Workhorse ADCP has four acoustic beams arranged symmetrically around the
instrument axis at a fixed beam angle of 20$^\circ$ from vertical. Each beam
measures the radial velocity of acoustic scatterers in the water column using
the Doppler effect. The four beams are labeled clockwise in the order 3-1-4-2
relative to the instrument frame, with the X-axis pointing from transducer 1
toward transducer 2, and the Y-axis pointing from transducer 4 toward
transducer 3 (TRDI, 2010a, p. 16-17).

Raw L0 velocity data (VELPROF-B1 through VELPROF-B4) are output in the
instrument's PD0 binary format in units of mm/s. For instruments programmed in
beam coordinates (ADCPS-I, ADCPS-K, ADCPT-B, ADCPT-D, ADCPT-E), the
processing chain is: beam to instrument coordinates, instrument to Earth
coordinates, then magnetic declination correction. For instruments programmed
in Earth coordinates (ADCPA, ADCPS-J, ADCPS-L, ADCPS-N, ADCPT-C, ADCPT-F,
ADCPT-G, ADCPT-M), only the magnetic declination correction is applied. The
final L1 products are scaled from mm/s to m/s.

**ADCPA (glider):** The ADCPA produces L1 VELPROF directly onboard; only
the magnetic declination correction is applied by ion-functions.

#### Beam to Instrument Coordinates

For beam-coordinate instruments, the beam velocity components (b1, b2, b3, b4)
are converted to instrument-frame components (x, y, z, e) using the
transformation matrix $A$ defined in DPS 1341-00750, Equation 1:

$$A = \begin{pmatrix}
c \times a & -c \times a & 0 & 0 \\
0 & 0 & -c \times a & c \times a \\
b & b & b & b \\
d & d & -d & -d
\end{pmatrix}$$

where the constants are:

$$a = \frac{1}{2 \sin\theta}, \quad
b = \frac{1}{4 \cos\theta}, \quad
c = +1\ \text{(convex transducer head)}, \quad
d = \frac{a}{\sqrt{2}}$$

and $\theta = 20^\circ$ is the fixed beam angle for both the Workhorse Long
Ranger and Workhorse Quartermaster. Applying $A$ yields:

$$x = c \times a \times (b1 - b2)$$
$$y = c \times a \times (b4 - b3)$$
$$z = b \times (b1 + b2 + b3 + b4)$$
$$e = d \times (b1 + b2 - b3 - b4)$$

#### Three-Beam Solution

When one beam fails, as determined by the percent good variable for that beam
falling below 25%, a three-beam solution is applied if the remaining three
beams exceed 25% good. The failed beam's velocity is reconstructed by forcing
the error velocity to zero, substituting the reconstructed value into the
standard transformation (DPS 1341-00750, Section 5; TRDI, 2010a, p. 14). If
more than one beam fails per depth cell, the result is set to NaN.

#### Instrument to Earth Coordinates

The instrument-frame velocities (x, y, z) are rotated to Earth coordinates
(uu, vv, ww) using the rotation matrix $M$ (DPS 1341-00750, Equation 2):

$$M = \begin{pmatrix}
\cos H & \sin H & 0 \\
-\sin H & \cos H & 0 \\
0 & 0 & 1
\end{pmatrix}
\times
\begin{pmatrix}
1 & 0 & 0 \\
0 & \cos P & -\sin P \\
0 & \sin P & \cos P
\end{pmatrix}
\times
\begin{pmatrix}
\cos R & 0 & \sin R \\
0 & 1 & 0 \\
-\sin R & 0 & \cos R
\end{pmatrix}$$

where $H$, $P$, and $R$ are the heading, pitch, and roll angles. For
upward-looking instruments, 180$^\circ$ is added to the measured roll angle
before applying $M$ (DPS 1341-00750, Equations 3-4):

$$R = \text{Tilt2} + 180^\circ \quad \text{(upward-looking)}$$
$$P = \arctan[\tan(\text{Tilt1}) \times \cos(\text{Tilt2})]$$

Heading, pitch, and roll are recorded in the PD0 variable leader in
centidegrees and converted to degrees before the rotation is applied.

#### Magnetic Declination Correction

The final processing step rotates the horizontal velocity components (uu, vv)
clockwise by the magnetic declination angle $\theta$ (DPS 1341-00750,
Equation 5):

$$\begin{pmatrix} U \\ V \end{pmatrix}
= \begin{pmatrix}
\cos\theta & \sin\theta \\
-\sin\theta & \cos\theta
\end{pmatrix}
\times
\begin{pmatrix} uu \\ vv \end{pmatrix}$$

Magnetic declination is computed using the World Magnetic Model via
`magnetic_declination` in `generic_functions.py`, using the deployment
latitude, longitude, and timestamp.

#### Echo Intensity

Raw echo intensity (ECHOINT-B1 through ECHOINT-B4) is recorded as a 1-byte
integer (counts) per beam per depth cell. ECHOINT_L1 is computed by scaling
the raw counts by a factory-supplied scale factor (nominally 0.45 dB/count
for the Workhorse family):

$$\text{ECHOINT\_L1} = \text{raw counts} \times \text{sfactor}$$

#### Bin Depths

Center depths of the ADCP measurement bins are computed from the sensor
depth, blanking distance, bin size, and instrument orientation. The sensor
depth is derived from the instrument pressure reading using the TEOS-10
`z_from_p` function. Three pressure input unit variants are supported:
bar (`adcp_bin_depths_bar`), decaPascal (`adcp_bin_depths_dapa`), and
meters directly (`adcp_bin_depths_meters`). The unified `adcp_bin_depths`
function accepts depth in meters directly.

### VELTURB_L1 and ECHOINT_L1 -- Turbulent Velocity Profile

VELTURB_L1 is the turbulent velocity profile (m/s) in Earth coordinates with
magnetic variation correction applied, produced by the VADCP and VADCP-B
instruments deployed on RSN 200 m platforms. ECHOINT_L1 is the echo intensity
(dB) for each beam.

The original VADCP (Teledyne RDI Workhorse Sentinel, 600 kHz) is a 5-beam
instrument operated as two synchronized ADCPs in a master/slave configuration.
The primary (master) 4-beam unit provides heading, tilt, and compass data; the
secondary (slave) unit provides only beam 5 (the vertical beam). The 5th beam
points directly upward and provides an independent estimate of the vertical
velocity component.

The VADCP-B (Nortek Signature 55) uses a different beam geometry, beam
numbering convention, beam velocity units (m/s directly, not mm/s), and a
vendor-supplied 4x4 transformation matrix (`tm`) in place of the fixed
analytical matrix used for TRDI instruments.

#### VELTURB -- Original VADCP (TRDI)

The beam-to-instrument and instrument-to-Earth transforms for the original
VADCP use the same matrix formulation as VELPROF (DPS 1341-00760, Equations
1-4), with the beam angle fixed at 20$^\circ$ for the 4-beam unit. The
5th beam points at 0$^\circ$ (directly vertical) and is used only to
compute the true vertical velocity product (VELTURB-VLU-5BM_L1).

Two vertical velocity products are defined:

- **VELTURB-VLU-4BM_L1** -- estimated vertical velocity derived from the
  standard 4-beam transform; the z component of the instrument-to-Earth
  rotation applied to the 4 tilted beams.
- **VELTURB-VLU-5BM_L1** -- true vertical velocity taken directly from
  beam 5 (the nadir-pointing beam), passed through the instrument-to-Earth
  rotation using the beam 5 velocity in place of the z component. Depth
  cells where beam 5 percent good falls below 25% are set to the fill value.

The deprecated `vadcp_beam_eastward`, `vadcp_beam_northward`, and
`vadcp_beam_error` functions use the standard 4-beam `adcp_beam2ins` and
`adcp_ins2earth` pipeline. `vadcp_beam_vertical_est` and
`vadcp_beam_vertical_true` are not deprecated because no VADCP-B replacements
exist for the 5-beam vertical products; see
[Deprecated Functions](#deprecated-functions).

#### VELTURB -- VADCP-B (Nortek Signature 55)

The VADCP-B beam-to-instrument transform uses the vendor-supplied
transformation matrix `tm` (a 4x4 matrix) applied via Einstein summation.
The beam velocities are already in m/s (no mm/s-to-m/s scaling required).
The transform produces four instrument-frame components (x, y, z1, z2), where
z1 and z2 are two estimates of the vertical component from beam pairs 1-3 and
2-4, respectively. This two-z output replicates the Nortek-provided MATLAB
reference code structure.

The instrument-to-Earth rotation for the VADCP-B uses a Nortek-convention
heading offset (heading minus 90$^\circ$) and constructs an extended 4x4
rotation matrix to accommodate the two z components. The final vertical
velocity is the average of the Earth-rotated z1 and z2.

The three-beam solution threshold for the VADCP-B is 50% (versus 25% for
TRDI instruments), following the vendor-specified percent good floor.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.adcp_functions.adcp_beam_velocity

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2014-02-03 | Christopher Wingard | Updated to use WMM 2010 magnetic declination |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2014-06-25 | Christopher Wingard | Corrected units for heading, pitch, roll, depth |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions; removed depth dependence from magnetic declination |
| 2019-03-11 | Christopher Wingard | Removed multiple wrapper functions; made function standalone returning all velocity components |
| 2023-08-15 | Samuel Dahlberg | Ported from Pyseas for CGSN compatibility |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_backscatter

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-21 | Christopher Wingard | Initial code |
| 2015-06-25 | Russell Desiderio | Added fill value to NaN conversion |
| 2023-08-15 | Samuel Dahlberg | Added default value to sfactor |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_beam_vertical_est

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions |
| 2015-06-22 | Russell Desiderio | Renamed data product |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_beam_vertical_true

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions |
| 2015-06-22 | Russell Desiderio | Renamed data product |
| 2015-06-25 | Russell Desiderio | Added beam 5 fill value to NaN conversion |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

### VELPROF -- Beam Coordinates

::: ion_functions.data.adcp_functions.adcp_beam_eastward

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2014-02-03 | Christopher Wingard | Updated to use WMM 2010 magnetic declination |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2014-06-25 | Christopher Wingard | Corrected units for heading, pitch, roll, depth |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions; removed depth dependence from magnetic declination |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2023-08-15 | Samuel Dahlberg | Refactored as wrapper calling adcp_beam_velocity; added NTP to Unix timestamp conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_beam_northward

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2014-02-03 | Christopher Wingard | Updated to use WMM 2010 magnetic declination |
| 2014-03-28 | Russell Desiderio | Documentation correction |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2014-06-25 | Christopher Wingard | Corrected units for heading, pitch, roll, depth |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions; removed depth dependence from magnetic declination |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2023-08-15 | Samuel Dahlberg | Refactored as wrapper calling adcp_beam_velocity; added NTP to Unix timestamp conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_beam_vertical

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2014-02-03 | Christopher Wingard | Updated to use WMM 2010 magnetic declination |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2014-06-25 | Christopher Wingard | Corrected units for heading, pitch, roll, depth |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_beam_error

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2015-06-10 | Russell Desiderio | Moved beam conditioning to helper function |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

### VELPROF -- Earth Coordinates

::: ion_functions.data.adcp_functions.adcp_earth_eastward

#### Additional Notes
The `z` parameter (instrument pressure, daPa) is accepted by this function
but not used in the computation. The depth-to-dbar conversion lines in the
function body are commented out. This is a known issue in the current
implementation.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2014-02-03 | Christopher Wingard | Updated to use WMM 2010 magnetic declination |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2014-06-25 | Christopher Wingard | Corrected units for heading, pitch, roll, depth |
| 2015-06-10 | Russell Desiderio | Removed depth dependence from magnetic declination |
| 2015-06-25 | Russell Desiderio | Added fill value to NaN conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_earth_northward

#### Additional Notes
The `z` parameter (instrument pressure, daPa) is accepted by this function
but not used in the computation. The depth-to-dbar conversion lines in the
function body are commented out. This is a known issue in the current
implementation.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2014-02-03 | Christopher Wingard | Updated to use WMM 2010 magnetic declination |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2014-06-25 | Christopher Wingard | Corrected units for heading, pitch, roll, depth |
| 2015-06-10 | Russell Desiderio | Removed depth dependence from magnetic declination |
| 2015-06-25 | Russell Desiderio | Added fill value to NaN conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_earth_vertical

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2015-06-25 | Russell Desiderio | Added fill value to NaN conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_earth_error

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2015-06-25 | Russell Desiderio | Added fill value to NaN conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

### VELTURB -- VADCP-B (Nortek Signature 55)

::: ion_functions.data.adcp_functions.vadcp_b_beam_eastward

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code (as adcp_beam_eastward) |
| 2014-02-03 | Christopher Wingard | Updated to use WMM 2010 magnetic declination |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2014-06-25 | Christopher Wingard | Corrected units for heading, pitch, roll, depth |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions; removed depth dependence from magnetic declination |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2025-02-19 | Samuel Dahlberg | Duplicated and modified for VADCP-B beam geometry, units, and transformation matrix |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_b_beam_northward

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code (as adcp_beam_northward) |
| 2014-02-03 | Christopher Wingard | Updated to use WMM 2010 magnetic declination |
| 2014-03-28 | Russell Desiderio | Documentation correction |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2014-06-25 | Christopher Wingard | Corrected units for heading, pitch, roll, depth |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions; removed depth dependence from magnetic declination |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2025-02-19 | Samuel Dahlberg | Duplicated and modified for VADCP-B beam geometry, units, and transformation matrix |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_b_beam_vertical_est

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code (as vadcp_beam_vertical_est) |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions |
| 2015-06-22 | Russell Desiderio | Renamed data product |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2025-02-19 | Samuel Dahlberg | Duplicated and modified for VADCP-B beam geometry, units, and transformation matrix |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_b_beam_vertical_true

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code (as vadcp_beam_vertical_true) |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions |
| 2015-06-22 | Russell Desiderio | Renamed data product |
| 2015-06-25 | Russell Desiderio | Added beam 5 fill value to NaN conversion |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2025-02-19 | Samuel Dahlberg | Duplicated and modified for VADCP-B beam geometry, units, and transformation matrix |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.adcp_functions.adcp_beam2ins

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2015-06-24 | Russell Desiderio | Added fill value to NaN conversion |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_ins2earth

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2015-06-24 | Russell Desiderio | Corrected vertical fill value propagation; added fill value to NaN conversion |
| 2023-08-15 | Samuel Dahlberg | Updated local variable naming convention |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.magnetic_correction

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-04 | Russell Desiderio | Initial code |
| 2015-04-10 | Russell Desiderio | Corrected typo in array initialization |
| 2023-08-15 | Samuel Dahlberg | Updated local variable naming convention |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_b_beam2ins

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code (as adcp_beam2ins) |
| 2015-06-24 | Russell Desiderio | Added fill value to NaN conversion |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2025-02-19 | Samuel Dahlberg | Duplicated and modified for VADCP-B beam geometry, units, and transformation matrix |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_b_ins2earth

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-10 | Christopher Wingard | Initial code (as adcp_ins2earth) |
| 2014-04-04 | Russell Desiderio | Replaced for loops with np.einsum |
| 2015-06-24 | Russell Desiderio | Corrected vertical fill value propagation; added fill value to NaN conversion |
| 2025-02-19 | Samuel Dahlberg | Duplicated and modified for VADCP-B heading offset, extended rotation matrix, and averaged vertical velocity |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_bin_depths

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-29 | Craig Risien | Initial code |
| 2015-06-26 | Russell Desiderio | Fixed pressure variable handling; time-vectorized code |
| 2015-06-30 | Russell Desiderio | Added fill value to NaN conversion |
| 2019-03-11 | Christopher Wingard | Removed OOI CI constraints; made depth input type-agnostic |
| 2023-08-15 | Samuel Dahlberg | Ported from Pyseas for CGSN compatibility |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_bin_depths_bar

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-29 | Craig Risien | Initial code |
| 2015-06-26 | Russell Desiderio | Fixed pressure variable handling; time-vectorized code |
| 2015-06-30 | Russell Desiderio | Added fill value to NaN conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_bin_depths_dapa

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-29 | Craig Risien | Initial code |
| 2015-06-26 | Russell Desiderio | Fixed pressure variable handling; time-vectorized code |
| 2015-06-30 | Russell Desiderio | Added fill value to NaN conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.adcp_bin_depths_meters

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-30 | Craig Risien | Initial code |
| 2015-06-26 | Russell Desiderio | Time-vectorized code |
| 2015-06-30 | Russell Desiderio | Added fill value to NaN conversion |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.depth_from_dbar

#### History
| Date | Author | Change |
|---|---|---|
| 2025-02-19 | Samuel Dahlberg | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.depth_from_pressure_dbar

#### History
| Date | Author | Change |
|---|---|---|
| 2019-06-29 | Mark Steiner | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_b_bin_depths

#### History
| Date | Author | Change |
|---|---|---|
| 2025-02-19 | Samuel Dahlberg | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.z_from_p

#### Additional Notes
This function is a local implementation of the TEOS-10 `z_from_p` algorithm,
included in this module before the `gsw` library was available in the OOI
processing environment. It is a candidate for replacement with `gsw.z_from_p`
in a future refactor. See `enthalpy_SSO_0_p` below.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-07-01 | Russell Desiderio | Updated to TEOS-10 ver 3.05 |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.enthalpy_SSO_0_p

#### Additional Notes
This function is a local implementation of a TEOS-10 subroutine called by
`z_from_p`. Both are candidates for replacement with `gsw.z_from_p` in a
future refactor.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-07-01 | Russell Desiderio | Updated to TEOS-10 ver 3.05 |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Deprecated Functions

!!! warning "Deprecated"
    The functions in this section are deprecated. They implement the original
    VADCP processing pipeline (Teledyne RDI Workhorse Sentinel) using
    `adcp_beam2ins` and `adcp_ins2earth`. They have been superseded by the
    VADCP-B wrapper functions (`vadcp_b_beam_eastward`,
    `vadcp_b_beam_northward`) for new deployments. `vadcp_beam_vertical_est`
    and `vadcp_beam_vertical_true` are **not** deprecated; see
    [Core Functions](#core-functions).

::: ion_functions.data.adcp_functions.vadcp_beam_eastward

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions; removed depth dependence from magnetic declination |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_beam_northward

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2015-06-10 | Russell Desiderio | Moved beam and compass conditioning to helper functions; removed depth dependence from magnetic declination |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.adcp_functions.vadcp_beam_error

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2015-06-10 | Russell Desiderio | Moved beam conditioning to helper function |
| 2019-08-13 | Christopher Wingard | Added 3-beam solution support |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

OOI (2013). Data Product Specification for Turbulent Velocity Profile and Echo
Intensity. Document Control Number 1341-00760.
https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00760_Data_Product_SPEC_VELTURB_ECHOINT_OOI.pdf

OOI (2020). Data Product Specification for Velocity Profile and Echo
Intensity. Document Control Number 1341-00750.
https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00750_Data_Product_SPEC_VELPROF_ECHOINT_OOI.pdf

Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate polynomial
expressions for the density and specific volume of seawater using the TEOS-10
standard. Ocean Modelling.

Teledyne RD Instruments. *ADCP Coordinate Transformation, Formulas and
Calculations*. P/N 951-6079-00. San Diego: Teledyne RD Instruments, 2010.

Teledyne RD Instruments. *Workhorse Commands and Output Data Format*.
P/N 957-6156-00. San Diego: Teledyne RD Instruments, 2010.
