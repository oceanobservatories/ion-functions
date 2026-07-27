# FDC Functions

!!! note "Inactive Module"
    These functions were used to process data from the FDCHP instrument. They
    are no longer in use because the OOI data system assembles instrument
    bursts in a way that is incompatible with this processing pipeline. They
    are retained for reference but are not used in OOI data production.

## Background

The Flux Direct Covariance High Power (FDCHP) system directly computes air-sea fluxes of
momentum and buoyancy using the motion-corrected direct covariance (MCDC) method. The system
combines a Gill Windmaster Pro (Model 1561-PK-020) 3-axis sonic anemometer/thermometer with a
Lord MicroStrain 3DM-GX3-25 Inertial Measurement Unit (IMU, integrating 3-axis linear
accelerometers, 3-axis angular rate sensors, and a 3-axis magnetometer), forming the FDCHP
instrument class deployed on select OOI surface buoys. MicroStrain was acquired by Hottinger
Bruel and Kjaer (HBK) in September 2023 and now operates as MicroStrain by HBK.

FDCHP collects 10 Hz data for 20 minutes out of every hour, yielding roughly 12000 records per
dataset. `fdc_functions.py` parses the raw 1D input streams into discrete 12000-record dataset
packets, computes the motion-corrected L1 wind velocity (WINDTUR) and sonic temperature
(TMPATUR) products, and derives the L2 momentum (FLUXMOM) and buoyancy (FLUXHOT) flux products
from those L1 products. Two auxiliary time-base products, not specified in the DPS, provide
timestamps for the L1 and L2 arrays.

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00280 | [OOI (2015). Data Product Specification for FDCHP Data Products. Document Control Number 1341-00280.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00280_Data_Product_SPEC_DCHPFLX.pdf) |

### WINDTUR_L1 -- Motion-Corrected Wind Velocity

WINDTUR_L1 is the motion-corrected wind velocity (m/s) relative to Earth, reported as three
directional components: WINDTUR-VLN_L1 (positive North), WINDTUR-VLW_L1 (positive West), and
WINDTUR-VLU_L1 (positive up). The sign convention follows a right-handed (x, y, z) buoy frame
with x positive toward the buoy vane, y positive to port, z positive upward, and yaw positive
for counter-clockwise rotation (the opposite of the typical left-handed compass convention).

#### Instrument Frame Data and Calibration

Raw sonic anemometer velocities and speed of sound are converted from counts to physical units
with a gain of 0.01 (cm/s to m/s). Raw angular rates and heading/pitch/roll are stored as
decimal counts equal to the physical value directly (gain of 1, no unit conversion).

The DPS Table 1 (Section 4.2) states a gain of 1 for the accelerometers, implying raw counts
are already in m/s$^2$. This conflicts with the DPS's own Appendix A MATLAB reference code,
which multiplies the raw counts by G = 9.80665 (`platform = dcfsdata(13:15,1:L)*G;`), and with
the module's own docstrings, which label the raw units as counts of 9.80665 m/s$^2$ (i.e. units
of g). Since the code and the DPS's own reference implementation agree with each other, and
only the summary table disagrees, `ion-functions` follows the code and Appendix A: raw
accelerometer counts are multiplied by G = 9.80665 to obtain m/s$^2$. Table 1 is the outlier
here.

Before use, the IMU angular rate (y, z), acceleration (y, z), and heading are negated to
convert from the IMU's North-East-Down convention to the sonic anemometer's North-West-Up
convention.

#### Euler Angle Estimation

Platform orientation is described by Euler angles roll ($\phi$), pitch ($\theta$), and yaw
($\psi$), combined into the coordinate transformation matrix:

$$T(\phi, \theta, \psi) = A(\psi) A(\theta) A(\phi)$$

$$A(\psi) = \begin{pmatrix}
\cos\psi & -\sin\psi & 0 \\
\sin\psi & \cos\psi & 0 \\
0 & 0 & 1
\end{pmatrix}, \quad
A(\theta) = \begin{pmatrix}
\cos\theta & 0 & \sin\theta \\
0 & 1 & 0 \\
-\sin\theta & 0 & \cos\theta
\end{pmatrix}, \quad
A(\phi) = \begin{pmatrix}
1 & 0 & 0 \\
0 & \cos\phi & -\sin\phi \\
0 & \sin\phi & \cos\phi
\end{pmatrix}$$

Because the FDCHP is a strapped-down system (the motion sensors are fixed to the buoy frame,
not gyro-stabilized), the Euler angles are not measured directly and must be estimated from the
angular rate sensors. The time derivatives of the Euler angles are related to the
strapped-down angular rate vector by:

$$\begin{pmatrix} \dot\phi \\ \dot\theta \\ \dot\psi \end{pmatrix} =
\begin{pmatrix}
\dot\phi_{obs} + [\dot\psi_{obs}\cos\phi + \dot\theta_{obs}\sin\phi]\tan\theta \\
\dot\theta_{obs}\cos\phi - \dot\psi_{obs}\sin\phi \\
[\dot\psi_{obs}\cos\phi + \dot\theta_{obs}\sin\phi] / \cos\theta
\end{pmatrix}$$

implemented in `fdc_update`. Integrating the angular rates directly is subject to low-frequency
sensor drift, so a complementary filter combines the integrated (fast) angle estimate with a
low-frequency (slow) reference angle derived from accelerometer tilt, using a 4th-order
Butterworth high-pass filter (`fdc_filtcoef`) applied forward and backward (zero-phase) to
avoid distorting the phase:

$$\phi \approx \left[\frac{\ddot y_{obs}}{g} -
HP\left(\frac{\ddot y_{obs}}{g}\right)\right] + HP(\phi_{obs})$$

where $HP$ denotes the zero-phase high-pass filter operator and $g$ is the gravitational
acceleration from `fdc_grv`. The slow yaw angle is derived from the unwrapped gyro (compass)
signal when the compass is judged reliable (`goodcompass`, from `fdc_process_compass_data`);
otherwise it is held at the median gyro value. Starting from the slow angles as a first guess,
`fdc_anglesclimodeyaw` refines the Euler angle estimate over 5 iterations, each time
integrating and filtering the updated angular rates and adding the result to the slow angles.
The high-pass filter cutoff frequency is hardcoded to 1/12 Hz per e-mail guidance from DPS
author Jim Edson, consistent with the DPS's discussion of theoretical limitations (Section
3.3), which notes the initial code implementation used the same cutoff.

#### Platform Velocity

The Euler angle estimate is used to rotate the measured linear accelerations into the earth
frame, remove the gravity vector, and integrate to platform velocity, high-pass filtering to
remove drift (`fdc_accelsclimode`):

$$\vec{V}_{hp} = HP\left[\int \left(T\ddot{\vec x}_{obs} + \vec g\right) dt\right]$$

where $\vec g = -\hat k g$.

#### True Wind Velocity

The measured wind velocity is combined with the angular-rate cross-product correction (which
accounts for the offset between the IMU and the sonic anemometer sampling volume) and rotated
into the earth frame, then added to the platform velocity to give the wind velocity relative to
the buoy (`fdc_sonic`, `fdc_trans`):

$$\vec{V}^{buoy}_{true} = T\left(\vec{V}_{obs} + \vec\Omega_{obs} \times \vec R\right) +
\vec{V}_{hp}$$

$\vec R$ is the fixed distance vector between the IMU and the sonic anemometer sampling volume,
hardcoded to 0.753 m vertical separation, matching the DPS Appendix A code. The first and last
30 seconds of $\vec{V}^{buoy}_{true}$ are discarded to remove filter edge effects, yielding the
19-minute WINDTUR-VLN_L1, WINDTUR-VLW_L1, and WINDTUR-VLU_L1 time series.

#### Streamwise Rotation

For the L2 flux calculations, the WINDTUR_L1 velocity components are rotated into the
longitudinal (streamwise) wind direction (`fdc_alignwind`), forcing the mean cross-wind and
vertical components to zero. This rotation has been shown to reduce the effect of flow
distortion on the fluxes (Oost et al., 1994).

### TMPATUR_L1 -- Sonic Temperature

TMPATUR_L1 is the sonic temperature ($^\circ$C) derived from the speed of sound measured by
the sonic anemometer:

$$T_{s}\ [^\circ C] = \frac{C_s^2}{403} - 273.15$$

where $C_s$ is the calibrated speed of sound (m/s). This conversion is applied both by the
standalone `fdc_tmpatur` wrapper and internally within `fdc_flux_and_wind` in the course of
computing FLUXHOT_L2; the two implementations are not linked by a function call and must be
kept in sync manually if the conversion changes.

### FLUXMOM_L2 -- Momentum Flux

FLUXMOM_L2 is the kinematic form of the momentum flux (m$^2$/s$^2$) (Stull, 1988), computed
using the motion-corrected direct covariance method. After streamwise rotation, the along-wind
and cross-wind velocity components are linearly detrended (removing a least-squares fit) to
give the fluctuating components $u'$, $v'$, $w'$, from which the along-wind and cross-wind
momentum flux components are computed over each 19-minute dataset:

$$\text{FLUXMOM-U\_L2} = \overline{u'w'}, \qquad \text{FLUXMOM-V\_L2} = \overline{v'w'}$$

The along-wind component (FLUXMOM-U_L2) generally carries most of the momentum flux; the
cross-wind component (FLUXMOM-V_L2) is usually smaller but can approach or exceed the along-wind
component near the ocean surface in the presence of waves, or in light winds where the wind and
stress vectors are poorly defined.

### FLUXHOT_L2 -- Buoyancy Flux

FLUXHOT_L2 is the kinematic form of the buoyancy flux (m/s K), computed as the covariance of
the detrended vertical wind fluctuation $w'$ with the detrended sonic temperature fluctuation
$T_s'$:

$$\text{FLUXHOT\_L2} = \overline{w'T_s'}$$

The buoyancy flux is properly defined in terms of fluctuations in virtual temperature, which the
sonic temperature closely approximates; the small correction needed to convert from sonic to
virtual temperature (using an independent estimate of the latent heat flux) is left to
post-processing and is not computed by `ion-functions`.

For a discussion of output accuracy, the DPS points to Bigorre et al. (2013) and Bradley and
Fairall (2006) (Appendix B); neither is cited with enough bibliographic detail elsewhere in the
DPS to reproduce a full reference here.

Full algorithm derivations, calibration procedures, and source references are listed in the
[References](#references) section.

---

## Core Functions

::: ion_functions.data.fdc_functions.fdc_flux_and_wind

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-20 | Russell Desiderio | Initial code |
| 2014-11-06 | Russell Desiderio | Incorporated fdc_quantize_data routine |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

::: ion_functions.data.fdc_functions.fdc_tmpatur

#### Additional Notes
Independently reproduces the sonic temperature conversion also computed internally by
`fdc_flux_and_wind`; the two are not linked by a function call, so a change to one requires a
matching change to the other.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2015-01-29 | Russell Desiderio | Changed units from kelvin to degrees Celsius |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_windtur_north

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2015-01-29 | Russell Desiderio | Removed temperature from calling arguments |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_windtur_up

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2015-01-29 | Russell Desiderio | Removed temperature from calling arguments |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_windtur_west

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2015-01-29 | Russell Desiderio | Removed temperature from calling arguments |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_fluxhot

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_fluxmom_alongwind

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2015-01-29 | Russell Desiderio | Removed temperature from calling arguments |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_fluxmom_crosswind

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2015-01-29 | Russell Desiderio | Removed temperature from calling arguments |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.fdc_functions.fdc_accelsclimode

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-19 | Russell Desiderio | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_alignwind

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-19 | Russell Desiderio | Initial code; converted arithmetic to matrix multiplication |
| 2014-05-29 | Russell Desiderio | Programmed original arithmetic, faster than matrix operations |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_anglesclimodeyaw

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-19 | Russell Desiderio | Initial code |
| 2014-09-25 | Russell Desiderio | Incorporated latest changes to pitch and roll calculation |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_despikesimple

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-19 | Craig Risien | Initial code |
| 2014-05-30 | Russell Desiderio | Vectorized some of the code |
| 2014-11-19 | Russell Desiderio | Made code more robust to avoid runtime errors on out-of-range interpolations |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_detrend

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-19 | Russell Desiderio | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_filtcoef

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-25 | Russell Desiderio | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_grv

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-15 | Russell Desiderio | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_process_compass_data

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-25 | Russell Desiderio | Separated out from main code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_quantize_data

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-06 | Russell Desiderio | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_sonic

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-16 | Craig Risien | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_time_L1

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_time_L2

#### History
| Date | Author | Change |
|---|---|---|
| 2014-11-17 | Russell Desiderio | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_trans

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-16 | Craig Risien | Initial code |
| 2014-05-30 | Russell Desiderio | Removed conditional iflag (always True) |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.fdc_functions.fdc_update

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-16 | Craig Risien | Initial code |
| 2026-07-24 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

Edson, J.B., Hinton, A.A., Prada, K.E., Hare, J.E., and Fairall, C.W. (1998). Direct covariance
flux estimates from mobile platforms at sea. J. Atmos. Oceanic Tech., 15, 547-562.

Fujitani, T. (1981). Direct measurement of turbulent fluxes over the sea during AMTEX. Pap.
Meteor. Geophys., 32, 119-134.

Fujitani, T. (1985). Method of turbulent flux measurement on a ship by using a stable platform
system. Pap. Meteor. Geophys., 36, 157-170.

Goldstein, H. (1965). Classical Mechanics. Addison-Wesley, 398 pp.

Miller, S., Friehe, C., Hristov, T., and Edson, J. (2008). Platform motion effects on
measurements of turbulence and air-sea exchange over the open ocean. J. Atmos. Oceanic Tech.,
25, 1683-1694.

[OOI (2015). Data Product Specification for FDCHP Data Products. Document Control Number
1341-00280.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00280_Data_Product_SPEC_DCHPFLX.pdf)

Oost, W.A., Fairall, C.W., Edson, J.B., Smith, S.D., Anderson, R.J., Wills, J.A.B., Katsaros,
K.B., and DeCosmo, J. (1994). Flow distortion calculations and their application in HEXMAX. J.
Atmos. Oceanic Technol., 11, 366-386.

Stull, R.B. (1988). An Introduction to Boundary Layer Meteorology. Kluwer Academic Publishers,
666 pp.

Ware, J. (2014). Interface Document for the Flux Direct Covariance High Power System. Ocean
Observatories Initiative internal document.
