# WAV Functions

## Background

The WAV instrument family consists of the TRIAXYS OEM Directional Wave Sensor
(WAVSS), manufactured by AXYS Technologies, Inc. The WAVSS measures
three-dimensional acceleration, rotation rate, and orientation in the Earth's
magnetic field. From these measurements the instrument reconstructs the motion
of the sea surface in the wave field. All wave statistics computations are
performed internally by the instrument's proprietary firmware; ion-functions
does not re-implement those algorithms.

The WAVSS acquires accelerations over a sampling interval (typically 20--30
minutes) and produces a single set of wave property values representing that
interval, along with a time series of platform displacement over the course of
the acquisition.

### Primary Sources

| DCN | Document |
|---|---|
| [1341-00450](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00450_Data_Product_SPEC_WAVSTAT_OOI.pdf) | OOI (2012). Data Product Specification for Wave Statistics. Document Control Number 1341-00450. |

### WAVSTAT Data Products

The WAVSS reports a large set of wave statistics under the collective product
name WAVSTAT. Most of these are extracted directly from the instrument's NMEA
output sentences without further computation by ion-functions. The table below
lists all WAVSTAT sub-products, their sources, and the ion-functions function
responsible where one exists.

| Product ID | Alt. | Description | Units | Source |
|---|---|---|---|---|
| WAVSTAT-N0 | N0 | Number of zero crossings in displacement data (QC use) | -- | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-HMAX | Hmax | Maximum wave height | m | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-HAVG | Havg | Average wave height | m | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-TAVG | Tavg | Average wave period | s | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-HSIG | Hsig | Significant wave height (average of highest 1/3) | m | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-HMO | Hmo | Significant wave height from spectral moments | m | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-TSIG | Tsig | Significant wave period (average period of Hsig waves) | s | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-H10 | H10 | Average height of highest 1/10 of waves | m | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-T10 | T10 | Average period of H10 waves | s | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-TP | TP | Peak wave period | s | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-TP5 | TP5 | Peak wave period via Read method (alternative to TP) | s | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-D_L0 | D | Mean wave direction (magnetic) | deg | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-D_L2 | D | Mean wave direction (true north) | deg | `wav_triaxys_correct_mean_wave_direction` |
| WAVSTAT-DS | DS | Mean directional spread of wave field | deg | Instrument firmware; extracted from $TSPWA |
| WAVSTAT-PND | PND | Power spectral density, non-directional spectra | m$^2$ Hz$^{-1}$ | Instrument firmware; extracted from $TSPNA |
| WAVSTAT-FND_L1 | FND | Frequency values for non-directional spectral bins | Hz | `wav_triaxys_nondir_freq` |
| WAVSTAT-FDS_L1 | FDS | Frequency values for directional spectral bins | Hz | `wav_triaxys_dir_freq` |
| WAVSTAT-PDS | PDS | Power spectral density, directional spectra | m$^2$ Hz$^{-1}$ | Instrument firmware; extracted from $TSPMA |
| WAVSTAT-DDS_L0 | DDS | Wave directions from directional spectra (magnetic) | deg | Instrument firmware; extracted from $TSPMA |
| WAVSTAT-DDS_L2 | DDS | Wave directions from directional spectra (true north) | deg | `wav_triaxys_correct_directional_wave_direction` |
| WAVSTAT-SDS | SDS | Directional spread from directional spectra | deg | Instrument firmware; extracted from $TSPMA |
| WAVSTAT-MOTX_L0 | X | Eastward buoy displacement (magnetic frame) | m | Instrument firmware; extracted from $TSPHA |
| WAVSTAT-MOTX_L1 | X | Eastward buoy displacement (true east) | m | `wav_triaxys_magcor_buoymotion_x` |
| WAVSTAT-MOTY_L0 | Y | Northward buoy displacement (magnetic frame) | m | Instrument firmware; extracted from $TSPHA |
| WAVSTAT-MOTY_L1 | Y | Northward buoy displacement (true north) | m | `wav_triaxys_magcor_buoymotion_y` |
| WAVSTAT-MOTZ | Z | Vertical buoy displacement | m | Instrument firmware; extracted from $TSPHA |
| WAVSTAT-MOTT_L1 | t | Time of each buoy displacement measurement | s since 1900-01-01 | `wav_triaxys_buoymotion_time` |

### Frequency Vector Reconstruction

WAVSTAT-FND_L1 and WAVSTAT-FDS_L1 are the frequency axis vectors for the
non-directional and directional power spectral density arrays, respectively.
Rather than transmitting every frequency value, the instrument reports only
the number of bins, the initial frequency, and the uniform frequency spacing.
ion-functions reconstructs the full frequency vector from these three values.

For the non-directional spectrum, all data packets share the same number of
frequency bins, so the reconstruction is vectorized across packets. For the
directional spectrum, the number of active bins (nfreq_dir) can vary between
data packets as a function of ocean conditions -- it is always less than or
equal to the number of non-directional bins. The output array is therefore
sized to the (fixed) non-directional bin count, with fill values occupying
unused positions.

The reconstruction for both products follows:

$$f_i = f_0 + (i - 1) \times \Delta f, \quad i = 1, 2, \ldots, N$$

where $f_0$ is the initial frequency [Hz], $\Delta f$ is the frequency spacing
[Hz], and $N$ is the number of active bins for that data packet.

### Platform Motion Time Reconstruction

WAVSTAT-MOTT_L1 is the time axis corresponding to the WAVSTAT-MOTX,
WAVSTAT-MOTY, and WAVSTAT-MOTZ displacement time series. The instrument
reports the NTP timestamp of the data sentence, the elapsed time to the first
displacement sample, and the uniform sampling interval. ion-functions
reconstructs the absolute NTP time of each sample:

$$t_i = t_0 + t_{init} + (i - 1) \times \Delta t, \quad i = 1, 2, \ldots, N$$

where $t_0$ is the NTP timestamp from the $TSPHA sentence [s since 1900-01-01],
$t_{init}$ is the initial time offset to the first sample [s], $\Delta t$ is
the time spacing between samples [s], and $N$ is the number of displacement
measurements.

### Magnetic Declination Correction

The WAVSS contains an internal compass. Wave directions (WAVSTAT-D and
WAVSTAT-DDS) and buoy displacement components (WAVSTAT-MOTX and WAVSTAT-MOTY)
are reported in the magnetic reference frame. ion-functions corrects these to
true north using the WMM2010 magnetic declination model via
`generic_functions.magnetic_declination`.

The direction correction applies a rotation modulo 360 deg:

$$D_{true} = (D_{mag} + \theta + 360) \bmod 360$$

where $\theta$ is the magnetic declination [deg] at the instrument's location
and time. The same formula applies element-wise to the WAVSTAT-DDS array;
fill-valued positions are preserved through the operation.

The buoy displacement correction rotates the (X, Y) coordinate pair by
$\theta$ using a 2-D rotation matrix:

$$\begin{align}
X_{true} &= X_{mag} \times \cos\theta + Y_{mag} \times \sin\theta \\
Y_{true} &= -X_{mag} \times \sin\theta + Y_{mag} \times \cos\theta
\end{align}$$

The WAVSS is a surface sensor; the depth parameter for the declination
calculation defaults to 0 (sea level).

### ADCP Wave Statistics

Wave statistics are also measured by RDI/Teledyne acoustic Doppler current
profilers (ADCPs) mounted on the seafloor at the two shallowest OOI Endurance
Inshore mooring locations. Those instruments process wave data using RDI's
WavesMon software, which produces outputs equivalent to the WAVSS WAVSTAT
products. ion-functions is not involved in that processing pipeline.

Full algorithm derivations, processing references, and source documentation
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.wav_functions.wav_triaxys_nondir_freq

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-03 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_dir_freq

#### Additional Notes

The number of active directional frequency bins (nfreq_dir) can vary between
data packets as a function of measured ocean conditions at fixed instrument
settings. As a result, the size of the WAVSTAT-FDS_L1 output array, and
correspondingly of the WAVSTAT-PDS and WAVSTAT-SDS arrays and the
WAVSTAT-DDS_L2 product, is not fixed across packets. The output is sized to
the (fixed) non-directional bin count, with fill values in unused positions.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-03 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_buoymotion_time

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-07 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_correct_mean_wave_direction

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-08 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_correct_directional_wave_direction

#### Additional Notes

The variable-length behavior of nfreq_dir described in the Additional Notes
for `wav_triaxys_dir_freq` applies equally here: the WAVSTAT-DDS_L2 output
array is sized to the non-directional bin count, with fill values in unused
positions. The correction preserves fill-valued positions by converting them
to NaN before applying the rotation and restoring fills afterward.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-09 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.wav_functions.magnetic_correction_einsum

#### Additional Notes

The executable code in this function is identical to
`magnetic_correction_vctrzd` in `adcp_functions.py`. It was written
separately to handle the vectorized (i, j) case of one magnetic declination
per ensemble of (u, v) pairs without a for loop. The ADCP DPS citations in
the docstring reflect the origin of the rotation algorithm; the function
operates identically for (X, Y) displacement coordinates as for (u, v)
velocity components.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-04 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

::: ion_functions.data.wav_functions.wav_triaxys_magcor_buoymotion_x

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-10 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.wav_functions.wav_triaxys_magcor_buoymotion_y

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-10 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

[OOI (2012). Data Product Specification for Wave Statistics. Document Control
Number 1341-00450.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00450_Data_Product_SPEC_WAVSTAT_OOI.pdf)

[OOI (2020). Data Product Specification for Velocity Profile and Echo
Intensity. Document Control Number 1341-00750.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00750_Data_Product_SPEC_VELPROF_ECHOINT_OOI.pdf)

[OOI (2013). Data Product Specification for Turbulent Velocity Profile and
Echo Intensity. Document Control Number 1341-00760.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00760_Data_Product_SPEC_VELTURB_ECHOINT_OOI.pdf)

AXYS Technologies Inc. *TRIAXYS OEM Directional Wave Sensor, User's Manual.*
Sidney, BC: AXYS Technologies Inc. May 2005.

Teledyne RD Instruments. *WavesMon v3.08 User's Guide.* P/N 957-6232-00.
Teledyne RD Instruments.

Teledyne RD Instruments. *Waves Primer: Wave Measurements and the RDI ADCP
Waves Array Technique.* Teledyne RD Instruments.
