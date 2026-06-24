# PRS Functions

## Background

The Seafloor Pressure (PRS) family covers the Bottom Pressure and Tilt
(BOTPT) instrument class deployed by OOI on the Regional Cabled Array
(RCA) at Axial Seamount. The BOTPT instrument integrates three sensor
subsystems on a single seafloor frame: a Paroscientific Digiquartz
nano-resolution pressure transducer (Model 42.4K-265) paired with a
Paroscientific Intelligent Interface Board, an Applied Geomechanics LILY
Self-Leveling Borehole Tiltmeter, and two engineering-grade tiltmeters
(IRIS low-resolution and ADXL327 coarse-resolution). All PRS processing
functions are implemented in `prs_functions.py`.

| Class | Hardware | Platform | Designator Meaning |
|---|---|---|---|
| BOTPT | Paroscientific 42.4K-265 + LILY tiltmeter | Seafloor (RCA) | Bottom pressure and tilt |

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00060 | [OOI (2013). Data Product Specification for Seafloor High-Resolution Tilt.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00060_Data_Product_Spec_BOTTILT_OOI.pdf) |
| 1341-00070 | [OOI (2013). Data Product Specification for Nano-resolution Bottom Pressure.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00070_Data_Product_SPEC_BOTPRES_OOI.pdf) |
| 1341-00080 | OOI (2015). Data Product Specification for Seafloor Uplift and Subsidence (BOTSFLU) from the BOTPT Instrument. Document Control Number 1341-00080. (Not released.) |

---

### BOTPRES_L1 -- Nano-resolution Bottom Pressure

BOTPRES_L1 is the nano-resolution bottom pressure (psi) at the seafloor,
produced by the Paroscientific Digiquartz Pressure Transducer (Model
42.4K-265) integrated in the BOTPT instrument. As specified in DPS
1341-00070, the L1 product is computed entirely onboard by the paired
pressure transducer and Intelligent Interface Board using
temperature-compensated frequency-to-pressure conversion; no
off-instrument computation is required. The instrument transmits the
result directly as the L1 BOTPRES data product in psi (as absolute pressure,
seawater and atmospheric pressure) at 40 Hz. ion-functions is not invoked 
for this product.

---

### BOTTILT -- Seafloor High-Resolution Tilt

BOTTILT_L1 is the seafloor high-resolution tilt data product (L1) derived
from the Applied Geomechanics LILY Self-Leveling Borehole Tiltmeter
integrated in the BOTPT instrument. The LILY sensor outputs X-tilt and
Y-tilt components (BOTTILT-XTLT_L0, BOTTILT-YTLT_L0, microradians) and an
uncorrected sensor compass direction (BOTTILT-SCMP_L0, degrees) at 1 Hz.
Three L1 sub-products are computed from these L0 inputs as specified in
DPS 1341-00060.

**BOTTILT-CCMP_L1** (Corrected Compass Direction, integer degrees CW from
north) is derived by applying a sensor-specific lookup table to the rounded
L0 compass reading. The lookup table, stored in `prs_functions_ccmp.py`,
corrects for calibration offsets and for the 17-degree east magnetic
declination at Axial Seamount. The L0 sensor compass direction is the
heading of the negative Y-tilt axis measured CCW from north; the lookup
table converts it to the azimuth of the positive Y-tilt axis measured CW
from north.

**BOTTILT-TMAG_L1** (Tilt Magnitude, microradians) is the resultant tilt
magnitude computed from the X- and Y-tilt L0 components:

$$TMAG = \sqrt{XTLT^2 + YTLT^2}$$

**BOTTILT-TDIR_L1** (Tilt Direction, integer degrees CW from north) is
the azimuth of the resultant downward tilt vector. The algorithm, specified
in DPS 1341-00060, uses `arctan2` to compute the angle between the tilt
vector and the corrected compass reference:

$$TDIR = (450 - \arctan2(YTLT,\ XTLT)_{deg} + CCMP) \bmod 360$$

where the 450-degree addend ensures a positive result before the modulo
operation. The output is rounded to the nearest integer degree.

---

### BOTSFLU -- Seafloor Uplift and Subsidence

The BOTSFLU data products characterize seafloor vertical displacement and
its rate of change at Axial Seamount, derived from the BOTPRES_L1 pressure
time series after de-tiding. All BOTSFLU products are computed in two
binning stages: 20 Hz raw pressure data are first binned into 15-second
mean values, and those 15-second records are then processed into daily mean
depths and multi-week rate estimates.

**Tide prediction.** Predicted tides at Axial Seamount caldera center
(lat = 45.95547, lon = -130.00957) are stored in a precomputed lookup
table (`prs_functions_tides_2014_thru_2019.mat`) at 15-second resolution
from 2014-01-01 through 2020-01-01. Tide values were computed using the
Tide Model Driver software with the TPXO7.2 global tidal model. All three
BOTPT sites at Axial Seamount are sufficiently close together that a single
caldera-center location is used for all. Tide values are stored as signed
4-byte integers in units of 0.001 mm and scaled to meters on read.

**15-second binning.** Raw pressure data are binned into 15-second
time-centered means anchored at the quarter-minute (0, 15, 30, 45 seconds
past each minute). Bins are identified by floored elapsed time in units of
15 seconds using `numpy.bincount`. NaN values and non-physical readings
(pressure <= 0) are discarded before binning. Empty bins are not
represented; the boolean mask `mask_nonzero` records the positions of
non-empty bins within the full time span for use by downstream products.

**Depth conversion and de-tiding.** Binned mean pressure (MEANPRES_L2,
psi) is converted to depth and de-tided to produce MEANDEPTH_L2 (m). 
Atmospheric pressure is not subtracted from the L1 pressure data despite 
its units of psi. The conversion factor is -0.67 m/psi, and depth follows
a negative convention, so the predicted tide (m) is added to de-tide the 
record:

$$MEANDEPTH = (MEANPRES \times -0.67) + PREDTIDE$$

**Daily binning.** The 15-second MEANDEPTH record is binned into 24-hour
noon-to-noon bins anchored at midnight to produce DAYDEPTH_L2 (m). A
fractional coverage threshold (default 0.90) controls quality: bins with
fewer than 90% of the maximum possible 5760 values per day are assigned
NaN. All days in the record span are represented, including data-gap days.

**Rate products.** Short-term rates are computed directly from the
15-second MEANDEPTH record; long-term rates use backwards-looking linear
regression on the daily DAYDEPTH record.

- **5MINRATE_L2** (cm/min): instantaneous rate using a 20-interval
  (5-minute) backward difference of MEANDEPTH.
- **10MINRATE_L2** (cm/hr): mean uplift rate using a 40-interval
  (10-minute) sliding window mean of MEANDEPTH, differenced over 40
  intervals.
- **4WKRATE_L2** (cm/yr): mean rate from a 29-day backwards-looking linear
  regression on DAYDEPTH; slope converted from m/day to cm/yr by
  multiplying by 36500.
- **8WKRATE_L2** (cm/yr): mean rate from a 57-day backwards-looking linear
  regression on DAYDEPTH; converted identically.

The regression algorithm (used for 4WKRATE and 8WKRATE) applies the
normal linear regression equations from Press et al. (1986) in a
backwards-looking sliding window. Windows with fewer non-NaN values than
specified by a fractional coverage threshold (default 0.75) are assigned
NaN. The data vector is front-padded with NaN so that windows near the
start of the record that satisfy the coverage criterion produce non-NaN
values.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.prs_functions.prs_bottilt_ccmp

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-10 | Christopher Wingard | Initial code |
| 2014-03-20 | Russell Desiderio | Alternate code: faster, but less direct |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_bottilt_tmag

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-10 | Christopher Wingard | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_bottilt_tdir

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-10 | Christopher Wingard | Initial code |
| 2014-03-20 | Russell Desiderio | Replaced initial code with arctan2 implementation |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_time15s

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-13 | Russell Desiderio | Initial code |
| 2017-05-15 | Russell Desiderio | Included botpres as input argument to delete timestamps of bad raw values |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_meanpres

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-13 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_predtide

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-13 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_meandepth

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_5minrate

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_10minrate

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_time24h

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2017-05-05 | Russell Desiderio | Changed time24h time base to span entire dataset including data gaps |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_daydepth

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2017-05-05 | Russell Desiderio | Changed time24h time base to span entire dataset including data gaps |
| 2017-05-11 | Russell Desiderio | Incorporated daydepth coverage threshold |
| 2020-06-01 | Mark Steiner | Required predtide as argument; refactored daydepth calculation to expand API |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_4wkrate

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2017-05-05 | Russell Desiderio | Changed time24h time base to span entire dataset including data gaps |
| 2017-05-12 | Russell Desiderio | Incorporated daydepth and rate coverage thresholds |
| 2020-06-01 | Mark Steiner | Required predtide as argument; refactored rate calculation to expand API |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_8wkrate

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2017-05-05 | Russell Desiderio | Changed time24h time base to span entire dataset including data gaps |
| 2017-05-12 | Russell Desiderio | Incorporated daydepth and rate coverage thresholds |
| 2020-06-01 | Mark Steiner | Required predtide as argument; refactored rate calculation to expand API |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.prs_functions.prs_botsflu_daydepth_from_15s_meandepth

#### History
| Date | Author | Change |
|---|---|---|
| 2020-06-01 | Mark Steiner | Extracted from prs_botsflu_daydepth to expand API |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_4wkrate_from_daydepth

#### History
| Date | Author | Change |
|---|---|---|
| 2020-06-01 | Mark Steiner | Extracted from prs_botsflu_4wkrate to expand API |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_botsflu_8wkrate_from_daydepth

#### History
| Date | Author | Change |
|---|---|---|
| 2020-06-01 | Mark Steiner | Extracted from prs_botsflu_8wkrate to expand API |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.anchor_bin_raw_data_to_15s

#### History
| Date | Author | Change |
|---|---|---|
| 2017-05-05 | Russell Desiderio | Initial code from modifying anchor_bin; added nan and bad-value trapping |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.anchor_bin_detided_data_to_24h

#### History
| Date | Author | Change |
|---|---|---|
| 2017-05-05 | Russell Desiderio | Initial code from modifying anchor_bin |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.calc_meandepth_plus

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.calculate_sliding_means

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-13 | Russell Desiderio | Initial code |
| 2017-05-06 | Russell Desiderio | Updated to work with odd window sizes |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.calculate_sliding_slopes

#### History
| Date | Author | Change |
|---|---|---|
| 2017-05-03 | Russell Desiderio | Initial code; replaced Moore-Penrose method to support nan-masking |
| 2017-05-08 | Russell Desiderio | Added fractional coverage criterion |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Deprecated Functions

!!! note "Deprecated"
    The functions in this section are retained for reference only and
    should not be used in new code. They have been superseded by the
    functions listed above.

::: ion_functions.data.prs_functions.prs_tsunami_detection

#### Additional Notes
The DPS document for BOTSFLU (1341-00080) was never publicly released.
This function was coded from pseudocode in that DPS. Its robustness has
not been verified with actual data.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_eruption_imminent

#### Additional Notes
The DPS document for BOTSFLU (1341-00080) was never publicly released.
This function was coded from pseudocode in that DPS. Its robustness has
not been verified with actual data.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.prs_eruption_occurred

#### Additional Notes
The DPS document for BOTSFLU (1341-00080) was never publicly released.
This function was coded from pseudocode in that DPS. Its robustness has
not been verified with actual data.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.anchor_bin

#### Additional Notes
Deprecated May 2017. Superseded by `anchor_bin_raw_data_to_15s` and
`anchor_bin_detided_data_to_24h`, which handle the raw-data bad-value
check and extended 24-hour timestamp records respectively.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-13 | Russell Desiderio | Initial code |
| 2015-01-14 | Russell Desiderio | Changed output arguments; added conditionals for efficiency |
| 2017-05-05 | Russell Desiderio | Deprecated |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.calc_daydepth_plus

#### Additional Notes
Deprecated May 2017. Superseded by `calc_meandepth_plus` and
`prs_botsflu_daydepth`.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-14 | Russell Desiderio | Initial code |
| 2017-05-05 | Russell Desiderio | Deprecated |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.calculate_all_sliding_slopes_then_Nan

#### Additional Notes
Deprecated May 2017. Superseded by `calculate_sliding_slopes`, which
pre-filters windows by fractional coverage before computing slopes rather
than computing all slopes and then applying NaN masking.

#### History
| Date | Author | Change |
|---|---|---|
| 2017-05-03 | Russell Desiderio | Initial code |
| 2017-05-08 | Russell Desiderio | Added coverage criterion |
| 2017-05-08 | Russell Desiderio | Deprecated |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.prs_functions.calculate_sliding_slopes__MoorePenrose

#### Additional Notes
Deprecated May 2017. Superseded by `calculate_sliding_slopes`. The
Moore-Penrose pseudoinverse method cannot accommodate NaN values in the
data window.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-01-13 | Russell Desiderio | Initial code |
| 2017-05-03 | Russell Desiderio | Deprecated; nan-masking cannot be implemented with this algorithm |
| 2025-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

[OOI (2013). Data Product Specification for Seafloor High-Resolution
Tilt. Document Control Number
1341-00060.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00060_Data_Product_Spec_BOTTILT_OOI.pdf)

[OOI (2013). Data Product Specification for Nano-resolution Bottom
Pressure. Document Control Number
1341-00070.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00070_Data_Product_SPEC_BOTPRES_OOI.pdf)

OOI (2015). Data Product Specification for Seafloor Uplift and
Subsidence (BOTSFLU) from the BOTPT Instrument. Document Control Number
1341-00080. (Not released.)

Polster, A., Fabian, M., and Villinger, H. (2009). Effective resolution
and drift of Paroscientific pressure sensors derived from long-term
seafloor measurements. *Geochemistry, Geophysics, Geosystems*, 10,
Q08009. doi:10.1029/2009GC002532.

Press, W.H., Flannery, B.P., Teukolsky, S.A., and Vetterling, W.T.
(1986). *Numerical Recipes*. Cambridge University Press. p. 507.
