# QC Functions

## Background

`qc_functions.py` implements the automated quality control (QC) algorithms specified in the OOI Data Product
Specifications (DPS) for quality control, ported from the MATLAB example code in those documents. Each test examines a
data product and returns a flag for every data point. The tests are configured with parameters (valid ranges,
accuracies, window lengths, thresholds) taken from QC lookup tables, which the DPS documents keep separate from the
algorithms themselves.

### Primary Sources

| DCN | Document |
|---|---|
| 1341-10001 | [OOI (2012). Data Product Specification for Modulus (MODULUS).](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10001_Data_Product_SPEC_MODULUS_OOI.pdf) |
| 1341-10002 | [OOI (2012). Data Product Specification for 1-D Interpolation (INTERP1).](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10002_Data_Product_SPEC_INTERP1_OOI.pdf) |
| 1341-10003 | [OOI (2012). Data Product Specification for Evaluate Polynomial (POLYVAL).](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10003_Data_Product_SPEC_POLYVAL_OOI.pdf) |
| 1341-10004 | [OOI (2012). Data Product Specification for Global Range Test.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10004_Data_Product_SPEC_GLBLRNG_OOI.pdf) |
| 1341-10005 | [OOI (2012). Data Product Specification for Local Range Test.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10005_Data_Product_SPEC_LOCLRNG_OOI.pdf) |
| 1341-10006 | [OOI (2012). Data Product Specification for Spike Test.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10006_Data_Product_SPEC_SPKETST_OOI.pdf) |
| 1341-10007 | [OOI (2012). Data Product Specification for Trend Test.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10007_Data_Product_SPEC_TRNDTST_OOI.pdf) |
| 1341-10008 | [OOI (2012). Data Product Specification for Stuck Value Test.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10008_Data_Product_SPEC_STUCKVL_OOI.pdf) |
| 1341-10010 | [OOI (2012). Data Product Specification for Gradient Test.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10010_Data_Product_SPEC_GRADTST_OOI.pdf) |
| 1341-10011 | [OOI (2012). Data Product Specification for Solar Elevation.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10011_Data_Product_SPEC_SOLAREL_OOI.pdf) |
| 1341-10012 | [OOI (2012). Data Product Specification for Combined QC Flags.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10012_Data_Product_SPEC_CMBNFLG_OOI.pdf) |
| 1341-10030 | [OOI (2012). Data Product Specification for Conductivity Compressibility Compensation.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10030_Data_Product_SPEC_CONDCMP_OOI.pdf) |

### Flag Convention

Following the DPS documents, a flag of 1 marks data that passed a test (presumed good) and 0 marks data that failed it
(presumed bad). Two further cases return -99:

- The OOI wrapper functions return -99 for every data point when a required test parameter is missing or holds the
  QC lookup table fill value of -9999, meaning the test was not run.
- The gradient test returns -99 for data points it skipped because they are closer than the minimum separation in the
  coordinate (`mindx`) to the previous point.

-99 does not appear in OOI data sets. The OOI stream engine converts any result other than 0 or 1 to 0 (and logs an
error), so a test that returns -99 is recorded as executed and failed. For each tested parameter, the delivered data
carry two bitmap variables, `<parameter>_qc_executed` and `<parameter>_qc_results`. Each bit corresponds to one test: in `qc_executed` a set bit
means the test was run, and in `qc_results` a set bit means the data passed (1) and a clear bit that they failed (0).

| Bit | Test |
|---|---|
| 0 | `dataqc_globalrangetest` |
| 1 | `dataqc_localrangetest` |
| 2 | `dataqc_spiketest` |
| 3 | `dataqc_polytrendtest` |
| 4 | `dataqc_stuckvaluetest` |
| 5 | `dataqc_gradienttest` |
| 6 | Not assigned |
| 7 | `dataqc_propagateflags` |

See OOI, Interpreting QC Variables and Results, in the [References](#references) for more on reading these variables.

### Tests Applied by the OOI Data System

The OOI stream engine calls the wrapper functions (not the core functions), passing the test parameters for each data
product from the QC lookup tables maintained in the
[oceanobservatories/qc-lookup](https://github.com/oceanobservatories/qc-lookup) repository. A test runs for a data
product only when the lookup tables supply parameters for it; the `qc_executed` bits show which tests were run on a
given data set.

- The lookup tables define parameters for the global range, spike, stuck value, and trend tests for most data
  products, and for the gradient test for a small subset.
- The local range test and the combined QC flags are not applied. The local range test wrapper requires a callback
  argument (`pval_callback`) that the stream engine does not supply, and there is no lookup table for the combined
  flags, whose wrapper expects flag vectors as positional arguments.

### Compiled Routines

The stuck value, spike, and gradient tests run in compiled C code (`extensions/stuck.c`, `spike.c`, and
`gradient.c`), called through the Cython module `ion_functions/qc/qc_extensions.pyx` as `stuckvalues`, `spikevalues`,
and `gradientvalues`. The same module provides `ntp_to_month`, which the local range test wrapper uses to convert NTP
time stamps into the month of the year. Because these routines are compiled, they are described here rather than
listed as function entries below; `qc_functions.py` requires the compiled `qc_extensions` module to import.

### GLBLRNG -- Global Range Test

The global range test flags data points (1 or 0) according to whether they fall within a fixed valid range, thought
of as universally valid -- hence "global". Data flagged 0 should be assumed to be of bad quality and not used for
scientific analysis.

### LOCLRNG -- Local Range Test

The local range test uses valid ranges that are site-specific -- hence "local" -- and can vary as a function of one or
more independent coordinates, such as the depth of the measurement or the month of the year. The minimum and maximum
limits, given in a lookup table at fixed coordinates, are linearly interpolated to the coordinates of each data point:
with `numpy.interp` for one coordinate, or with `scipy.interpolate.LinearNDInterpolator` (Delaunay triangulation) for
more than one. Data points outside the coordinates covered by the table are flagged 0. The validity of the output is
limited by the quality of the ranges in the lookup table.

### SPKETST -- Spike Test

The spike test flags single values that deviate significantly from the surrounding values. Each value is compared with
its $L-1$ neighbors within a window of length $L$ (an odd integer). With $R$ the range (maximum minus minimum) of the
neighbors and $ACC$ the measurement accuracy, the value is presumed good if it deviates from the mean of the neighbors
by less than

$$N \times \max(R, ACC)$$

Away from the ends of the series the neighbors lie symmetrically before and after the test value; within $(L-1)/2$
values of either end, the neighbors are the first or last $L$ values. $ACC$ restricts spike detection to deviations
larger than $N \times ACC$ even when the data have little variability; setting it to 0 disables this. The test can only
detect single outliers, and it assumes regularly spaced measurements, although OOI runs it on data with gaps or
changes in sampling.

### TRNDTST -- Trend Test

The trend test checks whether a significant fraction of the variability in a time series can be explained by a
polynomial drift, possibly a sensor drift. A polynomial of order $n$ (1 for a linear drift) is fitted to the series;
the series is taken as containing a significant trend if

$$nstd \times \sigma(x - P) < \sigma(x)$$

where $P$ is the fitted polynomial, $\sigma$ the standard deviation, and $nstd$ a threshold factor. The result applies
to the whole series, so every point receives the same flag (0 if a trend is detected). The DPS specifies that the input
is aggregated over a measurement interval set in the lookup table.

### STUCKVL -- Stuck Value Test

The stuck value test flags repeated occurrences of one value. When at least $N$ successive values lie within the
resolution $R$ of each other, all of them are flagged 0. Like the spike test, it assumes regularly spaced measurements
but is run on data with gaps or changes in sampling.

### GRADTST -- Gradient Test

The gradient test scans the data point by point. Starting from a point assumed to be good, the next point is good if
the change between the two, divided by the change in the coordinate $x$ (e.g. time or pressure), lies within the valid
range `ddatdx`; otherwise it is bad. Once a point is bad, the following points stay bad until one returns to within
`toldat` of the last good value. By default the first point is presumed good; a start value (`startdat`) can be given
instead, in which case the first point within `toldat` of it is the first good point. Because a large change between
two points does not show which of them is wrong, the result depends on this starting assumption, and a natural trend
during a stretch of bad data can prevent a return to good values. Points closer together in $x$ than `mindx` can be
excluded to avoid exaggerated gradients from small coordinate steps. The implementation differs from the DPS example
code in several respects; see the Additional Notes for `dataqc_gradienttest`.

### CMBNFLG -- Combined QC Flags

The combined flag merges several sets of flags into one with a logical "and": the result is 1 only where every input
flag is 1. One use is propagating the flags of input parameters (e.g. conductivity, temperature, and pressure) to a
derived product (e.g. salinity); another is merging the results of several tests on a single data product.

### Specifications Without Functions

DPS 1341-10001 (MODULUS), 1341-10002 (INTERP1), and 1341-10003 (POLYVAL) specify single-step operations used to apply
calibrations or adjustments: wrapping values such as compass directions back into a nominal range, one-dimensional
linear interpolation, and evaluating a polynomial. They correspond to the NumPy functions `numpy.mod`,
`numpy.interp`, and `numpy.polyval` and are not implemented as functions in this module.

Full algorithm derivations and source references are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.qc.qc_functions.dataqc_globalrangetest

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-06 | Christopher Wingard | Initial code |
| 2013-05-30 | Christopher Mueller | Performance improvements; added strict_validation flag |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_localrangetest

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-06 | Christopher Wingard | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_spiketest

#### Additional Notes

If the input has fewer values than the window length $L$, the test is not applied and every point is flagged 1.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-06 | Christopher Wingard | Initial code |
| 2013-05-30 | Christopher Mueller | Performance improvements |
| 2013-06-03 | Luke Campbell | Moved the test into compiled C code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_polytrendtest

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-06 | Christopher Wingard | Initial code |
| 2013-05-30 | Christopher Mueller | Performance improvements |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_stuckvaluetest

#### Additional Notes

As in the DPS example code, every point is flagged 0 when the input has fewer values than `num`. Within each window
of `num` values, the compiled code compares each value with the last value of the window, where the DPS example code
compares with the first; a `num` of 0 is treated as the default of 10.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-06 | Christopher Wingard | Initial code |
| 2013-05-31 | Luke Campbell | Moved the test into compiled C code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_gradienttest

#### Additional Notes

The implementation differs from the DPS example code:

| Behavior | DPS 1341-10010 example code | This implementation |
|---|---|---|
| Points closer than `mindx` | Removed; the output is shorter than the input | Kept and flagged -99, so the output matches the input length. If the next gradient fails, the skipped points are set to 0. |
| No start value | NaN or empty `startdat` presumes the first point good | NaN is converted to 0, and 0 means "not set"; a start value of exactly 0 cannot be used. |
| NaN data values | Removed before testing | Not removed |
| Too few points | Warns when one point or fewer remains | Returns 1 for every point, as floating-point values, when the full range of `x` is less than `mindx` |

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-06 | Christopher Wingard | Initial code |
| 2013-06-05 | Luke Campbell | Moved the test into compiled C code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_propagateflags

#### Additional Notes

An `inflags` array with fewer than two rows raises a `NameError` (the code calls an undefined `error` function) rather
than a descriptive error.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-06 | Christopher Wingard | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

The OOI data system calls these wrappers. Each extracts scalar test parameters from the arrays supplied by the system
(using the last element), returns -99 for every point if a parameter is missing or holds the -9999 fill value, and
otherwise calls the corresponding core function. External users should call the core functions directly.

::: ion_functions.qc.qc_functions.dataqc_globalrangetest_minmax

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-02 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_localrangetest_wrapper

#### Additional Notes

The second fill-value check tests `datlim` again where it should test `datlimz`, so a -9999 fill value in `datlimz`
is not detected.

The OOI stream engine does not supply `pval_callback`, so this wrapper is not run in the OOI data system.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-14 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_spiketest_wrapper

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-05 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_polytrendtest_wrapper

#### Additional Notes

The fill-value check tests `ord_n` twice and never tests `nstd`, so a -9999 fill value in `nstd` is not detected.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-06 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_stuckvaluetest_wrapper

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-05 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_gradienttest_wrapper

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-05 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.dataqc_propagateflags_wrapper

#### Additional Notes

The check that the flag vectors have equal lengths runs only when `strict_validation` is False.

The OOI stream engine passes only keyword arguments and has no lookup table for this test, so this wrapper is not run
in the OOI data system.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-04 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.qc.qc_functions.is_fill

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-05 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.qc.qc_functions.is_none

#### History
| Date | Author | Change |
|---|---|---|
| 2013-06-05 | Luke Campbell | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Deprecated Functions

!!! note "Deprecated"
    `dataqc_solarelevation` and `dataqc_condcompress` were never used in the OOI data system and are deprecated.

::: ion_functions.qc.qc_functions.dataqc_solarelevation

#### Additional Notes

Adapted from `soradna1.m` in the Air-Sea Toolbox, which takes its expressions from Appendix E of the 1978 Almanac for
Computers (U.S. Naval Observatory) and its solar constant (1368 W m$^{-2}$) from Coffey et al. (1995). The inputs must
be array-like: `lon`, `lat`, and `dt` are measured with `len()`, so scalar inputs raise a `TypeError`.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-07 | Christopher Wingard | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation; decorated with `@deprecated` (never used) |

---

::: ion_functions.qc.qc_functions.dataqc_condcompress

#### Additional Notes

Implements the compressibility compensation from Sea-Bird Application Note 10, as given in DPS 1341-10030:

$$c_{new} = c_{orig} \times \frac{1 + CPcor \times p_{orig}}{1 + CPcor \times p_{new}}$$

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-07 | Christopher Wingard | Initial code |
| 2026-10-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation; decorated with `@deprecated` (never used) |

## References

Coffey et al. (1995). Earth System Monitor, 6, 6-10.

OOI. [Interpreting QC Variables and Results.](https://oceanobservatories.org/knowledgebase/interpreting-qc-variables-and-results/)

[OOI (2012). Data Product Specification for Modulus (MODULUS). Document Control Number 1341-10001.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10001_Data_Product_SPEC_MODULUS_OOI.pdf)

[OOI (2012). Data Product Specification for 1-D Interpolation (INTERP1). Document Control Number 1341-10002.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10002_Data_Product_SPEC_INTERP1_OOI.pdf)

[OOI (2012). Data Product Specification for Evaluate Polynomial (POLYVAL). Document Control Number 1341-10003.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10003_Data_Product_SPEC_POLYVAL_OOI.pdf)

[OOI (2012). Data Product Specification for Global Range Test. Document Control Number 1341-10004.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10004_Data_Product_SPEC_GLBLRNG_OOI.pdf)

[OOI (2012). Data Product Specification for Local Range Test. Document Control Number 1341-10005.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10005_Data_Product_SPEC_LOCLRNG_OOI.pdf)

[OOI (2012). Data Product Specification for Spike Test. Document Control Number 1341-10006.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10006_Data_Product_SPEC_SPKETST_OOI.pdf)

[OOI (2012). Data Product Specification for Trend Test. Document Control Number 1341-10007.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10007_Data_Product_SPEC_TRNDTST_OOI.pdf)

[OOI (2012). Data Product Specification for Stuck Value Test. Document Control Number 1341-10008.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10008_Data_Product_SPEC_STUCKVL_OOI.pdf)

[OOI (2012). Data Product Specification for Gradient Test. Document Control Number 1341-10010.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10010_Data_Product_SPEC_GRADTST_OOI.pdf)

[OOI (2012). Data Product Specification for Solar Elevation. Document Control Number 1341-10011.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10011_Data_Product_SPEC_SOLAREL_OOI.pdf)

[OOI (2012). Data Product Specification for Combined QC Flags. Document Control Number 1341-10012.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10012_Data_Product_SPEC_CMBNFLG_OOI.pdf)

[OOI (2012). Data Product Specification for Conductivity Compressibility Compensation. Document Control Number 1341-10030.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-10030_Data_Product_SPEC_CONDCMP_OOI.pdf)

Sea-Bird Scientific. *Compressibility Compensation of Sea-Bird Conductivity Sensors*. Application Note 10.
Bellevue, WA: Sea-Bird Scientific.

U.S. Naval Observatory, Nautical Almanac Office (1978). Almanac for Computers.
