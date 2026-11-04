# Generic Functions

## Background

The Generic Functions module has no dedicated instrument family or Data Product
Specification. Rather, this module holds functions needed by more than one
instrument family's module -- fill-value handling, magnetic declination and
velocity correction, array utilities, and small comparison helpers used by
several data products across the codebase.

Of particular importance are the magnetic correction algorithms. Magnetic declination is computed from the
International Geomagnetic Reference Field (IGRF), 14th generation, via the `ppigrf` Python package.
`magnetic_declination` is the array entry point used throughout the codebase. Because `ppigrf` evaluates every
combination of the dates and positions it is given, the function first rounds the samples (time to the day, evaluated
at mid-day; position to four decimal places, about 10 m), removes duplicates, and evaluates the IGRF once on the grid of
unique days by unique positions, in chunks to bound memory, before mapping results back to the original samples. The
rounding adds at most about 0.0003 degrees of error, far below the uncertainty of the IGRF itself. `igrf_declination`
computes the exact, unrounded declination for a single location and date and serves as the reference implementation;
it rejects array inputs.

`magnetic_correction` rotates horizontal velocity components from magnetic to true Earth coordinates using one
declination value per sample. Velocities may be one value per sample (e.g. VELPTMN, WINDAVG) or a profile per sample.
`adcp_functions.magnetic_correction` and `wav_functions.magnetic_correction_einsum` implement the same rotation
separately; both obtain declination from `magnetic_declination`.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.generic_functions.magnetic_declination

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |
| 2026-10-16 | Christopher Wingard | Replaced `np.vectorize` over `igrf_declination` with rounded, de-duplicated grid evaluation of the IGRF. |

---

::: ion_functions.data.generic_functions.magnetic_correction

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |
| 2026-10-16 | Christopher Wingard | Accepts array inputs directly, one declination per sample. |

---

::: ion_functions.data.generic_functions.igrf_declination

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; migrated from pyIGRF to ppigrf (IGRF-14 model); fixed z/zflag sign-check bug; corrected docstring Example, which had `z` and `timestamp` swapped in the call. |
| 2026-10-16 | Christopher Wingard | Restricted to single-valued inputs; raises `ValueError` for arrays. |

---

## Helper Functions

::: ion_functions.data.generic_functions.replace_fill_with_nan

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.extract_parameter

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.select_non_zero_arg

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.select_arg_within_tolerance_of_std

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.error

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; added a docstring (previously undocumented). |

---

::: ion_functions.data.generic_functions.bilinear_interpolation

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

## Deprecated Functions

!!! note "Deprecated"
    `ntp_to_unix_time` is deprecated. Use the fixed offset of 2208988800
    seconds between the NTP and Unix epochs directly instead of calling
    this function.

::: ion_functions.data.generic_functions.ntp_to_unix_time

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; decorated with `@deprecated` and added to `DEPRECATIONS.md` (previously marked deprecated in docstring prose only). |

---

## References

Beggan, C.D., et al. (2026). International Geomagnetic Reference Field: the fourteenth
generation. *Earth Planets Space*, 78, 127. <https://doi.org/10.1186/s40623-025-02360-0>

International Association of Geomagnetism and Aeronomy (2024). [IGRF-14.
Zenodo.](https://doi.org/10.5281/zenodo.14012302)

OOI (2020). [Data Product Specification for Velocity Profile and Echo
Intensity. Document Control Number
1341-00750.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00750_Data_Product_SPEC_VELPROF_ECHOINT_OOI.pdf)

OOI (2013). [Data Product Specification for Turbulent Velocity Profile and
Echo Intensity. Document Control Number
1341-00760.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00760_Data_Product_SPEC_VELTURB_ECHOINT_OOI.pdf)

Strom, K.M., and Reistad, H. (2024). [ppigrf: Python package for computing
the International Geomagnetic Reference Field
(IGRF).](https://github.com/IAGA-VMOD/ppigrf)
