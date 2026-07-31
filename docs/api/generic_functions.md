# Generic Functions

## Background

The Generic Functions module has no dedicated instrument family or Data Product
Specification. Rather, this module holds functions needed by more than one
instrument family's module -- fill-value handling, magnetic declination and
velocity correction, array utilities, and small comparison helpers used by
several data products across the codebase.

Of particular importance are the magnetic correction algorithms. Magnetic 
declination is computed from the International Geomagnetic Reference Field 
(IGRF), 14th generation, via the `ppigrf` Python package. `igrf_declination` 
computes declination for a single location and date; `magnetic_declination` 
vectorizes it over arrays of samples. `magnetic_correction` applies that 
declination to rotate a velocity vector from magnetic to true compass coordinates, 
and is shared by several data products (e.g. VELPROF, WINDAVG) across multiple 
instrument classes -- it is not used for the ADCP's
velocity profiles, which have their own implementation in 
`adcp_functions.magnetic_correction` (adjusted to work with velocity profiles, but
still uses the `magnetic_declination` function defined here).

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.generic_functions.igrf_declination

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; migrated from pyIGRF to ppigrf (IGRF-14 model); fixed z/zflag sign-check bug; corrected docstring Example, which had `z` and `timestamp` swapped in the call. |

---

::: ion_functions.data.generic_functions.magnetic_declination

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.magnetic_correction

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

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

International Association of Geomagnetism and Aeronomy (2024). [IGRF-14.
Zenodo.](https://doi.org/10.5281/zenodo.14012302)

OOI (2012). [Data Product Specification for Velocity Profile and Echo
Intensity. Document Control Number
1341-00750.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

OOI (2013). [Data Product Specification for Turbulent Velocity Profile and
Echo Intensity. Document Control Number
1341-00760.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00760_Data_Product_SPEC_VELPROF_OOI.pdf)

Strom, K.M., and Reistad, H. (2024). [ppigrf: Python package for computing
the International Geomagnetic Reference Field
(IGRF).](https://github.com/IAGA-VMOD/ppigrf)
