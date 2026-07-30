# Generic Functions

## Background

`generic_functions.py` has no dedicated instrument family or Data Product
Specification of its own. It holds functions needed by more than one
instrument family's module -- fill-value handling, magnetic declination and
velocity correction, array utilities, and small comparison helpers used by
several data products across the codebase.

Magnetic declination is computed from the International Geomagnetic
Reference Field (IGRF), 14th generation, via the `ppigrf` Python package.
`igrf_declination` computes declination for a single location and date;
`magnetic_declination` vectorizes it over arrays of samples. `magnetic_correction`
applies that declination to rotate a velocity vector from magnetic to true
compass coordinates, and is shared by several data products (e.g. VELPROF,
WINDAVG) across multiple instrument classes -- it is not used for the ADCP's
velocity profiles, which have their own implementation in
`adcp_functions.adcp_magvar`.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.generic_functions.replace_fill_with_nan

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.magnetic_declination

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.magnetic_correction

#### Additional Notes

Several callers wrap `magnetic_correction` with `np.vectorize` to apply it
across arrays of samples (e.g. `met_functions.met_wind_mag_corr`). On the
current numpy version, this raises `ValueError: setting an array element
with a sequence`, because `magnetic_correction` returns `cor[0], cor[1]`,
each a length-1 array rather than a true scalar. This is a pre-existing
issue, not introduced by any change on this page, and is not fixed here;
see the ion-functions `CLAUDE.md` Known Issues section.

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.igrf_declination

#### Additional Notes

Prior to 2026-07-30 this function used the `pyIGRF` package, which is
unmaintained and is not installed in the current `ion` conda environment at
all -- the module could not be imported. It was rewritten to use `ppigrf`,
which defaults to the IGRF-14 coefficient set (valid 1900-01-01 to
2030-01-01). The rewrite also fixed a pre-existing sign-check bug: the
depth/height flip (`if z > 0 & zflag == -1`) used the bitwise `&` operator,
which binds tighter than the comparison operators here, so the condition
was always `False` and the sign flip never actually ran regardless of
input.

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; migrated from pyIGRF to ppigrf (IGRF-14 model); fixed z/zflag sign-check bug; corrected docstring Example, which had `z` and `timestamp` swapped in the call. |

---

::: ion_functions.data.generic_functions.extract_parameter

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.select_non_zero_arg

#### Additional Notes

Not called elsewhere in `ion_functions` and has no test coverage in this
repository; called directly by the external CI stream engine.

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

::: ion_functions.data.generic_functions.select_arg_within_tolerance_of_std

#### Additional Notes

Not called elsewhere in `ion_functions` and has no test coverage in this
repository; called directly by the external CI stream engine.

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---

## Helper Functions

::: ion_functions.data.generic_functions.error

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-30 | Christopher Wingard | Converted to NumPy docstring format; added a docstring (previously undocumented). |

---

## Utility Functions

Functions in this section do not fit the Core/Helper/Wrapper classification
in [Function Types](../function_types.md) -- they are not called by any
other function in this module or elsewhere in `ion_functions`, and have no
known external caller.

::: ion_functions.data.generic_functions.bilinear_interpolation

#### Additional Notes

Not called elsewhere in `ion_functions` and has no test coverage; no known
external caller has been identified either.

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
