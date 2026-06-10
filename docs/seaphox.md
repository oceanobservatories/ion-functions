# Deep SeapHOx V2

The Sea-Bird Scientific Deep SeapHOx V2 is an integrated moored instrument
that combines three sensors into a single, integrated instrument: an SBE
37 MicroCAT CTD, an SBE 63 optical dissolved oxygen sensor, and a Deep
SeaFET V2 ISFET pH sensor. OOI deploys the Deep SeapHOx V2 under two
instrument class designators that differ only in pressure rating:

| OOI Class | Pressure Rating | Depth Limit |
|---|---|---|
| PHSEN-G | 350 m | Shallow moored platforms |
| PHSEN-H | 2000 m | Deep moored platforms |

Both classes produce the same data products using the same processing
functions.

## Instrument Architecture

The Deep SeapHOx V2 is described by Sea-Bird Scientific as an [Ocean
CT(D)-pH-DO Sensor](https://www.seabird.com/products/deep-seaphox-v2-ph-sensor). 
From the user perspective it behaves as a single instrument, but internally it 
comprises three sensor subsystems that each require their own calibration coefficients 
and processing functions:

| Sub-instrument | Measurements | OOI Data Products                            |
|---|---|----------------------------------------------|
| SBE 37 MicroCAT | Temperature, pressure, conductivity | TEMPWAT, PRESWAT, CONDWAT, PRACSAL, DENSITY  |
| SBE 63 Optical DO | Dissolved oxygen, thermistor temperature | DOCONCS, DOXYGEN |
| Deep SeaFET V2 | ISFET voltage | PHWATER                                      |

The SBE 37 provides the temperature, pressure, and salinity context required
by both the SBE 63 oxygen algorithm and the Deep SeaFET pH algorithm. The
SBE 63 uses its own internal thermistor as the primary temperature input for
the oxygen calculation; the SBE 37 provides the salinity and pressure
corrections. The Deep SeaFET uses the SBE 37 temperature, pressure, and
salinity as inputs to the pH algorithm.

## Processing Functions by Sub-instrument

Because the Deep SeapHOx V2 contains three sensor types, its processing
functions are split across three instrument family documentation pages
rather than collected on a single page (although they are all collected in a 
single python module, `phsen_h_functions.py`). Each sub-instrument's functions 
are documented with the family they belong to scientifically:

**SBE 37 MicroCAT (CTD)** — Raw A/D counts to calibrated temperature,
pressure, and conductivity. The SBE 37 uses the same calibration polynomial
families as other Sea-Bird MicroCAT instruments. See
[CTD Functions — SeapHOx CTD Functions](api/ctd_functions.md#seaphox-ctd-functions).

Functions: `temperature_raw_conversion`, `pressure_raw_conversion`,
`conductivity_raw_conversion` (all in `phsen_h_functions.py`).

**SBE 63 Optical DO Sensor** — Raw phase and thermistor counts to dissolved
oxygen concentration in $\mu$mol kg$^{-1}$. Uses a Stern-Volmer fluorescence
quenching model with salinity and pressure corrections. See
[DO Functions — PHSEN-H Dissolved Oxygen](api/do2_functions.md#phsen-h-dissolved-oxygen-sbe-63-optical-sensor).

Functions: `dissolved_oxygen`, `convert_sbe63_thermistor`
(both in `phsen_h_functions.py`).

**Deep SeaFET V2 ISFET (pH)** — Raw ISFET voltage counts to pH on the total
hydrogen ion scale (PHWATER_L2). See
[pH Functions — Sea-Bird Scientific Deep SeapHOx V2](api/ph_functions.md#sea-bird-scientific-deep-seaphox-v2).

Functions: `convert_ph_voltage_counts`, `ph_total`
(both in `phsen_h_functions.py`).

## Internal Sensor Functions

The Deep SeapHOx V2 housing includes an internal SHT-series sensor that
monitors the instrument's internal environment. These functions process that
sensor's raw output and are not associated with any OOI science data product.

::: ion_functions.data.phsen_h_functions.internal_temperature

#### History
| Date | Author | Change |
|---|---|---|
| 2025-05-15 | Christopher Wingard | Initial NumPy docstring |

---

::: ion_functions.data.phsen_h_functions.internal_humidity

#### History
| Date | Author | Change |
|---|---|---|
| 2025-05-15 | Christopher Wingard | Initial NumPy docstring |

---

## References

[Sea-Bird Scientific. *Deep SeapHOx V2. pH, Conductivity, Temperature, Pressure, Dissolved
Oxygen*. User Manual, 2025-04-29, Rev E. Bellevue, WA.](https://www.seabird.com/products/deep-seaphox-v2-ph-sensor)
