# CO2 Functions

## Background

The Ocean Observatories Initiative deploys two families of CO2 instruments
across its surface and subsurface platforms. The table below lists the OOI
instrument classes covered by this module.

| Class | Hardware | Platform | Designator meaning |
|---|---|---|--------------------|
| PCO2W | Sunburst SAMI2-CO2 | Moored (fixed depth) | pCO2, Seawater     |
| PCO2A | Pro-Oceanus CO2-Pro | Surface buoys | pCO2, Air-Sea      |

`co2_functions.py` processes data from both instrument classes and computes
a third derived product using auxiliary bulk meteorology data from the METBK
instrument class.

---

### PCO2WAT — Partial Pressure of CO2 in Seawater (L0 to L1)

PCO2WAT is the partial pressure of CO2 ($\text{pCO}_2$) in seawater at
depth, computed from the Sunburst SAMI2-CO2 (PCO2W). The SAMI2-CO2 uses a
bromothymol blue (BTB) pH indicator solution that equilibrates with the
ambient seawater through a permeable silicone membrane. The equilibrated
indicator solution passes through an optical cell where absorbance is
measured at 434 nm and 620 nm — the peak absorbance wavelengths for the
protonated and deprotonated forms of the indicator. Periodic blank
measurements (every 3.5 days) correct for drift in the electro-optical
system (DPS 1341-00490, §3.1).

#### L0 Inputs and Blank Normalization

The SAMI2-CO2 outputs two record types in hexadecimal serial format:

- **Type 4** — regular measurement (used to compute PCO2WAT_L1)
- **Type 5** — blank measurement (used to correct Type 4 records)

Each record contains a 434 nm ratio (`Ratio434`, column index 6 of the
light array), a 620 nm ratio (`Ratio620`, column index 7), and a raw
thermistor count (`Therm`). Raw blank counts from Type 5 records are
normalized to dimensionless ratios by dividing by 16384 ($2^{14}$):

$$A_{434,\text{blank}} = \frac{\text{Ratio434}_\text{blank}}{16384}$$

$$A_{620,\text{blank}} = \frac{\text{Ratio620}_\text{blank}}{16384}$$

#### Thermistor Temperature

The raw thermistor count is converted to degrees C. The conversion depends
on the SAMI hardware generation (12-bit or 14-bit ADC). For 12-bit hardware
(full-scale count $= 4096$):

$$r_t = \ln\!\left(\frac{\text{Therm}}{4096 - \text{Therm}} \times 17400\right)$$

For 14-bit hardware (full-scale count $= 16384$):

$$r_t = \ln\!\left(\frac{\text{Therm}}{16384 - \text{Therm}} \times 17400\right)$$

In both cases temperature is then:

$$T_C = \frac{1}{0.0010183 + 0.000241\, r_t + 1.5 \times 10^{-7}\, r_t^3} - 273.15$$

#### PCO2WAT_L1 Calculation

The L1 pCO2 is computed from the blank-corrected absorbances and a
temperature-corrected SAMI response. The blank-corrected absorbances at
each wavelength are:

$$A_{434} = -\log_{10}\!\left(\frac{\text{Ratio434}}{A_{434,\text{blank}}}\right)$$

$$A_{620} = -\log_{10}\!\left(\frac{\text{Ratio620}}{A_{620,\text{blank}}}\right)$$

The absorbance ratio is:

$$R = \frac{A_{620}}{A_{434}}$$

The SAMI response $\text{RCO2}$ and intermediate temperature-correction
terms are:

$$\text{RCO2} = -\log_{10}\!\left(\frac{R - e_1}{e_2 - e_3 R}\right)$$

$$\text{RCO2}_2 = \text{RCO2} + 0.007\,(T_C - \text{CalT})$$

$$T_\text{coeff} = 0.0075778 + 0.0012389\,\text{RCO2}_2 - 0.00048757\,\text{RCO2}_2^2$$

$$\text{Tcor\_RCO2} = \text{RCO2} + T_\text{coeff}\,(T_C - \text{CalT})$$

The final pCO2 is obtained from the quadratic calibration curve:

$$\text{pCO2} = 10^{\,\dfrac{-\text{CalB} + \sqrt{\text{CalB}^2 - 4\,\text{CalA}\,(\text{CalC} - \text{Tcor\_RCO2})}}{2\,\text{CalA}}}$$

The equilibration constants $e_1 = 0.0043$, $e_2 = 2.136$, $e_3 = 0.2105$
are fixed values provided by Sunburst Sensors based on laboratory
determinations with the BTB indicator batch and are not recalculated per
calibration cycle (DPS 1341-00490, §2.5.3).

The instrument-specific calibration coefficients CalT, CalA, CalB, and CalC
are provided by the manufacturer after each calibration or refurbishment
cycle and are stored as metadata (DPS 1341-00490, §4.2). Instrument results
are valid between 0 and 35 °C (DPS 1341-00490, §3.3).

**Algorithm note:** The Python implementation corrects errors present in the
vendor-supplied Matlab example code distributed with DPS 1341-00490
(Appendix A, 2018).

Output accuracy: $\pm 3\,\mu\text{atm}$; precision $<1\,\mu\text{atm}$
(DPS 1341-00490, Appendix B).

---

### PCO2ATM and PCO2SSW — pCO2 in Air and Surface Seawater (L0 to L1)

PCO2ATM and PCO2SSW are the partial pressures of $\text{CO}_2$ in
atmosphere and surface seawater, respectively. Both are computed from the
Pro-Oceanus CO2-Pro (PCO2A), which determines the CO2 mole fraction
($x_{\text{CO2}}$, in ppm) internally by measuring the infrared absorbance
of $\text{CO}_2$ and compensating for pressure, temperature, and humidity
using onboard firmware (DPS 1341-00260, §3.1). The instrument alternates
between air sampling (producing XCO2ATM_L0, labeled "A" in the serial
output) and water sampling (producing XCO2SSW_L0, labeled "W").

The L0 CO2 mole fraction and the L0 internal gas stream pressure PRESAIR
are converted to $\mu\text{atm}$ using:

$$\text{pCO2} = \frac{x_{\text{CO2}} \times p}{\text{STD}}$$

where $p$ is PRESAIR in mbar and $\text{STD} = 1013.25$ mbar/atm is
standard atmospheric pressure. Because $x_{\text{CO2}}$ is in ppm
($10^{-6}$), the result is directly in $\mu\text{atm}$ without further
scaling (DPS 1341-00260, §3.2).

The same function `pco2_ppressure` computes both PCO2ATM_L1 and PCO2SSW_L1;
the distinction between air and surface seawater measurements is determined
by the L0 input supplied (XCO2ATM or XCO2SSW) and the record type flag in
the serial data stream.

Output accuracy: $\pm 1\,\mu\text{atm}$; precision $\pm 0.01\,\mu\text{atm}$
(DPS 1341-00260, Appendix B).

---

### CO2FLUX — Air-Sea Flux of CO2 (L1 to L2)

CO2FLUX is the OOI Level 2 estimate of the flux of $\text{CO}_2$ across
the air-sea interface. It is computed from L1 PCO2SSW and PCO2ATM (from
PCO2A) and L1/L2 bulk meteorology inputs (from METBK): sea surface
temperature (TEMPSRF_L1), sea surface salinity (SALSURF_L2), and wind
speed at 10 m height (WIND10M_L2).

The flux is positive when directed from the ocean to the atmosphere
(DPS 1341-00270, §2.2.7). The computation follows these steps:

**Step 1** — Convert pCO2 from $\mu\text{atm}$ to atm.

**Step 2** — Compute the Schmidt number $S_c$ from sea surface temperature
$t$ ($^\circ\text{C}$) (Wanninkhof, 1992, Table A1):

$$S_c = 2073.1 - 125.62\,t + 3.6276\,t^2 - 0.043219\,t^3$$

**Step 3** — Compute the gas transfer velocity $k$ in cm h$^{-1}$ and
convert to m s$^{-1}$ (Sweeney et al., 2007, Fig. 3 and Table 1):

$$k = \frac{0.27\, u_{10}^2 \sqrt{660 / S_c}}{100 \times 3600}$$

**Step 4** — Compute absolute temperature $T = t + 273.15$ (K), then
compute CO2 solubility $K_0$ in mol atm$^{-1}$ m$^{-3}$ using the volume
formulation (Weiss, 1974, Eqn. 12 and Table I):

$$K_0 = 1000 \exp\!\left(-58.0931 + \frac{90.5069 \times 100}{T} + 22.2940\,\ln\!\frac{T}{100}\right.$$
$$\left.+ S\!\left(0.027766 - 0.025888\,\frac{T}{100} + 0.0050578\,\left(\frac{T}{100}\right)^{\!2}\right)\right)$$

**Step 5** — Compute the flux (Wanninkhof, 1992, Eqn. A2):

$$F = k\, K_0\,(\text{pCO2}_w - \text{pCO2}_a)$$

The inherent uncertainty of the flux estimate is approximately 10%
(DPS 1341-00270, Appendix B).

---

## Core functions

::: ion_functions.data.co2_functions.pco2_abs434_ratio

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_abs620_ratio

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_blank

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2014-02-28 | Christopher Wingard | Updated to accept raw blank values from a sparse array |
| 2018-03-04 | Christopher Wingard | Updated blank normalization per vendor-corrected algorithm |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_thermistor

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2023-01-12 | Mark Steiner | Added sami_bits argument to handle 14-bit hardware |
| 2023-08-15 | Samuel Dahlberg | Renamed local variables; replaced numexpr with numpy |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_battery

#### History
| Date | Author | Change |
|---|---|---|
| 2023-02-23 | Mark Steiner | Initial code |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_calc_pco2

#### Additional Notes
The `ea434`, `eb434`, `ea620`, and `eb620` calibration coefficient
parameters are accepted as inputs for backward compatibility but are
not used in the calculation. The original formulation derived $e_1$,
$e_2$, and $e_3$ from these per-instrument values; the current
algorithm uses fixed constants supplied by Sunburst Sensors
($e_1 = 0.0043$, $e_2 = 2.136$, $e_3 = 0.2105$) for all instruments.
The retained parameters and the commented-out original formulation in
the code document this evolution.

#### History
| Date | Author | Change |
|---|---|---|
| unknown | J. Newton (Sunburst Sensors, LLC) | Original Matlab code |
| 2013-04-20 | Christopher Wingard | Initial Python code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2014-03-19 | Christopher Wingard | Optimized |
| 2018-03-04 | Christopher Wingard | Corrected blank normalization and temperature correction per vendor-revised algorithm |
| 2023-01-12 | Mark Steiner | Changed therm argument to degrees C input |
| 2023-08-15 | Samuel Dahlberg | Renamed local variables to follow naming convention |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_ppressure

#### History
| Date | Author | Change |
|---|---|---|
| 2014-10-27 | Christopher Wingard | Initial Python code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_co2flux

#### History
| Date | Author | Change |
|---|---|---|
| 2012-03-28 | Matthias Lankhorst | Original Matlab code |
| 2013-04-20 | Christopher Wingard | Initial Python code |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

`pco2_pco2wat` is the OOI single-output wrapper for the PCO2WAT_L1 data
product. It calls `pco2_calc_pco2` for all records and then sets blank
measurement records (mtype == 5) to the system fill value, satisfying the
OOI single-output data product requirement. External users should call
`pco2_calc_pco2` directly.

::: ion_functions.data.co2_functions.pco2_pco2wat

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2014-03-19 | Christopher Wingard | Optimized |
| 2017-04-04 | Pete Cable | Updated to use thermistor/blank counts per DPS |
| 2025-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

DeGrandpre, M. D., Baehr, M. M., and Hammar, T. R. (1999). Calibration-free
optical chemical sensors. Analytical Chemistry, 71(6), 1152-1159.

Sweeney, C., Gloor, E., Jacobson, A. R., Key, R. M., McKinley, G.,
Sarmiento, J. L., and Wanninkhof, R. (2007). Constraining global air-sea
gas exchange for CO2 with recent bomb 14C measurements. Global
Biogeochemical Cycles, 21, GB2015.

Wanninkhof, R. (1992). Relationship between wind speed and gas exchange
over the ocean. Journal of Geophysical Research, 97(C5), 7373-7382.

Weiss, R. F. (1974). Carbon dioxide in water and seawater: the solubility
of a non-ideal gas. Marine Chemistry, 2, 203-215.

[OOI (2012). Data Product Specification for Partial Pressure of CO2 in Air
and Surface Seawater. Document Control Number 1341-00260.](https://oceanobservatories.org/wp-content/uploads/2015/09/1341-00260_Data_Product_SPEC_PCO2ATM_PCO2SSW_OOI.pdf)

[OOI (2018). Data Product Specification for Partial Pressure of CO2 in
Seawater. Document Control Number 1341-00490.](https://oceanobservatories.org/wp-content/uploads/2015/09/1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)

[OOI (2012). Data Product Specification for Flux of CO2 from the Ocean into
the Atmosphere. Document Control Number 1341-00270.](https://oceanobservatories.org/wp-content/uploads/2015/09/1341-00270_Data_Product_SPEC_CO2FLUX_OOI.pdf)
