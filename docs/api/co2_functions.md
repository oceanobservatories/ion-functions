# CO2 Functions

## Background

The Ocean Observatories Initiative deploys two families of CO$_2$ instruments
across its surface and subsurface platforms. The table below lists the OOI
instrument classes covered by this module.

| Class | Hardware | Platform                           | Designator meaning |
|---|---|------------------------------------|--------------------|
| PCO2W | Sunburst SAMI-CO2 | Moored (fixed depth) and Profiling | pCO$_2$, Seawater  |
| PCO2A | Pro-Oceanus CO2-Pro | Surface buoys                      | pCO$_2$, Air-Sea   |

`co2_functions.py` processes data from both instrument classes and computes
a third derived product, CO$_2$ flux from the ocean to the atmosphere 
(CO2FLUX_L2), using auxiliary bulk meteorology data from the METBK instrument 
class. Within the OOI data system, PCO2W instruments fall under the Water
Column science regime and the Dissolved CO$_2$ category. PCO2A instruments
fall under the Air-Sea Interface and Surface Water regimes and the same
Dissolved CO$_2$ category.

### Primary Sources

| DCN | Document                                                                                                                                                                                                          |
|---|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1341-00490 | [OOI (2018). Data Product Specification for Partial Pressure of CO2 in Seawater.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00490_Data_Product_Spec_PCO2WAT_OOI-2.pdf)                       |
| 1341-00260 | [OOI (2012). Data Product Specification for Partial Pressure of CO2 in Air and Surface Seawater.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00260_Data_Product_SPEC_PCO2SSW_PCO2ATM_OOI.pdf) |
| 1341-00270 | [OOI (2012). Data Product Specification for Flux of CO2 from the Ocean into the Atmosphere.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00270_Data_Product_SPEC_CO2FLUX_OOI.pdf)              |

---

### PCO2WAT_L1 — Partial Pressure of CO<sub>2</sub> in Seawater

PCO2WAT_L1 is the partial pressure of CO$_2$ (pCO$_2$, $\mu$atm) in 
seawater at depth, computed from the [Sunburst SAMI-CO2](https://www.sunburstsensors.com/products/oceanographic-carbon-dioxide-sensor.html) (PCO2W). The SAMI-CO2 
uses a bromothymol blue (BTB) pH indicator solution that equilibrates with the
ambient seawater through a permeable silicone membrane. The equilibrated
indicator solution passes through an optical cell where absorbance is
measured at 434 nm and 620 nm — the peak absorbance wavelengths for the
protonated and deprotonated forms of the indicator. Periodic blank
measurements (every 3.5 days) correct for drift in the electro-optical
system.

#### L0 Inputs and Blank Normalization

The SAMI-CO2 outputs two record types:

- **Type 4** — regular measurement (used to compute PCO2WAT_L1)
- **Type 5** — blank measurement (used to correct Type 4 records)

Each record contains a 434 nm ratio (`Ratio434`), a 620 nm ratio (`Ratio620`), 
and a raw thermistor count (`Therm`). Raw blank counts from Type 5 records are
normalized to dimensionless ratios by dividing by 16384:

$$A_{434,blank} = \frac{\text{Ratio434}_{blank}}{16384}$$

$$A_{620,blank} = \frac{\text{Ratio620}_{blank}}{16384}$$

#### Thermistor Temperature

The raw thermistor count is converted to degrees C. The conversion depends
on the SAMI hardware generation (12-bit or 14-bit ADC). For 12-bit hardware
(full-scale count $= 4096$):

$$r_t = \ln\left(\frac{\text{Therm}}{4096 - \text{Therm}} \times 17400\right)$$

For 14-bit hardware (full-scale count $= 16384$):

$$r_t = \ln\left(\frac{\text{Therm}}{16384 - \text{Therm}} \times 17400\right)$$

In both cases, temperature is then:

$$T_C = \frac{1}{0.0010183 + 0.000241 \times r_t + 1.5 \times 10^{-7} \times {r_t}^3} - 273.15$$

#### PCO2WAT_L1 Calculation

The L1 pCO$_2$ is computed from the blank-corrected absorbances and a
temperature-corrected SAMI response. The blank-corrected absorbances at
each wavelength are:

$$A_{434} = -\log_{10}\left(\frac{\text{Ratio434}_{meas}}{A_{434,blank}}\right)$$

$$A_{620} = -\log_{10}\left(\frac{\text{Ratio620}_{meas}}{A_{620,blank}}\right)$$

The absorbance ratio is:

$$R = \frac{A_{620}}{A_{434}}$$

The SAMI response ($RCO_2$) and intermediate temperature-correction
terms are:

$$RCO_2 = -\log_{10}\left(\frac{R - e_1}{e_2 - e_3 \times R}\right)$$

$$RCO_2 = RCO_2 + 0.007 \times (T_C - \text{CalT})$$

$$T_{coeff} = 0.0075778 + 0.0012389 \times RCO_2 - 0.00048757 \times {RCO_2}^2$$

$$Tcor\_RCO_2 = RCO_2 + T_{coeff} \times (T_C - \text{CalT})$$

The final pCO$_2$ is obtained from the quadratic calibration curve:

$$qcc = \frac{-\text{CalB} + \sqrt{\text{CalB}^2 - 4 \times \text{CalA} \times (\text{CalC} - Tcor\_RCO_2)}}{2 \times \text{CalA}}$$

$$pCO_2 = 10^{qcc}$$

The equilibration constants $e_1 = 0.0043$, $e_2 = 2.136$, $e_3 = 0.2105$
are fixed values provided by Sunburst Sensors based on laboratory
determinations with the BTB indicator batch and are not recalculated per
calibration cycle.

The instrument-specific calibration coefficients CalT, CalA, CalB, and CalC
are provided by the manufacturer after each calibration or refurbishment
cycle and are made available to the [OOI Asset Management database via CSV files 
uploaded to GitHub](https://github.com/oceanobservatories/asset-management/).

**Algorithm note:** The Python implementation corrects errors present in the
vendor-supplied Matlab example code distributed with DPS 1341-00490
(Appendix A, 2018).

Output accuracy: $\pm 3 \times \mu\text{atm}$; precision $<1 \times \mu\text{atm}$
(DPS 1341-00490, Appendix B).

---

### PCO2ATM_L1 and PCO2SSW_L1 — pCO<sub>2</sub> in Air and Surface Seawater

PCO2ATM_L1 and PCO2SSW_L1 are the partial pressures of CO$_2$ (pCO$_2$, $\mu$atm) in
atmosphere and surface seawater, respectively. Both are computed from the
[Pro-Oceanus CO2-Pro Atmosphere](https://pro-oceanus.com/products/pro-series/co2-pro-atm) (PCO2A) instrument, which determines 
the CO$_2$ mole fraction (xCO$_2$, in ppm) internally by measuring the infrared 
absorbance of CO$_2$ and compensating for pressure, temperature, and humidity
using onboard firmware. The instrument alternates between air sampling 
(producing XCO2ATM_L0) and water sampling (producing XCO2SSW_L0).

The xCO$_2$ mole fractions and the internal gas stream pressure (PRESAIR_L0)
are converted to pCO$_2$ ($\mu$atm) using:

$$pCO_2 = \frac{\text{xCO}_2 \times p}{\text{STD}}$$

where $p$ is PRESAIR_L0 in mbar and $\text{STD} = 1013.25$ mbar/atm, the
standard atmospheric pressure. Because xCO$_2$ is in ppm
($10^{-6}$), the result is directly in $\mu$atm without further
scaling.

The same function `pco2_ppressure` computes both PCO2ATM_L1 and PCO2SSW_L1;
the distinction between air and surface seawater measurements is determined
by the L0 input supplied (XCO2ATM or XCO2SSW) and the record type flag in
the serial data stream.

Output accuracy: $\pm 1 \times \mu\text{atm}$; precision $\pm 0.01 \times \mu\text{atm}$
(DPS 1341-00260, Appendix B).

---

### CO2FLUX_L2 — Air-Sea Flux of CO<sub>2</sub>

CO2FLUX_L2 is the estimate of the flux of CO$_2$ across the air-sea 
interface. It is computed from PCO2SSW_L1 and PCO2ATM_L1 and select 
L1/L2 bulk meteorology inputs (from the METBK): sea surface
temperature (TEMPSRF_L1), sea surface salinity (SALSURF_L2), and wind
speed at 10 m height (WIND10M_L2).

The flux is positive when directed from the ocean to the atmosphere. The 
computation follows these steps:

**Step 1** — Convert pCO$_2$ (for both air and water) from $\mu$atm to atm.

$$pCO_2,air = pCO_2,air \div 10^6$$

$$pCO_2,water = pCO_2,water \div 10^6$$

**Step 2** — Compute the Schmidt number $S_c$ from sea surface temperature
$t$ ($^\circ$C) (Wanninkhof, 1992, Table A1):

$$S_c = 2073.1 - 125.62 \times t + 3.6276 \times t^2 - 0.043219 \times t^3$$

**Step 3** — Compute the gas transfer velocity $k$ in cm h$^{-1}$ and
convert to m s$^{-1}$ (Sweeney et al., 2007, Fig. 3 and Table 1):

$$k = \frac{0.27 \times {u_{10}}^2 \times \sqrt{660 \div S_c}}{100 \times 3600}$$

where $u_{10}$ is the 10 m wind speed in m s$^{-1}$.

**Step 4** — Compute absolute temperature $T = t + 273.15$ (K), then
compute CO$_2$ solubility $K_0$ in mol atm$^{-1}$ m$^{-3}$ using the volume
formulation (Weiss, 1974, Eqn. 12 and Table I):

$$T100 = T \div 100$$

$$\begin{align}
K_0 &= 1000 \times \exp(-58.0931 + 90.5069 \times (100 / T) + 22.2940 \times \ln(T100) \\
&\quad + s \times (0.027766 - 0.025888 \times T100 + 0.0050578 \times T100^2))
\end{align}$$

where $s$ is the sea surface salinity.

**Step 5** — Compute the flux (Wanninkhof, 1992, Eqn. A2):

$$F = k \times K_0 \times (pCO_2,water - pCO_2,air)$$

The inherent uncertainty of the flux estimate is approximately 10%
(DPS 1341-00270, Appendix B).

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

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
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_ppressure

#### History
| Date | Author | Change |
|---|---|---|
| 2014-10-27 | Christopher Wingard | Initial Python code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_co2flux

#### History
| Date | Author | Change |
|---|---|---|
| 2012-03-28 | Matthias Lankhorst | Original Matlab code |
| 2013-04-20 | Christopher Wingard | Initial Python code |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.co2_functions.pco2_abs434_ratio

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_abs620_ratio

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_blank

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2014-02-28 | Christopher Wingard | Updated to accept raw blank values from a sparse array |
| 2018-03-04 | Christopher Wingard | Updated blank normalization per vendor-corrected algorithm |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_thermistor

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-20 | Christopher Wingard | Initial code |
| 2014-02-19 | Christopher Wingard | Updated comments |
| 2023-01-12 | Mark Steiner | Added sami_bits argument to handle 14-bit hardware |
| 2023-08-15 | Samuel Dahlberg | Renamed local variables; replaced numexpr with numpy |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.co2_functions.pco2_battery

#### History
| Date | Author | Change |
|---|---|---|
| 2023-02-23 | Mark Steiner | Initial code |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

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
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

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
and Surface Seawater. Document Control Number 1341-00260.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00260_Data_Product_SPEC_PCO2SSW_PCO2ATM_OOI.pdf)

[OOI (2018). Data Product Specification for Partial Pressure of CO2 in
Seawater. Document Control Number 1341-00490.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00490_Data_Product_Spec_PCO2WAT_OOI-2.pdf)

[OOI (2012). Data Product Specification for Flux of CO2 from the Ocean into
the Atmosphere. Document Control Number 1341-00270.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00270_Data_Product_SPEC_CO2FLUX_OOI.pdf)
