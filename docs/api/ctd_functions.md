# CTD Functions

## Background

The Ocean Observatories Initiative deploys the Sea-Bird Electronics
conductivity-temperature-depth (CTD) family of instruments across its moored,
profiling, and mobile platforms. The table below lists the OOI instrument
classes covered by this module.

| Class | Hardware | Platform type | Designator meaning |
|---|---|---|---------------|
| CTDBP | SBE 16Plus V2 | Moored (fixed depth) | CTD, Bottom Pumped |
| CTDMO | SBE 37IM | Moored (fixed depth) | CTD, Modem (Inductive) |
| CTDPF | SBE 16Plus V2 (A/B) or SBE 52MP (C/K/L) | Profiling | CTD, Profiler |
| CTDGV | SBE GPCTD (Seabird Payload CTD) | Gliders | CTD, Glider Vehicle |
| PHSEN | SBE Deep SeapHOx V2 | Moored (fixed depth) | pH Sensor |

`ctd_functions.py` converts raw Sea-Bird Electronics CTD data into calibrated 
L1 engineering products (TEMPWAT, PRESWAT, CONDWAT) and computes the L2 
derived products Practical Salinity (PRACSAL) and in-situ Density (DENSITY) 
using the [TEOS-10 GSW library](https://teos-10.github.io/GSW-Python/). All 
calibration coefficients are from factory calibration values supplied with 
individual instruments. Within the OOI data system, CTD functions fall under
the Water Column science regime and the Conductivity, Temperature, Pressure,
Salinity, and Density categories.

The L1 engineering products (TEMPWAT, PRESWAT, CONDWAT) may be reported directly
by the instrument (computed onboard the sensor using vendor firmware) or they 
may be reported as L0 values in raw units (counts, frequency, or linearly scaled 
integers). Depending on the instrument, telemetered data may be reported directly 
in L1 units while the recovered instrument data is reported in L0 units 
requiring the conversions below.  

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00010 | [OOI (2013). Data Product Specification for Water Temperature.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf) |
| 1341-00020 | [OOI (2013). Data Product Specification for Pressure (Depth).](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf) |
| 1341-00030 | [OOI (2013). Data Product Specification for Conductivity.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf) |
| 1341-00040 | [OOI (2013). Data Product Specification for Salinity.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf) |
| 1341-00050 | [OOI (2012). Data Product Specification for Density.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00050_Data_Product_SPEC_DENSITY_OOI.pdf) |

---

### CONDWAT_L1 — Seawater Conductivity

The raw seawater conductivity (CONDWAT_L0, $c_0$), is reported in either 
counts or scaled L1 units. The converted seawater conductivity, 
CONDWAT_L1, is reported in S m$^{-1}$. As with temperature, conversions from 
L0 to L1 depend on the instrument class and/or data delivery method.

**SBE 16Plus** (CTDBP, CTDPF-A/B): The raw count is converted to a frequency
$f$ in kHz, then evaluated with a polynomial corrected for temperature and
pressure:

$$f = \frac{c_0 \div  256}{1000}$$

$$C_{L1} = \frac{g + h f^2 + i f^3 + j f^4}{1 + \text{CTcor}\times T + \text{CPcor}\times P}$$

where $T$ is TEMPWAT_L1 ($^\circ\text{C}$), $P$ is PRESWAT_L1 (dbar), and
$g$, $h$, $i$, $j$, CTcor, CPcor are factory calibration coefficients.

**SBE 37IM instrument-recovered** (CTDMO): Same polynomial as the SBE 16Plus
path, but includes an additional wbotc correction to the frequency:

$$f = \frac{c_0 \div  256}{1000} \times \sqrt{1 + \text{wbotc}\times T}$$

**SBE 37IM telemetered and recovered_host** (CTDMO):

$$C_{L1} = c_0 \div 100000 - 0.5$$

**SBE 37 Deep SeapHOx V2 streamed and recovered_instrument** (PHSEN):
#TODO

**SBE 52MP** (CTDPF-C/K/L): Linear scaling with a unit conversion from
mS cm$^{-1}$ to S m$^{-1}$:

$$C_{L1} = (c_0 \div 10000 - 0.5) \times 0.1$$

**CTDGV (glider)**: The SBE GPCTD computes conductivity onboard the vehicle
using vendor software and transmits the result already in S m$^{-1}$;
`ion-functions` is not invoked for those deployments.

Output accuracy: SBE 16Plus V2 $\pm 0.0005$ S m$^{-1}$; SBE 37IM
$\pm 0.0003$ S m$^{-1}$.

---

### TEMPWAT_L1 — Seawater Temperature

The raw seawater temperature (TEMPWAT_L0, $t_0$), is reported in either 
counts or in L1 units ($^\circ$C), but scaled to an integer to compress the 
data for internal storage or transmission. The converted seawater temperature, 
TEMPWAT_L1, is reported in $^\circ$C. Conversions from L0 to L1 depend on the 
instrument class and/or data delivery method.

**SBE 16Plus** (CTDBP, CTDPF-A/B): Raw counts $t_0$ are converted to seawater
temperature ($^\circ$C) via an intermediate resistance quantity $R$:

$$MV = \frac{t_0 - 524288}{1.6 \times 10^7}$$

$$R = \frac{MV \times 2.9 \times 10^9 + 1.024 \times 10^8}{2.048 \times 10^4 - MV \times 2.0 \times 10^5}$$

$$T_{L1} = \frac{1}{a_0 + a_1 \ln R + a_2 \ln^2 R + a_3 \ln^3 R} - 273.15$$

where $a_0$, $a_1$, $a_2$, $a_3$ are factory calibration coefficients.

**SBE 37IM telemetered and recovered_host** (CTDMO): The instrument outputs
engineering units scaled to an integer. The L1 conversion is:

$$T_{L1} = t_0 / 10000 - 10$$

**SBE 37IM instrument-recovered** (CTDMO): Raw decimal counts $t_0$ are
retained in the instrument-stored file and require the calibration-coefficient
equation:

$$T_{L1} = \frac{1}{a_0 + a_1 \ln t_0 + a_2 \ln^2 t_0 + a_3 \ln^3 t_0} - 273.15$$

**SBE 37 Deep SeapHOx V2 streamed and recovered_instrument** (PHSEN):
#TODO

**SBE 52MP** (CTDPF-C/K/L): The instrument outputs engineering units scaled 
to an integer. The L1 conversion is:

$$T_{L1} = t_0 / 10000 - 5$$

**CTDGV (glider)**: The SBE GPCTD computes temperature onboard the vehicle
using vendor software and transmits the result already in $^\circ$C;
`ion-functions` is not invoked for those deployments.

Output accuracy: SBE 16Plus V2 $\pm 0.005\ ^\circ$C; SBE 37IM
$\pm 0.002\ ^\circ$C.

---

### PRESWAT_L1 — Seawater Pressure

The raw seawater pressure (PRESWAT_L0, $p_0$), is reported in either 
counts, frequency or in scaled L1 units (dbar). The converted seawater 
pressure, PRESWAT_L1, is reported in dbar relative to one standard atmosphere
(10.1325 dbar). As with temperature, conversions from L0 to L1 depend on the 
instrument class and/or data delivery method.

**SBE 16Plus — strain-gauge pressure sensor** (CTDBP except N/O, CTDPF-A/B):

$$t_v = $pt_0$ \div  13107$$

$$t = \text{PTEMPA0} + \text{PTEMPA1} \times t_v + \text{PTEMPA2} \times t_v^2$$

$$x = p_0 - \text{PTCA0} - \text{PTCA1 }\times t - \text{PTCA2} \times t^2$$

$$n = \frac{x \times \text{PTCB0}}{\text{PTCB0} + \text{PTCB1} \times t + \text{PTCB2} \times t^2}$$

$$p_{\text{psi}} = \text{PA0} + \text{PA1} \times n + \text{PA2} \times n^2$$

$$P_{L1}\ = p_{\text{psi}} \times 0.689475729 - 10.1325 + \delta$$

where $pt_0$ is the pressure sensor thermistor in counts, $\delta$ is an 
optional Druck sensor offset correction (default 0 dbar) and all 
calibration coefficients are from factory calibration sheets.

**SBE 16Plus — digiquartz pressure sensor** (CTDBP-N and CTDBP-O only):

$$p_f = p_0 \div  256$$

$$t_v = t_0 \div  13107$$

$$U = 23.7(t_v + 9.7917) - 273.15$$

$$C = C_1 + C_2 U + C_3 U^2$$

$$D = D_1 + D_2 U$$

$$T_0 = T_1 + T_2 U + T_3 U^2 + T_4 U^3 + T_5 U^4$$

$$\tau = (1\div p_f) \times 10^6$$

$$p_{\text{psi}} = C\left(1 - \frac{T_0^2}{\tau^2}\right)\left(1 - D\left(1 - \frac{T_0^2}{\tau^2}\right)\right)$$

$$P_{L1} = p_{\text{psi}} \times 0.689475729 - 10.1325$$

All calibration coefficients are from factory calibration sheets.

**SBE 37IM telemetered and recovered_host** (CTDMO): Raw pressure is a scaled 
integer relative to a factory-set full-scale pressure range $P_{rng}$ (in psi):

$$P_{rng,dbar} = (P_{rng,psi} - 14.7) \times 0.6894757$$

$$P_{L1} = \frac{p_0 \times P_{rng,dbar}}{0.85 \times 65536} - 0.05 \times P_{rng,dbar}$$

**SBE 37IM instrument-recovered** (CTDMO): Uses the same strain-gauge
polynomial as the SBE 16Plus strain-gauge path, but with the raw thermistor
count $pt_0$ used directly in the polynomial rather than being first
converted to voltage.

**SBE 37 Deep SeapHOx V2 streamed and recovered_instrument** (PHSEN):
#TODO

**SBE 52MP** (CTDPF-C/K/L):

$$P_{L1} = p_0 \div 100 - 10$$

**CTDGV (glider)**: The SBE GPCTD reports pressure in bar; `ion-functions`
converts to dbar:

$$P_{L1}\ [\text{dbar}] = p_{\text{bar}} \times 10$$

Output accuracy: SBE 16Plus V2 and SBE 37IM 0.1 % of full-scale range.

---

### PRACSAL_L2 — Seawater Practical Salinity

Seawater practical salinity is computed from L1 conductivity, temperature, and
pressure using the TEOS-10 GSW library function `gsw.SP_from_C`, which
implements the Practical Salinity Scale 1978 (PSS-78) algorithm.
Conductivity must be converted from S m$^{-1}$ to mS cm$^{-1}$ (multiply by 10) 
before calling the GSW function.

Seawater practical salinity is dimensionless and reported without units on 
the PSS-78 scale.

$$PS_{L2} =\text{gsw.SP_from_C}(C_{L1} * 10, T_{L1}, P_{L1})$$

---

### DENSITY_L2 — In-situ Density

In-situ seawater density is computed via a three-step chain using the TEOS-10
GSW library:

**Step 1 — Absolute Salinity:**

$$S_A = \text{gsw.SA_from_SP}(PS_{L2}, P_{L1}, \text{lon}, \text{lat})$$

Absolute Salinity $S_A$ (g kg$^{-1}$) is derived from Practical Salinity
using a lookup table of the Absolute Salinity Anomaly (SAAR) as a function of
location and pressure. For moored instruments, latitude and longitude are the
mooring position metadata; for gliders, they are the vehicle position at each
sample.

**Step 2 — Conservative Temperature:**

$$\Theta = \text{gsw.CT_from_t}(S_A, T_{L1}, P_{L1})$$

**Step 3 — In-situ density:**

$$\rho = \text{gsw.rho}(S_A, \Theta, p)$$

Density is computed using the computationally-efficient 48-term expression
described in McDougall et al. (2011).

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core functions

::: ion_functions.data.ctd_functions.ctd_sbe16plus_tempwat

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-12 | Luke Campbell | Initial code |
| 2013-04-12 | Christopher Wingard | Minor edits |
| 2013-05-10 | Christopher Wingard | Minor edits to comments |
| 2014-01-31 | Russell Desiderio | Standardized comment format |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe37im_tempwat_instrument_recovered

#### History
| Date | Author | Change |
|---|---|---|
| 2016-06-16 | Russell Desiderio | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe37im_tempwat

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-05 | Russell Desiderio | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe52mp_tempwat

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-17 | Russell Desiderio | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe16plus_preswat

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-12 | Christopher Wingard | Initial code |
| 2013-05-10 | Christopher Wingard | Minor edits to comments |
| 2014-01-31 | Russell Desiderio | Standardized comment format |
| 2017-03-31 | Dan Mergens | Added Druck offset correction |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe16digi_preswat

#### Additional Notes
Per the SBE 16Plus V2 User Manual, the raw pressure input `p0` is in units
of Hz (counts divided by 256), not raw A/D counts; the code implements this
correctly.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-05-10 | Christopher Wingard | Initial code |
| 2014-01-31 | Russell Desiderio | Standardized comment format; corrected pressure period calculation to use Hz input |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe37im_preswat_instrument_recovered

#### History
| Date | Author | Change |
|---|---|---|
| 2016-06-16 | Russell Desiderio | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe37im_preswat

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-05 | Russell Desiderio | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_glider_preswat

#### History
| Date | Author | Change |
|---|---|---|
| 2015-10-28 | Russell Desiderio | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe52mp_preswat

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-17 | Russell Desiderio | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe16plus_condwat

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-12 | Christopher Wingard | Initial code |
| 2013-05-10 | Christopher Wingard | Minor edits to comments |
| 2014-01-31 | Russell Desiderio | Standardized comment format |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe37im_condwat_instrument_recovered

#### History
| Date | Author | Change |
|---|---|---|
| 2016-06-16 | Russell Desiderio | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe37im_condwat

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-05 | Russell Desiderio | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_sbe52mp_condwat

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-17 | Russell Desiderio | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_pracsal

#### History
| Date | Author | Change |
|---|---|---|
| 2013-03-13 | Christopher Wingard | Initial code |
| 2013-05-10 | Christopher Wingard | Minor edits to comments |
| 2014-01-31 | Russell Desiderio | Standardized comment format |
| 2023-08-15 | Samuel Dahlberg | Replaced incompatible pygsw with GSW library |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ctd_functions.ctd_density

#### History
| Date | Author | Change |
|---|---|---|
| 2013-03-11 | Christopher Mueller | Initial code |
| 2013-03-13 | Christopher Wingard | Added commenting; moved to ctd_functions |
| 2013-05-10 | Christopher Wingard | Minor edits to comments |
| 2014-01-31 | Russell Desiderio | Standardized comment format |
| 2023-08-15 | Samuel Dahlberg | Replaced incompatible pygsw with GSW library |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## SeapHOx CTD Functions

The Sea-Bird Scientific Deep SeapHOx V2 (PHSEN-G and PHSEN-H) integrates
an SBE 37 MicroCAT CTD alongside an SBE 63 dissolved oxygen sensor and a
Deep SeaFET pH sensor. The three functions below convert raw SBE 37 A/D
counts to calibrated temperature, pressure, and conductivity for use by
the SeapHOx dissolved oxygen and pH algorithms. They live in
`phsen_h_functions.py`, but are documented here because they implement the 
SBE 37 calibration equations already covered by this page. 
See [Deep SeapHOx V2](../seaphox.md) for instrument architecture context.

::: ion_functions.data.phsen_h_functions.temperature_raw_conversion

#### History
| Date | Author | Change |
|---|---|---|
| 2025-05-15 | Christopher Wingard | Initial NumPy docstring; added to CTD family documentation |

---

::: ion_functions.data.phsen_h_functions.pressure_raw_conversion

#### History
| Date | Author | Change |
|---|---|---|
| 2025-05-15 | Christopher Wingard | Initial NumPy docstring; added to CTD family documentation |

---

::: ion_functions.data.phsen_h_functions.conductivity_raw_conversion

#### History
| Date | Author | Change |
|---|---|---|
| 2025-05-15 | Christopher Wingard | Initial NumPy docstring; added to CTD family documentation |

---

## References

Feistel, R. (2008). A Gibbs function for seawater thermodynamics for -6 to
80 degrees C and salinity up to 120 g/kg. Deep Sea Research I, 55,
1639-1671.

Fofonoff, N. P. and Millard, R. C. (1983). Algorithms for computation of
fundamental properties of seawater. UNESCO Technical Papers in Marine
Science, 44, 1-53.

Hill, K. D., Dauphinee, T. M., and Woods, D. J. (1986). The extension of the
Practical Salinity Scale 1978 to low salinities. IEEE Journal of Oceanic
Engineering, OE-11(1), 109-112.

IOC, SCOR and IAPSO (2010). The international thermodynamic equation of
seawater - 2010: Calculation and use of thermodynamic properties.
Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
UNESCO, 196 pp.

McDougall, T. J., Barker, P. M., Feistel, R., and Jackett, D. R. (2011). A
computationally efficient 48-term expression for the density of seawater in
terms of Conservative Temperature, and related properties of seawater.
Journal of Atmospheric and Oceanic Technology, 28, 1464-1477.

McDougall, T. J., Jackett, D. R., and Millero, F. J. (2010). An algorithm
for estimating Absolute Salinity in the global ocean. Ocean Science
Discussions, 6, 215-242.

Millero, F. J., et al. (2008). The composition of Standard Seawater and the
definition of the Reference-Composition Salinity Scale. Deep Sea Research I,
55, 50-72.

Pawlowicz, R. (2010). What every oceanographer needs to know about TEOS-10
(The TEOS-10 Primer). Thermodynamic Equation Of Seawater - 2010 (TEOS-10)
website: http://www.teos-10.org/.

Sea-Bird Electronics (2008). Application Note 10. Compressibility
Compensation of Sea-Bird Conductivity Sensors. Revision March 2008. Sea-Bird
Electronics, Inc.

Sea-Bird Electronics (2009). SBE 16plus V2 SEACAT User's Manual. Manual
Version 005. Sea-Bird Electronics, Inc.

Sea-Bird Electronics (2010). Application Note 42. ITS-90 Temperature Scale.
Revision February 2010. Sea-Bird Electronics, Inc.

Sea-Bird Electronics (2011). SBE 37-IM MicroCAT User's Manual. Manual
Version 027. Sea-Bird Electronics, Inc.

[OOI (2013). Data Product Specification for Water Temperature. Document
Control Number 1341-00010.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)

[OOI (2013). Data Product Specification for Pressure (Depth). Document
Control Number 1341-00020.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)

[OOI (2013). Data Product Specification for Conductivity. Document Control
Number 1341-00030.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)

[OOI (2013). Data Product Specification for Salinity. Document Control
Number 1341-00040.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf)

[OOI (2012). Data Product Specification for Density. Document Control Number
1341-00050.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00050_Data_Product_SPEC_DENSITY_OOI.pdf)

Sea-Bird Scientific. *Deep SeapHOx V2 Ocean CT(D)-pH-DO Sensor*. Data
sheet DS.53.May25. Bellevue, WA: Sea-Bird Scientific.
