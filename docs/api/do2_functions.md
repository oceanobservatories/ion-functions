# DO2 Functions

## Background

The Ocean Observatories Initiative deploys two families of dissolved oxygen
instruments across its moored, profiling, and mobile platforms. A third
dissolved oxygen sensor is integrated into the Deep SeapHOx V2 pH instrument
(PHSEN-G and PHSEN-H). The table below lists the OOI instrument classes
covered by this module.

| Class | Hardware | Platform | Designator Meaning |
|---|---|---|---|
| DOSTA | [Aanderaa Optode 4831](#dosta-functions-aanderaa-optode-4831) | Moored (fixed depth) and gliders | DO, Stable |
| DOFST-A | [Sea-Bird Scientific SBE 43](#dofst-functions-sbe-43-sbe-43f) | RCA cabled shallow profilers and platforms (CTDPF-A/B) | DO, Fast, SBE 43 |
| DOFST-K | [Sea-Bird Scientific SBE 43F](#dofst-functions-sbe-43-sbe-43f) | Wire-following profilers (CTDPF-K/L) | DO, Fast, SBE 43F |
| PHSEN-G/H | [Sea-Bird Scientific Deep SeapHOx V2 (SBE 63)](#deep-seaphox-v2-do-functions-sbe-63) ([instrument notes](../seaphox.md)) | Moored (fixed depth) | pH Sensor (integrated DO) |

`do2_functions.py` processes data from the DOSTA and DOFST instrument classes.
The dissolved oxygen functions for the PHSEN-G/H integrated SBE 63 sensor live
in `phsen_h_functions.py`; they are documented here because the SBE 63 is a
dissolved oxygen sensor. See [Deep SeapHOx V2](../seaphox.md) for instrument
architecture context.

Within the OOI data system, dissolved oxygen data products fall under the Water Column [sampling
regime](https://oceanobservatories.org/data-products/) (see the [complete list of OOI data
products](https://oceanobservatories.org/data-products/data-products-complete-list/)).

**Units identify the processing level.** All OOI dissolved oxygen products in $\mu$mol L$^{-1}$ or mL L$^{-1}$
(DOCONCS_L1) are Level 1 products that have not been corrected for salinity or pressure. All products in
$\mu$mol kg$^{-1}$ (DOXYGEN_L2, DOCONCF_L2) are Level 2 products that have been corrected for salinity and pressure;
converting from per-volume to per-mass units requires salinity, which is what makes the result a Level 2 product.

Calibration coefficients are from the factory calibration sheets supplied with individual instruments, except where
noted in the Calibration Coefficients tables below.

This page documents the algorithms used to compute the dissolved oxygen data products. For general guidance on working
with OOI biogeochemical sensor data, see the [OOI Biogeochemical Sensor Data Best Practices and User
Guide](https://doi.org/10.25607/OBP-1865.2); for locating calibration coefficients and raw inputs, see [Calibration
Coefficients](../calibrations.md) and [Using These Functions with OOI Data](../using-ooi-data.md).

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00520 | [OOI (2014). Data Product Specification for Oxygen Concentration from "Stable" Instruments.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf) |
| 1341-00521 | [OOI (2013). Data Product Specification for Fast Dissolved Oxygen.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00521_Data_Product_Spec_DOCONCF_OOI.pdf) |

---

### DOCONCS_L1 and DOXYGEN_L2 -- Stable Dissolved Oxygen (DOSTA)

The Aanderaa Optode 4831 (DOSTA) measures dissolved oxygen using a
luminescence quenching principle. The optode's fluorescent foil emits light
at a frequency that depends on the oxygen partial pressure at the foil surface.
The instrument measures the phase shift between the excitation and emission
signals, producing a calibrated phase value (CalPhase, $P_t$) as its L0
output. The optode also reports the foil temperature ($T_{opt}$) from an
internal thermistor positioned directly at the foil.

DOCONCS_L1 is the dissolved oxygen concentration ($\mu$mol L$^{-1}$)
uncorrected for salinity and pressure. DOXYGEN_L2 is the final
salinity- and pressure-corrected concentration ($\mu$mol kg$^{-1}$). The
Level 2 product was named DOCONCS_L2 in early OOI documentation and renamed DOXYGEN_L2
in DPS 1341-00520 version 1-02 (2014-04-11); `do2_salinity_correction`
produces DOXYGEN_L2.

#### Calibration Coefficients

| Function Argument | Calibration Coefficient | Source |
|---|---|---|
| `csv` | `CC_csv` | Factory calibration sheet (SVU foil calibration, $\text{csv}_1$-$\text{csv}_7$) |
| `conc_coef` | `CC_conc_coef` | Factory calibration sheet (secondary calibration; Form 710 "PhaseCoef") |
| `sref` | `CC_ref_salinity` | Instrument configuration (optode firmware reference salinity; glider DOSTA only) |

OOI calibration coefficient names are explained in [Using These Functions with OOI Data](../using-ooi-data.md).

#### DOSTA L0 Data Streams

The DOSTA reaches the OOI data system through three distinct hardware
configurations, each producing a different L0 format:

**Analog output to CTD voltage channel** -- The optode outputs $P_t$ and
$T_{opt}$ as 0-5 V analog signals routed into CTD auxiliary voltage channels.
The voltage signals are parsed as DOCONCS-VLT_L0 (phase) and a separate
temperature voltage, then converted to degrees and $^\circ$C:

$$P_t\ [\text{deg}] = 10 + 12 \times P_t\ [\text{V}]$$

$$T_{opt}\ [^\circ\text{C}] = -5 + 8 \times T_{opt}\ [\text{V}]$$

These fixed conversion coefficients are universal for all Aanderaa optodes
(Shawn Sneddon, Xylem-Aanderaa).

**Digital output via RS-232 to CTD** -- The optode reports an oxygen
concentration directly in internal units transmitted over RS-232 to an SBE
16Plus V2 CTD, parsed as DOCONCS-CNT_L0:

$$\text{DOCONCS_L1}\ [\mu\text{mol/L}] = \frac{\text{DOCONCS-CNT_L0}}{10000} - 10$$

**Autonomous digital output** -- The optode operates independently and reports
$P_t$ (DOCONCS-DEG_L0) and $T_{opt}$ directly in degrees and $^\circ$C,
requiring no voltage conversion before the SVU equation.

#### DOSTA DOCONCS_L1 -- Stern-Volmer-Uchida Equation

Where DOCONCS-DEG_L0 is available (analog or autonomous configurations),
DOCONCS_L1 is computed using the Stern-Volmer-Uchida (SVU) equation
(Uchida et al. 2008). The SVU equation relates the calibrated phase to
oxygen concentration with a temperature-dependent sensitivity:

$$K_{SV} = \text{csv}_1 + \text{csv}_2 \times T_{opt} + \text{csv}_3 \times {T_{opt}}^2$$

$$P_0 = \text{csv}_4 + \text{csv}_5 \times T_{opt}$$

$$P_c = \text{csv}_6 + \text{csv}_7 \times P_t$$

$$O_2 = \frac{(P_0 / P_c) - 1}{K_{SV}}$$

where $T_{opt}$ is the optode foil temperature ($^\circ$C), $P_t$ is the
calibrated phase (degrees), $P_0$ is the zero-oxygen phase shift, $P_c$ is
the phase at the current oxygen level, and $\text{csv}_1$ through
$\text{csv}_7$ (`csv`) are from the factory calibration sheet.

A secondary calibration (`conc_coef` = [offset, slope], from the factory calibration sheet) is applied after the SVU
equation:

$$O_2 = \text{conc_coef}[0] + \text{conc_coef}[1] \times O_2$$

Aanderaa applies this two-point correction after refurbishment; for new
optodes or new SVU foil determinations conc_coef defaults to [0, 1].

The result is DOCONCS_L1 in $\mu$mol L$^{-1}$.

The optode foil thermistor temperature ($T_{opt}$) must be used as the
temperature input to the SVU equation, not CTD temperature. The thermistor is
situated directly at the foil; the SVU calibration coefficients are derived to
compensate for changes in oxygen permeability through the foil as a function
of its temperature. On gliders, differences between CTD and optode temperature
of 1$^\circ$C correspond to approximately 5% differences in calculated oxygen
concentration.

#### DOSTA DOXYGEN_L2 -- Salinity and Pressure Correction

DOCONCS_L1 ($\mu$mol L$^{-1}$) is corrected for salinity and pressure using
CTD data from the co-located instrument to produce DOXYGEN_L2 ($\mu$mol
kg$^{-1}$). The correction has three sequential steps.

**Step 1 -- Volume to mass conversion.** Potential density $\rho$ is computed
at reference pressure $p_{ref} = 0$ dbar from absolute salinity and
conservative temperature (TEOS-10 GSW library):

$$S_A = \text{gsw.SA_from_SP}(S_P, p, \text{lon}, \text{lat})$$

$$\Theta = \text{gsw.CT_from_t}(S_A, T, p)$$

$$\rho = \text{gsw.rho}(S_A, \Theta, p_{ref})$$

$$O_{2,\text{mass}} = \frac{1000 \times O_2 [\mu\text{mol/L}]}{\rho\ [\text{kg/m}^3]}$$

**Step 2 -- Pressure correction** (Uchida et al. 2008):

$$O_{2,\text{pres}} = \left(1 + \frac{0.032 \times p}{1000}\right) \times O_{2,\text{mass}}$$

**Step 3 -- Salinity correction** (Garcia and Gordon 1992, Table 1,
combined fit):

$$t_s = \ln\!\left(\frac{298.15 - T}{273.15 + T}\right)$$

$$B(t_s) = B_0 + B_1 t_s + B_2 {t_s}^2 + B_3 {t_s}^3$$

$$\text{DOXYGEN_L2} = \exp\!\left[(S_P - S_{ref}) \times B(t_s) + C_0 \times ({S_P}^2 - {S_{ref}}^2)\right] \times O_{2,\text{pres}}$$

where $B_0 = -6.24097 \times 10^{-3}$, $B_1 = -6.93498 \times 10^{-3}$,
$B_2 = -6.90358 \times 10^{-3}$, $B_3 = -4.29155 \times 10^{-3}$,
$C_0 = -3.11680 \times 10^{-7}$, and $S_{ref}$ is the reference salinity configured in the Aanderaa optode firmware
(typically 0 or 35). For the glider DOSTA it is passed as `sref`; the other DOSTA streams use the function default of 0.

---

### DOCONCF_L2 -- Fast Dissolved Oxygen (DOFST)

The Sea-Bird Scientific SBE 43 (DOFST-A) and SBE 43F (DOFST-K) are Clark-type
polarographic dissolved oxygen sensors integrated into profiling CTD systems.
The SBE 43 is mounted on the SBE 16Plus V2 CTD (CTDPF-A/B) and reports
a voltage signal; the SBE 43F is mounted on the SBE 52-MP profiling CTD
(CTDPF-K/L) and reports a frequency signal.

DOCONCF_L2 is the dissolved oxygen concentration ($\mu$mol kg$^{-1}$)
corrected for temperature, salinity, and pressure. There is no intermediate
L1 product; DOCONCF_L0 (counts representing voltage or frequency) is converted
directly to DOCONCF_L2.

#### Calibration Coefficients

| Function Argument | Calibration Coefficient | Source |
|---|---|---|
| `Soc` | `CC_oxygen_signal_slope` | Factory calibration sheet |
| `Foffset` | `CC_frequency_offset` | Factory calibration sheet (SBE 43F) |
| `Voffset` | `CC_frequency_offset` | Factory calibration sheet (SBE 43; a voltage offset despite the name) |
| `A` | `CC_residual_temperature_correction_factor_a` | Factory calibration sheet |
| `B` | `CC_residual_temperature_correction_factor_b` | Factory calibration sheet |
| `C` | `CC_residual_temperature_correction_factor_c` | Factory calibration sheet |
| `E` | `CC_residual_temperature_correction_factor_e` | Factory calibration sheet |

OOI calibration coefficient names are explained in [Using These Functions with OOI Data](../using-ooi-data.md).

#### DOFST L0 Input Conversion

**SBE 43 (DOFST-A):** The raw hex voltage counts are converted to volts:

$$V = \frac{\text{DOCONCF-CNT_L0}}{13107}$$

**SBE 43F (DOFST-K):** The raw hex frequency value is used directly as
$F$ [Hz].

#### DOFST DOCONCF_L2 Algorithm

Sea-Bird uses an algorithm based on Owens and Millard (1985). With the
derivative term removed by setting tau = 0 (recommended by Sea-Bird and the
OOI DPS to prevent amplification of residual noise in deep water):

**Step 1 -- Oxygen solubility** (Garcia and Gordon 1992, Table 1, 1st column):

$$T_s = \ln\!\left(\frac{298.15 - T}{273.15 + T}\right)$$

$$\begin{align}
\text{Oxsol}(T, S) &= \exp\!\big[A_0 + A_1 T_s + A_2 {T_s}^2 + A_3 {T_s}^3
+ A_4 {T_s}^4 + A_5 {T_s}^5 \\
&\quad + S(B_0 + B_1 T_s + B_2 {T_s}^2 + B_3 {T_s}^3) + C_0 S^2\big]
\end{align}$$

where $A_0 = 2.00907$, $A_1 = 3.22014$, $A_2 = 4.0501$, $A_3 = 4.94457$,
$A_4 = -0.256847$, $A_5 = 3.88767$, $B_0 = -0.00624523$,
$B_1 = -0.00737614$, $B_2 = -0.010341$, $B_3 = -0.00817083$,
$C_0 = -4.88682 \times 10^{-7}$.

**Step 2 -- Intermediate oxygen concentration** [mL L$^{-1}$]:

$$\text{DO_int} = \text{Soc} \times (x + \text{offset}) \times \text{Oxsol}(T, S)
\times (1 + AT + BT^2 + CT^3) \times e^{EP/K}$$

where $x$ is $V$ (SBE 43) or $F$ (SBE 43F), $K = T + 273.15$ is absolute
temperature, and Soc, offset, A, B, C, E (`Soc`, `Foffset` or `Voffset`, `A`, `B`, `C`, `E`) are from the factory
calibration sheet. The SBE 43 voltage offset (`Voffset`, V) is stored under the name `CC_frequency_offset`, the same
name as the SBE 43F frequency offset.

**Step 3 -- Convert to $\mu$mol kg$^{-1}$:**

Potential density $\rho_\theta$ is computed using `gsw.pot_rho_t_exact` with
a reference pressure of 0 dbar:

$$\text{DOCONCF_L2} = \frac{\text{DO_int} \times 44660}{\rho_\theta}$$

---

### DOXYGEN_L2 -- SeapHOx Dissolved Oxygen

The Sea-Bird Scientific Deep SeapHOx V2 (PHSEN-G and PHSEN-H; see the [Deep SeapHOx V2 instrument notes](../seaphox.md)) integrates an
SBE 63 optical dissolved oxygen sensor alongside an SBE 37 MicroCAT CTD and
a Deep SeaFET pH sensor. The SBE 63 uses a Stern-Volmer fluorescence quenching
approach with a calibration model distinct from the DOSTA Aanderaa optode.
Processing depends on the data delivery method.

**Streamed and recovered_instrument**: Raw phase output ($\mu$s) is converted
directly to DOXYGEN_L2 ($\mu$mol kg$^{-1}$) via `dissolved_oxygen` in
`phsen_h_functions.py`. Salinity and pressure corrections and the conversion
from mL L$^{-1}$ to $\mu$mol kg$^{-1}$ using potential density are all
applied within a single function call.

#### Calibration Coefficients

| Function Argument | Calibration Coefficient | Source |
|---|---|---|
| `a0` | `CC_set_a0` | Factory calibration sheet |
| `a1` | `CC_set_a1` | Factory calibration sheet |
| `a2` | `CC_set_a2` | Factory calibration sheet |
| `b0` | `CC_set_b0` | Factory calibration sheet |
| `b1` | `CC_set_b1` | Factory calibration sheet |
| `c0` | `CC_set_c0` | Factory calibration sheet |
| `c1` | `CC_set_c1` | Factory calibration sheet |
| `c2` | `CC_set_c2` | Factory calibration sheet |
| `coeff_e` | `CC_set_e` | Factory calibration sheet |
| `therm_ta0` | `CC_set_ta0` | Factory calibration sheet (thermistor) |
| `therm_ta1` | `CC_set_ta1` | Factory calibration sheet (thermistor) |
| `therm_ta2` | `CC_set_ta2` | Factory calibration sheet (thermistor) |
| `therm_ta3` | `CC_set_ta3` | Factory calibration sheet (thermistor) |

OOI calibration coefficient names are explained in [Using These Functions with OOI Data](../using-ooi-data.md).

**Telemetered**: The instrument reports dissolved oxygen in mL L$^{-1}$ (DOCONCS_L1). No ion-functions function
currently converts this value to DOXYGEN_L2. DOXYGEN_L2 is produced only from the raw phase in the streamed data
via (`dissolved_oxygen`).

---

## DOSTA Functions (Aanderaa Optode 4831)

### Core Functions

::: ion_functions.data.do2_functions.do2_SVU

#### Additional Notes
The model 4831 thermistor time constant is less than 2 seconds, while the foil (and therefore calphase) response
time is 8 or 24 seconds depending on the optode, so there is little or no advantage in using a faster temperature
sensor such as a CTD; it is more important that the temperature reflect the foil temperature.

For standard refurbishment Aanderaa recommends the secondary calibration rather than a new multi-point SVU foil
calibration. The conc_coef values from the secondary calibration are stored in the optode firmware and listed on the
Aanderaa Form No. 710 calibration certificate, where they are mislabelled "PhaseCoef". (As of 2015 Aanderaa was
moving the secondary calibration to a one-point, 100% oxygen calibration, still expressed as offset and slope.)

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-26 | Stuart Pearce | Initial code |
| 2015-04-10 | Russell Desiderio | Revised for CI calibration coefficient implementation |
| 2015-08-04 | Russell Desiderio | Added documentation |
| 2015-08-10 | Russell Desiderio | Added conc_coef calibration array to argument list |
| 2015-10-28 | Russell Desiderio | Added atleast_2d handling for 1D conc_coef arrays |
| 2023-08-15 | Samuel Dahlberg | Renamed local variables to follow naming convention |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.do2_salinity_correction

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-26 | Stuart Pearce | Initial code |
| 2015-08-04 | Russell Desiderio | Added Garcia-Gordon reference |
| 2021-12-16 | Stuart Pearce | Added salinity reference parameter |
| 2023-08-15 | Samuel Dahlberg | Removed numexpr; replaced pygsw with gsw |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

### Helper Functions

::: ion_functions.data.do2_functions.dosta_phase_volt_to_degree

#### History
| Date | Author | Change |
|---|---|---|
| 2015-08-04 | Russell Desiderio | Initial code |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.dosta_Topt_volt_to_degC

#### History
| Date | Author | Change |
|---|---|---|
| 2015-08-04 | Russell Desiderio | Initial code |
| 2023-08-15 | Samuel Dahlberg | Renamed local variables to follow naming convention |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.o2_counts_to_uM

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-26 | Stuart Pearce | Initial code |
| 2015-04-10 | Russell Desiderio | Added documentation and fill value handling |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## DOFST Functions (SBE 43 / SBE 43F)

### Core Functions

::: ion_functions.data.do2_functions.dofst_calc

#### History
| Date | Author | Change |
|---|---|---|
| 2013-08-20 | Stuart Pearce | Initial code |
| 2015-08-04 | Russell Desiderio | Added Garcia-Gordon reference |
| 2023-08-15 | Samuel Dahlberg | Added freq variable for CGSN compatibility |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

### Wrapper Functions

::: ion_functions.data.do2_functions.do2_dofst_volt

#### History
| Date | Author | Change |
|---|---|---|
| 2013-08-20 | Stuart Pearce | Initial code |
| 2015-08-05 | Russell Desiderio | Added fill value conversion to NaN |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.do2_dofst_frequency

#### History
| Date | Author | Change |
|---|---|---|
| 2013-08-20 | Stuart Pearce | Initial code |
| 2015-08-05 | Russell Desiderio | Added fill value conversion to NaN |
| 2026-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Deep SeapHOx V2 DO Functions (SBE 63)

### Core Functions

See the [Deep SeapHOx V2 instrument notes](../seaphox.md) for how the pH, CTD, and DO subsystems and their data streams fit together.

::: ion_functions.data.phsen_h_functions.dissolved_oxygen

#### History
| Date | Author | Change |
|---|---|---|
| 2026-05-15 | Christopher Wingard | Initial NumPy docstring; added to DO family documentation |

---

### Helper Functions

::: ion_functions.data.phsen_h_functions.convert_sbe63_thermistor

#### History
| Date | Author | Change |
|---|---|---|
| 2026-05-15 | Christopher Wingard | Initial NumPy docstring; added to DO family documentation |

---

## References

Aanderaa Data Instruments. *TD 269 Operating Manual Oxygen Optode 4330,
4831, 4835*. Bergen: Aanderaa Data Instruments.

Garcia, H. E. and Gordon, L. I. (1992). Oxygen solubility in seawater:
Better fitting equations. Limnology and Oceanography, 37(6), 1307-1312.

[OOI (2014). Data Product Specification for Oxygen Concentration from
"Stable" Instruments. Document Control Number 1341-00520.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf)

[OOI (2013). Data Product Specification for Fast Dissolved Oxygen.
Document Control Number 1341-00521.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00521_Data_Product_Spec_DOCONCF_OOI.pdf)

Owens, W. B. and Millard, R. C. (1985). A new algorithm for CTD oxygen
calibration. Journal of Physical Oceanography, 15, 621-631.

Sea-Bird Scientific. *SBE 43 Dissolved Oxygen Sensor User's Manual*.
Bellevue, WA: Sea-Bird Scientific.

Sea-Bird Scientific. *SBE 43F Dissolved Oxygen Sensor User's Manual*.
Bellevue, WA: Sea-Bird Scientific.

Sea-Bird Scientific. *SBE 63 Optical Dissolved Oxygen Sensor User's
Manual*. Bellevue, WA: Sea-Bird Scientific.

Sea-Bird Scientific. *Deep SeapHOx V2 Ocean CT(D)-pH-DO Sensor*.
Data sheet DS.53.May25. Bellevue, WA: Sea-Bird Scientific.

Uchida, H., Kawano, T., Kaneko, I., and Fukasawa, M. (2008). In situ
calibration of optode-based oxygen sensors. Journal of Atmospheric and
Oceanic Technology, 25, 2271-2281.
