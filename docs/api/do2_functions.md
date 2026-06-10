# Dissolved Oxygen Functions

## Background

The Ocean Observatories Initiative deploys two families of dissolved oxygen
instruments across its moored, profiling, and mobile platforms. A third
dissolved oxygen sensor is integrated into the Deep SeapHOx V2 pH instrument
(PHSEN-G and PHSEN-H). The table below lists the OOI instrument classes
covered by this module.

| Class | Hardware | Platform type | Designator meaning |
|---|---|---|---|
| DOSTA | Aanderaa Optode 4831 | Moored (fixed depth) and gliders | DO, Stable |
| DOFST-A | Sea-Bird Scientific SBE 43 | Moored profiler (CTDPF-A/B) | DO, Fast, SBE 43 |
| DOFST-K | Sea-Bird Scientific SBE 43F | Moored profiler (CTDPF-C/K/L) | DO, Fast, SBE 43F |
| PHSEN-G/H | Sea-Bird Scientific Deep SeapHOx V2 (SBE 63) | Moored (fixed depth) | pH Sensor (integrated DO) |

`do2_functions.py` processes data from the DOSTA and DOFST instrument classes.
The dissolved oxygen functions for the PHSEN-G/H integrated SBE 63 sensor live
in `phsen_h_functions.py`; they are documented here because the SBE 63 is a
dissolved oxygen sensor. See [Deep SeapHOx V2](../seaphox.md) for instrument
architecture context.

All calibration coefficients are from factory calibration sheets supplied with
individual instruments. Within the OOI data system, dissolved oxygen functions
fall under the Water Column science regime and the Dissolved O$_2$ category.

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00520 | [OOI (2014). Data Product Specification for Oxygen Concentration from "Stable" Instruments.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00520_Data_Product_SPEC_DOCONCS_OOI.pdf) |
| 1341-00521 | [OOI (2013). Data Product Specification for Fast Dissolved Oxygen.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00521_Data_Product_Spec_DOCONCF_OOI.pdf) |

---

### DOCONCS_L1 and DOXYGEN_L2 — Stable Dissolved Oxygen (DOSTA)

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
product was named DOCONCS_L2 in early OOI documentation and renamed DOXYGEN
in DPS 1341-00520 version 1-02 (2014-04-11); `do2_salinity_correction`
produces DOXYGEN_L2.

#### L0 Data Streams

The DOSTA reaches the OOI data system through three distinct hardware
configurations, each producing a different L0 format:

**Analog output to CTD voltage channel** — The optode outputs $P_t$ and
$T_{opt}$ as 0-5 V analog signals routed into CTD auxiliary voltage channels.
The voltage signals are parsed as DOCONCS-VLT_L0 (phase) and a separate
temperature voltage, then converted to degrees and $^\circ$C:

$$P_t\ [\text{deg}] = 10 + 12 \times P_t\ [\text{V}]$$

$$T_{opt}\ [^\circ\text{C}] = -5 + 8 \times T_{opt}\ [\text{V}]$$

These fixed conversion coefficients are universal for all Aanderaa optodes
(Shawn Sneddon, Xylem-Aanderaa).

**Digital output via RS-232 to CTD** — The optode reports an oxygen
concentration directly in internal units transmitted over RS-232 to an SBE
16Plus V2 CTD, parsed as DOCONCS-CNT_L0:

$$\text{DOCONCS_L1}\ [\mu\text{mol/L}] = \frac{\text{DOCONCS-CNT_L0}}{10000} - 10$$

**Autonomous digital output** — The optode operates independently and reports
$P_t$ (DOCONCS-DEG_L0) and $T_{opt}$ directly in degrees and $^\circ$C,
requiring no voltage conversion before the SVU equation.

#### DOCONCS_L1 — Stern-Volmer-Uchida Equation

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
$\text{csv}_7$ are factory calibration coefficients.

A secondary calibration (conc_coef = [offset, slope]) is applied after the
SVU equation:

$$O_2 = \text{conc_coef}[0] + \text{conc_coef}[1] \times O_2$$

Aanderaa applies this two-point correction after refurbishment; for new
optodes or new SVU foil determinations conc_coef defaults to [0, 1].

The result is DOCONCS_L1 in $\mu$mol L$^{-1}$. Note: the OOI documentation
incorrectly lists the units for this product as $\mu$mol kg$^{-1}$.

The optode foil thermistor temperature ($T_{opt}$) must be used as the
temperature input to the SVU equation, not CTD temperature. The thermistor is
situated directly at the foil; the SVU calibration coefficients are derived to
compensate for changes in oxygen permeability through the foil as a function
of its temperature. On gliders, differences between CTD and optode temperature
of 1$^\circ$C correspond to approximately 5% differences in calculated oxygen
concentration.

#### DOXYGEN_L2 — Salinity and Pressure Correction

DOCONCS_L1 ($\mu$mol L$^{-1}$) is corrected for salinity and pressure using
CTD data from the co-located instrument to produce DOXYGEN_L2 ($\mu$mol
kg$^{-1}$). The correction has three sequential steps.

**Step 1 — Volume to mass conversion.** Potential density $\rho$ is computed
at reference pressure $p_{ref} = 0$ dbar from absolute salinity and
conservative temperature (TEOS-10 GSW library):

$$S_A = \text{gsw.SA_from_SP}(S_P, p, \text{lon}, \text{lat})$$

$$\Theta = \text{gsw.CT_from_t}(S_A, T, p)$$

$$\rho = \text{gsw.rho}(S_A, \Theta, p_{ref})$$

$$O_{2,\text{mass}} = \frac{1000 \times O_2 [\mu\text{mol/L}]}{\rho\ [\text{kg/m}^3]}$$

**Step 2 — Pressure correction** (Uchida et al. 2008):

$$O_{2,\text{pres}} = \left(1 + \frac{0.032 \times p}{1000}\right) \times O_{2,\text{mass}}$$

**Step 3 — Salinity correction** (Garcia and Gordon 1992, Table 1,
combined fit):

$$t_s = \ln\!\left(\frac{298.15 - T}{273.15 + T}\right)$$

$$B(t_s) = B_0 + B_1 t_s + B_2 {t_s}^2 + B_3 {t_s}^3$$

$$\text{DOXYGEN_L2} = \exp\!\left[(S_P - S_{ref}) \times B(t_s) + C_0 \times ({S_P}^2 - {S_{ref}}^2)\right] \times O_{2,\text{pres}}$$

where $B_0 = -6.24097 \times 10^{-3}$, $B_1 = -6.93498 \times 10^{-3}$,
$B_2 = -6.90358 \times 10^{-3}$, $B_3 = -4.29155 \times 10^{-3}$,
$C_0 = -3.11680 \times 10^{-7}$, and $S_{ref}$ is the reference salinity
configured in the Aanderaa optode firmware (typically 0 or 35).

---

### DOCONCF_L2 — Fast Dissolved Oxygen (DOFST)

The Sea-Bird Scientific SBE 43 (DOFST-A) and SBE 43F (DOFST-K) are Clark-type
polarographic dissolved oxygen sensors integrated into profiling CTD systems.
The SBE 43 is mounted on the SBE 16Plus V2 CTD (CTDPF-A/B) and reports
a voltage signal; the SBE 43F is mounted on the SBE 52-MP profiling CTD
(CTDPF-C/K/L) and reports a frequency signal.

DOCONCF_L2 is the dissolved oxygen concentration ($\mu$mol kg$^{-1}$)
corrected for temperature, salinity, and pressure. There is no intermediate
L1 product; DOCONCF_L0 (counts representing voltage or frequency) is converted
directly to DOCONCF_L2.

#### L0 Input Conversion

**SBE 43 (DOFST-A):** The raw hex voltage counts are converted to volts:

$$V = \frac{\text{DOCONCF-CNT_L0}}{13107}$$

**SBE 43F (DOFST-K):** The raw hex frequency value is used directly as
$F$ [Hz].

#### DOCONCF_L2 Algorithm

Sea-Bird uses an algorithm based on Owens and Millard (1985). With the
derivative term removed by setting tau = 0 (recommended by Sea-Bird and the
OOI DPS to prevent amplification of residual noise in deep water):

**Step 1 — Oxygen solubility** (Garcia and Gordon 1992, Table 1, 1st column):

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

**Step 2 — Intermediate oxygen concentration** [mL L$^{-1}$]:

$$\text{DO_int} = \text{Soc} \times (x + \text{offset}) \times \text{Oxsol}(T, S)
\times (1 + AT + BT^2 + CT^3) \times e^{EP/K}$$

where $x$ is $V$ (SBE 43) or $F$ (SBE 43F), $K = T + 273.15$ is absolute
temperature, and Soc, offset, A, B, C, E are factory calibration coefficients.

**Step 3 — Convert to $\mu$mol kg$^{-1}$:**

Potential density $\rho_\theta$ is computed using `gsw.pot_rho_t_exact` with
a reference pressure of 0 dbar:

$$\text{DOCONCF_L2} = \frac{\text{DO_int} \times 44660}{\rho_\theta}$$

---

### DOXYGEN_L2 — SeapHOx Dissolved Oxygen

The Sea-Bird Scientific Deep SeapHOx V2 (PHSEN-G and PHSEN-H) integrates an
SBE 63 optical dissolved oxygen sensor alongside an SBE 37 MicroCAT CTD and
a Deep SeaFET pH sensor. The SBE 63 uses a Stern-Volmer fluorescence quenching
approach with a calibration model distinct from the DOSTA Aanderaa optode.
Processing depends on the data delivery method.

**Streamed and recovered_instrument**: Raw phase output ($\mu$s) is converted
directly to DOXYGEN_L2 ($\mu$mol kg$^{-1}$) via `dissolved_oxygen` in
`phsen_h_functions.py`. Salinity and pressure corrections and the conversion
from mL L$^{-1}$ to $\mu$mol kg$^{-1}$ using potential density are all
applied within a single function call.

**Telemetered**: The instrument reports dissolved oxygen in mL L$^{-1}$
(DOCONCS_L1). The salinity and pressure correction to produce DOXYGEN_L2
follows the same `do2_salinity_correction` path as the DOSTA optode. This
processing path is not yet implemented in `do2_functions.py`.

See [Deep SeapHOx V2](../seaphox.md) for instrument architecture context.

---

## Core Functions

::: ion_functions.data.do2_functions.do2_SVU

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-26 | Stuart Pearce | Initial code |
| 2015-04-10 | Russell Desiderio | Revised for CI calibration coefficient implementation |
| 2015-08-04 | Russell Desiderio | Added documentation |
| 2015-08-10 | Russell Desiderio | Added conc_coef calibration array to argument list |
| 2015-10-28 | Russell Desiderio | Added atleast_2d handling for 1D conc_coef arrays |
| 2023-08-15 | Samuel Dahlberg | Renamed local variables to follow naming convention |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.do2_salinity_correction

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-26 | Stuart Pearce | Initial code |
| 2015-08-04 | Russell Desiderio | Added Garcia-Gordon reference |
| 2021-12-16 | Stuart Pearce | Added salinity reference parameter |
| 2023-08-15 | Samuel Dahlberg | Removed numexpr; replaced pygsw with gsw |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.dofst_calc

#### History
| Date | Author | Change |
|---|---|---|
| 2013-08-20 | Stuart Pearce | Initial code |
| 2015-08-04 | Russell Desiderio | Added Garcia-Gordon reference |
| 2023-08-15 | Samuel Dahlberg | Added freq variable for CGSN compatibility |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.phsen_h_functions.dissolved_oxygen

#### History
| Date | Author | Change |
|---|---|---|
| 2025-05-15 | Christopher Wingard | Initial NumPy docstring; added to DO family documentation |

---

## Helper Functions

::: ion_functions.data.do2_functions.dosta_phase_volt_to_degree

#### History
| Date | Author | Change |
|---|---|---|
| 2015-08-04 | Russell Desiderio | Initial code |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.dosta_Topt_volt_to_degC

#### History
| Date | Author | Change |
|---|---|---|
| 2015-08-04 | Russell Desiderio | Initial code |
| 2023-08-15 | Samuel Dahlberg | Renamed local variables to follow naming convention |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.o2_counts_to_uM

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-26 | Stuart Pearce | Initial code |
| 2015-04-10 | Russell Desiderio | Added documentation and fill value handling |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.phsen_h_functions.convert_sbe63_thermistor

#### History
| Date | Author | Change |
|---|---|---|
| 2025-05-15 | Christopher Wingard | Initial NumPy docstring; added to DO family documentation |

---

## Wrapper Functions

::: ion_functions.data.do2_functions.do2_dofst_volt

#### History
| Date | Author | Change |
|---|---|---|
| 2013-08-20 | Stuart Pearce | Initial code |
| 2015-08-05 | Russell Desiderio | Added fill value conversion to NaN |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.do2_functions.do2_dofst_frequency

#### History
| Date | Author | Change |
|---|---|---|
| 2013-08-20 | Stuart Pearce | Initial code |
| 2015-08-05 | Russell Desiderio | Added fill value conversion to NaN |
| 2025-05-15 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

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
