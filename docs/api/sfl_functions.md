# SFL Functions

## Background

The Seafloor Properties (SFL) family covers four distinct instrument classes
deployed by OOI to characterize hydrothermal vent chemistry and seafloor
pressure: the Hydrothermal Vent Fluid In-situ Chemistry instrument (THSPH),
the Temperature-Resistivity Probe (TRHPH), the Seafloor Pressure instrument
(PRESF, Sea-Bird SBE 26plus), and the Seafloor Pressure instrument (PREST,
Sea-Bird SBE 54 Tsunameter). All SFL processing functions are implemented in
`sfl_functions.py`.

| Class | Hardware | Platform | Designator Meaning |
|---|---|---|---|
| THSPH | Custom multi-electrode chemistry sonde | Seafloor | Hydrothermal vent fluid in-situ chemistry |
| TRHPH | Custom temperature-resistivity probe | Seafloor | Temperature-resistivity probe for hydrothermal vents |
| PRESF | Sea-Bird SBE 26plus | Seafloor | Seafloor pressure (SBE 26plus) |
| PREST | Sea-Bird SBE 54 Tsunameter | Seafloor | Seafloor pressure (SBE 54) |

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00120 | [OOI (2014). Data Product Specification for Vent Fluid Temperature from THSPH.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00120_Data_Product_SPEC_THSPHTE_OOI.pdf) |
| 1341-00150 | [OOI (2014). Data Product Specification for Vent Fluid Temperature from TRHPH.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00150_Data_Product_Spec_TRHPHTE_OOI.pdf) |
| 1341-00160 | [OOI (2013). Data Product Specification for Vent Fluid Chloride Concentration.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00160_Data_Product_Spec_TRHPHCC_OOI.pdf) |
| 1341-00170 | [OOI (2013). Data Product Specification for Vent Fluid Oxidation-Reduction Potential (ORP).](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00170_Data_Product_Spec_TRHPHEH_OOI.pdf) |
| 1341-00190 | OOI (2014). Data Product Specification for Vent Fluid pH. Document Control Number 1341-00190. (Not released.) |
| 1341-00200 | OOI (2014). Data Product Specification for Vent Fluid Hydrogen Sulfide Concentration. Document Control Number 1341-00200. (Not released.) |
| 1341-00210 | OOI (2014). Data Product Specification for Vent Fluid Hydrogen Concentration. Document Control Number 1341-00210. (Not released.) |
| 1341-00230 | [OOI (2013). Data Product Specification for Seafloor Pressure from Sea-Bird SBE 26PLUS.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00230_Data_Product_Spec_SFLPRES_OOI.pdf) |
| 1341-00231 | [OOI (2012). Data Product Specification for Seafloor Pressure from SeaBird 54 Tsunameter.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00231_Data_Product_SPEC_SLFPRES_OOI.pdf) |

---

### THSPH -- Hydrothermal Vent Fluid In-situ Chemistry

The THSPH instrument is a custom multi-electrode chemistry sonde deployed
directly in hydrothermal vent fluids on the Regional Cabled Array. It
measures temperature and chemical concentrations in high-temperature
vent fluids using a combination of thermocouples, thermistors, and
electrochemical electrodes. The instrument outputs voltage and resistance 
data that are converted off-instrument to the L1 and L2 data products 
described below.

The THSPH contains eight channels. The temperature-relevant channels are:
a high-temperature thermocouple (TCH), a low-temperature thermocouple
(TCL), a reference thermistor at the cold-junction inside the sensor wand
(REF), and a board thermistor inside the electronics housing (INT). The
chemical measurement channels are: a YSZ (yttria-stabilized zirconia) pH
electrode, an AgCl reference electrode, a hydrogen electrode, and a sulfide
electrode.

All THSPH processing uses 5th-degree polynomial evaluation via `eval_poly`,
which implements Horner's algorithm for numerical stability. Calibration
coefficients for each sensor are supplied as 6-element arrays in
descending-degree order. All calibration coefficients are from factory
calibration sheets.

#### THSPHTE_L1 -- Vent Fluid Temperature from THSPH

THSPHTE_L1 is the vent fluid temperature ($^\circ$C) measured by the THSPH
instrument. Two "final" temperature products are produced -- THSPHTE-TH
(temperature at the high-temperature sample inlet) and THSPHTE-TL
(temperature near the vent) -- along with four intermediate temperature
sub-products used for instrument health monitoring and quality control:
THSPHTE-TCH, THSPHTE-TCL, THSPHTE-REF, and THSPHTE-INT.

The processing chain for the final temperature products (TH and TL) follows
seven steps, as specified in DPS 1341-00120:

**Step 1 -- Raw to engineering units.** Decimal thermocouple counts are
converted to engineering voltages [mV], and decimal thermistor counts are
converted to engineering resistances [$\Omega$]:

$$\begin{align}
V_{tc,eng} &= \frac{tc_{rawdec} \times 0.25 - 1024}{61606} \\
R_{ts,eng} &= \frac{10000 \times ts_{rawdec} \times 0.125}
              {2048 - ts_{rawdec} \times 0.125}
\end{align}$$

**Step 2 -- Engineering to lab-calibrated values.** A 5th-degree polynomial
correction (coefficients `e2l`) converts engineering values to lab-calibrated
voltages or resistances.

**Step 3 -- Lab-calibrated to scientific units.** For the thermocouples, a
5th-degree polynomial (coefficients `l2s`) converts calibrated voltage [mV] to
temperature [$^\circ$C]. For the thermistors, a 4th-degree polynomial (coefficients 
`l2s`) is evaluated at $\ln(R_{ts,actual})$ to produce an intermediate value
$pval$, and temperature is recovered as:

$$T_{ts} = \frac{1}{pval} - 273.15$$

**Step 4 -- Cold-junction correction.** The reference thermistor temperature
is converted to a thermocouple-equivalent voltage using `s2v` coefficients,
which is added to the calibrated thermocouple voltage before the final
polynomial evaluation to yield the absolute (cold-junction-corrected)
temperature at each thermocouple site.

**Step 5 -- Final linear calibration.** A linear correction (coefficients `s2f`)
is applied to the combined thermocouple-plus-thermistor temperature to yield the 
final L1 products TH and TL.

The four intermediate sub-products (TCH, TCL, REF, INT) are produced at
earlier stages of this chain and are useful for instrument diagnostics.
`sfl_thsph_temp_labcal_h` and `sfl_thsph_temp_labcal_l` implement the
engineering-to-lab-calibrated conversion for each thermocouple channel and
are called internally by the core temperature functions.

#### THSPHHC_L2, THSPHHS_L2, and THSPHPH_L2 -- Vent Fluid Chemistry

The THSPH chemical data products are computed from electrode potential
measurements referenced against the YSZ electrode and corrected for
temperature using the Nernst equation. All three products share the
temperature-dependent Nernst factor:

$$E_{nernst} = 1.9842 \times 10^{-4} \times (T + 273.15)$$

where $T$ is the THSPHTE-TH temperature in $^\circ$C. Raw electrode counts
are first converted to lab-calibrated voltages via `v_labcal`, which applies
a counts-to-engineering-values conversion followed by a 5th-degree
polynomial correction using `e2l` calibration coefficients.

**THSPHHC_L2** is the hydrogen concentration (mmol/kg) at the vent. The
measured H$_2$ electrode potential is computed as the difference between the
YSZ and H$_2$ lab-calibrated voltages. The log of hydrogen fugacity is then
derived from this potential and a temperature-dependent correction for the
electrode material (`arr_hgo`), and hydrogen concentration is recovered
from the log fugacity using equilibrium calibration coefficients
(`arr_logkfh2g`).

**THSPHHS_L2** is the hydrogen sulfide concentration (mmol/kg) at the vent.
The algorithm follows the same structure as THSPHHC, with the addition of
temperature-dependent coefficients for the H$_2$S electrode (`arr_eh2sg`)
and a fugacity-to-concentration quotient (`arr_yh2sg`).

**THSPHPH_L2** is the vent fluid pH (dimensionless). Four sub-products are
produced depending on data availability:

- **THSPHPH-PH**: measured AgCl reference electrode and chloride from
  TRHPHCC_L2.
- **THSPHPH-PH-ACL**: measured AgCl reference electrode, assumed chloride
  (default 250 mmol/kg).
- **THSPHPH-PH-NOREF**: theoretical reference electrode potential computed
  from temperature (`arr_agclref`), chloride from TRHPHCC_L2.
- **THSPHPH-PH-NOREF-ACL**: theoretical reference electrode potential,
  assumed chloride.

All four share the `calculate_vent_pH` helper. The pH is computed from the
calibrated electrode potential difference, temperature-dependent electrode
material corrections (`arr_hgo`, `arr_agcl`), and chloride activity. Chloride
activity is computed by `chloride_activity` as a polynomial function of
temperature, using four sets of 5th-degree coefficients (`arr_tac`,
`arr_tbc1`, `arr_tbc2`, `arr_tbc3`). When no chloride measurement is
available, a default of 250 mmol/kg is used. Out-of-range electrode
potentials ($E_{ph}$ outside [-0.7, 0.0] V) and out-of-range pH values
(outside [3.0, 7.0]) are set to NaN.

---

### TRHPH -- Temperature-Resistivity Probe

The TRHPH instrument is a custom probe placed directly in high-temperature
hydrothermal vent fluid to measure temperature, oxidation-reduction potential 
(ORP), and resistivity. It uses a Type K thermocouple for high-temperature 
measurement and a thermistor as a cold-junction reference, along with a Pt-Ag/AgCl 
ORP electrode pair and three resistivity circuits scaled to different concentration ranges.

#### TRHPHTE_L1 -- Vent Fluid Temperature from TRHPH

TRHPHTE_L1 is the vent fluid temperature ($^\circ$C) measured by the TRHPH
instrument. The algorithm, specified in DPS 1341-00150, first computes
thermistor temperature $T_{ts}$ from the thermistor voltage $V_{ts}$:

$$T_{ts} = 27.50133 - 17.2658 \times V_{ts} + \frac{15.83424}{V_{ts}}$$

The final temperature $T$ is then determined by three cases:

- When $V_{tc} \leq 0$: 

$$T = T_{ts}$$

- When $V_{tc} > 0$ and $T_{ts} > 10$ $^\circ$C:

$$T = (V_{tc} + c_3 {T_{ts}}^3 + c_2 {T_{ts}}^2 + c_1 T_{ts} + c_0)
      \times 244.97$$

- When $V_{tc} > 0$ and $0 < T_{ts} \leq 10$ $^\circ$C:

$$T = (V_{tc} + V_{tc} \times 244.97 \times tc_{slope}
       + T_{ts} \times ts_{slope}) \times 244.97$$

The polynomial coefficients $(c_3, c_2, c_1, c_0) = (-10^{-6},\ 7 \times
10^{-5},\ 0.0024,\ 0.015)$ are fixed constants specified in the DPS.
`tc_slope` and `ts_slope` are instrument-specific factory calibration
coefficients. The auxiliary product TRHPHTE-T\_TS-AUX (thermistor
temperature only) is computed by `sfl_trhph_vfl_thermistor_temp` and is
useful as an instrument diagnostic.

#### TRHPHEH_L1 -- Vent Fluid Oxidation-Reduction Potential

TRHPHEH_L1 is the vent fluid oxidation-reduction potential (ORP, mV)
measured by the Pt-Ag/AgCl electrode pair on the TRHPH instrument. As
specified in DPS 1341-00170, the conversion from raw sensor voltage to mV
corrects for the electronic offset and gain introduced by the A/D board:

$$ORP = \frac{V \times 1000 - offset}{gain}$$

where $V$ is the raw ORP voltage (TRHPHVO_L0) in volts (V), $offset$ is the
electronic offset calibration coefficient in mV, and $gain$ is the gain
multiplier calibration coefficient. The result is rounded to the nearest
integer mV. Because the reference electrode is not a Standard Hydrogen Electrode, 
ORP values from this instrument cannot be directly compared with standard in-situ 
ORP measurements; the primary use of this product is to quantify change in ORP 
with respect to time.

#### TRHPHCC_L2 -- Vent Fluid Chloride Concentration

TRHPHCC_L2 is the vent fluid chloride concentration (mmol/kg) computed from
TRHPH resistivity and temperature data. The algorithm, specified in DPS
1341-00160, uses a three-dimensional temperature-conductivity-chloride
calibration surface developed by Larson et al. (2007) and stored as the
`Larson_2007surface.mat` reference dataset (imported from
`sfl_functions_surface.py` as `tdat`, `sdat`, `cdat`).

The three resistivity voltage channels (V\_R1, V\_R2, V\_R3) are scaled to
different measurement ranges. The algorithm selects the optimal channel:

$$V_R = \begin{cases}
V_{R3} / 5 & \text{if } V_{R2} < 0.75 \\
V_{R2}     & \text{if } 0.75 \leq V_{R2} < 3.90 \\
V_{R1} \times 5 & \text{if } V_{R2} \geq 3.90
\end{cases}$$

Conductivity is computed as $C = 1 / V_R$. For each observation where
temperature falls within the bounds of the calibration surface (103 to
382 $^\circ$C), a conductivity isotherm is interpolated from the surface at
the observed temperature using `RectBivariateSpline`, and the measured
conductivity is mapped to chloride concentration. Values outside the
calibration surface bounds are returned as NaN. The result is converted from
mol/kg to mmol/kg and rounded to the nearest integer.

---

### PRESF/PREST -- Seafloor Pressure

The PRESF and PREST instrument classes both produce the SFLPRES_L1 Seafloor
Pressure data product (dbar). PRESF uses the Sea-Bird SBE 26plus, which
measures absolute pressure (hydrostatic plus atmospheric) using an internal
quartz crystal resonator and supports three sampling modes. PREST uses the
Sea-Bird SBE 54 Tsunameter, which streams real-time absolute pressure data.
Both instruments report the absolute pressure in psi and only require unit 
conversions.

#### SFLPRES_L1 -- Seafloor Pressure

SFLPRES_L1 is the absolute seafloor pressure (dbar) including both
hydrostatic and atmospheric pressure. Three sub-products correspond to the
three SBE 26plus operating modes:

**SFLPRES-RTIME** is produced from real-time ASCII output. The instrument
internally converts raw counts to pressure in psi using onboard calibration
coefficients. The L1 product is a direct unit conversion from psi to dbar:

$$p_{dbar} = p_{psi} \times 0.689475728$$

**SFLPRES-TIDE** is produced from post-recovery tide data. The raw pressure 
measurement `p_dec_tide` is converted to psi using instrument calibration 
coefficients $M$, $B$, and optional slope and offset corrections, then scaled
to dbar:

$$p_{psi} = slope \times \frac{p_{dec\_tide} - B}{M} + offset$$

$$p_{dbar} = p_{psi} \times 0.689475728$$

**SFLPRES-WAVE** is produced from post-recovery wave burst data. The pressure 
temperature compensation frequency (PTCF) and pressure frequency (PF) are derived 
from raw inputs, then pressure is computed using quartz transducer calibration 
coefficients as specified in DPS 1341-00230:

$$\begin{align}
U &= \frac{10^6}{PTCF} - U_0 \\
C &= C_1 + C_2 U + C_3 U^2 \\
D &= D_1 + D_2 \\
T_0 &= \frac{T_1 + T_2 U + T_3 U^2 + T_4 U^3}{10^6} \\
W &= 1 - {T_0}^2 \times PF^2 \\
p_{psi} &= slope \times [C \times W \times (1 - D \times W) + POffset]
            + offset
\end{align}$$

The temperature metadata product PRESTMP_L1 (seawater temperature in
$^\circ$C) is decoded from the tide record temperature number $t_0$ as:

$$t = \frac{t_0}{1000} - 10$$

The SBE 54 Tsunameter (PREST) transmits pressure already converted to psi
using onboard calibration coefficients; the L1 product is computed with the
same unit conversion as SFLPRES-RTIME and also represents the absolute seafloor
pressure..

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.sfl_functions.sfl_thsph_temp_th

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-01 | Russell Desiderio | Initial code |
| 2014-06-30 | Russell Desiderio | DPS modifications to cal equations implemented |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_temp_tl

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-01 | Russell Desiderio | Initial code |
| 2014-06-30 | Russell Desiderio | DPS modifications to cal equations implemented |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_temp_tch

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-01 | Russell Desiderio | Initial code |
| 2014-06-30 | Russell Desiderio | DPS modifications to cal equations implemented |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_temp_tcl

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-01 | Russell Desiderio | Initial code |
| 2014-06-30 | Russell Desiderio | DPS modifications to cal equations implemented |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_temp_ref

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-01 | Russell Desiderio | Initial code |
| 2014-06-30 | Russell Desiderio | DPS modifications to cal equations implemented |
| 2015-07-24 | Russell Desiderio | Added call to replace_fill_with_nan; cleaned up error-checking |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_temp_int

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-01 | Russell Desiderio | Initial code |
| 2014-06-30 | Russell Desiderio | DPS modifications to cal equations implemented |
| 2015-07-24 | Russell Desiderio | Added call to replace_fill_with_nan; cleaned up error-checking |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_hydrogen

#### Additional Notes
The DPS document for THSPHHC (1341-00210) was never publicly released.
The algorithm is documented from the code and code comments only.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_sulfide

#### Additional Notes
The DPS document for THSPHHS (1341-00200) was never publicly released.
The algorithm is documented from the code and code comments only.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_ph

#### Additional Notes
The DPS document for THSPHPH (1341-00190) was never publicly released.
The algorithm is documented from the code and code comments only.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2015-07-24 | Russell Desiderio | Incorporated calculate_vent_pH function |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_ph_acl

#### Additional Notes
The DPS document for THSPHPH (1341-00190) was never publicly released.
The algorithm is documented from the code and code comments only.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2015-07-24 | Russell Desiderio | Incorporated calculate_vent_pH function |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_ph_noref

#### Additional Notes
The DPS document for THSPHPH (1341-00190) was never publicly released.
The algorithm is documented from the code and code comments only.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2015-07-24 | Russell Desiderio | Incorporated calculate_vent_pH function |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_ph_noref_acl

#### Additional Notes
The DPS document for THSPHPH (1341-00190) was never publicly released.
The algorithm is documented from the code and code comments only.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2015-07-24 | Russell Desiderio | Incorporated calculate_vent_pH function |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_trhph_vfltemp

#### History
| Date | Author | Change |
|---|---|---|
| 2013-05-01 | Christopher Wingard | Initial code |
| 2014-02-27 | Russell Desiderio | Added documentation; implemented Horner's method |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_trhph_vflorp

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-28 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_trhph_chloride

#### History
| Date | Author | Change |
|---|---|---|
| 2013-05-01 | Christopher Wingard | Initial code |
| 2014-02-28 | Russell Desiderio | Modified code to better handle nans and fill values; added documentation |
| 2014-03-10 | Russell Desiderio | Removed unnecessary np.vectorized wrapper; improved speed |
| 2014-03-26 | Russell Desiderio | Incorporated optimization from Chris Fortin; sped up execution by factor of 5 |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_sflpres_rtime

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-31 | Craig Risien | Initial code |
| 2014-09-23 | Christopher Wingard | Minor edits |
| 2015-07-22 | Russell Desiderio | Removed replace_fill_with_nan call (no integer inputs) |
| 2023-08-15 | Samuel Dahlberg | Removed use of Numexpr library |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_sflpres_tide

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-23 | Christopher Wingard | Initial code |
| 2015-07-22 | Russell Desiderio | Added call to replace_fill_with_nan |
| 2023-08-15 | Samuel Dahlberg | Removed use of Numexpr library |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_sflpres_wave

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-23 | Christopher Wingard | Initial code |
| 2015-07-20 | Russell Desiderio | Modified code to accept p_dec_wave as 2D array |
| 2015-07-22 | Russell Desiderio | Added call to replace_fill_with_nan |
| 2023-08-15 | Samuel Dahlberg | Removed use of Numexpr library |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_sbe26plus_prestmp

#### History
| Date | Author | Change |
|---|---|---|
| 2015-10-28 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.sfl_functions.sfl_thsph_temp_labcal_h

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-30 | Russell Desiderio | Initial code |
| 2015-07-24 | Russell Desiderio | Added call to replace_fill_with_nan |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_thsph_temp_labcal_l

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-30 | Russell Desiderio | Initial code |
| 2015-07-24 | Russell Desiderio | Added call to replace_fill_with_nan |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.calculate_vent_pH

#### Additional Notes
The DPS document for THSPHPH (1341-00190) was never publicly released.
The algorithm is documented from the code and code comments only.

#### History
| Date | Author | Change |
|---|---|---|
| 2015-07-24 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.chloride_activity

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.v_labcal

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2015-07-22 | Russell Desiderio | Added call to replace_fill_with_nan |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.nernst

#### History
| Date | Author | Change |
|---|---|---|
| 2014-07-08 | Russell Desiderio | Initial code |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.eval_poly

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-01 | Russell Desiderio | Initial code (no arrays) |
| 2014-07-02 | Russell Desiderio | 2D calibration coefficient array implementation |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.sfl_functions.sfl_trhph_vfl_thermistor_temp

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-28 | Russell Desiderio | Initial code |
| 2015-01-06 | Russell Desiderio | Documented product as TRHPHTE-T_TS-AUX |
| 2026-06-10 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

Larson, B.I., Olson, E.J., and Lilley, M.D. (2007). In situ measurement of
dissolved chloride in high temperature hydrothermal vent fluids with a
diamond anvil cell. *Geochimica et Cosmochimica Acta*, 71(11), 2673-2683.

[OOI (2014). Data Product Specification for Vent Fluid Temperature from
THSPH. Document Control Number 1341-00120.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00120_Data_Product_SPEC_THSPHTE_OOI.pdf)

OOI (2014). Data Product Specification for Vent Fluid pH. Document Control
Number 1341-00190. (Not released.)

OOI (2014). Data Product Specification for Vent Fluid Hydrogen Sulfide
Concentration. Document Control Number 1341-00200. (Not released.)

OOI (2014). Data Product Specification for Vent Fluid Hydrogen
Concentration. Document Control Number 1341-00210. (Not released.)

[OOI (2014). Data Product Specification for Vent Fluid Temperature from
TRHPH. Document Control Number 1341-00150.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00150_Data_Product_Spec_TRHPHTE_OOI.pdf)

[OOI (2013). Data Product Specification for Vent Fluid Chloride
Concentration. Document Control Number 1341-00160.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00160_Data_Product_Spec_TRHPHCC_OOI.pdf)

[OOI (2013). Data Product Specification for Vent Fluid Oxidation-Reduction
Potential (ORP). Document Control Number 1341-00170.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00170_Data_Product_Spec_TRHPHEH_OOI.pdf)

[OOI (2013). Data Product Specification for Seafloor Pressure from Sea-Bird
SBE 26PLUS. Document Control Number 1341-00230.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00230_Data_Product_Spec_SFLPRES_OOI.pdf)

[OOI (2012). Data Product Specification for Seafloor Pressure from SeaBird
54 Tsunameter. Document Control Number 1341-00231.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00231_Data_Product_SPEC_SLFPRES_OOI.pdf)
