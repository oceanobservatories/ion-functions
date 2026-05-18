# pH Functions

## Background

The Ocean Observatories Initiative deploys two families of pH instruments
across its moored platforms. The table below lists the OOI instrument classes
covered by this module.

| Class | Hardware | Platform | Designator meaning |
|---|---|---|---|
| PHSEN-A/D/E/F | Sunburst SAMI-pH | Moored (fixed depth) | pH Sensor |
| PHSEN-G/H | Sea-Bird Scientific Deep SeapHOx V2 | Moored (fixed depth) | pH Sensor |

`ph_functions.py` processes data from the [Sunburst Sensors SAMI-pH](https://www.sunburstsensors.com/products/oceanographic-ph-sensor.html)
(PHSEN-A/D/E/F) and computes the L2 pH of seawater data product 
(PHWATER_L2). `phsen_h_functions.py` processes data from the [Sea-Bird 
Scientific Deep SeapHOx V2](https://www.seabird.com/products/deep-seaphox-v2-ph-sensor) 
(PHSEN-G and PHSEN-H) and computes the same PHWATER_L2 data product. Both 
instruments measure pH on the total hydrogen ion scale ($\text{pH}_T$). All 
calibration coefficients are from factory calibration sheets supplied with 
individual instruments. Within the OOI data system, pH functions fall under
the Water Column science regime and the pH category.

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00510 | [OOI (2012). Data Product Specification for pH of Seawater.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00510_Data_Product_SPEC_PHWATER_OOI.pdf) |

---

### PHWATER_L2 — pH of Seawater

PHWATER_L2 is the pH of seawater on the total hydrogen ion scale
($\text{pH}_T$). Both PHSEN-A/D/E/F and PHSEN-G/H produce PHWATER_L2 directly 
from their raw L0 data; there is no intermediate L1 product.

#### Sunburst Sensors SAMI-pH

The SAMI-pH measures pH using a colorimetric reaction with the pH indicator
meta-Cresol Purple (mCP). A seawater sample is pumped through a flow cell and
injected with a pulse of indicator solution. Two LEDs illuminate the
indicator-sample mixture at 434 nm and 578 nm — the peak absorbance wavelengths
of the protonated ($\text{HI}^-$) and deprotonated ($\text{I}^{2-}$) forms of
mCP, respectively. The ratio of the absorbances at these two wavelengths is used
to compute pH on the total scale.

**L0 Inputs**

Each SAMI-II data record contains two types of raw light measurements:

- **Blank measurements** — 4 sets of 4 interleaved readings
  $[\text{ref}_{434},\ \text{sig}_{434},\ \text{ref}_{578},\ \text{sig}_{578}]$
  collected while pumping pure seawater (16 values total).
- **Sample measurements** — 23 sets of 4 interleaved readings in the same
  arrangement, collected while pumping the seawater–indicator mixture
  (92 values total).

Raw signal intensities range from 0 to 4096 counts.

**Thermistor Temperature**

The raw thermistor count is converted to temperature in $^\circ$C. The
conversion depends on the SAMI hardware generation (12-bit or 14-bit ADC).
For 12-bit hardware (full-scale count $= 4096$):

$$r_t = \ln\left(\frac{\text{Therm}}{4096 - \text{Therm}} \times 17400\right)$$

For 14-bit hardware (full-scale count $= 16384$):

$$r_t = \ln\left(\frac{\text{Therm}}{16384 - \text{Therm}} \times 17400\right)$$

In both cases temperature is then:

$$T_C = \frac{1}{0.0010183 + 0.000241 \times r_t + 1.5 \times 10^{-7} \times {r_t}^3} - 273.15$$

The 14-bit ADC is a post-DPS hardware change not described in DPS 1341-00510.

**Blank Normalization**

The blank signal intensity ratio at each wavelength is computed from the 4
blank measurement sets and averaged:

$$B_{434} = \frac{1}{4}\sum_{k=1}^{4} \frac{\text{sig}_{434,k}}{\text{ref}_{434,k}}$$

$$B_{578} = \frac{1}{4}\sum_{k=1}^{4} \frac{\text{sig}_{578,k}}{\text{ref}_{578,k}}$$

**Absorbance**

The blank-corrected absorbances at each wavelength are computed using
Beer's Law:

$$A_{434} = -\log_{10}\left(\frac{\text{sig}_{434}}{\text{ref}_{434}}\right) - \left(-\log_{10} B_{434}\right)$$

$$A_{578} = -\log_{10}\left(\frac{\text{sig}_{578}}{\text{ref}_{578}}\right) - \left(-\log_{10} B_{578}\right)$$

The absorbance ratio is:

$$R = \frac{A_{578}}{A_{434}}$$

**Temperature-Dependent Molar Absorptivities**

The molar absorptivities of the two indicator forms are adjusted for
temperature using the factory-supplied reference values $e_{a434}$,
$e_{b434}$, $e_{a578}$, $e_{b578}$ at a reference temperature of
24.788 $^\circ$C:

$$E_{a434} = e_{a434} - 26 \times (T_C - 24.788)$$

$$E_{a578} = e_{a578} + (T_C - 24.788)$$

$$E_{b434} = e_{b434} + 12 \times (T_C - 24.788)$$

$$E_{b578} = e_{b578} - 71 \times (T_C - 24.788)$$

The absorptivity ratios used in the pH equation are:

$$e_1 = E_{a578} / E_{a434}, e_2 = E_{b578} / E_{a434}, e_3 = E_{b434} / E_{a434}$$

**pKa**

The apparent dissociation constant of mCP is computed from temperature and
salinity (Clayton and Byrne, 1993):

$$pK'_a = \frac{1245.69}{T_C + 273.15} + 3.8275 + 0.0021 \times (35 - S)$$

where $S$ is the seawater practical salinity from a co-located CTD (default 
35.0 if no CTD data are available).

**PHWATER_L2 Calculation**

Point-by-point pH values are computed for each of the 23 sample measurement
sets:

$$\text{pH}_\text{point} = pK'_a + \log_{10}\left(\frac{R - e_1}{e_2 - R \times e_3}\right)$$

Indicator concentrations for the protonated and deprotonated forms are:

$$[\text{HI}^-] = \frac{A_{434} \times E_{b578} - A_{578} \times E_{b434}}{E_{a434} \times E_{b578} - E_{b434} \times E_{a578}}$$

$$[\text{I}^{2-}] = \frac{A_{578} \times E_{a434} - A_{434} \times E_{a578}}{E_{a434} \times E_{b578} - E_{b434} \times E_{a578}}$$

$$C_\text{ind} = [\text{HI}^-] + [\text{I}^{2-}]$$

The first 5 of the 23 measurement sets are discarded. From the remaining 18,
the 8 consecutive points with the highest linear $R^2$ between $C_\text{ind}$
and $\text{pH}_\text{point}$ are selected. The final PHWATER_L2 value is
the y-intercept of a linear regression of $\text{pH}_\text{point}$ on
$C_\text{ind}$ through those 8 points — i.e., the pH extrapolated to zero
indicator concentration.

An impurity correction is applied when the extrapolated pH exceeds 8.2:

$$\text{pH}_\text{final} = \text{pH} \times \text{ind_slp} + \text{ind_off}$$

where `ind_slp` and `ind_off` are instrument-specific correction factors
not described in DPS 1341-00510.

Output accuracy: $\pm 0.01$ pH units; precision $\pm 0.005$ pH units
(DPS 1341-00510, §4.4). Algorithm results are valid between 0 and 35
$^\circ$C and at salinities of $35 \pm 1$; salinity corrections
from a co-located CTD extend the valid salinity range (DPS 1341-00510,
§3.3).

---

#### Sea-Bird Scientific Deep SeapHOx V2

The Deep SeapHOx V2 combines the Deep SeaFET V2 ISFET pH sensor with the
Sea-Bird Electronics SBE 37-SMP-ODO MicroCAT CTD+DO sensor. PHSEN-G and
PHSEN-H differ only in the pressure rating of their strain-gauge pressure
sensor (Sea-Bird Scientific Deep SeapHOx V2 data sheet, DS53, May 2025).

The ISFET external electrochemical cell exhibits a Nernstian response to pH
and is sensitive to chloride activity. The raw ISFET voltage is digitized
by a 23-bit ADC with a 2.5 V reference and unity gain and converted to
volts before the pH calculation (Sea-Bird Scientific Application Note 99).

**PHWATER_L2 Calculation**

The Nernst factor is computed from temperature and fundamental constants
(Application Note 99):

$$S_\text{nernst} = \frac{R \times T \times \ln(10)}{F}$$

where $R = 8.3144621\ \text{J} \times (\text{mol} \times \text{K})^{-1}$,
$T$ is temperature in K, and $F = 96485.365\ \text{C} \times \text{mol}^{-1}$.

The pressure response of the sensor is modeled by a 6th-order polynomial
in pressure $P$ (dbar) using factory coefficients $f_1$ through $f_6$
(Application Note 99):

$$f(P) = f_1 P + f_2 P^2 + f_3 P^3 + f_4 P^4 + f_5 P^5 + f_6 P^6$$

The coefficient $f_0$ is captured in $k_0$ and is not used separately.

Total pH on the total hydrogen ion scale is then (Application Note 99,
Johnson et al. 2016, Johnson et al. 2017):

$$\begin{align}
pH_T &= \frac{V_\text{FET/REF} - k_0 - k_2 \times t - f(P)}{S_\text{nernst}} \\
&\quad + \log_{10}(Cl_T) + 2 \times \log_{10}(\gamma_{\pm\text{HCl}})_{T\&P} \\
&\quad - \log_{10}\left(1 + \frac{S_T}{K_{S,T\&P}}\right) \\
&\quad - \log_{10}\left(\frac{1000 - 1.005 \times S}{1000}\right)
\end{align}$$

where $t$ is temperature in $^\circ$C, $S$ is practical salinity,
$Cl_T$ is total chloride, $(\gamma_{\pm\text{HCl}})_{T\&P}$ is the HCl
activity coefficient corrected for temperature and pressure, $S_T$ is
total sulfate, and $K_{S,T\&P}$ is the acid dissociation constant of
$\text{HSO}_4^-$ corrected for temperature and pressure.

The intermediate quantities are computed as follows (Application Note 99):

**Total chloride** (Dickson et al. 2007):

$$Cl_T = \frac{0.99889}{35.453} \times \frac{S}{1.80655} \times \frac{1000}{1000 - 1.005 \times S}$$

**Sample ionic strength** (Dickson et al. 2007):

$$I = \frac{19.924 \times S}{1000 - 1.005 \times S}$$

**Debye-Hückel constant** (Khoo et al. 1977):

$$A_{DH} = 3.4286 \times 10^{-6} \times t^2 + 6.7503 \times 10^{-4} \times t + 0.49172143$$

**HCl activity coefficient — temperature only** (Khoo et al. 1977):

$$\log(\gamma_{\pm\text{HCl}})_T = \frac{-A_{DH} \times \sqrt{I}}{1 + 1.394 \times \sqrt{I}} + (0.08885 - 0.000111 \times t) \times I$$

**Partial molal volume of HCl** (Millero 1983):

$$\bar{V}_\text{HCl} = 17.85 + 0.1044 \times t - 0.0001316 \times t^2$$

**HCl activity coefficient — temperature and pressure** (Johnson et al. 2017):

$$\log(\gamma_{\pm\text{HCl}})_{T\&P} = \log(\gamma_{\pm\text{HCl}})_T + \frac{\bar{V}_\text{HCl} \times p}{2 \times \ln(10) \times R \times T \times 10}$$

where $p$ is pressure in bar.

**Total sulfate** (Dickson et al. 2007):

$$S_T = \frac{0.1400}{96.062} \times \frac{S}{1.80655}$$

**Acid dissociation constant of HSO$_4$$^-$ at $T$** (Dickson et al. 2007):

$$K_S = (1 - 0.001005 \times S) \times \exp(adc)$$

where $adc$ is:

$$\begin{align}
adc &= (\frac{-4276.1}{T} + 141.328 - 23.093 \times \ln(T) \\
&\quad + \left(\frac{-13856}{T} + 324.57 - 47.986\ln(T)\right) \\
&\quad \times \sqrt{I} + \left(\frac{35474}{T} - 771.54 + 114.723 \times \ln(T)\right) \\
&\quad \times I - \frac{2698}{T} \times I^{1.5} + \frac{1776}{T} \times I^2
\end{align}$$

**Partial molal volume of $\text{HSO}_4^-$** (Millero 1983):

$$\bar{V}_S = -18.03 + 0.0466 \times t + 0.000316 \times t^2$$

**Compressibility of $\text{HSO}_4^-$** (Millero 1983):

$$\bar{K}_S = \frac{-4.53 + 0.09 \times t}{1000}$$

**Acid dissociation constant corrected for $T$ and $P$** (Millero 1982):

$$K_{S,T\&P} = K_S \times \exp\left(\frac{-\bar{V}_S \times p + 0.5 \times \bar{K}_S \times p^2}{R \times T \times 10}\right)$$

Calibration coefficients $k_0$, $k_2$, and $f_1$–$f_6$ are from factory
calibration sheets supplied by Sea-Bird Scientific with each instrument.

Output accuracy: $\pm 0.05$ pH units; resolution 0.004 pH units; typical
stability 0.003 pH units per month (Sea-Bird Scientific Deep SeapHOx V2
data sheet, DS53, May 2025).

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.ph_functions.ph_calc_phwater

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-19 | Christopher Wingard | Initial code |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.phsen_h_functions.ph_total

#### History
| Date | Author | Change |
|---|---|---|
| 2026-01-21 | Samuel Dahlberg | Initial code, adapted from Christopher Wingard's ph_total in cgsn-processing |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.ph_functions.ph_434_intensity

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-19 | Christopher Wingard | Initial code |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ph_functions.ph_578_intensity

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-19 | Christopher Wingard | Initial code |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ph_functions.ph_thermistor

#### History
| Date | Author | Change |
|---|---|---|
| 2014-05-01 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr; added default for sami_bits |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.ph_functions.ph_battery

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-19 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Added 14-bit hardware support |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.phsen_h_functions.convert_ph_voltage_counts

#### History
| Date | Author | Change |
|---|---|---|
| 2026-01-21 | Samuel Dahlberg | Initial code, adapted from Sea-Bird Scientific processing library |
| 2026-04-20 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

Byrne, R. H., Robertbaldo, G., Thompson, S. W. and Chen, C. T. A. (1988).
Seawater pH measurements — an at-sea comparison of spectrophotometric and
potentiometric methods. Deep-Sea Research Part A, 35(8): 1405–1410.

Clayton, T. D. and Byrne, R. H. (1993). Spectrophotometric seawater pH
measurements — total hydrogen-ion concentration scale calibration of
m-Cresol Purple and at-sea results. Deep-Sea Research Part I, 40(10):
2115–2129.

Dickson, A. G., Sabine, C. L., and Christian, J. R. (2007). Guide to Best
Practices for Ocean CO2 Measurements. PICES Special Publication 3, IOCCP
Report No. 8.

Johnson, K. S., Jannasch, H. W., Coletti, L. J., Elrod, V. A., Martz, T. R.,
Takeshita, Y., Carlson, R. J., and Connery, J. G. (2016). Deep-Sea DuraFET:
A pressure tolerant pH sensor designed for global sensor networks. Analytical
Chemistry, 88: 3249–3256.

Johnson, K. S., Plant, J. N., and Maurer, T. L. (2017). Processing BGC-Argo
pH data at the DAC level. BGC-Argo document.

Khoo, K. H., Ramette, R. W., Culberson, C. H., and Bates, R. G. (1977).
Determination of hydrogen ion concentrations in seawater from 5 C to 40 C:
standard potentials at salinities 20 to 45%. Analytical Chemistry, 49:
29–34.

Liu, X., Patsavas, M. C., and Byrne, R. H. (2011). Purification and
characterization of meta-cresol purple for spectrophotometric seawater pH
measurements. Environmental Science and Technology, 45: 4862–4868.

Martz, T. R., Connery, J. G., and Johnson, K. S. (2010). Testing the
Honeywell Durafet for seawater pH applications. Limnology and Oceanography:
Methods, 8: 172–184.

Martz, T. R., Carr, J. J., French, C. R., and DeGrandpre, M. D. (2003). A
submersible autonomous sensor for spectrophotometric pH measurements of
natural waters. Analytical Chemistry, 75: 1844–1850.

Millero, F. J. (1982). The effect of pressure on the solubility of minerals
in water and seawater. Geochimica et Cosmochimica Acta, 46: 11–22.

Millero, F. J. (1983). In Chemical Oceanography; Riley, J. P., Chester, R.,
Eds.; Academic Press: London, Vol. 8, pp 1–88.

[OOI (2012). Data Product Specification for pH of Seawater. Document Control
Number 1341-00510.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00510_Data_Product_SPEC_PHWATER_OOI.pdf)

Sea-Bird Scientific (2025). Deep SeapHOx V2 Data Sheet. Document DS53,
May 2025.

Sea-Bird Scientific. Application Note 99: Calculating pH from ISFET pH
Sensors. SeaFET V2, Shallow SeapHOx V2, Deep SeapHOx V2, Floats.

Seidel, M. P., DeGrandpre, M. D., and Dickson, A. G. (2008). A sensor for
in situ indicator-based measurements of seawater pH. Marine Chemistry,
109: 18–28.
