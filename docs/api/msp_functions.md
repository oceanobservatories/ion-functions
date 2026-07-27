# MSP Functions

!!! note "Inactive Module"
    These functions were written to process data from the MASSP instrument.
    As far as can be determined from the codebase and available records,
    this pipeline was never placed into operational use -- no MASSP
    instrument deployment has been confirmed to have run through this
    module. They are retained for reference but are not used in OOI data
    production.

## Background

The Multi-species Mass Spectrometer (MASSP) measures dissolved gas
concentrations in hydrothermal vent and cold seep fluids using an
integrated Residual Gas Analyzer (RGA). `msp_functions.py` (originally
implemented by Craig Risien) computes the OOI Level 1 Dissolved Gas
Concentrations (DISSGAS) core data product from the L0 Mass Spectral
Intensities and the sample temperature, both measured by the MASSP
instrument. DISSGAS is composed of the dissolved concentrations (uM) of
seven gases -- methane, ethane, hydrogen, argon, hydrogen sulfide, oxygen,
and carbon dioxide -- in up to four fluid sources: the in situ sample
water, a background water reference, and two shipboard calibration
solutions. For methane, Nafion mode scan data is used; for the other six
gases, Direct mode scan data is used.

A Level 2 Total Dissolved Gas Concentration product (TOTLGAS) is derived
from the L1 DISSGAS hydrogen sulfide and carbon dioxide concentrations,
applying a temperature-, salinity-, and pressure-dependent speciation
correction using the equilibrated water pH.

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00240 | Data Product Specification for Dissolved Gas Concentrations (DISSGAS). (Not released.) |

No DCN was ever assigned to a DPS for the L2 TOTLGAS product; the module's
own header comment carries a placeholder ("1341-00XXX") rather than a real
document number. TOTLGAS is documented here from the code only.

### DISSGAS_L1 -- Dissolved Gas Concentrations

DISSGAS_L1 is the in situ concentration (uM) of a dissolved gas in one of
four fluid sources, computed by `gas_concentration` from a mass-spectral
intensity, a mode-averaged inlet temperature, and the in situ pressure
(derived from `sensor_depth`). The function selects the calibration-table
column(s) bracketing the averaged temperature, evaluates a
pressure-dependent polynomial-plus-exponential fit at each, and
interpolates between two columns when the temperature falls between
calibration points:

$$conc_T = \alpha \times x^2 + \beta \times x + \delta \times e^{\zeta x} + \gamma$$

where $x$ is the deconvolution-corrected intensity from
`deconvolution_correction` and $\alpha$, $\beta$, $\gamma$, $\delta$,
$\zeta$ are pressure-dependent polynomial coefficients drawn from the
calibration table. The final concentration is a linear rescaling of
$conc_T$ using two additional calibration table coefficients.

Twenty-two DISSGAS_L1 products are produced: all seven gases for the
sample and background water, and methane plus carbon dioxide only for
each of the two calibration solutions.

### CALRANG_L1 -- MASSP Calibration Range

CALRANG_L1 is the quality status (dimensionless flag) for each DISSGAS_L1
concentration, computed alongside it by `gas_concentration`. A value of 0
indicates both the intensity and temperature used are within the
calibration range; -1 indicates the intensity was below the calibration
minimum; 1 indicates the intensity was above the calibration maximum but
the temperature was in range; 2 indicates the intensity was in range but
the temperature was above the calibration maximum; 3 indicates both the
intensity and temperature were above range.

### TSTAMP_L1 -- Scan Timestamps

TSTAMP_L1 is the mean port timestamp (seconds) over the Direct- or
Nafion-mode averaging window used to compute the corresponding DISSGAS_L1,
CALRANG_L1, MSINLET_L1, or NAFEFF auxiliary product. A separate TSTAMP_L1
wrapper exists for each DISSGAS_L1/CALRANG_L1 product plus the pH
intensity and Nafion Drier Efficiency products.

### MSINLET_L1 -- Sample Inlet pH Intensity

MSINLET_L1 is the raw pH signal intensity (dimensionless) of the sample,
background, or calibration solution fluid at the time of dissolved gas
measurement, averaged over the last minute of the relevant mode. Four
MSINLET_L1 products exist, one per fluid source.

### NAFEFF -- Nafion Drier Efficiency

NAFEFF is an indicator (percent) of the drying efficiency of the Nafion
drier, computed as 100 times the ratio of the mz 18 (water) signal seen in
Nafion mode to that seen in Direct mode, for the sample water only.

### GASMODE and SMPMODE -- Operating Mode Indicators

GASMODE indicates the MASSP operating mode (-1 for another mode, 0 for
Direct, 1 for Nafion) from four sample valve statuses (`GasModeDetermination`).
SMPMODE indicates the sample measurement mode (2, 1, -1, or -2, depending
on the valve combination) from five external valve statuses
(`SmpModeDetermination`). Both are computed directly from valve status
inputs with no intermediate preprocessing step.

### TOTLGAS_L2 -- Total Dissolved Gas Concentration

TOTLGAS_L2 is the total dissolved concentration (uM) of hydrogen sulfide
or carbon dioxide in the sample or background water, computed from the
corresponding L1 DISSGAS concentration and a speciation correction factor
that accounts for the fraction of the total dissolved species present in
the form measured by the RGA. The correction depends on the equilibrated
water pH (`calc_l2_mswater_smpphval` / `calc_l2_mswater_bkgphval`), the
MSINLET-TEMP inlet temperature, an assumed or measured salinity, and the
depth-derived pressure.

For hydrogen sulfide, the speciation factor is:

$$\beta = 1 + \frac{K_1}{10^{-pH}}$$

For carbon dioxide, the speciation factor is:

$$\alpha = 1 + \frac{K_1}{10^{-pH}} + \frac{K_1 K_2}{\left(10^{-pH}\right)^2}$$

$K_1$ (and $K_2$ for carbon dioxide) are computed from temperature-,
salinity-, and pressure-dependent empirical expressions in the code;
neither the expressions nor their source is named beyond what appears in
the code, since no DPS or code comment attributes them.

### MSWATER_L2 -- Equilibrated Water pH

MSWATER_L2 is the mass-spectrometer-equilibrated water pH (dimensionless),
computed from MSINLET-TEMP and the raw MSINLET pH intensity using a
six-coefficient L2 pH calibration table:

$$pH = \left(A_0 + A_1 \times T + A_2 \times T^2\right) \times
\left(a_2 \times V^2 + a_1 \times V + a_0 + 7\right)$$

where $T$ is MSINLET-TEMP and $V$ is the raw pH intensity. A computed
value outside the range 2-12 is flagged as -9999999.0 rather than
reported.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.msp_functions.gas_concentration

#### Additional Notes
The `concT2_flag` variable used later in the function is only assigned inside the branch that selects two bracketing
calibration-temperature columns; the other two branches (an exact temperature match, or a temperature at or above the
highest calibration point) never set it, so those code paths would raise a `NameError` if reached with a fresh Python
interpreter state. Documented as-is; not corrected, per scope.

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.GasModeDetermination

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.SmpModeDetermination

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_l2_mswater_smpphval

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_l2_mswater_bkgphval

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_l2_totlgas_smph2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_l2_totlgas_smpco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_l2_totlgas_bkgh2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_l2_totlgas_bkgco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
## Helper Functions

::: ion_functions.data.msp_functions.average_mz

#### Additional Notes
The code comment describes taking "the median of the three highest values," but the implementation sorts each scan and
takes the single second-highest value (index -2), not a median of three. Documented from the code; the comment is not
corrected, per scope.

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---
::: ion_functions.data.msp_functions.deconvolution_correction

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---
::: ion_functions.data.msp_functions.rga_status_process

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---
::: ion_functions.data.msp_functions.SamplePreProcess

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---
::: ion_functions.data.msp_functions.BackgroundPreProcess

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---
::: ion_functions.data.msp_functions.Cal1PreProcess

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---
::: ion_functions.data.msp_functions.Cal2PreProcess

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation. |

---
## Wrapper Functions

::: ion_functions.data.msp_functions.calc_dissgas_smpmethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_smpethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_smph2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_smparcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_smph2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_smpo2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_smpco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_bkgmethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_bkgethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_bkgh2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_bkgarcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_bkgh2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_bkgo2con

#### Additional Notes
The module's own header comment table lists this function's OOI product name as DISSGAS-BKGCO2CON and
`calc_dissgas_bkgco2con`'s as DISSGAS-BKGO2CON -- the reverse of the pattern followed by every other gas/water-source
pair. The function body (calibration column range, source mz value) is internally consistent with the name used here;
the header comment appears to be the error and is not propagated into this documentation.

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_bkgco2con

#### Additional Notes
The module's own header comment table lists this function's OOI product name as DISSGAS-BKGCO2CON and
`calc_dissgas_bkgco2con`'s as DISSGAS-BKGO2CON -- the reverse of the pattern followed by every other gas/water-source
pair. The function body (calibration column range, source mz value) is internally consistent with the name used here;
the header comment appears to be the error and is not propagated into this documentation.

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_cal1methcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_cal1co2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_cal2methcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_dissgas_cal2co2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_smpmethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_smpethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_smph2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_smparcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_smph2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_smpo2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_smpco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_bkgmethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_bkgethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_bkgh2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_bkgarcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_bkgh2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_bkgo2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_bkgco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_cal1methcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_cal1co2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_cal2methcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_calrang_cal2co2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_smpmethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_smpethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_smph2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_smparcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_smph2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_smpo2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_smpco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_bkgmethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_bkgethcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_bkgh2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_bkgarcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_bkgh2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_bkgo2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_bkgco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_cal1methcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_cal1co2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_cal2methcon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_cal2co2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_msinlet_smpphint

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_msinlet_smpphint_timestamp

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_msinlet_bkgphint

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_msinlet_bkgphint_timestamp

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_msinlet_cal1phint

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_msinlet_cal1phint_timestamp

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_msinlet_cal2phint

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_msinlet_cal2phint_timestamp

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_smpnafeff

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_smpnafeff_timestamp

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_totlgas_smph2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_totlgas_smpco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_totlgas_bkgh2scon

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
::: ion_functions.data.msp_functions.calc_timestamp_totlgas_bkgco2con

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; added @deprecated (never implemented). |

---
## References

OOI Data Product Specification for Dissolved Gas Concentrations (DISSGAS),
Document Control Number 1341-00240. (Not released.)
