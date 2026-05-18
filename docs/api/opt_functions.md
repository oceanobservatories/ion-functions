# OPT Functions

## Background

The OPT instrument family covers three distinct optical measurement systems
deployed across the OOI arrays: the Sea-Bird Scientific ac-s spectral
absorption and attenuation meter (OPTAA), scalar PAR sensors from multiple
manufacturers (PARAD), and the Sea-Bird Scientific OCR-507 multispectral
downwelling irradiance radiometer (SPKIR). All three produce L1 or L2 data products from raw
counts or voltage, using instrument-specific factory calibration
coefficients.

`opt_functions.py` computes five data products: OPTATTN_L2 and OPTABSN_L2
from the OPTAA instrument, OPTPARW_L1 from three PARAD variants, and
SPECTIR_L1 from the SPKIR instrument. A companion module,
`opt_functions_tscor.py`, provides the wavelength-indexed dictionary of
temperature and salinity correction coefficients used in the OPTAA
processing chain.

| Instrument Class | Hardware | Data Product |
|---|---|---|
| OPTAA | Sea-Bird Scientific ac-s | OPTATTN_L2, OPTABSN_L2 |
| PARAD (RSN Shallow Profiler) | Sea-Bird Scientific PAR LIN 600m | OPTPARW_L1 |
| PARAD (CSPP, PARAD-J) | Sea-Bird Scientific ECO PAR (ECOPARS) | OPTPARW_L1 |
| PARAD (glider, QSP-2100) | Biospherical Instruments QSP-2100 | OPTPARW_L1 |
| PARAD (WFP, QSP-2200) | Biospherical Instruments QSP-2200 | OPTPARW_L1 |
| SPKIR | Sea-Bird Scientific OCR-507 | SPECTIR_L1 |

### Primary Sources

| Document | DCN |
|---|---|
| [OOI (2014). Data Product Specification for Optical Beam Attenuation Coefficient. Document Control Number 1341-00690.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf) | 1341-00690 |
| [OOI (2014). Data Product Specification for Optical Absorption Coefficient. Document Control Number 1341-00700.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf) | 1341-00700 |
| [OOI (2012). Data Product Specification for Photosynthetically Active Radiation (PAR) from Satlantic Instrument on RSN Shallow Profiler. Document Control Number 1341-00720.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00720_Data_Product_SPEC_OPTPARW_Satl_OOI.pdf) | 1341-00720 |
| [OOI (2014). Data Product Specification for Photosynthetically Active Radiation (PAR) from Biospherical Instruments on CGSN Profilers and Mobile Assets. Document Control Number 1341-00721.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00721_Data_Product_SPEC_OPTPARW_Bios_OOI.pdf) | 1341-00721 |
| [OOI (2014). Data Product Specification for Photosynthetically Active Radiation (PAR) from WET Labs Instrument on Coastal Surface Piercing Profiler. Document Control Number 1341-00722.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00722_Data_Product_SPEC_OPTPARW_WETLabs.pdf) | 1341-00722 |
| [OOI (2014). Data Product Specification for Downwelling Spectral Irradiance. Document Control Number 1341-00730.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00730_DATA_PRODUCT_SPEC_SPECTIR_OOI.pdf) | 1341-00730 |

---

### OPTATTN_L2 — Optical Beam Attenuation Coefficient

OPTATTN_L2 is the spectral beam attenuation coefficient for all water
impurities ($c_{pd}(\lambda)$, m$^{-1}$), computed from raw measurements
made by the Sea-Bird Scientific ac-s (OPTAA). It represents the difference between
the total beam attenuation of the water mixture and that of pure seawater:
$c_{pd} = c - c_w$.

The ac-s instrument provides approximately 75 wavelength channels spanning
roughly 400–750 nm in approximately 4 nm steps. It performs concurrent
attenuation and absorption measurements via a dual optical path in a single
instrument, with an optical path length of 0.25 m. The attenuation channel
uses a non-reflecting tube; the absorption channel uses a highly polished
reflecting tube.

OPTATTN_L2 is an L2 product because it requires co-located and synchronized
measurements of water temperature (TEMPWAT_L1) and practical salinity
(PRACSAL_L2) from a CTD in addition to the raw ac-s signals.

#### Calibration Coefficients

Each OPTAA instrument ships with a device file (`*.dev`) containing four
types of instrument-specific calibration data:

- `coff` — pure water offsets for the attenuation ('c') channel (m$^{-1}$),
  one per wavelength
- `tcal` — the factory calibration reference temperature (pure water,
  $^\circ$C) at which the offsets were determined
- `tbins` — internal temperature calibration bin values ($^\circ$C)
- `tc_arr` — internal temperature calibration correction coefficients for
  the 'c' channel (m$^{-1}$), a 2-D array indexed by wavelength and
  temperature bin

In addition, a community-standard file of wavelength-dependent temperature
and salinity correction coefficients ($\Psi_t$, $\Psi_{sc}$) from the
Sea-Bird Scientific ac Meter Protocol Document (Revision P) is tabulated in
`opt_functions_tscor.py`.

#### OPTATTN_L2 Algorithm

The computation proceeds in four steps (DPS 1341-00690, Section 4.3):

**Step 1 — Compute internal instrument temperature** from the raw thermistor
count (OPTTEMP_L0):

$$T_{intr} = \frac{1}{a + b \ln R + c (\ln R)^3} - 273.15$$

where $R$ is the thermistor resistance derived from the raw count, and
$a$, $b$, $c$ are fixed constants.

**Step 2 — Convert raw counts to uncorrected beam attenuation** and apply
the internal temperature correction:

$$c_{pd}(\lambda) = \left[c_{off}(\lambda) - \frac{1}{r}
\ln\!\left(\frac{C_{sig}}{C_{ref}}\right)\right] - \Delta T(\lambda)$$

where $r = 0.25$ m is the path length, $C_{sig}$ and $C_{ref}$ are the
raw signal and reference counts (OPTCSIG_L0 and OPTCREF_L0), $c_{off}$ is
the pure water offset from the device file, and $\Delta T(\lambda)$ is the
linearly interpolated internal temperature correction:

$$\Delta T(\lambda) = \Delta T_1(\lambda) +
\frac{T_{intr} - T_1}{T_2 - T_1}
\times \left[\Delta T_2(\lambda) - \Delta T_1(\lambda)\right]$$

where $T_1$ and $T_2$ are the bracketing temperature bin values and
$\Delta T_1$, $\Delta T_2$ are the corresponding correction coefficients
from the device file.

**Step 3 — Apply the water temperature and salinity correction:**

$$c_{pd;ts}(\lambda) = c_{pd}(\lambda)
- \Psi_t(\lambda) \times (t - t_r)
- \Psi_{sc}(\lambda) \times s$$

where $t$ is the in situ temperature (TEMPWAT_L1), $t_r$ is the factory
calibration reference temperature, $s$ is the in situ practical salinity
(PRACSAL_L2), and $\Psi_t$ and $\Psi_{sc}$ are the wavelength-dependent
temperature and salinity correction constants from the Sea-Bird Scientific
ac Meter Protocol Document (Revision P). The resulting $c_{pd;ts}(\lambda)$
is OPTATTN_L2.

---

### OPTABSN_L2 — Optical Absorption Coefficient

OPTABSN_L2 is the spectral optical absorption coefficient for all water
impurities ($a_{pd}(\lambda)$, m$^{-1}$), computed from raw measurements
made by the Sea-Bird Scientific ac-s (OPTAA). It represents the difference between
the total absorption of the water mixture and that of pure seawater:
$a_{pd} = a - a_w$.

OPTABSN_L2 is an L2 product requiring co-located CTD temperature
(TEMPWAT_L1) and salinity (PRACSAL_L2), and the simultaneously computed
OPTATTN_L2 for the scattering correction.

#### Calibration Coefficients

The OPTABSN processing uses the same device file as OPTATTN, but draws on
the absorption ('a') channel coefficients:

- `aoff` — pure water offsets for the absorption ('a') channel (m$^{-1}$)
- `ta_arr` — internal temperature correction coefficients for the 'a'
  channel (m$^{-1}$), indexed by wavelength and temperature bin
- `tcal`, `tbins` — shared with the attenuation channel (above)

#### OPTABSN_L2 Algorithm

The computation follows four steps (DPS 1341-00700, Section 4.3), using
the same raw-to-scientific conversion and internal temperature correction as
OPTATTN_L2 (Steps 1 and 2 above), but applied to the absorption channel
(OPTASIG_L0 and OPTAREF_L0 in place of OPTCSIG_L0 and OPTCREF_L0, and
`aoff` / `ta_arr` in place of `coff` / `tc_arr`).

**Step 3 — Apply the water temperature and salinity correction:**

$$a_{pd;ts}(\lambda) = a_{pd}(\lambda)
- \Psi_t(\lambda) \times (t - t_r)
- \Psi_{sa}(\lambda) \times s$$

where $\Psi_{sa}$ is the absorption-channel salinity correction constant
(distinct from the attenuation-channel $\Psi_{sc}$ used in OPTATTN).

**Step 4 — Apply the scattering correction.** The ac-s absorption tube does
not detect light scattered into the backward direction; the measured signal
therefore overestimates absorption. The correction parameterizes the
fraction of undetected scatter $\tilde{\beta}$ using a reference wavelength
(default 715 nm) where absorption by water impurities is assumed negligible:

$$\tilde{\beta} = \frac{a_{pd;ts}(\lambda_{ref})}
{c_{pd;ts}(\lambda_{ref}) - a_{pd;ts}(\lambda_{ref})}$$

$$a_{pd;ts;s}(\lambda) = a_{pd;ts}(\lambda)
- \tilde{\beta} \times \left[c_{pd;ts}(\lambda) - a_{pd;ts}(\lambda)\right]$$

where $c_{pd;ts}(\lambda)$ values are interpolated from the attenuation
channel wavelength grid onto the absorption channel wavelength grid before
the subtraction. The correction is suppressed when
$c_{pd;ts}(\lambda_{ref}) - a_{pd;ts}(\lambda_{ref}) \leq 0.02$ or when
$a_{pd;ts}(\lambda_{ref}) \leq 0$. The resulting
$a_{pd;ts;s}(\lambda)$ is OPTABSN_L2.

---

### OPTPARW_L1 — Photosynthetically Active Radiation

OPTPARW_L1 is the photosynthetically active radiation ($\mu$mol photons
m$^{-2}$ s$^{-1}$) within the spectral range 400–700 nm. Three instrument
variants are deployed across the OOI arrays, each using a different
calibration equation.

**Satlantic PAR LIN 600m** (PARAD, RSN Shallow Profiler): applies a linear
calibration to 32-bit ADC counts (DPS 1341-00720):

$$OPTPARW = I_m \times a_1 \times (x - a_0)$$

where $x$ is the L0 count output, $a_0$ is the voltage offset (counts),
$a_1$ is the scaling factor ($\mu$mol photons m$^{-2}$ s$^{-1}$ count$^{-1}$),
and $I_m$ is the immersion coefficient.

**Sea-Bird Scientific ECO PAR / ECOPARS** (PARAD-J, CSPP): applies an exponential
calibration to 14-bit ADC counts (DPS 1341-00722):

$$OPTPARW = I_m \times 10^{(x - a_0) / a_1}$$

where $x$, $a_0$, $a_1$, and $I_m$ carry the same definitions as above.

**Biospherical QSP-2100** (PARAD, glider / mobile assets): applies a
linear calibration to a voltage signal (DPS 1341-00721):

$$OPTPARW = \frac{output - dark\_offset}{scale\_wet}$$

where $output$ is the sensor reading in volts, $dark\_offset$ is the dark
reading in volts, and $scale\_wet$ is the wet calibration scale factor
(V per $\mu$mol photons m$^{-2}$ s$^{-1}$).

**Biospherical QSP-2200** (PARAD, WFP): uses the same equation as the
QSP-2100, but the input is in millivolts and $scale\_wet$ is in
V per quanta cm$^{-2}$ s$^{-1}$. The function converts millivolt inputs to
volts and the scale factor to SI units
(1 $\mu$mol photons m$^{-2}$ s$^{-1}$ = $6.02 \times 10^{13}$ quanta
cm$^{-2}$ s$^{-1}$) before applying the linear equation.

Calibration coefficients for all four variants come from factory calibration
sheets supplied by the manufacturer.

---

### SPECTIR_L1 — Downwelling Spectral Irradiance

SPECTIR_L1 is the downwelling spectral "vector" irradiance ($E_d(\lambda)$,
$\mu$W cm$^{-2}$ nm$^{-1}$) measured by the Sea-Bird Scientific OCR-507
multispectral radiometer (SPKIR). The OCR-507 provides seven user-defined wavelength
channels, typically spanning 400–700 nm in bandwidths of nominally 10 or
20 nm, with center wavelengths specified at instrument purchase.

SPECTIR_L1 is an L1 product computed from 32-bit ADC counts using three
instrument-specific calibration coefficients per channel, supplied by the
manufacturer:

$$E_d(\lambda) = I_m \times a_1 \times (x - a_0)$$

where $x$ is the raw ADC count (SPECTIR_L0), $a_0$ is the dark offset
(counts), $a_1$ is the scale factor ($\mu$W cm$^{-2}$ nm$^{-1}$ count$^{-1}$),
and $I_m$ is the immersion coefficient (set to 1.0 when the sensor is not
immersed). The calibration equation type is designated `OPTICS2` in the
Sea-Bird Scientific calibration file format (DPS 1341-00730, Section 4.2).

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core Functions

::: ion_functions.data.opt_functions.opt_internal_temp

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-25 | Christopher Wingard | Initial implementation |
| 2014-03-07 | Russell Desiderio | Reduced calls to np.log |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_pd_calc

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-25 | Christopher Wingard | Initial implementation |
| 2014-02-19 | Russell Desiderio | Expanded documentation |
| 2015-04-21 | Russell Desiderio | Added diagnostics to ValueError exceptions |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_tempsal_corr

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-25 | Christopher Wingard | Initial implementation |
| 2014-02-19 | Russell Desiderio | Expanded documentation; removed incorrect vector length requirement |
| 2014-03-21 | Russell Desiderio | Added dictionary comprehension to vectorize correction calculation |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_scatter_corr

#### Additional Notes
The scattering correction is suppressed when
$c_{pd;ts}(\lambda_{ref}) - a_{pd;ts}(\lambda_{ref}) \leq 0.02$ or when
$a_{pd;ts}(\lambda_{ref}) \leq 0$, to avoid numerical instability from
near-zero denominators and to avoid applying corrections when scattering
is too small to be meaningfully resolved. The threshold value of 0.02 was
chosen empirically. This behavior is not specified in DPS 1341-00700 but
is documented in the source code (2015-12-08 revision).

The code comments also note that the accepted deployment protocol (separate
intake tubes for the 'a' and 'c' flow paths) means the two paths do not
sample exactly the same water, introducing uncertainty in the
$c_{pd;ts}(\lambda_{ref}) - a_{pd;ts}(\lambda_{ref})$ quantity used to
compute the scatter ratio.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-25 | Christopher Wingard | Initial implementation |
| 2014-02-19 | Russell Desiderio | Trapped out potential problems in scatter ratio calculation |
| 2015-12-08 | Russell Desiderio | Made the scatter ratio calculation more robust; added minimum threshold |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_par_satlantic

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-31 | Craig Risien | Initial code |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_par_wetlabs

#### History
| Date | Author | Change |
|---|---|---|
| 2014-12-10 | Craig Risien | Initial code |
| 2015-04-09 | Russell Desiderio | Fixed bug so that the function runs correctly on time-vectorized arguments |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_par_biospherical_mobile

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-31 | Craig Risien | Initial code |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_par_biospherical_wfp

#### History
| Date | Author | Change |
|---|---|---|
| 2014-03-07 | Craig Risien | Initial code |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_ocr507_irradiance

#### History
| Date | Author | Change |
|---|---|---|
| 2014-03-14 | Russell Desiderio | Initial code |
| 2014-03-25 | Russell Desiderio | Changed code to require data inputs to be arrays with 7 columns |
| 2015-04-09 | Russell Desiderio | Revised for time-vectorized calibration coefficient arrays |
| 2015-04-21 | Russell Desiderio | Revised for 2-D array input and output |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

Two helper functions support optional auxiliary sensors that may be fitted
to some OPTAA units. Neither is used in the computation of OPTATTN_L2 or
OPTABSN_L2.

::: ion_functions.data.opt_functions.opt_pressure

#### Additional Notes
This function is not used in the computation of OPTATTN_L2 or OPTABSN_L2.
It supports the optional auxiliary pressure sensor available on some OPTAA
units.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-25 | Christopher Wingard | Initial implementation |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_external_temp

#### Additional Notes
This function is not used in the computation of OPTATTN_L2 or OPTABSN_L2.
It supports the optional auxiliary temperature sensor available on some
OPTAA units.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-25 | Christopher Wingard | Initial implementation |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

::: ion_functions.data.opt_functions.opt_beam_attenuation

#### Additional Notes
All input arrays are assumed to be vectorized over data packets such that
the first dimension iterates over packet number. The function loops over
packets and calls `opt_internal_temp`, `opt_pd_calc`, and
`opt_tempsal_corr` for each.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-25 | Christopher Wingard | Initial implementation |
| 2014-03-06 | Russell Desiderio | Reset dimensions of arguments; implemented loop for vectorized input |
| 2014-03-07 | Russell Desiderio | Added documentation |
| 2014-03-21 | Russell Desiderio | Added wavelength rounding to match tscor dictionary keys |
| 2014-05-29 | Russell Desiderio | Added handling for wavelengths outside the T/S correction table range |
| 2015-04-17 | Russell Desiderio | Use np.nan instead of fill value |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.opt_functions.opt_optical_absorption

#### Additional Notes
All input arrays are assumed to be vectorized over data packets such that
the first dimension iterates over packet number. The function loops over
packets and calls `opt_internal_temp`, `opt_pd_calc`, `opt_tempsal_corr`,
and `opt_scatter_corr` for each.

The `rwlngth` parameter (default 715 nm) sets the reference wavelength for
the proportional scattering correction. The default of 715 nm is specified
in DPS 1341-00700; the parameter is exposed to allow overriding when the
instrument's channel closest to 715 nm differs significantly.

#### History
| Date | Author | Change |
|---|---|---|
| 2013-04-25 | Christopher Wingard | Initial implementation |
| 2014-02-19 | Russell Desiderio | Added rwlngth to argument list |
| 2014-03-06 | Russell Desiderio | Reset dimensions of arguments; implemented loop for vectorized input |
| 2014-03-07 | Russell Desiderio | Added documentation |
| 2014-03-21 | Russell Desiderio | Added wavelength rounding to match tscor dictionary keys |
| 2014-05-29 | Russell Desiderio | Added handling for wavelengths outside the T/S correction table range |
| 2015-04-17 | Russell Desiderio | Use np.nan instead of fill value |
| 2025-05-18 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

Ackleson, S. G. and O'Donnell, J. (2011). Small-scale variability in
suspended matter associated with the Connecticut River plume front.
J. Geophys. Res., 116(C10013), doi:10.1029/2011JC007053.

Boss, E. et al. (2001). Spectral particulate attenuation and particle size
distribution in the bottom boundary layer of a continental shelf.
J. Geophys. Res., 106(C5), 9509–9516, doi:10.1029/2000JC900077.

Kitchen, J. C., Zaneveld, J. R. V., and Pak, H. (1982). Effects of particle
size distribution and chlorophyll content on beam attenuation spectra.
Appl. Opt., 21, 3913–3918.

Mobley, C. D. (1994). Light and Water. San Diego, CA: Academic Press.

Mueller, J. L. et al. (2003). Ocean Optics Protocols for Satellite Ocean
Color Sensor Validation, Revision 4, Volume IV: Inherent Optical Properties.
NASA/TM-2003-211621, Vol. IV.

[OOI (2014). Data Product Specification for Optical Beam Attenuation Coefficient. Document Control Number 1341-00690.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00690_Data_Product_SPEC_OPTATTN_OOI.pdf)

[OOI (2014). Data Product Specification for Optical Absorption Coefficient. Document Control Number 1341-00700.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00700_Data_Product_SPEC_OPTABSN_OOI.pdf)

[OOI (2012). Data Product Specification for Photosynthetically Active Radiation (PAR) from Satlantic Instrument on RSN Shallow Profiler. Document Control Number 1341-00720.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00720_Data_Product_SPEC_OPTPARW_Satl_OOI.pdf)

[OOI (2014). Data Product Specification for Photosynthetically Active Radiation (PAR) from Biospherical Instruments on CGSN Profilers and Mobile Assets. Document Control Number 1341-00721.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00721_Data_Product_SPEC_OPTPARW_Bios_OOI.pdf)

[OOI (2014). Data Product Specification for Photosynthetically Active Radiation (PAR) from WET Labs Instrument on Coastal Surface Piercing Profiler. Document Control Number 1341-00722.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00722_Data_Product_SPEC_OPTPARW_WETLabs.pdf)

[OOI (2014). Data Product Specification for Downwelling Spectral Irradiance. Document Control Number 1341-00730.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00730_DATA_PRODUCT_SPEC_SPECTIR_OOI.pdf)

Pegau, W. S., Gray, D., and Zaneveld, J. R. V. (1997). Absorption and
attenuation of visible and near-infrared light in water: dependence on
temperature and salinity. Appl. Opt., 36, 6035–6046.

Sea-Bird Scientific. *Spectral Absorption and Attenuation Meter ac-s
User's Guide, Revision J.* Philomath, OR: Sea-Bird Scientific.

Sea-Bird Scientific. *ac Meter Protocol Document, Revision P.* Philomath,
OR: Sea-Bird Scientific.

Sea-Bird Scientific. *ECO PAR Sensor User Manual, Edition 5.* Philomath,
OR: Sea-Bird Scientific.

Sea-Bird Scientific. *OCR-507 Operation Manual, SAT-DN-0027, Rev E.*
Philomath, OR: Sea-Bird Scientific.

Sea-Bird Scientific. *PAR Sensor Operation Manual, SAT-DN-00462, Rev B.*
Philomath, OR: Sea-Bird Scientific.

Sullivan, J. M., Twardowski, M. S., Zaneveld, J. R. V., Moore, C. M.,
Barnard, A. H., Donaghay, P. L., and Rhoades, B. (2006). The hyperspectral
temperature and salt dependencies of absorption by water and heavy water in
the 400–750 nm spectral range. Appl. Opt., 45, 5294–5309.
