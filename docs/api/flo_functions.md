# FLO Functions

## Background

The Ocean Observatories Initiative deploys the WET Labs ECO fluorometer
family across its moored, profiling, and mobile platforms to measure
fluorescence and optical backscatter. Two instrument classes are covered
by this module.

| Class | Hardware | Channels | Designator meaning |
|---|---|---|---|
| FLORT | WET Labs ECO Triplet | 3 | Fluorometer, Three Wavelength |
| FLORD | WET Labs ECO Dual | 2 | Fluorometer, Two Wavelength |

Fixed-platform instruments include a wiper to actively limit biofouling.
Mobile-platform instruments (profilers, gliders, AUVs) use only passive
mitigation (copper faceplates).

`flo_functions.py` converts raw L0 count data from these instruments into
three L1 data products — fluorometric chlorophyll-a concentration
(CHLAFLO_L1), fluorometric CDOM concentration (CDOMFLO_L1), and the
volume scattering function (FLUBSCT_L1) — and one L2 product, the total
optical backscatter coefficient (FLUBSCT_L2). All calibration coefficients
are from factory calibration sheets. Within the OOI data system, FLO
functions fall under the Water Column science regime. FLORD instruments
fall under the Chlorophyll a and Inherent Optical Properties categories;
FLORT instruments additionally include the Colored Dissolved Organic Matter
category.

### Primary Sources

| DCN | Document                                                                                                                                                                                                |
|---|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1341-00530 | [OOI (2012). Data Product Specification for Fluorometric Chlorophyll-a Concentration.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00530_Data_Product_Specification_CHLAFLO_OOI.pdf) |
| 1341-00550 | [OOI (2012). Data Product Specification for Fluorometric CDOM Concentration.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00550_Data_Product_Specification_CDOMFLO_OOI.pdf)                   |
| 1341-00540 | [OOI (2014). Data Product Specification for Optical Backscatter (Red Wavelengths).](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00540_Data_Product_SPEC_FLUBSCT_OOI.pdf)             |

---

### CHLAFLO and CDOMFLO — Fluorescence Concentrations

Raw output from the FLORD/FLORT instrument is in counts, ranging from 0
to approximately 4096. The L1 fluorescence concentrations are computed
from a linear scale-and-offset equation using factory-supplied calibration
coefficients:

$$XX = (C_\text{output} - C_\text{dc}) \times \text{SF}$$

where $C_\text{output}$ is the raw count output, $C_\text{dc}$ is the
dark count from the factory calibration sheet, and SF is the scale factor
from the factory calibration sheet. The same equation applies to both
CHLAFLO and CDOMFLO; only the calibration coefficients and units differ.

For **CHLAFLO** (chlorophyll-a), SF has units of $\mu g\ L^{-1}\ \text{count}^{-1}$
and the output is in $\mu g\ L^{-1}$. The instrument measures at
excitation/emission wavelengths of 470/695 nm.

For **CDOMFLO**, SF has units of $\text{ppb}\ \text{count}^{-1}$ and
the output is in parts per billion (ppb). The instrument measures at
excitation/emission wavelengths of 370/460 nm.

The dark count is the instrument's signal output in clean water with
black tape over the detector, measured at the factory. The scale factor
is calculated at the factory by obtaining a consistent output in a
solution of known concentration and subtracting the dark count.

**FLORT/FLORD gliders and AUVs**: Raw data from fluorometers on gliders
and AUVs are processed onboard the vehicle using vendor software and
transmitted already in decimal engineering units; `ion-functions` is not
invoked for those deployments.

---

### FLUBSCT_L1: Total Volume Scattering Function

The L1 optical backscatter product is the total volume scattering function
$\beta(\theta, \lambda)$ ($\text{m}^{-1}\ \text{sr}^{-1}$) at the
instrument's effective (centroid) backscatter angle $\theta$ and measurement
wavelength $\lambda$. It is computed from raw counts using the same
scale-and-offset equation as the fluorescence products:

$$\beta(\theta, \lambda) = (C_\text{output} - C_\text{dc}) \times \text{SF}$$

where SF has units of $(\text{m}^{-1}\ \text{sr}^{-1})\ \text{count}^{-1}$.
Raw counts range from 0 to approximately 4210. The dark count and scale
factor are from the factory calibration sheet.

The centroid angle $\theta$ and chi factor $\chi$ (used below) are 
instrument-dependent:

| Instrument type | OOI classes | $\theta$ | $\chi$ |
|---|---|---|---|
| ECO 3-channel (FLORT D/J/K/M/N/O, FLORD D) | FLORT, FLORD | 124 deg | 1.076 |
| ECO 2-channel (FLORD G/L/M, FLNTU A) | FLORD, FLNTU | 140 deg | 1.096 |

All FLORD and FLORT backscatter channels use a nominal measurement
wavelength of 700 nm. The ECO-BB3 (FLORT series O) is an exception with
three backscatter channels at different visible wavelengths.

### FLUBSCT_L2: Total Optical Backscatter Coefficient

The L2 product is the total optical backscatter coefficient $b_b$
($\text{m}^{-1}$), computed in four steps:

**Step 1** — Compute the seawater volume scattering function $\beta_\text{sw}$
and total seawater scattering coefficient $b_\text{sw}$ using the Zhang
et al. (2009) model with co-located CTD temperature and salinity:

$$[\beta_\text{sw}, b_\text{sw}] = \text{flo_zhang_scatter_coeffs}(\text{degC}, \text{psu}, \theta, \lambda)$$

**Step 2** — Subtract the seawater contribution to obtain the particulate
volume scattering function:

$$\beta_p(\theta, \lambda) = \beta(\theta, \lambda) - \beta_\text{sw}(\theta, \lambda)$$

**Step 3** — Convert to the particulate backscattering coefficient using
the chi factor $\chi$:

$$b_{bp} = \chi \times 2\pi \times \beta_p(\theta, \lambda)$$

The factor of $2\pi$ arises from integration over the polar angle. The
chi factor relates the particulate volume scattering at angle $\theta$ to
the total particulate backscattering coefficient integrated over all
backward angles.

**Step 4** — Add the seawater backscattering contribution. Because
seawater scattering is symmetric in forward and backward directions, the
backward component is half the total:

$$b_b = b_{bp} + \frac{b_\text{sw}}{2}$$

#### Zhang et al. (2009) Seawater Scattering Model

The seawater scattering model (`flo_zhang_scatter_coeffs`) computes the
volume scattering function of pure seawater at a specified angle and the
total seawater scattering coefficient. The model sums contributions from
density fluctuations and concentration fluctuations:

$$\beta_{90,\text{sw}} = \beta_\text{df} + \beta_\text{cf}$$

The density fluctuation term $\beta_\text{df}$ depends on the isothermal
compressibility (Lepple and Millero, 1971), absolute temperature, and the
density derivative of the refractive index (PMH model). The concentration
fluctuation term $\beta_\text{cf}$ depends on the partial derivative of
the refractive index with respect to salinity (Quan and Fry, 1994) and
the water activity (Millero and Leung, 1976). The total scattering
coefficient is obtained by integrating $\beta_{90,\text{sw}}$ over all
angles. The volume scattering function at angle $\theta$ is then:

$$\beta_\text{sw}(\theta) = \beta_{90,\text{sw}} \left(1 + \frac{1 - \delta}{1 + \delta} \cos^2\theta \right)$$

where $\delta = 0.039$ is the default depolarization ratio (Farinato and
Roswell, 1976). Seawater density uses the UNESCO (1981) formulation. The
refractive index of air uses the Ciddor (1996) formula.

Full algorithm derivations, calibration procedures, and source references
are listed in the [References](#references) section.

---

## Core functions

::: ion_functions.data.flo_functions.flo_scale_and_offset

#### Additional Notes
This is the shared core algorithm for all three L1 FLORT/FLORD data
products. The named OOI wrappers `flo_chla`, `flo_cdom`, and `flo_beta`
each call this function with product-specific calibration coefficients.
External users should call this function directly rather than the wrappers
when a single computation is needed for multiple products.

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-30 | Craig Risien | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_bback_total

#### History
| Date | Author | Change |
|---|---|---|
| 2013-07-16 | Christopher Wingard | Initial code |
| 2014-04-23 | Christopher Wingard | Revised to address integration issues and meet DPS intent |
| 2015-10-26 | Russell Desiderio | Removed default argument values; revised documentation |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.flo_functions.flo_zhang_scatter_coeffs

#### History
| Date | Author | Change |
|---|---|---|
| 2013-07-15 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_refractive_index

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-21 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_isotherm_compress

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-21 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_density_seawater

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-21 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

::: ion_functions.data.flo_functions.flo_scat_seawater

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-24 | Christopher Wingard | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

The three functions below are named-product wrappers that apply
`flo_scale_and_offset` under instrument- and product-specific names to
satisfy the OOI single-output data product requirement. External users
should call `flo_scale_and_offset` directly.

::: ion_functions.data.flo_functions.flo_chla

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-30 | Craig Risien | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_cdom

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-30 | Craig Risien | Initial code |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_beta

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-30 | Craig Risien | Initial code |
| 2015-10-23 | Russell Desiderio | Revised documentation |
| 2025-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

Ciddor, P. E. (1996). Refractive index of air: New equations for the
visible and near infrared. Applied Optics, 35(9), 1566-1573.

Farinato, R. S. and Roswell, R. L. (1976). New values of the light
scattering depolarization and anisotropy of water. Journal of Chemical
Physics, 65(2), 593-595.

Lepple, F. K. and Millero, F. J. (1971). The isothermal compressibility
of seawater near one atmosphere. Deep-Sea Research, 18(12), 1233-1254.

Millero, F. J. and Leung, W. H. (1976). The thermodynamics of seawater
at one atmosphere. American Journal of Science, 276, 1035-1077.

Millero, F. J. (1980). The equation of state of seawater. Deep-Sea
Research, 27(3-4), 255-274.

Quan, X. and Fry, E. S. (1994). Empirical equation for the index of
refraction of seawater. Applied Optics, 33(15), 3241-3243.

Sullivan, J. M., Twardowski, M. S., Zaneveld, J. R. V., and Moore, C. C.
(2013). Measuring optical backscattering in water. In A. A. Kokhanovsky
(Ed.), Light Scattering Reviews 7 (pp. 189-224). Springer.

UNESCO (1981). Tenth report of the joint panel on oceanographic tables
and standards. UNESCO Technical Papers in Marine Science No. 38. UNESCO,
Paris.

Zhang, X., Hu, L., and He, M. (2009). Scattering by pure seawater:
Effect of salinity. Optics Express, 17(7), 5698-5710.

[OOI (2012). Data Product Specification for Fluorometric Chlorophyll-a
Concentration. Document Control Number 1341-00530.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00530_Data_Product_Specification_CHLAFLO_OOI.pdf)

[OOI (2012). Data Product Specification for Fluorometric CDOM
Concentration. Document Control Number 1341-00550.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00550_Data_Product_Specification_CDOMFLO_OOI.pdf)

[OOI (2014). Data Product Specification for Optical Backscatter (Red
Wavelengths). Document Control Number 1341-00540.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00540_Data_Product_SPEC_FLUBSCT_OOI.pdf)
