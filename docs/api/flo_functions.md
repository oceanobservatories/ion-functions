# FLO Functions

## Background

The Ocean Observatories Initiative deploys the Sea-Bird Scientific (formerly WET Labs) ECO fluorometer
family across its moored, profiling, and mobile platforms to measure
fluorescence and optical backscatter. Two instrument classes are covered
by this module.

| Class | Hardware | Platform | Designator Meaning |
|---|---|---|---|
| FLORT | Sea-Bird Scientific (formerly WET Labs) ECO Triplet | Moored (fixed depth), cabled shallow profilers, profilers, gliders, and AUVs | Fluorometer, Three Wavelength |
| FLORD | Sea-Bird Scientific (formerly WET Labs) ECO Dual | Moored (fixed depth), cabled 200 m platforms, wire-following profilers, and gliders | Fluorometer, Two Wavelength |

Fixed-platform instruments include a wiper to actively limit biofouling.
Mobile-platform instruments (profilers, gliders, AUVs) use only passive
mitigation (copper faceplates).

`flo_functions.py` converts raw L0 count data from these instruments into
three L1 data products -- fluorometric chlorophyll-a concentration
(CHLAFLO_L1), fluorometric CDOM concentration (CDOMFLO_L1), and the
volume scattering function (FLUBSCT_L1) -- and one L2 product, the total
optical backscatter coefficient (FLUBSCT_L2).

Within the OOI data system, FLO data products fall under the Water Column [sampling
regime](https://oceanobservatories.org/data-products/) (see the [complete list of OOI data
products](https://oceanobservatories.org/data-products/data-products-complete-list/)).

Calibration coefficients are from the factory calibration sheets supplied with individual instruments, except where
noted in the Calibration Coefficients tables below.

This page documents the algorithms used to compute the FLO data products. For general guidance on working with OOI
biogeochemical sensor data, see the [OOI Biogeochemical Sensor Data Best Practices and User
Guide](https://doi.org/10.25607/OBP-1865.2); for locating calibration coefficients and raw inputs, see [Calibration
Coefficients](../calibrations.md) and [Using These Functions with OOI Data](../using-ooi-data.md).

### Primary Sources

| DCN | Document                                                                                                                                                                                                |
|---|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1341-00530 | [OOI (2012). Data Product Specification for Fluorometric Chlorophyll-a Concentration.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00530_Data_Product_Specification_CHLAFLO_OOI.pdf) |
| 1341-00550 | [OOI (2012). Data Product Specification for Fluorometric CDOM Concentration.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00550_Data_Product_Specification_CDOMFLO_OOI.pdf)                   |
| 1341-00540 | [OOI (2014). Data Product Specification for Optical Backscatter (Red Wavelengths).](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00540_Data_Product_SPEC_FLUBSCT_OOI.pdf)             |

---

### CHLAFLO and CDOMFLO -- Fluorescence Concentrations

Raw output from the FLORD/FLORT instrument is in counts, ranging from 0 to approximately 4096.

#### Calibration Coefficients

| Function Argument | Calibration Coefficient | Source |
|---|---|---|
| `counts_dark` | `CC_dark_counts_chlorophyll_a` | Factory calibration sheet (CHLAFLO) |
| `scale_factor` | `CC_scale_factor_chlorophyll_a` | Factory calibration sheet (CHLAFLO) |
| `counts_dark` | `CC_dark_counts_cdom` | Factory calibration sheet (CDOMFLO) |
| `scale_factor` | `CC_scale_factor_cdom` | Factory calibration sheet (CDOMFLO) |

OOI calibration coefficient names are explained in [Using These Functions with OOI Data](../using-ooi-data.md).

#### CHLAFLO and CDOMFLO Algorithm

The L1 fluorescence concentrations are computed
from a linear scale-and-offset equation using factory-supplied calibration
coefficients:

$$XX = (C_\text{output} - C_\text{dc}) \times \text{SF}$$

where $C_\text{output}$ is the raw count output, $C_\text{dc}$ is the dark count (`counts_dark`) and SF is the scale
factor (`scale_factor`), both from the factory calibration sheet. The same equation applies to both
CHLAFLO and CDOMFLO; only the calibration coefficients and units differ.

Chlorophyll-a absorbs light in the visible spectrum and re-emits a portion as fluorescence; the intensity of that
fluorescence is used as a proxy for chlorophyll-a concentration, and therefore for phytoplankton biomass. CDOM is the
fraction of dissolved organic matter that absorbs ultraviolet light and fluoresces in the visible blue; it is largely
made up of tannins and lignins from decaying plant material and gives some water masses their tea-like color.

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

**FLORT/FLORD gliders and AUVs**: Fluorescence (CHLAFLO, CDOMFLO) from gliders and AUVs is computed onboard the
vehicle using vendor software and transmitted already in engineering units; `ion-functions` is not invoked for those
products. Glider backscatter is the exception (see [FLUBSCT_L1](#flubsct_l1-total-volume-scattering-function) and
[FLUBSCT_L2](#flubsct_l2-total-optical-backscatter-coefficient)).

---

### FLUBSCT_L1: Total Volume Scattering Function

The L1 optical backscatter product is the total volume scattering function $\beta(\theta, \lambda)$
($\text{m}^{-1}\ \text{sr}^{-1}$) at the instrument's effective (centroid) backscatter angle $\theta$ and measurement
wavelength $\lambda$.

#### Calibration Coefficients

| Function Argument | Calibration Coefficient | Source |
|---|---|---|
| `counts_dark` | `CC_dark_counts_volume_scatter` | Factory calibration sheet |
| `scale_factor` | `CC_scale_factor_volume_scatter` | Factory calibration sheet |
| `counts_dark` | `CC_1_`, `CC_2_`, `CC_3_dark_counts_volume_scatter` | Factory calibration sheet (ECO-BB3, 660/532/470 nm channels) |
| `scale_factor` | `CC_1_`, `CC_2_`, `CC_3_scale_factor_volume_scatter` | Factory calibration sheet (ECO-BB3, 660/532/470 nm channels) |

OOI calibration coefficient names are explained in [Using These Functions with OOI Data](../using-ooi-data.md).

#### FLUBSCT_L1 Algorithm

The product is computed from raw counts using the same
scale-and-offset equation as the fluorescence products:

$$\beta(\theta, \lambda) = (C_\text{output} - C_\text{dc}) \times \text{SF}$$

where SF has units of $(\text{m}^{-1}\ \text{sr}^{-1})\ \text{count}^{-1}$.
Raw counts range from 0 to approximately 4210. The dark count (`counts_dark`) and scale factor (`scale_factor`) are
from the factory calibration sheet.

On gliders, the FLORT-M and FLORD-M volume scattering function is computed onboard the vehicle. The ECO-BB3 (FLORT
series O) on coastal gliders reports raw counts for each of its three channels, which `flo_beta` converts using the
`CC_1_`, `CC_2_`, and `CC_3_` coefficients for the 660, 532, and 470 nm channels.

The centroid angle $\theta$ and chi factor $\chi$ (used below) are 
instrument-dependent:

| Instrument type | OOI classes | $\theta$ | $\chi$ |
|---|---|---|---|
| ECO 3-channel (FLORT D/J/K/M/N/O, FLORD D) | FLORT, FLORD | 124 deg | 1.076 |
| ECO 2-channel (FLORD G/L/M, FLNTU A) | FLORD, FLNTU | 140 deg | 1.096 |

All FLORD and FLORT backscatter channels use a nominal measurement
wavelength of 700 nm. The ECO-BB3 (FLORT series O) is an exception with
three backscatter channels at different visible wavelengths.

The 3-channel values ($\theta$ = 124 deg, $\chi$ = 1.076) are from Sullivan et al. (2013), Table 6.2b ("ECO-BB");
older WET Labs documentation lists the centroid angle as 117 deg. The 2-channel values are from M. Twardowski
(personal communication). The chi factor is thought to be at most weakly dependent on wavelength, so these values
apply at any visible wavelength. The ECO-BB3 was originally classified as OPTAA series M and is now FLORT series O.

The chi factor is a scaling factor, not an angular resolution. In the OOI calibration records it is stored as
`CC_angular_resolution`, a historical mislabel.

In this documentation "total" can mean seawater plus particulate scattering, forward plus backward scattering, or
backscatter integrated over all backward angles, depending on context.

### FLUBSCT_L2: Total Optical Backscatter Coefficient

The L2 product is the total optical backscatter coefficient $b_b$
($\text{m}^{-1}$), computed in four steps. It is computed by `ion-functions` for all FLORT and FLORD backscatter,
including gliders, using the co-located (or glider) CTD temperature and practical salinity.

#### Calibration Coefficients

| Function Argument | Calibration Coefficient | Source |
|---|---|---|
| `theta` | `CC_scattering_angle` | OOI processing setting (centroid angle; see the table in FLUBSCT_L1) |
| `xfactor` | `CC_angular_resolution` | OOI processing setting (chi factor; the name is a historical mislabel) |
| `wlngth` | `CC_measurement_wavelength` | Instrument specification (nominal 700 nm) |
| `wlngth` | `CC_1_`, `CC_2_`, `CC_3_measurement_wavelength` | Instrument specification (ECO-BB3 channel wavelengths) |
| `delta` | `CC_depolarization_ratio` | OOI processing setting (0.039; Farinato and Roswell, 1976) |

OOI calibration coefficient names are explained in [Using These Functions with OOI Data](../using-ooi-data.md).

#### FLUBSCT_L2 Algorithm

**Step 1** -- Compute the seawater volume scattering function $\beta_\text{sw}$
and total seawater scattering coefficient $b_\text{sw}$ using the Zhang
et al. (2009) model with co-located CTD temperature and salinity:

$$[\beta_\text{sw}, b_\text{sw}] = \text{flo_zhang_scatter_coeffs}(\text{degC}, \text{psu}, \theta, \lambda)$$

where $\theta$ (`theta`) is the instrument centroid angle and $\lambda$ (`wlngth`) the measurement wavelength, both
OOI processing settings from the table in [FLUBSCT_L1](#flubsct_l1-total-volume-scattering-function).

**Step 2** -- Subtract the seawater contribution to obtain the particulate
volume scattering function:

$$\beta_p(\theta, \lambda) = \beta(\theta, \lambda) - \beta_\text{sw}(\theta, \lambda)$$

**Step 3** -- Convert to the particulate backscattering coefficient using the chi factor $\chi$ (`xfactor`, an OOI
processing setting from the table in [FLUBSCT_L1](#flubsct_l1-total-volume-scattering-function)):

$$b_{bp} = \chi \times 2\pi \times \beta_p(\theta, \lambda)$$

The factor of $2\pi$ arises from integration over the polar angle. The
chi factor relates the particulate volume scattering at angle $\theta$ to
the total particulate backscattering coefficient integrated over all
backward angles.

**Step 4** -- Add the seawater backscattering contribution. Because
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

where $\delta = 0.039$ (`delta`) is the depolarization ratio (Farinato and Roswell, 1976), an OOI processing
setting. Seawater density uses the UNESCO (1981) formulation. The
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
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_bback_total

#### History
| Date | Author | Change |
|---|---|---|
| 2013-07-16 | Christopher Wingard | Initial code |
| 2014-04-23 | Christopher Wingard | Revised to address integration issues and meet DPS intent |
| 2015-10-26 | Russell Desiderio | Removed default argument values; revised documentation |
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.flo_functions.flo_zhang_scatter_coeffs

#### History
| Date | Author | Change |
|---|---|---|
| unknown | Xiaodong Zhang (University of North Dakota) | Original Matlab code |
| 2013-07-15 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_refractive_index

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-21 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_isotherm_compress

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-21 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_density_seawater

#### History
| Date | Author | Change |
|---|---|---|
| 2014-02-21 | Christopher Wingard | Initial code |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

::: ion_functions.data.flo_functions.flo_scat_seawater

#### History
| Date | Author | Change |
|---|---|---|
| 2014-04-24 | Christopher Wingard | Initial code |
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

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
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_cdom

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-30 | Craig Risien | Initial code |
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.flo_functions.flo_beta

#### History
| Date | Author | Change |
|---|---|---|
| 2014-01-30 | Craig Risien | Initial code |
| 2015-10-23 | Russell Desiderio | Revised documentation |
| 2026-04-17 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

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
