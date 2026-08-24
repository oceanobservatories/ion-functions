# MET Functions

## Background

The Surface Mooring [Bulk Meteorology (METBK) instrument package](https://www.whoi.edu/what-we-do/explore/instruments/instruments-sensors-samplers/air-sea-interaction-meteorology-the-asimet-system/) integrates a suite of sensors
that measure the surface meteorological and near-surface oceanographic variables needed to
compute air-sea fluxes at OOI surface moorings. 

| Data Product | Sensor | Units |
|---|---|---|
| RELHUMI (relative humidity) | Rotronic MP-101A | % |
| TEMPAIR (air temperature) | Rotronic MP-101A | $^\circ$C |
| BARPRES (barometric pressure) | Heise DXD | mbar |
| WINDAVG (mean wind velocity) | Gill Windobserver II | m/s |
| PRECIPM (precipitation) | RM Young 50202 | mm |
| SHRTIRR (downwelling shortwave irradiance) | Kipp & Zonen CMP 21 | W/m$^2$ |
| LONGIRR (downwelling longwave irradiance) | Eppley PIR | W/m$^2$ |
| CONDSRF (sea surface conductivity) | Sea-Bird SBE-37 | S/m |
| TEMPSRF (sea surface temperature) | Sea-Bird SBE-37 | $^\circ$C |

`met_functions.py` covers two distinct processing tiers. A small set of functions compute simple 
L1 conversions and metadata products directly from the sensor readings above recorded at 1-minute 
resolution; the larger remainder of the module implements the TOGA-COARE bulk flux algorithm (warmlayer 
and coolskin corrections feeding the COARE 3.5 iterative solver) used to compute the L2 BULKFLX 
air-sea flux data products. The bulk flux engine and every data product built on it have been deprecated
(see [Deprecated Functions](#deprecated-functions)) in favor of the official TOGA-COARE reference implementation maintained
at [github.com/NOAA-PSL/COARE-algorithm](https://github.com/NOAA-PSL/COARE-algorithm); users requiring bulk air-sea flux products should use that 
implementation directly starting with the 1-minute resolution data defined above.

### Primary Sources

| DCN | Document |
|---|---|
| 1341-00360 | [OOI (2013). Data Product Specification for L1 Bulk Meteorological Data Products.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00360_Data_Product_Spec_BULKMET_OOI.pdf) |
| 1341-00040 | [OOI (2012). Data Product Specification for Salinity.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf) |
| 1341-00370 | OOI. Data Product Specification for L2 BULKFLX Data Products. Document Control Number 1341-00370. (Not released.) |

DPS 1341-00370 (BULKFLX), which specifies the L2 bulk flux algorithm and is cited throughout the
legacy code comments for nearly every product in this module, was never released. Background content 
for CURRENT_DIR, CURRENT_SPD, RELWIND_DIR-AUX, RELWIND_SPD-AUX, NETSIRR_L2, and SPECHUM_L2 below is 
therefore derived from the code and its comments only.

---

### BARPRES_L1 -- Barometric Pressure

BARPRES_L1 is barometric pressure (Pa) at the mooring, computed by converting the L0 Heise DXD
reading from millibars to Pascals (1 mbar = 100 Pa). Per DPS 1341-00360, no other processing is
required for this conversion.

### WINDAVG_L1 -- Mean Wind Velocity

WINDAVG_L1 is the mean wind velocity (m/s), reported as northward and eastward components
relative to true north, representing wind direction toward (not from) the given heading. The
Gill Windobserver II ultrasonic anemometer measures the L0 eastward and northward wind
components in the reference frame of the buoy compass; DPS 1341-00360 specifies rotating these
components by the local magnetic declination to reference true north and east.

The magnetic declination is computed via `generic_functions.magnetic_declination`, which uses
the IGRF-14 model (the DPS specifies the World Magnetic Model; this repository-wide substitution
was made during the WAV and ADCP family sessions and applies here as well). A second correction,
not present in the original DPS, is also applied: a wind-speed-dependent linear or piecewise-
linear correction (offset, slope, and an optional slope change above a threshold) addresses a
known under-reporting of wind speed at higher values. The correction coefficients are supplied
as calibration coefficients and default to values that apply no correction.

### CURRENT_DIR and CURRENT_SPD -- Surface Current Direction and Speed

CURRENT_DIR (degrees, [0, 360)) and CURRENT_SPD (m/s) are metadata products describing the
direction and magnitude of the surface current, computed from the eastward and northward
surface current velocity components (VELPTMN-VLE_L1, VELPTMN-VLN_L1) reported by a co-located
VELPT current meter. Neither product feeds into any other METBK calculation in this module; the
vector difference of wind and current velocity (RELWIND_SPD-AUX, below) is the quantity actually
used by the flux algorithms.

Suspect VELPTMN readings are excluded via a time-vectorized quality flag
(`use_velptmn_with_metbk`), set per-deployment to indicate whether a mooring's VELPT compass is
considered reliable. Earlier deployments of Endurance Array moorings have VELPT current meters mounted
close enough to stainless steel mooring ballast plates that their compass readings were aliased; the 
flag exists to exclude those current data from the calculations. The stainless steel mooring ballast
plates were subsequently replaced with lead plates. Pioneer and Global Array surface moorings do not
have a current meter on the mooring. The closest is located on the midwater platform at 7 or 10 m, 
respectively.

### RELWIND_DIR-AUX and RELWIND_SPD-AUX -- Relative Wind Direction and Speed

RELWIND_DIR-AUX (degrees, [0, 360)) and RELWIND_SPD-AUX (m/s) describe the wind velocity relative
to the surface current: the vector difference of the METBK wind velocity (WINDAVG_L1) and the
VELPT surface current velocity. RELWIND_SPD-AUX is the fundamental wind speed input to every
TOGA-COARE bulk flux calculation in this module.

The two products differ in how they treat unavailable or suspect current data.
RELWIND_DIR-AUX, like CURRENT_DIR and CURRENT_SPD, returns nan for suspect current values.
RELWIND_SPD-AUX instead substitutes 0 for suspect or missing current, so that the calculation
falls back to using the wind speed alone -- matching the DPS-specified behavior for the bulk flux
inputs and reflecting the fact that, at the time this code was written, most OOI surface moorings
had no independent surface current measurements.

### NETSIRR_L2 -- Net Downward Shortwave Radiation

NETSIRR_L2 is the net downward shortwave radiation (W/m$^2$, 0.3 to 3.0 $\mu$m wavelengths),
computed by subtracting the fraction of the measured downward shortwave irradiance (SHRTIRR_L1)
reflected at the sea surface, using a fixed albedo of 0.055 (the reflection coefficient used in
the original TOGA-COARE code). NETSIRR_HOURLY_L2 is the same calculation binned to an hourly
timebase.

### SALSURF_L2 -- Sea Surface Practical Salinity

SALSURF_L2 is the sea surface practical salinity (PSS-78, unitless), computed from the METBK sea
surface conductivity (CONDSRF_L1) and temperature (TEMPSRF_L1) measurements using TEOS-10 (the
`gsw` library), per DPS 1341-00040. The measurement depth serves as a proxy for pressure in the
salinity calculation, consistent with the CTD family's practical salinity computation.

### SPECHUM_L2 -- Air Specific Humidity

SPECHUM_L2 is the air specific humidity (g/kg), computed from air temperature (TEMPAIR_L1),
barometric pressure (BARPRES_L0, not the Pa-converted BARPRES_L1), and relative humidity
(RELHUMI_L1) via a saturation-vapor-pressure formulation. SPECHUM_L2 is not to be confused with
SPHUM2M_L2, the modelled specific humidity at a 2 m reference height computed by the (now
deprecated) bulk flux engine.

Full algorithm derivations, calibration procedures, and source references are listed in the
[References](#references) section.

---

## Core Functions

::: ion_functions.data.met_functions.met_barpres

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_wind_mag_corr

#### History
| Date | Author | Change |
|---|---|---|
| 2025-10-10 | Christopher Wingard | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_current_direction

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-27 | Russell Desiderio | Initial code |
| 2015-07-10 | Russell Desiderio | Added data quality flags (use_velptmn_with_metbk) to argument list |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_current_speed

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-26 | Christopher Wingard | Initial code |
| 2014-08-27 | Russell Desiderio | Added documentation, changed variable names |
| 2015-07-10 | Russell Desiderio | Added data quality flags (use_velptmn_with_metbk) to argument list |
| 2023-08-15 | Samuel Dahlberg | Removed use of numexpr |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_relwind_direction

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-26 | Russell Desiderio | Initial code |
| 2015-07-10 | Russell Desiderio | Set default calling water velocity values and implemented use_velptmn_with_metbk switch |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_relwind_speed

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-26 | Russell Desiderio | Initial code |
| 2015-07-10 | Russell Desiderio | Set default calling water velocity values; added use_velptmn_with_metbk switch and code to vet surface current data |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_netsirr

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-27 | Russell Desiderio | Initial code |
| 2017-02-03 | Russell Desiderio | Added timebase documentation |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_salsurf

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2014-08-26 | Russell Desiderio | Changed variable names |
| 2023-08-15 | Samuel Dahlberg | Replaced incompatible pygsw with GSW library |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_spechum

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-26 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Helper Functions

::: ion_functions.data.met_functions.vet_velptmn_data

#### History
| Date | Author | Change |
|---|---|---|
| 2015-07-10 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Wrapper Functions

::: ion_functions.data.met_functions.met_windavg_mag_corr_east

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2014-08-26 | Russell Desiderio | Added documentation |
| 2025-10-10 | Christopher Wingard | Converted to wrapper function and added wind speed correction factors |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_windavg_mag_corr_north

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-25 | Christopher Wingard | Initial code |
| 2014-08-26 | Russell Desiderio | Added documentation |
| 2025-10-10 | Christopher Wingard | Converted to wrapper function and added wind speed correction factors |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## Deprecated Functions

!!! note "Deprecated"
    The functions in this section implement the TOGA-COARE bulk flux algorithm (warmlayer,
    coolskin, and the L2 BULKFLX data products built on them). They are superseded by the
    official TOGA-COARE reference implementation at
    [github.com/NOAA-PSL/COARE-algorithm](https://github.com/NOAA-PSL/COARE-algorithm) and
    should not be used in new code.

### Core Functions

::: ion_functions.data.met_functions.seasurface_skintemp_correct

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.warmlayer

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2015-07-08 | Russell Desiderio | Added array subscripts to sensor height arrays so these parameters can be either 1-element or time-vectorized arrays |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.coare35vn

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

### Helper Functions

::: ion_functions.data.met_functions.air_density

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-29 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.airtemp_at_refheight

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.calc_rain_rate

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-19 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.gravity

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-26 | Christopher Wingard | Initial code |
| 2014-08-26 | Russell Desiderio | Optimized |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.latent_heat_vaporization_pure_water

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-29 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.net_longwave_up

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.psit_26

#### History
| Date | Author | Change |
|---|---|---|
| 2014-06-26 | Christopher Wingard | Initial code |
| 2014-09-02 | Russell Desiderio | Prevented raising a negative number to a fractional power |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.psiu_26

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-03 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.rain_heat_flux

#### History
| Date | Author | Change |
|---|---|---|
| 2014-10-28 | Russell Desiderio | Initial code (new derivation) |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.sea_spechum

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.spechum_at_refheight

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.water_thermal_expansion

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-29 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.windspeed_at_refheight

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.charnock_wind

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.coolskin_parameters

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.effective_relwind

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.obukhov_for_init

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.obukhov_length_scale

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.roughness_lengths

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.roughness_lengths_for_init

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.scaling_parameters

#### History
| Date | Author | Change |
|---|---|---|
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.condition_data

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-19 | Russell Desiderio | Initial code |
| 2015-07-08 | Russell Desiderio | Added sensor height indices for ztmpwat, zwindsp, ztmpair, zhumair so these parameters can be either 1-element or time-vectorized arrays; conditioned jcool/jwarm to 1-element switches |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.make_hourly_data

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-18 | Russell Desiderio | Initial code |
| 2014-10-22 | Russell Desiderio | Added capability to process timestamps as the sole input argument |
| 2015-07-08 | Russell Desiderio | Deleted sensor height indices from the skip list so these parameters can be either 1-element or time-vectorized arrays |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.warmlayer_time_keys

#### History
| Date | Author | Change |
|---|---|---|
| 2014-10-22 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

### Wrapper Functions

::: ion_functions.data.met_functions.met_buoyfls

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_buoyflx

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_frshflx

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_heatflx

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_latnflx

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_mommflx

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_netlirr

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_rainflx

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2014-10-28 | Russell Desiderio | Incorporated new subroutine for rain heat flux |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_sensflx

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_sphum2m

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_stablty

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_tempa2m

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_tempskn

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_wind10m

#### History
| Date | Author | Change |
|---|---|---|
| 2014-09-01 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_heatflx_minute

#### History
| Date | Author | Change |
|---|---|---|
| 2017-02-13 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_latnflx_minute

#### History
| Date | Author | Change |
|---|---|---|
| 2017-02-13 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_netlirr_minute

#### History
| Date | Author | Change |
|---|---|---|
| 2017-02-13 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_sensflx_minute

#### History
| Date | Author | Change |
|---|---|---|
| 2017-02-13 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_timeflx

#### History
| Date | Author | Change |
|---|---|---|
| 2014-10-22 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_netsirr_hourly

#### History
| Date | Author | Change |
|---|---|---|
| 2017-02-03 | Russell Desiderio | Initial code |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

::: ion_functions.data.met_functions.met_rainrte

#### History
| Date | Author | Change |
|---|---|---|
| 2014-08-27 | Russell Desiderio | Initial code |
| 2014-09-19 | Russell Desiderio | Added front end to convert each-minute data to hourly |
| 2026-07-27 | Christopher Wingard | Converted to NumPy docstring format; updated documentation |

---

## References

Fairall, C.W., E.F. Bradley, D.P. Rogers, J.B. Edson, and G.S. Young (1996). Bulk parameterization
of air-sea fluxes for Tropical Ocean-Global Atmosphere Coupled-Ocean Atmosphere Response
Experiment. *Journal of Geophysical Research*, 101(C2), 3747-3764.

Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B. Edson, and G.S. Young (1996). Cool-skin
and warm-layer effects on sea surface temperature. *Journal of Geophysical Research*, 101(C1),
1295-1308.

Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003). Bulk parameterization
of air-sea fluxes: updates and verification for the COARE algorithm. *Journal of Climate*, 16,
571-590.

Gosnell, R., C.W. Fairall, and P.J. Webster (1995). The sensible heat of rainfall in the tropical
ocean. *Journal of Geophysical Research*, 100(C9), 18437-18442.

NOAA Physical Sciences Laboratory. [COARE-algorithm.](https://github.com/NOAA-PSL/COARE-algorithm)
GitHub repository.

[OOI (2012). Data Product Specification for Salinity. Document Control Number
1341-00040.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf)

[OOI (2013). Data Product Specification for L1 Bulk Meteorological Data Products. Document
Control Number
1341-00360.](https://oceanobservatories.org/wp-content/uploads/2023/09/1341-00360_Data_Product_Spec_BULKMET_OOI.pdf)

OOI. Data Product Specification for L2 BULKFLX Data Products. Document Control Number
1341-00370. (Not released.)
