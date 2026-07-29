#!/usr/bin/env python
"""
@package ion_functions.data.met_functions
@file ion_functions/data/met_functions.py
@author Russell Desiderio
@brief Module containing functions for the met family of instruments
"""

import numpy as np
import gsw

from ion_functions.data.generic_functions import magnetic_declination, magnetic_correction
from ion_functions import deprecated

### July 2015
""" use_velptmn_with_metbk """
# Until VELPT current meters are
#     (1) co-located with the METBK instrumentation to measure surface currents and
#     (2) without ferrous interference aliasing the VELPT compass measurements,
# VELPT current measurements should not be used to calculate relative wind speed (with
# respect to water), which is the fundamental windspeed variable used in the METBK
# calculations. Almost all of the METBK L2 data products require the relative wind speed
# as a calling argument.
#
# The DPA to calculate relative wind speed over water is currently set to return actual
# wind speed (as if the current velocities were measured to be 0) for all cases of input
# current velocity values (absent, present, nan).
#
# A subset of the Endurance moorings are the only ones that have VELPT instruments
# mounted to measure surface currents. However, their compass readings are inaccurate due
# to the mounted instruments' proximity to iron ballast in the mooring.
#
# It is anticipated that these moorings will be modified to eliminate this magnetic
# interference. To use the velptmn measurements for METBK calculations on these moorings,
# a 5th variable, use_velptmn_with_metbk, has been added to the argument list of the
# function met_relwind_speed. Implementation of the use_velptmn_with_metbk variable will
# require that it be treated as a "platform/instrument instance specific metadata parameter
# changeable in time". It has been coded to accept time-vectorized input.
#
### Further documentation is contained in the Notes to function met_relwind_speed in this module.

"""
    METBK SENSOR HEIGHTS

    Note that the sensor heights may depend on the type of mooring:
"""
# these 4 sensor height variables were time-vectorized in the July 2015 revision.
#     zwindsp = height of the wind measurement [m]
#     ztmpair = height of air temperature measurement [m]
#     zhumair = height of air humidity measurement [m]
#     ztmpwat = depth of bulk sea surface water measurements [m]

#     zvelptm = depth of surface current measurement [m]:
#         this parameter is specified as metadata in the DPS;
#         however, it is not used in the code.

#     zinvpbl = planetary boundary layer/inversion height: this is
#               set to a default value of 600m in the code, as is
#               used in the DPS. this variable was written to accept
#               time-vectorized input in the initial METBK code.

"""
    Set algorithm switches used in METBK bulk flux calculations
"""
# The jcool and jwarm switches should be set to 1, always!
JCOOLFL = 1      # 1=do coolskin calc
JWARMFL = 1      # 1=do warmlayer calc
#JWAVEFL         # only the windspeed parametrization of the charnok
                 # variable is coded; therefore this switch is not used.

"""
    LISTING OF SUBROUTINES BY ORDER IN THIS MODULE
        Grouped by sections; alphabetical within each section.

        The functions which directly calculate data products, both formal and meta,
        are listed here as the actual data product name in upper case, rather than
        by the name of the function; the functions themselves are named as
        "met_prdname" except as noted.

        All other functions are listed by function name.
#...................................................................................
    Functions to compute the L1 BULKMET (METBK) data products:
    these do not require the 'warmlayer/coolskin' iteration algorithm:
        BARPRES
        WINDAVG-VLE  (function name: met_windavg_mag_corr_east)
        WINDAVG-VLN  (function name: met_windavg_mag_corr_north)
    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute).
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) metadata products that do not require
    the 'warmlayer/coolskin' iteration algorithms:
        CURRENT_DIR
        CURRENT_SPD
        RELWIND_DIR-AUX
        RELWIND_SPD-AUX
        TIMEFLX-AUX
    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute), EXCEPT for TIMEFLX-AUX (hourly).
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) L2 METBK data products that do not require
    the 'warmlayer/coolskin' iteration algorithm:
        NETSIRR (this may operationally be an L1 product)
        NETSIRR_HOURLY (on an hourly time base)
        RAINRTE
        SALSURF
        SPECHUM
    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute) except as noted.
#...................................................................................
#...................................................................................
    Functions to compute the L2 METBK data products that do require
    the 'warmlayer/coolskin' iteration algorithm:
        BUOYFLS:  added DPA to match FDCHP, not in original DPS
        BUOYFLX:  added DPA to match FDCHP, not in original DPS
        FRSHFLX
        HEATFLX
        LATNFLX
        MOMMFLX
        NETLIRR
        RAINFLX
        SENSFLX
        SPHUM2M
        STABLTY:  metadata
        TEMPA2M
        TEMPSKN:  metadata
        WIND10M
    These products are calculated on hourly averages.
#...................................................................................
#...................................................................................
    Functions to compute the L2 METBK data products (that do require
    the 'warmlayer/coolskin' iteration algorithm) at the native temporal
    resolution of the instrument suite (roughly per minute).

    These functions were not specified in the DPS.
        HEATFLX_MINUTE
        LATNFLX_MINUTE
        NETLIRR_MINUTE
        SENSFLX_MINUTE
#...................................................................................
#...................................................................................
    Simple subroutines used in the routines in the sections above.
        air_density
        airtemp_at_refheight
        calc_rain_rate
        gravity
        latent_heat_vaporization_pure_water
        net_longwave_up
        psit_26
        psiu_26
        rain_heat_flux
        sea_spechum
        spechum_at_refheight
        water_thermal_expansion
        windspeed_at_refheight
#...................................................................................
#...................................................................................
    seasurface_skintemp_correct  (wrapper; calls warmlayer and coare35vn)
#...................................................................................
#...................................................................................
    warmlayer ('warmlayer' toga-coare routine)
#...................................................................................
#...................................................................................
    coare35vn (bulk calculation + 'coolskin' toga-coare routines; plus subroutines)
        charnock_wind
        coolskin_parameters
        effective_relwind
        obukhov_for_init
        obukhov_length_scale
        roughness_lengths
        roughness_lengths_for_init
        scaling_parameters
#...................................................................................
#...................................................................................
    Data conditioning and averaging routines
        vet_velptmn_data
        condition_data
        make_hourly_data
        warmlayer_time_keys
#...................................................................................

"""

####################################################################################
####################################################################################
####################################################################################

"""
#...................................................................................
#...................................................................................
    Functions to compute the L1 BULKMET (METBK) data products:
    do not require the 'warmlayer/coolskin' iteration algorithm:

        BARPRES
        WINDAVG-VLE  (met_windavg_mag_corr_east)
        WINDAVG-VLN  (met_windavg_mag_corr_north)

    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute).
#...................................................................................
#...................................................................................

"""


def met_barpres(mbar):
    """
    Computes BARPRES_L1, the OOI Level 1 barometric pressure core data
    product, by scaling the measured barometric pressure from mbar to
    Pascals.

    Parameters
    ----------
    mbar : array_like
        Barometric pressure (BARPRES_L0) [mbar].

    Returns
    -------
    Pa : array_like
        Barometric pressure (BARPRES_L1) [Pa].
    """
    Pa = mbar * 100.
    return Pa


def met_windavg_mag_corr_east(uu, vv, lat, lon, timestamp, spd_corr=[0.0, 1.0], zwindsp=0.0):
    """
    OOI single-output wrapper for WINDAVG-VLE_L1. Returns the METBK
    eastward wind speed [m/s] corrected for magnetic declination and
    wind speed under-reporting at higher wind speeds.

    See Also
    --------
    met_wind_mag_corr : Core algorithm; use directly for multi-output
        access to both corrected wind components.
    """
    # calculate the magnetic declination using the WMM model and rotate the vectors
    # from the magnetic to the true compass frame using met_wind_mag_corr
    uu_cor, vv_cor = met_wind_mag_corr(uu, vv, lat, lon, timestamp, spd_corr, zwindsp)
    return uu_cor


def met_windavg_mag_corr_north(uu, vv, lat, lon, timestamp, spd_corr=[0.0, 1.0], zwindsp=0.0):
    """
    OOI single-output wrapper for WINDAVG-VLN_L1. Returns the METBK
    northward wind speed [m/s] corrected for magnetic declination and
    wind speed under-reporting at higher wind speeds.

    See Also
    --------
    met_wind_mag_corr : Core algorithm; use directly for multi-output
        access to both corrected wind components.
    """
    # calculate the magnetic declination using the WMM model and rotate the vectors
    # from the magnetic to the true compass frame using met_wind_mag_corr
    uu_cor, vv_cor = met_wind_mag_corr(uu, vv, lat, lon, timestamp, spd_corr, zwindsp)
    return vv_cor


def met_wind_mag_corr(uu, vv, lat, lon, timestamp, spd_corr=[0.0, 1.0], zwindsp=0.0):
    """
    Computes WINDAVG_L1 (eastward and northward), the OOI Level 1 core
    wind speed data product, by correcting the METBK wind measurement
    for magnetic declination and for wind speed under-reporting at
    higher wind speeds.

    Parameters
    ----------
    uu : array_like
        Eastward wind speed (WINDAVG-VLE_L0), uncorrected [m/s].
    vv : array_like
        Northward wind speed (WINDAVG-VLN_L0), uncorrected [m/s].
    lat : float
        Instrument deployment latitude [decimal degrees].
    lon : float
        Instrument deployment longitude [decimal degrees].
    timestamp : array_like
        Sample date and time [seconds since 1900-01-01].
    spd_corr : array_like, optional
        Wind speed correction coefficients, shape (n, 2) for a linear
        offset and slope, or (n, 4) for a piecewise-linear correction
        (offset, slope, slope change, threshold). Default [0.0, 1.0]
        applies no correction.
    zwindsp : float, optional
        Height of the wind speed sensor above sea level [m]. Default
        0.0.

    Returns
    -------
    uu_cor : array_like
        Eastward wind speed (WINDAVG-VLE_L1), corrected for magnetic
        declination and wind speed under-reporting [m/s].
    vv_cor : array_like
        Northward wind speed (WINDAVG-VLN_L1), corrected for magnetic
        declination and wind speed under-reporting [m/s].

    Notes
    -----
    Magnetic declination is calculated using the IGRF-14 model via
    generic_functions.magnetic_declination.
    """
    # calculate the magnetic declination for the site
    zflag = 1  # denotes that z is a height above sea level.
    mag_dec = magnetic_declination(lat, lon, timestamp, zwindsp, zflag)

    # rotate the vectors from the magnetic to the true compass frame
    magvar = np.vectorize(magnetic_correction)
    uu_cor, vv_cor = magvar(mag_dec, uu, vv)

    # convert the wind components to speed and direction
    wspd = np.sqrt(uu_cor**2 + vv_cor**2)
    wdir = np.arctan2(uu_cor, vv_cor)
    wdir = np.where(wdir < 0, wdir + np.pi * 2, wdir)  # 0 to 360 degrees, but still in radians

    # apply the wind speed correction factors (array specific and provided as calibration coefficients).
    # calibration coefficients are provided as an array of (n, 2) or (n, 4), where n is the number of records.
    if spd_corr.shape[1] == 2:
        # linear correction only (offset and slope, where the offset can be zero for a pure slope correction)
        wspd_cor = spd_corr[:, 0] + spd_corr[:, 1] * wspd
    elif spd_corr.shape[1] == 4:
        # piecewise linear correction (offset, slope, slope change, threshold)
         wspd_cor = spd_corr[:, 0] + spd_corr[:, 1] * wspd + spd_corr[:, 2] * np.maximum(wspd - spd_corr[:, 3], 0)
    else:
        raise ValueError('spd_corr must be a list of 2 or 4 values.')

    # ensure there are no negative wind speeds after correction
    wspd_cor = np.where(wspd_cor < 0, 0.0, wspd_cor)

    # with the adjusted wind speed, re-calculate the eastward and northward wind components
    uu_cor = wspd_cor * np.sin(wdir)  # eastward wind component
    vv_cor = wspd_cor * np.cos(wdir)  # northward wind component
    return uu_cor, vv_cor

"""
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) metadata products that do not require
    the warmlayer/coolskin iteration algorithms:
        CURRENT_DIR
        CURRENT_SPD
        RELWIND_DIR-AUX
        RELWIND_SPD-AUX
        TIMEFLX-AUX

    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute), EXCEPT for TIMEFLX-AUX (hourly).
#...................................................................................
#...................................................................................
"""
def met_current_direction(vle_water, vln_water, use_velptmn_with_metbk=0):
    """
    Computes CURRENT_DIR, the direction of the METBK surface current,
    from the VELPT eastward and northward velocity components.

    Parameters
    ----------
    vle_water : array_like
        Eastward surface current (VELPTMN-VLE_L1) [m/s].
    vln_water : array_like
        Northward surface current (VELPTMN-VLN_L1) [m/s].
    use_velptmn_with_metbk : array_like, optional
        Time-vectorized data quality flag: 0 for bad VELPTMN current
        data, 1 for good VELPTMN current data. Default 0.

    Returns
    -------
    current_dir : array_like
        Direction of the surface current (CURRENT_DIR) [0, 360)
        degrees.

    Notes
    -----
    Not used by any other function in this module; calculated for
    its scientific interest. See Notes to met_relwind_speed.
    """
    # replace aliased current values with nans.
    vle_water, vln_water = vet_velptmn_data(vle_water, vln_water, use_velptmn_with_metbk)

    # use arctan2, which will properly handle arguments in all 4 cartesian quadrants,
    # and gives answers [-180 180] after conversion to degrees.
    cartesian_dir = np.degrees(np.arctan2(vln_water, vle_water))

    # the angle above is defined for a cartesian (x,y) coordinate system:
    # an angle of 0 points along the positive x-axis ("east" instead of
    # "north"), and increasingly positive angles indicate a ccw rotation
    # (instead of clockwise, which is the compass heading convention).

    # to convert to the compass convention, flip the sign and then add 90 degrees.
    # to change the range of values to [0 360), add 360 and 'mod' the result.
    current_dir = np.mod(450 - cartesian_dir, 360)

    return current_dir


def met_current_speed(vle_water, vln_water, use_velptmn_with_metbk=0):
    """
    Computes CURRENT_SPD, the magnitude of the METBK surface current,
    from the VELPT eastward and northward velocity components.

    Parameters
    ----------
    vle_water : array_like
        Eastward surface current (VELPTMN-VLE_L1) [m/s].
    vln_water : array_like
        Northward surface current (VELPTMN-VLN_L1) [m/s].
    use_velptmn_with_metbk : array_like, optional
        Time-vectorized data quality flag: 0 for bad VELPTMN current
        data, 1 for good VELPTMN current data. Default 0.

    Returns
    -------
    current_spd : array_like
        Magnitude of the surface current (CURRENT_SPD) [m/s].

    Notes
    -----
    Not used elsewhere in this module; the vector difference of wind
    and current is used instead (see RELWIND_SPD-AUX). See Notes to
    met_relwind_speed.
    """
    # replace aliased current values with nans.
    vle_water, vln_water = vet_velptmn_data(vle_water, vln_water, use_velptmn_with_metbk)

    current_spd = np.sqrt(vle_water**2 + vln_water**2)
    return current_spd


def met_relwind_direction(vle_wind, vln_wind, vle_water=None, vln_water=None, use_velptmn_with_metbk=0):
    """
    Computes RELWIND_DIR-AUX, the direction of the vector difference
    between METBK wind velocity and VELPT surface current velocity.

    Parameters
    ----------
    vle_wind : array_like
        Eastward wind speed (WINDAVG-VLE_L1) [m/s].
    vln_wind : array_like
        Northward wind speed (WINDAVG-VLN_L1) [m/s].
    vle_water : array_like, optional
        Eastward surface current (VELPTMN-VLE_L1) [m/s]. If not
        supplied, current is treated as unavailable and nan is
        returned.
    vln_water : array_like, optional
        Northward surface current (VELPTMN-VLN_L1) [m/s]. If not
        supplied, current is treated as unavailable and nan is
        returned.
    use_velptmn_with_metbk : array_like, optional
        Time-vectorized data quality flag: 0 for bad VELPTMN current
        data, 1 for good VELPTMN current data. Default 0.

    Returns
    -------
    u_dir : array_like
        Direction of relative wind (RELWIND_DIR-AUX) [0, 360) degrees.

    Notes
    -----
    Not used elsewhere in this module; calculated for its scientific
    interest. See Notes to met_relwind_speed.
    """
    # if this function is called without using surface current data, return nan
    if vle_water is None or vln_water is None:
        u_dir = vle_wind * np.nan
        return u_dir

    # replace aliased current values with nans.
    vle_water, vln_water = vet_velptmn_data(vle_water, vln_water, use_velptmn_with_metbk)

    # use arctan2, which will properly handle arguments in all 4 cartesian quadrants,
    # and gives answers [-180 180] after conversion to degrees.
    cartesian_dir = np.degrees(np.arctan2(vln_wind - vln_water, vle_wind - vle_water))

    # the angle above is defined for a cartesian (x,y) coordinate system:
    # an angle of 0 points along the positive x-axis ("east" instead of
    # "north"), and increasingly positive angles indicate a ccw rotation
    # (instead of clockwise, which is the compass heading convention).

    # to convert to the compass convention, flip the sign and then add 90 degrees.
    # to change the range of values to [0 360), add 360 and 'mod' the result.
    u_dir = np.mod(450 - cartesian_dir, 360)

    return u_dir


def met_relwind_speed(vle_wind, vln_wind, vle_water=None, vln_water=None, use_velptmn_with_metbk=0):
    """
    Computes RELWIND_SPD-AUX, wind speed relative to water, as the
    magnitude of the vector difference between METBK wind velocity and
    VELPT surface current velocity. This is the fundamental wind speed
    variable used by the METBK toga-coare algorithms.

    Parameters
    ----------
    vle_wind : array_like
        Eastward wind speed (WINDAVG-VLE_L1) [m/s].
    vln_wind : array_like
        Northward wind speed (WINDAVG-VLN_L1) [m/s].
    vle_water : array_like, optional
        Eastward surface current (VELPTMN-VLE_L1) [m/s]. If not
        supplied, current is treated as 0.
    vln_water : array_like, optional
        Northward surface current (VELPTMN-VLN_L1) [m/s]. If not
        supplied, current is treated as 0.
    use_velptmn_with_metbk : array_like, optional
        Time-vectorized data quality flag: 0 for bad VELPTMN current
        data, 1 for good VELPTMN current data. Default 0.

    Returns
    -------
    u_rel : array_like
        Magnitude of wind speed relative to water (RELWIND_SPD-AUX)
        [m/s].

    Notes
    -----
    Unlike met_relwind_direction, met_current_direction, and
    met_current_speed, bad or missing current data are set to 0 (not
    nan) here so actual wind speed is used in place of relative wind
    speed, per the DPS. See docs site for full discussion.
    """
    # If the surface current velocities are missing or invalid, the actual windspeed
    # will be used in place of the relative windspeed over water to calculate the METBK
    # data products.
    #
    # if this function is called without using surface current data, set current data to 0.
    if vle_water is None or vln_water is None:
        vle_water = np.zeros(vle_wind.shape[0])
        vln_water = np.zeros(vle_wind.shape[0])

    # find nans in the current record and if found replace both the east
    # and north components with values of 0.
    nanmask = np.isnan(vle_water * vln_water)
    vle_water[nanmask] = 0.0
    vln_water[nanmask] = 0.0

    # expand use_velptmn_with_metbk if it is called as a scalar
    if np.atleast_1d(use_velptmn_with_metbk).shape[0] == 1:
        use_velptmn_with_metbk = np.tile(use_velptmn_with_metbk, vle_wind.shape[0])

    # vet the surface current data - but don't use Nans
    #   when use_velptmn_with_metbk=0, set surface current velocities to 0.
    #   when use_velptmn_with_metbk=1, use surface current velocities as received.
    vle_water = vle_water * use_velptmn_with_metbk
    vln_water = vln_water * use_velptmn_with_metbk

    u_rel = np.sqrt((vle_water - vle_wind)**2 + (vln_water - vln_wind)**2)
    return u_rel


@deprecated
def met_timeflx(timestamp):
    """
    Computes TIMEFLX-AUX, the UTC timestamps corresponding to the
    hourly averaged METBK data products [seconds since 1900-01-01].

    Parameters
    ----------
    timestamp : array_like
        Sample date and time [seconds since 1900-01-01].

    Returns
    -------
    fluxtime_hourly : array_like
        UTC timestamp for hourly data [seconds since 1900-01-01].

    Notes
    -----
    Timestamps mark the midpoint of each hourly bin, starting a half
    hour after the first data record. See make_hourly_data.
    """
    # here, the output of make_hourly_data is a list,
    # the only element of which is the desired rank 1 np.array
    fluxtime_hourly = make_hourly_data(timestamp)[0]
    return fluxtime_hourly


"""
#...................................................................................
#...................................................................................
    Functions to compute the (simpler) L2 METBK data products that do not require
    the 'warmlayer/coolskin' iteration algorithm:

        NETSIRR (this may operationally be an L1 product)
        NETSIRR_HOURLY (hourly time base)
        RAINRTE
        SALSURF
        SPECHUM

    These products are calculated at the native temporal resolution of the
    instrument suite (roughly each minute) except as noted.
#...................................................................................
#...................................................................................

"""


def met_netsirr(shortwave_down):
    """
    Computes NETSIRR_L2, the net downward shortwave radiation (0.3 to
    3.0 um wavelengths), by subtracting the reflected component from
    the measured downward shortwave irradiance. Calculated on the
    native METBK timebase (roughly per minute).

    Parameters
    ----------
    shortwave_down : array_like
        Measured downward shortwave radiation (SHRTIRR_L1) [W/m^2].

    Returns
    -------
    net_shortwave_down : array_like
        Net downward shortwave radiation (NETSIRR_L2) [W/m^2].

    Notes
    -----
    Uses a fixed albedo of 0.055, the reflection coefficient used in
    the original toga-coare code.
    """
    # net down = total down - reflected up
    albedo = 0.055    # value for reflection coefficient used in toga-coare code
    net_shortwave_down = (1.0 - albedo) * shortwave_down

    return net_shortwave_down


@deprecated
def met_netsirr_hourly(shortwave_down, timestamp):
    """
    Computes NETSIRR_HOURLY_L2, the net downward shortwave radiation
    (0.3 to 3.0 um wavelengths), binned to an hourly timebase. Does
    not require the coolskin/warmlayer algorithms.

    Parameters
    ----------
    shortwave_down : array_like
        Measured downward shortwave radiation (SHRTIRR_L1) [W/m^2].
    timestamp : array_like
        Sample date and time [seconds since 1900-01-01].

    Returns
    -------
    net_shortwave_down_hourly : array_like
        Net downward shortwave radiation, hourly (NETSIRR_HOURLY_L2)
        [W/m^2].
    """
    shortwave_down, timestamp = condition_data(shortwave_down, timestamp)

    (shortwave_down_hourly, _) = make_hourly_data(shortwave_down, timestamp)

    net_shortwave_down_hourly = met_netsirr(shortwave_down_hourly)

    return net_shortwave_down_hourly


@deprecated
def met_rainrte(cumulative_precipitation, timestamp):
    """
    Computes RAINRTE_L2, rain rate, binned to an hourly timebase.

    Parameters
    ----------
    cumulative_precipitation : array_like
        Measured cumulative rain level (PRECIPM_L1) [mm].
    timestamp : array_like
        Sample date and time [seconds since 1900-01-01].

    Returns
    -------
    rainrte : array_like
        Rain rate, hourly (RAINRTE_L2) [mm/hr].

    Notes
    -----
    Likely an L1 product despite the L2 designation in the DPS.
    """
    cumulative_precipitation, timestamp = condition_data(cumulative_precipitation, timestamp)

    # trap out scalar case; return a value of 0 as does the DPS code.
    if cumulative_precipitation.size == 1:
        return 0.0

    cumu_prcp_hrly, time_hrly = make_hourly_data(cumulative_precipitation, timestamp)

    rainrte = calc_rain_rate(cumu_prcp_hrly, time_hrly)

    return rainrte


def met_salsurf(cond, tC_sea, ztmpwat):
    """
    Computes SALSURF_L2, the OOI Level 2 sea surface practical
    salinity (PSS-78), using TEOS-10 (GSW) with METBK conductivity
    and temperature measurements.

    Parameters
    ----------
    cond : array_like
        Sea surface conductivity (CONDSRF_L1) [S/m].
    tC_sea : array_like
        Sea surface temperature (TEMPSRF_L1) [deg_C].
    ztmpwat : array_like
        Depth of the conductivity and temperature measurements,
        serving as a proxy for pressure [m].

    Returns
    -------
    SP : array_like
        Practical salinity, PSS-78 (SALSURF_L2) [unitless].
    """
    # Convert L1 Conductivity from S/m to mS/cm
    C10 = cond * 10.0

    # Calculate the Practical Salinity (PSS-78) [unitless]
    SP = gsw.SP_from_C(C10, tC_sea, ztmpwat)
    return SP


def met_spechum(tC_air, pr_air, relhum):
    """
    Computes SPECHUM_L2, the OOI air specific humidity core data
    product. Not to be confused with SPHUM2M_L2.

    Parameters
    ----------
    tC_air : array_like
        Air temperature (TEMPAIR_L1) [deg_C].
    pr_air : array_like
        Air pressure (BARPRES_L0) [mbar]. Note this is BARPRES_L0,
        not BARPRES_L1 [Pa].
    relhum : array_like
        Relative humidity (RELHUMI_L1) [%].

    Returns
    -------
    q_air : array_like
        Air specific humidity (SPECHUM_L2) [g/kg].
    """
    # calculate saturated vapor pressure es in mbar
    es = 6.1121 * np.exp(17.502 * tC_air/(tC_air+240.97)) * (1.0007 + 3.46e-6 * pr_air)
    # calculate vapor pressure em from definition of relative humidity
    em = 0.01 * relhum * es
    # specific humidity q [g/kg] is then:
    q_air = 621.97 * em/(pr_air - 0.378 * em)
    return q_air


"""
#...................................................................................
#...................................................................................
    Functions to compute the L2 METBK data products that do require
    the 'warmlayer/coolskin' iteration algorithm:

        BUOYFLS:  added DPA to match FDCHP, not in original DPS
        BUOYFLX:  added DPA to match FDCHP, not in original DPS
        FRSHFLX
        HEATFLX
        LATNFLX
        MOMMFLX
        NETLIRR
        RAINFLX
        SENSFLX
        SPHUM2M
        STABLTY:  metadata
        TEMPA2M
        TEMPSKN:  metadata
        WIND10M

    These products are calculated on hourly averages.
#...................................................................................
#...................................................................................

"""


@deprecated
def met_buoyfls(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for BUOYFLS_L2. Returns the sonic
    buoyancy flux [W/m^2], computed using sonic temperature rather
    than virtual temperature. Not specified in the original DPS.

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    c2k = 273.15   # celsius to kelvin temperature constant
    tssr = tsr + 0.51 * (tC_air + c2k) * qsr
    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    # sonic buoyancy flux
    hsbb = -rhoa * cpa * usr * tssr

    return hsbb


@deprecated
def met_buoyflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for BUOYFLX_L2. Returns the buoyancy
    flux [W/m^2], computed using virtual temperature. Not specified
    in the original DPS.

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    c2k = 273.15   # celsius to kelvin temperature constant
    tvsr = tsr + 0.61 * (tC_air + c2k) * qsr
    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    # buoyancy flux
    hbb = -rhoa * cpa * usr * tvsr

    return hbb


@deprecated
def met_frshflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for FRSHFLX_L2. Returns the upward
    freshwater flux [mm/hr].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, _, qsr, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    # jim edson uses freshwater density; whoi dps uses seawater density;
    # perhaps the w/v concentration of pure water in seawater should be used
    # (which would be < 1000 kg/m^3).
    rho_purewater = 1000.0  # kg/m^3
    # the factor of 1000 converts from m -> mm; 3600, per sec -> per hr.
    evap = -rhoa * usr * qsr / rho_purewater * 1000.0 * 3600.0    # [mm/hr]
    frshflx = evap - rain_rate

    return frshflx


@deprecated
def met_heatflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for HEATFLX_L2. Returns the total net
    upward heat flux [W/m^2].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    met_heatflx_minute : Same calculation on the native per-minute
        timebase, not binned to hourly.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, qsr, _, dter, dqer, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, Rshort_down,
        Rlong_down, _, _, _, _) = args

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    rhoa = air_density(tC_air, pr_air, relhum)
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr                                              # positive up
    hsb = -rhoa * cpa * usr * tsr                                             # positive up
    Rns_down = met_netsirr(Rshort_down)                                       # positive down
    Rnl_up = net_longwave_up(tC_sea + dsea - dter, Rlong_down)                # positive up
    rainflx = rain_heat_flux(rain_rate, tC_sea+dsea, tC_air, relhum, pr_air)  # positive up

    heatflx = hlb + hsb - Rns_down + Rnl_up + rainflx

    return heatflx


@deprecated
def met_latnflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for LATNFLX_L2. Returns the upward
    latent heat flux [W/m^2].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    met_latnflx_minute : Same calculation on the native per-minute
        timebase, not binned to hourly.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dsea is the warmlayer correction to the sea surface temperature
    (usr, _, qsr, _, _, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    # note that the original (coare ver. 3.5) code:
    #    (a) uses Le for pure water, not seawater.
    #    (b) does not include the coolskin correction to sea surface temperature.
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr

    return hlb


@deprecated
def met_mommflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for MOMMFLX_L2. Returns the absolute
    value of the momentum flux, also called the wind stress tau
    [N/m^2].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, _, _, ut, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, wnd, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    # the wind stress tau is the magnitude of the momentum flux.
    tau = rhoa * usr * usr * wnd / ut

    return tau


@deprecated
def met_netlirr(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for NETLIRR_L2. Returns the net upward
    longwave irradiance [W/m^2].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    met_netlirr_minute : Same calculation on the native per-minute
        timebase, not binned to hourly.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dter is the coolskin temperature depression [degC]
    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, dter, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, _, _, _, _, _, _, Rlong_down, _, _, _, _) = args

    Rnl = net_longwave_up(tC_sea + dsea - dter, Rlong_down)

    return Rnl


@deprecated
def met_rainflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for RAINFLX_L2. Returns the net upward
    rain heat flux [W/m^2].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    rain_heat_flux : Underlying rain heat flux calculation.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, _, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air,
        _, _, _, _, _, _) = args

    # the raindrops penetrate the sea surface on the order of cm, which is where the
    # heat is 'exchanged'. therefore, use the warmlayer correction but not the
    # coolskin correction (which is order microns (?) thick) to the sea temperature.
    rainflx = rain_heat_flux(rain_rate, tC_sea+dsea, tC_air, relhum, pr_air)

    return rainflx


@deprecated
def met_sensflx(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for SENSFLX_L2. Returns the net upward
    sensible heat flux [W/m^2].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    met_sensflx_minute : Same calculation on the native per-minute
        timebase, not binned to hourly.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, _, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    hsb = -rhoa * cpa * usr * tsr

    return hsb


@deprecated
def met_sphum2m(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for SPHUM2M_L2. Returns the modelled
    specific humidity at a reference height of 2 m [g/kg].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    zrefht = 2.0  # [m]

    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # L is the Obukhov length scale [m]
    (_, _, qsr, _, _, _, _, L, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, zhumair, pr_air, _, _, _, _, _, _) = args

    sphum2m = spechum_at_refheight(tC_air, pr_air, relhum, qsr, zrefht, zhumair, L)

    return sphum2m


@deprecated
def met_stablty(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for STABLTY_L2. Returns the
    Monin-Obukhov stability parameter [unitless].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # L is the Obukhov length scale [m]
    (_, _, _, _, _, _, _, L, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, zwindsp, _, _, _, _, _, _, _, _, _, _, _) = args

    return zwindsp / L


@deprecated
def met_tempa2m(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for TEMPA2M_L2. Returns the modelled
    air temperature at a reference height of 2 m [degC].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    zrefht = 2.0  # [m]

    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # L is the Obukhov length scale [m]
    (_, tsr, _, _, _, _, _, L, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, ztmpair, _, _, _, _, _, lat, _, _, _) = args

    tempa2m = airtemp_at_refheight(tC_air, tsr, zrefht, ztmpair, L, lat)

    return tempa2m


@deprecated
def met_tempskn(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for TEMPSKN_L2. Returns the skin sea
    surface temperature [degC], from the warmlayer and coolskin
    (coare35vn) model.

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dter is the coolskin temperature depression [degC]
    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, dter, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, _, _, _, _, _, _, _, _, _, _, _) = args

    # warmlayer corrections are added; coolskin corrections are subtracted
    tempskn = tC_sea + dsea - dter

    return tempskn


@deprecated
def met_wind10m(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for WIND10M_L2. Returns the modelled
    wind speed at a reference height of 10 m [m/s].

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    """
    zrefht = 10.0  # [m]

    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args = make_hourly_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # L is the Obukhov length scale [m]
    (usr, _, _, ut, _, _, _, L, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, wnd, zwindsp, _, _, _, _, _, _, _, _, _, _, _) = args

    wind10m = windspeed_at_refheight(wnd, usr, ut, zrefht, zwindsp, L)

    return wind10m

"""
#...................................................................................
#...................................................................................
    Functions to compute the L2 METBK data products that do require
    the 'warmlayer/coolskin' iteration algorithm, but at the native
    temporal resolution of the instrument suite (roughly per minute).

    These functions were not specified in the DPS.
        HEATFLX_MINUTE
        LATNFLX_MINUTE
        NETLIRR_MINUTE
        SENSFLX_MINUTE
#...................................................................................
#...................................................................................
"""


@deprecated
def met_heatflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                       zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                       Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                       zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for HEATFLX_MINUTE_L2. Returns the
    total net upward heat flux [W/m^2] on the native METBK per-minute
    timebase. Not specified in the DPS.

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    met_heatflx : Same calculation binned to hourly averages.

    Notes
    -----
    Differs from met_heatflx only in that make_hourly_data is not
    called to bin the input data into hourly bins.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, qsr, _, dter, dqer, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (rain_rate, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, Rshort_down,
        Rlong_down, _, _, _, _) = args

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    rhoa = air_density(tC_air, pr_air, relhum)
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr                                              # positive up
    hsb = -rhoa * cpa * usr * tsr                                             # positive up
    Rns_down = met_netsirr(Rshort_down)                                       # positive down
    Rnl_up = net_longwave_up(tC_sea + dsea - dter, Rlong_down)                # positive up
    rainflx = rain_heat_flux(rain_rate, tC_sea+dsea, tC_air, relhum, pr_air)  # positive up

    heatflx = hlb + hsb - Rns_down + Rnl_up + rainflx

    return heatflx


@deprecated
def met_latnflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                       zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                       Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                       zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for LATNFLX_MINUTE_L2. Returns the
    upward latent heat flux [W/m^2] on the native METBK per-minute
    timebase. Not specified in the DPS.

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    met_latnflx : Same calculation binned to hourly averages.

    Notes
    -----
    Differs from met_latnflx only in that make_hourly_data is not
    called to bin the input data into hourly bins.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dsea is the warmlayer correction to the sea surface temperature
    (usr, _, qsr, _, _, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    # note that the original (coare ver. 3.5) code:
    #    (a) uses Le for pure water, not seawater.
    #    (b) does not include the coolskin correction to sea surface temperature.
    Le = latent_heat_vaporization_pure_water(tC_sea + dsea)

    hlb = -rhoa * Le * usr * qsr

    return hlb


@deprecated
def met_netlirr_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                       zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                       Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                       zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for NETLIRR_MINUTE_L2. Returns the net
    upward longwave irradiance [W/m^2] on the native METBK per-minute
    timebase. Not specified in the DPS.

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    met_netlirr : Same calculation binned to hourly averages.

    Notes
    -----
    Differs from met_netlirr only in that make_hourly_data is not
    called to bin the input data into hourly bins.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    # dter is the coolskin temperature depression [degC]
    # dsea is the warmlayer correction to the sea surface temperature [degC]
    (_, _, _, _, dter, _, _, _, _, _, _, _, _, dsea) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, tC_sea, _, _, _, _, _, _, _, _, Rlong_down, _, _, _, _) = args

    Rnl_native = net_longwave_up(tC_sea + dsea - dter, Rlong_down)

    return Rnl_native


@deprecated
def met_sensflx_minute(tC_sea, wnd, tC_air, relhum, timestamp, lon, ztmpwat,
                       zwindsp, ztmpair, zhumair, lat=45.0, pr_air=1013.0,
                       Rshort_down=150.0, Rlong_down=370.0, cumu_prcp=0.0,
                       zinvpbl=600.0, jwarm=JWARMFL, jcool=JCOOLFL):
    """
    OOI single-output wrapper for SENSFLX_MINUTE_L2. Returns the net
    upward sensible heat flux [W/m^2] on the native METBK per-minute
    timebase. Not specified in the DPS.

    See Also
    --------
    seasurface_skintemp_correct : Core algorithm; use directly for
        multi-output access.
    met_sensflx : Same calculation binned to hourly averages.

    Notes
    -----
    Differs from met_sensflx only in that make_hourly_data is not
    called to bin the input data into hourly bins.
    """
    # package input arguments.
    # 1st 4 arguments are warmlayer, followed by coolskin, then switches.
    args = [cumu_prcp, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp,
            tC_air, ztmpair, relhum, zhumair, pr_air, Rshort_down,
            Rlong_down, lat, zinvpbl, jcool, jwarm]

    args = condition_data(*args)

    args[0] = calc_rain_rate(*args[0:2])

    (usr, tsr, _, _, _, _, _, _, _, _, _, _, _, _) = seasurface_skintemp_correct(*args)

    # make the necessary processed hourly data available for the final calculation
    (_, _, _, _, _, _, _, tC_air, _, relhum, _, pr_air, _, _, _, _, _, _) = args

    rhoa = air_density(tC_air, pr_air, relhum)

    cpa = 1004.67  # specific heat capacity of (dry) air [J/kg/K]
    hsb = -rhoa * cpa * usr * tsr

    return hsb


"""
#...................................................................................
#...................................................................................
    Simple subroutines used in the routines in the sections above.
        Does NOT include:
            the routines condition_data and make_hourly_data
            warmlayer and coolskin (coare35vn) routines

    air_density
    airtemp_at_refheight
    calc_rain_rate
    gravity
    latent_heat_vaporization_pure_water
    net_longwave_up
    psit_26
    psiu_26
    rain_heat_flux
    sea_spechum
    spechum_at_refheight
    water_thermal_expansion
    windspeed_at_refheight
#...................................................................................
#...................................................................................

"""


def air_density(tC_air, air_pressure_mbar, relhum):
    """
    Returns the density of air.

    Parameters
    ----------
    tC_air : array_like
        Air temperature (TEMPAIR_L1) [deg_C].
    air_pressure_mbar : array_like
        Barometric pressure (BARPRES_L0) [mbar].
    relhum : array_like
        Relative humidity of air (RELHUMI_L1) [%].

    Returns
    -------
    rho_air : array_like
        Density of air [kg/m^3].
    """
    Rgas = 287.05  # gas constant [J/kg/K] for dry(!) air
    c2k = 273.15  # celsius to kelvin temperature constant
    sp_hum_air = met_spechum(tC_air, air_pressure_mbar, relhum)/1000.0   # units kg/kg
    # the factor of 100 converts pressure from mbar to pascal [N/m^2]
    rho_air = air_pressure_mbar * 100.0 / (1.0 + 0.61 * sp_hum_air) / (Rgas * (tC_air + c2k))

    return rho_air


def airtemp_at_refheight(tC_air, tsr, zrefht, ztmpair, L, lat):
    """
    Computes air temperature at an arbitrary reference height using
    the Monin-Obukhov similarity profile.

    Parameters
    ----------
    tC_air : array_like
        Air temperature at the measurement height [deg_C].
    tsr : array_like
        Temperature scaling parameter [K].
    zrefht : float
        Reference height at which to compute air temperature [m].
    ztmpair : array_like
        Height of the air temperature measurement [m].
    L : array_like
        Obukhov length scale [m].
    lat : float
        Latitude, used to compute the lapse rate correction [deg].

    Returns
    -------
    temp : array_like
        Air temperature at zrefht [deg_C].

    Notes
    -----
    Equivalent to DPS section A.3.3, with the addition of a lapse-rate
    term present in J. Edson's coare35vn v3.5 code.
    """
    von = 0.4      # von Karman constant
    cpa = 1004.67  # heat capacity of dry air [J/kg/K]

    lapse = gravity(lat)/cpa
    temp = tC_air + tsr / von * (np.log(zrefht/ztmpair) - psit_26(zrefht/L) +
                                 psit_26(ztmpair/L)) + lapse * (ztmpair-zrefht)

    return temp


def calc_rain_rate(cumulative_precipitation, timestamp):
    """
    Computes rain rate from cumulative precipitation and timestamp.
    Not itself the formal RAINRTE data product, but used by
    met_rainrte and by the warmlayer routine.

    Parameters
    ----------
    cumulative_precipitation : array_like
        Measured cumulative rain level (PRECIPM_L1) [mm].
    timestamp : array_like
        Sample date and time [seconds since 1900-01-01].

    Returns
    -------
    rain_rate : array_like
        Rain rate [mm/hr].

    Notes
    -----
    Negative rates (from sensor resets or evaporation exceeding
    precipitation) are set to 0, per the DPS.
    """
    # trap out scalar case; return a value of 0
    if cumulative_precipitation.size == 1:
        return 0.0

    # calculate the amount of rain fallen in each time interval
    rainfall = np.diff(cumulative_precipitation)   # [mm]
    # calculate each time interval and convert to hours
    delta_time = np.diff(timestamp) / 3600.0
    # calculate rainrate
    rain_rate = rainfall / delta_time
    # prepend a 0 to represent the first (unknown) rainrate value.
    rain_rate = np.hstack((0.0, rain_rate))
    # follow the DPS and set all negative values to 0.
    #     values of 0 could arise when (1) the sensor resets its rainlevel
    #     to 0 by draining the accumulated rainwater, or when (2) the
    #     evaporation is greater than the precipitation.
    rain_rate[np.less(rain_rate, 0.0)] = 0.0

    return rain_rate


def gravity(lat):
    """
    Returns acceleration due to gravity as a function of latitude.

    Parameters
    ----------
    lat : array_like
        Latitude [degrees].

    Returns
    -------
    g : array_like
        Acceleration due to gravity [m/s/s].

    Notes
    -----
    From grv.m in the TOGA COARE 3.0 MATLAB toolbox.
    """
    gamma = 9.7803267715
    c1 = 0.0052790414
    c2 = 0.0000232718
    c3 = 0.0000001262
    c4 = 0.0000000007

    phi = np.radians(lat)
    x = np.sin(phi)
    xsq = x * x
    g = gamma * (1.0 + xsq * (c1 + xsq * (c2 + xsq * (c3 + xsq * c4))))

    return g


def latent_heat_vaporization_pure_water(tC_water):
    """
    Returns the latent heat of vaporization of pure water.

    Parameters
    ----------
    tC_water : array_like
        Water temperature [deg_C].

    Returns
    -------
    Le_water : array_like
        Latent heat of vaporization of pure water [J/kg].
    """
    return (2500.8 - 2.37 * tC_water) * 1000.0


def net_longwave_up(tC_water, total_longwave_down):
    """
    Computes the net upward longwave radiation flux from water
    temperature and downward longwave radiation.

    Parameters
    ----------
    tC_water : array_like
        Water temperature [deg_C].
    total_longwave_down : array_like
        Measured downward longwave radiation (LONGIRR_L1), positive
        downward [W/m^2].

    Returns
    -------
    Rnl : array_like
        Net longwave radiation, positive upward [W/m^2].

    Notes
    -----
    Does not by itself compute NETLIRR_L2, which requires the sea
    surface skin temperature corrected for warmlayer and coolskin
    effects. Uses a fixed blackbody emissivity of 0.97.
    """
    sigma = 5.67e-8  # Stefan-Boltzmann constant [W/(m^2 K^4)]
    eps = 0.97
    c2k = 273.15  # degC to kelvin conversion constant
    Rnl = eps * (sigma * (tC_water + c2k) ** 4 - total_longwave_down)

    return Rnl


def psit_26(zet):
    """
    Computes the temperature structure function, used to calculate
    air temperature and specific humidity at reference heights above
    sea level.

    Parameters
    ----------
    zet : array_like
        Monin-Obukhov stability parameter (zu/L) [dimensionless].

    Returns
    -------
    psit : array_like
        Temperature structure function value at zet.
    """
    # force the shape of the input to a 1D array.
    zet = np.atleast_1d(zet)

    # stable case: zet > 0.
    # calculate for all zet, and overwrite zet < 0 cases in last section, as
    # done in the original code.
    #
    # because negative zet cases will be overwritten, can take the absolute value
    # of zet underneath the fractional exponent (1.5) without affecting the end
    # result. this is done to avoid the possibility of raising a negative number
    # to a fractional power which will result in complex values.
    dzet = np.minimum(50.0, 0.35 * zet)
    psit = -((1.0 + 0.6667 * np.abs(zet))**1.5 +
             0.6667 * (zet - 14.28) * np.exp(-dzet) + 8.525)

    # overwrite psit for zet < 0 values (unstable case).
    # trap out nans; psit already has nans where we want them
    zet[np.isnan(zet)] = 1.0
    k = zet < 0.0
    x = (1.0 - 16.0 * zet[k])**0.5
    psik = 2.0 * np.log((1.0 + x) / 2.0)
    x = (1.0 - 34.15 * zet[k])**0.3333
    psic = (1.5 * np.log((1.0 + x + x**2) / 3.0) -
            np.sqrt(3) * np.arctan((1.0 + 2.0 * x) / np.sqrt(3)) +
            4.0 * np.arctan(1) / np.sqrt(3))
    f = zet[k]**2 / (1.0 + zet[k]**2)
    psit[k] = (1.0 - f) * psik + f * psic

    return psit


def psiu_26(zet):
    """
    Computes the velocity structure function, used to calculate wind
    speed at reference heights above sea level.

    Parameters
    ----------
    zet : array_like
        Monin-Obukhov stability parameter (zu/L) [dimensionless].

    Returns
    -------
    psiu : array_like
        Velocity structure function value at zet.
    """
    # force the shape of the input to a 1D array.
    zet = np.atleast_1d(zet)

    # stable case: zet > 0.
    # calculate for all zet, and overwrite zet < 0 cases in last section, as
    # done in the original code.
    dzet = np.minimum(50.0, 0.35 * zet)
    (a, b, c, d) = (0.7, 0.75, 5.0, 0.35)
    psiu = -(a * zet + b * (zet - c / d) * np.exp(-dzet) + b * c / d)

    # overwrite psiu for zet < 0 values (unstable case).
    # trap out nans; psiu already has nans where we want them
    zet[np.isnan(zet)] = 1.0
    k = zet < 0.0
    x = (1.0 - 16.0 * zet[k])**0.25
    psik = (2.0 * np.log((1.0 + x) / 2.0) + np.log((1.0 + x**2) / 2.0) -
            2.0 * np.arctan(x) + 2.0 * np.arctan(1))
    x = (1.0 - 10.15 * zet[k])**0.3333
    psic = (1.5 * np.log((1.0 + x + x**2) / 3.0) -
            np.sqrt(3) * np.arctan((1.0 + 2.0 * x) / np.sqrt(3)) +
            4.0 * np.arctan(1) / np.sqrt(3))
    f = zet[k]**2 / (1.0 + zet[k]**2)
    psiu[k] = (1.0 - f) * psik + f * psic

    return psiu


def rain_heat_flux(rainrate, Tsea, Tair, relhum, pr_air):
    """
    Computes the rain heat flux, the heat flux due to rain falling
    into the ocean; positive values indicate heat flowing from the
    ocean to the atmosphere (rain cools the ocean).

    Parameters
    ----------
    rainrate : array_like
        Rain fall rate [mm/hr].
    Tsea : array_like
        Bulk sea surface temperature, corrected for warmlayer only
        [deg_C].
    Tair : array_like
        Air temperature [deg_C].
    relhum : array_like
        Relative humidity [%].
    pr_air : array_like
        Air pressure [mb].

    Returns
    -------
    RHF : array_like
        Rain heat flux [W/m^2].

    Notes
    -----
    Follows Gosnell, Fairall, and Webster (1995), with the raindrop
    wetbulb temperature approximated by a Taylor series expansion at
    the air temperature rather than the sea temperature, per a
    correction suggested by Simon de Szoeke. See
    rain_heat_flux_FLAWED [removed] for the superseded derivation this
    replaced.
    """
    c2k = 273.15       # celsius to kelvin temperature constant
    Rgas_air = 287.05  # gas constant [J/kg/K] for air
    Rgas_wtr = Rgas_air / 0.622  # gas constant for water wapor [J/kg/K]
    cp_air = 1004.67   # specific heat capacity of dry air [J/kg/K]
    cp_rain = 4186.0   # specific heat capacity of freshwater at T=25 degC [J/kg/K]
    rho_rain = 1000.0  # density of freshwater [kg/m^3]

    rho_air = air_density(Tair, pr_air, relhum)               # units kg/m^3
    Lv_at_Tair = latent_heat_vaporization_pure_water(Tair)    # units J/kg
    # specific humidity of air
    qair = met_spechum(Tair, pr_air, relhum)/1000.0           # units kg/kg
    # saturation humidity of air at a temperature of Tair
    qsat_at_Tair = met_spechum(Tair, pr_air, 100.0)/1000.0    # units kg/kg

    # diffusivity expressions are taken from original code;
    # these equations were not checked.
    vapor_diffusivity = 2.11e-5 * ((Tair + c2k) / c2k) ** 1.94         # water vapour diffusivity
    heat_diffusivity = (1.0 + 3.309e-3 * Tair -
                        1.44e-6 * Tair * Tair) * 0.02411 / (rho_air * cp_air)  # heat diffusivity

    psi = Lv_at_Tair / cp_air * vapor_diffusivity / heat_diffusivity

    # Clausius-Clapeyron; the factor of 0.622 is included in Rgas_wtr
    CC_at_Tair = Lv_at_Tair * qsat_at_Tair / (Rgas_wtr * (Tair + c2k)**2)

    # rain temperature is assumed to be at the wet bulb temperature.
    Train = Tair - psi * (qsat_at_Tair - qair) / (1.0 + psi * CC_at_Tair)

    # most code versions omit rho_rain and the factor of 1/1000.
    RHF = rainrate * rho_rain * cp_rain * (Tsea - Train) / 1000.0 / 3600.0

    return RHF


def sea_spechum(tC_sea, p_air):
    """
    Computes the sea surface value of specific humidity, Qsea. Not to
    be confused with any of the specific humidity in air variables
    (e.g. Qair).

    Parameters
    ----------
    tC_sea : array_like
        Seawater temperature [deg_C].
    p_air : array_like
        Air pressure (BARPRES_L0) [mbar].

    Returns
    -------
    q_sea : array_like
        Sea surface specific humidity [g/kg].

    Notes
    -----
    Named qsat26sea in J. Edson's coare35vn v3.5 MATLAB code and qsee
    in the v3.0 DPS code.
    """
    # calculate saturated vapor pressure es in mbar
    es = 6.1121 * np.exp(17.502 * tC_sea/(tC_sea + 240.97)) * (1.0007 + 3.46e-6 * p_air)
    # the factor of 0.98 arises from the fact that the saturation vapor pressure above a
    # salt solution is less than that for pure water; Raoult's law (applicable for the
    # case of an ideal solution) specifies this factor as the mole fraction of water,
    # which for a salinity of 30 is 0.99. Empirically the correct value to use is 0.98.
    esat_sea = 0.98 * es
    # specific humidity q_sea [g/kg] is then:
    q_sea = 621.97 * esat_sea/(p_air - 0.378 * esat_sea)
    return q_sea


def spechum_at_refheight(tC_air, pr_air, relhum, qsr, zrefht, zhumair, L):
    """
    Computes specific humidity at an arbitrary reference height using
    the Monin-Obukhov similarity profile.

    Parameters
    ----------
    tC_air : array_like
        Air temperature at the measurement height [deg_C].
    pr_air : array_like
        Air pressure (BARPRES_L0) [mbar].
    relhum : array_like
        Relative humidity (RELHUMI_L1) [%].
    qsr : array_like
        Specific humidity scaling parameter [kg/kg].
    zrefht : float
        Reference height at which to compute specific humidity [m].
    zhumair : array_like
        Height of the humidity measurement [m].
    L : array_like
        Obukhov length scale [m].

    Returns
    -------
    spechum : array_like
        Specific humidity at zrefht [g/kg].

    Notes
    -----
    Equivalent to DPS section A.3.3, per J. Edson's coare35vn v3.5
    code.
    """
    qsr = qsr * 1000.0  # change units to g/kg for this calculation
    von = 0.4      # von Karman constant
    Q_air = met_spechum(tC_air, pr_air, relhum)
    spechum = Q_air + qsr / von * (np.log(zrefht/zhumair) - psit_26(zrefht/L) +
                                   psit_26(zhumair/L))

    return spechum


def water_thermal_expansion(tC_water):
    """
    Returns the thermal expansion coefficient of water. Used in both
    the warmlayer and coolskin algorithms.

    Parameters
    ----------
    tC_water : array_like
        Water temperature [deg_C].

    Returns
    -------
    Al : array_like
        Water thermal expansion coefficient. Units not documented in
        the original code.
    """

    return 2.1e-5 * (tC_water + 3.2)**0.79


def windspeed_at_refheight(wnd, usr, ut, zrefht, zwindsp, L):
    """
    Computes wind speed at an arbitrary reference height using the
    Monin-Obukhov similarity profile.

    Parameters
    ----------
    wnd : array_like
        Wind speed relative to current (RELWIND_SPD-AUX) [m/s].
    usr : array_like
        Friction velocity including gustiness [m/s].
    ut : array_like
        Effective relative wind speed including gustiness [m/s].
    zrefht : float
        Reference height at which to compute wind speed [m].
    zwindsp : array_like
        Height of the wind speed measurement [m].
    L : array_like
        Obukhov length scale [m].

    Returns
    -------
    windspeed : array_like
        Wind speed at zrefht [m/s].

    Notes
    -----
    Equivalent to DPS section A.3.3, per J. Edson's coare35vn v3.5
    code.
    """
    von = 0.4      # von Karman constant

    windspeed = wnd + usr / von / ut * wnd * (np.log(zrefht/zwindsp) -
                                              psiu_26(zrefht/L) +
                                              psiu_26(zwindsp/L))

    return windspeed


"""
#...................................................................................
#...................................................................................

    Wrapper function which calls the warmlayer and coolskin (coare35vn) routines:

        seasurface_skintemp_correct

#...................................................................................
#...................................................................................
"""


@deprecated
def seasurface_skintemp_correct(*args):
    """
    Wrapper that applies the METBK sea surface skin temperature
    correction algorithms -- warmlayer and coolskin (coare35vn) -- to
    compute the friction velocity, temperature and humidity scaling
    parameters, and other fundamental bulk parameters needed by the
    L2 BULKFLX data products. Applies both corrections by OOI default
    (JWARMFL = JCOOLFL = 1); the switch construct is retained for
    generality. Warmlayer corrections (dsea) are added to, and
    coolskin corrections (dter, dqer) are subtracted from, the bulk
    sea temperature.

    Parameters
    ----------
    args : tuple
        Positional arguments, in order:

        rain_rate : array_like
            Rain rate [mm/hr].
        timestamp : array_like
            Sample date and time [seconds since 1900-01-01].
        lon : array_like
            Longitude [deg].
        ztmpwat : array_like
            Depth of the bulk sea temperature measurement [m].
        tC_sea : array_like
            Bulk sea surface temperature [deg_C].
        wnd : array_like
            Wind speed relative to current [m/s].
        zwindsp : array_like
            Height of the wind speed measurement [m].
        tC_air : array_like
            Air temperature [deg_C].
        ztmpair : array_like
            Height of the air temperature measurement [m].
        relhum : array_like
            Relative humidity [%].
        zhumair : array_like
            Height of the air humidity measurement [m].
        pr_air : array_like
            Air pressure [mb].
        Rshort_down : array_like
            Downwelling shortwave irradiation [W/m^2].
        Rlong_down : array_like
            Downwelling longwave irradiation [W/m^2].
        lat : array_like
            Latitude [deg].
        zinvpbl : array_like
            Planetary boundary layer inversion height, default 600 m
            [m].
        jcool : array_like
            Coolskin algorithm switch (hardwired to 1 for OOI).
        jwarm : array_like
            Warmlayer algorithm switch (hardwired to 1 for OOI).

    Returns
    -------
    usr : array_like
        Friction velocity including gustiness [m/s].
    tsr : array_like
        Temperature scaling parameter [K].
    qsr : array_like
        Specific humidity scaling parameter [kg/kg].
    ut : array_like
        Effective relative wind speed including gustiness [m/s].
    dter : array_like
        Coolskin temperature depression [deg_C].
    dqer : array_like
        Coolskin humidity depression [kg/kg].
    tkt : array_like
        Coolskin thickness [m].
    L : array_like
        Obukhov length scale [m].
    zou : array_like
        Wind roughness length [m].
    zot : array_like
        Thermal roughness length [m].
    zoq : array_like
        Moisture roughness length [m].
    dt_wrm : array_like
        Warming across the entire warmlayer [deg_C].
    tk_pwp : array_like
        Warmlayer thickness [m].
    dsea : array_like
        Additive warmlayer temperature correction [deg_C].

    See Also
    --------
    warmlayer : Computes dt_wrm, tk_pwp, dsea.
    coare35vn : Computes usr, tsr, qsr, ut, dter, dqer, tkt, L, zou,
        zot, zoq.

    Notes
    -----
    When jwarm = jcool = 0, input sea temperatures are treated as
    already-measured skin temperatures and neither correction is
    applied; when jwarm = jcool = 1 (the OOI case), the input is
    treated as a bulk measurement requiring both corrections. Each
    calling wrapper on this page explicitly unpacks this function's
    output tuple with underscore placeholders, to make clear which
    corrections apply to which intermediate variable in a given data
    product's calculation.
    """
    jwarm = args[-1]    # jwarm (and jcool) are scalars
    if jwarm:
        (dt_wrm, tk_pwp, dsea) = warmlayer(*args[0:-1])  # does not pass jwarm
    else:
        # the tk_pwp parameter is often used as a divisor in warmlayer calculations to
        # compare the warmlayer depth with the depth of the bulk temperature sensor.
        # when the warmlayer code is not run, the desired results will be obtained if
        # dt_warm and dsea are set to 0 where tk_pwp is nonzero so that a divide by
        # zero error does not result. the value chosen is the default value specified
        # in the warmlayer code itself.
        (dt_wrm, tk_pwp, dsea) = (0.0, 19.0, 0.0)

    # construct tuple containing coolskin input arguments;
    # add the warmlayer temperature correction to the msrd bulk sea temp.
    coolskin_args = (args[4]+dsea,) + args[5:-1]    # does not pass jwarm
    # append results of warmlayer calculation to output,
    # as is also done in original coare35vn warmlayer matlab code.
    return coare35vn(*coolskin_args) + (dt_wrm, tk_pwp, dsea)

"""
#...................................................................................
#...................................................................................

    warmlayer

#...................................................................................
#...................................................................................
"""


@deprecated
def warmlayer(rain_rate, timestamp, lon, ztmpwat, tC_sea, wnd, zwindsp, tC_air, ztmpair, relhum,
              zhumair, pr_air, Rshort_down, Rlong_down, lat, zinvpbl, jcool):
    """
    Computes the warmlayer correction to bulk sea surface temperature,
    accounting for solar heating between the sub-surface temperature
    sensor and the air-sea interface. Refactored from
    coare35vnWarm.m.

    Parameters
    ----------
    rain_rate : array_like
        Rain rate, hourly [mm/hr].
    timestamp : array_like
        Sample date and time, hourly [seconds since 1900-01-01].
    lon : array_like
        Longitude [deg].
    ztmpwat : array_like
        Depth of the bulk sea temperature measurement [m].
    tC_sea : array_like
        Bulk sea surface temperature [deg_C].
    wnd : array_like
        Wind speed relative to current [m/s].
    zwindsp : array_like
        Height of the wind speed measurement [m].
    tC_air : array_like
        Air temperature [deg_C].
    ztmpair : array_like
        Height of the air temperature measurement [m].
    relhum : array_like
        Relative humidity [%].
    zhumair : array_like
        Height of the air humidity measurement [m].
    pr_air : array_like
        Air pressure [mb].
    Rshort_down : array_like
        Downwelling shortwave irradiation [W/m^2].
    Rlong_down : array_like
        Downwelling longwave irradiation [W/m^2].
    lat : array_like
        Latitude [deg].
    zinvpbl : array_like
        Planetary boundary layer inversion height, default 600 m [m].
    jcool : array_like
        Coolskin algorithm switch (hardwired to 1 for OOI).

    Returns
    -------
    dt_wrm : array_like
        Warming across the entire warmlayer [deg_C].
    tk_pwp : array_like
        Warmlayer thickness [m].
    dsea : array_like
        Additive warmlayer temperature correction [deg_C]; the key
        warmlayer output.

    Notes
    -----
    Only processes days that have data before a 6 AM local threshold
    (equatorial sunrise); other days return nan (see
    warmlayer_time_keys). Refactored from the original code: local
    time is no longer allowed to go negative, and a spurious
    half-hour offset present in the original local-time calculation
    was removed.
    """
    # set constants
    c2k = 273.15        # Converts degC to Kelvin
    cpw = 4000.0        # specific heat capacity of sw at T=20 degC, S=35 [J/kg/K]
    rhow = 1022.0       # density of seawater at T=20C, S=31 kg/m^3.
    cpa = 1004.67       # specific heat capacity of (dry) air [J/kg/K]

    #.. hardcoded warmlayer parameters
    rich = 0.65         # critical Richardson number

    #.. initialize warmlayer variables
    fxp = 0.5           # initial value of solar flux absorption
    max_pwp = 19.0      # maximum depth of warm layer (adjustable)
    jamset = 0          # warmlayer threshold indicator
    qcol_ac = 0.0       # accumulates heat from integral
    tau_ac = 0.0        # accumulates stress from integral

    #.. vector calculation of variables used in loop.
    rhoa = air_density(tC_air, pr_air, relhum)
    Rns = met_netsirr(Rshort_down)
    Al = water_thermal_expansion(tC_sea)   # original code does not use dter nor dsea
    grav = gravity(lat)
    ctd1 = np.sqrt(2.0 * rich * cpw / (Al * grav * rhow))         # mess-o-constants 1
    ctd2 = np.sqrt(2.0 * Al * grav / (rich * rhow)) / (cpw**1.5)  # mess-o-constants 2

    #**********************************************************
    nx = timestamp.size        # number of lines of data

    #.. initialize warmlayer products with the default values applicable to
    #.. the case where there is no warmlayer correction to the coare35vn
    #.. coolskin calculation.
    dt_wrm = np.zeros(nx)            # warming across entire warm layer deg.C
    tk_pwp = np.zeros(nx) + max_pwp  # warm layer thickness m
    dsea = np.zeros(nx)              # correction to get to interface

    # local solar time adjustment is a function of longitude:
    #..  360 degrees = 24 *3600 seconds,
    #..  so each degree is worth 240 seconds of time.
    # the OOI timestamp is seconds since midnight 01-jan-1900; therefore
    # local time will still be positive for the case of lon = -180deg.
    local_date_time = timestamp + lon * 240.0

    # and calculate all the delta times [sec] for the integrals' abscissae.
    #.. prepend a zero to line up delta_time with iteration number.
    #.. values at newday records are not used in the calculations.
    delta_time = np.hstack((0.0, np.diff(local_date_time)))

    # determine:
    #    idx_warm: indices of data for days that have data before 6AM
    #    newday_bool: boolean mask, true for the first record of each day.
    #    nanmask: boolean, true for records of days that do not have data before 6AM.
    #             the output at these indices will be changed from initialized to nan.
    idx_warm, newday_bool, nanmask = warmlayer_time_keys(local_date_time)

    #.. the original code has been changed to show the explicit dependence
    #.. of the variables upon iteration count (data record number).
    for ii in idx_warm:   # step through each timepoint

        #.. warmlayer values for the following case are just the initialized
        #.. values. so, instead of using if-then-else, simplify indentation
        #.. by using 'if' only, reset variables, and jump to next iteration.
        if newday_bool[ii]:  # re-zero when starting a new day
            # dt_wrm[ii] = 0.0;
            # tk_pwp[ii] = max_pwp;
            # dsea[ii]   = 0.0;
            jamset = 0
            fxp = 0.5
            tau_ac = 0.0
            qcol_ac = 0.0
            continue  # go to next time (data) record
            # end midnight reset

        #*****  dependent variables for the [ii]th warm layer calculation
        #*****  of dsea are fluxes, coolskin correction dter, and dsea itself,
        #*****  which are derived from the previous ([ii-1]th) data record.
        #
        # because of the dependence on the previous value of dsea, this calculation
        # cannot be vectorized.
        tsea_corr = tC_sea[ii-1] + dsea[ii-1]

        # slicing 1D arrays with [ii-1:ii] returns a 1-element nd.array variable which
        # can be indexed, whereas slicing with [ii-1] returns a variable which cannot be
        # indexed. [ii-1:ii] slicing is used so that coare35vn can be run with both
        # 'scalar' and 'vector' input.
        args = (tsea_corr, wnd[ii-1:ii], zwindsp[ii-1:ii], tC_air[ii-1:ii], ztmpair[ii-1:ii], relhum[ii-1:ii], zhumair[ii-1:ii],
                pr_air[ii-1:ii], Rshort_down[ii-1:ii], Rlong_down[ii-1:ii], lat[ii-1:ii], zinvpbl[ii-1:ii],
                jcool)

        (usr, tsr, qsr, ut, dter, dqer, _, _, _, _, _) = coare35vn(*args)

        # in the original matlab code, Le was calculated inside of the coare35vn
        # subroutine, which was called using tC_sea+dsea for seawater temperature:
        Le = latent_heat_vaporization_pure_water(tsea_corr)
        tau_old = rhoa[ii-1] * usr * usr * wnd[ii-1] / ut  # stress
        hs_old = -rhoa[ii-1] * cpa * usr * tsr              # sensible heat flux
        hl_old = -rhoa[ii-1] * Le * usr * qsr                      # latent heat flux

        # note:
        #     the source v3.5 matlab code is followed here: it does not use dsea
        #     in the Rnl expression used in the warmlayer calculation, although
        #     dsea is used in the expression for RF_old.
        Rnl = net_longwave_up(tC_sea[ii]-dter, Rlong_down[ii])
        RF_old = rain_heat_flux(rain_rate[ii-1], tC_sea[ii-1]+dsea[ii-1], tC_air[ii-1],
                                relhum[ii-1], pr_air[ii-1])

        #********************************************************
        #****  Compute warm layer correction *******************
        #********************************************************
        qr_out = Rnl + hs_old + hl_old + RF_old  # total cooling at surface
        q_pwp = fxp * Rns[ii] - qr_out          # tot heat abs in warm layer

        # calculate dt_wrm and tk_pwp for this iteration.
        if q_pwp >= 50.0 or jamset == 1:         # Check for threshold
            jamset = 1			         # indicates threshold crossed
            tau_ac = tau_ac + np.maximum(.002, tau_old) * delta_time[ii]  # momentum integral

            # check threshold for warm layer existence
            if qcol_ac + q_pwp * delta_time[ii] > 0.0:
                #******************************************
                # Compute the absorption profile
                #******************************************
                #.. tk_pwp can iteratively change value in the following loop,
                #.. requiring the creation of the variable tkpwp.
                tkpwp = tk_pwp[ii-1]
                for i in range(5):               # loop 5 times for fxp
                    fxp = 1.0 - (0.28 * 0.014 * (1.0 - np.exp(-tkpwp / 0.014)) +
                                 0.27 * 0.357 * (1.0 - np.exp(-tkpwp / 0.357)) +
                                 0.45 * 12.82 * (1.0 - np.exp(-tkpwp / 12.82))) / tkpwp
                    qjoule = (fxp * Rns[ii] - qr_out) * delta_time[ii]
                    if qcol_ac + qjoule > 0.0:   # Compute warm-layer depth
                        tkpwp = np.minimum(max_pwp,
                                           ctd1[ii] * tau_ac / np.sqrt(qcol_ac + qjoule))
                tk_pwp[ii] = tkpwp
            else:                                # warm layer wiped out
                fxp = 0.75
                tk_pwp[ii] = max_pwp
                qjoule = (fxp * Rns[ii] - qr_out) * delta_time[ii]

            qcol_ac = qcol_ac + qjoule           # heat integral

            #*******  compute dt_warm  ******
            if qcol_ac > 0.0:
                dt_wrm[ii] = ctd2[ii] * (qcol_ac)**1.5 / tau_ac
            else:
                dt_wrm[ii] = 0.0

        else:   # propagate dt_wrm and tk_pwp values
            dt_wrm[ii] = dt_wrm[ii-1]
            tk_pwp[ii] = tk_pwp[ii-1]

        # Compute warm layer correction dsea
        if tk_pwp[ii] < ztmpwat[ii]:
            dsea[ii] = dt_wrm[ii]
        else:
            dsea[ii] = dt_wrm[ii] * ztmpwat[ii] / tk_pwp[ii]

    # for all days that did not begin before 6AM, return NaNs
    dt_wrm[nanmask] = np.nan
    tk_pwp[nanmask] = np.nan
    dsea[nanmask] = np.nan

    return dt_wrm, tk_pwp, dsea

"""
#...................................................................................
#...................................................................................

    coare35vn coolskin code and the subroutines unique to it.
        coare35vn also calls subroutines located elsewhere in this module.

    coare35vn

        charnock_wind
        coolskin_parameters
        effective_relwind
        obukhov_for_init
        obukhov_length_scale
        roughness_lengths
        roughness_lengths_for_init
        scaling_parameters
#...................................................................................
#...................................................................................
"""


@deprecated
def coare35vn(tC_sea, wnd, zwindsp, tC_air, ztmpair, relhum, zhumair, pr_air,
              Rshort_down, Rlong_down, lat, zinvpbl, jcool):
    """
    Iteratively computes the fundamental bulk parameters, with
    coolskin correction, from which METBK air-sea flux data products
    are calculated. Transliterated from version 3.5 of coare35vn.m.
    Unlike the original fortran and MATLAB versions, this function
    does not directly compute data products; see the wrapper
    functions on this docs page for the individual L2 products.

    Parameters
    ----------
    tC_sea : array_like
        Bulk sea surface temperature [deg_C].
    wnd : array_like
        Wind speed relative to current [m/s].
    zwindsp : array_like
        Height of the wind speed measurement [m].
    tC_air : array_like
        Air temperature [deg_C].
    ztmpair : array_like
        Height of the air temperature measurement [m].
    relhum : array_like
        Relative humidity [%].
    zhumair : array_like
        Height of the air humidity measurement [m].
    pr_air : array_like
        Air pressure [mb].
    Rshort_down : array_like
        Downwelling shortwave irradiation [W/m^2].
    Rlong_down : array_like
        Downwelling longwave irradiation [W/m^2].
    lat : array_like
        Latitude [deg].
    zinvpbl : array_like
        Planetary boundary layer inversion height, default 600 m [m].
    jcool : array_like
        Coolskin algorithm switch (hardwired to 1 for OOI).

    Returns
    -------
    usr : array_like
        Friction velocity including gustiness [m/s].
    tsr : array_like
        Temperature scaling parameter [K].
    qsr : array_like
        Specific humidity scaling parameter [kg/kg].
    ut : array_like
        Effective relative wind speed including gustiness [m/s]; not
        an output of the original code.
    dter : array_like
        Coolskin temperature depression [deg_C].
    dqer : array_like
        Coolskin humidity depression [kg/kg].
    tkt : array_like
        Coolskin thickness [m].
    L : array_like
        Obukhov length scale [m].
    zou : array_like
        Wind roughness length [m].
    zot : array_like
        Thermal roughness length [m].
    zoq : array_like
        Moisture roughness length [m].

    See Also
    --------
    charnock_wind, coolskin_parameters, effective_relwind,
    obukhov_for_init, obukhov_length_scale, roughness_lengths,
    roughness_lengths_for_init, scaling_parameters : Per-iteration
        subroutines called by the loop in this function.

    Notes
    -----
    Iterates 6 times (hardwired). The per-iteration calculation was
    split into named subroutines to clarify which intermediate
    variables are recomputed each pass; the original code carried
    little documentation of this structure. dter and dqer are forced
    to 0 when jcool=0, so downstream callers need no jcool
    multiplicative factor.
    """
    # convert relative humidity to specific humidity [kg/kg]
    Qsea = sea_spechum(tC_sea, pr_air) / 1000.0          # surface water specific humidity
    Qair = met_spechum(tC_air, pr_air, relhum) / 1000.0  # specific humidity of air

    #***********  set constants **********************************************
    Beta = 1.2
    von = 0.4     # von karman constant
    fdg = 1.00    # Turbulent Prandtl number
    c2k = 273.15  # degC to kelvin
    grav = gravity(lat)

    #***********  air constants **********************************************
    Rgas = 287.05
    Le = latent_heat_vaporization_pure_water(tC_sea)
    cpa = 1004.67   # specific heat capacity of dry air, J/kg/K
    rhoa = air_density(tC_air, pr_air, relhum)
    visa = 1.326e-5 * (1.0 + 6.542e-3 * tC_air + 8.301e-6 * tC_air**2 -
                       4.84e-9 * tC_air**3)
    lapse = grav / cpa

    #***********  cool skin constants  ***************************************
    Al = water_thermal_expansion(tC_sea)
    be = 0.026     # salinity expansion coeff.
    cpw = 4000.0   # specific heat of seawater at T=20C, S=35 J/kg/K.
    rhow = 1022.0  # density of seawater at T=20C, S=31 kg/m^3.
    visw = 1.0e-6
    tcw = 0.6
    bigc = 16.0 * grav * cpw * (rhow * visw)**3 / (tcw**2 * rhoa**2)

    #.. derived variables unchanged by loop
    dq = Qsea - Qair
    dt = tC_sea - tC_air - lapse * ztmpair
    tK_air = tC_air + c2k
    tv = tK_air * (1.0 + 0.61 * Qair)   # virtual temperature
    ug = 0.5

    #*********** initialization **********************************************

    #.. coolskin parameters (changed inside loop and used after loop)
    #
    #.. if jcool is set to 0, then all dter and dqer values in this code will
    #.. also be 0 (and therefore not require a multiplicative factor of jcool
    #.. as in the original matlab code) *except* for when they are directly
    #.. calculated inside the loop.
    dter = 0.3 * jcool
    wetc = 0.622 * Le * Qsea / (Rgas * (tC_sea + c2k)**2)
    dqer = dter * wetc
    tkt = 0.001 + np.zeros(wnd.size)
    ut = np.sqrt(wnd**2 + ug**2)
    # for initialization of usr, tsr, qsr, charnC
    # original code does use (10.0/1e-4)
    u10 = ut * np.log(10.0/1.0e-4) / np.log(zwindsp/1.0e-4)

    #.. scaling parameters usr,qsr,tsr
    zo10, zot10 = roughness_lengths_for_init(0.035*u10, grav, visa, von)
    #.. k50 is used inside the loop to save loop variables for ii=0
    L10, k50 = obukhov_for_init(von, grav, tK_air, dt, dter, dq, ut,
                                zwindsp, ztmpair, zo10, zot10, zinvpbl, Beta)
    usr, qsr, tsr = scaling_parameters(dter, dqer, von, fdg, zwindsp,
                                       zhumair, ztmpair, zo10, zot10, zot10,
                                       L10, ut, dq, dt)

    #***********  net radiation fluxes ***************************************
    Rns = met_netsirr(Rshort_down)                  # net shortwave radiation DOWN
    Rnl = net_longwave_up(tC_sea-dter, Rlong_down)  # net longwave radiation UP

    #**********************************************************
    #  The following gives the new formulation for the
    #  Charnock variable in COARE 3.5
    #**********************************************************

    charnC = charnock_wind(u10)

    nits = 6  # hardwired number of iterations

    #**************  bulk loop ***********************************************
    for ii in range(nits):

        L = obukhov_length_scale(von, grav, tK_air, Qair, usr, tsr, qsr)
        zou, zoq, zot = roughness_lengths(charnC, usr, grav, visa)

        usr, qsr, tsr = scaling_parameters(dter, dqer, von, fdg, zwindsp,
                                           zhumair, ztmpair, zou, zoq, zot, L, ut,
                                           dq, dt)

        dter, dqer, tkt = coolskin_parameters(usr, qsr, tsr, Rnl, Rns, rhoa, cpa,
                                              Le, tkt, Al, be, cpw, visw, rhow,
                                              bigc, tcw, tC_sea, Qsea, pr_air)

        # these coolskin parameters must be reset to 0 if coolskin is off, so:
        dter = dter * jcool
        dqer = dqer * jcool

        Rnl = net_longwave_up(tC_sea-dter, Rlong_down)

        ut = effective_relwind(tsr, tK_air, qsr, grav, tv, usr, Beta, zinvpbl, wnd)

        #.. update charnock variable.
        u10N_for_charnC = usr / von / ut * wnd * np.log(10.0/zou)
        charnC = charnock_wind(u10N_for_charnC)

        # for stable cases designated by k50 save the results from the
        # first iteration as the algorithm output. this construction also
        # works as desired if k50 is empty.
        if ii == 0:
            stable_cases = (usr[k50], tsr[k50], qsr[k50], ut[k50],
                            dter[k50], dqer[k50], tkt[k50], L[k50],
                            zou[k50], zot[k50], zoq[k50])

    # loop is finished:
    # insert first iteration solution for stable cases.
    (usr[k50], tsr[k50], qsr[k50], ut[k50], dter[k50], dqer[k50],
        tkt[k50], L[k50], zou[k50], zot[k50], zoq[k50]) = stable_cases

    # the whoi v3.0 and jbe v3.5 qsr units differ.
    #    whoi (DPS): [kg/kg] (same as [g/g])
    #    jbe - output units are [g/kg]
    #
    #    because the calculations of fluxes are natively done with qsr in units
    #    of [kg/kg], this is what will be used.
    #
    #    HOWEVER - when calculating specific humidity at a reference height,
    #              qsr must have units of g/kg:
    #              qsr = qsr * 1000    # changes units from kg/kg to g/kg

    # note: original DPS code returns nonzero values for dter and dqer even
    # when jcool = 0; the edson v3.5 code returns a nonzero dter value, but
    # a zero value for dqer.
    #
    # this code gives dter = dqer = 0 when jcool = 0. this also obviates the
    # need for multiplicative factors of jcool when this output is used in
    # subsequent routines.

    return (usr, tsr, qsr, ut, dter, dqer, tkt, L, zou, zot, zoq)


#-------------------------------------------------------------------------
#---------------- subroutines unique to coare35vn-------------------------
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
def charnock_wind(u10N):
    """
    Computes the Charnock coefficient as a function of the 10 m
    neutral wind speed, following the COARE 3.5 wind-speed
    parameterization.

    Parameters
    ----------
    u10N : array_like
        10 m neutral wind speed [m/s].

    Returns
    -------
    charnC : array_like
        Charnock coefficient, capped at its value for u10N = 19 m/s.
    """
    umax = 19
    a1 = 0.0017
    a2 = -0.0050
    charnC = a1 * u10N + a2
    # trap out nans after nans have propagated into charnC
    u10N[np.isnan(u10N)] = -1.0  # anything less than umax
    charnC[u10N > umax] = a1 * umax + a2
    return charnC


#-------------------------------------------------------------------------
def coolskin_parameters(usr, qsr, tsr, Rnl, Rns, rhoa, cpa, Le, tkt, Al, be,
                        cpw, visw, rhow, bigc, tcw, tC_sea, Qsea, pr_air):
    """
    Computes the coolskin temperature depression, humidity depression,
    and layer thickness for one iteration of the coare35vn coolskin
    loop.

    Parameters
    ----------
    usr : array_like
        Friction velocity including gustiness [m/s].
    qsr : array_like
        Specific humidity scaling parameter [kg/kg].
    tsr : array_like
        Temperature scaling parameter [K].
    Rnl : array_like
        Net upward longwave radiation [W/m^2].
    Rns : array_like
        Net downward shortwave radiation [W/m^2].
    rhoa : array_like
        Air density [kg/m^3].
    cpa : float
        Specific heat capacity of dry air [J/kg/K].
    Le : array_like
        Latent heat of vaporization of seawater [J/kg].
    tkt : array_like
        Coolskin thickness from the previous iteration [m].
    Al : array_like
        Water thermal expansion coefficient.
    be : float
        Salinity expansion coefficient.
    cpw : float
        Specific heat capacity of seawater [J/kg/K].
    visw : float
        Kinematic viscosity of seawater [m^2/s].
    rhow : float
        Density of seawater [kg/m^3].
    bigc : array_like
        Coolskin scaling constant.
    tcw : float
        Thermal conductivity of seawater.
    tC_sea : array_like
        Bulk sea surface temperature [deg_C].
    Qsea : array_like
        Sea surface specific humidity [kg/kg].
    pr_air : array_like
        Air pressure [mb].

    Returns
    -------
    dter : array_like
        Coolskin temperature depression [deg_C].
    dqer : array_like
        Coolskin humidity depression [kg/kg].
    tkt : array_like
        Coolskin thickness [m].
    """
    #.. dter: coolskin temperature depression
    #.. dqer: coolskin humidity depression
    #.. tkt: coolskin thickness
    N = usr.shape[0]
    hsb = -rhoa * cpa * usr * tsr
    hlb = -rhoa * Le * usr * qsr
    qout = Rnl + hsb + hlb
    dels = Rns * (0.065 + 11.0 * tkt - 6.6e-5 / tkt *
                  (1.0 - np.exp(-tkt / 8.0e-4)))
    qcol = qout - dels
    alq = Al * qcol + be * hlb * cpw / Le
    xlamx = 6.0 * np.ones(N)
    tkt = np.minimum(0.01, xlamx * visw / (np.sqrt(rhoa / rhow) * usr))
    # trap out nans; tkt already has nans where we want them
    alq[np.isnan(alq)] = -1.0
    k = np.where(alq > 0)
    xlamx[k] = 6.0 / (1.0 + (bigc[k] * alq[k] / usr[k]**4)**0.75)**0.333
    tkt[k] = xlamx[k] * visw / (np.sqrt(rhoa[k] / rhow) * usr[k])
    dter = qcol * tkt / tcw
    # formerly, dqer = wetc * dter
    dqer = Qsea - sea_spechum(tC_sea - dter, pr_air) / 1000.0
    return dter, dqer, tkt


#-------------------------------------------------------------------------
def effective_relwind(tsr, tK_air, qsr, grav, tv, usr, Beta, zinvpbl, wnd):
    """
    Computes the effective relative wind speed, adding a
    gustiness contribution to the mean relative wind speed under
    convective (unstable) conditions.

    Parameters
    ----------
    tsr : array_like
        Temperature scaling parameter [K].
    tK_air : array_like
        Air temperature [K].
    qsr : array_like
        Specific humidity scaling parameter [kg/kg].
    grav : array_like
        Acceleration due to gravity [m/s/s].
    tv : array_like
        Virtual air temperature [K].
    usr : array_like
        Friction velocity including gustiness [m/s].
    Beta : float
        Gustiness coefficient.
    zinvpbl : array_like
        Planetary boundary layer inversion height [m].
    wnd : array_like
        Wind speed relative to current [m/s].

    Returns
    -------
    ut : array_like
        Effective relative wind speed including gustiness [m/s].
    """
    N = wnd.shape[0]
    tvsr = tsr + 0.61 * tK_air * qsr
    Bf = -grav / tv * usr * tvsr
    ug = 0.2 * np.ones(N)
    nanmask = np.isnan(Bf)
    Bf[nanmask] = -1.0
    ug[Bf > 0] = np.maximum(0.2,
                            Beta * (Bf[Bf > 0] * zinvpbl[Bf > 0]) ** 0.333)
    ug[nanmask] = np.nan
    ut = np.sqrt(wnd**2 + ug**2)
    return ut


#-------------------------------------------------------------------------
def obukhov_for_init(von, grav, tK_air, dt, dter, dq, ut, zwindsp,
                     ztmpair, zo10, zot10, zinvpbl, Beta):
    """
    Computes an initial estimate of the Obukhov length scale for the
    coare35vn iteration loop, and identifies stable-case indices with
    a Monin-Obukhov length very thin relative to zwindsp.

    Parameters
    ----------
    von : float
        Von Karman constant.
    grav : array_like
        Acceleration due to gravity [m/s/s].
    tK_air : array_like
        Air temperature [K].
    dt : array_like
        Sea-air temperature difference [deg_C].
    dter : array_like
        Coolskin temperature depression [deg_C].
    dq : array_like
        Sea-air specific humidity difference [kg/kg].
    ut : array_like
        Effective relative wind speed including gustiness [m/s].
    zwindsp : array_like
        Height of the wind speed measurement [m].
    ztmpair : array_like
        Height of the air temperature measurement [m].
    zo10 : array_like
        Wind roughness length estimate at 10 m [m].
    zot10 : array_like
        Thermal roughness length estimate at 10 m [m].
    zinvpbl : array_like
        Planetary boundary layer inversion height [m].
    Beta : float
        Gustiness coefficient.

    Returns
    -------
    L_init : array_like
        Initial Obukhov length scale estimate [m].
    k50 : array_like
        Indices where the stability parameter exceeds 50 (stable
        cases with a very thin Monin-Obukhov length).
    """
    #.. calculates an obukhov length scale for loop variable initializations.
    #.. also finds indices (k50) of zetu stability values with very thin
    #..    Monin-Obukhov lengths relative to zwindsp.
    Ribu = -grav * zwindsp / tK_air * (dt - dter + 0.61 * tK_air * dq) / ut**2
    #.. CC calculation is left unchanged from earlier code.
    Cd = (von / np.log(zwindsp/zo10))**2
    Ct = von / np.log(ztmpair/zot10)
    CC = von * Ct / Cd
    zetu = CC * Ribu * (1.0 + 27.0 / 9.0 * Ribu / CC)
    # trap out nans
    nanmask = np.isnan(zetu)
    zetu[nanmask] = 0.0
    k50 = np.where(zetu > 50)  # stable with very thin M-O length relative to zwindsp
    # restore nans to zetu
    zetu[nanmask] = np.nan

    Ribcu = -zwindsp / zinvpbl / 0.004 / Beta**3
    # trap out nans; Ribu>0 values are not further processed, so:
    Ribu[np.isnan(Ribu)] = 1.0
    k = Ribu < 0
    zetu[k] = CC[k] * Ribu[k] / (1.0 + Ribu[k] / Ribcu[k])
    L_init = zwindsp / zetu
    return L_init, k50


#-------------------------------------------------------------------------
def obukhov_length_scale(von, grav, tK_air, Qair, usr, tsr, qsr):
    """
    Computes the Obukhov length scale from average air temperature
    and humidity, and the friction velocity, temperature, and
    humidity scaling parameters.

    Parameters
    ----------
    von : float
        Von Karman constant.
    grav : array_like
        Acceleration due to gravity [m/s/s].
    tK_air : array_like
        Air temperature [K].
    Qair : array_like
        Air specific humidity [kg/kg].
    usr : array_like
        Friction velocity including gustiness [m/s].
    tsr : array_like
        Temperature scaling parameter [K].
    qsr : array_like
        Specific humidity scaling parameter [kg/kg].

    Returns
    -------
    L : array_like
        Obukhov length scale [m].

    Notes
    -----
    Per Liu et al. (1979); DPS documentation cites this as the source
    formulation.
    """
    tv = tK_air * (1.0 + 0.61 * Qair)
    tvsr = tsr * (1.0 + 0.61 * Qair) + 0.61 * tK_air * qsr
    # add tvsr adjustment to avoid program failure when tvsr very small
    #.. (and evidently, negative).
    #.. note also that tvsr=0 is not trapped out in original code.
    nanmask = np.logical_or(np.isnan(tvsr), tvsr == 0)
    tvsr[nanmask] = 100.0
    mask = np.abs(tvsr) < 1.e-3
    tvsr[mask] = np.abs(tvsr[mask])
    tvsr[nanmask] = np.nan
    L = tv * usr * usr / (grav * von * tvsr)

    return L


#-------------------------------------------------------------------------
def roughness_lengths(charn, usr, grav, visa):
    """
    Computes the wind, moisture, and thermal roughness lengths for
    one iteration of the coare35vn loop.

    Parameters
    ----------
    charn : array_like
        Charnock coefficient.
    usr : array_like
        Friction velocity including gustiness [m/s].
    grav : array_like
        Acceleration due to gravity [m/s/s].
    visa : array_like
        Kinematic viscosity of air [m^2/s].

    Returns
    -------
    zo : array_like
        Wind roughness length [m].
    zoq : array_like
        Moisture roughness length [m].
    zot : array_like
        Thermal roughness length [m].

    Notes
    -----
    The thermal and moisture roughness lengths are set equal, chosen
    to give Stanton and Dalton numbers that closely approximate
    COARE 3.0.
    """
    zo = charn * usr**2 / grav + 0.11 * visa / usr  # surface roughness
    rr = zo * usr / visa
    # These thermal roughness lengths give Stanton and
    # Dalton numbers that closely approximate COARE 3.0
    zoq = np.minimum(1.6e-4, 5.8e-5 / rr**0.72)
    zot = zoq
    return zo, zoq, zot


#-------------------------------------------------------------------------
def roughness_lengths_for_init(usr, grav, visa, von):
    """
    Computes initial wind and thermal roughness length estimates at
    10 m, used to seed the coare35vn iteration loop.

    Parameters
    ----------
    usr : array_like
        Friction velocity estimate [m/s].
    grav : array_like
        Acceleration due to gravity [m/s/s].
    visa : array_like
        Kinematic viscosity of air [m^2/s].
    von : float
        Von Karman constant.

    Returns
    -------
    zo10 : array_like
        Wind roughness length estimate at 10 m [m].
    zot10 : array_like
        Thermal roughness length estimate at 10 m [m].

    Notes
    -----
    Uses a fixed stand-in Charnock coefficient of 0.011 for this
    initial estimate only.
    """
    charn_standin = 0.011
    zo10 = charn_standin * usr**2 / grav + 0.11 * visa / usr
    #.. the following code is unchanged from the original version,
    #.. in case there's ever a question of where zot10 came from.
    Cd10 = (von / np.log(10./zo10))**2
    Ch10 = 0.00115
    Ct10 = Ch10 / np.sqrt(Cd10)
    zot10 = 10. / np.exp(von/Ct10)
    return zo10, zot10


#-------------------------------------------------------------------------
def scaling_parameters(dter, dqer, von, fdg, zwindsp, zhumair,
                       ztmpair, zou, zoq, zot, L, ut, dq, dt):
    """
    Computes the friction velocity, specific humidity, and
    temperature scaling parameters for one iteration of the
    coare35vn loop.

    Parameters
    ----------
    dter : array_like
        Coolskin temperature depression [deg_C].
    dqer : array_like
        Coolskin humidity depression [kg/kg].
    von : float
        Von Karman constant.
    fdg : float
        Turbulent Prandtl number.
    zwindsp : array_like
        Height of the wind speed measurement [m].
    zhumair : array_like
        Height of the humidity measurement [m].
    ztmpair : array_like
        Height of the air temperature measurement [m].
    zou : array_like
        Wind roughness length [m].
    zoq : array_like
        Moisture roughness length [m].
    zot : array_like
        Thermal roughness length [m].
    L : array_like
        Obukhov length scale [m].
    ut : array_like
        Effective relative wind speed including gustiness [m/s].
    dq : array_like
        Sea-air specific humidity difference [kg/kg].
    dt : array_like
        Sea-air temperature difference [deg_C].

    Returns
    -------
    usr : array_like
        Friction velocity including gustiness [m/s].
    qsr : array_like
        Specific humidity scaling parameter [kg/kg].
    tsr : array_like
        Temperature scaling parameter [K].
    """
    cdhf = von / (np.log(zwindsp / zou) - psiu_26(zwindsp / L))
    cqhf = von * fdg / (np.log(zhumair / zoq) - psit_26(zhumair / L))
    cthf = von * fdg / (np.log(ztmpair / zot) - psit_26(ztmpair / L))
    usr = ut * cdhf
    qsr = -(dq - dqer) * cqhf
    tsr = -(dt - dter) * cthf
    return usr, qsr, tsr


"""
#...................................................................................
#...................................................................................
    Data conditioning and averaging routines

        vet_velptmn_data
        condition_data
        make_hourly_data
        warmlayer_time_keys
#...................................................................................
#...................................................................................
"""


def vet_velptmn_data(vle, vln, use_velptmn):
    """
    Replaces suspect VELPTMN current values with nan where the data
    quality flag indicates bad data.

    Parameters
    ----------
    vle : array_like
        Eastward surface current speed [m/s].
    vln : array_like
        Northward surface current speed [m/s].
    use_velptmn : array_like
        Time-vectorized data quality flag: 0 for bad current data, 1
        for good current data.

    Returns
    -------
    vle_out : array_like
        Eastward surface current speed, suspect values set to nan
        [m/s].
    vln_out : array_like
        Northward surface current speed, suspect values set to nan
        [m/s].

    Notes
    -----
    Not for use with met_relwind_speed, which replaces suspect values
    with 0 instead of nan.
    """
    # expand use_velptmn_with_metbk if it is called as a scalar
    if np.atleast_1d(use_velptmn).shape[0] == 1:
        use_velptmn = np.tile(use_velptmn, vle.shape[0])

    # to prevent what turned out to be very much unanticipated "call-by-reference"-ish
    # ramifications in the unit tests
    vle_out = np.copy(vle)
    vln_out = np.copy(vln)
    # replace aliased current values with nans.
    nanmask = use_velptmn == 0
    vle_out[nanmask] = np.nan
    vln_out[nanmask] = np.nan

    return vle_out, vln_out


def condition_data(*args):
    """
    Conditions the input argument list for the warmlayer/coolskin
    algorithm: coerces all arguments to at least 1D arrays, expands
    scalar-default arguments (including sensor heights ztmpwat,
    zwindsp, ztmpair, zhumair) to match the length of the
    time-vectorized arguments, and reduces the jcool/jwarm switches
    to single-element arrays.

    Parameters
    ----------
    *args : array_like
        Argument list of input data.

    Returns
    -------
    args_out : list of array_like
        Argument list of conditioned output data.

    Notes
    -----
    Called at the front end of every function that requires the
    warmlayer/coolskin algorithm. May also be called with fewer than
    the full 18-argument set (e.g. for met_rainrte), in which case
    only array-shape coercion is performed.
    """
    # to enable modification of the input arguments in place,
    # make sure that the args are in a list
    args = list(args)

    number_of_bulk_vars = 18

    # zinvpbl [15] must always be expanded
    idx_of_args_to_expand = [0, 3, 6, 8, 10, 11, 12, 13, 14, 15]

    nargs = len(args)
    for ii in range(nargs):
        args[ii] = np.atleast_1d(args[ii])

    # for rainrte and testing
    if nargs < number_of_bulk_vars:
        return args

    # condition the jcool and jwarm switches.
        # (1) these should be type integer, so there is no need to trap out Nans.
        # (2) these switches should not be time-vectorized in the code. therefore,
        #     the code is made compatible with time-vectorized inputs of these
        #     switches by using only the 1st value; and, if it is not zero, change
        #     it to 1. This will also trap out -99999999 system fillvalues.
    for ii in [16, 17]:
        if args[ii][0] != 0:
            args[ii][0] = 1
        args[ii] = args[ii][0:1]  # return 1 element as an ndarray

    # expand if only one item in argument
    n_records = args[1].size
    for ii in idx_of_args_to_expand:
        if args[ii].size == 1:
            args[ii] = np.zeros(n_records) + args[ii]

    return args


def make_hourly_data(*args):
    """
    Computes hourly averages of the variables passed in the argument
    list, binned by timestamp. Works for sporadically spaced data and
    for data with time gaps; no records are produced for missing time
    bins.

    Parameters
    ----------
    *args : array_like
        Argument list of each-minute data. All arguments must be 1D
        arrays of the same length. Timestamp (seconds) must be the
        second element unless it is the only argument. If the
        argument count is 18 (the warmlayer/coolskin input set),
        constant and switch arguments are not averaged.

    Returns
    -------
    args_out : list of array_like
        Argument list of hourly-averaged data.

    Notes
    -----
    Each hourly timestamp marks the midpoint of its bin; the first
    hourly timestamp is a half hour after the first each-minute
    datum, so hourly timestamps do not in general fall on the top of
    the hour. Uses np.bincount, weighted, to sum and average values
    per bin.
    """
    args = list(args)

    number_of_bulk_vars = 18

    # prep all variables for rainrte (nargs<18).
    # for nargs >= 18, skip arguments that are constants (except for zinvpbl)
    # and switches.
    idx_to_skip = [17, 16]

    nargs = len(args)
    idx = list(range(nargs))
    if nargs >= number_of_bulk_vars:
        for ii in idx_to_skip:
            del idx[ii]

    # timestamps must be the 2nd variable in the input argument list,
    # unless there is only 1 variable.
    index_timedata = np.sign(nargs-1)
    time_sec = args[index_timedata]
    time_elapsed_hr = (time_sec - time_sec[0])/3600.0

    # assign each timestamp a bin number index based on its elapsed time in hrs.
    bin_number = np.floor(time_elapsed_hr).astype(int)
    # the number of elements in each hourly bin is given by
    bin_count = np.bincount(bin_number).astype(float)
    # create a logical mask of non-zero bin_count values
    mask = (bin_count != 0)
    # and keep only the non-zero values for calculating the average
    bin_count_no_zeros = bin_count[mask]

    # average the values in each hourly bin for each input variable;
    # the np.bincount function only works on 1D arrays, so
    for ivar in idx:
        # sum the values in each hourly bin for the ivar[th] variable
        args[ivar] = np.bincount(bin_number, args[ivar])
        # discard trivial bin sums of 0 where there were no bin elements
        args[ivar] = args[ivar][mask]
        # divide the bin sums by the number of elements in each bin
        args[ivar] = args[ivar] / bin_count_no_zeros

    # hourly timestamp calculation:
    #     note that the midpoint of the data interval is used, not the timestamp
    #     of the first nor last point.
    #
    #     use the midpoint of the bins as the timestamp, instead of the average
    #     of the timestamps within the bin as calculated in the above loop; this
    #     would give significantly different values only if there are missing data.
    bin_time_sec = time_sec[0] + 1800.0 + 3600.0 * np.array(list(range(len(mask))))
    # delete bins with no entries as before
    bin_time_sec = bin_time_sec[mask]

    args[index_timedata] = bin_time_sec

    #for ii in idx:
    #    print args[ii]

    return args


def warmlayer_time_keys(localdate):
    """
    Computes indices and flags used by the warmlayer routine to
    identify which days of local-time data qualify for warmlayer
    processing.

    Parameters
    ----------
    localdate : array_like
        Local (not UTC) date and time [seconds since 1900-01-01].

    Returns
    -------
    idx_warm : array_like
        Indices of data records to be processed by the warmlayer
        routine; these are records for days that have data before a
        6 AM local threshold.
    newday : array_like of bool
        True for the first record of each day, False otherwise.
    nanmask : array_like of bool
        True for indices of data records not to be processed by the
        warmlayer routine (days without data before the threshold).

    Notes
    -----
    The 6 AM (equatorial sunrise) threshold matches the original
    fortran and MATLAB code. Applies a quality-control check the
    original code lacked: a day is only warmlayer-processed if it has
    data before the threshold, preventing e.g. a data gap that
    straddles midnight from being misattributed to the wrong day.
    """
    warmlayer_threshold_OOI = 21600.0  # equatorial sunrise: 6AM

    # in the original code, no checks were included to trap out the kind
    # of situation in which there are data for a given day from local
    # times 0500-1200 immediately followed by data for the following day
    # from 1300-1800.
    #
    # finding the start of each day when the timestamps have units of seconds
    # since 01-jan-1900 is straightforward:
    newday = np.diff(np.floor(localdate/86400.0)) > 0
    # prepend a True value to get the index count correct and to start data at a new day.
    newday = np.hstack((True, newday))
    # find the indices of the start of newdays
    idx_nd = np.nonzero(newday)[0]
    # append a bracketing end index for the for loop to follow
    idx_nd = np.hstack((idx_nd, newday.size))

    # the warmlayer routine is to be run only on days which start earlier than threshold.
    time_of_day = np.mod(localdate, 86400.0)
    earlier_than_threshold = time_of_day <= warmlayer_threshold_OOI

    # initialize warmmask to all False (do not process any data records with warmlayer)
    warmmask = np.zeros(time_of_day.size, dtype=bool)
    # every warmmask value will be overwritten by the loop.
    #..  each iteration processes one day's worth of data records whose indices are
    #..     [idx_nd[ii]:idx_nd[ii+1]] (python does not use last index in its range).
    #..  the values for that day's records are determined by whether the first local
    #..     time for that day is earlier than the threshold (T) or not (F).
    for ii in range(idx_nd.size-1):
        warmmask[idx_nd[ii]:idx_nd[ii+1]] = earlier_than_threshold[idx_nd[ii]]

    idx_warm = np.nonzero(warmmask)[0]
    nanmask = ~warmmask

    return idx_warm, newday, nanmask
