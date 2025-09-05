#!/usr/bin/env python

# @package ion_functions.data.adcp_functions
# @file ion_functions/data/adcp_functions.py
# @author Christopher Wingard, Russell Desiderio, Craig Risien
# @brief Module containing ADCP related data-calculations.

"""
Module containing ADCP related data-calculations.
This module implements the algorithms for calculating ADCP velocity profiles 
and echo intensity from ADCP beam coordinate transformed velocity profiles.
This module also implements the algorithms for calculating ADCP velocity bin 
depths for the pd0 and pd8 output formats. This module is used by the OOI 
Cyberinfrastructure to calculate the L1 and L2 data products from the L0 data 
products for the ADCP family of instruments.

## Overview of functions
### For instruments programmed in beam coordinates:
    (ADCPS-I,K;  ADCPT-B,D,E)

* adcp_beam_eastward -- calculates VELPROF-VLE_L1
* adcp_beam_northward -- calculates VELPROF-VLN_L1
* adcp_beam_vertical -- calculates VELPROF-VLU_L1
* adcp_beam_error -- calculates VELPROF-ERR_L1

### For instruments programmed in earth coordinates:
    (ADCPA;  ADCPS-J,L,N; ADCPT-C,F,G,M)

* adcp_earth_eastward -- calculates VELPROF-VLE_L1
* adcp_earth_northward -- calculates VELPROF-VLN_L1
* adcp_earth_vertical -- calculates VELPROF-VLU_L1
* adcp_earth_error -- calculates VELPROF-ERR_L1

### For the VADCP programmed in beam coordinates:
* vadcp_beam_eastward -- calculates VELTURB-VLE_L1
* vadcp_beam_northward -- calculates VELTURB-VLN_L1
* vadcp_beam_vertical_true -- calculates VELTURB-VLU-5BM_L1
* vadcp_beam_vertical_est -- calculates VELTURB-VLU-4BM_L1
* vadcp_beam_error -- calculates VELTURB-ERR_L1

### For all tRDI ADCP instruments:
* adcp_backscatter -- calculates ECHOINT-B1_L1, ECHOINT-B2_L1, ECHOINT-B3_L1, ECHOINT-B4_L1.

### Base functions used by above functions
* adcp_beam2ins -- applies the beam to instrument transform using either a 4
    or 3 beam solution for instruments programmed in beam coordinates
* adcp_ins2earth -- applies the instrument to Earth transform for all
    instruments originally programmed in beam coordinates.
* magnetic_correction -- corrects horizontal velocities for the magnetic
        variation (declination) at the measurement location.

### Supplementary functions to calculate velocity bin depths:
* adcp_bin_depths -- calculates bin depths for the pd0 output format
                        (virtually all tRDI ADCPs deployed by OOI); uses
                        TEOS-10 functions p_from_z and enthalpy_SSO_0_p.
* adcp_bin_depths_pd8 -- calculates bin depths for the pd8 output format,
                            assuming that (1) the ADCP operator recorded the
                            necessary input variables and (2) these are somehow
                            entered into the CI system.
"""

import numpy as np

from ion_functions import deprecated
from ion_functions.data.generic_functions import (
    magnetic_declination,
    replace_fill_with_nan,
)

# instrument fill value unprocessed by CI
# (bad beam velocity sentinel output by tRDI ADCP instruments)
ADCP_FILLVALUE = -32768

def adcp_beam_velocity(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt):
    """
    Compute Earth referenced velocity data from beam coordinate transformed velocity profiles.

    This function returns all velocity components in Earth coordinates, corrected for magnetic declination.
    It can be used for all VADCP processing, with the vertical velocity component representing an estimate
    of the vertical velocity.

    Parameters
    ----------
    b1 : array_like
        "beam 1" velocity profiles in beam coordinates [mm s-1].
    b2 : array_like
        "beam 2" velocity profiles in beam coordinates [mm s-1].
    b3 : array_like
        "beam 3" velocity profiles in beam coordinates [mm s-1].
    b4 : array_like
        "beam 4" velocity profiles in beam coordinates [mm s-1].
    pg1 : float or array_like
        Percent good estimate for beam 1 [%].
    pg2 : float or array_like
        Percent good estimate for beam 2 [%].
    pg3 : float or array_like
        Percent good estimate for beam 3 [%].
    pg4 : float or array_like
        Percent good estimate for beam 4 [%].
    h : float or array_like
        Instrument's uncorrected magnetic heading [cdegrees].
    p : float or array_like
        Instrument pitch [cdegrees].
    r : float or array_like
        Instrument roll [cdegrees].
    vf : int or array_like
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking).
    lat : float or array_like
        Instrument's deployment latitude [decimal degrees].
    lon : float or array_like
        Instrument's deployment longitude [decimal degrees].
    dt : float or array_like
        Sample date and time value [seconds since 1970-01-01] (Unix Time Format).

    Returns
    -------
    u : ndarray
        East velocity profiles in Earth coordinates, corrected for magnetic declination [m s-1].
    v : ndarray
        North velocity profiles in Earth coordinates, corrected for magnetic declination [m s-1].
    w : ndarray
        Vertical velocity profiles in Earth coordinates [m s-1].
    e : ndarray
        Error velocity profiles in Earth coordinates [m s-1].

    Examples
    --------
    >>> u, v, w, e = adcp_beam_velocity(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt)

    Notes
    -----
    Relies on the base functions [adcp_beam2ins][ion_functions.data.adcp_functions.adcp_beam2ins], [adcp_ins2earth][ion_functions.data.adcp_functions.adcp_ins2earth], and [magnetic_correction][ion_functions.data.adcp_functions.magnetic_correction].

    References
    ----------
    .. [1] OOI (2012). Data Product Specification for Velocity Profile and Echo Intensity. Document Control Number 1341-00750.
    .. [2] OOI (2013). Data Product Specification for Turbulent Velocity Profile and Echo Intensity. Document Control Number 1341-00760.
    """
    # force shapes of inputs to arrays of the correct dimensions
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    x, y, z, e = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # compute the instrument to earth beam transform
    u, v, w = adcp_ins2earth(x, y, z, h, p, r, vf)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it.
    u_cor, v_cor = magnetic_correction(theta, u, v)

    # scale velocities to m/s from mm/s
    u_cor = u_cor / 1000.
    v_cor = v_cor / 1000.
    w = w / 1000.
    e = e / 1000.

    # return the velocity profiles
    return u_cor, v_cor, w, e


# Wrapper functions to create the VELPROF L1 data products for instruments
# programmed in beam coordinates by RSN (ADCPS-I,K and ADCPT-B,D,E)
def adcp_beam_eastward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt):
    """
    Wrapper function to compute the Eastward Velocity Profile (VELPROF-VLE) 
    from beam coordinate transformed velocity profiles.

    Parameters
    ----------
    b1 : array_like
        "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1].
    b2 : array_like
        "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1].
    b3 : array_like
        "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1].
    b4 : array_like
        "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1].
    pg1 : float or array_like
        Percent good estimate for beam 1 [%].
    pg2 : float or array_like
        Percent good estimate for beam 2 [%].
    pg3 : float or array_like
        Percent good estimate for beam 3 [%].
    pg4 : float or array_like
        Percent good estimate for beam 4 [%].
    h : float or array_like
        Instrument's uncorrected magnetic heading [cdegrees]
    p : float or array_like
        Instrument pitch [cdegrees]
    r : float or array_like
        Instrument roll [cdegrees]
    vf : int or array_like
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking)
    lat : float or array_like
        Instrument's deployment latitude [decimal degrees]
    lon : float or array_like
        Instrument's deployment longitude [decimal degrees]
    dt : float or array_like
        Sample date and time value [seconds since 1900-01-01] (NTP Time Format)

    Notes
    -----
    Uses the [adcp_beam_velocity][ion_functions.data.adcp_functions.adcp_beam_velocity] 
    function to compute the eastward velocity component.
    
    Returns
    -------
    u_cor : array_like
        Eastward velocity profiles in Earth coordinates corrected for the 
        magnetic declination (VELPROF-VLE_L1) [m s-1]
    """

    # Convert the given ntp epoch timestamp in unix epoch timestamp.
    dt_unix = dt - 2208988800

    # call central adcp_beam_velocity function to get specific eastward velocity profile
    u_cor, _, _, _ = adcp_beam_velocity(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt_unix)

    # return the eastward velocity profile
    return u_cor


def adcp_beam_northward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt):
    """
    Wrapper function to compute the Northward Velocity Profile (VELPROF-VLN) 
    from beam coordinate transformed velocity profiles.

    Parameters
    ----------
    b1 : array_like
        "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1].
    b2 : array_like
        "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1].
    b3 : array_like
        "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1].
    b4 : array_like
        "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1].
    pg1 : float or array_like
        Percent good estimate for beam 1 [%].
    pg2 : float or array_like
        Percent good estimate for beam 2 [%].
    pg3 : float or array_like
        Percent good estimate for beam 3 [%].
    pg4 : float or array_like
        Percent good estimate for beam 4 [%].
    h : float or array_like
        Instrument's uncorrected magnetic heading [cdegrees]
    p : float or array_like
        Instrument pitch [cdegrees]
    r : float or array_like
        Instrument roll [cdegrees]
    vf : int or array_like
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking)
    lat : float or array_like
        Instrument's deployment latitude [decimal degrees]
    lon : float or array_like
        Instrument's deployment longitude [decimal degrees]
    dt : float or array_like
        Sample date and time value [seconds since 1900-01-01] (NTP Time Format)

    Notes
    -----
    Uses the [adcp_beam_velocity][ion_functions.data.adcp_functions.adcp_beam_velocity] 
    function to compute the northward velocity component.

    Returns
    -------
    v_cor : array_like
        Northward velocity profiles in Earth coordinates corrected for the 
        magnetic declination (VELPROF-VLN_L1) [m s-1]
    """

    # Convert the given ntp epoch timestamp in unix epoch timestamp.
    dt_unix = dt - 2208988800

    # call central adcp_beam_velocity function to get specific eastward velocity profile
    _, v_cor, _, _ = adcp_beam_velocity(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt_unix)

    # return the northward velocity profile
    return v_cor


def adcp_beam_vertical(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf):
    """
    Wrapper function to compute the Upward Velocity Profile (VELPROF-VLU) 
    from beam coordinate transformed velocity profiles.

    Parameters
    ----------
    b1 : array_like
        "beam 1" velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1].
    b2 : array_like
        "beam 2" velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1].
    b3 : array_like
        "beam 3" velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1].
    b4 : array_like
        "beam 4" velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1].
    pg1 : float or array_like
        Percent good estimate for beam 1 [%].
    pg2 : float or array_like
        Percent good estimate for beam 2 [%].
    pg3 : float or array_like
        Percent good estimate for beam 3 [%].
    pg4 : float or array_like
        Percent good estimate for beam 4 [%].
    h : float or array_like
        Instrument's uncorrected magnetic heading [cdegrees]
    p : float or array_like
        Instrument pitch [cdegrees]
    r : float or array_like
        Instrument roll [cdegrees]
    vf : int or array_like
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking)

    Returns
    -------
    w : array_like
        Vertical velocity profiles (VELPROF-VLU_L1) [m s-1]
    """
    # compute the beam to instrument transform
    x, y, z, _ = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # compute the instrument to earth beam transform
    _, _, w = adcp_ins2earth(x, y, z, h, p, r, vf)

    # scale the vertical velocity to m/s
    w = w / 1000.  # mm/s -> m/s

    # return the vertical velocity profile
    return w


def adcp_beam_error(b1, b2, b3, b4, pg1, pg2, pg3, pg4):
    """
    Wrapper function to compute the Error Velocity Profile (VELPROF-ERR) 
    from beam coordinate transformed velocity profiles.

    Parameters
    ----------
    b1 : array_like
        beam 1 velocity profiles in beam coordinates (VELPROF-B1_L0) [mm s-1].
    b2 : array_like
        beam 2 velocity profiles in beam coordinates (VELPROF-B2_L0) [mm s-1].
    b3 : array_like
        beam 3 velocity profiles in beam coordinates (VELPROF-B3_L0) [mm s-1].
    b4 : array_like
        beam 4 velocity profiles in beam coordinates (VELPROF-B4_L0) [mm s-1].
    pg1 : float or array_like
        Percent good estimate for beam 1 [%].
    pg2 : float or array_like
        Percent good estimate for beam 2 [%].
    pg3 : float or array_like
        Percent good estimate for beam 3 [%].
    pg4 : float or array_like
        Percent good estimate for beam 4 [%].

    Returns
    -------
    e : array_like
        Error velocity profiles (VELPROF-ERR_L1) [m s-1]
    """
    # compute the beam to instrument transform
    _, _, _, e = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # scale error velocity to m/s
    e = e / 1000.   # mm/s

    # return the Error Velocity Profile
    return e


# Wrapper functions to create the VELPROF L1 data products for instruments
# programmed in Earth coordinates by CGSN (Pioneer and Endurance) (ADCPA,
# ADCPS-J,L,N and ADCPT-C,F,G,M)
def adcp_earth_eastward(u, v, z, lat, lon, dt):
    """
    Wrapper function to compute the Eastward Velocity Profile (VELPROF-VLE)
    from Earth coordinate transformed velocity profiles.

    Parameters
    ----------
    u : array_like
        Eastward velocity profiles (VELPROF-VLE_L0) [mm s-1]
    v : array_like
        Northward velocity profiles (VELPROF-VLN_L0) [mm s-1]
    z : array_like
        Instrument's pressure sensor reading (depth) [daPa]
    lat : float or array_like
        Instrument's deployment latitude [decimal degrees]
    lon : float or array_like
        Instrument's deployment longitude [decimal degrees]
    dt : float or array_like
        Sample date and time value [seconds since 1900-01-01]

    Returns
    -------
    uu_cor : array_like
        Eastward velocity profiles in Earth coordinates corrected for the 
        magnetic declination (VELPROF-VLE_L1) [m s-1]
    """
    # force shapes of inputs to arrays
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)

    # on input, the elements of u and v are of type int.
    u, v = replace_fill_with_nan(ADCP_FILLVALUE, u, v)

    #z = np.atleast_1d(z) / 1000.  # scale daPa depth input to dbar
    #z = z * 1.019716  # use a simple approximation to calculate depth in m
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    uu_cor, _ = magnetic_correction(theta, u, v)

    # scale velocity to m/s
    uu_cor = uu_cor / 1000.  # mm/s -> m/s

    # return the Eastward Velocity Profile
    return uu_cor


def adcp_earth_northward(u, v, z, lat, lon, dt):
    """
    Wrapper function to compute the Northward Velocity Profile (VELPROF-VLN) 
    from Earth coordinate transformed velocity profiles.

    Parameters
    ----------
    u : array_like
        Eastward velocity profiles (VELPROF-VLE_L0) [mm s-1]
    v : array_like
        Northward velocity profiles (VELPROF-VLN_L0) [mm s-1]
    z : array_like
        Instrument's pressure sensor reading (depth) [daPa]
    lat : float or array_like
        Instrument's deployment latitude [decimal degrees]
    lon : float or array_like
        Instrument's deployment longitude [decimal degrees]
    dt : float or array_like
        Sample date and time value [seconds since 1900-01-01]

    Returns
    -------
    vv_cor : array_like
        Northward velocity profiles in Earth coordinates corrected for the 
        magnetic declination (VELPROF-VLN_L1) [m s-1]
    """
    # force shapes of inputs to arrays
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)

    # on input, the elements of u and v are of type int.
    u, v = replace_fill_with_nan(ADCP_FILLVALUE, u, v)

    #z = np.atleast_1d(z) / 1000.  # scale daPa depth input to dbar
    #z = z * 1.019716  # use a simple approximation to calculate depth in m
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    _, vv_cor = magnetic_correction(theta, u, v)

    # scale velocity to m/s
    vv_cor = vv_cor / 1000.  # mm/s -> m/s

    # return the Northward Velocity Profile
    return vv_cor


def adcp_earth_vertical(w):
    """
    Wrapper function to compute the Upward Velocity Profile (VELPROF-VLU) 
    from Earth coordinate transformed velocity profiles.

    Parameters
    ----------
    w : array_like
        Upward velocity profiles (VELPROF-VLU_L0) [mm s-1]

    Returns
    -------
    w_scl : array_like
        Scaled upward velocity profiles in Earth coordinates (VELPROF-VLN_L1) [m s-1]
    """
    w = replace_fill_with_nan(ADCP_FILLVALUE, w)

    # scale velocity to m/s
    w_scl = w / 1000.  # mm/s -> m/s

    # return the Upward Velocity Profile
    return w_scl


def adcp_earth_error(e):
    """
    Wrapper function to compute the Error Velocity Profile (VELPROF-ERR) 
    from Earth coordinate transformed velocity profiles.

    Parameters
    ----------
    e : array_like
        Error velocity profiles (VELPROF-ERR_L0) [mm s-1]

    Returns
    -------
    e_scl : array_like
        Scaled error velocity profiles in Earth coordinates (VELPROF-ERR_L1) [m s-1]
    """
    e = replace_fill_with_nan(ADCP_FILLVALUE, e)

    # scale velocity to m/s
    e_scl = e / 1000.  # mm/s -> m/s

    # return the scaled Error Velocity Profile
    return e_scl


# Compute the VELTURB_L1 data products for the VADCP instrument deployed by RSN.
@deprecated
def vadcp_beam_eastward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt):
    """
    Wrapper function to compute the Eastward Velocity Profile (VELTURB-VLE)
    from beam coordinate transformed velocity profiles.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument's uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like or int
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking).
    lat : array_like
        Instrument's deployment latitude [decimal degrees].
    lon : array_like
        Instrument's deployment longitude [decimal degrees].
    dt : array_like
        Sample date and time value [seconds since 1900-01-01].

    Returns
    -------
    u_cor : ndarray
        Eastward velocity profiles in Earth coordinates corrected for the 
        magnetic declination (VELTURB-VLE_L1) [m s-1].

    Notes
    -----
    - Input velocities are expected in mm/s and output is in m/s.

    References
    ----------
    - Data Product Specification for Turbulent Velocity Profile and Echo Intensity - DCN 1341-00760.

    Examples
    --------
    >>> u_cor = vadcp_beam_eastward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt)
    """
    # force shapes of some inputs to arrays of the correct dimensions
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    x, y, z, _ = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # compute the instrument to earth beam transform
    u, v, _ = adcp_ins2earth(x, y, z, h, p, r, vf)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    u_cor, _ = magnetic_correction(theta, u, v)

    # scale velocity to m/s
    u_cor = u_cor / 1000.  # mm/s -> m/s

    # return the eastward velocity profile
    return u_cor


@deprecated
def vadcp_beam_northward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt):
    """
    Wrapper function to compute the Northward Velocity Profile (VELTURB-VLN)
    from beam coordinate transformed velocity profiles.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument's uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like or int
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking).
    lat : array_like
        Instrument's deployment latitude [decimal degrees].
    lon : array_like
        Instrument's deployment longitude [decimal degrees].
    dt : array_like
        Sample date and time value [seconds since 1900-01-01].

    Returns
    -------
    v_cor : ndarray
        Northward velocity profiles in Earth coordinates corrected for the magnetic declination (VELTURB-VLN_L1) [m s-1].

    Notes
    -----
    - Input velocities are expected in mm/s and output is in m/s.

    References
    ----------
    - Data Product Specification for Turbulent Velocity Profile and Echo Intensity - DCN 1341-00760.

    Examples
    --------
    >>> v_cor = vadcp_beam_northward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt)
    """
    # force shapes of some inputs to arrays of the correct dimensions
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # compute the beam to instrument transform
    x, y, z, _ = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # compute the instrument to earth beam transform
    u, v, _ = adcp_ins2earth(x, y, z, h, p, r, vf)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    _, v_cor = magnetic_correction(theta, u, v)

    # scale velocity to m/s
    v_cor = v_cor / 1000.  # mm/s -> m/s

    # return the northward velocity profile
    return v_cor


def vadcp_beam_vertical_est(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf):
    """
     Wrapper function to compute the estimated upward velocity profile 
     (VELTURB-VLU-4BM) from beam coordinate transformed velocity profiles
    using a 4- or 3-beam solution, where each beam is oriented facing outward 
    at 20 degrees relative to vertical.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument's uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like or int
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking).

    Returns
    -------
    w : ndarray
        Estimated vertical velocity profiles in Earth coordinates (VELTURB-VLU-4BM_L1) [m s-1].

    Notes
    -----
    - Input velocities are expected in mm/s and output is in m/s.

    References
    ----------
    Data Product Specification for Turbulent Velocity Profile and Echo Intensity - DCN 1341-00760.
    """
    # compute the beam to instrument transform
    x, y, z, _ = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # compute the instrument to earth beam transform
    _, _, w = adcp_ins2earth(x, y, z, h, p, r, vf)

    # scale upward velocity to m/s
    w = w / 1000.  # mm/s -> m/s

    # return the estimated Upward Velocity Profile
    return w


def vadcp_beam_vertical_true(b1, b2, b3, b4, b5, pg1, pg2, pg3, pg4, pg5, h, p, r, vf):
    """Computes the "true" Upward Velocity Profile (VELTURB-VLU-5BM) from the 
    beam coordinate transformed velocity profiles.
    This provides a better vertical velocity estimate since beam 5 is oriented vertically.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1].
    b5 : array_like
        Beam 5 velocity profiles in beam coordinates (VELTURB-B5_L0) [mm s-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    pg5 : array_like
        Percent good estimate for beam 5 [percent].
    h : array_like
        Instrument's uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like or int
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking).
    Returns
    -------
    w : ndarray
        True vertical velocity profiles in Earth coordinates (VELTURB-VLU-5BM_L1) [m s-1].

    Notes
    -----
    - Input velocities are expected in mm/s and output is in m/s.

    References
    ----------
    - Data Product Specification for Turbulent Velocity Profile and Echo Intensity - DCN 1341-00760.
    """
    # compute the beam to instrument transform
    x, y, _, _ = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # check percent good data for beam 5, reset to fill if less than 25%
    b5 = np.ma.filled(np.ma.masked_where(pg5 < 25, b5), ADCP_FILLVALUE)

    # check b5 for the presence of fill values and replace with NaN
    b5 = replace_fill_with_nan(ADCP_FILLVALUE, b5)

    # compute the instrument to earth beam transform
    _, _, w = adcp_ins2earth(x, y, b5, h, p, r, vf)

    # scale upward velocity to m/s
    w = w / 1000.  # mm/s -> m/s

    # return the true Upward Velocity Profile
    return w


@deprecated
def vadcp_beam_error(b1, b2, b3, b4, pg1, pg2, pg3, pg4):
    """
    Wrapper function to compute the Error Velocity Profile (VELTURB-ERR)
    from the beam coordinate transformed velocity profiles.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1].
    pg1 : array_like
        Percent good estimate for beam 1 [%].
    pg2 : array_like
        Percent good estimate for beam 2 [%].
    pg3 : array_like
        Percent good estimate for beam 3 [%].
    pg4 : array_like
        Percent good estimate for beam 4 [%].

    Returns
    -------
    e : ndarray
        Error velocity profiles (VELTURB-ERR_L1) [m s-1].

    Examples
    --------
    >>> e = vadcp_beam_error(b1, b2, b3, b4, pg1, pg2, pg3, pg4)
    """
    # compute the beam to instrument transform
    _, _, _, e = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # scale error velocity to m/s
    e = e / 1000.   # mm/s

    # return the Error Velocity Profile
    return e


# Calculates ECHOINT_L1 for all tRDI ADCPs
def adcp_backscatter(raw, sfactor=0.45):
    """
    Converts the echo intensity data from counts to dB using a factory
    specified scale factor.

    Parameters
    ----------
    raw : array_like
        Raw echo intensity (ECHOINT_L0) [count].
    sfactor : float or array_like, optional
        Factory supplied scale factor, instrument and beam specific [dB/count].
        Default is 0.45.

    Returns
    -------
    dB : array_like
        Relative Echo Intensity (ECHOINT_L1) [dB].

    Notes
    -----
    * The ADCP outputs the raw echo intensity as a 1-byte integer, so the ADCP_FILLVALUE
    cannot apply (requires 2 bytes).
    * The default scale factor is nominally 0.45 dB/count for the Workhorse
        family of ADCPs and 0.61 dB/count for the ExplorerDVL family.

    References
    ----------
    OOI (2012). Data Product Specification for Velocity Profile and Echo Intensity.
        Document Control Number 1341-00750.
        https://alfresco.oceanobservatories.org/
    """
    if np.isscalar(sfactor) is False:
        sfactor = sfactor.reshape(sfactor.shape[0], 1)

    # check raw for the presence of system fill values
    raw = replace_fill_with_nan(None, raw)

    dB = raw * sfactor
    return dB


##### ADCP Beam to Earth Transforms and Magnetic Variation Corrections
def adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4):
    """
    Converts the Beam Coordinate transformed velocity profiles to the instrument coordinate system.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0) [mm s-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0) [mm s-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0) [mm s-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0) [mm s-1].
    pg1 : array_like
        Percent good estimate for beam 1 [%].
    pg2 : array_like
        Percent good estimate for beam 2 [%].
    pg3 : array_like
        Percent good estimate for beam 3 [%].
    pg4 : array_like
        Percent good estimate for beam 4 [%].

    Returns
    -------
    x : array_like
        x axis velocity profiles in instrument coordinates [mm s-1].
    y : array_like
        y axis velocity profiles in instrument coordinates [mm s-1].
    z : array_like
        z axis velocity profiles in instrument coordinates [mm s-1].
    e : array_like
        Error velocity profiles [mm s-1].

    References
    ----------
    OOI (2012). Data Product Specification for Velocity Profile and Echo Intensity. Document Control Number
        1341-00750. https://alfresco.oceanobservatories.org/
    OOI (2013). Data Product Specification for Turbulent Velocity Profile and Echo Intensity. Document Control
        Number 1341-00760. https://alfresco.oceanobservatories.org/
    Teledyne RD Instruments (2008). ADCP Coordinate Transformation, Formulas and Calculations.
    """
    # raw beam velocities, set to correct shape
    b1 = np.atleast_2d(b1)
    b2 = np.atleast_2d(b2)
    b3 = np.atleast_2d(b3)
    b4 = np.atleast_2d(b4)

    # percentage of good pings for each beam per depth cell, set to correct shape
    pg1 = np.atleast_2d(pg1)
    pg2 = np.atleast_2d(pg2)
    pg3 = np.atleast_2d(pg3)
    pg4 = np.atleast_2d(pg4)

    # using the vendor specified percent good floor of 25%, create masked arrays with fill values set to compute
    # a 3-beam solution, if applicable.
    ma1 = np.ma.masked_where(pg1 < 25, b1)
    bm1 = ma1.filled((b2 - b3 - b4) * -1)
    ma2 = np.ma.masked_where(pg2 < 25, b2)
    bm2 = ma2.filled((b1 - b3 - b4) * -1)
    ma3 = np.ma.masked_where(pg3 < 25, b3)
    bm3 = ma3.filled(b1 + b2 - b4)
    ma4 = np.ma.masked_where(pg4 < 25, b4)
    bm4 = ma4.filled(b1 + b2 - b3)

    # sum across the masked arrays to determine if more than 1 beam is bad per depth cell, if so we cannot compute a
    # 3-beam solution and need to set the fill value to a NaN.
    mad = np.ma.dstack((ma1, ma2, ma3, ma4))    # stack the masked arrays in depth
    mas = np.ma.count_masked(mad, axis=2)       # count the number of depth cells masked in the depth stacked array

    # using the above, reset the raw beams. fill with 3-beam if applicable, otherwise use a NaN
    bm1 = np.ma.filled(np.ma.masked_where(mas > 1, bm1), ADCP_FILLVALUE)
    bm2 = np.ma.filled(np.ma.masked_where(mas > 1, bm2), ADCP_FILLVALUE)
    bm3 = np.ma.filled(np.ma.masked_where(mas > 1, bm3), ADCP_FILLVALUE)
    bm4 = np.ma.filled(np.ma.masked_where(mas > 1, bm4), ADCP_FILLVALUE)

    bm1, bm2, bm3, bm4 = replace_fill_with_nan(ADCP_FILLVALUE, bm1, bm2, bm3, bm4)

    theta = 20.0 / 180.0 * np.pi
    a = 1.0 / (2.0 * np.sin(theta))
    b = 1.0 / (4.0 * np.cos(theta))
    c = 1.0   # +1.0 for convex transducer head, -1 for concave
    d = a / np.sqrt(2.0)

    x = c * a * (bm1 - bm2)
    y = c * a * (bm4 - bm3)
    z = b * (bm1 + bm2 + bm3 + bm4)
    e = d * (bm1 + bm2 - bm3 - bm4)

    return x, y, z, e


def adcp_ins2earth(u, v, w, heading, pitch, roll, vertical):
    """
    Converts the Instrument Coordinate transformed velocity profiles to the Earth coordinate system.

    Parameters
    ----------
    u : array_like
        East velocity profiles in instrument coordinates [mm s-1]
    v : array_like
        North velocity profiles in instrument coordinates [mm s-1]
    w : array_like
        Vertical velocity profiles in instrument coordinates [mm s-1]
    heading : float or array_like
        Instrument's uncorrected magnetic heading [centidegrees]
    pitch : float or array_like
        Instrument pitch [centidegrees]
    roll : float or array_like
        Instrument roll [centidegrees]
    vertical : int or array_like
        Instrument's vertical orientation (0 = downward looking, 1 = upward looking)

    Returns
    -------
    uu : array_like
        "East" velocity profiles in earth coordinates [mm s-1]
    vv : array_like
        "North" velocity profiles in earth coordinates [mm s-1]
    ww : array_like
        "Vertical" velocity profiles in earth coordinates [mm s-1]

    References
    ----------
    OOI (2012). Data Product Specification for Velocity Profile and Echo Intensity. Document Control Number
        1341-00750. https://alfresco.oceanobservatories.org/
    """
    ### the input beam data for adcp_ins2earth are always called using the output
    ### of adcp_beam2ins, so the following lines are not needed.
    # insure we are dealing with array inputs
    #u = np.atleast_2d(u)
    #v = np.atleast_2d(v)
    #w = np.atleast_2d(w)

    # check for CI fill values before changing units.
    # this function 'conditions' (np.atleast_1d) its inputs.
    # TRDI does not apply its ADCP fill/bad value sentinels to compass data.
    heading, pitch, roll, vertical = replace_fill_with_nan(None, heading, pitch, roll, vertical)

    # change units from centidegrees to degrees
    heading = heading / 100.0
    pitch = pitch / 100.0
    roll = roll / 100.0

    # better way to calculate roll from the vertical orientation toggle;
    # this will propagate R as nans if the vertical variable is missing from the data.
    r = roll + vertical * 180.0

    # roll
    r_rad = np.radians(r)
    cos_r = np.cos(r_rad)
    sin_r = np.sin(r_rad)
    # heading
    h_rad = np.radians(heading)
    cos_h = np.cos(h_rad)
    sin_h = np.sin(h_rad)
    # pitch
    t1rad = np.radians(pitch)
    t2rad = np.radians(roll)
    p_rad = np.arctan(np.tan(t1rad) * np.cos(t2rad))
    cos_p = np.cos(p_rad)
    sin_p = np.sin(p_rad)

    # determine array size
    n_packets = u.shape[0]
    n_uvw = u.shape[1]

    # initialize vectors to be used as matrix elements
    ones = np.ones(n_packets)
    zeros = ones * 0.0

    # the rollaxis calls reorient the matrices so that their lead index is
    # the data packet index
    m1 = np.array([[cos_h, sin_h, zeros],
                   [-sin_h, cos_h, zeros],
                   [zeros, zeros, ones]])
    m1 = np.rollaxis(m1, 2)
    m2 = np.array([[ones, zeros, zeros],
                   [zeros, cos_p, -sin_p],
                   [zeros, sin_p, cos_p]])
    m2 = np.rollaxis(m2, 2)
    m3 = np.array([[cos_r, zeros, sin_r],
                   [zeros, ones, zeros],
                   [-sin_r, zeros, cos_r]])
    m3 = np.rollaxis(m3, 2)

    # construct input array of coordinates (velocities) to be transformed.
    # the basis set is 3D (E,N,U) so that the middle dimension is sized at 3.
    uvw = np.zeros((n_packets, 3, n_uvw))

    # pack the coordinates (velocities) to be transformed into the appropriate
    # slices.
    uvw[:, 0, :] = u
    uvw[:, 1, :] = v
    uvw[:, 2, :] = w

    # the Einstein summation is here configured to do the matrix
    # multiplication MM(i,l) = M1(i,j) * M2(j,k) * M3(k,l) on each slice h.
    mm = np.einsum('hij,hjk,hkl->hil', m1, m2, m3)

    # the Einstein summation is here configured to do the matrix
    # multiplication uvw_earth(i,m) = MM(i,l) * uvw(l,m) on each slice h.
    uvw_earth = np.einsum('hil,hlm->him', mm, uvw)

    # NOTE:
    # these last two executable statements run about a factor of 2
    # faster in the 10000 data packet performance tests versus combining
    # these operations into the one statement:
    #     uvw_earth = np.einsum('hij,hjk,hkl,hlm->him', M1, M2, M3, uvw)

    # break out the coordinate slices and return them
    uu = uvw_earth[:, 0, :]
    vv = uvw_earth[:, 1, :]
    ww = uvw_earth[:, 2, :]

    return uu, vv, ww


def magnetic_correction(theta, u, v):
    
    """
    Corrects velocity profiles for the magnetic variation (declination) at the 
    measurement location.
    The magnetic declination is obtained from the 2010 World Magnetic Model 
    (WMM2010) provided by NOAA (see wmm_declination).

    Parameters
    ----------
    theta : float or array_like
        Magnetic variation based on location (latitude, longitude, altitude) 
        and date [degrees]
    u : array_like
        Uncorrected eastward velocity profiles in earth coordinates
    v : array_like
        Uncorrected northward velocity profiles in earth coordinates

    Returns
    -------
    u_cor : array_like
        Eastward velocity profiles, in earth coordinates, with the correction 
        for magnetic variation applied.
    v_cor : array_like
        Northward velocity profiles, in earth coordinates, with the correction 
        for magnetic variation applied.

    Notes
    -----
    This version handles 'vectorized' input variables without using for
    loops. It was specifically written to handle the case of a 1D array of
    theta values, theta=f(i), with corresponding sets of 'u' and 'v' values
    such that u=f(i,j) and v=f(i,j), where there are j 'u' and 'v' values
    for each theta(i).
    
    References
    ----------
    OOI (2012). Data Product Specification for Velocity Profile and Echo Intensity. Document Control Number
        1341-00750. https://alfresco.oceanobservatories.org/
    OOI (2013). Data Product Specification for Turbulent Velocity Profile and Echo Intensity. Document Control
        Number 1341-00760. https://alfresco.oceanobservatories.org/
    """
    # force shapes of inputs to arrays
    theta = np.atleast_1d(theta)
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)

    r_theta = np.radians(theta)
    cos_t = np.cos(r_theta)
    sin_t = np.sin(r_theta)

    m = np.array([[cos_t, sin_t],
                  [-sin_t, cos_t]])

    # roll axes so that the lead index represents data packet #.
    m = np.rollaxis(m, 2)

    # the coordinate system is 2D, so the middle dimension is sized at 2.
    uv = np.zeros((u.shape[0], 2, u.shape[1]))

    # pack the coordinates to be rotated into the appropriate slices
    uv[:, 0, :] = u
    uv[:, 1, :] = v

    # the Einstein summation is here configured to do the matrix
    # multiplication uv_cor(i,k) = M(i,j) * uv(j,k) on each slice h.
    uv_cor = np.einsum('hij,hjk->hik', m, uv)

    # the magnetically corrected u values are:
    u_cor = uv_cor[:, 0, :]

    # the magnetically corrected v values are:
    v_cor = uv_cor[:, 1, :]

    # return corrected u and v values
    return (u_cor, v_cor)


def adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude):
    """
    Calculates the center bin depths for PD0 and PD12 ADCP data.

    Parameters
    ----------
    dist_first_bin : float or array_like
        Distance to the first ADCP bin [centimeters]
    bin_size : float or array_like
        Depth of each ADCP bin [centimeters]
    num_bins : int or array_like
        Number of ADCP bins [unitless]
    pressure : float or array_like
        Pressure at the sensor head [bar]
    adcp_orientation : int or array_like
        1=upward looking or 0=downward looking [unitless]
    latitude : float or array_like
        Latitude of the instrument [degrees]

    Returns
    -------
    bin_depths : array_like
        Bin depths [meters]

    References
    ----------
    OOI (2012). Data Product Specification for Velocity Profile and Echo Intensity. Document Control Number
        1341-00750. https://alfresco.oceanobservatories.org/
    """
    # check for CI fill values.
    pressure = replace_fill_with_nan(None, pressure)

    # Convert pressure from bar to decibar
    pressure_dbar = pressure * 10.0

    # Calculate sensor depth using TEOS-10 toolbox z_from_p function
    # note change of sign to make the sensor_depth variable positive
    sensor_depth = -z_from_p(pressure_dbar, latitude)

    return adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)


def adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude):
    """
    Calculates the center bin depths for PD0 and PD12 ADCP data.

    Parameters
    ----------
    dist_first_bin : float or array_like
        Distance to the first ADCP bin [centimeters]
    bin_size : float or array_like
        Depth of each ADCP bin [centimeters]
    num_bins : int or array_like
        Number of ADCP bins [unitless]
    pressure : float or array_like
        Pressure at the sensor head [daPa]
    adcp_orientation : int or array_like
        1=upward looking or 0=downward looking [unitless]
    latitude : float or array_like
        Latitude of the instrument [degrees]

    Returns
    -------
    bin_depths : array_like
        Bin depths [meters]

    References
    ----------
    OOI (2012). Data Product Specification for Velocity Profile and Echo Intensity. Document Control Number
        1341-00750. https://alfresco.oceanobservatories.org/
    """
    # check for CI fill values.
    pressure = replace_fill_with_nan(None, pressure)

    # Convert pressure from decaPascal to decibar
    pressure_dbar = pressure / 1000.0

    # Calculate sensor depth using TEOS-10 toolbox z_from_p function
    # note change of sign to make the sensor_depth variable positive
    sensor_depth = -z_from_p(pressure_dbar, latitude)

    return adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)


def adcp_bin_depths(blanking_distance, bin_size, number_bins, orientation, depth):
    """
    Calculates the center bin depths for ADCP data.

    Parameters
    ----------
    blanking_distance : float or array_like
        Distance to the first ADCP bin [centimeters]
    bin_size : float or array_like
        Size, or cell length, of each ADCP bin [centimeters]
    number_bins : int or array_like
        Number of ADCP bins [unitless]
    orientation : int or array_like
        1=upward looking or 0=downward looking [unitless]
    depth : float or array_like
        Depth of the sensor head [m]

    Returns
    -------
    bin_depths : array_like
        Bin depths [meters]

    References
    ----------
    OOI (2012). Data Product Specification for Velocity Profile and Echo Intensity. Document Control Number
        1341-00750. https://alfresco.oceanobservatories.org/
    """
    # Convert from cm to meters
    blanking_distance = blanking_distance / 100.0
    bin_size = bin_size / 100.0

    # Following the PD0 convention, where:
    #     orientation = 0 is downward looking, bin depths are added to sensor depth
    #                 = 1 is upward looking, bin depths are subtracted from sensor depth
    if orientation == 0:
        z_sign = 1
    elif orientation == 1:
        z_sign = -1

    # Calculate bin depths
    depth = np.atleast_2d(depth).T
    number_bins = np.atleast_2d(number_bins)
    bin_depths = depth + z_sign * (blanking_distance + bin_size * number_bins)

    return bin_depths


def z_from_p(p, lat, geo_strf_dyn_height=0, sea_surface_geopotential=0):
    """
    Calculates height from sea pressure using the computationally-efficient 
    75-term expression for density.

    Parameters
    ----------
    p : float or array_like
        Pressure [dbar]
    lat : float or array_like
        Latitude in decimal degrees north [-90..+90]
    geo_strf_dyn_height : float, optional
        Dynamic height anomaly [m^2/s^2]
    sea_surface_geopotential : float, optional
        Geopotential at zero sea pressure [m^2/s^2]

    Returns
    -------
    z : array_like
        TEOS-10 height [m] (negative below sea surface)

    Notes
    -----
    - Dynamic height anomaly, geo_strf_dyn_height, if provided, must be
    computed with its p_ref=0 (the surface). 
    - Calls a function which calculates enthalpy assuming standard ocean salinity
    and 0 degrees celsius.

    Examples
    --------

    >>> p = [10, 50, 125, 250, 600, 1000]
    >>> lat = 4
    >>> z_from_p(p, lat) =
    [  -9.9445834469453,  -49.7180897012550, -124.2726219409978,
     -248.4700576548589, -595.8253480356214, -992.0919060719987]

    References
    ----------
    IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of seawater - 2010.
    McDougall, T.J., et al., 2003: Accurate and computationally efficient algorithms for potential temperature and density of seawater.J. Atmosph. Ocean. Tech., 20, pp. 730-741.
    Moritz, 2000: Goedetic reference system 1980. J. Geodesy, 74, 128-133.
    Roquet, F., et al., 2015: Accurate polynomial expressions for the density and specifc volume of seawater using the TEOS-10 standard. Ocean Modelling.
    Saunders, P. M., 1981: Practical conversion of pressure to depth. Journal of Physical Oceanography, 11, 573-574.
    """
    x = np.sin(np.deg2rad(lat))
    sin2 = x ** 2
    b = 9.780327 * (1.0 + (5.2792e-3 + (2.32e-5 * sin2)) * sin2)
    gamma = 2.26e-07
    a = -0.5 * gamma * b
    c = enthalpy_SSO_0_p(p) - geo_strf_dyn_height

    return -2 * c / (b + np.sqrt(b ** 2 - 4 * a * c))


def enthalpy_SSO_0_p(p):
    """
    Calculates enthalpy at the Standard Ocean Salinity, SSO, and at a Conservative Temperature of zero degrees C, as a function of pressure.

    Parameters
    ----------
    p : float or array_like
        Pressure [dbar]

    Returns
    -------
    enthalpy_SSO_0 : array_like
        Enthalpy at SSO and 0 deg C

    Notes
    -----
    - The python implementation comes directly from the matlab coding of this function:
        VERSION NUMBER: 3.05 (27th January 2015)

    References
    ----------
    Roquet, F., et al., 2015: Accurate polynomial expressions for the density and specifc volume of seawater using the TEOS-10 standard. Ocean Modelling.
    """
    z = p * 1e-4

    h006 = -2.1078768810e-9
    h007 = 2.8019291329e-10

    dynamic_enthalpy_SSO_0_p = z * (9.726613854843870e-4 + z * (-2.252956605630465e-5 + z * (
        2.376909655387404e-6 + z * (-1.664294869986011e-7 + z * (
            -5.988108894465758e-9 + z * (h006 + h007 * z))))))

    enthalpy_SSO_0 = dynamic_enthalpy_SSO_0_p * 1.e8  # Note. 1e8 = db2Pa*1e4

    return enthalpy_SSO_0


def adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation):
    """
    Calculates the center bin depths for PD0, PD8 and PD12 ADCP data.

    Parameters
    ----------
    dist_first_bin : float or array_like
        Distance to the first ADCP bin [centimeters]
    bin_size : float or array_like
        Depth of each ADCP bin [centimeters]
    num_bins : int or array_like
        Number of ADCP bins [unitless]
    sensor_depth : float or array_like
        Estimated depth at the sensor head [meters]
    adcp_orientation : int or array_like
        1=upward looking or 0=downward looking [unitless]

    Returns
    -------
    bin_depths_pd8 : array_like
        Bin depths [meters]

    Notes
    -----
    The PD8 output format is a very sparse format. Other than num_bins, it does *not* record
    any of the other input variables required by this DPA. Those must somehow be supplied "by hand".
    """
    # check for CI fill values.
    #
    # Note that these input parameters will not come from an IDD driver (except for possibly
    # (num_bins) because the PD8 output format does not output them. Therefore, I don't know
    # if they will be of type integer or not. However, ndarrays composed of float types are
    # passed through the check-code unchanged, so run the inputs through in case they are of
    # type int and in case -999999999 fill values are somehow present.
    dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation = replace_fill_with_nan(
        None, dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)

    # note, there is a CI problem not yet addressed if the time-vectorized values
    # in num_bins are not all the same!! For now, assume they are all the same:
    num_bins_constant = num_bins[0]
    # make bin_numbers a row vector
    bin_numbers = np.array([np.arange(num_bins_constant)])

    # Convert from cm to meters
    # the input variables are type integer, so divide by a real number
    # to avoid truncation errors.
    dist_first_bin = dist_first_bin / 100.0
    bin_size = bin_size / 100.0

    # make sure sensor depth is positive
    sensor_depth = np.fabs(sensor_depth)

    # Following the PD0 convention where
    #     adcp_orientation = 0 is downward looking, bindepths are added to sensor depth
    #                      = 1 is upward looking, bindepths are subtracted from sensor depth
    z_sign = 1.0 - 2.0 * adcp_orientation

    # to broadcast the vertical time dimension correctly with the horizontal bin_numbers dimension,
    # make all the 1D time arrays into column vectors to be processed with the bin_numbers row vector.
    sensor_depth = sensor_depth.reshape(-1, 1)
    z_sign = z_sign.reshape(-1, 1)
    dist_first_bin = dist_first_bin.reshape(-1, 1)
    bin_size = bin_size.reshape(-1, 1)

    # Calculate bin depths
    bin_depths_pd8 = sensor_depth + z_sign * (dist_first_bin + bin_size * bin_numbers)

    return bin_depths_pd8


def depth_from_pressure_dbar(pressure, latitude, pressure_scale_factor=None):
    """
    Calculates depth from pressure.

    Parameters
    ----------
    pressure : float or array_like
        Pressure [dbar]
    latitude : float or array_like
        Latitude of the instrument [degrees]
    pressure_scale_factor : float, optional
        Scale factor to convert pressure to dbar [unitless]

    Returns
    -------
    depth : array_like
        Depths [meters]
    """
    # check for CI fill values.
    pressure = replace_fill_with_nan(None, pressure)

    # Apply scale factor to convert pressure to decibar
    pressure_dbar = pressure * pressure_scale_factor if pressure_scale_factor else pressure;

    # Calculate sensor depth using TEOS-10 toolbox z_from_p function
    # note change of sign to make the sensor_depth variable positive
    return -z_from_p(pressure_dbar, latitude)
