#!/usr/bin/env python
"""
Functions for computing OOI ADCP family data products.

Covers VELPROF (Velocity Profile) and ECHOINT (Echo Intensity) for Teledyne
RDI Workhorse ADCP instruments (ADCPA, ADCPS, ADCPT), VELTURB (Turbulent
Velocity Profile) and ECHOINT for the original VADCP (Teledyne RDI Workhorse
Sentinel 5-beam), and VELTURB and ECHOINT for the VADCP-B (Nortek Signature
55).
"""
import numpy as np

from ion_functions import deprecated
from ion_functions.data.generic_functions import magnetic_declination
from ion_functions.data.generic_functions import replace_fill_with_nan

# instrument fill value unprocessed by CI
# (bad beam velocity sentinel output by tRDI ADCP instruments)
ADCP_FILLVALUE = -32768
VADCP_B_FILLVALUE = np.nan


def adcp_beam_velocity(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf,
                       lat, lon, dt):
    """
    Compute Earth-referenced velocity profiles from beam coordinate data.

    Applies the beam-to-instrument transform, instrument-to-Earth transform,
    and magnetic declination correction to produce VELPROF-VLE_L1,
    VELPROF-VLN_L1, VELPROF-VLU_L1, and VELPROF-ERR_L1 for instruments
    programmed in beam coordinates (ADCPS-I, ADCPS-K, ADCPT-B, ADCPT-D,
    ADCPT-E) as defined in DPS 1341-00750.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELPROF-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELPROF-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELPROF-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELPROF-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1970-01-01] (Unix time).

    Returns
    -------
    u_cor : ndarray
        Eastward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELPROF-VLE_L1) [m s^-1].
    v_cor : ndarray
        Northward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELPROF-VLN_L1) [m s^-1].
    w : ndarray
        Upward velocity profiles in Earth coordinates
        (VELPROF-VLU_L1) [m s^-1].
    e : ndarray
        Error velocity profiles in Earth coordinates
        (VELPROF-ERR_L1) [m s^-1].

    Notes
    -----
    Beam velocities are input in mm s^-1 and scaled to m s^-1 on output.
    A 3-beam solution is applied when one beam falls below 25% good; cells
    with more than one bad beam are set to NaN.
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
def adcp_beam_eastward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf,
                       lat, lon, dt):
    """
    OOI wrapper for VELPROF-VLE_L1 from beam-coordinate ADCP data.

    Returns the eastward velocity profile (VELPROF-VLE_L1) for instruments
    programmed in beam coordinates (ADCPS-I, ADCPS-K, ADCPT-B, ADCPT-D,
    ADCPT-E) as defined in DPS 1341-00750.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELPROF-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELPROF-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELPROF-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELPROF-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1900-01-01] (NTP time).

    Returns
    -------
    u_cor : ndarray
        Eastward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELPROF-VLE_L1) [m s^-1].

    See Also
    --------
    adcp_beam_velocity : Core function; use directly for multi-output access.
    """
    # Convert the given ntp epoch timestamp in unix epoch timestamp.
    dt_unix = dt - 2208988800

    # call central adcp_beam_velocity function to get specific eastward
    # velocity profile
    u_cor, _, _, _ = adcp_beam_velocity(
        b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt_unix)

    # return the eastward velocity profile
    return u_cor


def adcp_beam_northward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf,
                        lat, lon, dt):
    """
    OOI wrapper for VELPROF-VLN_L1 from beam-coordinate ADCP data.

    Returns the northward velocity profile (VELPROF-VLN_L1) for instruments
    programmed in beam coordinates (ADCPS-I, ADCPS-K, ADCPT-B, ADCPT-D,
    ADCPT-E) as defined in DPS 1341-00750.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELPROF-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELPROF-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELPROF-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELPROF-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1900-01-01] (NTP time).

    Returns
    -------
    v_cor : ndarray
        Northward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELPROF-VLN_L1) [m s^-1].

    See Also
    --------
    adcp_beam_velocity : Core function; use directly for multi-output access.
    """
    # Convert the given ntp epoch timestamp in unix epoch timestamp.
    dt_unix = dt - 2208988800

    # call central adcp_beam_velocity function to get specific northward
    # velocity profile
    _, v_cor, _, _ = adcp_beam_velocity(
        b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt_unix)

    # return the northward velocity profile
    return v_cor


def adcp_beam_vertical(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf):
    """
    OOI wrapper for VELPROF-VLU_L1 from beam-coordinate ADCP data.

    Returns the upward velocity profile (VELPROF-VLU_L1) for instruments
    programmed in beam coordinates (ADCPS-I, ADCPS-K, ADCPT-B, ADCPT-D,
    ADCPT-E) as defined in DPS 1341-00750.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELPROF-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELPROF-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELPROF-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELPROF-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).

    Returns
    -------
    w : ndarray
        Upward velocity profiles in Earth coordinates
        (VELPROF-VLU_L1) [m s^-1].

    See Also
    --------
    adcp_beam_velocity : Core function; use directly for multi-output access.
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
    OOI wrapper for VELPROF-ERR_L1 from beam-coordinate ADCP data.

    Returns the error velocity profile (VELPROF-ERR_L1) for instruments
    programmed in beam coordinates (ADCPS-I, ADCPS-K, ADCPT-B, ADCPT-D,
    ADCPT-E) as defined in DPS 1341-00750.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELPROF-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELPROF-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELPROF-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELPROF-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].

    Returns
    -------
    e : ndarray
        Error velocity profiles in Earth coordinates
        (VELPROF-ERR_L1) [m s^-1].

    See Also
    --------
    adcp_beam_velocity : Core function; use directly for multi-output access.
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
    OOI wrapper for VELPROF-VLE_L1 from Earth-coordinate ADCP data.

    Returns the eastward velocity profile (VELPROF-VLE_L1) for instruments
    programmed in Earth coordinates (ADCPA, ADCPS-J, ADCPS-L, ADCPS-N,
    ADCPT-C, ADCPT-F, ADCPT-G, ADCPT-M) as defined in DPS 1341-00750.

    Parameters
    ----------
    u : array_like
        Eastward velocity profiles (VELPROF-VLE_L0) [mm s^-1].
    v : array_like
        Northward velocity profiles (VELPROF-VLN_L0) [mm s^-1].
    z : array_like
        Instrument pressure sensor reading [daPa]. Accepted but not used
        in the current implementation; see Additional Notes on the docs
        page.
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1900-01-01] (NTP time).

    Returns
    -------
    uu_cor : ndarray
        Eastward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELPROF-VLE_L1) [m s^-1].
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
    OOI wrapper for VELPROF-VLN_L1 from Earth-coordinate ADCP data.

    Returns the northward velocity profile (VELPROF-VLN_L1) for instruments
    programmed in Earth coordinates (ADCPA, ADCPS-J, ADCPS-L, ADCPS-N,
    ADCPT-C, ADCPT-F, ADCPT-G, ADCPT-M) as defined in DPS 1341-00750.

    Parameters
    ----------
    u : array_like
        Eastward velocity profiles (VELPROF-VLE_L0) [mm s^-1].
    v : array_like
        Northward velocity profiles (VELPROF-VLN_L0) [mm s^-1].
    z : array_like
        Instrument pressure sensor reading [daPa]. Accepted but not used
        in the current implementation; see Additional Notes on the docs
        page.
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1900-01-01] (NTP time).

    Returns
    -------
    vv_cor : ndarray
        Northward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELPROF-VLN_L1) [m s^-1].
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
    OOI wrapper for VELPROF-VLU_L1 from Earth-coordinate ADCP data.

    Returns the upward velocity profile (VELPROF-VLU_L1) for instruments
    programmed in Earth coordinates as defined in DPS 1341-00750.

    Parameters
    ----------
    w : array_like
        Upward velocity profiles (VELPROF-VLU_L0) [mm s^-1].

    Returns
    -------
    w_scl : ndarray
        Upward velocity profiles in Earth coordinates
        (VELPROF-VLU_L1) [m s^-1].
    """
    w = replace_fill_with_nan(ADCP_FILLVALUE, w)

    # scale velocity to m/s
    w_scl = w / 1000.  # mm/s -> m/s

    # return the Upward Velocity Profile
    return w_scl


def adcp_earth_error(e):
    """
    OOI wrapper for VELPROF-ERR_L1 from Earth-coordinate ADCP data.

    Returns the error velocity profile (VELPROF-ERR_L1) for instruments
    programmed in Earth coordinates as defined in DPS 1341-00750.

    Parameters
    ----------
    e : array_like
        Error velocity profiles (VELPROF-ERR_L0) [mm s^-1].

    Returns
    -------
    e_scl : ndarray
        Error velocity profiles in Earth coordinates
        (VELPROF-ERR_L1) [m s^-1].
    """
    e = replace_fill_with_nan(ADCP_FILLVALUE, e)

    # scale velocity to m/s
    e_scl = e / 1000.  # mm/s -> m/s

    # return the scaled Error Velocity Profile
    return e_scl


def depth_from_dbar(pressure, ctd_pressure, deployment_depth, latitude):
    """
    Compute depth from instrument pressure with CTD and static fallbacks.

    If the onboard pressure sensor returns valid data, depth is computed
    from that pressure. If not, the nearest CTD pressure is used. If
    neither is available, the deployment depth calibration coefficient is
    returned.

    Parameters
    ----------
    pressure : array_like
        Onboard pressure [dbar].
    ctd_pressure : array_like
        Nearest CTD pressure [dbar].
    deployment_depth : float
        Deployment sheet CC_depth value [m].
    latitude : array_like
        Deployment latitude [decimal degrees].

    Returns
    -------
    depth : ndarray
        Depth [m].
    """
    if pressure.all() > 0:
        depth = -z_from_p(pressure, latitude)
    elif ctd_pressure.all() > 0:
        depth = -z_from_p(ctd_pressure, latitude)
    else:
        depth = deployment_depth

    return depth


def vadcp_b_bin_depths(depth, cell_positions, orientation):
    """
    Compute bin depths for VADCP-B from sensor depth and cell positions.

    Computes bin depths using the sensor depth, the Nortek cell_position
    parameter, and the CC_vadcpb_orientation calibration coefficient. Returns
    an empty array if orientation is not provided or is not a valid value.

    Parameters
    ----------
    depth : array_like
        Sensor depth [m].
    cell_positions : array_like
        Array of vertical distances from the instrument to each bin [m].
    orientation : array_like
        Instrument orientation (-1 = upward facing, 1 = downward facing).

    Returns
    -------
    bin_depths : ndarray
        VADCP-B bin depths [m]. Returns an empty list if orientation is
        invalid.
    """
    if not orientation.all() or orientation.all() not in [-1, 1]:
        return []

    bin_depths = []
    for i, bin in enumerate(cell_positions):
        result = depth[i] + orientation[i]*bin
        bin_depths.append(result)
    return np.array(bin_depths)


def vadcp_b_beam_eastward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r,
                          lat, lon, dt, tm):
    """
    OOI wrapper for VELPROF-VLE_L1 from VADCP-B beam-coordinate data.

    Returns the eastward velocity profile (VELPROF-VLE_L1) for the VADCP-B
    (Nortek Signature 55) as defined in DPS 1341-00750.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELPROF-B1_L0)
        [m s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELPROF-B2_L0)
        [m s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELPROF-B3_L0)
        [m s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELPROF-B4_L0)
        [m s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [degrees].
    p : array_like
        Instrument pitch [degrees].
    r : array_like
        Instrument roll [degrees].
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1900-01-01] (NTP time).
    tm : array_like
        Instrument-specific transformation matrix [4 x 4].

    Returns
    -------
    u_cor : ndarray
        Eastward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELPROF-VLE_L1) [m s^-1].

    See Also
    --------
    vadcp_b_beam2ins : Beam-to-instrument transform for VADCP-B.
    vadcp_b_ins2earth : Instrument-to-Earth transform for VADCP-B.
    """
    # force shapes of some inputs to arrays of the correct dimensions
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # Compute the beam to instrument transform. Output two identical "Z"
    # beams to replicate the Nortek provided matlab code for transforming
    # beams, allowing the use of the Nortek provided transformation matrix
    # without any changes
    x, y, z1, z2 = vadcp_b_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4, tm)

    # compute the instrument to earth beam transform
    u, v, _ = vadcp_b_ins2earth(x, y, z1, z2, h, p, r)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    u_cor, _ = magnetic_correction(theta, u, v)

    # return the eastward velocity profile
    return u_cor


def vadcp_b_beam_northward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r,
                           lat, lon, dt, tm):
    """
    OOI wrapper for VELPROF-VLN_L1 from VADCP-B beam-coordinate data.

    Returns the northward velocity profile (VELPROF-VLN_L1) for the VADCP-B
    (Nortek Signature 55) as defined in DPS 1341-00750.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELPROF-B1_L0)
        [m s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELPROF-B2_L0)
        [m s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELPROF-B3_L0)
        [m s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELPROF-B4_L0)
        [m s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [degrees].
    p : array_like
        Instrument pitch [degrees].
    r : array_like
        Instrument roll [degrees].
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1900-01-01] (NTP time).
    tm : array_like
        Instrument-specific transformation matrix [4 x 4].

    Returns
    -------
    v_cor : ndarray
        Northward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELPROF-VLN_L1) [m s^-1].

    See Also
    --------
    vadcp_b_beam2ins : Beam-to-instrument transform for VADCP-B.
    vadcp_b_ins2earth : Instrument-to-Earth transform for VADCP-B.
    """
    # force shapes of some inputs to arrays of the correct dimensions
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    dt = np.atleast_1d(dt)

    # Compute the beam to instrument transform. Output two identical "Z"
    # beams to replicate the Nortek provided matlab code for transforming
    # beams, allowing the use of the Nortek provided transformation matrix
    # without any changes
    x, y, z1, z2 = vadcp_b_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4, tm)

    # compute the instrument to earth beam transform
    u, v, _ = vadcp_b_ins2earth(x, y, z1, z2, h, p, r)

    # compute the magnetic variation, and ...
    theta = magnetic_declination(lat, lon, dt)

    # ... correct for it
    _, v_cor = magnetic_correction(theta, u, v)

    # return the northward velocity profile
    return v_cor


def vadcp_b_beam_vertical_est(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r,
                              tm):
    """
    OOI wrapper for VELTURB-VLU-4BM_L1 from VADCP-B beam-coordinate data.

    Returns the estimated upward velocity profile (VELTURB-VLU-4BM_L1) for
    the VADCP-B (Nortek Signature 55) as defined in DPS 1341-00760. This
    product is the traditional estimate of vertical velocity derived from the
    4-beam transform, where each beam is angled at 25 degrees from vertical.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0)
        [m s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0)
        [m s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0)
        [m s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0)
        [m s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [degrees].
    p : array_like
        Instrument pitch [degrees].
    r : array_like
        Instrument roll [degrees].
    tm : array_like
        Instrument-specific transformation matrix [4 x 4].

    Returns
    -------
    w : ndarray
        Estimated upward velocity profiles in Earth coordinates
        (VELTURB-VLU-4BM_L1) [m s^-1].

    See Also
    --------
    vadcp_b_beam2ins : Beam-to-instrument transform for VADCP-B.
    vadcp_b_ins2earth : Instrument-to-Earth transform for VADCP-B.
    """
    # Compute the beam to instrument transform. Output two identical "Z"
    # beams to replicate the Nortek provided matlab code for transforming
    # beams, allowing the use of the Nortek provided transformation matrix
    # without any changes
    x, y, z1, z2 = vadcp_b_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4, tm)

    # compute the instrument to earth beam transform
    _, _, w = vadcp_b_ins2earth(x, y, z1, z2, h, p, r)

    # return the estimated Upward Velocity Profile
    return w


def vadcp_b_beam_vertical_true(b1, b2, b3, b4, b5, pg1, pg2, pg3, pg4, pg5,
                               h, p, r, tm):
    """
    OOI wrapper for VELTURB-VLU-5BM_L1 from VADCP-B beam-coordinate data.

    Returns the true upward velocity profile (VELTURB-VLU-5BM_L1) for the
    VADCP-B (Nortek Signature 55) as defined in DPS 1341-00760. This product
    uses beam 5 (pointing directly upward) for a better estimate of the true
    vertical velocity component. Depth cells where beam 5 percent good falls
    below 50% are set to the fill value.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0)
        [m s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0)
        [m s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0)
        [m s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0)
        [m s^-1].
    b5 : array_like
        Beam 5 velocity profiles in beam coordinates (VELTURB-B5_L0)
        [m s^-1].
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
        Instrument uncorrected magnetic heading [degrees].
    p : array_like
        Instrument pitch [degrees].
    r : array_like
        Instrument roll [degrees].
    tm : array_like
        Instrument-specific transformation matrix [4 x 4].

    Returns
    -------
    w : ndarray
        True upward velocity profiles in Earth coordinates
        (VELTURB-VLU-5BM_L1) [m s^-1].

    See Also
    --------
    vadcp_b_beam2ins : Beam-to-instrument transform for VADCP-B.
    vadcp_b_ins2earth : Instrument-to-Earth transform for VADCP-B.
    """
    # compute the beam to instrument transform
    x, y, _, _ = vadcp_b_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4, tm)

    # check percent good data for beam 5, reset to fill if less than 50%
    # (vendor provided "percent good" value)
    b5 = np.ma.filled(np.ma.masked_where(pg5 < 50, b5), VADCP_B_FILLVALUE)

    # compute the instrument to earth beam transform
    _, _, w = vadcp_b_ins2earth(x, y, b5, b5, h, p, r)

    # return the true Upward Velocity Profile
    return w


def vadcp_b_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4, tm):
    """
    Convert VADCP-B beam velocities to instrument coordinates.

    Applies the vendor-supplied transformation matrix to convert beam
    coordinate velocity profiles to the instrument coordinate system for the
    VADCP-B (Nortek Signature 55) as defined in DPS 1341-00750. Returns two
    vertical components (z1, z2) to replicate the Nortek MATLAB reference
    implementation.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates [m s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates [m s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates [m s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates [m s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    tm : array_like
        Instrument-specific transformation matrix [4 x 4].

    Returns
    -------
    x : ndarray
        X-axis velocity profiles in instrument coordinates [m s^-1].
    y : ndarray
        Y-axis velocity profiles in instrument coordinates [m s^-1].
    z1 : ndarray
        Z-axis velocity profiles derived from beams 1 and 3 [m s^-1].
    z2 : ndarray
        Z-axis velocity profiles derived from beams 2 and 4 [m s^-1].

    Notes
    -----
    The three-beam solution threshold for the VADCP-B is 50% (versus 25%
    for TRDI instruments), following the vendor-specified percent good
    floor.
    """
    # raw beam velocities, set to correct shape
    b1 = np.atleast_2d(b1)
    b2 = np.atleast_2d(b2)
    b3 = np.atleast_2d(b3)
    b4 = np.atleast_2d(b4)

    # adjust tm array to be 4 x 4 array
    tm = np.reshape(tm, (-1, 4, 4))

    # get number of elements and number of bins
    n_elements = b1.shape[0]
    n_bins = b1.shape[1]

    # percentage of good pings for each beam per depth cell, set to correct
    # shape
    pg1 = np.atleast_2d(pg1)
    pg2 = np.atleast_2d(pg2)
    pg3 = np.atleast_2d(pg3)
    pg4 = np.atleast_2d(pg4)

    # using the vendor specified percent good floor of 50%, create masked
    # arrays with fill values set to compute a 3-beam solution, if applicable.
    ma1 = np.ma.masked_where(pg1 < 50, b1)
    bm1 = ma1.filled((b3 - b4 - b2) * -1)
    ma2 = np.ma.masked_where(pg2 < 50, b2)
    bm2 = ma2.filled(b1 + b3 - b4)
    ma3 = np.ma.masked_where(pg3 < 50, b3)
    bm3 = ma3.filled((b1 - b4 - b2) * -1)
    ma4 = np.ma.masked_where(pg4 < 50, b4)
    bm4 = ma4.filled(b1 + b3 - b2)

    # sum across the masked arrays to determine if more than 1 beam is bad
    # per depth cell, if so we cannot compute a 3-beam solution and need to
    # set the fill value to a NaN.
    mad = np.ma.dstack((ma1, ma2, ma3, ma4))  # stack the masked arrays
    mas = np.ma.count_masked(mad, axis=2)     # count masked depth cells

    # using the above, reset the raw beams. fill with 3-beam if applicable,
    # otherwise use a NaN
    bm1 = np.ma.filled(np.ma.masked_where(mas > 1, bm1), VADCP_B_FILLVALUE)
    bm2 = np.ma.filled(np.ma.masked_where(mas > 1, bm2), VADCP_B_FILLVALUE)
    bm3 = np.ma.filled(np.ma.masked_where(mas > 1, bm3), VADCP_B_FILLVALUE)
    bm4 = np.ma.filled(np.ma.masked_where(mas > 1, bm4), VADCP_B_FILLVALUE)

    # Create a numpy array with dimensions of number of elements, number of
    # beams (4), and number of bins
    beams = np.zeros((n_elements, 4, n_bins))

    # pack the beam velocities into the appropriate slices.
    beams[:, 0, :] = bm1
    beams[:, 1, :] = bm2
    beams[:, 2, :] = bm3
    beams[:, 3, :] = bm4

    # the Einstein summation is here configured to do the matrix
    # multiplication xyz(h,i,k) = tm(i,j) * beams(h,j,k) on each slice h.
    xyz = np.einsum('hij,hjk->hik', tm, beams)

    # break out the xyz velocities and return them
    x = xyz[:, 0, :]
    y = xyz[:, 1, :]
    z1 = xyz[:, 2, :]
    z2 = xyz[:, 3, :]

    return x, y, z1, z2


def vadcp_b_ins2earth(u, v, w1, w2, heading, pitch, roll):
    """
    Convert VADCP-B instrument velocities to Earth coordinates.

    Applies the Nortek-convention heading offset and an extended 4x4 rotation
    matrix to convert instrument-frame velocity profiles to Earth coordinates
    for the VADCP-B (Nortek Signature 55) as defined in DPS 1341-00750. The
    final vertical velocity is the average of the two Earth-rotated vertical
    components.

    Parameters
    ----------
    u : array_like
        X-axis velocity profiles in instrument coordinates [m s^-1].
    v : array_like
        Y-axis velocity profiles in instrument coordinates [m s^-1].
    w1 : array_like
        Vertical velocity profiles derived from beams 1 and 3 [m s^-1].
    w2 : array_like
        Vertical velocity profiles derived from beams 2 and 4 [m s^-1].
    heading : array_like
        Instrument uncorrected magnetic heading [degrees].
    pitch : array_like
        Instrument pitch [degrees].
    roll : array_like
        Instrument roll [degrees].

    Returns
    -------
    uu : ndarray
        Eastward velocity profiles in Earth coordinates [m s^-1].
    vv : ndarray
        Northward velocity profiles in Earth coordinates [m s^-1].
    ww : ndarray
        Upward velocity profiles in Earth coordinates, averaged from
        the two vertical components [m s^-1].

    Notes
    -----
    The heading is offset by -90 degrees following the Nortek reference
    implementation.
    """
    # check for CI fill values before changing units.
    # this function 'conditions' (np.atleast_1d) its inputs.
    # TRDI does not apply its ADCP fill/bad value sentinels to compass data.
    heading, pitch, roll = replace_fill_with_nan(None, heading, pitch, roll)

    # heading, subtracting 90 degrees from heading based on Nortek provided
    # code
    H = heading - 90

    # roll
    Rrad = np.radians(roll)
    cos_R = np.cos(Rrad)
    sin_R = np.sin(Rrad)
    # heading
    Hrad = np.radians(H)
    cos_H = np.cos(Hrad)
    sin_H = np.sin(Hrad)
    # pitch
    Prad = np.radians(pitch)
    cos_P = np.cos(Prad)
    sin_P = np.sin(Prad)

    # determine array size
    n_packets = u.shape[0]
    n_uvw = u.shape[1]

    # initialize vectors to be used as matrix elements
    ones = np.ones(n_packets)
    zeros = ones * 0.0

    # heading and pitch matrices, from Nortek provided code
    Hmat = np.array([[cos_H, sin_H, zeros],
                     [-sin_H, cos_H, zeros],
                     [zeros, zeros, ones]])
    HmatRoll = np.rollaxis(Hmat, 2)

    Pmat = np.array([[cos_P, -sin_P*sin_R, -cos_R*sin_P],
                     [zeros, cos_R, -sin_R],
                     [sin_P, sin_R*cos_P, cos_P*cos_R]])
    PmatRoll = np.rollaxis(Pmat, 2)

    # Create a rotation matrix using the heading and pitch matrices
    R = np.einsum('hij,hjk->hik', HmatRoll, PmatRoll)
    # Extend the rotation matrix to a 4 x 4 matrix, to properly fit with the
    # 4 velocities
    R_ext = np.pad(R, ((0, 0), (0, 1), (0, 1)), 'constant',
                   constant_values=0)
    R_ext[:, -1, 0:3] = R[:, -1, :]  # Fill new row with values from last row
    R_ext[:, 0:3, -1] = R[:, :, -1]  # Fill new col with values from last col
    R_ext[:, 3, 2] = 0
    R_ext[:, 2, 3] = 0
    R_ext[:, 3, 3] = R_ext[:, 2, 2]

    # construct input array of coordinates (velocities) to be transformed.
    # the basis set is 3D (E,N,U) so that the middle dimension is sized at 3.
    uvw = np.zeros((n_packets, 4, n_uvw))

    # pack the coordinates (velocities) to be transformed into the appropriate
    # slices.
    uvw[:, 0, :] = u
    uvw[:, 1, :] = v
    uvw[:, 2, :] = w1
    uvw[:, 3, :] = w2

    # the Einstein summation is here configured to do the matrix
    # multiplication uvw_earth(i,m) = R_ext(h,i,l) * uvw(h,l,m) on each
    # slice h.
    uvw_earth = np.einsum('hil,hlm->him', R_ext, uvw)

    # break out the coordinate slices and return them
    uu = uvw_earth[:, 0, :]
    vv = uvw_earth[:, 1, :]
    ww1 = uvw_earth[:, 2, :]
    ww2 = uvw_earth[:, 3, :]
    ww = (ww1 + ww2) / 2

    return uu, vv, ww


# Compute the VELTURB_L1 data products for the VADCP instrument deployed by
# RSN.
def vadcp_beam_eastward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf,
                        lat, lon, dt):
    """
    OOI wrapper for VELTURB-VLE_L1 from original VADCP (TRDI) beam data.

    Returns the eastward velocity profile (VELTURB-VLE_L1) for the original
    VADCP (Teledyne RDI Workhorse Sentinel) as defined in DPS 1341-00760.
    For VADCP-B deployments, use vadcp_b_beam_eastward instead.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1900-01-01] (NTP time).

    Returns
    -------
    u_cor : ndarray
        Eastward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELTURB-VLE_L1) [m s^-1].

    See Also
    --------
    vadcp_b_beam_eastward : Equivalent function for VADCP-B deployments.
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


def vadcp_beam_northward(b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf,
                         lat, lon, dt):
    """
    OOI wrapper for VELTURB-VLN_L1 from original VADCP (TRDI) beam data.

    Returns the northward velocity profile (VELTURB-VLN_L1) for the original
    VADCP (Teledyne RDI Workhorse Sentinel) as defined in DPS 1341-00760.
    For VADCP-B deployments, use vadcp_b_beam_northward instead.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).
    lat : array_like
        Deployment latitude [decimal degrees].
    lon : array_like
        Deployment longitude [decimal degrees].
    dt : array_like
        Sample timestamp [seconds since 1900-01-01] (NTP time).

    Returns
    -------
    v_cor : ndarray
        Northward velocity profiles in Earth coordinates, corrected for
        magnetic declination (VELTURB-VLN_L1) [m s^-1].

    See Also
    --------
    vadcp_b_beam_northward : Equivalent function for VADCP-B deployments.
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
    Compute VELTURB-VLU-4BM_L1 for the original VADCP (TRDI 5-beam).

    Returns the estimated upward velocity profile (VELTURB-VLU-4BM_L1) from
    the 4-beam transform for the original VADCP (Teledyne RDI Workhorse
    Sentinel), where each beam is angled at 20 degrees from vertical. This
    product is the traditional estimate of vertical velocity using the
    standard beam-to-Earth transform as defined in DPS 1341-00760.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].
    h : array_like
        Instrument uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).

    Returns
    -------
    w : ndarray
        Estimated upward velocity profiles in Earth coordinates
        (VELTURB-VLU-4BM_L1) [m s^-1].
    """
    # compute the beam to instrument transform
    x, y, z, _ = adcp_beam2ins(b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    # compute the instrument to earth beam transform
    _, _, w = adcp_ins2earth(x, y, z, h, p, r, vf)

    # scale upward velocity to m/s
    w = w / 1000.  # mm/s -> m/s

    # return the estimated Upward Velocity Profile
    return w


def vadcp_beam_vertical_true(b1, b2, b3, b4, b5, pg1, pg2, pg3, pg4, pg5,
                             h, p, r, vf):
    """
    Compute VELTURB-VLU-5BM_L1 for the original VADCP (TRDI 5-beam).

    Returns the true upward velocity profile (VELTURB-VLU-5BM_L1) using
    beam 5 (pointing directly upward) for the original VADCP (Teledyne RDI
    Workhorse Sentinel), as defined in DPS 1341-00760. Depth cells where
    beam 5 percent good falls below 25% are set to the fill value.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0)
        [mm s^-1].
    b5 : array_like
        Beam 5 velocity profiles in beam coordinates (VELTURB-B5_L0)
        [mm s^-1].
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
        Instrument uncorrected magnetic heading [cdegrees].
    p : array_like
        Instrument pitch [cdegrees].
    r : array_like
        Instrument roll [cdegrees].
    vf : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).

    Returns
    -------
    w : ndarray
        True upward velocity profiles in Earth coordinates
        (VELTURB-VLU-5BM_L1) [m s^-1].
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


def vadcp_beam_error(b1, b2, b3, b4, pg1, pg2, pg3, pg4):
    """
    OOI wrapper for VELTURB-ERR_L1 from original VADCP (TRDI) beam data.

    Returns the error velocity profile (VELTURB-ERR_L1) for the original
    VADCP (Teledyne RDI Workhorse Sentinel) as defined in DPS 1341-00760.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates (VELTURB-B1_L0)
        [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates (VELTURB-B2_L0)
        [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates (VELTURB-B3_L0)
        [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates (VELTURB-B4_L0)
        [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].

    Returns
    -------
    e : ndarray
        Error velocity profiles (VELTURB-ERR_L1) [m s^-1].
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
    Compute relative echo intensity (ECHOINT_L1) from raw counts.

    Converts raw echo intensity from counts to dB using a factory-supplied
    scale factor as defined in DPS 1341-00750. The nominal scale factor is
    0.45 dB/count for the Workhorse family.

    Parameters
    ----------
    raw : array_like
        Raw echo intensity (ECHOINT_L0) [count].
    sfactor : float or array_like, optional
        Factory-supplied scale factor, instrument and beam specific
        [dB count^-1]. Default is 0.45.

    Returns
    -------
    dB : ndarray
        Relative echo intensity (ECHOINT_L1) [dB].

    Notes
    -----
    The ADCP outputs raw echo intensity as a 1-byte integer, so the
    ADCP_FILLVALUE (-32768, requiring 2 bytes) cannot appear in this data.
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
    Convert TRDI ADCP beam velocities to instrument coordinates.

    Applies the beam-to-instrument transformation matrix to convert beam
    coordinate velocity profiles to the instrument coordinate system for
    TRDI Workhorse ADCP instruments as defined in DPS 1341-00750.

    Parameters
    ----------
    b1 : array_like
        Beam 1 velocity profiles in beam coordinates [mm s^-1].
    b2 : array_like
        Beam 2 velocity profiles in beam coordinates [mm s^-1].
    b3 : array_like
        Beam 3 velocity profiles in beam coordinates [mm s^-1].
    b4 : array_like
        Beam 4 velocity profiles in beam coordinates [mm s^-1].
    pg1 : array_like
        Percent good estimate for beam 1 [percent].
    pg2 : array_like
        Percent good estimate for beam 2 [percent].
    pg3 : array_like
        Percent good estimate for beam 3 [percent].
    pg4 : array_like
        Percent good estimate for beam 4 [percent].

    Returns
    -------
    x : ndarray
        X-axis velocity profiles in instrument coordinates [mm s^-1].
    y : ndarray
        Y-axis velocity profiles in instrument coordinates [mm s^-1].
    z : ndarray
        Z-axis velocity profiles in instrument coordinates [mm s^-1].
    e : ndarray
        Error velocity profiles [mm s^-1].

    Notes
    -----
    A 3-beam solution is applied when one beam falls below 25% good,
    following the TRDI vendor specification (TRDI, 2010a, p. 14). Depth
    cells with more than one bad beam are set to NaN.
    """
    # raw beam velocities, set to correct shape
    b1 = np.atleast_2d(b1)
    b2 = np.atleast_2d(b2)
    b3 = np.atleast_2d(b3)
    b4 = np.atleast_2d(b4)

    # percentage of good pings for each beam per depth cell, set to correct
    # shape
    pg1 = np.atleast_2d(pg1)
    pg2 = np.atleast_2d(pg2)
    pg3 = np.atleast_2d(pg3)
    pg4 = np.atleast_2d(pg4)

    # using the vendor specified percent good floor of 25%, create masked
    # arrays with fill values set to compute a 3-beam solution, if applicable.
    ma1 = np.ma.masked_where(pg1 < 25, b1)
    bm1 = ma1.filled((b2 - b3 - b4) * -1)
    ma2 = np.ma.masked_where(pg2 < 25, b2)
    bm2 = ma2.filled((b1 - b3 - b4) * -1)
    ma3 = np.ma.masked_where(pg3 < 25, b3)
    bm3 = ma3.filled(b1 + b2 - b4)
    ma4 = np.ma.masked_where(pg4 < 25, b4)
    bm4 = ma4.filled(b1 + b2 - b3)

    # sum across the masked arrays to determine if more than 1 beam is bad
    # per depth cell, if so we cannot compute a 3-beam solution and need to
    # set the fill value to a NaN.
    mad = np.ma.dstack((ma1, ma2, ma3, ma4))  # stack the masked arrays
    mas = np.ma.count_masked(mad, axis=2)     # count masked depth cells

    # using the above, reset the raw beams. fill with 3-beam if applicable,
    # otherwise use a NaN
    bm1 = np.ma.filled(np.ma.masked_where(mas > 1, bm1), ADCP_FILLVALUE)
    bm2 = np.ma.filled(np.ma.masked_where(mas > 1, bm2), ADCP_FILLVALUE)
    bm3 = np.ma.filled(np.ma.masked_where(mas > 1, bm3), ADCP_FILLVALUE)
    bm4 = np.ma.filled(np.ma.masked_where(mas > 1, bm4), ADCP_FILLVALUE)

    bm1, bm2, bm3, bm4 = replace_fill_with_nan(
        ADCP_FILLVALUE, bm1, bm2, bm3, bm4)

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
    Convert TRDI ADCP instrument velocities to Earth coordinates.

    Applies the instrument-to-Earth rotation matrix to convert
    instrument-frame velocity profiles to Earth coordinates for TRDI
    Workhorse ADCP instruments as defined in DPS 1341-00750.

    Parameters
    ----------
    u : array_like
        X-axis velocity profiles in instrument coordinates [mm s^-1].
    v : array_like
        Y-axis velocity profiles in instrument coordinates [mm s^-1].
    w : array_like
        Z-axis velocity profiles in instrument coordinates [mm s^-1].
    heading : array_like
        Instrument uncorrected magnetic heading [cdegrees].
    pitch : array_like
        Instrument pitch [cdegrees].
    roll : array_like
        Instrument roll [cdegrees].
    vertical : array_like
        Instrument vertical orientation (0 = downward looking,
        1 = upward looking).

    Returns
    -------
    uu : ndarray
        Eastward velocity profiles in Earth coordinates [mm s^-1].
    vv : ndarray
        Northward velocity profiles in Earth coordinates [mm s^-1].
    ww : ndarray
        Upward velocity profiles in Earth coordinates [mm s^-1].

    Notes
    -----
    For upward-looking instruments, 180 degrees is added to the measured
    roll angle before applying the rotation matrix (DPS 1341-00750,
    Equation 3). Heading, pitch, and roll are converted from cdegrees to
    degrees internally.
    """
    # check for CI fill values before changing units.
    # this function 'conditions' (np.atleast_1d) its inputs.
    # TRDI does not apply its ADCP fill/bad value sentinels to compass data.
    heading, pitch, roll, vertical = replace_fill_with_nan(
        None, heading, pitch, roll, vertical)

    # change units from centidegrees to degrees
    heading = heading / 100.0
    pitch = pitch / 100.0
    roll = roll / 100.0

    # better way to calculate roll from the vertical orientation toggle;
    # this will propagate R as nans if the vertical variable is missing from
    # the data.
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

    # break out the coordinate slices and return them
    uu = uvw_earth[:, 0, :]
    vv = uvw_earth[:, 1, :]
    ww = uvw_earth[:, 2, :]

    return uu, vv, ww


def magnetic_correction(theta, u, v):
    """
    Correct velocity profiles for magnetic declination.

    Rotates horizontal velocity profiles by the magnetic declination angle
    to convert from magnetic to true Earth coordinates as defined in DPS
    1341-00750.

    Parameters
    ----------
    theta : array_like
        Magnetic declination at the measurement location [degrees].
        Positive values indicate magnetic north is east of true north.
    u : array_like
        Uncorrected eastward velocity profiles in Earth coordinates.
    v : array_like
        Uncorrected northward velocity profiles in Earth coordinates.

    Returns
    -------
    u_cor : ndarray
        Eastward velocity profiles corrected for magnetic declination.
    v_cor : ndarray
        Northward velocity profiles corrected for magnetic declination.

    Notes
    -----
    This function handles vectorized input where theta is a 1D array of
    shape (i,) and u, v are 2D arrays of shape (i, j), applying one
    declination value per data packet.
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


def adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure,
                        adcp_orientation, latitude):
    """
    Compute ADCP bin center depths from pressure in bar.

    Calculates center bin depths for PD0 and PD12 ADCP data using pressure
    input in bar as defined in DPS 1341-00750.

    Parameters
    ----------
    dist_first_bin : array_like
        Distance to the first ADCP bin [centimeters].
    bin_size : array_like
        Depth of each ADCP bin [centimeters].
    num_bins : array_like
        Number of ADCP bins [unitless].
    pressure : array_like
        Pressure at the sensor head [bar].
    adcp_orientation : array_like
        Instrument orientation (1 = upward looking, 0 = downward looking).
    latitude : array_like
        Deployment latitude [decimal degrees].

    Returns
    -------
    bin_depths : ndarray
        Center depths of each ADCP bin [meters].

    See Also
    --------
    adcp_bin_depths_meters : Core bin depth calculation from depth in meters.
    """
    # check for CI fill values.
    pressure = replace_fill_with_nan(None, pressure)

    # Convert pressure from bar to decibar
    pressure_dbar = pressure * 10.0

    # Calculate sensor depth using TEOS-10 toolbox z_from_p function
    # note change of sign to make the sensor_depth variable positive
    sensor_depth = -z_from_p(pressure_dbar, latitude)

    return adcp_bin_depths_meters(
        dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)


def adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure,
                         adcp_orientation, latitude):
    """
    Compute ADCP bin center depths from pressure in decaPascals.

    Calculates center bin depths for PD0 and PD12 ADCP data using pressure
    input in decaPascals as defined in DPS 1341-00750.

    Parameters
    ----------
    dist_first_bin : array_like
        Distance to the first ADCP bin [centimeters].
    bin_size : array_like
        Depth of each ADCP bin [centimeters].
    num_bins : array_like
        Number of ADCP bins [unitless].
    pressure : array_like
        Pressure at the sensor head [daPa].
    adcp_orientation : array_like
        Instrument orientation (1 = upward looking, 0 = downward looking).
    latitude : array_like
        Deployment latitude [decimal degrees].

    Returns
    -------
    bin_depths : ndarray
        Center depths of each ADCP bin [meters].

    See Also
    --------
    adcp_bin_depths_meters : Core bin depth calculation from depth in meters.
    """
    # check for CI fill values.
    pressure = replace_fill_with_nan(None, pressure)

    # Convert pressure from decaPascal to decibar
    pressure_dbar = pressure / 1000.0

    # Calculate sensor depth using TEOS-10 toolbox z_from_p function
    # note change of sign to make the sensor_depth variable positive
    sensor_depth = -z_from_p(pressure_dbar, latitude)

    return adcp_bin_depths_meters(
        dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)


def adcp_bin_depths(blanking_distance, bin_size, number_bins, orientation,
                    depth):
    """
    Compute ADCP bin center depths from sensor depth in meters.

    Calculates center bin depths for ADCP data from sensor depth in meters
    as defined in DPS 1341-00750.

    Parameters
    ----------
    blanking_distance : array_like
        Distance to the first ADCP bin [centimeters].
    bin_size : array_like
        Size (cell length) of each ADCP bin [centimeters].
    number_bins : array_like
        Number of ADCP bins [unitless].
    orientation : int
        Instrument orientation (1 = upward looking, 0 = downward looking).
    depth : array_like
        Depth of the sensor head [meters].

    Returns
    -------
    bin_depths : ndarray
        Center depths of each ADCP bin [meters].
    """
    # Convert from cm to meters
    blanking_distance = blanking_distance / 100.0
    bin_size = bin_size / 100.0

    # Following the PD0 convention, where:
    #     orientation = 0 is downward looking, bin depths are added to sensor
    #                     depth
    #                 = 1 is upward looking, bin depths are subtracted from
    #                     sensor depth
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
    Compute height from sea pressure using a 75-term TEOS-10 expression.

    Calculates height from sea pressure using the computationally-efficient
    75-term expression for density (Roquet et al., 2015). Calls
    enthalpy_SSO_0_p assuming standard ocean salinity and 0 degrees Celsius.

    Parameters
    ----------
    p : array_like
        Sea pressure [dbar].
    lat : array_like
        Latitude [decimal degrees north, -90 to +90].
    geo_strf_dyn_height : array_like, optional
        Dynamic height anomaly computed with p_ref = 0 [m^2 s^-2].
        Default is 0.
    sea_surface_geopotential : array_like, optional
        Geopotential at zero sea pressure [m^2 s^-2]. Default is 0.

    Returns
    -------
    z : ndarray
        TEOS-10 height [m]. Height is negative in the ocean (depth is the
        absolute value).

    Notes
    -----
    This function is a local implementation of the TEOS-10 z_from_p
    algorithm included before the gsw library was available in the OOI
    processing environment. It is a candidate for replacement with
    gsw.z_from_p in a future refactor.
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
    Compute enthalpy at standard ocean salinity and 0 degC as a function of p.

    Calculates enthalpy at the Standard Ocean Salinity (SSO) and at a
    Conservative Temperature of zero degrees Celsius as a function of
    pressure, using a streamlined version of the 76-term expression for
    specific volume (TEOS-10 ver 3.05).

    Parameters
    ----------
    p : array_like
        Sea pressure [dbar].

    Returns
    -------
    enthalpy_SSO_0 : ndarray
        Enthalpy at SSO and 0 degC [J kg^-1].

    Notes
    -----
    This function is a local implementation of a TEOS-10 subroutine
    called by z_from_p. Both are candidates for replacement with
    gsw.z_from_p in a future refactor.
    """
    z = p * 1e-4

    h006 = -2.1078768810e-9
    h007 = 2.8019291329e-10

    dynamic_enthalpy_SSO_0_p = z * (
        9.726613854843870e-4 + z * (
            -2.252956605630465e-5 + z * (
                2.376909655387404e-6 + z * (
                    -1.664294869986011e-7 + z * (
                        -5.988108894465758e-9 + z * (h006 + h007 * z))))))

    enthalpy_SSO_0 = dynamic_enthalpy_SSO_0_p * 1.e8  # db2Pa*1e4

    return enthalpy_SSO_0


def adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth,
                           adcp_orientation):
    """
    Compute ADCP bin center depths from sensor depth in meters.

    Calculates center bin depths for PD0, PD8, and PD12 ADCP data from
    sensor depth in meters as defined in DPS 1341-00750.

    Parameters
    ----------
    dist_first_bin : array_like
        Distance to the first ADCP bin [centimeters].
    bin_size : array_like
        Depth of each ADCP bin [centimeters].
    num_bins : array_like
        Number of ADCP bins [unitless].
    sensor_depth : array_like
        Estimated depth at the sensor head [meters].
    adcp_orientation : array_like
        Instrument orientation (1 = upward looking, 0 = downward looking).

    Returns
    -------
    bin_depths_pd8 : ndarray
        Center depths of each ADCP bin [meters].

    Notes
    -----
    The PD8 output format does not record dist_first_bin, bin_size, or
    num_bins; those values must be supplied externally for PD8 data.
    """
    # check for CI fill values.
    dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation = (
        replace_fill_with_nan(
            None, dist_first_bin, bin_size, num_bins, sensor_depth,
            adcp_orientation))

    # note, there is a CI problem not yet addressed if the time-vectorized
    # values in num_bins are not all the same!! For now, assume they are all
    # the same:
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
    #     adcp_orientation = 0 is downward looking, bindepths are added to
    #                          sensor depth
    #                      = 1 is upward looking, bindepths are subtracted
    #                          from sensor depth
    z_sign = 1.0 - 2.0 * adcp_orientation

    # to broadcast the vertical time dimension correctly with the horizontal
    # bin_numbers dimension, make all the 1D time arrays into column vectors
    # to be processed with the bin_numbers row vector.
    sensor_depth = sensor_depth.reshape(-1, 1)
    z_sign = z_sign.reshape(-1, 1)
    dist_first_bin = dist_first_bin.reshape(-1, 1)
    bin_size = bin_size.reshape(-1, 1)

    # Calculate bin depths
    bin_depths_pd8 = (sensor_depth
                      + z_sign * (dist_first_bin + bin_size * bin_numbers))

    return bin_depths_pd8


def depth_from_pressure_dbar(pressure, latitude, pressure_scale_factor=None):
    """
    Compute depth from pressure in dbar with optional scale factor.

    Calculates depth from pressure in decibar using the TEOS-10 z_from_p
    function, with an optional scale factor to convert from other pressure
    units to dbar.

    Parameters
    ----------
    pressure : array_like
        Pressure [dbar].
    latitude : array_like
        Deployment latitude [decimal degrees].
    pressure_scale_factor : float, optional
        Scale factor to convert pressure to dbar [unitless]. If None,
        pressure is used as-is.

    Returns
    -------
    depth : ndarray
        Depth [meters].
    """
    # check for CI fill values.
    pressure = replace_fill_with_nan(None, pressure)

    # Apply scale factor to convert pressure to decibar
    pressure_dbar = (pressure * pressure_scale_factor
                     if pressure_scale_factor else pressure)

    # Calculate sensor depth using TEOS-10 toolbox z_from_p function
    # note change of sign to make the sensor_depth variable positive
    return -z_from_p(pressure_dbar, latitude)
