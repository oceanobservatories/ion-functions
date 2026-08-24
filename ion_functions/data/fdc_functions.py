#!/usr/bin/env python
"""
ion_functions.data.fdc_functions

Module containing data processing functions for the FDCHP (Flux
Direct Covariance High Power) instrument family.
"""

import numpy as np
import scipy as sp
from scipy import integrate
from scipy import interpolate
from scipy import signal

from ion_functions import deprecated


"""
#...................................................................................
#...................................................................................
    The FDCHP instrument outputs a dataset of 10 Hz data for approximately 20 minutes
        every hour, resulting in around 12000 data records for each dataset.
    The FDCHP data product algorithms parse incoming 1D array data into separate
        datasets and operate on these datasets individually.
    Each of the four L1 products (wind and temperature) consists of 11400 values for
        each 20 minute dataset.
    Each of the three L2 flux products consists of one value per 20 minute dataset.

    Two additional auxiliary data products, not specified in the DPS, have been coded
        to provide time bases for the L1 and L2 products.
#...................................................................................
#...................................................................................

    LISTING OF SUBROUTINES BY ORDER IN THIS MODULE
        Grouped by sections; alphabetical within each section.

#...................................................................................
#...................................................................................
    Functions to compute the L1 FDCHP data products:

        fdc_tmpatur:        TMPATUR
        fdc_windtur_north:  WINDTUR-VLN
        fdc_windtur_up:     WINDTUR-VLU
        fdc_windtur_west:   WINDTUR-VLW
#...................................................................................
#...................................................................................
    Functions to compute the L2 FDCHP data products:

        fdc_fluxhot:            FLUXHOT
        fdc_fluxmom_alongwind:  FLUXMOM-U
        fdc_fluxmom_crosswind:  FLUXMOM-V
#...................................................................................
#...................................................................................
    Functions to compute the auxiliary time base data products:

        fdc_time_L1:  TIME_L1-AUX
        fdc_time_L2:  TIME_L2-AUX
#...................................................................................
#...................................................................................
    Primary routine to directly compute L1 wind products and L2 flux products:

        fdc_flux_and_wind
#...................................................................................
#...................................................................................
    Subroutines called by the primary routine fdc_flux_and_wind and its subroutines:

        fdc_accelsclimode
        fdc_alignwind
        fdc_anglesclimodeyaw
        fdc_despikesimple
        fdc_detrend
        fdc_filtcoef
        fdc_grv
        fdc_process_compass_data
        fdc_quantize_data
        fdc_sonic
        fdc_trans
        fdc_update
#...................................................................................
#...................................................................................

"""
####################################################################################
####################################################################################
####################################################################################
"""
#...................................................................................
#...................................................................................
    Functions to compute the L1 FDCHP data products:

        fdc_tmpatur:        TMPATUR
        fdc_windtur_north:  WINDTUR-VLN
        fdc_windtur_up:     WINDTUR-VLU
        fdc_windtur_west:   WINDTUR-VLW
#...................................................................................
#...................................................................................
"""


@deprecated
def fdc_tmpatur(timestamp, sonicT):
    """
    OOI single-product wrapper for TMPATUR_L1, sonic temperature [degC].

    Independently reproduces the sonic temperature conversion also
    computed internally by `fdc_flux_and_wind`; does not call that
    function directly.

    Parameters
    ----------
    timestamp : ndarray
        Data date and time values, TIME [seconds since 1900-01-01].
    sonicT : ndarray
        Raw speed of sound measured by the sonic anemometer,
        TMPATUR_L0 [counts].

    Returns
    -------
    Ts : ndarray
        Sonic temperature, TMPATUR_L1 [degC].

    See Also
    --------
    fdc_flux_and_wind : Core algorithm; computes TMPATUR_L1 internally
        as part of the full motion-correction pipeline.
    """
    # condition data and parse it into discrete datasets.
    data = fdc_quantize_data(timestamp, sonicT)
    # the shape of the data array is [n_var, n_packets, pts per packet].

    # there is no despiking or filtering to be done on these data, so truncate now.
    # number of seconds of data to remove from the beginning and end of the processed
    # data before calculating the mean of the products of the elements of two vectors.
    edge_sec = 30
    # sampling frequency
    fs = 10
    # number of edge data values to remove, based on sampling frequency fs in Hz
    edge = fs * edge_sec
    data = data[:, :, edge:-edge]

    # for clarity in following the DPS code, unpack the data into its constituent
    # variables as 2D arrays so that the index of the lead dimension of each
    # variable array indicates the dataset number; sonicT[0, :] will be a vector
    # containing the L0 temperature data for the first dataset packet.
    sonicT = data[1, :, :]

    # process L0 temperature data
    Ts = 0.01 * sonicT
    Ts = Ts * Ts / 403.0 - 273.15
    Ts = Ts.flatten()

    return Ts


@deprecated
def fdc_windtur_north(timestamp, sonicU, sonicV, sonicW, heading,
                      rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    OOI single-output wrapper for WINDTUR-VLN_L1, motion-corrected
    northward wind speed [m/s], uncorrected for magnetic variation.

    Parameters
    ----------
    timestamp : ndarray
        Data date and time values, TIME [seconds since 1900-01-01].
    sonicU : ndarray
        U-component windspeed in the buoy frame, WINDTUR-U_L0 [cm/s].
    sonicV : ndarray
        V-component windspeed in the buoy frame, WINDTUR-V_L0 [cm/s].
    sonicW : ndarray
        W-component windspeed in the buoy frame, WINDTUR-W_L0 [cm/s].
    heading : ndarray
        Magnetometer heading (yaw), MOTFLUX-YAW_L0 [radians].
    rateX : ndarray
        Roll angular rate from the gyro, MOTFLUX-ROLL_RATEX_L0
        [radians/s].
    rateY : ndarray
        Pitch angular rate from the gyro, MOTFLUX-PITCH_RATEY_L0
        [radians/s].
    rateZ : ndarray
        Yaw angular rate from the gyro, MOTFLUX-YAW_RATEZ_L0
        [radians/s].
    accX : ndarray
        X-component platform linear acceleration, MOTFLUX-ACX_L0
        [counts; 1 count = 9.80665 m/s^2].
    accY : ndarray
        Y-component platform linear acceleration, MOTFLUX-ACY_L0
        [counts; 1 count = 9.80665 m/s^2].
    accZ : ndarray
        Z-component platform linear acceleration, MOTFLUX-ACZ_L0
        [counts; 1 count = 9.80665 m/s^2].
    lat : float or ndarray
        Latitude of the instrument [decimal degrees].

    Returns
    -------
    wind_north : ndarray
        Motion-corrected northward wind speed, WINDTUR-VLN_L1 [m/s].

    Notes
    -----
    Pitch and roll are not used in any FDCHP data product calculation;
    only heading (yaw) and the gyro angular rates are used, per the
    DPS.

    See Also
    --------
    fdc_flux_and_wind : Core algorithm; computes all three WINDTUR_L1
        components and all three L2 flux products in a single call.
    """
    # this data product is temperature independent
    sonicT = sonicW * np.nan
    _, windspeeds = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                      heading, rateX, rateY, rateZ, accX, accY,
                                      accZ, lat)

    wind_north = np.asarray(windspeeds[0]).flatten()

    return wind_north


@deprecated
def fdc_windtur_up(timestamp, sonicU, sonicV, sonicW, heading,
                   rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    OOI single-output wrapper for WINDTUR-VLU_L1, motion-corrected
    upward wind speed [m/s].

    Parameters
    ----------
    timestamp : ndarray
        Data date and time values, TIME [seconds since 1900-01-01].
    sonicU : ndarray
        U-component windspeed in the buoy frame, WINDTUR-U_L0 [cm/s].
    sonicV : ndarray
        V-component windspeed in the buoy frame, WINDTUR-V_L0 [cm/s].
    sonicW : ndarray
        W-component windspeed in the buoy frame, WINDTUR-W_L0 [cm/s].
    heading : ndarray
        Magnetometer heading (yaw), MOTFLUX-YAW_L0 [radians].
    rateX : ndarray
        Roll angular rate from the gyro, MOTFLUX-ROLL_RATEX_L0
        [radians/s].
    rateY : ndarray
        Pitch angular rate from the gyro, MOTFLUX-PITCH_RATEY_L0
        [radians/s].
    rateZ : ndarray
        Yaw angular rate from the gyro, MOTFLUX-YAW_RATEZ_L0
        [radians/s].
    accX : ndarray
        X-component platform linear acceleration, MOTFLUX-ACX_L0
        [counts; 1 count = 9.80665 m/s^2].
    accY : ndarray
        Y-component platform linear acceleration, MOTFLUX-ACY_L0
        [counts; 1 count = 9.80665 m/s^2].
    accZ : ndarray
        Z-component platform linear acceleration, MOTFLUX-ACZ_L0
        [counts; 1 count = 9.80665 m/s^2].
    lat : float or ndarray
        Latitude of the instrument [decimal degrees].

    Returns
    -------
    wind_up : ndarray
        Motion-corrected upward wind speed, WINDTUR-VLU_L1 [m/s].

    Notes
    -----
    Pitch and roll are not used in any FDCHP data product calculation;
    only heading (yaw) and the gyro angular rates are used, per the
    DPS.

    See Also
    --------
    fdc_flux_and_wind : Core algorithm; computes all three WINDTUR_L1
        components and all three L2 flux products in a single call.
    """
    # this data product is temperature independent
    sonicT = sonicW * np.nan
    _, windspeeds = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                      heading, rateX, rateY, rateZ, accX, accY,
                                      accZ, lat)

    wind_up = np.asarray(windspeeds[2]).flatten()

    return wind_up


@deprecated
def fdc_windtur_west(timestamp, sonicU, sonicV, sonicW, heading,
                     rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    OOI single-output wrapper for WINDTUR-VLW_L1, motion-corrected
    westward wind speed [m/s], uncorrected for magnetic variation.

    Parameters
    ----------
    timestamp : ndarray
        Data date and time values, TIME [seconds since 1900-01-01].
    sonicU : ndarray
        U-component windspeed in the buoy frame, WINDTUR-U_L0 [cm/s].
    sonicV : ndarray
        V-component windspeed in the buoy frame, WINDTUR-V_L0 [cm/s].
    sonicW : ndarray
        W-component windspeed in the buoy frame, WINDTUR-W_L0 [cm/s].
    heading : ndarray
        Magnetometer heading (yaw), MOTFLUX-YAW_L0 [radians].
    rateX : ndarray
        Roll angular rate from the gyro, MOTFLUX-ROLL_RATEX_L0
        [radians/s].
    rateY : ndarray
        Pitch angular rate from the gyro, MOTFLUX-PITCH_RATEY_L0
        [radians/s].
    rateZ : ndarray
        Yaw angular rate from the gyro, MOTFLUX-YAW_RATEZ_L0
        [radians/s].
    accX : ndarray
        X-component platform linear acceleration, MOTFLUX-ACX_L0
        [counts; 1 count = 9.80665 m/s^2].
    accY : ndarray
        Y-component platform linear acceleration, MOTFLUX-ACY_L0
        [counts; 1 count = 9.80665 m/s^2].
    accZ : ndarray
        Z-component platform linear acceleration, MOTFLUX-ACZ_L0
        [counts; 1 count = 9.80665 m/s^2].
    lat : float or ndarray
        Latitude of the instrument [decimal degrees].

    Returns
    -------
    wind_west : ndarray
        Motion-corrected westward wind speed, WINDTUR-VLW_L1 [m/s].

    Notes
    -----
    Pitch and roll are not used in any FDCHP data product calculation;
    only heading (yaw) and the gyro angular rates are used, per the
    DPS.

    See Also
    --------
    fdc_flux_and_wind : Core algorithm; computes all three WINDTUR_L1
        components and all three L2 flux products in a single call.
    """
    # this data product is temperature independent
    sonicT = sonicW * np.nan
    _, windspeeds = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                      heading, rateX, rateY, rateZ, accX, accY,
                                      accZ, lat)

    wind_west = np.asarray(windspeeds[1]).flatten()

    return wind_west


"""
#...................................................................................
#...................................................................................
    Functions to compute the L2 FDCHP data products:

        fdc_fluxhot:            FLUXHOT
        fdc_fluxmom_alongwind:  FLUXMOM-U
        fdc_fluxmom_crosswind:  FLUXMOM-V
#...................................................................................
#...................................................................................
"""


@deprecated
def fdc_fluxhot(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    OOI single-output wrapper for FLUXHOT_L2, the sonic buoyancy flux
    [m/s K].

    Parameters
    ----------
    timestamp : ndarray
        Data date and time values, TIME [seconds since 1900-01-01].
    sonicU : ndarray
        U-component windspeed in the buoy frame, WINDTUR-U_L0 [cm/s].
    sonicV : ndarray
        V-component windspeed in the buoy frame, WINDTUR-V_L0 [cm/s].
    sonicW : ndarray
        W-component windspeed in the buoy frame, WINDTUR-W_L0 [cm/s].
    sonicT : ndarray
        Raw speed of sound measured by the sonic anemometer,
        TMPATUR_L0 [counts].
    heading : ndarray
        Magnetometer heading (yaw), MOTFLUX-YAW_L0 [radians].
    rateX : ndarray
        Roll angular rate from the gyro, MOTFLUX-ROLL_RATEX_L0
        [radians/s].
    rateY : ndarray
        Pitch angular rate from the gyro, MOTFLUX-PITCH_RATEY_L0
        [radians/s].
    rateZ : ndarray
        Yaw angular rate from the gyro, MOTFLUX-YAW_RATEZ_L0
        [radians/s].
    accX : ndarray
        X-component platform linear acceleration, MOTFLUX-ACX_L0
        [counts; 1 count = 9.80665 m/s^2].
    accY : ndarray
        Y-component platform linear acceleration, MOTFLUX-ACY_L0
        [counts; 1 count = 9.80665 m/s^2].
    accZ : ndarray
        Z-component platform linear acceleration, MOTFLUX-ACZ_L0
        [counts; 1 count = 9.80665 m/s^2].
    lat : float or ndarray
        Latitude of the instrument [decimal degrees].

    Returns
    -------
    fluxhot : ndarray
        Sonic buoyancy flux, FLUXHOT_L2, one value per 20-minute
        dataset [m/s K].

    See Also
    --------
    fdc_flux_and_wind : Core algorithm; computes all three L2 flux
        products and all three WINDTUR_L1 components in a single
        call.
    """
    fluxes, _ = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                  heading, rateX, rateY, rateZ, accX, accY,
                                  accZ, lat)

    fluxhot = fluxes[2]

    return fluxhot


@deprecated
def fdc_fluxmom_alongwind(timestamp, sonicU, sonicV, sonicW, heading,
                          rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    OOI single-output wrapper for FLUXMOM-U_L2, the along-wind
    component of the momentum flux [m^2/s^2].

    Parameters
    ----------
    timestamp : ndarray
        Data date and time values, TIME [seconds since 1900-01-01].
    sonicU : ndarray
        U-component windspeed in the buoy frame, WINDTUR-U_L0 [cm/s].
    sonicV : ndarray
        V-component windspeed in the buoy frame, WINDTUR-V_L0 [cm/s].
    sonicW : ndarray
        W-component windspeed in the buoy frame, WINDTUR-W_L0 [cm/s].
    heading : ndarray
        Magnetometer heading (yaw), MOTFLUX-YAW_L0 [radians].
    rateX : ndarray
        Roll angular rate from the gyro, MOTFLUX-ROLL_RATEX_L0
        [radians/s].
    rateY : ndarray
        Pitch angular rate from the gyro, MOTFLUX-PITCH_RATEY_L0
        [radians/s].
    rateZ : ndarray
        Yaw angular rate from the gyro, MOTFLUX-YAW_RATEZ_L0
        [radians/s].
    accX : ndarray
        X-component platform linear acceleration, MOTFLUX-ACX_L0
        [counts; 1 count = 9.80665 m/s^2].
    accY : ndarray
        Y-component platform linear acceleration, MOTFLUX-ACY_L0
        [counts; 1 count = 9.80665 m/s^2].
    accZ : ndarray
        Z-component platform linear acceleration, MOTFLUX-ACZ_L0
        [counts; 1 count = 9.80665 m/s^2].
    lat : float or ndarray
        Latitude of the instrument [decimal degrees].

    Returns
    -------
    fluxmom_along : ndarray
        Along-wind component of momentum flux, FLUXMOM-U_L2, one
        value per 20-minute dataset [m^2/s^2].

    See Also
    --------
    fdc_flux_and_wind : Core algorithm; computes all three L2 flux
        products and all three WINDTUR_L1 components in a single
        call.
    """
    # this data product is temperature independent
    sonicT = sonicW * np.nan
    fluxes, _ = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                  heading, rateX, rateY, rateZ, accX, accY,
                                  accZ, lat)

    fluxmom_along = fluxes[0]

    return fluxmom_along


@deprecated
def fdc_fluxmom_crosswind(timestamp, sonicU, sonicV, sonicW, heading,
                          rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    OOI single-output wrapper for FLUXMOM-V_L2, the cross-wind
    component of the momentum flux [m^2/s^2].

    Parameters
    ----------
    timestamp : ndarray
        Data date and time values, TIME [seconds since 1900-01-01].
    sonicU : ndarray
        U-component windspeed in the buoy frame, WINDTUR-U_L0 [cm/s].
    sonicV : ndarray
        V-component windspeed in the buoy frame, WINDTUR-V_L0 [cm/s].
    sonicW : ndarray
        W-component windspeed in the buoy frame, WINDTUR-W_L0 [cm/s].
    heading : ndarray
        Magnetometer heading (yaw), MOTFLUX-YAW_L0 [radians].
    rateX : ndarray
        Roll angular rate from the gyro, MOTFLUX-ROLL_RATEX_L0
        [radians/s].
    rateY : ndarray
        Pitch angular rate from the gyro, MOTFLUX-PITCH_RATEY_L0
        [radians/s].
    rateZ : ndarray
        Yaw angular rate from the gyro, MOTFLUX-YAW_RATEZ_L0
        [radians/s].
    accX : ndarray
        X-component platform linear acceleration, MOTFLUX-ACX_L0
        [counts; 1 count = 9.80665 m/s^2].
    accY : ndarray
        Y-component platform linear acceleration, MOTFLUX-ACY_L0
        [counts; 1 count = 9.80665 m/s^2].
    accZ : ndarray
        Z-component platform linear acceleration, MOTFLUX-ACZ_L0
        [counts; 1 count = 9.80665 m/s^2].
    lat : float or ndarray
        Latitude of the instrument [decimal degrees].

    Returns
    -------
    fluxmom_cross : ndarray
        Cross-wind component of momentum flux, FLUXMOM-V_L2, one
        value per 20-minute dataset [m^2/s^2].

    See Also
    --------
    fdc_flux_and_wind : Core algorithm; computes all three L2 flux
        products and all three WINDTUR_L1 components in a single
        call.
    """
    # this data product is temperature independent
    sonicT = sonicW * np.nan
    fluxes, _ = fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT,
                                  heading, rateX, rateY, rateZ, accX, accY,
                                  accZ, lat)

    fluxmom_cross = fluxes[1]

    return fluxmom_cross


"""
#...................................................................................
#...................................................................................
    Functions to compute the auxiliary time base data products:

        fdc_time_L1:  TIME_L1-AUX
        fdc_time_L2:  TIME_L2-AUX
#...................................................................................
#...................................................................................
"""


@deprecated
def fdc_time_L1(timestamp):
    """
    Computes TIME_L1-AUX, the timestamps associated with the L1 wind
    and temperature data products.

    Parameters
    ----------
    timestamp : ndarray
        Input data date and time values, TIME [seconds since
        1900-01-01].

    Returns
    -------
    time_L1 : ndarray
        Timestamps associated with the L1 data products,
        TIME_L1-AUX [seconds since 1900-01-01].

    Notes
    -----
    Not specified in the DPS. The L1 data products discard 30 seconds
    of data from the beginning and end of each 20-minute, 12000-record
    dataset; TIME_L1-AUX provides the corresponding timestamps for the
    remaining records, since the DPS does not otherwise document how
    timestamps are associated with the L1 arrays.
    """
    # condition data and parse it into discrete datasets.
    data = fdc_quantize_data(timestamp)
    # the shape of the data array is [n_var, n_packets, pts per packet].

    # number of seconds of data to remove from the beginning and end of the processed data
    edge_sec = 30
    # sampling frequency
    fs = 10
    # number of edge data values to remove, based on sampling frequency fs in Hz
    edge = fs * edge_sec
    data = data[:, :, edge:-edge]

    # for clarity in following the DPS code, unpack the data into its constituent
    # variables as 2D arrays so that the index of the lead dimension of each
    # variable array indicates the dataset number.
    tmstmp = data[0, :, :]

    time_L1 = tmstmp.flatten()

    return time_L1


@deprecated
def fdc_time_L2(timestamp):
    """
    Computes TIME_L2-AUX, the timestamps associated with the L2 flux
    data products.

    Parameters
    ----------
    timestamp : ndarray
        Input data date and time values, TIME [seconds since
        1900-01-01].

    Returns
    -------
    time_L2 : ndarray
        Median timestamp of each 20-minute dataset, one value per
        dataset packet, TIME_L2-AUX [seconds since 1900-01-01].

    Notes
    -----
    Not specified in the DPS. Each L2 flux product is a single value
    per 20-minute dataset; TIME_L2-AUX provides the corresponding
    (rounded, median) timestamp for each dataset, since the DPS does
    not otherwise document how timestamps are associated with the L2
    values.
    """
    # condition data and parse it into discrete datasets.
    data = fdc_quantize_data(timestamp)
    # the shape of the data array is [n_var, n_packets, pts per packet].

    # number of seconds of data to remove from the beginning and end of the processed
    # data before calculating the mean of the products of the elements of two vectors.
    edge_sec = 30
    # sampling frequency
    fs = 10
    # number of edge data values to remove, based on sampling frequency fs in Hz
    edge = fs * edge_sec
    data = data[:, :, edge:-edge]

    # for clarity in following the DPS code, unpack the data into its constituent
    # variables as 2D arrays so that the index of the lead dimension of each
    # variable array indicates the dataset number.
    tmstmp = data[0, :, :]

    # round the median L1 time values to get the L2 flux timestamps
    time_L2 = np.around(np.median(tmstmp, axis=-1))

    return time_L2


"""
#...................................................................................
#...................................................................................
    Primary routine to directly compute L1 wind products and L2 flux products:

        fdc_flux_and_wind
#...................................................................................
#...................................................................................
"""


@deprecated
def fdc_flux_and_wind(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                      rateX, rateY, rateZ, accX, accY, accZ, lat):
    """
    Computes the three L2 flux data products and the three L1 wind
    velocity components for the FDCHP instrument using the
    motion-corrected direct covariance (MCDC) method.

    Parses the input 1D arrays into discrete 20-minute, 12000-record
    datasets, estimates platform Euler angles from a complementary
    filter combining gyro-integrated angular rate and accelerometer
    tilt reference, computes platform velocity by rotating and
    integrating the linear accelerations, assembles the true wind
    velocity relative to Earth, rotates it into the streamwise
    coordinate system, and computes the kinematic momentum and
    buoyancy fluxes from the detrended fluctuations.

    Parameters
    ----------
    timestamp : ndarray
        Data date and time values, TIME [seconds since 1900-01-01].
    sonicU : ndarray
        U-component windspeed in the buoy frame, WINDTUR-U_L0 [cm/s].
    sonicV : ndarray
        V-component windspeed in the buoy frame, WINDTUR-V_L0 [cm/s].
    sonicW : ndarray
        W-component windspeed in the buoy frame, WINDTUR-W_L0 [cm/s].
    sonicT : ndarray
        Raw speed of sound measured by the sonic anemometer,
        TMPATUR_L0 [counts].
    heading : ndarray
        Magnetometer heading (yaw), MOTFLUX-YAW_L0 [radians].
    rateX : ndarray
        Roll angular rate from the gyro, MOTFLUX-ROLL_RATEX_L0
        [radians/s].
    rateY : ndarray
        Pitch angular rate from the gyro, MOTFLUX-PITCH_RATEY_L0
        [radians/s].
    rateZ : ndarray
        Yaw angular rate from the gyro, MOTFLUX-YAW_RATEZ_L0
        [radians/s].
    accX : ndarray
        X-component platform linear acceleration, MOTFLUX-ACX_L0
        [counts; 1 count = 9.80665 m/s^2].
    accY : ndarray
        Y-component platform linear acceleration, MOTFLUX-ACY_L0
        [counts; 1 count = 9.80665 m/s^2].
    accZ : ndarray
        Z-component platform linear acceleration, MOTFLUX-ACZ_L0
        [counts; 1 count = 9.80665 m/s^2].
    lat : float or ndarray
        Latitude of the instrument [decimal degrees].

    Returns
    -------
    fluxes : tuple of ndarray
        Three-element tuple of L2 flux products, one value per
        dataset: (fluxmom_u, fluxmom_v, fluxhot). fluxmom_u is the
        along-wind momentum flux, FLUXMOM-U_L2 [m^2/s^2]; fluxmom_v
        is the cross-wind momentum flux, FLUXMOM-V_L2 [m^2/s^2];
        fluxhot is the buoyancy flux, FLUXHOT_L2 [m/s K].
    windspeeds : tuple of list of ndarray
        Three-element tuple of L1 windspeed products, one array per
        dataset per element: (windtur_vln, windtur_vlw, windtur_vlu),
        the motion-corrected North, West, and Up wind velocity
        components, WINDTUR-VLN_L1, WINDTUR-VLW_L1, WINDTUR-VLU_L1
        [m/s].

    Notes
    -----
    Pitch and roll are accepted as inputs but not used in any FDCHP
    data product calculation, per the DPS; the machinery to process
    them is retained for possible future use. The high-pass filter
    cutoff frequency is hardcoded to 1/12 Hz per e-mail guidance from
    DPS author Jim Edson, consistent with the DPS's own discussion of
    theoretical limitations (Section 3.3).
    """
    # uncertainty in how latitude will be broadcasted. so.
    lat = np.atleast_1d(lat)
    if lat.size == 1:
        lat = np.repeat(lat, sonicT.size)

    # condition data and parse it into discrete datasets.
    # the heading data is passed 2 extra times, once to take the place of pitch
    # data, once for roll, which in the original matlab code are processed but
    # not used to calculate any data products.
    data = fdc_quantize_data(timestamp, sonicU, sonicV, sonicW, sonicT, heading,
                             heading, heading, rateX, rateY, rateZ, accX, accY,
                             accZ, lat)

    # the shape of the data array is [n_var, n_packets, pts per packet].
    #print data.shape
    n_pack = data.shape[1]

    # for clarity in following the DPS code, unpack the data into its constituent
    # variables as 2D arrays so that the index of the lead dimension of each
    # variable array indicates the dataset number; sonicU[0, :] will be a vector
    # containing the U wind data for the first dataset packet.
    sonicU = data[1, :, :]
    sonicV = data[2, :, :]
    sonicW = data[3, :, :]
    sonicT = data[4, :, :]
    heading = data[5, :, :]
    roll = data[6, :, :]      # not used to calculate any data products
    pitch = data[7, :, :]     # not used to calculate any data products
    rateX = data[8, :, :]
    rateY = data[9, :, :]
    rateZ = data[10, :, :]
    accX = data[11, :, :]
    accY = data[12, :, :]
    accZ = data[13, :, :]
    lat = data[14, :, :]

    # pitch and roll aren't currently used. to emphasize this:
    roll = np.nan
    pitch = np.nan

    # calculate the gravitational acceleration for each dataset
    gv = fdc_grv(np.median(lat, axis=-1))

    # process L0 data
    sonicU = 0.01 * sonicU
    sonicV = 0.01 * sonicV
    sonicW = 0.01 * sonicW
    sonicT = 0.01 * sonicT
    sonicT = sonicT * sonicT / 403.0

    # convert IMU from N,E,Down to match Sonic N,W,Up coordinate system
    rateY = -rateY
    rateZ = -rateZ
    accY = -accY
    accZ = -accZ
    pitch = -pitch
    heading = -heading

    # hardcoded variables in the DPS code
    G = 9.80665  # units of accelerometer values
    roffset = 0.0
    poffset = 0.0
    # z distance between IMU and sonic sampling volume: ok to hardcode
    z_imu_2_smplvol = 0.753
    # distance vector between IMU and sonic sampling volume
    Rvec = np.array([0.0, 0.0, z_imu_2_smplvol])

    fs = 10                         # sampling frequency, Hz
    #fltr_cutoff_freq = 10.0        # cutoff frequency to generate filter coeffs
    # 16-oct-2014 e-mail from Jim Edson (DPS author): use
    fltr_cutoff_freq = 12.0

    bhi, ahi = fdc_filtcoef(fs, 1.0/fltr_cutoff_freq)

    # gyro is the processed compass data, and,
    # goodcompass is a switch signifying whether these data are good.
    # this subroutine is vectorized.
    gyro, goodcompass = fdc_process_compass_data(heading)

    # number of seconds of data to remove from the beginning and end of the processed
    # data before calculating the mean of the products of the elements of two vectors.
    edge_sec = 30
    # number of edge data values to remove, based on sampling frequency
    edge = int(fs * edge_sec)
    # set up sonic temperature for buoyancy flux calculation;
    # the temperature processing can be vectorized outside the loop
    Ts_L1 = sonicT[:, edge:-edge]
    Ts = fdc_detrend(Ts_L1, -1, 'linear')

    # initialize L2 dataproduct arrays
    fluxmom_u = np.zeros(n_pack)
    fluxmom_v = np.zeros(n_pack)
    fluxhot = np.zeros(n_pack)

    # and lists to contain the L1 dataproducts
    vln = [None] * n_pack
    vlw = [None] * n_pack
    vlu = [None] * n_pack
    Tmp = [None] * n_pack

    # process one datapacket at a time
    for ii in range(n_pack):
        # wind speeds
        sonics = np.vstack((sonicU[ii, :], sonicV[ii, :], sonicW[ii, :]))

        # process angular rate data; already in radians
        deg_rate = np.vstack((rateX[ii, :], rateY[ii, :], rateZ[ii, :]))
        deg_rate = fdc_despikesimple(deg_rate)

        # process the linear accelerometer data:
        platform = np.vstack((accX[ii, :], accY[ii, :], accZ[ii, :])) * G
        platform = fdc_despikesimple(platform)
        gcomp = np.mean(platform, axis=-1)
        g = np.array([np.sqrt(np.sum(gcomp*gcomp))])
        platform = platform * gv[ii]/g

        platform[0, :] = platform[0, :] + poffset
        platform[1, :] = platform[1, :] + roffset

        gcomp = np.mean(platform, axis=-1)
        g = np.array([np.sqrt(np.sum(gcomp*gcomp))])
        platform = platform * gv[ii] / g

        euler, dr = fdc_anglesclimodeyaw(ahi, bhi, fs, platform, deg_rate,
                                         gyro[ii, :], goodcompass[ii, 0])

        # euler angles are right-handed
        _, uvwplat, _ = fdc_accelsclimode(bhi, ahi, fs, platform, euler)

        uvw, _, _ = fdc_sonic(sonics, dr, euler, uvwplat, Rvec)

        UVW_L1 = uvw[:, edge:-edge]

        # rotate wind velocity components into windstream
        u = fdc_alignwind(UVW_L1)

        u = fdc_detrend(u, -1, 'linear')

        # calculate flux products
        fluxmom_u[ii] = np.mean(u[2, :] * u[0, :])
        fluxmom_v[ii] = np.mean(u[2, :] * u[1, :])
        fluxhot[ii] = np.mean(u[2, :] * Ts[ii])

        # save the L1 wind data products
        (vln[ii], vlw[ii], vlu[ii]) = (UVW_L1[0, :], UVW_L1[1, :], UVW_L1[2, :])

    fluxes = (fluxmom_u, fluxmom_v, fluxhot)

    windspeeds = (vln, vlw, vlu)

    return fluxes, windspeeds


"""
#...................................................................................
#...................................................................................
    Subroutines called by the primary routine fdc_flux_and_wind and its subroutines:

        fdc_accelsclimode
        fdc_alignwind
        fdc_anglesclimodeyaw
        fdc_despikesimple
        fdc_detrend
        fdc_filtcoef
        fdc_grv
        fdc_process_compass_data
        fdc_quantize_data
        fdc_sonic
        fdc_trans
        fdc_update
#...................................................................................
#...................................................................................
"""


def fdc_accelsclimode(bhi, ahi, sf, accm, euler):
    """
    Rotates linear accelerations measured on the FDCHP platform into
    an earth reference system, then integrates to get platform
    velocity and displacement.

    Parameters
    ----------
    bhi : ndarray
        Numerator coefficients for the high-pass filter.
    ahi : ndarray
        Denominator coefficients for the high-pass filter.
    sf : float
        Sampling frequency [Hz].
    accm : ndarray
        Measured platform linear accelerations, shape (3, n_rec)
        [m/s^2].
    euler : ndarray
        Euler angles phi, theta, psi, shape (3, n_rec) [radians].

    Returns
    -------
    acc : ndarray
        Linear accelerations in the earth reference frame, shape
        (3, n_rec) [m/s^2]. Not used downstream in the data product
        algorithm.
    uvwplat : ndarray
        Platform linear velocities at the point of measurement, shape
        (3, n_rec) [m/s].
    xyzplat : ndarray
        Platform displacements from the mean position, shape
        (3, n_rec) [m]. Not used downstream in the data product
        algorithm.

    Notes
    -----
    Follows the DPS Appendix A MATLAB code directly. Gravity is
    removed from the rotated vertical acceleration before integrating
    to velocity and displacement; both integration stages are high-
    pass filtered forward and backward (filtfilt) to remove drift.
    """
    # keep DPS variable names and comments
    gravxyz = np.mean(accm, axis=-1)
    gravity = np.sqrt(gravxyz.dot(gravxyz))

    # rotate measured accelerations into earth frame
    acc = fdc_trans(accm, euler)
    acc[2, :] = acc[2, :] - gravity

    # integrate accelerations to get velocities
    uvwplat = integrate.cumtrapz(acc, axis=-1, initial=0.0) / sf
    uvwplat = signal.filtfilt(bhi, ahi, uvwplat, axis=-1, padtype='odd', padlen=15)

    # integrate again to get displacements
    xyzplat = integrate.cumtrapz(uvwplat, axis=-1, initial=0.0) / sf
    xyzplat = signal.filtfilt(bhi, ahi, xyzplat, axis=-1, padtype='odd', padlen=15)

    return acc, uvwplat, xyzplat


def fdc_alignwind(u):
    """
    Rotates wind velocity components into the streamwise wind.

    Parameters
    ----------
    u : ndarray
        Wind velocity in the earth frame, uncorrected for magnetic
        declination, shape (3, n_rec) [m/s].

    Returns
    -------
    u_rot : ndarray
        Wind velocity rotated into the streamwise coordinate system,
        shape (3, n_rec) [m/s].

    Notes
    -----
    The rotation angles are computed from the mean wind vector so
    that the mean cross-wind and vertical components are forced to
    zero, per the DPS Appendix A MATLAB code. Implemented here as
    direct arithmetic rather than matrix multiplication for speed.
    """
    # mean wind velocity components
    u_mean = np.mean(u, axis=-1)

    # calculate angles for coordinate rotation
    u_hor = np.sqrt(u_mean[0] * u_mean[0] + u_mean[1] * u_mean[1])
    beta = np.arctan2(u_mean[2], u_hor)
    alpha = np.arctan2(u_mean[1], u_mean[0])

    # populate rotation matrix
    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    sin_b = np.sin(beta)
    cos_b = np.cos(beta)

    ur = u[0, :] * cos_a * cos_b + u[1, :] * sin_a * cos_b + u[2, :] * sin_b
    vr = -u[0, :] * sin_a + u[1, :] * cos_a
    wr = -u[0, :] * cos_a * sin_b - u[1, :] * sin_a * sin_b + u[2, :] * cos_b

    u_rot = np.vstack((ur, vr, wr))

    return u_rot


def fdc_anglesclimodeyaw(ahi, bhi, sf, accm, ratem, gyro, goodcompass):
    """
    Computes the Euler angles for the FDCHP instrument platform using
    complementary filtering of accelerometer tilt and integrated
    gyro angular rate.

    Parameters
    ----------
    ahi : ndarray
        Denominator coefficients for the high-pass filter.
    bhi : ndarray
        Numerator coefficients for the high-pass filter.
    sf : float
        Sampling frequency [Hz].
    accm : ndarray
        Recalibrated platform linear accelerations, shape (3, n_rec)
        [m/s^2].
    ratem : ndarray
        Recalibrated angular rates, shape (3, n_rec) [radians/s].
    gyro : ndarray
        Gyro signal (heading/compass), shape (1, n_rec) [radians].
    goodcompass : bool
        Whether the gyro (compass) measurements are reliable enough
        to use.

    Returns
    -------
    euler : ndarray
        Euler angles phi, theta, psi, shape (3, n_rec) [radians].
    dr : ndarray
        Detrended angular rate velocities, shape (3, n_rec)
        [radians/s].

    Notes
    -----
    Follows the DPS Appendix A MATLAB code. Slow (low-frequency) roll
    and pitch angles are derived from the accelerometer tilt; slow
    yaw is derived from the unwrapped gyro signal if goodcompass is
    True, otherwise held at the median gyro value. The slow angles
    are used as a first guess, then refined over 5 iterations by
    integrating and high-pass filtering the angular rates and adding
    the result to the slow angles.
    """
    # hardcoded number of iterations for integrations
    n_iterations = 5

    # keep DPS variable names and comments
    gravxyz = np.mean(accm, axis=-1)
    gravity = np.sqrt(np.sum(gravxyz*gravxyz))

    # unwrap compass
    gyro = np.unwrap(gyro)

    # DPS comment: "remove mean from rate sensors".
    # However, the DPS code removes the linear trend.
    # note that the matlab function detrend used with the 'linear'
    # option specified subtracts piecewise linear values, whereas
    # the scipy version acts to subtract least fit values (as does
    # the matlab version with no options specified).
    ratem = fdc_detrend(ratem, -1, 'linear')

    ### documentation from DPS code verbatim:
    # low frequency angles from accelerometers and gyro
    # slow roll from gravity effects on horizontal accelerations. low pass
    # filter since high frequency horizontal accelerations may be 'real'

    ## avoid imaginary numbers by making sure that arcsin operates on [-1 1].
    # RAD version of revised DPS pitch
    # also, trap out runtime warnings when theta has nans
    theta = np.array(-accm[0, :] / gravity, ndmin=2)
    nanmask = np.isnan(theta)
    theta[nanmask] = 0.0
    theta[theta <= -1.0] = -np.pi/2.0
    theta[theta >= 1.0] = np.pi/2.0
    mask = np.absolute(theta) < 1.0
    theta[mask] = np.arcsin(theta[mask])
    theta[nanmask] = np.nan
    theta_slow = theta - signal.filtfilt(bhi, ahi, theta, axis=-1, padtype='odd', padlen=15)

    # RAD version of revised DPS roll
    # nanmask used for the same reason as above
    phi = np.array(accm[1, :] / gravity, ndmin=2) / np.cos(theta_slow)
    nanmask = np.isnan(phi)
    phi[nanmask] = 0.0
    phi[phi <= -1.0] = -np.pi/2.0
    phi[phi >= 1.0] = np.pi/2.0
    mask = np.absolute(phi) < 1.0
    phi[mask] = np.arcsin(phi[mask])
    phi[nanmask] = np.nan
    phi_slow = phi - signal.filtfilt(bhi, ahi, phi, axis=-1, padtype='odd', padlen=15)

    ### documentation from DPS code verbatim:
    # yaw
    # here, we estimate the slow heading. the 'fast heading' is not needed
    # for the euler angle update matrix. the negative sign puts the gyro
    # signal into a right handed system.

    # these are from revised DPS code: fs = 10, cutoff_freq = 1/240
    #ahi2 = np.array([1.0, -4.989457431527359, 9.957885277614746, -9.936911062858101,
    #                 4.957996019344925, -0.989512802573847])
    #bhi2 = np.array([0.994742581059968, -4.973712905299841, 9.947425810599682, -9.947425810599682,
    #                 4.973712905299841, -0.994742581059968])

    # Using 1.0/240.0 as the second argument results in matrices approaching singularity
    # (computation on the edge of robustness).
    bhi2, ahi2 = fdc_filtcoef(sf, 1.0/240.0)

    if goodcompass:
        psi_slow = -gyro - signal.filtfilt(bhi2, ahi2, -gyro, axis=-1, padtype='odd', padlen=15)
    else:
        psi_slow = -np.median(gyro)*np.ones(phi.shape)

    # use slow angles as first guess
    euler = np.vstack((phi_slow, theta_slow, psi_slow))
    rates = fdc_update(ratem, euler)

    # "i will use this filter with a lower cutoff for yaw
    #  since the compass is having issues"

    # integrate and filter angle rates, and add to slow angles
    for ii in range(n_iterations):
        phi_int = integrate.cumtrapz(rates[0, :], axis=-1, initial=0.0) / sf
        phi = phi_slow + signal.filtfilt(bhi, ahi, phi_int, axis=-1, padtype='odd', padlen=15)
        theta_int = integrate.cumtrapz(rates[1, :], axis=-1, initial=0.0) / sf
        theta = theta_slow + signal.filtfilt(bhi, ahi, theta_int, axis=-1, padtype='odd', padlen=15)
        psi_int = integrate.cumtrapz(rates[2, :], axis=-1, initial=0.0) / sf
        # rad: note that psi_slow values are also a function of the goodcompass value
        if goodcompass:
            psi = psi_slow + signal.filtfilt(bhi2, ahi2, psi_int, axis=-1, padtype='odd', padlen=15)
        else:
            psi = psi_slow + psi_int

        euler = np.vstack((phi, theta, psi))
        rates = fdc_update(ratem, euler)
        rates = fdc_detrend(rates, -1, 'constant')

    dr = ratem

    return euler, dr


def fdc_despikesimple(data):
    """
    Removes outliers from FDCHP platform motion and angular rate
    data streams by iterative median/standard-deviation thresholding
    and nearest-neighbor interpolation.

    Parameters
    ----------
    data : ndarray
        Data values, shape (3, N).

    Returns
    -------
    data : ndarray
        Despiked data values, shape (3, N).

    Notes
    -----
    Matches the DPS Appendix A MATLAB code: 3 iterations, points
    outside median +/- 4 standard deviations are replaced by nearest-
    neighbor interpolation. Python's nearest-neighbor interpolation
    replaces from the left where MATLAB's interp1 replaces from the
    right for equidistant abscissae; the input is flipped before and
    after despiking to better match the MATLAB reference output.
    bounds_error is set to False (with fill_value=nan) to avoid
    ValueErrors when replacing an outlier would require extrapolation,
    which requires nanmedian/nanstd in place of median/std so that
    interp1d is never called with an empty array. This is a
    documented weak point in the original algorithm: it does not
    exclude low-frequency trend/variability from the spike threshold,
    so a Tukey-style approach may be more robust.
    """
    # at first the matlab test code flux results for the goodcompass=0 case differed
    # in the 4th decimal place compared to the python results; this agreement extends
    # out to the 11th decimal place if the data is flipped before and after this routine
    # is run. (the rest of the discrepancy results because filtfilt.m and scipy.filtfilt
    # do not give exactly the same results).

    # the reason flipping the python data increases agreement is because matlab
    # nearest neighbor interpolation replaces from the right while scipy replaces
    # from the left.
    data = np.fliplr(data)  # to match the matlab results

    # number of times to run each vector stream through the despiking routine
    n_iterations = 3
    # standard deviation span; this was 6 in the original DPS and revised code;
    # Jim Edson (DPS author) says to set this at 4
    n_std = 4
    # interpolation method; 'nearest' in test code
    ntrpmeth = 'nearest'
    #print "interpolation method", ntrpmeth

    array_size = np.shape(data)
    t = np.arange(0, array_size[1])

    for jj in range(n_iterations):

        # vectorize the median, stdev, and masking operations outside of the inner
        # loop, which should help program efficiency when processing large numbers
        # of datasets.

        # calculate the median and stdev as column vectors for broadcasting
        M = np.atleast_2d(np.nanmedian(data, axis=-1)).T
        # original code used matlab std function, which has a "N-1" in denominator -
        # so, ddof=1.
        Sn = np.nanstd(data, axis=-1, ddof=1, keepdims=True) * n_std
        mask = np.logical_and(data <= M + Sn, data >= M - Sn)
        # the interp1d function is vectorized for the second argument 2D array ONLY IF
        # the first argument 1D array is unchanging - which it's not.
        # therefore, here a for loop is required.
        for ii in range(array_size[0]):
            f = interpolate.interp1d(t[mask[ii, :]], data[ii, mask[ii, :]], kind=ntrpmeth,
                                     bounds_error=False, fill_value=np.nan)
            data[ii, :] = f(t)

        ## as coded in DPS
        #for tot in range(array_size[0]):
        #    M = np.median(data[tot, :])
        #    # original code used matlab std function, which has a "N-1" in denom, so ddof=1
        #    stan = np.std(data[tot, :], ddof=1)
        #    Sn = stan * n_std
        #    mask = np.logical_and(data[tot, :] < M + Sn, data[tot, :] > M - Sn)
        #    f = interpolate.interp1d(t[mask], data[tot, mask], kind=ntrpmeth)
        #    data[tot, :] = f(t)

    # re-orient the data to the way it came into the routine
    data = np.fliplr(data)

    return data


def fdc_detrend(data, axis=-1, type='linear'):
    """
    Detrends data using scipy.signal.detrend, trapping NaN values
    that would otherwise raise a runtime execution error.

    Parameters
    ----------
    data : ndarray
        Data to detrend.
    axis : int, optional
        Data axis along which detrend is applied. Default is the
        last axis (-1).
    type : str, optional
        Detrend method, either 'linear' or 'constant'. Default is
        'linear'.

    Returns
    -------
    detrended : ndarray
        Detrended data, same shape as `data`. If any input value is
        non-finite (NaN or inf), all output values are set to NaN,
        regardless of the dimensionality of `data` or the axis
        designated, since scipy.signal.detrend cannot process NaNs.

    Notes
    -----
    Not present in the DPS; a Python-only wrapper around
    scipy.signal.detrend needed to avoid runtime errors on NaN-
    containing input.
    """
    # trap out nans and infs
    if np.any(~np.isfinite(data)):
        detrended = np.zeros(data.shape) + np.nan
        return detrended
    # else run the data through the scipy function
    detrended = signal.detrend(data, axis=axis, type=type)
    return detrended


def fdc_filtcoef(fs, fc):
    """
    Computes Butterworth high-pass filter coefficients as a function
    of sampling frequency and cutoff frequency.

    Parameters
    ----------
    fs : float
        Sampling frequency [Hz].
    fc : float
        Inverse of the cutoff frequency used to derive the filter
        (e.g. 1.0 / fltr_cutoff_freq) [s].

    Returns
    -------
    bhigh : ndarray
        Numerator (b) Butterworth filter coefficients, for use with
        scipy.signal.filtfilt.
    ahigh : ndarray
        Denominator (a) Butterworth filter coefficients, for use with
        scipy.signal.filtfilt.

    Notes
    -----
    Not present in the DPS, which hardcodes the filter coefficients
    directly in Appendix A rather than computing them from a cutoff
    frequency; this function derives equivalent coefficients via
    scipy.signal.buttord and scipy.signal.butter. High-pass filter
    intended to retain real acceleration signal while removing drift;
    the current filter design (chosen so the doubly-integrated
    acceleration spectrum matches the frequency-domain integrated
    power spectrum in the pass band) supersedes an earlier 1998
    design.
    """
    nfreq = fs / 2.0
    wp = fc / nfreq
    ws = .7 * wp
    n, wn = signal.buttord(wp, ws, 10.0, 25.0)
    bhigh, ahigh = signal.butter(n, wn, 'high')
    return bhigh, ahigh


def fdc_grv(lat):
    """
    Computes gravitational acceleration as a function of latitude.

    Parameters
    ----------
    lat : float or ndarray
        Latitude of the instrument [decimal degrees].

    Returns
    -------
    g : ndarray
        Acceleration due to earth's gravitational field [m/s^2].

    Notes
    -----
    Evaluated as a polynomial in sin(lat)^2 using Horner's method,
    following the DPS Appendix A MATLAB code exactly (equatorial
    gravity gamma = 9.7803267715 m/s^2, plus latitude-dependent
    correction terms).
    """
    # constants from the DPS:
    # equatorial value for 'g'
    gamma = 9.7803267715
    # coefficients of polynomial in (sin(lat))^2
    c1 = 0.0052790414
    c2 = 0.0000232718
    c3 = 0.0000001262
    c4 = 0.0000000007

    x = sp.sin(np.radians(lat))
    xsq = x * x

    # Horner's method for calculating polynomials
    g = gamma * (1.0 + xsq * (c1 + xsq * (c2 + xsq * (c3 + xsq * c4))))

    ## straightforward powers method
    #g=gamma*(1.0+c1*x**2+c2*x**4+c3*x**6+c4*x**8)

    return g


def fdc_process_compass_data(heading):
    """
    Vectorized routine that processes raw compass (heading) data from
    the FDCHP instrument and determines whether it is reliable enough
    to use.

    Parameters
    ----------
    heading : ndarray
        Heading in a (N, W, Up) coordinate system, shape
        (n_packets, n_pts) [radians].

    Returns
    -------
    gyro : ndarray
        Processed compass data, same shape as `heading` [radians].
    goodcompass : ndarray of bool
        Per-packet switch denoting reliability of the heading data:
        False if the heading range exceeds 120 degrees or the
        standard deviation exceeds 45 degrees, True otherwise.

    Notes
    -----
    Not a standalone function in the DPS Appendix A MATLAB code,
    where this processing is inline in the main routine; factored out
    here for clarity and to vectorize across dataset packets. The
    first and last 10 samples of each packet are overwritten with
    their nearest interior neighbor before despiking, since edge
    compass readings are often bad.
    """
    # this routine is vectorized as noted above.
    # number of values on either edge of the compass readings to overwrite:
    edge_compass = 10

    # gyro = 'heading' = compass (may be in N,W,Up cpoordinate sytem); already in radians
    gyro = heading
    # overwrite edge values
    gyro[:, 0:edge_compass] = gyro[:, [edge_compass]]     # right side is a column vector
    gyro[:, -edge_compass:] = gyro[:, [-edge_compass-1]]  # right side is a column vector

    # process gyro values
    gx = np.cos(gyro)
    gy = np.sin(gyro)
    gx = fdc_despikesimple(gx)
    gy = fdc_despikesimple(gy)
    gyro = np.arctan2(gy, gx)
    gyro[gyro < 0] = gyro[gyro < 0] + 2.0 * np.pi

    # determine whether gyro data is good
    gchk = np.unwrap(gyro)
    # matlab std uses (N-1) in the denominator, so set ddof=1
    stdhdg = np.std(gchk, axis=-1, ddof=1, keepdims=True)
    hdg_range = np.amax(gchk, axis=-1, keepdims=True) - np.amin(gchk, axis=-1, keepdims=True)

    # set the goodcompass vector:
    #    if ( hdg_range>(120/180*pi) or stdhdg>(45/180*pi) )
    #        goodcompass = 0
    #    else
    #        goodcompass = 1
    #    end
    #
    # do the same thing without conditional
    goodcompass = np.logical_not(np.logical_or(
        hdg_range > (120.0/180.0*np.pi), stdhdg > (45.0/180.0*np.pi)))

    #print 'goodcompass value: ', goodcompass
    return gyro, goodcompass


def fdc_quantize_data(*args):
    """
    Groups FDCHP data into discrete 12000-record (20-minute) dataset
    packets by parsing record timestamps.

    Parameters
    ----------
    *args : ndarray
        Variable-length argument list of 1D input data arrays, all
        the same length (on the order of npts * n_packets). The
        first array must be the timestamps.

    Returns
    -------
    data : ndarray
        3D array of the input variables grouped into datasets, shape
        (n_variables, n_packets, 12000).

    Notes
    -----
    Not present in the DPS Appendix A code, which operates on a
    single pre-chunked 12000-point dataset. FDCHP collects 10 Hz data
    for 20 minutes every hour, so consecutive dataset packets are
    separated by roughly 40 minutes; packets are identified from gaps
    greater than 1800 seconds in the timestamp array. Per guidance
    from DPS author Jim Edson, packets with more than 12000 records
    are truncated to 12000, and packets with fewer than 12000 records
    are padded to 12000 using the last data record. Also works if the
    only input variable is a 1D array of timestamps.
    """
    # target number of datapoints per dataset: 10Hz for 20 minutes.
    npts = 12000
    # time in between dataset chunks (last of n_th and first of n+1_th)
    # is expected to be 40 minutes = 2400 seconds; use a lower value as
    # the time discriminant.
    time_gap = 1800.0

    data = np.atleast_2d(args)  # data is a 2D array

    # parse into discrete datasets by finding the number of datapoints in each chunk.
    # data are expected to come in chunks of 20 minutes duration, separated by 40 min.
    # the first row of data must be the timestamps.
    # first find the indices at these timegaps
    idx_at_gap = np.where(np.diff(data[0, :]) > time_gap)[0]

    # prepend and append values to get accurate counts for 1st and last dataset
    idx_at_gap = np.hstack((-1, idx_at_gap, data.shape[1]-1))

    # difference to get answer
    chunklengths = np.diff(idx_at_gap)
    #print chunklengths

    # process one dataset chunk at a time.
    # process all data streams for each dataset at the same time,
    # inspecting each chunk to see if it has (less than), (equal to),
    # or (greater than) npts and then processing it accordingly.
    for ii in range(chunklengths.size):
        if chunklengths[ii] < npts:
            # pad dataset with last set of datapoints.
            # here idx is the insertion point index (1st element = 1);
            # pad data and insert after this point and before start of next dataset
            idx = ii * npts + chunklengths[ii]
            # number of rows to insert
            n_nsrt = npts-chunklengths[ii]
            # now go to python array indexing conventions;
            # tile the column vector data[:,idx-1:idx]
            filldata = np.tile(data[:, idx-1:idx], (1, n_nsrt))
            # correct the timestamps of the filldata.
            delta_time = np.median(np.diff(data[0, idx-chunklengths[ii]:idx]))
            filldata[0, :] = filldata[0, :] + np.arange(1.0, n_nsrt+1) * delta_time
            # and insert filldata into the data array
            data = np.hstack((data[:, 0:idx], filldata, data[:, idx:]))
        elif chunklengths[ii] == npts:
            continue  # no action needed
        else:  # chunklengths[ii] > npts
            # delete data after npts in this dataset and before next dataset.
            idx_del_beg = (ii + 1) * npts
            idx_del_end = ii * npts + chunklengths[ii]
            data = np.delete(data, np.s_[idx_del_beg:idx_del_end], 1)

    # convert data to a 3D array so that calling program can parse its shape
    # to figure out dataset dimensions (n_var, n_dataset_packets, npts per dataset)
    data = np.reshape(data, (data.shape[0], -1, npts))

    return data


def fdc_sonic(sonics, omegam, euler, uvwplat, dist_vec):
    """
    Corrects the sonic anemometer wind velocity components for
    platform motion and orientation.

    Parameters
    ----------
    sonics : ndarray
        Sonic anemometer wind velocity components in the buoy frame,
        shape (3, N) [m/s].
    omegam : ndarray
        Measured angular rate in the platform frame, shape (3, N)
        [radians/s].
    euler : ndarray
        Euler angles phi, theta, psi, shape (3, N) [radians].
    uvwplat : ndarray
        Platform velocities, shape (3, N) [m/s].
    dist_vec : ndarray
        Distance vector between the IMU and the sonic anemometer
        sampling volume, shape (3,) [m].

    Returns
    -------
    uvw : ndarray
        Corrected sonic anemometer components in the fixed earth
        (North, West, Up) frame, shape (3, N) [m/s].
    uvwr : ndarray
        Intermediate rotated sonic components prior to adding
        platform velocity, shape (3, N) [m/s]. Not used downstream in
        the data product algorithm.
    uvwrot : ndarray
        Angular-rate cross-product correction term, shape (3, N)
        [m/s]. Not used downstream in the data product algorithm.

    Notes
    -----
    From the EDDYCORR toolbox, per the DPS Appendix A MATLAB code.
    """

    n_rec = euler.shape[-1]
    Rvec = np.transpose(np.tile(dist_vec, (n_rec, 1)))
    # override default cross product vector axis definition, which is -1
    uvwrot = np.cross(omegam, Rvec, axis=-2)

    uvwr = fdc_trans(sonics + uvwrot, euler)
    uvw = uvwr + uvwplat

    return uvw, uvwr, uvwrot


def fdc_trans(ang_rates, angles):
    """
    Rotates a 3-component vector (linear acceleration or wind
    velocity) from the platform frame into the earth frame of
    reference using the Euler angle rotation matrix.

    Parameters
    ----------
    ang_rates : ndarray
        3-component vector to rotate (e.g. linear accelerations or
        sonic wind velocity plus angular-rate cross-product term),
        shape (3, N).
    angles : ndarray
        Euler angles phi, theta, psi, shape (3, N) [radians].

    Returns
    -------
    values : ndarray
        Rotated vector in the earth frame, shape (3, N).

    Notes
    -----
    Implements the coordinate transformation matrix T(phi, theta,
    psi) = A(psi) A(theta) A(phi) defined in the DPS (Section 3.2.1,
    Equation 5).
    """

    p = angles[0, :]
    t = angles[1, :]
    ps = angles[2, :]

    up = ang_rates[0, :]
    vp = ang_rates[1, :]
    wp = ang_rates[2, :]

    u = (up * sp.cos(t) * sp.cos(ps) + vp * (sp.sin(p) * sp.sin(t) *
         sp.cos(ps) - sp.cos(p) * sp.sin(ps)) + wp * (sp.cos(p) *
         sp.sin(t) * sp.cos(ps) + sp.sin(p) * sp.sin(ps)))
    v = (up * sp.cos(t) * sp.sin(ps) + vp * (sp.sin(p) * sp.sin(t) *
         sp.sin(ps) + sp.cos(p) * sp.cos(ps)) + wp * (sp.cos(p) *
         sp.sin(t) * sp.sin(ps) - sp.sin(p) * sp.cos(ps)))
    w = (up * (-sp.sin(t)) + vp * (sp.cos(t) * sp.sin(p)) + wp *
         (sp.cos(t) * sp.cos(p)))

    values = np.vstack((u, v, w))

    return values


def fdc_update(ang_rates, angles):
    """
    Computes the angular update matrix relating measured
    strapped-down angular rates to the time derivatives of the Euler
    angles.

    Parameters
    ----------
    ang_rates : ndarray
        Angular rates in the platform frame, shape (3, N)
        [radians/s].
    angles : ndarray
        Euler angles phi, theta, psi, shape (3, N) [radians].

    Returns
    -------
    values : ndarray
        Euler angle rates (phi_dot, theta_dot, psi_dot), shape
        (3, N) [radians/s].

    Notes
    -----
    From the EDDYCORR toolbox. Implements the strapped-down angular
    rate relation given in the DPS (Section 3.2.2, Equation 8);
    matches Edson et al. (1998).
    """

    p = angles[0, :]
    t = angles[1, :]
    ps = angles[2, :]

    up = ang_rates[0, :]
    vp = ang_rates[1, :]
    wp = ang_rates[2, :]

    u = up + vp * sp.sin(p) * np.tan(t) + wp * sp.cos(p) * np.tan(t)
    v = 0 + vp * sp.cos(p) - wp * sp.sin(p)
    w = 0 + vp * sp.sin(p) / sp.cos(t) + wp * sp.cos(p) / sp.cos(t)

    values = np.vstack((u, v, w))
    return values
