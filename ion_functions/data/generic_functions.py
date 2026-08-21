#!/usr/bin/env python
"""
Generic data-calculation functions used across multiple OOI instrument
families. These functions have no dedicated instrument family or Data
Product Specification of their own; they are organized here because
more than one family's module needs them.

Author: Christopher Wingard, Stuart Pearce, Russell Desiderio
"""

# Common imports
import datetime
import pandas as pd
import ppigrf
import numpy as np
import numexpr as ne
import time
from numbers import Integral

from ion_functions import deprecated

# CyberInfrastructure fill value for all integer data types
SYSTEM_FILLVALUE = -999999999


def replace_fill_with_nan(instrument_fillvalue, *args):
    """
    Converts integer data types to float and replaces both system and
    instrument integer fill values with nan.

    Parameters
    ----------
    instrument_fillvalue : int, list, or ndarray, optional
        Instrument-specific integer fill value(s). Use None, an empty
        list, or an empty array if there are none.
    *args : array_like
        One or more variables to check for fill values. Integer arrays
        are converted to float with fill values replaced by nan; float
        arrays are passed through unchanged.

    Returns
    -------
    out : array_like or list of array_like
        The input variable(s) with fill values replaced by nan. A
        single array if only one variable was passed, otherwise a
        list of arrays.

    Notes
    -----
    Float fill values are np.nan, but integer types cannot hold nan;
    integers are stored as signed 32-bit values with a single system
    fill value (SYSTEM_FILLVALUE). Instrument fill values vary by
    instrument and are passed in via instrument_fillvalue. The
    'beams' variable in vel3dk is used as a dictionary key and should
    not be processed with this function.
    """
    # while all input arguments presented to the DPA functions by CI will be np.ndarrays of
    # at least one dimension, many of the unit tests were written before this policy was
    # established. as a result many use core python scalars to represent single-valued
    # parameters ( for example: lat = 45.0 instead of lat = np.array([45.0]) ). because of
    # this, for compatibility input arguments that are core python scalars will be supported.
    args = np.atleast_1d(*args)
    # the output of this last statement is a list, UNLESS there is only one element,
    # in which case instead of a list of one element the element itself is returned.

    # if args consists of a single element make it into a list so that when it is indexed,
    # the index refers to the element itself instead of the element's subelements.
    if not isinstance(args, list):
        args = [args]

    # the first fill value to be checked will be the system fillvalue
    all_fillvalues = np.array([SYSTEM_FILLVALUE])
    # in case instrument_fillvalue is passed as an ndarray, turn it into a flattened list;
    # testing a list for whether it is empty (F) or not (T) is more straightforward than
    # dealing with ndarrays of sizes from empty to multiple elements - some of my tests
    # with these revealed non-intuitive results when subjected to boolean testing.
    if isinstance(instrument_fillvalue, np.ndarray):
        instrument_fillvalue = list(instrument_fillvalue.flatten())
    if instrument_fillvalue is not None:
        all_fillvalues = np.hstack((all_fillvalues, instrument_fillvalue))

    ## original coding loops
    #for ii in range(len(args)):
    #    # check to see if the first element is an integer datatype
    #    if isinstance(np.ndarray.flatten(args[ii])[0], (long, int)):
    #        mask = np.zeros_like(args[ii], dtype=bool)
    #        for jj in range(len(all_fillvalues)):
    #            mask = np.logical_or(mask, args[ii] == all_fillvalues[jj])
    #        args[ii] = args[ii].astype('float')
    #        args[ii][mask] = np.nan

    # more pythonic, perhaps faster
    for ii, val in enumerate(args):
        # check to see if the first element is an integer datatype
        if isinstance(val.flatten()[0], Integral):
            mask = np.zeros_like(val, dtype=bool)
            for jj, fil in enumerate(all_fillvalues):
                mask = np.logical_or(mask, val == fil)
            val = val.astype('float')
            val[mask] = np.nan
            args[ii] = val

    # if args has only one element, unit tests fail if it is passed out as a list.
    if len(args) == 1:
        args = args[0]
    return args


def magnetic_declination(lat, lon, timestamp, z=0.0, zflag=-1, ntp=1):
    """
    Computes magnetic declination for a platform location and date,
    vectorizing per-sample calls to igrf_declination.

    Parameters
    ----------
    lat : array_like
        Latitude of the instrument [decimal degrees]. East is
        positive, West negative.
    lon : array_like
        Longitude of the instrument [decimal degrees]. North is
        positive, South negative.
    timestamp : array_like
        Time stamp from a data particle, either NTP [secs since
        1900-01-01] or Unix [secs since 1970-01-01].
    z : array_like, optional
        Depth or height of the instrument relative to sea level
        [meters]. Positive values only. Default 0.
    zflag : array_like, optional
        -1 to treat z as a depth (i.e. -z), 1 to treat z as a height
        (i.e. +z). Default -1.
    ntp : array_like, optional
        1 if timestamp is NTP time, 0 if Unix time. Default 1.

    Returns
    -------
    mag_dec : array_like
        Magnetic declination (magnetic variation) [degrees from N].
        Positive values are eastward, negative westward of North.

    See Also
    --------
    igrf_declination : Core per-sample algorithm.
    """
    decln = np.vectorize(igrf_declination)
    mag_dec = decln(lat, lon, timestamp, z, zflag, ntp)
    return mag_dec


def magnetic_correction(theta, u, v):
    """
    Corrects velocity vectors for magnetic variation (declination) at
    the measurement location. Used by several data products (e.g.
    VELPROF, WINDAVG) across multiple instrument classes.

    Parameters
    ----------
    theta : array_like
        Magnetic variation based on location and date [degrees]. See
        igrf_declination.
    u : array_like
        Uncorrected eastward velocity in Earth coordinates.
    v : array_like
        Uncorrected northward velocity in Earth coordinates.

    Returns
    -------
    u_cor : array_like
        Eastward velocity, in Earth coordinates, corrected for
        magnetic variation.
    v_cor : array_like
        Northward velocity, in Earth coordinates, corrected for
        magnetic variation.
    """
    theta_rad = np.radians(theta)
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)

    m = np.array([
        [cos_t, sin_t],
        [-1*sin_t, cos_t]
    ])

    u = np.atleast_1d(u)
    v = np.atleast_1d(v)
    cor = np.dot(m, np.array([u, v]))

    return cor[0], cor[1]


def igrf_declination(lat, lon, timestamp, z=0.0, zflag=-1, ntp=1):
    """
    Computes magnetic declination (magnetic variation) for a single
    location and date from the International Geomagnetic Reference
    Field (IGRF).

    Parameters
    ----------
    lat : float
        Latitude of the instrument [decimal degrees]. East is
        positive, West negative.
    lon : float
        Longitude of the instrument [decimal degrees]. North is
        positive, South negative.
    timestamp : float
        Time stamp from a data particle, either NTP [secs since
        1900-01-01] or Unix [secs since 1970-01-01].
    z : float, optional
        Depth or height of the instrument relative to sea level
        [meters]. Positive values only. Default 0.
    zflag : float, optional
        -1 to treat z as a depth (i.e. -z), 1 to treat z as a height
        (i.e. +z). Default -1.
    ntp : float, optional
        1 if timestamp is NTP time, 0 if Unix time. Default 1.

    Returns
    -------
    mag_dec : float
        Magnetic declination (magnetic variation) [degrees from N].
        Positive values are eastward, negative westward of North.

    Examples
    --------
    >>> igrf_declination(45.0, -128, 3574792037.958, 1000)
    16.45715213214582

    See Also
    --------
    magnetic_declination : Public, array-vectorized entry point.
    """

    if ntp == 1:
        timestamp = timestamp - 2208988800.

    date = pd.to_datetime(timestamp, unit='s')

    # set the depth to negative for below sealevel (if needed) and convert from
    # meters to kilometers.
    z = z / 1000.  # m -> km
    if z > 0 and zflag == -1:   # check that depth is a positive number first
        z = zflag * z    # convert z to indicate depth

    # calculate the magnetic declination
    be, bn, bu = ppigrf.igrf(lon, lat, z, date)
    _, mag_dec = ppigrf.get_inclination_declination(be, bn, bu)
    return mag_dec[0]


@deprecated
def ntp_to_unix_time(ntp_timestamp):
    """
    Converts an NTP time stamp (epoch 1900-01-01) to a Unix time
    stamp (epoch 1970-01-01).

    Parameters
    ----------
    ntp_timestamp : array_like
        NTP timestamp(s) [seconds since 1900-01-01].

    Returns
    -------
    unix_timestamp : array_like
        Unix timestamp(s) [seconds since 1970-01-01].

    Notes
    -----
    Deprecated: use the fixed offset of 2208988800 seconds between
    the NTP and Unix epochs directly instead of calling this
    function.
    """
    SYSTEM_EPOCH = datetime.date(*time.gmtime(0)[0:3])
    NTP_EPOCH = datetime.date(1900, 1, 1)
    NTP_DELTA = (SYSTEM_EPOCH - NTP_EPOCH).total_seconds()

    unix_timestamp = ne.evaluate('ntp_timestamp - NTP_DELTA')
    return unix_timestamp


def extract_parameter(in_array, index):
    """
    Extracts a single value from an array. Used, for example, to
    extract the L0 PH434SI value from the array holding the 24 sets
    of 4 light measurements made by the Sunburst SAMI-II pH
    instrument (PHSEN).

    Parameters
    ----------
    in_array : array_like
        The input array holding the value.
    index : int
        0-based index into in_array.

    Returns
    -------
    out_value : object
        The extracted value at index.
    """
    out_value = in_array[index]
    return out_value


def bilinear_interpolation(x, y, points):
    """
    Interpolates the value at (x, y) from four points forming a
    rectangle.

    Parameters
    ----------
    x : float
        x-coordinate to interpolate at.
    y : float
        y-coordinate to interpolate at.
    points : array_like
        Four (x, y, value) triplets, in any order, forming the
        corners of a rectangle.

    Returns
    -------
    value : float
        The bilinearly interpolated value at (x, y).

    Examples
    --------
    >>> bilinear_interpolation(12, 5.5,
    ...                        [(10, 4, 100),
    ...                         (20, 4, 200),
    ...                         (10, 6, 150),
    ...                         (20, 6, 300)])
    165.0

    Notes
    -----
    See https://en.wikipedia.org/wiki/Bilinear_interpolation for the
    formula.
    """
    # order points by x, then by y
    pts = np.sort(points.view('f8,f8,f8'), order=['f0', 'f1'])
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1) + 0.0)


def error(x, y):
    """
    Relative error of x with respect to y.

    Parameters
    ----------
    x : array_like
        Value(s) to compare.
    y : array_like
        Reference value(s).

    Returns
    -------
    out : array_like
        abs(x - y) / abs(y).

    Notes
    -----
    Internal helper used only by select_arg_within_tolerance_of_std.
    """
    return np.abs(x - y) / np.abs(y)


def select_non_zero_arg(a1, a2=None, a1_scale_factor=None, a2_scale_factor=None):
    """
    Tests a1 and a2 for non-zero values and returns the non-zero
    array, scaled if a scale factor is supplied.

    Parameters
    ----------
    a1 : array_like
        An input array to test for non-zero elements.
    a2 : array_like, optional
        An input array to test for non-zero elements.
    a1_scale_factor : array_like, optional
        Scale factor to apply to a1.
    a2_scale_factor : array_like, optional
        Scale factor to apply to a2.

    Returns
    -------
    out : array_like
        The scaled array containing at least one non-zero element.
    """
    if np.any(a1):
        return a1 * a1_scale_factor if np.any(a1_scale_factor) else a1
    if np.any(a2):
        return a2 * a2_scale_factor if np.any(a2_scale_factor) else a2
    return a1


def select_arg_within_tolerance_of_std(a1, a2=None, std=None, tol=0.25,
                                       a1_scale_factor=1, a2_scale_factor=1,
                                       std_scale_factor=1):
    """
    Tests a1 and a2 for non-zero values and returns the non-zero
    array whose average is within tolerance of the standard, scaling
    inputs as needed for comparability.

    Parameters
    ----------
    a1 : array_like
        An input array to test for non-zero elements and tolerance.
    a2 : array_like, optional
        An input array to test for non-zero elements and tolerance.
    std : array_like, optional
        The standard that a1 and a2 are compared against.
    tol : array_like, optional
        Relative tolerance the returned array's average must be
        within, compared to the standard. Default 0.25.
    a1_scale_factor : array_like, optional
        Scale factor to apply to a1. Default 1.
    a2_scale_factor : array_like, optional
        Scale factor to apply to a2. Default 1.
    std_scale_factor : array_like, optional
        Scale factor to apply to std. Default 1.

    Returns
    -------
    out : array_like
        The scaled array containing at least one non-zero element
        with an average within tolerance of the standard.

    See Also
    --------
    error : Relative-error helper used internally.
    """
    a1s = a1 * a1_scale_factor if np.any(a1) and np.any(a1_scale_factor) else a1
    a2s = a2 * a2_scale_factor if np.any(a2) and np.any(a2_scale_factor) else a2

    if not np.any(a2s) or not np.any(std) or not np.any(tol):
        return a1s
    if not np.any(a1s):
        return a2s

    stds = std * std_scale_factor if np.any(std_scale_factor) else std
    avg_tol = np.average(tol)

    if np.average(error(a1s, stds)) < avg_tol:
        return a1s
    if np.average(error(a2s, stds)) < avg_tol:
        return a2s
    return a1s
