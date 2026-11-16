#!/usr/bin/env python

"""
Automated quality control (QC) tests for OOI data products, ported from
the MATLAB example code in the OOI Data Product Specifications. The stuck
value, spike, and gradient tests, and the month conversion used by the
local range test wrapper, call compiled routines in qc_extensions.

Author: Christopher Mueller
"""
from ion_functions.qc.qc_extensions import stuckvalues, spikevalues, gradientvalues, ntp_to_month

import time
import numpy as np
import numexpr as ne
from scipy.interpolate import LinearNDInterpolator
from ion_functions import deprecated, utils
from ion_functions.utils import fill_value

# try to load the OOI logging module, using default Python logging module if
# unavailable
try:
    from ooi.logging import log
except ImportError:
    import logging
    log = logging.getLogger('ion-functions')

def is_fill(arr):
    """
    Tests whether a QC parameter holds the -9999 fill value.

    Parameters
    ----------
    arr : array_like
        QC parameter value(s) from the QC lookup table. Only the last
        element is tested.

    Returns
    -------
    out : bool
        True if the last element equals -9999.

    Notes
    -----
    -9999 is the fill value used in the QC parameter tables; it is not the
    system fill value used for integer data products.
    """
    return np.atleast_1d(arr)[-1] == -9999. # Not the normal fill value, it's hardcoded to the QC params
def is_none(arr):
    """
    Tests whether a QC parameter is missing.

    Parameters
    ----------
    arr : array_like or None
        QC parameter value(s) from the QC lookup table. Only the last
        element is tested.

    Returns
    -------
    out : bool
        True if arr is None or its last element is None.
    """
    return arr is None or (np.atleast_1d(arr).flatten()[-1] == None)

def dataqc_globalrangetest_minmax(dat, dat_min, dat_max, strict_validation=False):
    """
    OOI wrapper for the Global Range Test (GLBLRNG). Combines separate
    minimum and maximum limits into the two-element range expected by
    dataqc_globalrangetest.

    Parameters
    ----------
    dat : ndarray
        Input data.
    dat_min : array_like
        Minimum valid value; the last element is used.
    dat_max : array_like
        Maximum valid value; the last element is used.
    strict_validation : bool, optional
        Passed to dataqc_globalrangetest. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        1 where dat is within the range, 0 elsewhere. -99 for every point
        if either limit is missing or holds the -9999 fill value.

    See Also
    --------
    dataqc_globalrangetest : Core algorithm.
    """
    if is_none(dat_min) or is_none(dat_max) or is_fill(dat_min) or is_fill(dat_max):
        out = np.empty(dat.shape, dtype=np.int8)
        out.fill(-99)
        return out
    return dataqc_globalrangetest(dat, [np.atleast_1d(dat_min)[-1], np.atleast_1d(dat_max)[-1]], strict_validation=strict_validation)


def dataqc_globalrangetest(dat, datlim, strict_validation=False):
    """
    Global Range Test (GLBLRNG): flags data that fall outside a fixed
    valid range.

    Parameters
    ----------
    dat : array_like
        Input data; numeric and real.
    datlim : array_like
        Valid range. The minimum and maximum of datlim are used as the
        lower and upper limits.
    strict_validation : bool, optional
        If True, raise ValueError unless dat and datlim are numeric and
        real and datlim has at least 2 elements. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        1 where min(datlim) <= dat <= max(datlim), 0 elsewhere.
    """
    dat = np.atleast_1d(dat)
    datlim = np.atleast_1d(datlim)

    if strict_validation:
        if not utils.isnumeric(dat).all():
            raise ValueError('\'dat\' must be numeric')

        if not utils.isreal(dat).all():
            raise ValueError('\'dat\' must be real')

        if not utils.isnumeric(datlim).all():
            raise ValueError('\'datlim\' must be numeric')

        if not utils.isreal(datlim).all():
            raise ValueError('\'datlim\' must be real')

        if len(datlim) < 2:  # Must have at least 2 elements
            raise ValueError('\'datlim\' must have at least 2 elements')

    return (datlim.min() <= dat) & (dat <= datlim.max()).astype('int8')

def dataqc_localrangetest_wrapper(dat, datlim, datlimz, dims, pval_callback):
    """
    OOI wrapper for the Local Range Test (LOCLRNG). Builds the coordinate
    array z from the named dimensions and calls dataqc_localrangetest.

    Parameters
    ----------
    dat : ndarray
        Input data.
    datlim : ndarray
        Valid ranges, one [min, max] row per row of datlimz.
    datlimz : ndarray
        Coordinates at which datlim is given, one column per dimension.
    dims : list of str
        Names of the coordinate dimensions, in the column order of datlimz.
        'month' is computed from the 'time' parameter (NTP seconds); any
        other name is fetched with pval_callback.
    pval_callback : callable
        Function returning the values of a named parameter.

    Returns
    -------
    qcflag : ndarray of int8
        Output of dataqc_localrangetest. -99 for every point if an input
        is missing or datlim holds the -9999 fill value.

    See Also
    --------
    dataqc_localrangetest : Core algorithm.
    """
    if is_none(datlim) or np.all(np.atleast_1d(datlim).flatten() == -9999):
        out = np.empty(dat.shape, dtype=np.int8)
        out.fill(-99)
        return out
    if is_none(datlimz) or np.all(np.atleast_1d(datlim).flatten() == -9999):
        out = np.empty(dat.shape, dtype=np.int8)
        out.fill(-99)
        return out
    if is_none(dims):
        out = np.empty(dat.shape, dtype=np.int8)
        out.fill(-99)
        return out
    if is_none(pval_callback):
        out = np.empty(dat.shape, dtype=np.int8)
        out.fill(-99)
        return out

    z = []
    for dim in dims:
        if dim == 'month':
            # Convert time vector to vector of months
            v = pval_callback('time')
            v = np.asanyarray(v, dtype=float)
            v = ntp_to_month(v)
            z.append(v)
        else:
            # Fetch the dimension from the callback method
            v = pval_callback(dim)
            z.append(v)
    if len(dims)>1:
        z = np.column_stack(z)
    else:
        z = z[0]
        datlimz = datlimz[:,0]
    return dataqc_localrangetest(dat, z, datlim, datlimz)


def dataqc_localrangetest(dat, z, datlim, datlimz, strict_validation=False):
    """
    Local Range Test (LOCLRNG): flags data that fall outside a valid range
    that varies with one or more coordinates (e.g. depth or month).

    Parameters
    ----------
    dat : ndarray
        Input data; a numeric, real vector.
    z : ndarray
        Coordinates of each value in dat: one row per value and the same
        number of columns as datlimz.
    datlim : ndarray
        Valid ranges: a 2-column array of [min, max], one row per row of
        datlimz.
    datlimz : ndarray
        Coordinates at which datlim is given.
    strict_validation : bool, optional
        If True, check that dat is a vector, datlim is a matrix, and all
        inputs are numeric and real. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        1 where dat is within the interpolated limits, 0 elsewhere.

    Raises
    ------
    ValueError
        If the shapes of dat, z, datlim, and datlimz are inconsistent, or
        a minimum in datlim is not less than its maximum.

    Notes
    -----
    The limits are linearly interpolated to z. Points outside the
    coordinates covered by datlimz are flagged 0.
    """
    if strict_validation:
        # check if dat and datlim are matrices
        if not utils.isvector(dat):
            raise ValueError('\'dat\' must be a matrix')

        if not utils.ismatrix(datlim):
            raise ValueError('\'datlim\' must be a matrix')

        # check if all inputs are numeric and real
        for k, arg in {'dat': dat, 'z': z, 'datlim': datlim,
                       'datlimz': datlimz}.items(): # Might not be a proper 2to3 port, keep an eye on it. https://stackoverflow.com/questions/13998492/when-should-iteritems-be-used-instead-of-items
            if not utils.isnumeric(arg).all():
                raise ValueError('\'{0}\' must be numeric'.format(k))

            if not utils.isreal(arg).all():
                raise ValueError('\'{0}\' must be real'.format(k))

    if len(datlim.shape) == 3 and datlim.shape[0] == 1:
        datlim = datlim.reshape(datlim.shape[1:])

    if len(datlimz.shape) == 3 and datlimz.shape[0] == 1:
        datlimz = datlimz.reshape(datlimz.shape[1:])

    # test size and shape of the input arrays datlimz and datlim, setting test
    # variables.
    array_size = datlimz.shape
    if len(array_size) == 1:
        numlim = array_size[0]
        ndim = 1
    else:
        numlim = array_size[0]
        ndim = array_size[1]

    array_size = datlim.shape
    tmp1 = array_size[0]
    tmp2 = array_size[1]
    if tmp1 != numlim:
        raise ValueError('\'datlim\' and \'datlimz\' must '
                         'have the same number of rows.')

    if tmp2 != 2:
        raise ValueError('\'datlim\' must be structured as 2-D array '
                         'with exactly 2 columns and 1 through N rows.')

    # test the size and shape of the z input array
    array_size = z.shape
    if len(array_size) == 1:
        num = array_size[0]
        tmp2 = 1
    else:
        num = array_size[0]
        tmp2 = array_size[1]

    if tmp2 != ndim:
        raise ValueError('\'z\' must have the same number of columns '
                         'as \'datlimz\'.')

    if num != dat.size:
        raise ValueError('Len of \'dat\' must match number of '
                         'rows in \'z\'')

    # test datlim, values in column 2 must be greater than those in column 1
    if not all(datlim[:, 1] > datlim[:, 0]):
        raise ValueError('Second column values of \'datlim\' should be '
                         'greater than first column values.')

    # calculate the upper and lower limits for the data set
    if ndim == 1:
        # determine the lower limits using linear interpolation
        lim1 = np.interp(z, datlimz, datlim[:, 0], left=np.nan, right=np.nan)
        # determine the upper limits using linear interpolation
        lim2 = np.interp(z, datlimz, datlim[:, 1], left=np.nan, right=np.nan)
    else:
        # Compute Delaunay Triangulation and use linear interpolation to
        # determine the N-dimensional lower limits
        F = LinearNDInterpolator(datlimz, datlim[:, 0].reshape(numlim, 1))
        lim1 = F(z).reshape(dat.size)

        # Compute Delaunay Triangulation and use linear interpolation to
        # determine the N-dimensional upper limits
        F = LinearNDInterpolator(datlimz, datlim[:, 1].reshape(numlim, 1))
        lim2 = F(z).reshape(dat.size)

    # replace NaNs from above interpolations
    ff = (np.isnan(lim1)) | (np.isnan(lim2))
    lim1[ff] = np.max(datlim[:, 1])
    lim2[ff] = np.min(datlim[:, 0])

    # compute the qcflags
    qcflag = (dat >= lim1) & (dat <= lim2)
    return qcflag.astype('int8')


def dataqc_spiketest_wrapper(dat, acc, N, L, strict_validation=False):
    """
    OOI wrapper for the Spike Test (SPKETST). Extracts scalar test
    parameters and calls dataqc_spiketest.

    Parameters
    ----------
    dat : ndarray
        Input data.
    acc : array_like
        Measurement accuracy; the last element is used.
    N : array_like
        Range multiplier; the last element is used.
    L : array_like
        Window length; the last element is used.
    strict_validation : bool, optional
        Passed to dataqc_spiketest. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        Output of dataqc_spiketest. -99 for every point if a parameter is
        missing or holds the -9999 fill value.

    See Also
    --------
    dataqc_spiketest : Core algorithm.
    """
    if is_none(acc) or is_fill(acc) or is_none(N) or is_fill(N) or is_none(L) or is_fill(L):
        out = np.empty(dat.shape, dtype=np.int8)
        out.fill(-99)
        return out
    return dataqc_spiketest(dat, np.atleast_1d(acc)[-1], np.atleast_1d(N)[-1], np.atleast_1d(L)[-1], strict_validation=strict_validation)

def dataqc_spiketest(dat, acc, N=5, L=5, strict_validation=False):
    """
    Spike Test (SPKETST): flags single values that deviate strongly from
    their neighbors.

    Parameters
    ----------
    dat : array_like
        Input data; a numeric, real vector.
    acc : float
        Measurement accuracy. Deviations smaller than N * acc are never
        flagged; use 0 to disable this threshold.
    N : float, optional
        Range multiplier. Default 5.
    L : int, optional
        Window length; an odd integer. Default 5.
    strict_validation : bool, optional
        If True, check that dat is a numeric, real vector and that acc, N,
        and L are numeric and real. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        1 for good data, 0 for spikes.

    Notes
    -----
    Each value is compared with the mean of its L-1 neighbors and flagged
    if it deviates by N * max(R, acc) or more, where R is the range of the
    neighbors. The test runs in compiled code (qc_extensions).
    """
    dat = np.atleast_1d(dat)

    if strict_validation:
        if not utils.isnumeric(dat).all():
            raise ValueError('\'dat\' must be numeric')

        if not utils.isreal(dat).all():
            raise ValueError('\'dat\' must be real')

        if not utils.isvector(dat):
            raise ValueError('\'dat\' must be a vector')

        for k, arg in {'acc': acc, 'N': N, 'L': L}.items(): # Might not be a proper 2to3 port, keep an eye on it. https://stackoverflow.com/questions/13998492/when-should-iteritems-be-used-instead-of-items
            if not utils.isnumeric(arg).all():
                raise ValueError('\'{0}\' must be numeric'.format(k))

            if not utils.isreal(arg).all():
                raise ValueError('\'{0}\' must be real'.format(k))
    dat = np.asanyarray(dat, dtype=float)
    
    out = spikevalues(dat, L, N, acc)
    return out


def dataqc_polytrendtest_wrapper(dat, t, ord_n, nstd, strict_validation=False):
    """
    OOI wrapper for the Trend Test (TRNDTST). Extracts scalar test
    parameters and calls dataqc_polytrendtest.

    Parameters
    ----------
    dat : ndarray
        Input data.
    t : ndarray
        Time associated with dat.
    ord_n : array_like
        Polynomial order; the last element is used.
    nstd : array_like
        Threshold factor; the last element is used.
    strict_validation : bool, optional
        Passed to dataqc_polytrendtest. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        Output of dataqc_polytrendtest. -99 for every point if a parameter
        is missing or ord_n holds the -9999 fill value.

    See Also
    --------
    dataqc_polytrendtest : Core algorithm.
    """
    if is_none(ord_n) or is_fill(ord_n) or is_none(nstd) or is_fill(ord_n):
        out = np.empty(dat.shape, dtype=np.int8)
        out.fill(-99)
        return out
    return dataqc_polytrendtest(dat, t, np.atleast_1d(ord_n)[-1], np.atleast_1d(nstd)[-1], strict_validation=strict_validation)

def dataqc_polytrendtest(dat, t, ord_n=1, nstd=3, strict_validation=False):
    """
    Trend Test (TRNDTST): flags a time series that contains a significant
    polynomial trend, such as sensor drift.

    Parameters
    ----------
    dat : array_like
        Input data; a numeric, real vector.
    t : array_like
        Time associated with dat.
    ord_n : int, optional
        Polynomial order; 1 tests for a linear trend. Default 1.
    nstd : int, optional
        Threshold factor. Default 3.
    strict_validation : bool, optional
        If True, check that dat and t are numeric, real vectors and that
        ord_n and nstd are numeric, real scalars. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        The same value for every point: 0 if a trend is detected, 1
        otherwise.

    Notes
    -----
    A trend is detected if nstd times the standard deviation of the
    residuals from the polynomial fit is less than the standard deviation
    of dat. ord_n is rounded to the nearest integer and nstd is truncated
    to an integer.
    """
    dat = np.atleast_1d(dat)
    t = np.atleast_1d(t)

    if strict_validation:
        for k, arg in {'dat': dat, 't': t, 'ord_n': ord_n, 'nstd': nstd}.items(): # Might not be a proper 2to3 port, keep an eye on it. https://stackoverflow.com/questions/13998492/when-should-iteritems-be-used-instead-of-items
            if not utils.isnumeric(arg).all():
                raise ValueError('\'{0}\' must be numeric'.format(k))

            if not utils.isreal(arg).all():
                raise ValueError('\'{0}\' must be real'.format(k))

        for k, arg in {'dat': dat, 't': t}.items(): # Might not be a proper 2to3 port, keep an eye on it. https://stackoverflow.com/questions/13998492/when-should-iteritems-be-used-instead-of-items
            if not utils.isvector(arg):
                raise ValueError('\'{0}\' must be a vector'.format(k))

        for k, arg in {'ord_n': ord_n, 'nstd': nstd}.items(): # Might not be a proper 2to3 port, keep an eye on it. https://stackoverflow.com/questions/13998492/when-should-iteritems-be-used-instead-of-items
            if not utils.isscalar(arg):
                raise ValueError('\'{0}\' must be a scalar'.format(k))

    ord_n = int(round(abs(ord_n)))
    nstd = int(abs(nstd))

    ll = len(dat)
    # Not needed because time is incorporated as 't'
    # t = range(ll)

    pp = np.polyfit(t, dat, ord_n)
    datpp = np.polyval(pp, t)

    # test for a trend
    if np.atleast_1d((np.std(dat - datpp) * nstd) < np.std(dat)).all():
        trndtst = 0
    else:
        trndtst = 1

    # insure output size equals input, even though test yields a single value.
    qcflag = np.ones(dat.shape).astype('int8') * trndtst
    return qcflag


def dataqc_stuckvaluetest_wrapper(x, reso, num, strict_validation=False):
    """
    OOI wrapper for the Stuck Value Test (STUCKVL). Extracts scalar test
    parameters and calls dataqc_stuckvaluetest.

    Parameters
    ----------
    x : ndarray
        Input time series.
    reso : array_like
        Resolution; the last element is used.
    num : array_like
        Minimum number of repeats; the last element is used.
    strict_validation : bool, optional
        Passed to dataqc_stuckvaluetest. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        Output of dataqc_stuckvaluetest. -99 for every point if a parameter
        is missing or holds the -9999 fill value.

    See Also
    --------
    dataqc_stuckvaluetest : Core algorithm.
    """
    if is_none(reso) or is_fill(reso) or is_none(num) or is_fill(num):
        out = np.empty(x.shape, np.int8)
        out.fill(-99)
        return out
    return  dataqc_stuckvaluetest(x, np.atleast_1d(reso)[-1], np.atleast_1d(num)[-1], strict_validation=strict_validation)

def dataqc_stuckvaluetest(x, reso, num=10, strict_validation=False):
    """
    Stuck Value Test (STUCKVL): flags runs of repeated values in a time
    series.

    Parameters
    ----------
    x : array_like
        Input time series; a numeric, real vector.
    reso : float
        Resolution; successive values less than reso apart count as
        repeats.
    num : int, optional
        Minimum number of successive repeats that marks values as stuck.
        Default 10.
    strict_validation : bool, optional
        If True, check that x is a numeric, real vector and that reso and
        num are numeric, real scalars. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        0 where stuck values are found, 1 elsewhere. Every point is 0 if x
        has fewer than num values.

    Notes
    -----
    The test runs in compiled code (qc_extensions).
    """
    dat = np.atleast_1d(x)

    if strict_validation:
        if not utils.isnumeric(dat).all():
            raise ValueError('\'x\' must be numeric')

        if not utils.isvector(dat):
            raise ValueError('\'x\' must be a vector')

        if not utils.isreal(dat).all():
            raise ValueError('\'x\' must be real')

        for k, arg in {'reso': reso, 'num': num}.items(): # Might not be a proper 2to3 port, keep an eye on it. https://stackoverflow.com/questions/13998492/when-should-iteritems-be-used-instead-of-items
            if not utils.isnumeric(arg).all():
                raise ValueError('\'{0}\' must be numeric'.format(k))

            if not utils.isscalar(arg):
                raise ValueError('\'{0}\' must be a scalar'.format(k))

            if not utils.isreal(arg).all():
                raise ValueError('\'{0}\' must be real'.format(k))

    num = np.abs(num)
    dat = np.asanyarray(dat, dtype=float)
    ll = len(x)
    if ll < num:
        # Warn - 'num' is greater than len(x), returning zeros
        out = np.zeros(dat.size, dtype='int8')
    else:
        out = stuckvalues(dat, reso, num)

    return out

def dataqc_gradienttest_wrapper(dat, x, ddatdx, mindx, startdat, toldat, strict_validation=False):
    """
    OOI wrapper for the Gradient Test (GRADTST). Builds the symmetric valid
    range [-ddatdx, ddatdx] from a single gradient limit and calls
    dataqc_gradienttest.

    Parameters
    ----------
    dat : ndarray
        Input data.
    x : ndarray
        Coordinate of each value in dat.
    ddatdx : array_like
        Maximum absolute change in dat per unit x; the last element is
        used.
    mindx : array_like
        Minimum separation in x; the last element is used.
    startdat : array_like
        Start value presumed good; the last element is used.
    toldat : array_like
        Tolerance for returning to good data; the last element is used.
    strict_validation : bool, optional
        Passed to dataqc_gradienttest. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        Output of dataqc_gradienttest. -99 for every point if a parameter
        is missing or holds the -9999 fill value.

    See Also
    --------
    dataqc_gradienttest : Core algorithm.
    """
    if is_none(ddatdx) or is_fill(ddatdx) or is_none(mindx) or is_fill(mindx) or is_none(startdat) or is_fill(startdat) or is_none(toldat) or is_fill(toldat):
        out = np.empty(dat.shape, dtype=np.int8)
        out.fill(-99)
        return out
    outqc = dataqc_gradienttest(dat, x, [-np.atleast_1d(ddatdx)[-1], np.atleast_1d(ddatdx)[-1]], np.atleast_1d(mindx)[-1], np.atleast_1d(startdat)[-1], np.atleast_1d(toldat)[-1], strict_validation=strict_validation)
    return outqc


def dataqc_gradienttest(dat, x, ddatdx, mindx, startdat, toldat, strict_validation=False):
    """
    Gradient Test (GRADTST): flags data where the change between successive
    points falls outside a valid range.

    Parameters
    ----------
    dat : array_like
        Input data; a numeric, real vector.
    x : array_like
        Coordinate (e.g. time or distance) of each value in dat; the same
        length as dat and strictly increasing.
    ddatdx : array_like
        Two-element [min, max] valid range of the change in dat per unit x.
    mindx : float
        Minimum separation in x for the test to be applied; points closer
        than mindx to the last tested point are skipped. NaN or 0 disables
        this.
    startdat : float
        Value presumed good, used to test the first point. NaN or 0 means
        the first point is presumed good.
    toldat : float
        Tolerance. After a failure, data are flagged 0 until a value comes
        within toldat of the last good value.
    strict_validation : bool, optional
        If True, check that dat and x are vectors of equal length and that
        x is strictly increasing. Default False.

    Returns
    -------
    qcflag : ndarray of int8
        1 for good data, 0 for bad data, and -99 for points skipped because
        of mindx.

    Notes
    -----
    The test runs in compiled code (qc_extensions).
    """

    if strict_validation:
        if not utils.isvector(dat) or not utils.isvector(x):
            raise ValueError('\'dat\' and \'x\' must be vectors')

        if len(dat) != len(x):
            raise ValueError('\'dat\' and \'x\' must be of equal len')

        if not all(np.diff(x) > 0):
            raise ValueError('\'x\' must be montonically increasing')

    dat = np.asanyarray(dat, dtype=float).flatten()
    x = np.asanyarray(x, dtype=float).flatten()

    if np.isnan(mindx):
        mindx = 0
    mindx = mindx or 0
    if np.isnan(startdat):
        startdat = 0
    startdat = startdat or 0
    
    # No strict validation here, they are scalards and they must be validated
    # before going into the C-layer
    if not utils.isscalar(mindx):
        raise ValueError("'mindx' must be scalar, NaN, or empty.")
    if not utils.isscalar(startdat):
        raise ValueError("'startdat' must be scalar, NaN, or empty.")


    # Confirm that there are still data points left, else abort:
    if np.abs(x[0] - x[-1]) < mindx:
        out = np.zeros(x.shape)
        out.fill(1)
        log.warn('Too few values to inspect')
        return out


    grad_min = ddatdx[0]
    grad_max = ddatdx[1]
    out = gradientvalues(dat, x, grad_min, grad_max, mindx, startdat, toldat)

    return out


@deprecated
def dataqc_solarelevation(lon, lat, dt):
    """
    Solar Elevation (SOLAREL): computes the solar elevation and the solar
    radiation in the absence of an atmosphere for a given location and
    time.

    Parameters
    ----------
    lon : array_like
        Longitude [decimal degrees]; east is positive.
    lat : array_like
        Latitude [decimal degrees]; north is positive.
    dt : array_like
        Time stamp [seconds since 1970-01-01 UTC].

    Returns
    -------
    z : ndarray
        Solar elevation (altitude) [degrees].
    sorad : ndarray
        Solar radiation with no atmosphere [W m^-2]; 0 when the sun is
        below the horizon.

    Notes
    -----
    Python implementation of soradna1.m from the Air-Sea Toolbox, using
    reduced-accuracy expressions valid for the years 1800-2100.
    Deprecated: never used in the OOI data system.
    """
    # Test lengths and types of inputs. Latitude and longitude must be the same
    # size and can either be a scalar or a vecotr. The date and time stamp
    # can also be either a scalar or a vector. If all three inputs are vectors,
    # they must be of the same length.
    if len(lon) != len(lat):
        raise ValueError('\'lon\' and \'lat\' must be the same size')

    if utils.isvector(lon) and utils.isvector(lat) and utils.isvector(dt):
        # test their lengths
        if not len(lon) == len(lat) == len(dt):
            raise ValueError('If all inputs are vectors, these must all '
                             'be of the same length')

    # set constants (using values from as_consts.m)
    # ------ short-wave flux calculations
    #   the solar constant [W m^-2] represents a mean of satellite measurements
    #   made over the last sunspot cycle (1979-1995), taken from Coffey et al.
    #   (1995), Earth System Monitor, 6, 6-10.
    solar_const = 1368.0

    # Create a time tuple in UTC from the Epoch time input, and then create
    # scalars or numpy arrays of time elements for subsequent calculations.
    ldt = len(dt)
    yy = np.zeros(ldt, dtype=int)
    mn = np.zeros(ldt, dtype=int)
    dd = np.zeros(ldt, dtype=int)
    hh = np.zeros(ldt, dtype=int)
    mm = np.zeros(ldt, dtype=int)
    ss = np.zeros(ldt, dtype=int)
    for i in range(ldt):
        # create time tuple in UTC
        gtime = time.gmtime(dt[i])
        # create scalar elements
        yy[i] = gtime[0]
        mn[i] = gtime[1]
        dd[i] = gtime[2]
        hh[i] = gtime[3]
        mm[i] = gtime[4]
        ss[i] = gtime[5]

    #constants used in function
    deg2rad = np.pi / 180.0
    rad2deg = 1 / deg2rad

    # compute Universal Time in hours
    utime = hh + (mm + ss / 60.0) / 60.0

    # compute Julian ephemeris date in days (Day 1 is 1 Jan 4713 B.C. which
    # equals -4712 Jan 1)
    jed = (367.0 * yy - np.fix(7.0*(yy+np.fix((mn+9)/12.0))/4.0)
           + np.fix(275.0*mn/9.0) + dd + 1721013 + utime / 24.0)

    # compute interval in Julian centuries since 1900
    jc_int = (jed - 2415020.0) / 36525.0

    # compute mean anomaly of the sun
    ma_sun = 358.475833 + 35999.049750 * jc_int - 0.000150 * jc_int**2
    ma_sun = (ma_sun - np.fix(ma_sun/360.0) * 360.0) * deg2rad

    # compute mean longitude of sun
    ml_sun = 279.696678 + 36000.768920 * jc_int + 0.000303 * jc_int**2
    ml_sun = (ml_sun - np.fix(ml_sun/360.0) * 360.0) * deg2rad

    # compute mean anomaly of Jupiter
    ma_jup = 225.444651 + 2880.0 * jc_int + 154.906654 * jc_int
    ma_jup = (ma_jup - np.fix(ma_jup/360.0) * 360.0) * deg2rad

    # compute longitude of the ascending node of the moon's orbit
    an_moon = (259.183275 - 1800 * jc_int - 134.142008 * jc_int
               + 0.002078 * jc_int**2)
    an_moon = (an_moon - np.fix(an_moon/360.0) * 360.0 + 360.0) * deg2rad

    # compute mean anomaly of Venus
    ma_ven = (212.603219 + 58320 * jc_int + 197.803875 * jc_int
              + 0.001286 * jc_int**2)
    ma_ven = (ma_ven - np.fix(ma_ven/360.0) * 360.0) * deg2rad

    # compute sun theta
    theta = (0.397930 * np.sin(ml_sun) + 0.009999 * np.sin(ma_sun-ml_sun)
             + 0.003334 * np.sin(ma_sun+ml_sun) - 0.000208 * jc_int
             * np.sin(ml_sun) + 0.000042 * np.sin(2*ma_sun+ml_sun) - 0.000040
             * np.cos(ml_sun) - 0.000039 * np.sin(an_moon-ml_sun) - 0.000030
             * jc_int * np.sin(ma_sun-ml_sun) - 0.000014
             * np.sin(2*ma_sun-ml_sun) - 0.000010
             * np.cos(ma_sun-ml_sun-ma_jup) - 0.000010 * jc_int
             * np.sin(ma_sun+ml_sun))

    # compute sun rho
    rho = (1.000421 - 0.033503 * np.cos(ma_sun) - 0.000140 * np.cos(2*ma_sun)
           + 0.000084 * jc_int * np.cos(ma_sun) - 0.000033
           * np.sin(ma_sun-ma_jup) + 0.000027 * np.sin(2.*ma_sun-2.*ma_ven))

    # compute declination
    decln = np.arcsin(theta/np.sqrt(rho))

    # compute equation of time (in seconds of time)
    l = 276.697 + 0.98564734 * (jed-2415020.0)
    l = (l - 360.0 * np.fix(l/360.0)) * deg2rad
    eqt = (-97.8 * np.sin(l) - 431.3 * np.cos(l) + 596.6 * np.sin(2*l)
           - 1.9 * np.cos(2*l) + 4.0 * np.sin(3*l) + 19.3 * np.cos(3*l)
           - 12.7 * np.sin(4*l))
    eqt = eqt / 60.0

    # compute local hour angle from global hour angle
    gha = 15.0 * (utime-12) + 15.0 * eqt / 60.0
    lha = gha - lon

    # compute radius vector
    rv = np.sqrt(rho)

    # compute solar altitude
    sz = (np.sin(deg2rad*lat) * np.sin(decln) + np.cos(deg2rad*lat)
          * np.cos(decln) * np.cos(deg2rad*lha))
    z = rad2deg * np.arcsin(sz)

    # compute solar radiation outside atmosphere (defaults to 0 when solar
    # altitude is below the horizon)
    sorad = (solar_const / rv**2) * np.sin(deg2rad * z)
    sorad[z < 0] = 0

    return (z, sorad)

def dataqc_propagateflags_wrapper(strict_validation=False, *args):
    """
    OOI wrapper for Combined QC Flags (CMBNFLG). Stacks a variable number of
    flag vectors and calls dataqc_propagateflags.

    Parameters
    ----------
    strict_validation : bool, optional
        Passed to dataqc_propagateflags. Default False.
    *args : ndarray
        Flag vectors of the same length.

    Returns
    -------
    outflag : ndarray of int8
        1 where all inputs are 1, 0 elsewhere.

    Raises
    ------
    ValueError
        If the vectors differ in length (checked only when
        strict_validation is False).

    See Also
    --------
    dataqc_propagateflags : Core algorithm.
    """

    if not strict_validation:

        shapes = np.array([i.shape[0] for i in args])
        if not (shapes == shapes[0]).all():
            raise ValueError('Input vectors are not the same shape')

    return dataqc_propagateflags(np.array(args), strict_validation=strict_validation)



def dataqc_propagateflags(inflags, strict_validation=False):
    """
    Combined QC Flags (CMBNFLG): merges several sets of QC flags into one,
    e.g. to flag a derived product such as salinity from the flags of its
    inputs (conductivity, temperature, and pressure).

    Parameters
    ----------
    inflags : ndarray
        M-by-N array of flags, one row per input data set; 1 is good and 0
        is bad.
    strict_validation : bool, optional
        If True, raise ValueError unless inflags contains only 0 and 1.
        Default False.

    Returns
    -------
    outflag : ndarray of int8
        Length-N flags: 1 where every input is 1, 0 otherwise.
    """

    if strict_validation:
        if not utils.islogical(inflags):
            raise ValueError('\'inflags\' must be \'0\' or \'1\' '
                             'integer flag array')

    array_size = inflags.shape
    nrows = array_size[0]
    if nrows < 2:
        error('\'inflags\' must be at least a two-dimensional array')

    outflag = np.all(inflags, 0)
    return outflag.astype('int8')


@deprecated
def dataqc_condcompress(p_orig, p_new, c_orig, cpcor=-9.57e-8):
    """
    Conductivity Compressibility Compensation (CONDCMP): rescales
    conductivity computed with one pressure to a corrected pressure.

    Parameters
    ----------
    p_orig : array_like
        Pressure used to compute the original conductivity, typically L1a
        PRESWAT [dbar].
    p_new : array_like
        Updated pressure, typically L1b PRESWAT [dbar].
    c_orig : array_like
        Original conductivity, typically L1a CONDWAT [S m^-1].
    cpcor : float, optional
        Pressure correction coefficient used to compute c_orig. Default
        -9.57e-8.

    Returns
    -------
    c_new : array_like
        Conductivity corrected to p_new [S m^-1].

    Notes
    -----
    Deprecated: never used in the OOI data system.
    """
    c_new = c_orig * (1 + cpcor * p_orig) / (1 + cpcor * p_new)
    return c_new
