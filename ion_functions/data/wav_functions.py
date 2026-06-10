#!/usr/bin/env python
"""
Module containing WAV (WAVSS) wave statistics data-calculation functions
for the Ocean Observatories Initiative.
"""
import numpy as np

from ion_functions.data.generic_functions import magnetic_declination
from ion_functions.utils import fill_value


def wav_triaxys_nondir_freq(nfreq, freq0, delta_freq):
    """
    Compute frequency values for non-directional wave spectral bins
    (WAVSTAT-FND_L1) for the WAVSS instrument class.

    Parameters
    ----------
    nfreq : array_like
        Number of non-directional frequency bins, from the value specified
        in the WAVSS $TSPNA data sentence. All data packets must have the
        same value.
    freq0 : array_like
        Initial frequency [Hz], from the value specified in the WAVSS
        $TSPNA data sentence.
    delta_freq : array_like
        Frequency spacing [Hz], from the value specified in the WAVSS
        $TSPNA data sentence.

    Returns
    -------
    fnd : ndarray
        Frequency values for non-directional wave spectral bins
        (WAVSTAT-FND_L1) [Hz]. Shape is (npackets, nfreq).

    Notes
    -----
    The frequency vector is reconstructed from the initial frequency and
    uniform spacing reported by the instrument. All data packets share the
    same number of non-directional frequency bins, so the reconstruction
    is fully vectorized across packets.
    """
    # condition input variables
    nfreq = np.array(nfreq, ndmin=1)
    freq0 = np.array(freq0, ndmin=1)
    delta_freq = np.array(delta_freq, ndmin=1)

    # each set of data inputs will call for the same number of frequency values nfreq.
    # therefore can vectorize (without using forloop as had to be done for the
    # directional frequencies case) by setting up all variables as 2D arrays
    # of size(npackets, nfreq).
    npackets = nfreq.shape[0]
    n_freqs = nfreq[0]

    # orient arrays such that lead index indexes each set of data inputs
    freq0_2d = np.tile(freq0, (n_freqs, 1)).transpose()
    delta_freq_2d = np.tile(delta_freq, (n_freqs, 1)).transpose()
    steps_2d = np.tile(np.arange(n_freqs), (npackets, 1))
    fnd = freq0_2d + steps_2d * delta_freq_2d

    return fnd


def wav_triaxys_dir_freq(nfreq_nondir, nfreq_dir, freq0, delta_freq):
    """
    Compute frequency values for directional wave spectral bins
    (WAVSTAT-FDS_L1) for the WAVSS instrument class.

    Parameters
    ----------
    nfreq_nondir : array_like
        Number of non-directional wave frequency bins, from the value
        specified in the WAVSS $TSPNA data sentence. All data packets must
        have the same value; this sets the output array width.
    nfreq_dir : array_like
        Number of directional wave frequency bins, from the value specified
        in the WAVSS $TSPMA data sentence. May vary between data packets.
    freq0 : array_like
        Initial frequency [Hz], from the value specified in the WAVSS
        $TSPMA data sentence.
    delta_freq : array_like
        Frequency spacing [Hz], from the value specified in the WAVSS
        $TSPMA data sentence.

    Returns
    -------
    fds : ndarray
        Frequency values for directional wave spectral bins
        (WAVSTAT-FDS_L1) [Hz]. Shape is (npackets, nfreq_nondir), with
        fill values in positions beyond nfreq_dir for each packet.

    Notes
    -----
    The number of active directional frequency bins (nfreq_dir) can vary
    between data packets as a function of measured ocean conditions at
    fixed instrument settings. The output array is sized to the fixed
    non-directional bin count, with fill values in unused positions. This
    variable-length behavior also affects the WAVSTAT-PDS and WAVSTAT-SDS
    arrays and the WAVSTAT-DDS_L2 product.
    """
    # condition input variables.
    # all delta_freq and freq0 values will be floats.
    nfreq_nondir = np.array(nfreq_nondir, ndmin=1)
    nfreq_dir = np.array(nfreq_dir, ndmin=1)
    freq0 = np.array(freq0, ndmin=1)
    delta_freq = np.array(delta_freq, ndmin=1)

    # each data packet may call for a different number of directional frequency values nfreq_dir.
    # however, this number will always be <= nfreq_nondir, and all the nfreq_nondir values will be identical.
    npackets = nfreq_nondir.shape[0]
    fds = np.zeros((npackets, int(nfreq_nondir[0]))) + fill_value

    for ii in range(npackets):
        fds[ii, 0:nfreq_dir[ii]] = freq0[ii] + np.arange(nfreq_dir[ii]) * delta_freq[ii]

    return fds


def wav_triaxys_buoymotion_time(ntp_timestamp, ntime, time0, delta_time):
    """
    Compute NTP timestamps for buoy displacement measurements
    (WAVSTAT-MOTT_L1) for the WAVSS instrument class.

    Parameters
    ----------
    ntp_timestamp : array_like
        NTP timestamp of the $TSPHA data sentence [s since 1900-01-01].
    ntime : array_like
        Number of displacement time values, from the value specified in
        the WAVSS $TSPHA data sentence.
    time0 : array_like
        Elapsed time between ntp_timestamp and the first displacement
        sample, from the value specified in the WAVSS $TSPHA data sentence
        ("Initial Time") [s].
    delta_time : array_like
        Time interval between successive displacement measurements, from
        the value specified in the WAVSS $TSPHA data sentence
        ("Time Spacing") [s].

    Returns
    -------
    mott : ndarray
        NTP timestamps corresponding to buoy displacement measurements
        (WAVSTAT-MOTT_L1) [s since 1900-01-01]. Shape is
        (npackets, ntime).

    Notes
    -----
    The time vector is reconstructed from the sentence timestamp, initial
    time offset, and uniform sampling interval reported by the instrument.
    The resulting timestamps correspond element-wise to the
    WAVSTAT-MOTX, WAVSTAT-MOTY, and WAVSTAT-MOTZ displacement arrays.
    """
    # condition input variables;
    # make sure time interval is not type integer
    ntime = np.array(ntime, ndmin=1)
    time0 = np.array(time0, ndmin=1)
    delta_time = np.array(delta_time, dtype='float', ndmin=1)

    # this algorithm is almost identical to that contained in def wav_wavss_nondir_freq above.
    # these are the dimensions of all the 2D arrays used in the calculation
    npackets = ntime.shape[0]
    n_time_values = ntime[0]

    # orient the lead index to iterate over the data packet number
    ntp0_2d = np.tile(ntp_timestamp + time0, (n_time_values, 1)).transpose()
    delta_time_2d = np.tile(delta_time, (n_time_values, 1)).transpose()
    steps_2d = np.tile(np.arange(n_time_values), (npackets, 1))

    mott = ntp0_2d + steps_2d * delta_time_2d

    return mott


def wav_triaxys_correct_mean_wave_direction(dir_raw, lat, lon, ntp_ts):
    """
    Compute mean wave direction corrected for magnetic declination
    (WAVSTAT-D_L2) for the WAVSS instrument class.

    Parameters
    ----------
    dir_raw : array_like
        Uncorrected mean wave direction (WAVSTAT-D_L0) [deg, [0, 360)].
    lat : array_like
        Deployment latitude of the instrument [decimal deg]. North is
        positive, South is negative.
    lon : array_like
        Deployment longitude of the instrument [decimal deg]. East is
        positive, West is negative.
    ntp_ts : array_like
        NTP timestamp from the data particle [s since 1900-01-01].

    Returns
    -------
    dir_cor : ndarray
        Mean wave direction corrected for magnetic declination
        (WAVSTAT-D_L2) [deg, [0, 360)].

    Notes
    -----
    Magnetic declination is computed using the WMM2010 model via
    generic_functions.magnetic_declination. The WAVSS is a surface sensor;
    depth defaults to 0 (sea level) in the declination calculation.
    """
    # calculate the magnetic declination using the WWM2010 model
    # the WAVSS is a surface wave sensor, so that height above sealevel = 0,
    # which is the default value used in the magnetic_declination calculation.
    theta = magnetic_declination(lat, lon, ntp_ts)

    # directions are [0,360) degrees; and magnetic declinations can be positive or negative
    dir_cor = np.mod(dir_raw + theta + 360, 360)

    return dir_cor


def wav_triaxys_correct_directional_wave_direction(dir_raw, lat, lon, ntp_ts):
    """
    Compute directional wave directions corrected for magnetic declination
    (WAVSTAT-DDS_L2) for the WAVSS instrument class.

    Parameters
    ----------
    dir_raw : array_like
        Uncorrected directional wave directions (WAVSTAT-DDS_L0)
        [deg, [0, 360)]. Expected as a 2-D array of shape
        (npackets, nfreq_nondir) with fill values in unused positions.
    lat : array_like
        Deployment latitude of the instrument [decimal deg]. North is
        positive, South is negative.
    lon : array_like
        Deployment longitude of the instrument [decimal deg]. East is
        positive, West is negative.
    ntp_ts : array_like
        NTP timestamp from the data particle [s since 1900-01-01].

    Returns
    -------
    dir_cor : ndarray
        Directional wave directions corrected for magnetic declination
        (WAVSTAT-DDS_L2) [deg, [0, 360)]. Same shape as dir_raw; fill
        values in unused positions are preserved.

    Notes
    -----
    Magnetic declination is computed using the WMM2010 model via
    generic_functions.magnetic_declination. The WAVSS is a surface sensor;
    depth defaults to 0 (sea level) in the declination calculation.

    The number of active directional frequency bins can vary between data
    packets. Fill values mark unused positions; they are converted to NaN
    before the correction and restored afterward so that array operations
    do not alter those positions.
    """
    # assume that the dir_raw data product comes in as a 2D numpy array with fill values
    # appropriately placed to account for the cases in which the number of reported
    # directional wave frequency bins differs from data packet to data packet (and is
    # less than the number of reported non-directional frequency bins).
    dir_raw = np.array(dir_raw, ndmin=2)

    # change fill values to Nans, so that subsequent array operations will leave the
    # Nan entries unchanged.
    dir_raw[dir_raw == fill_value] = np.nan

    # calculate the magnetic declination using the WWM2010 model
    # the WAVSS is a surface wave sensor, so that height above sealevel = 0,
    # which is the default value used in the magnetic_declination calculation.
    theta = magnetic_declination(lat, lon, ntp_ts)

    # theta in general will be a vector, so replicate it into a matrix to match the dir_raw dimensions.
    theta = np.tile(theta, (dir_raw.shape[1], 1)).transpose()

    # directions are [0,360) degrees; and magnetic declinations can be positive or negative
    dir_cor = np.mod(dir_raw + theta + 360, 360)

    # replace Nans with fills
    dir_cor[np.isnan(dir_cor)] = fill_value

    return dir_cor


def wav_triaxys_magcor_buoymotion_x(x, y, lat, lon, ntp_timestamp):
    """
    OOI single-output wrapper for WAVSTAT-MOTX_L1. Returns eastward buoy
    displacement corrected for magnetic declination [m].

    Parameters
    ----------
    x : array_like
        Uncorrected eastward buoy displacement (WAVSTAT-MOTX_L0) [m].
    y : array_like
        Uncorrected northward buoy displacement (WAVSTAT-MOTY_L0) [m].
    lat : array_like
        Deployment latitude of the instrument [decimal deg]. North is
        positive, South is negative.
    lon : array_like
        Deployment longitude of the instrument [decimal deg]. East is
        positive, West is negative.
    ntp_timestamp : array_like
        NTP timestamp of the $TSPHA data sentence [s since 1900-01-01].
        One timestamp per displacement ensemble is sufficient; the maximum
        sampling period (35 min) is short enough that magnetic declination
        does not change meaningfully within an ensemble.

    Returns
    -------
    motx : ndarray
        Eastward buoy displacement corrected for magnetic declination
        (WAVSTAT-MOTX_L1) [m].

    See Also
    --------
    magnetic_correction_einsum : Core rotation algorithm; returns both
        corrected components.
    """
    # force shapes of inputs to arrays
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    ntp_timestamp = np.atleast_1d(ntp_timestamp)

    # calculate the magnetic declination using the WWM2010 model.
    # the WAVSS surface wave sensor is at sealevel, which is the default z value for mag dec.
    theta = magnetic_declination(lat, lon, ntp_timestamp)

    # correct for declination by rotating coordinates.
    # the function magnetic_correction_einsum was written to correct (u,v) velocities, but
    # it also applies to (E,N) coordinates.
    motx, _ = magnetic_correction_einsum(theta, x, y)

    return motx


def wav_triaxys_magcor_buoymotion_y(x, y, lat, lon, ntp_timestamp):
    """
    OOI single-output wrapper for WAVSTAT-MOTY_L1. Returns northward buoy
    displacement corrected for magnetic declination [m].

    Parameters
    ----------
    x : array_like
        Uncorrected eastward buoy displacement (WAVSTAT-MOTX_L0) [m].
    y : array_like
        Uncorrected northward buoy displacement (WAVSTAT-MOTY_L0) [m].
    lat : array_like
        Deployment latitude of the instrument [decimal deg]. North is
        positive, South is negative.
    lon : array_like
        Deployment longitude of the instrument [decimal deg]. East is
        positive, West is negative.
    ntp_timestamp : array_like
        NTP timestamp of the $TSPHA data sentence [s since 1900-01-01].
        One timestamp per displacement ensemble is sufficient; the maximum
        sampling period (35 min) is short enough that magnetic declination
        does not change meaningfully within an ensemble.

    Returns
    -------
    moty : ndarray
        Northward buoy displacement corrected for magnetic declination
        (WAVSTAT-MOTY_L1) [m].

    See Also
    --------
    magnetic_correction_einsum : Core rotation algorithm; returns both
        corrected components.
    """
    # force shapes of inputs to arrays
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    lat = np.atleast_1d(lat)
    lon = np.atleast_1d(lon)
    ntp_timestamp = np.atleast_1d(ntp_timestamp)

    # calculate the magnetic declination using the WWM2010 model.
    # the WAVSS surface wave sensor is at sealevel, which is the default z value for mag dec.
    theta = magnetic_declination(lat, lon, ntp_timestamp)

    # correct for declination by rotating coordinates.
    # the function magnetic_correction_einsum was written to correct (u,v) velocities, but
    # it also applies to (E,N) coordinates.
    _, moty = magnetic_correction_einsum(theta, x, y)

    return moty


def magnetic_correction_einsum(theta, u, v):
    """
    Correct velocity or displacement profiles for magnetic declination
    using a vectorized 2-D rotation.

    Parameters
    ----------
    theta : array_like
        Magnetic declination at the measurement location [deg]. One value
        per data packet (ensemble).
    u : array_like
        Uncorrected eastward velocity or displacement profiles. Shape must
        be (npackets, nsamples).
    v : array_like
        Uncorrected northward velocity or displacement profiles. Same
        shape as u.

    Returns
    -------
    u_cor : ndarray
        Eastward component corrected for magnetic declination.
    v_cor : ndarray
        Northward component corrected for magnetic declination.

    Notes
    -----
    This function handles the vectorized case of one magnetic declination
    value per ensemble of (u, v) pairs, i.e. theta=f(i) with u=f(i,j)
    and v=f(i,j). It uses numpy.einsum to apply the rotation matrix to
    all ensembles without a Python for loop.

    The executable code is identical to magnetic_correction_vctrzd in
    adcp_functions.py. It is reproduced here to avoid a cross-module
    dependency within wav_functions.py. The ADCP DPS citations below
    reflect the origin of the rotation algorithm.
    """
    # force shapes of inputs to arrays
    theta = np.atleast_1d(theta)
    u = np.atleast_2d(u)
    v = np.atleast_2d(v)

    theta_rad = np.radians(theta)
    cosT = np.cos(theta_rad)
    sinT = np.sin(theta_rad)

    # set up rotation matrix
    M = np.array([[cosT, sinT],
                  [-sinT, cosT]])

    # roll axes so that the lead index represents data packet #.
    M = np.rollaxis(M, 2)

    # construct the uncorrected velocity matrix.
    # the coordinate system is 2D, so the middle dimension is sized at 2.
    uv = np.zeros((u.shape[0], 2, u.shape[1]))

    # load the coordinates to be rotated into the appropriate slices
    uv[:, 0, :] = u
    uv[:, 1, :] = v

    # the Einstein summation is here configured to do the matrix
    # multiplication uv_cor(i,k) = M(i,j) * uv(j,k) on each slice h.
    uv_cor = np.einsum('hij,hjk->hik', M, uv)

    u_cor = uv_cor[:, 0, :]
    v_cor = uv_cor[:, 1, :]

    return (u_cor, v_cor)
