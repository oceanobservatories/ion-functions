#!/usr/bin/env python
"""
Module containing calculations related to instruments in the Seafloor
Pressure (PRS) family.

Functions are organized by instrument subsystem:

  BOTTILT -- Seafloor High-Resolution Tilt (Applied Geomechanics LILY):

    prs_bottilt_ccmp -- computes BOTTILT-CCMP_L1
    prs_bottilt_tmag -- computes BOTTILT-TMAG_L1
    prs_bottilt_tdir -- computes BOTTILT-TDIR_L1

  BOTSFLU -- Seafloor Uplift and Subsidence:

    prs_botsflu_time15s   -- computes TIME15S-AUX
    prs_botsflu_meanpres  -- computes BOTSFLU-MEANPRES_L2
    prs_botsflu_predtide  -- computes BOTSFLU-PREDTIDE_L2
    prs_botsflu_meandepth -- computes BOTSFLU-MEANDEPTH_L2
    prs_botsflu_5minrate  -- computes BOTSFLU-5MINRATE_L2
    prs_botsflu_10minrate -- computes BOTSFLU-10MINRATE_L2
    prs_botsflu_time24h   -- computes TIME24H-AUX
    prs_botsflu_daydepth  -- computes BOTSFLU-DAYDEPTH_L2
    prs_botsflu_4wkrate   -- computes BOTSFLU-4WKRATE_L2
    prs_botsflu_8wkrate   -- computes BOTSFLU-8WKRATE_L2

  BOTSFLU helper functions:

    prs_botsflu_daydepth_from_15s_meandepth
    prs_botsflu_4wkrate_from_daydepth
    prs_botsflu_8wkrate_from_daydepth
    anchor_bin_raw_data_to_15s
    anchor_bin_detided_data_to_24h
    calc_meandepth_plus
    calculate_sliding_means
    calculate_sliding_slopes

  Deprecated functions (retained for reference):

    prs_tsunami_detection
    prs_eruption_imminent
    prs_eruption_occurred
    anchor_bin
    calc_daydepth_plus
    calculate_all_sliding_slopes_then_Nan
    calculate_sliding_slopes__MoorePenrose

Authors: Russell Desiderio, Christopher Wingard
"""
import pkg_resources
import numpy as np
import scipy.io as spio
from scipy import signal

from ion_functions import deprecated


#**********************************************************************
#.. BOTTILT: Core functions
#**********************************************************************
def prs_bottilt_ccmp(scmp, sn):
    """
    Compute the corrected compass direction BOTTILT-CCMP_L1 [degrees].

    Applies a sensor-specific lookup table to correct the L0 compass
    reading for calibration offsets and magnetic declination at Axial
    Seamount, as specified in DPS 1341-00060.

    Parameters
    ----------
    scmp : array_like
        Uncorrected sensor compass direction (BOTTILT-SCMP_L0) [degrees].
    sn : array_like
        LILY tiltmeter serial number for each sample [unitless].

    Returns
    -------
    ccmp : ndarray
        Corrected compass direction (BOTTILT-CCMP_L1) [integer degrees
        CW from north].

    Notes
    -----
    The lookup table is stored in prs_functions_ccmp.py as cmp_lookup,
    keyed on (serial_number, rounded_scmp) pairs. The L0 compass value
    is rounded to the nearest integer before lookup, as specified in
    the DPS.
    """
    """
        Currently, there are two coded algorithms:
            (1) the straightforward original, which uses a two-element keyed
                dictionary;
            (2) a faster version, which uses serial number keys to the
                dictionary.

        Since each algorithm uses its own dictionary, the corresponding
        import statements are TEMPORARILY placed at the beginning of their
        respective code sections instead of at module top.
    """
    ###  Original coding, using a dictionary constructed with 2-element keys.

    # load the corrected compass directions table [(sn, scmp) keys]
    from ion_functions.data.prs_functions_ccmp import cmp_lookup

    # use the lookup table to get the ccmp
    ccmp = np.zeros(len(scmp))

    for i in range(len(scmp)):
        sn_key = sn[i].decode('utf-8') if isinstance(sn[i], bytes) else sn[i]
        ccmp[i] = cmp_lookup[(sn_key, int(round(scmp[i])))]
    return ccmp


    ####  Faster coding, using a dictionary constructed with 1-element keys.
    #
    ## load the corrected compass directions table [sn keys]
    #from ion_functions.data.prs_functions_ccmp_lily_compass_cals import cmp_cal
    #
    ## initialize output array for vectorized masking operations. this will 'break'
    ##    the code if an invalid serial number is specified in the argument list.
    #ccmp = np.zeros(len(scmp)) + np.nan
    #
    ## round the uncorrected compass values to the nearest integer as specified in the DPS,
    ##    which uses a lookup table consisting of integral values to do the correction.
    #scmp = np.round(scmp)
    #
    ## find the supported tilt sensor serial numbers, which are keys in the dictionary
    #sernum = cmp_cal.keys()
    #
    #for ii in range(len(sernum)):
    #    # get the cal coeffs as a function of the iterated serial number;
    #    #    x is the raw, uncorrected reading (scmp)
    #    #    y is the corrected reading (ccmp)
    #    [x, y] = cmp_cal[sernum[ii]]
    #
    #    # the boolean mask has 'true' entries where the elements of input vector sn
    #    #    agree with the iterated serial number.
    #    # np.core.defchararray.equal handles vector string comparisons.
    #    mask = np.core.defchararray.equal(sn, sernum[ii])
    #
    #    ## np.interp is used to do the 'lookup' for performance reasons (vectorized)
    #    ccmp[mask] = np.interp(scmp[mask], x, y)
    #
    ## round to make sure we get an integral value (but not int type)
    #return np.round(ccmp)


def prs_bottilt_tmag(x_tilt, y_tilt):
    """
    Compute the resultant tilt magnitude BOTTILT-TMAG_L1 [microradians].

    Computes the vector magnitude of the X- and Y-tilt L0 components as
    specified in DPS 1341-00060.

    Parameters
    ----------
    x_tilt : array_like
        Sensor X-tilt (BOTTILT-XTLT_L0) [microradians].
    y_tilt : array_like
        Sensor Y-tilt (BOTTILT-YTLT_L0) [microradians].

    Returns
    -------
    tmag : ndarray
        Resultant tilt magnitude (BOTTILT-TMAG_L1) [microradians].
    """
    tmag = np.sqrt(x_tilt**2 + y_tilt**2)
    return tmag


def prs_bottilt_tdir(x_tilt, y_tilt, ccmp):
    """
    Compute the resultant tilt direction BOTTILT-TDIR_L1 [degrees].

    Computes the azimuth of the downward tilt vector using arctan2 and
    the corrected compass direction, as specified in DPS 1341-00060.

    Parameters
    ----------
    x_tilt : array_like
        Sensor X-tilt (BOTTILT-XTLT_L0) [microradians].
    y_tilt : array_like
        Sensor Y-tilt (BOTTILT-YTLT_L0) [microradians].
    ccmp : array_like
        Corrected compass direction (BOTTILT-CCMP_L1) [degrees].

    Returns
    -------
    tdir : ndarray
        Resultant tilt direction (BOTTILT-TDIR_L1) [integer degrees
        CW from north].

    Notes
    -----
    The addend of 450 (= 90 + 360) ensures a positive argument to the
    modulo operation before rounding to integer degrees.
    """
    ### As originally coded, according to the algorithm specified in the DPS:

    ## Calculate the angle to use in the tilt direction formula
    ## default angle calculation -- in degrees
    #angle = ne.evaluate('arctan(y_tilt / x_tilt)')
    #angle = np.degrees(angle)
    #
    ## if X-Tilt == 0 and Y-Tilt > 0
    #mask = np.logical_and(x_tilt == 0, y_tilt > 0)
    #angle[mask] = 90.0
    #
    ## if X-Tilt == 0 and Y-Tilt < 0
    #mask = np.logical_and(x_tilt == 0, y_tilt < 0)
    #angle[mask] = -90.0
    #
    ## if Y-Tilt == 0
    #mask = np.equal(y_tilt, np.zeros(len(y_tilt)))
    #angle[mask] = 0.0
    #
    ### Calculate the tilt direction, using the X-Tilt to set the equation
    ## default tilt direction equation
    #tdir = ne.evaluate('(270 - angle + ccmp) % 360')
    #
    ## if X-Tilt >= 0
    #tmp = ne.evaluate('(90 - angle + ccmp) % 360')
    #mask = np.greater_equal(x_tilt, np.zeros(len(x_tilt)))
    #tdir[mask] = tmp[mask]
    #
    #return np.round(tdir)

    # The calculation is faster and simpler if the arctan2 function is used.
    # Use 450=90+360 as the addend in the first argument to the mod(x,360)
    # function to make sure the result is positive.
    return np.round(np.mod(450 - np.degrees(np.arctan2(y_tilt, x_tilt)) + ccmp, 360))


#**********************************************************************
#.. BOTSFLU: Core functions
#**********************************************************************
def prs_botsflu_time15s(timestamp, botpres):
    """
    Compute the auxiliary BOTSFLU timestamp product TIME15S-AUX.

    Returns timestamps anchored at multiples of 15 seconds past the
    minute, corresponding to the time base for the BOTSFLU data products
    binned on 15-second intervals.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia]. Used to discard timestamps
        associated with bad (NaN or non-positive) raw pressure values.

    Returns
    -------
    time15s : ndarray
        15-second-anchored timestamps (TIME15S-AUX) [sec since
        1900-01-01].

    Notes
    -----
    The BOTSFLU data products on this time base are: MEANPRES, PREDTIDE,
    MEANDEPTH, 5MINRATE, and 10MINRATE.
    """
    # botpres is required to eliminate timestamps of bad input values
    time15s, _, _ = anchor_bin_raw_data_to_15s(timestamp, botpres)

    return time15s


def prs_botsflu_meanpres(timestamp, botpres):
    """
    Compute the BOTSFLU mean pressure product BOTSFLU-MEANPRES_L2 [psi].

    Bins raw 20 Hz bottom pressure data into 15-second mean values.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].

    Returns
    -------
    meanpres : ndarray
        15-second mean bottom pressure (BOTSFLU-MEANPRES_L2) [psi].

    Notes
    -----
    The associated time base data product is TIME15S-AUX.
    """
    _, meanpres, _ = anchor_bin_raw_data_to_15s(timestamp, botpres)

    return meanpres


def prs_botsflu_predtide(time):
    """
    Compute the BOTSFLU predicted tide product BOTSFLU-PREDTIDE_L2 [m].

    Assigns predicted tide heights for the three BOTPT instrument sites
    at Axial Seamount by positional indexing into a precomputed lookup
    table.

    Parameters
    ----------
    time : array_like
        BOTSFLU 15-second timestamps (TIME15S-AUX) [sec since
        1900-01-01].

    Returns
    -------
    tide : ndarray
        Predicted tide height (BOTSFLU-PREDTIDE_L2) [m].

    Notes
    -----
    The lookup table is stored in
    prs_functions_tides_2014_thru_2019.mat and contains tide values
    every 15 seconds from 2014-01-01 to 2020-01-01 at lat = 45.95547,
    lon = -130.00957 (Axial Seamount caldera center). Tide values were
    computed using the Tide Model Driver software with the TPXO7.2
    global tidal model. The three BOTPT sites are sufficiently close
    together that a single location is used for all. Tide values are
    stored as signed 4-byte integers in units of 0.001 mm and scaled to
    meters on read.

    A separate unit-test table covering February-April 2011 is loaded
    automatically when the input timestamps predate 2014-01-01.
    """
    time0 = 3597523200.0  # midnight, 2014-01-01
    time_interval = 15.0  # seconds

    # for unit test data, only, feb-apr 2011
    if time[0] < time0:
        time0 = 3502828800.0  # midnight, 2011-01-01
        matpath = 'data/matlab_scripts/botpt/tides_15sec_2011_for_unit_tests.mat'
    else:
        # else, OOI data from 2014 onwards
        matpath = 'data/prs_functions_tides_2014_thru_2019.mat'

    matstream = pkg_resources.resource_stream('ion_functions', matpath)
    dict_tides = spio.loadmat(matstream)
    # tide values are signed 4 byte integers, units [0.001mm]
    tidevector = 0.000001 * dict_tides['tides_mat']
    tidevector = tidevector.reshape((-1))
    # calculate tide vector index as a function of timestamp
    idx = np.around((time - time0) / time_interval)
    tide = tidevector[idx.astype(int)]

    return tide


def prs_botsflu_meandepth(timestamp, botpres, predtide):
    """
    Compute the BOTSFLU de-tided depth product BOTSFLU-MEANDEPTH_L2 [m].

    Converts 15-second mean pressure to depth and removes the predicted
    tidal signal.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].
    predtide : array_like
        Predicted tide (BOTSFLU-PREDTIDE_L2) [m].

    Returns
    -------
    meandepth : ndarray
        De-tided bottom depth (BOTSFLU-MEANDEPTH_L2) [m].

    Notes
    -----
    The associated time base data product is TIME15S-AUX.

    Atmospheric pressure is not subtracted from the L1 pressure data
    even though its units are psia, as specified in DPS 1341-00080.
    """
    _, meandepth, _ = calc_meandepth_plus(timestamp, botpres, predtide)

    return meandepth


def prs_botsflu_5minrate(timestamp, botpres, predtide):
    """
    Compute the BOTSFLU 5-minute rate product BOTSFLU-5MINRATE_L2
    [cm/min].

    Computes the instantaneous rate of depth change using a 5-minute
    backwards-looking difference of the de-tided depth record.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].
    predtide : array_like
        Predicted tide (BOTSFLU-PREDTIDE_L2) [m].

    Returns
    -------
    botsflu_5minrate : ndarray
        5-minute instantaneous depth change rate
        (BOTSFLU-5MINRATE_L2) [cm/min].

    Notes
    -----
    The associated time base data product is TIME15S-AUX.

    For 15-second binned data, 5 minutes corresponds to a lag of 20
    intervals. The conversion from m/5min to cm/min is a factor of 20.
    NaN values arising from data gaps are propagated through the
    difference but are then removed at data dropout positions to maintain
    1:1 correspondence with the TIME15S time base.
    """
    # calculate de-tided depth and the positions of non-zero bins in the original data.
    _, meandepth, mask_nonzero = calc_meandepth_plus(timestamp, botpres, predtide)

    # initialize data product including elements representing data gap positions
    botsflu_5minrate = np.zeros(mask_nonzero.size) + np.nan

    # re-constitute the original data, with data gaps represented by nans.
    data_w_gaps = np.copy(botsflu_5minrate)
    data_w_gaps[mask_nonzero] = meandepth

    # for 15s binned data, 5 minutes comes out to (5 minutes)/(0.25 min) = 20 intervals
    shift = 20
    # units of the subtraction are meter/5min; to convert to cm/min,
    # multiply by 100cm/m and divide by 5 = 20.
    botsflu_5minrate[shift:] = 20.0 * (data_w_gaps[shift:] - data_w_gaps[:-shift])

    # this rate product now has potentially two sources of nans;
    # definitely those at the start of the data record, and any that might
    # have been propagated into the calculation because of the presence of
    # data gaps. remove those only at the data dropout positions (if present)
    # so that this data product will have a 1:1 correspondence with
    # its associated timestamp variable (TIME15S).
    botsflu_5minrate = botsflu_5minrate[mask_nonzero]

    return botsflu_5minrate


def prs_botsflu_10minrate(timestamp, botpres, predtide):
    """
    Compute the BOTSFLU 10-minute rate product BOTSFLU-10MINRATE_L2
    [cm/hr].

    Computes the mean seafloor uplift rate using a 10-minute
    backwards-looking sliding mean of the de-tided depth record,
    differenced over 10 minutes.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].
    predtide : array_like
        Predicted tide (BOTSFLU-PREDTIDE_L2) [m].

    Returns
    -------
    botsflu_10minrate : ndarray
        10-minute mean depth change rate (BOTSFLU-10MINRATE_L2) [cm/hr].

    Notes
    -----
    The associated time base data product is TIME15S-AUX.

    For 15-second binned data, 10 minutes corresponds to a window of 40
    intervals. Sliding means are computed by digital convolution. The
    conversion from m/10min to cm/hr is a factor of 600.
    """
    # calculate de-tided depth and the positions of non-zero bins in the original data.
    _, meandepth, mask_nonzero = calc_meandepth_plus(timestamp, botpres, predtide)

    # initialize data product including elements representing data gap positions
    botsflu_10minrate = np.zeros(mask_nonzero.size) + np.nan

    # re-constitute the original data, with data gaps represented by nans.
    data_w_gaps = np.copy(botsflu_10minrate)
    data_w_gaps[mask_nonzero] = meandepth

    # now calculate sliding 10 minute means.
    # the mean of the 1st 40 values will be located at timestamp position 20
    # (python index 19).
    window_size = 40  # 10min averages on 0.25min binned data
    means = calculate_sliding_means(data_w_gaps, window_size)

    # as above, 10 minutes = 40 intervals for 15sec binned data.
    shift = 40
    # units of the subtraction are meter/10min; to convert to cm/hr,
    # multiply by 100cm/m and multiply by 6 = 600.
    botsflu_10minrate[shift:] = 600.0 * (means[shift:] - means[:-shift])

    # this rate product now has potentially two sources of nans;
    # definitely those at the start of the data record, and any that might
    # have been propagated into the calculation because of the presence of
    # data gaps. remove those only at the data dropout positions (if present)
    # so that this data product will have a 1:1 correspondence with
    # its associated timestamp variable (TIME15S).
    botsflu_10minrate = botsflu_10minrate[mask_nonzero]

    return botsflu_10minrate


def prs_botsflu_time24h(time15s):
    """
    Compute the auxiliary BOTSFLU timestamp product TIME24H-AUX.

    Returns timestamps anchored at midnight, corresponding to the time
    base for BOTSFLU data products binned on 24-hour (noon-to-noon)
    intervals.

    Parameters
    ----------
    time15s : array_like
        15-second-anchored timestamps (TIME15S-AUX) [sec since
        1900-01-01].

    Returns
    -------
    time24h : ndarray
        Midnight-anchored daily timestamps (TIME24H-AUX) [sec since
        1900-01-01].

    Notes
    -----
    The BOTSFLU data products on this time base are: DAYDEPTH, 4WKRATE,
    and 8WKRATE. The time base spans the entire dataset including data
    gaps.
    """
    # the second and third calling arguments are placeholders
    time24h, _, _ = anchor_bin_detided_data_to_24h(time15s, None, None)

    return time24h


def prs_botsflu_daydepth(timestamp, botpres, predtide, dday_coverage=0.90):
    """
    Compute the BOTSFLU daily depth product BOTSFLU-DAYDEPTH_L2 [m].

    Bins 15-second de-tided depth data into 24-hour (noon-to-noon) mean
    values anchored at midnight.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].
    predtide : array_like
        Predicted tide (BOTSFLU-PREDTIDE_L2) [m].
    dday_coverage : float, optional
        Fractional coverage threshold. Daily bins with fewer than this
        fraction of the maximum 5760 possible 15-second values are set
        to NaN. Default is 0.90.

    Returns
    -------
    daydepth : ndarray
        Daily mean de-tided bottom depth (BOTSFLU-DAYDEPTH_L2) [m].

    Notes
    -----
    The associated time base data product is TIME24H-AUX. All days in
    the record span are represented, including data-gap days.
    """
    # calculate 15sec bin timestamps and de-tided depth.
    time15s, meandepth, _ = calc_meandepth_plus(timestamp, botpres, predtide)

    # downstream data products no longer require the mask_nonzero variable
    return prs_botsflu_daydepth_from_15s_meandepth(time15s, meandepth, dday_coverage)


def prs_botsflu_4wkrate(timestamp, botpres, predtide, dday_coverage=0.9, rate_coverage=0.75):
    """
    Compute the BOTSFLU 4-week rate product BOTSFLU-4WKRATE_L2 [cm/yr].

    Computes the mean rate of seafloor depth change using 4-week
    backwards-looking linear regressions on the daily depth record.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].
    predtide : array_like
        Predicted tide (BOTSFLU-PREDTIDE_L2) [m].
    dday_coverage : float, optional
        Fractional coverage threshold for daily depth bins. Default is
        0.9.
    rate_coverage : float, optional
        Fractional window fill threshold for rate calculation. Windows
        below this fraction are set to NaN. Default is 0.75.

    Returns
    -------
    botsflu_4wkrate : ndarray
        4-week mean seafloor depth change rate (BOTSFLU-4WKRATE_L2)
        [cm/yr].

    Notes
    -----
    The associated time base data product is TIME24H-AUX. Regression
    slopes in m/day are converted to cm/yr by multiplying by 36500.
    """
    # calculate daydepth
    daydepth = prs_botsflu_daydepth(timestamp, botpres, predtide, dday_coverage)

    return prs_botsflu_4wkrate_from_daydepth(daydepth, rate_coverage)


def prs_botsflu_8wkrate(timestamp, botpres, predtide, dday_coverage=0.9, rate_coverage=0.75):
    """
    Compute the BOTSFLU 8-week rate product BOTSFLU-8WKRATE_L2 [cm/yr].

    Computes the mean rate of seafloor depth change using 8-week
    backwards-looking linear regressions on the daily depth record.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].
    predtide : array_like
        Predicted tide (BOTSFLU-PREDTIDE_L2) [m].
    dday_coverage : float, optional
        Fractional coverage threshold for daily depth bins. Default is
        0.9.
    rate_coverage : float, optional
        Fractional window fill threshold for rate calculation. Windows
        below this fraction are set to NaN. Default is 0.75.

    Returns
    -------
    botsflu_8wkrate : ndarray
        8-week mean seafloor depth change rate (BOTSFLU-8WKRATE_L2)
        [cm/yr].

    Notes
    -----
    The associated time base data product is TIME24H-AUX. Regression
    slopes in m/day are converted to cm/yr by multiplying by 36500.
    """
    # calculate daydepth
    daydepth = prs_botsflu_daydepth(timestamp, botpres, predtide, dday_coverage)

    return prs_botsflu_8wkrate_from_daydepth(daydepth, rate_coverage)


#**********************************************************************
#.. BOTSFLU: Helper functions
#**********************************************************************
def prs_botsflu_daydepth_from_15s_meandepth(time15s, meandepth, dday_coverage=0.90):
    """
    Compute BOTSFLU-DAYDEPTH_L2 [m] from 15-second binned meandepth.

    Bins the 15-second de-tided depth record into 24-hour mean values
    anchored at midnight. Exposed as a public entry point to allow
    computation of DAYDEPTH directly from pre-computed MEANDEPTH.

    Parameters
    ----------
    time15s : array_like
        15-second-anchored timestamps (TIME15S-AUX) [sec since
        1900-01-01].
    meandepth : array_like
        De-tided 15-second mean depth (BOTSFLU-MEANDEPTH_L2) [m].
    dday_coverage : float, optional
        Fractional coverage threshold. Daily bins with fewer than this
        fraction of the maximum 5760 possible values are set to NaN.
        Default is 0.90.

    Returns
    -------
    daydepth : ndarray
        Daily mean de-tided bottom depth (BOTSFLU-DAYDEPTH_L2) [m].

    Notes
    -----
    The associated time base data product is TIME24H-AUX.
    """
    # bin the 15sec data into 24 hour bins so that the timestamps are at midnight.
    # to calculate daydepth, don't need the time24h timestamps.
    _, daydepth, _ = anchor_bin_detided_data_to_24h(time15s, meandepth, dday_coverage)

    return daydepth


def prs_botsflu_4wkrate_from_daydepth(daydepth, rate_coverage=0.75):
    """
    Compute BOTSFLU-4WKRATE_L2 [cm/yr] from daily depth.

    Applies 4-week backwards-looking linear regression to the daily
    depth record. Exposed as a public entry point to allow computation
    of 4WKRATE directly from pre-computed DAYDEPTH.

    Parameters
    ----------
    daydepth : array_like
        Daily mean de-tided bottom depth (BOTSFLU-DAYDEPTH_L2) [m].
    rate_coverage : float, optional
        Fractional window fill threshold. Windows below this fraction
        are set to NaN. Default is 0.75.

    Returns
    -------
    botsflu_4wkrate : ndarray
        4-week mean seafloor depth change rate (BOTSFLU-4WKRATE_L2)
        [cm/yr].

    Notes
    -----
    The associated time base data product is TIME24H-AUX. The window
    size of 29 days is used. Regression slopes in m/day are converted
    to cm/yr by multiplying by 36500.
    """
    # 4 weeks of data
    window_size = 29
    botsflu_4wkrate = calculate_sliding_slopes(daydepth, window_size, rate_coverage)
    #  convert units:
    #    the units of the slopes are [y]/[x] = meters/day;
    #    to get units of cm/yr, multiply by 100cm/m * 365 days/yr
    botsflu_4wkrate = 100.0 * 365.0 * botsflu_4wkrate

    return botsflu_4wkrate


def prs_botsflu_8wkrate_from_daydepth(daydepth, rate_coverage=0.75):
    """
    Compute BOTSFLU-8WKRATE_L2 [cm/yr] from daily depth.

    Applies 8-week backwards-looking linear regression to the daily
    depth record. Exposed as a public entry point to allow computation
    of 8WKRATE directly from pre-computed DAYDEPTH.

    Parameters
    ----------
    daydepth : array_like
        Daily mean de-tided bottom depth (BOTSFLU-DAYDEPTH_L2) [m].
    rate_coverage : float, optional
        Fractional window fill threshold. Windows below this fraction
        are set to NaN. Default is 0.75.

    Returns
    -------
    botsflu_8wkrate : ndarray
        8-week mean seafloor depth change rate (BOTSFLU-8WKRATE_L2)
        [cm/yr].

    Notes
    -----
    The associated time base data product is TIME24H-AUX. The window
    size of 57 days is used. Regression slopes in m/day are converted
    to cm/yr by multiplying by 36500.
    """
    # 8 weeks of data
    window_size = 57
    botsflu_8wkrate = calculate_sliding_slopes(daydepth, window_size, rate_coverage)
    #  convert units:
    #    the units of the slopes are [y]/[x] = meters/day;
    #    to get units of cm/yr, multiply by 100cm/m * 365 days/yr
    botsflu_8wkrate = 100.0 * 365.0 * botsflu_8wkrate

    return botsflu_8wkrate


def anchor_bin_raw_data_to_15s(time, data):
    """
    Bin raw 20 Hz BOTPT pressure data into 15-second anchored means.

    Calculates anchored bin timestamps and mean-binned data using
    numpy.bincount. Discards NaN values and non-physical (<=0) pressure
    readings before binning. Empty bins are not represented in the
    output; the boolean mask records the positions of non-empty bins
    within the full time span.

    Parameters
    ----------
    time : array_like
        1D array of system timestamps [sec since 1900-01-01].
    data : array_like
        1D array of data values to be binned (BOTPRES_L1 [psia]).

    Returns
    -------
    bin_timestamps : ndarray
        1D array of centered timestamps for non-empty bins [sec since
        1900-01-01].
    binned_data : ndarray
        1D array of mean-binned data; no empty bins are represented.
    mask_nonzero : ndarray
        Boolean array where True values indicate non-empty bin positions
        within the full time span (in 15-second units).

    Notes
    -----
    Bin duration is hard-coded to 15 seconds. Timestamps are anchored
    at the quarter-minute (0, 15, 30, 45 seconds past the minute); each
    bin encompasses data 7.5 seconds on either side of the center.

    The binning follows the numpy.bincount accumarray pattern: elapsed
    time is floored in units of bin duration to assign each sample to a
    bin index, then bincount sums the weighted values.
    """
    bin_duration = 15.0  # seconds
    half_bin = bin_duration/2.0

    # one nan value in a bin will nan out the sum of the values of that bin.
    # throw out nan values and their associated timestamps.
    # also throw out non-physical data values <= 0.
    data[np.isnan(data)] = -999.0  # use any negative value
    mask = (data > 0.0)
    time = time[mask]

    # anchor time-centered bins by determining the start time to be half a bin
    # before the first 'anchor timestamp', which will be an integral number of
    # bin_durations after midnight.
    start_time = np.floor((time[0] - half_bin)/bin_duration) * bin_duration + half_bin
    # calculate elapsed time from start in units of bin_duration.
    elapsed_time = time
    elapsed_time -= start_time
    elapsed_time /= bin_duration
    # assign each timestamp a bin number index based on its elapsed time.
    bin_number = np.floor(elapsed_time).astype(int)
    # the number of elements in each bin is given by
    bin_count = np.bincount(bin_number).astype(float)
    # create a logical mask of non-zero bin_count values
    mask_nonzero = (bin_count != 0)

    # directly calculate bin timestamp, units of [sec]:
    # the midpoint of the data interval is used.
    bin_timestamps = start_time + half_bin + bin_duration * np.arange(bin_count.size)
    # keep only the bins with values
    bin_timestamps = bin_timestamps[mask_nonzero]

    data = data[mask]
    # sum the values in each time bin, and put into the variable binned_data
    binned_data = np.bincount(bin_number, data)
    # divide the values in non-empty bins by the number of values in each bin
    binned_data = binned_data[mask_nonzero]/bin_count[mask_nonzero]

    return bin_timestamps, binned_data, mask_nonzero


def anchor_bin_detided_data_to_24h(time, data, dday_coverage):
    """
    Bin 15-second BOTSFLU depth data into 24-hour anchored means.

    Calculates midnight-anchored bin timestamps and mean-binned data
    using numpy.bincount. All days in the record span are represented,
    including days below the coverage threshold (which receive NaN).

    Parameters
    ----------
    time : array_like
        1D array of 15-second timestamps [sec since 1900-01-01].
    data : array_like
        1D array of de-tided depth values to be binned (MEANDEPTH [m]).
    dday_coverage : float or None
        Fractional coverage threshold. Bins with fewer than this
        fraction of the maximum 5760 possible values are set to NaN.
        Pass None to skip threshold masking (used when computing
        timestamps only).

    Returns
    -------
    bin_timestamps : ndarray
        1D array of midnight-anchored daily timestamps [sec since
        1900-01-01].
    daydepth : ndarray
        1D array of daily mean depth values (BOTSFLU-DAYDEPTH_L2) [m].
        Bins below the coverage threshold contain NaN.
    raw_bincount : ndarray
        1D array of the raw count of values in each bin (used as a
        diagnostic in unit tests).

    Notes
    -----
    Bin duration is hard-coded to 86400 seconds (one day). Each bin
    spans noon-to-noon, anchored at midnight. The maximum possible bin
    count for 15-second data within a 24-hour bin is 5760.

    The binning follows the numpy.bincount accumarray pattern.
    """
    bin_duration = 86400.0  # number of seconds in a day
    half_bin = bin_duration/2.0
    max_count = bin_duration/15.0  # maximum number of values in a day's bin

    # anchor time-centered bins by determining the start time to be half a bin
    # before the first 'anchor timestamp', which will an integral number of
    # bin_durations after midnight.
    start_time = np.floor((time[0] - half_bin)/bin_duration) * bin_duration + half_bin
    # calculate elapsed time from start in units of bin_duration.
    time_elapsed = (time - start_time)/bin_duration
    # assign each timestamp a bin number index based on its elapsed time.
    bin_number = np.floor(time_elapsed).astype(int)
    # the number of elements in each bin is given by
    bin_count = np.bincount(bin_number).astype(float)
    raw_bincount = np.copy(bin_count)  # for unit test

    # bin_count is used as a divisor to calculate mean values at each bin
    #    bins with bincounts below the threshold value will have a nan value.
    #    bins with bincounts equal to or above the threshold value are non_Nan.
    bin_mean = bin_count//max_count
    if np.isscalar(dday_coverage):
        bin_count[bin_mean < dday_coverage] = np.nan
    #    use nans to prevent dividing by zero in case dday_coverage=0
    bin_count[bin_count == 0] = np.nan

    # directly calculate bin timestamp, units of [sec]:
    # the midpoint of the data interval is used.
    bin_timestamps = start_time + half_bin + bin_duration * np.arange(bin_count.size)

    # sum the values in each time bin, and put into the variable binned_data
    binned_data = np.bincount(bin_number, data)
    # divide the values in each bin by the number of values in each bin
    daydepth = binned_data/bin_count

    return bin_timestamps, daydepth, raw_bincount


def calc_meandepth_plus(timestamp, botpres, predtide):
    """
    Compute BOTSFLU-MEANDEPTH_L2 [m] plus auxiliary binning outputs.

    Bins raw pressure data to 15-second means, converts to de-tided
    depth, and returns intermediate variables needed by downstream
    BOTSFLU data product functions.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].
    predtide : array_like
        Predicted tide (BOTSFLU-PREDTIDE_L2) [m].

    Returns
    -------
    time15s : ndarray
        15-second-anchored timestamps [sec since 1900-01-01].
    meandepth : ndarray
        De-tided 15-second mean bottom depth
        (BOTSFLU-MEANDEPTH_L2) [m].
    mask_nonzero : ndarray
        Boolean array marking positions of non-empty 15-second bins
        within the full time span.

    Notes
    -----
    Atmospheric pressure is not subtracted from the L1 pressure data
    even though its units are psia, as specified in DPS 1341-00080.
    The conversion factor from psi to depth is -0.67 m/psi. The
    negative sign reflects the convention that depths are negative,
    so the predicted tide (positive) is added to de-tide the record.
    """
    # The pressure values do have units of psia. However, historically at these sites
    # atmospheric pressure has *not* been subtracted when converting the pressure data
    # to depth. Therefore the DPS authors do not want atmospheric pressure subtracted
    # in the DPA. To emphasize this, I have created the variable atm_press_psi and set
    # it to 0.
    atm_press_psi = 0.0
    psi_2_depth = -0.67  # psi to depth in meters
    bin_duration = 15.0  # seconds

    time15s, meanpres, mask_nonzero = anchor_bin_raw_data_to_15s(timestamp, botpres)

    # de-tide
    meandepth = ((meanpres - atm_press_psi) * psi_2_depth) + predtide

    # downstream data products require the time15s and mask_nonzero variables,
    # so pass these as output arguments so that they won't have to be recalculated.
    return time15s, meandepth, mask_nonzero


def calculate_sliding_means(data, window_size):
    """
    Compute time-centered sliding means by digital convolution.

    Used internally for the BOTSFLU 10MINRATE data product.

    Parameters
    ----------
    data : array_like
        1D array of input data.
    window_size : int
        Number of samples in the sliding window.

    Returns
    -------
    means : ndarray
        1D array of sliding window means; boundary elements are set to
        NaN.

    Notes
    -----
    For even window sizes, the result is rolled by one sample to match
    the behavior of MATLAB's convolution, which was used to generate the
    unit test reference values.
    """
    kk = np.ones(window_size) / window_size
    means = np.convolve(data, kk, 'same')
    # nan out data with boundary effects at edges.
    # integer arithmetic will 'truncate' 5/2 to 2 and -5/2 to -3.
    means[:window_size//2] = np.nan
    means[-((window_size-1)//2):] = np.nan

    # matlab and numpy behave differently for even window sizes, so roll the python
    # result to mimic matlab
    means = np.roll(means, -np.mod(window_size+1, 2))  # roll only if window_size is even

    return means


def calculate_sliding_slopes(data, window_size, coverage_threshold=0.0):
    """
    Compute backwards-looking sliding linear regression slopes.

    Used internally for the BOTSFLU 4WKRATE and 8WKRATE data products.
    NaN values within a window are excluded from the regression; windows
    with insufficient valid data are set to NaN.

    Parameters
    ----------
    data : array_like
        1D array of input data (typically daily depth values [m]).
    window_size : int
        Nominal window size in samples. If even, incremented by one to
        enforce an odd window.
    coverage_threshold : float, optional
        Minimum fraction of non-NaN values required within a window to
        compute a slope. Default is 0.0 (no threshold).

    Returns
    -------
    slopes : ndarray
        1D array of backwards-looking regression slopes in units of
        [data units]/[sample]. Length equals len(data).

    Notes
    -----
    The regression equations are from Press et al. (1986), equations
    14.2.15-14.2.17. The data vector is front-padded with NaN so that
    windows near the start of the record that satisfy the coverage
    criterion produce non-NaN values. Coverage is assessed by
    calculate_sliding_means on a binary good/bad mask, then shifted to
    align with the backwards-looking window.
    """
    # ODD WINDOW SIZES are expected; if even, increment by one
    window_size = window_size + 1 - np.mod(window_size, 2)
    half_window = window_size // 2  # this will 'floor' when window_size is odd as desired

    # pad front end of data with Nans (because data product and 'for loop' are backwards-looking)
    npts = 2 * half_window + data.size
    padded_data = np.zeros(npts) + np.nan
    padded_data[window_size-1:] = data

    # first calculate coverage and compare it to threshold so that only windows containing the
    # prescribed number of 'good' data points will have slope values calculated in the for loop.

    y = np.copy(padded_data)
    # determine the fraction of good values per window (vectorized):
    # first change non-nan values to 1, then nan values to 0, and take the average
    y[~np.isnan(y)] = 1.0
    y[np.isnan(y)] = 0.0
    fraction_good = calculate_sliding_means(y, window_size)  # centered windows
    # convert to backwards-looking
    fraction_good = np.roll(fraction_good, half_window)

    # deal with round-off error and boundary considerations at threshold values
    machine_epsilon = np.finfo(float).eps
    eps_x_100 = machine_epsilon * 100

    # rather than deleting padding, keep it for index registration with padded data;
    # the data padding will be deleted at the end of the routine anyway.
    # set nans to no coverage (avoids warning message in conditional statements)
    fraction_good[np.isnan(fraction_good)] = -1.0
    # indices with fraction_good near 0 are to be disregarded
    fraction_good[fraction_good <= eps_x_100] = -1.0

    # there can be issues involving roundoff error when considering 100% coverage, so:
    fraction_good = fraction_good + eps_x_100
    # find indices that satisfy coverage criterion,
    # which is the first element of the 'np.where' tuple
    indices = np.where(fraction_good >= coverage_threshold)[0]

    # actual slope calculation

    slopes = np.zeros(npts) + np.nan            # initialize
    abscissa = np.arange(npts).astype('float')  # unity spacing
    for ii in indices:
        x = abscissa[(ii-2*half_window):(ii+1)]
        y = padded_data[(ii-2*half_window):(ii+1)]
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]
        if x.size < 2:
            continue  # trap out windows with only 0 or 1 valid datapoint
        tti = x - x.mean(axis=0)
        stt = np.sum(tti * tti)
        slopes[ii] = np.sum(tti * y) / stt
    # get rid of padding
    slopes = slopes[2*half_window:]

    return slopes


#**********************************************************************
#.. Deprecated functions
#**********************************************************************
@deprecated
def prs_tsunami_detection(botsflu_5minrate, tsunami_detection_threshold=1.0):
    """
    OOI wrapper returning True if a tsunami event is detected.

    Deprecated. This function was coded from pseudocode in DPS
    1341-00080, which was never publicly released. Its robustness has
    not been verified with actual data.

    Parameters
    ----------
    botsflu_5minrate : array_like
        5-minute instantaneous depth change rate
        (BOTSFLU-5MINRATE_L2) [cm/min].
    tsunami_detection_threshold : float, optional
        Detection threshold [cm/min]. Default is 1.0.

    Returns
    -------
    boolean_tsunami_detection : bool
        True if any value of botsflu_5minrate meets or exceeds the
        threshold in absolute value; False otherwise.
    """
    # units of variable and threshold are [cm/min]
    boolean_tsunami_detection = False
    # get rid of runtime warnings if nans are present
    botsflu_5minrate[np.isnan(botsflu_5minrate)] = 0.0
    if np.any(np.abs(botsflu_5minrate) >= tsunami_detection_threshold):
        boolean_tsunami_detection = True
    return boolean_tsunami_detection


@deprecated
def prs_eruption_imminent(botsflu_10minrate, eruption_imminent_threshold=5.0):
    """
    OOI wrapper returning True if a volcanic eruption is imminent.

    Deprecated. This function was coded from pseudocode in DPS
    1341-00080, which was never publicly released. Its robustness has
    not been verified with actual data.

    Parameters
    ----------
    botsflu_10minrate : array_like
        10-minute mean depth change rate (BOTSFLU-10MINRATE_L2) [cm/hr].
    eruption_imminent_threshold : float, optional
        Detection threshold [cm/hr]. Default is 5.0.

    Returns
    -------
    boolean_eruption_imminent : bool
        True if any value of botsflu_10minrate meets or exceeds the
        threshold; False otherwise.
    """
    # units of variable and threshold are [cm/hr]
    boolean_eruption_imminent = False
    # get rid of runtime warnings if nans are present
    botsflu_10minrate[np.isnan(botsflu_10minrate)] = 0.0
    if np.any(botsflu_10minrate >= eruption_imminent_threshold):
        boolean_eruption_imminent = True
    return boolean_eruption_imminent


@deprecated
def prs_eruption_occurred(botsflu_10minrate, eruption_occurred_threshold=-5.0):
    """
    OOI wrapper returning True if a volcanic eruption has occurred.

    Deprecated. This function was coded from pseudocode in DPS
    1341-00080, which was never publicly released. Its robustness has
    not been verified with actual data.

    Parameters
    ----------
    botsflu_10minrate : array_like
        10-minute mean depth change rate (BOTSFLU-10MINRATE_L2) [cm/hr].
    eruption_occurred_threshold : float, optional
        Detection threshold [cm/hr]. Default is -5.0.

    Returns
    -------
    boolean_eruption_occurred : bool
        True if any value of botsflu_10minrate is at or below the
        threshold; False otherwise.
    """
    # units of variable and threshold are [cm/hr]
    boolean_eruption_occurred = False
    # get rid of runtime warnings if nans are present
    botsflu_10minrate[np.isnan(botsflu_10minrate)] = 0.0
    if np.any(botsflu_10minrate <= eruption_occurred_threshold):
        boolean_eruption_occurred = True
    return boolean_eruption_occurred


@deprecated
def anchor_bin(time, data, bin_duration, mode):
    """
    Bin BOTPT data into anchored time bins (general-purpose).

    Deprecated May 2017. Superseded by anchor_bin_raw_data_to_15s and
    anchor_bin_detided_data_to_24h, which handle the raw-data bad-value
    check and extended 24-hour timestamp records respectively.

    Parameters
    ----------
    time : array_like
        1D array of timestamps [sec since 1900-01-01].
    data : array_like or None
        1D array of data to be binned, or None when mode is 'time'.
    bin_duration : float
        Bin size [s].
    mode : str
        Output mode. One of 'time' (timestamps only), 'data' (binned
        data and mask only), or 'both' (timestamps, binned data, and
        mask).

    Returns
    -------
    Depends on mode:
    'time'  : bin_timestamps
    'data'  : binned_data, mask_nonzero
    'both'  : bin_timestamps, binned_data, mask_nonzero

    Notes
    -----
    Timestamps are anchored so that bin centers fall at integral
    multiples of bin_duration after midnight (e.g., at 0, 15, 30, 45
    seconds past the minute for bin_duration = 15 s).
    """
    half_bin = bin_duration/2.0
    # anchor time-centered bins by determining the start time to be half a bin
    # before the first 'anchor timestamp', which will an integral number of
    # bin_durations after midnight.
    start_time = np.floor((time[0] - half_bin)/bin_duration) * bin_duration + half_bin
    # calculate elapsed time from start in units of bin_duration.
    time_elapsed = (time - start_time)/bin_duration
    # assign each timestamp a bin number index based on its elapsed time.
    bin_number = np.floor(time_elapsed).astype(int)
    # the number of elements in each bin is given by
    bin_count = np.bincount(bin_number).astype(float)
    # create a logical mask of non-zero bin_count values
    mask_nonzero = (bin_count != 0)

    # to calculate timestamps and to get tides, without also binning data.
    # mask_nonzero is not needed.
    if mode == 'time':
        # directly calculate bin timestamp, units of [sec]:
        # the midpoint of the data interval is used.
        bin_timestamps = start_time + half_bin + bin_duration * np.arange(bin_count.size)
        # keep only the bins with values
        bin_timestamps = bin_timestamps[mask_nonzero]
        return bin_timestamps

    # for binning data when the resultant timestamps are not explicitly required.
    # daydepth_plus also requires mask_nonzero for downstream products 4wkrate and 8wkrate.
    elif mode == 'data':
        # sum the values in each time bin, and put into the variable binned_data
        binned_data = np.bincount(bin_number, data)
        # divide the values in non-empty bins by the number of values in each bin
        binned_data = binned_data[mask_nonzero]/bin_count[mask_nonzero]
        return binned_data, mask_nonzero

    # for when both timestamps and binned data are required.
    elif mode == 'both':
        bin_timestamps = start_time + half_bin + bin_duration * np.arange(bin_count.size)
        bin_timestamps = bin_timestamps[mask_nonzero]
        binned_data = np.bincount(bin_number, data)
        binned_data = binned_data[mask_nonzero]/bin_count[mask_nonzero]
        return bin_timestamps, binned_data, mask_nonzero


@deprecated
def calc_daydepth_plus(timestamp, botpres, predtide):
    """
    Compute BOTSFLU-DAYDEPTH_L2 [m] plus a non-empty bin mask.

    Deprecated May 2017. Superseded by calc_meandepth_plus and
    prs_botsflu_daydepth.

    Parameters
    ----------
    timestamp : array_like
        OOI system timestamps [sec since 1900-01-01].
    botpres : array_like
        Bottom pressure (BOTPRES_L1) [psia].
    predtide : array_like
        Predicted tide [m].

    Returns
    -------
    daydepth : ndarray
        Daily mean de-tided bottom depth (BOTSFLU-DAYDEPTH_L2) [m].
    """
    # calculate 15sec bin timestamps and de-tided depth.
    time15s, meandepth, _ = calc_meandepth_plus(timestamp, botpres, predtide)

    # bin the 15sec data into 24 hour bins so that the timestamps are at midnight.
    # to calculate daydepth, don't need the time24h timestamps.

    _, daydepth = anchor_bin(time15s, meandepth)

    # downstream data products no longer require the mask_nonzero variable
    return daydepth


@deprecated
def calculate_all_sliding_slopes_then_Nan(data, window_size, coverage_threshold):
    """
    Compute backwards-looking sliding slopes, then NaN below coverage.

    Deprecated May 2017. Superseded by calculate_sliding_slopes, which
    pre-filters windows by fractional coverage before computing slopes
    rather than computing all slopes first and then applying NaN masking.

    Parameters
    ----------
    data : array_like
        1D array of input data.
    window_size : int
        Nominal window size in samples. If even, incremented by one.
    coverage_threshold : float
        Minimum fraction of non-NaN values required within a window.
        Windows below this threshold are set to NaN after computation.

    Returns
    -------
    slopes : ndarray
        1D array of backwards-looking regression slopes. Length equals
        len(data).

    Notes
    -----
    The regression equations are from Press et al. (1986), equations
    14.2.15-14.2.17.
    """
    # ODD WINDOW SIZES are expected; if even, increment by one
    window_size = window_size + 1 - np.mod(window_size, 2)
    half_window = window_size / 2  # this will 'floor' when window_size is odd as desired

    # pad front end of data with Nans (because 'for loop' is backwards-looking)
    npts = 2 * half_window + data.size
    padded_data = np.zeros(npts) + np.nan
    padded_data[window_size-1:] = data

    # first calculate values for all sliding windows
    slopes = np.zeros(npts) + np.nan
    abscissa = np.arange(npts).astype('float')
    for ii in range(window_size-1, npts):
        x = abscissa[(ii-2*half_window):(ii+1)]  # rather than np.arange-ing in each iteration
        y = padded_data[(ii-2*half_window):(ii+1)]
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]
        if x.size < 2:
            continue  # trap out windows with only 0 or 1 valid datapoint
        tti = x - x.mean(axis=0)
        stt = np.sum(tti * tti)
        slopes[ii] = np.sum(tti * y) / stt
    # get rid of padding
    slopes = slopes[2*half_window:]

    # now determine the fraction of good values per window (vectorized):
    # first change non-nan values to 1, then nan values to 0, and take the average
    padded_data[~np.isnan(padded_data)] = 1.0
    padded_data[np.isnan(padded_data)] = 0.0
    fraction_good = calculate_sliding_means(padded_data, window_size)  # centered windows
    # convert to backwards-looking
    fraction_good = np.roll(fraction_good, half_window)
    # get rid of padding
    fraction_good = fraction_good[2*half_window:]

    # nan out fractional values (means) less than coverage threshold
    # there can be issues involving roundoff error when considering 100% coverage, so:
    machine_epsilon = np.finfo(float).eps
    fraction_good = fraction_good + 100 * machine_epsilon
    # avoid a python warning message by trapping out nans in the conditional
    fraction_good[np.isnan(fraction_good)] = -999.0  # any negative number will work as intended
    slopes[fraction_good < coverage_threshold] = np.nan

    return slopes


@deprecated
def calculate_sliding_slopes__MoorePenrose(data, window_size):
    # DEPRECATED because nan-masking cannot be implemented with this algorithm #
    """
    Compute backwards-looking sliding slopes via Moore-Penrose pseudoinverse.

    Deprecated May 2017. Superseded by calculate_sliding_slopes. The
    Moore-Penrose pseudoinverse method cannot accommodate NaN values in
    the data window.

    Parameters
    ----------
    data : array_like
        1D array of input data.
    window_size : int
        Number of samples in the sliding window.

    Returns
    -------
    slopes : ndarray
        1D array of backwards-looking regression slopes. The first
        window_size - 1 elements are NaN.

    Notes
    -----
    Algorithm from John D'Errico's response on Matlab Central (thread
    49181). The first non-NaN value occurs at index window_size - 1 and
    is the slope of the regression of the first window_size points. If
    time-centered slopes are needed, circularly shift the result by half
    a window.
    """
    column1 = np.ones((window_size, 1))
    column2 = -np.arange(float(window_size)).reshape(-1, 1)
    X = np.hstack((column1, column2))
    filtercoef = np.linalg.pinv(X)
    slopes = signal.lfilter(filtercoef[1, :], 1, data)
    slopes[0:window_size-1] = np.nan

    # if time-centered slopes are desired, circularly shift this by half a window.

    return slopes
