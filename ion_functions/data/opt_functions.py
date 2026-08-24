#!/usr/bin/env python
"""
ion_functions.data.opt_functions

Functions supporting the OPTAA (Sea-Bird Scientific ac-s), PARAD (scalar
PAR sensors from multiple manufacturers), and SPKIR (Sea-Bird Scientific
OCR-507) optical instrument families. Computes OPTATTN_L2 and OPTABSN_L2
from the ac-s, OPTPARW_L1 from three PARAD variants, and SPECTIR_L1 from
the OCR-507.
"""

import numpy as np

# load the temperature and salinity correction coefficients table
from ion_functions.data.opt_functions_tscor import tscor


# wrapper function to calculate the beam attenuation coefficients (OPTATTN_L2)
#   from the WET Labs, Inc. ACS (OPTAA).
def opt_beam_attenuation(cref, csig, traw, cwl, coff, tcal, tbins, tc_arr,
                         T, PS):

    """
    OOI wrapper for OPTATTN_L2. Computes the beam attenuation coefficient
    corrected for temperature and salinity from the Sea-Bird Scientific
    ac-s (OPTAA).

    Parameters
    ----------
    cref : ndarray
        Raw reference light measurements, OPTCREF_L0 [counts].
    csig : ndarray
        Raw signal light transmission measurements, OPTCSIG_L0
        [counts].
    traw : ndarray
        Raw internal instrument temperature, OPTTEMP_L0 [counts].
    cwl : ndarray
        Attenuation channel wavelengths [nm], from the ac-s device
        file.
    coff : ndarray
        Pure water offsets for the attenuation channel, from the ac-s
        device file [m-1].
    tcal : ndarray
        Factory calibration reference (pure water) temperature
        [deg_C], supplied by the manufacturer.
    tbins : ndarray
        Internal temperature calibration bin values, from the ac-s
        device file [deg_C].
    tc_arr : ndarray
        Internal temperature calibration correction coefficients for
        the attenuation channel, from the ac-s device file [m-1].
    T : ndarray
        In situ temperature, TEMPWAT_L1, from a co-located CTD
        [deg_C].
    PS : ndarray
        In situ practical salinity, PRACSAL_L2, from a co-located CTD
        [unitless].

    Returns
    -------
    cpd_ts : ndarray
        Beam attenuation coefficient corrected for temperature and
        salinity, OPTATTN_L2 [m-1].

    See Also
    --------
    opt_internal_temp : Computes the internal instrument temperature
        used here.
    opt_pd_calc : Computes the uncorrected beam attenuation
        coefficient.
    opt_tempsal_corr : Applies the temperature and salinity
        correction.

    Notes
    -----
    All input arrays are assumed to be vectorized over data packets,
    with the first dimension iterating over packet number; the
    function loops over packets, calling opt_internal_temp,
    opt_pd_calc, and opt_tempsal_corr for each.
    """
    # reset shapes of input arguments
    #    using np.array([], ndmin=#) seems faster than using np.atleast_#d
    cref = np.array(cref, ndmin=2)
    csig = np.array(csig, ndmin=2)
    traw = np.array(traw, ndmin=1)
    cwl = np.around(np.array(cwl, ndmin=2), decimals=1)
    coff = np.array(coff, ndmin=2)
    tcal = np.array(tcal, ndmin=1)
    tbins = np.array(tbins, ndmin=2)
    T = np.array(T, ndmin=1)
    PS = np.array(PS, ndmin=1)
    # note, np.atleast_3d appends the extra dimension;
    # np.array using ndmin prepends the extra dimension.
    tc_arr = np.array(tc_arr, ndmin=3)

    # size up inputs
    npackets = cwl.shape[0]
    nwavelengths = cwl.shape[1]
    # initialize output array
    cpd_ts = np.zeros([npackets, nwavelengths])

    for ii in range(npackets):

        # calculate the internal instrument temperature [deg_C]
        tintrn = opt_internal_temp(traw[ii])

        # calculate the uncorrected beam attenuation coefficient [m^-1]
        cpd, _ = opt_pd_calc(cref[ii, :], csig[ii, :], coff[ii, :], tintrn,
                             tbins[ii, :], tc_arr[ii, :, :])

        # correct the beam attenuation coefficient for temperature and salinity.
        cpd_ts_row = opt_tempsal_corr('c', cpd, cwl[ii, :], tcal[ii], T[ii], PS[ii])
        cpd_ts[ii, :] = cpd_ts_row

    # return the temperature and salinity corrected beam attenuation
    # coefficient OPTATTN_L2 [m^-1]
    return cpd_ts


# wrapper function to calculate the optical absorption coefficients (OPTABSN_L2)
#   from the WET Labs, Inc. ACS (OPTAA).
def opt_optical_absorption(aref, asig, traw, awl, aoff, tcal, tbins, ta_arr,
                           cpd_ts, cwl, T, PS, rwlngth=715.):
    """
    OOI wrapper for OPTABSN_L2. Computes the optical absorption
    coefficient corrected for temperature, salinity, and scattering
    from the Sea-Bird Scientific ac-s (OPTAA).

    Parameters
    ----------
    aref : ndarray
        Raw reference light measurements, OPTAREF_L0 [counts].
    asig : ndarray
        Raw signal light transmission measurements, OPTASIG_L0
        [counts].
    traw : ndarray
        Raw internal instrument temperature, OPTTEMP_L0 [counts].
    awl : ndarray
        Absorption channel wavelengths [nm], from the ac-s device
        file.
    aoff : ndarray
        Pure water offsets for the absorption channel, from the ac-s
        device file [m-1].
    tcal : ndarray
        Factory calibration reference (pure water) temperature
        [deg_C], supplied by the manufacturer.
    tbins : ndarray
        Internal temperature calibration bin values, from the ac-s
        device file [deg_C].
    ta_arr : ndarray
        Internal temperature calibration correction coefficients for
        the absorption channel, from the ac-s device file [m-1].
    cpd_ts : ndarray
        Beam attenuation coefficient corrected for temperature and
        salinity, OPTATTN_L2 [m-1], from opt_beam_attenuation.
    cwl : ndarray
        Attenuation channel wavelengths [nm], from the ac-s device
        file.
    T : ndarray
        In situ temperature, TEMPWAT_L1, from a co-located CTD
        [deg_C].
    PS : ndarray
        In situ practical salinity, PRACSAL_L2, from a co-located CTD
        [unitless].
    rwlngth : float, optional
        Scattering correction reference wavelength [nm] (default
        715).

    Returns
    -------
    apd_ts_s : ndarray
        Absorption coefficient corrected for temperature, salinity,
        and scattering, OPTABSN_L2 [m-1].

    See Also
    --------
    opt_internal_temp : Computes the internal instrument temperature
        used here.
    opt_pd_calc : Computes the uncorrected absorption coefficient.
    opt_tempsal_corr : Applies the temperature and salinity
        correction.
    opt_scatter_corr : Applies the scattering correction.

    Notes
    -----
    All input arrays are assumed to be vectorized over data packets,
    with the first dimension iterating over packet number; the
    function loops over packets, calling opt_internal_temp,
    opt_pd_calc, opt_tempsal_corr, and opt_scatter_corr for each.
    rwlngth defaults to 715 nm per DPS 1341-00700, but is exposed for
    override if needed.
    """
    # reset shapes of input arguments
    #    using np.array ndmin=# seems faster than using np.atleast_#d
    aref = np.array(aref, ndmin=2)
    asig = np.array(asig, ndmin=2)
    traw = np.array(traw, ndmin=1)
    awl = np.around(np.array(awl, ndmin=2), decimals=1)
    aoff = np.array(aoff, ndmin=2)
    tcal = np.array(tcal, ndmin=1)
    tbins = np.array(tbins, ndmin=2)
    # note, np.atleast_3d appends the extra dimension;
    # np.array using ndmin prepends the extra dimension.
    ta_arr = np.array(ta_arr, ndmin=3)
    cpd_ts = np.array(cpd_ts, ndmin=2)
    cwl = np.array(cwl, ndmin=2)
    T = np.array(T, ndmin=1)
    PS = np.array(PS, ndmin=1)

    # size up inputs
    npackets = awl.shape[0]
    nwavelengths = awl.shape[1]
    # initialize output array
    apd_ts_s = np.zeros([npackets, nwavelengths])

    for ii in range(npackets):

        # calculate the internal instrument temperature [deg_C]
        tintrn = opt_internal_temp(traw[ii])

        # calculate the uncorrected optical absorption coefficient [m^-1]
        apd, _ = opt_pd_calc(aref[ii, :], asig[ii, :], aoff[ii, :], tintrn,
                             tbins[ii, :], ta_arr[ii, :, :])

        # correct the optical absorption coefficient for temperature and salinity.
        apd_ts = opt_tempsal_corr('a', apd, awl[ii, :], tcal[ii], T[ii], PS[ii])

        # correct the optical absorption coefficient for scattering effects
        apd_ts_s_row = opt_scatter_corr(apd_ts, awl[ii, :], cpd_ts[ii, :], cwl[ii, :], rwlngth)
        apd_ts_s[ii, :] = apd_ts_s_row

    # return the temperature, salinity and scattering corrected optical
    # absorption coefficient OPTABSN_L2 [m^-1]
    return apd_ts_s


# Functions used in calculating optical absorption and beam attenuation
# coefficients from the OPTAA family of instruments.
def opt_internal_temp(traw):
    """
    Calculates the internal ac-s instrument temperature from the raw
    thermistor count. Used internally by opt_pd_calc.

    Parameters
    ----------
    traw : ndarray
        Raw internal instrument temperature, OPTTEMP_L0 [counts].

    Returns
    -------
    tintrn : ndarray
        Calculated internal instrument temperature [deg_C].
    """
    # convert counts to volts
    volts = 5. * traw / 65535.

    # calculate the resistance of the thermistor
    res = 10000. * volts / (4.516 - volts)

    # convert resistance to temperature
    a = 0.00093135
    b = 0.000221631
    c = 0.000000125741

    log_res = np.log(res)
    degC = (1. / (a + b * log_res + c * log_res**3)) - 273.15
    return degC


def opt_pd_calc(ref, sig, offset, tintrn, tbins, tarray):
    """
    Converts raw reference and signal measurements to an uncorrected
    beam attenuation or optical absorption coefficient. The 'c' (beam
    attenuation) and 'a' (absorption) cases are isomorphic; they differ
    only in which calibration coefficients are supplied. The returned
    values are not final data products.

    Parameters
    ----------
    ref : ndarray
        Raw reference light measurements, OPTCREF_L0 or OPTAREF_L0 as
        appropriate [counts].
    sig : ndarray
        Raw signal light measurements, OPTCSIG_L0 or OPTASIG_L0 as
        appropriate [counts].
    offset : ndarray
        Pure water offsets from the ac-s device file; use the 'c' or
        'a' offsets as appropriate [m-1].
    tintrn : float
        Internal instrument temperature [deg_C], from
        opt_internal_temp.
    tbins : ndarray
        Internal temperature calibration bin values, from the ac-s
        device file [deg_C].
    tarray : ndarray
        Internal temperature calibration correction coefficients from
        the ac-s device file, indexed by wavelength and temperature
        bin; use the 'c' or 'a' array as appropriate [m-1].

    Returns
    -------
    pd : ndarray
        Uncorrected beam attenuation or optical absorption
        coefficient [m-1].
    deltaT : ndarray
        Internal temperature correction applied to pd [m-1]; returned
        for unit testing and not used in subsequent processing.
    """
    # Raw reference and signal values are imported as 1D arrays. They must be
    # the same length.
    ref = np.atleast_1d(ref).astype(float)
    sig = np.atleast_1d(sig).astype(float)
    lFlag = len(ref) != len(sig)
    if lFlag:
        raise ValueError('Reference and Signal arrays must be the same length')

    nValues = len(sig)

    # The offsets are imported as a 1D array. They must be the same length as
    # ref and sig.
    offset = np.atleast_1d(offset)
    lFlag = len(offset) != nValues
    if lFlag:
        str1 = 'The number of calibration offset channels ('
        str2 = ') from the cal devfile must match the number of Signal and Reference channels ('
        str3 = ') from the rawdata.'
        offsetErrorString = str1 + str(len(offset)) + str2 + str(nValues) + str3
        raise ValueError(offsetErrorString)

    # The temperature bins are imported as a 1D array
    tbins = np.atleast_1d(tbins)
    tValues = np.size(tbins)

    # The internal temperature compensation array 'tarray' is a 2D array. The # of
    # columns must equal the length of the tbins array. The number of rows must equal
    # the number of wavelengths.
    tarray = np.atleast_2d(tarray)
    r, c = tarray.shape

    if r != nValues:
        str1 = 'The number of rows in the internal temperature compensation calibration array ('
        str2 = ') must match the number of wavelength channels in the rawdata ('
        str3 = ').'
        rErrorString = str1 + str(r) + str2 + str(nValues) + str3
        raise ValueError(rErrorString)

    if c != tValues:
        # since both tbins and tarray (should) come from the same dev file, this exception is not likely.
        str1 = 'Number of columns in the internal temperature compensation calibration array ('
        str2 = ') must match the number of internal temp comp cal bin values = ('
        cErrorString = str1 + str(c) + str2 + str(tValues) + ').'
        raise ValueError(cErrorString)

    # find the indexes in the temperature bins corresponding to the values
    # bracketing the internal temperature.
    ind1 = np.nonzero(tbins-tintrn < 0)[0][-1]
    ind2 = np.nonzero(tintrn-tbins < 0)[0][0]
    T0 = tbins[ind1]    # set first bracketing temperature
    T1 = tbins[ind2]    # set second bracketing temperaure

    # Calculate the linear temperature correction.
    dT0 = tarray[:, ind1]
    dT1 = tarray[:, ind2]
    deltaT = dT0 + ((tintrn - T0) / (T1 - T0)) * (dT1 - dT0)

    # Calculate the uncorrected signal [m-1]; the pathlength is 0.25m.
    # Apply the corrections for the clean water offsets (offset) and
    # the instrument's internal temperature (deltaT).
    pd = (offset - (1./0.25) * np.log(sig/ref)) - deltaT

    return pd, deltaT


def opt_tempsal_corr(channel, pd, wlngth, tcal, T, PS):
    """
    Applies the wavelength- and channel-specific temperature and
    salinity correction to an uncorrected beam attenuation or
    absorption coefficient.

    Parameters
    ----------
    channel : str
        Measurement channel: 'c' for beam attenuation or 'a' for
        absorption.
    pd : ndarray
        Uncorrected beam attenuation or absorption coefficient [m-1],
        from opt_pd_calc.
    wlngth : ndarray
        Wavelengths at which the measurements were made [nm], from
        the ac-s device file; use the 'c' or 'a' wavelengths as
        appropriate.
    tcal : float
        Factory calibration reference (pure water) temperature
        [deg_C], supplied by the manufacturer.
    T : float
        In situ temperature, TEMPWAT_L1, from a co-located CTD
        [deg_C].
    PS : float
        In situ practical salinity, PRACSAL_L2, from a co-located CTD
        [unitless].

    Returns
    -------
    pd_ts : ndarray
        Temperature- and salinity-corrected data [m-1]: OPTATTN_L2 for
        channel 'c'; an intermediate absorption product still
        requiring the scattering correction for channel 'a'.
    """
    # Absorption/attenuation and the wavelength values are imported as 1D
    # arrays. They must be the same length.
    pd = np.atleast_1d(pd)
    wlngth = np.atleast_1d(wlngth)
    lFlag = len(pd) != len(wlngth)
    if lFlag:
        raise ValueError('pd and wavelength arrays must be the same length')

    nValues = np.size(pd)

    # apply the temperature and salinity corrections for each wavelength
    # use a dictionary comprehension to read in only those values required into a np array
    np_tscor = np.array([tscor[ii] for ii in wlngth])
    dT = T - tcal
    if channel == 'a':
        pd_ts = pd - dT * np_tscor[:, 0] - PS * np_tscor[:, 2]
    elif channel == 'c':
        pd_ts = pd - dT * np_tscor[:, 0] - PS * np_tscor[:, 1]
    else:
        raise ValueError('Channel must be either "a" or "c"')

    return pd_ts


def opt_scatter_corr(apd_ts, awlngth, cpd_ts, cwlngth, rwlngth=715.):
    """
    Applies the scattering correction to the temperature- and
    salinity-corrected optical absorption coefficient, producing
    OPTABSN_L2.

    Parameters
    ----------
    apd_ts : ndarray
        Absorption coefficient corrected for temperature and
        salinity [m-1], from opt_tempsal_corr.
    awlngth : ndarray
        Absorption channel wavelengths [nm], from the ac-s device
        file.
    cpd_ts : ndarray
        Beam attenuation coefficient corrected for temperature and
        salinity, OPTATTN_L2 [m-1], from opt_tempsal_corr.
    cwlngth : ndarray
        Attenuation channel wavelengths [nm], from the ac-s device
        file.
    rwlngth : float, optional
        Scattering correction reference wavelength [nm] (default
        715).

    Returns
    -------
    apd_ts_s : ndarray
        Absorption coefficient corrected for temperature, salinity,
        and scattering, OPTABSN_L2 [m-1].

    Notes
    -----
    Numpy 1.10.1 interpolation changes caused unphysical values
    (-10^15) instead of the expected null correction when cref =
    aref, due to a divide-by-near-zero in the scatter ratio
    calculation. The fix traps the divide-by-zero more robustly and
    also accounts for baseline variability in the optical signals
    (not specified in the DPS): since attenuation c = absorption a +
    scattering b, and a, b, c are all non-negative, c is always >= a,
    so a small or noisy c-a difference should not drive an unstable
    correction. A minimum reference-wavelength scattering leakage
    value (ref_scatter_leakage_min = 0.02, chosen empirically) is used
    below which the correction is not applied; see Additional Notes on
    the docs page for further detail, including known limitations of
    the two-intake-tube deployment configuration.
    """
    ref_scatter_leakage_min = 0.02  # see Notes

    # Absorption and the absorption wavelength values are imported as 1D
    # arrays. They must be the same length.
    apd_ts = np.atleast_1d(apd_ts)
    awlngth = np.atleast_1d(awlngth)
    lFlag = len(apd_ts) != len(awlngth)
    if lFlag:
        raise ValueError('Absorption and absorption wavelength arrays must ',
                         'be the same length')

    # Attenuation and the attenuation wavelength values are imported as 1D
    # arrays. They must be the same length.
    cpd_ts = np.atleast_1d(cpd_ts)
    cwlngth = np.atleast_1d(cwlngth)
    lFlag = len(cpd_ts) != len(cwlngth)
    if lFlag:
        raise ValueError('Attenuation and attenuation wavelength arrays must ',
                         'be the same length')

    # find the the 'a' channel wavelength closest to the reference wavelength
    # for scattering and set the 'a' scattering reference value.
    idx = (np.abs(awlngth-rwlngth)).argmin()
    aref = apd_ts[idx]

    # interpolate the 'c' channel cpd_ts values to match the 'a' channel
    # wavelengths and set the 'c' scattering reference value. 
    cintrp = np.interp(awlngth, cwlngth, cpd_ts)
    cref = cintrp[idx]

    # trap out potential problems in scat_ratio calculation:
    # scat_ratio = aref / (cref - aref).
    #     aref must be > 0 AND scat_ratio must be > a pre-determined minimum;
    #     else, scat_ratio = 0.
    if aref <= 0.0:
        scat_ratio = 0.
    elif cref - aref <= ref_scatter_leakage_min:
        scat_ratio = 0.
    else:
        scat_ratio = aref / (cref - aref)

    # apply the scattering corrections
    apd_ts_s = apd_ts - scat_ratio * (cintrp - apd_ts)
    return apd_ts_s


# The next 2 functions are not used in calculating optical absorption and beam attenuation
# coefficients from the OPTAA family of instruments. However, some of these instruments
# may be outfitted with an auxiliary pressure sensor and/or external temperature sensor.
#
# opt_pressure is not used in the calculation of the final OPTAA data products.
def opt_pressure(praw, offset, sfactor):
    """
    Calculates the depth of the ac-s from its optional auxiliary
    pressure sensor. Not used in computing OPTATTN_L2 or OPTABSN_L2.

    Parameters
    ----------
    praw : ndarray
        Raw pressure reading [counts].
    offset : float
        Depth offset from the instrument device file [m].
    sfactor : float
        Scale factor from the instrument device file [m counts-1].

    Returns
    -------
    depth : ndarray
        Depth of the instrument [m].
    """
    depth = praw * sfactor + offset
    return depth


# opt_external_temp is not used in the calculation of the final OPTAA data products.
def opt_external_temp(traw):
    """
    Calculates the external environmental temperature of the ac-s from
    its optional auxiliary temperature sensor. Not used in computing
    OPTATTN_L2 or OPTABSN_L2.

    Parameters
    ----------
    traw : ndarray
        Raw external temperature [counts].

    Returns
    -------
    degC : ndarray
        Calculated external environment temperature [deg_C].
    """
    # convert counts to degrees Centigrade
    a = -7.1023317e-13
    b = 7.09341920e-08
    c = -3.87065673e-03
    d = 95.8241397

    degC = a * traw**3 + b * traw**2 + c * traw + d
    return degC


def opt_par_satlantic(counts_output, a0, a1, Im):
    """
    Computes OPTPARW_L1 from the Satlantic PAR LIN 600m instrument
    (PARAD, RSN Shallow Profiler) using a linear calibration.

    Parameters
    ----------
    counts_output : ndarray
        OPTPARW_L0 raw ADC counts.
    a0 : float
        Voltage offset [counts].
    a1 : float
        Scaling factor [umol photons m-2 s-1 count-1].
    Im : float
        Immersion coefficient.

    Returns
    -------
    OPTPARW_L1 : ndarray
        Photosynthetically active radiation [umol photons m-2 s-1].
    """

    OPTPARW_L1 = np.atleast_1d(Im * a1 * (counts_output - a0))

    return OPTPARW_L1


def opt_par_wetlabs(counts_output, a0, a1, Im):
    """
    Computes OPTPARW_L1 from the Sea-Bird Scientific ECO PAR /
    ECOPARS instrument (PARAD-J, CSPP) using an exponential
    calibration.

    Parameters
    ----------
    counts_output : ndarray
        OPTPARW_L0 raw ADC counts.
    a0 : float
        Voltage offset [counts].
    a1 : float
        Scaling factor [umol photons m-2 s-1 count-1].
    Im : float
        Immersion coefficient.

    Returns
    -------
    OPTPARW_L1 : ndarray
        Photosynthetically active radiation [umol photons m-2 s-1].
    """

    counts_output = counts_output * 1.0  # type conversion

    OPTPARW_L1 = np.atleast_1d(Im * 10**((counts_output - a0) / a1))

    return OPTPARW_L1


def opt_par_biospherical_mobile(output, dark_offset, scale_wet):
    """
    Computes OPTPARW_L1 from the Biospherical QSP-2100 scalar PAR
    instrument (PARAD, glider and other mobile assets).

    Parameters
    ----------
    output : ndarray
        OPTPARW_L0 raw sensor output [V].
    dark_offset : float
        Dark reading [V].
    scale_wet : float
        Wet calibration scale factor [V per umol photons m-2 s-1].

    Returns
    -------
    OPTPARW_L1 : ndarray
        Photosynthetically active radiation [umol photons m-2 s-1].
    """

    OPTPARW_L1 = np.atleast_1d((output - dark_offset) / scale_wet)

    return OPTPARW_L1


def opt_par_biospherical_wfp(output, dark_offset, scale_wet):
    """
    Computes OPTPARW_L1 from the Biospherical QSP-2200 scalar PAR
    instrument (PARAD, wire-following profiler).

    Parameters
    ----------
    output : ndarray
        OPTPARW_L0 raw sensor output [mV].
    dark_offset : float
        Dark reading [mV].
    scale_wet : float
        Wet calibration scale factor [V per quanta cm-2 s-1].

    Returns
    -------
    OPTPARW_L1 : ndarray
        Photosynthetically active radiation [umol photons m-2 s-1].

    Notes
    -----
    Converts the millivolt inputs to volts and the scale factor to SI
    units (1 umol photons m-2 s-1 = 6.02e13 quanta cm-2 s-1) before
    applying the same linear equation used by
    opt_par_biospherical_mobile.
    """

    #Convert output from mvolts to volts
    output_volts = output / 1000.

    #Convert dark_offset from mvolts to volts
    dark_offset_volts = dark_offset / 1000.

    #Convert scale_wet from Volts/(quanta/cm^2.s^1) to Volts/(umol photons/m^2.s^1)
    #1uE/sec/m^2 PAR= 1umole/sec/m^2 PAR = 6.02*10**13 quanta/sec/cm^2 PAR
    scale_wet_converted = scale_wet * (6.02 * 10**13)

    OPTPARW_L1 = np.atleast_1d((output_volts - dark_offset_volts) / scale_wet_converted)

    return OPTPARW_L1


def opt_ocr507_irradiance(counts, offset, scale, immersion_factor):
    """
    Computes SPECTIR_L1 from the Sea-Bird Scientific OCR-507
    multispectral radiometer (SPKIR).

    Parameters
    ----------
    counts : ndarray
        SPECTIR_L0 raw ADC counts, one column per wavelength channel
        (7 channels).
    offset : ndarray
        Dark offset calibration coefficient, per channel.
    scale : ndarray
        Scale factor calibration coefficient, per channel.
    immersion_factor : ndarray
        Immersion coefficient, per channel.

    Returns
    -------
    Ed : ndarray
        Downwelling vector irradiance, SPECTIR_L1 [uW cm-2 nm-1].
    """
    # condition input to be arrays for error-checking, in case scalars are input
    counts = np.atleast_2d(counts*1.0)  # type conversion from fix to float in case needed.
    offset = np.atleast_2d(offset)
    scale = np.atleast_2d(scale)
    immersion_factor = np.atleast_2d(immersion_factor)

    # check to see that there are 7 columns (corresponding to 7 wavelength channels) ...
    lFlag1 = counts.shape[-1] != 7
    # ... and check that the shapes of all input arguments are identical
    lFlag2 = not (counts.shape == offset.shape == scale.shape == immersion_factor.shape)
    if lFlag1 or lFlag2:
        raise ValueError('counts, offset, scale, and immersion arrays must have the same shape and have 7 columns')

    # Apply cal coeffs to raw data
    Ed = (counts - offset) * scale * immersion_factor
    return Ed
