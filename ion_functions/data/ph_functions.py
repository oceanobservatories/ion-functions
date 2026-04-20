#!/usr/bin/env python
"""
Module containing pH family instrument data processing functions for the
Sunburst SAMI-II pH instrument (PHSEN).
"""

# imports
import numpy as np
import scipy as sp


def ph_434_intensity(light):
    """
    Extract the L0 signal intensity at 434 nm (PH434SI_L0) from the PHSEN
    light measurement array.

    Parameters
    ----------
    light : array_like
        Raw light measurement array from the PHSEN instrument. May be a
        single record or an array of records. The array is reshaped
        internally to (-1, 23, 4) where the second index selects the
        measurement set and the third index selects the channel [counts].

    Returns
    -------
    si434 : ndarray, shape (nRec, 23)
        Signal intensity at 434 nm (PH434SI_L0) [counts].

    Notes
    -----
    Each PHSEN data record contains 23 sets of 4 light measurements
    interleaved as [ref434, sig434, ref578, sig578]. Column index 1
    (0-based) of each set is the 434 nm signal intensity.
    """
    light = np.atleast_3d(light).astype(float)
    new = np.reshape(light, (-1, 23, 4))
    si434 = new[:, :, 1]
    return si434  # signal intensity, 434 nm (PH434SI_L0)


def ph_578_intensity(light):
    """
    Extract the L0 signal intensity at 578 nm (PH578SI_L0) from the PHSEN
    light measurement array.

    Parameters
    ----------
    light : array_like
        Raw light measurement array from the PHSEN instrument. May be a
        single record or an array of records. The array is reshaped
        internally to (-1, 23, 4) where the second index selects the
        measurement set and the third index selects the channel [counts].

    Returns
    -------
    si578 : ndarray, shape (nRec, 23)
        Signal intensity at 578 nm (PH578SI_L0) [counts].

    Notes
    -----
    Each PHSEN data record contains 23 sets of 4 light measurements
    interleaved as [ref434, sig434, ref578, sig578]. Column index 3
    (0-based) of each set is the 578 nm signal intensity.
    """
    light = np.atleast_3d(light).astype(float)
    new = np.reshape(light, (-1, 23, 4))
    si578 = new[:, :, 3]
    return si578  # signal intensity, 578 nm (PH578SI_L0)


def ph_thermistor(traw, sami_bits=12):
    """
    Convert the PHSEN thermistor counts (ABSTHRM_L0) to temperature in
    degrees C.

    Parameters
    ----------
    traw : array_like
        Raw thermistor counts from the PHSEN instrument (ABSTHRM_L0)
        [counts].
    sami_bits : int, optional
        ADC resolution of the SAMI hardware generation. Use 12 for
        original SAMI-II hardware (full-scale 4096 counts) and 14 for
        newer hardware (full-scale 16384 counts). Default is 12.

    Returns
    -------
    therm : ndarray
        Thermistor temperature [deg_C].

    Notes
    -----
    The conversion uses a three-term polynomial in the natural log of the
    thermistor resistance. The full-scale count differs between SAMI
    hardware generations: 4096 for 12-bit hardware and 16384 for 14-bit
    hardware. The reference resistance is 17400 ohms in both cases.
    """
    traw = np.atleast_1d(traw)
    sami_bits = np.atleast_1d(sami_bits)

    if sami_bits[0] == 14:
        rt = np.log((traw / (16384.0 - traw)) * 17400.0)
    else:
        rt = np.log((traw / (4096.0 - traw)) * 17400.0)
    inv = 0.0010183 + 0.000241 * rt + 0.00000015 * rt**3
    therm = (1.0 / inv) - 273.15

    return therm


def ph_battery(braw, sami_bits=12):
    """
    Convert the PHSEN battery counts to battery voltage in Volts.

    Parameters
    ----------
    braw : array_like
        Raw battery counts from the PHSEN instrument [counts].
    sami_bits : int, optional
        ADC resolution of the SAMI hardware generation. Use 12 for
        original SAMI-II hardware (full-scale 4096 counts) and 14 for
        newer hardware (full-scale 16384 counts). Default is 12.

    Returns
    -------
    volts : ndarray
        Battery voltage [Volts].

    Notes
    -----
    The full-scale voltage and count differ between SAMI hardware
    generations. For 12-bit hardware the full-scale is 15 V at 4096
    counts; for 14-bit hardware the full-scale is 3 V at 4000 counts.
    """
    braw = np.atleast_1d(braw)
    sami_bits = np.atleast_1d(sami_bits)

    if sami_bits[0] == 14:
        volts = braw * 3. / 4000.
    else:
        volts = braw * 15. / 4096.
    return volts


def ph_calc_phwater(ref, light, therm, ea434, eb434, ea578, eb578,
                    ind_slp, ind_off, psal=35.0):
    """
    Compute the OOI L2 pH of seawater (PHWATER_L2) from the Sunburst
    SAMI-II pH instrument (PHSEN).

    Parameters
    ----------
    ref : array_like, shape (nRec, 16)
        Raw blank reference and signal measurements from the PHSEN blank
        cycle. Contains 4 sets of 4 interleaved measurements:
        [ref434, sig434, ref578, sig578] [counts].
    light : array_like, shape (nRec, 92)
        Raw reference and signal measurements from the PHSEN measurement
        cycle. Contains 23 sets of 4 interleaved measurements:
        [ref434, sig434, ref578, sig578] [counts].
    therm : array_like, shape (nRec,)
        Thermistor temperature at the end of the measurement cycle,
        converted to degrees C via ph_thermistor (ABSTHRM_L0) [deg_C].
    ea434 : array_like, shape (nRec,)
        Calibration coefficient 1. Molar absorptivity of the acidic
        indicator form at 434 nm at the reference temperature [unitless].
    eb434 : array_like, shape (nRec,)
        Calibration coefficient 2. Molar absorptivity of the basic
        indicator form at 434 nm at the reference temperature [unitless].
    ea578 : array_like, shape (nRec,)
        Calibration coefficient 3. Molar absorptivity of the acidic
        indicator form at 578 nm at the reference temperature [unitless].
    eb578 : array_like, shape (nRec,)
        Calibration coefficient 4. Molar absorptivity of the basic
        indicator form at 578 nm at the reference temperature [unitless].
    ind_slp : float or array_like, shape (nRec,)
        Indicator impurity slope correction factor applied to pH values
        >= 8.2 [unitless].
    ind_off : float or array_like, shape (nRec,)
        Indicator impurity offset correction factor applied to pH values
        >= 8.2 [unitless].
    psal : float or array_like, shape (nRec,), optional
        Practical salinity from a co-located CTD. Default is 35.0 if
        CTD data are unavailable [unitless].

    Returns
    -------
    ph : ndarray, shape (nRec,)
        pH of seawater on the total hydrogen ion scale (PHWATER_L2)
        [unitless].

    Notes
    -----
    The algorithm selects the 8 most linearly consistent measurement
    points from 23 collected during each cycle (skipping the first 5)
    by finding the window of 8 consecutive points with the highest
    linear correlation coefficient (R^2) between indicator concentration
    and point pH. The final pH is extrapolated to zero indicator
    concentration from this best-fit region.

    An impurity correction (ind_slp, ind_off) is applied when the
    calculated pH exceeds 8.2. This correction is not described in DPS
    1341-00510 and was added post-publication.
    """
    ref = (np.atleast_2d(ref)).astype(float)
    nRec = ref.shape[0]

    light = np.atleast_3d(light).astype(float)
    light = np.reshape(light, (nRec, 23, 4))

    therm = np.reshape(therm, (nRec, 1)).astype(float)

    ea434 = np.reshape(ea434, (nRec, 1)).astype(float)
    eb434 = np.reshape(eb434, (nRec, 1)).astype(float)
    ea578 = np.reshape(ea578, (nRec, 1)).astype(float)
    eb578 = np.reshape(eb578, (nRec, 1)).astype(float)

    if np.isscalar(ind_slp) is True:
        ind_slp = np.tile(ind_slp, (nRec)).astype(float)
    else:
        ind_slp = np.reshape(ind_slp, (nRec)).astype(float)

    if np.isscalar(ind_off) is True:
        ind_off = np.tile(ind_off, (nRec)).astype(float)
    else:
        ind_off = np.reshape(ind_off, (nRec)).astype(float)

    if np.isscalar(psal) is True:
        psal = np.tile(psal, (nRec, 1)).astype(float)
    else:
        psal = np.reshape(psal, (nRec, 1)).astype(float)

    # Calculate blanks from the 16 sets of reference light measurements
    arr434 = np.array([
        (ref[:, 1] / ref[:, 0]),
        (ref[:, 5] / ref[:, 4]),
        (ref[:, 9] / ref[:, 8]),
        (ref[:, 13] / ref[:, 12]),
    ])
    blank434 = np.reshape(np.mean(arr434, axis=0), (nRec, 1))

    arr578 = np.array([
        (ref[:, 3] / ref[:, 2]),
        (ref[:, 7] / ref[:, 6]),
        (ref[:, 11] / ref[:, 10]),
        (ref[:, 15] / ref[:, 14]),
    ])
    blank578 = np.reshape(np.mean(arr578, axis=0), (nRec, 1))

    # Extract 23 sets of 4 light measurements into arrays corresponding
    # to the raw reference and signal measurements at 434 and 578 nm.
    ref434 = light[:, :, 0]   # reference signal, 434 nm
    int434 = light[:, :, 1]   # signal intensity, 434 nm (PH434SI_L0)
    ref578 = light[:, :, 2]   # reference signal, 578 nm
    int578 = light[:, :, 3]   # signal intensity, 578 nm (PH578SI_L0)

    # Absorbance
    A434 = -sp.log10(int434 / ref434)
    A434blank = -sp.log10(blank434)
    abs434 = A434 - A434blank

    A578 = -sp.log10(int578 / ref578)
    A578blank = -sp.log10(blank578)
    abs578 = A578 - A578blank

    R = abs578 / abs434

    # pKa from Clayton and Byrne, 1993
    pKa = (1245.69 / (therm + 273.15)) + 3.8275 + (0.0021 * (35. - psal))
    pKa = np.reshape(pKa, (-1, 1))

    # Molar absorptivities
    Ea434 = ea434 - (26. * (therm - 24.788))
    Ea578 = ea578 + (therm - 24.788)
    Eb434 = eb434 + (12. * (therm - 24.788))
    Eb578 = eb578 - (71. * (therm - 24.788))
    e1 = Ea578 / Ea434
    e2 = Eb578 / Ea434
    e3 = Eb434 / Ea434

    V1 = R - e1
    V2 = e2 - R * e3

    # indicator concentration calculations
    HI = (abs434 * Eb578 - abs578 * Eb434) / (Ea434 * Eb578 - Eb434 * Ea578)
    I = (abs578 * Ea434 - abs434 * Ea578) / (Ea434 * Eb578 - Eb434 * Ea578)
    IndConc = HI + I
    pointph = np.real(pKa + sp.log10(V1 / V2))

    # determine the most linear region of points for pH of seawater
    # calculation, skipping the first 5 points.
    IndConca = IndConc[:, 5:]
    Y = pointph[:, 5:]
    X = np.linspace(1, 18, 18)

    # create arrays for vectorized computations used in sum of squares.
    step = 7  # number of points to use
    count = step + 1
    nPts = np.size(X) - step
    x = np.zeros((nPts, count))
    y = np.zeros((nRec, nPts, count))
    for i in range(nPts):
        x[i, :] = X[i:i+count]
        for j in range(nRec):
            y[j, i, :] = Y[j, i:i+count]

    # compute correlation coefficient for each window of 8 points
    sumx = np.sum(x, axis=1)
    sumy = np.sum(y, axis=2)
    sumxy = np.sum(x * y, axis=2)
    sumx2 = np.sum(x**2, axis=1)
    sumy2 = np.sum(y**2, axis=2)
    sumxx = sumx * sumx
    sumyy = sumy * sumy
    ssxy = sumxy - (sumx * sumy) / count
    ssx = sumx2 - (sumxx / count)
    ssy = sumy2 - (sumyy / count)
    r2 = ssxy**2 / (ssx * ssy)

    # Range of seawater points to use
    cutoff1 = np.argmax(r2, axis=1)  # Find the first, best R-squared value
    cutoff2 = cutoff1 + count

    # Indicator and pH range limited to best points
    IndConcS = np.zeros((nRec, count))
    pointphS = np.zeros((nRec, count))
    for i in range(nRec):
        IndConcS[i, :] = IndConca[i, cutoff1[i]:cutoff2[i]]
        pointphS[i, :] = Y[i, cutoff1[i]:cutoff2[i]]

    # Final pH calculation: extrapolate to zero indicator concentration
    sumx = np.sum(IndConcS, axis=1)
    sumy = np.sum(pointphS, axis=1)
    sumxy = np.sum(pointphS * IndConcS, axis=1)
    sumx2 = np.sum(IndConcS**2, axis=1)
    sumy2 = np.sum(pointphS**2, axis=1)
    xbar = np.mean(IndConcS, axis=1)
    ybar = np.mean(pointphS, axis=1)
    sumxx = sumx * sumx
    sumyy = sumy * sumy
    ssxy = sumxy - (sumx * sumy) / count
    ssx = sumx2 - (sumxx / count)
    ssy = sumy2 - (sumyy / count)
    slope = ssxy / ssx
    ph = ybar - slope * xbar

    # pH corrections due to indicator impurity if the calculated pH is
    # greater than 8.2.
    phFlag = ph >= 8.2
    ph[phFlag] = ph[phFlag] * ind_slp[phFlag] + ind_off[phFlag]

    return ph
