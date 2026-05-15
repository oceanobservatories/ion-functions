#!/usr/bin/env python
"""
Module containing NIT (NUTNR) data processing functions for the Ocean
Observatories Initiative. Computes the temperature and salinity corrected
dissolved nitrate concentration (NITRTSC_L2) from raw UV absorption spectra
produced by the Sea-Bird Scientific SUNA V2 nitrate sensor.
"""

import numpy as np


def ts_corrected_nitrate(cal_temp, wl, eno3, eswa, di, dark_value, ctd_t,
                         ctd_sp, data_in, frame_type, wllower=217,
                         wlupper=240):
    """
    Compute NITRTSC_L2, dissolved nitrate concentration corrected for
    temperature and salinity.

    Parameters
    ----------
    cal_temp : array_like, shape (N,)
        Calibration water temperature [degC]. Scalar value tiled to length
        N by the OOI data management system.
    wl : array_like, shape (N, 256)
        Wavelength bins for each spectral channel [nm].
    eno3 : array_like, shape (N, 256)
        Nitrate molar absorptivity (extinction coefficients) at each
        wavelength, from factory calibration [1/(M cm)].
    eswa : array_like, shape (N, 256)
        Seawater extinction coefficients at each wavelength at reference
        salinity 35 and calibration temperature cal_temp, from factory
        calibration [1/(M cm)].
    di : array_like, shape (N, 256)
        Deionized water reference spectrum from factory calibration
        [counts].
    dark_value : array_like, shape (N,)
        Dark current scalar, averaged from dark frame measurements
        [counts].
    ctd_t : array_like, shape (N,)
        In situ water temperature from co-located CTD (TEMPWAT_L1)
        [degC].
    ctd_sp : array_like, shape (N,)
        Practical salinity from co-located CTD (PRACSAL_L2) [unitless].
    data_in : array_like, shape (N, 256)
        Raw UV absorption spectrum from the SUNA V2 (NITROPT_L0)
        [counts].
    frame_type : array_like, shape (N,)
        Frame type string for each data packet. Light frames ('SLB') are
        processed; dark frames ('SDB', 'SDF', 'NDF') are filled with NaN.
    wllower : float or array_like, shape (N,), optional
        Lower wavelength limit for spectral fitting window [nm]. Default
        is 217 nm (1-cm pathlength probe tip). Use 220 nm for 4-cm
        pathlength probe tips.
    wlupper : float or array_like, shape (N,), optional
        Upper wavelength limit for spectral fitting window [nm]. Default
        is 240 nm (1-cm pathlength probe tip). Use 245 nm for 4-cm
        pathlength probe tips.

    Returns
    -------
    NO3_conc : ndarray, shape (N,)
        Temperature and salinity corrected dissolved nitrate concentration
        (NITRTSC_L2) [uM]. Dark frame records are filled with NaN.

    Notes
    -----
    Implements the Sakamoto et al. (2009) TS-corrected nitrate algorithm.
    For each light frame, absorbance is computed from the ratio of the
    deionized water reference to the dark-corrected seawater spectrum.
    The seawater bromide contribution is removed using temperature-
    corrected seawater extinction coefficients scaled by practical
    salinity. Nitrate concentration is then obtained by linear least
    squares using a model matrix that includes the nitrate extinction
    spectrum plus a linear baseline in wavelength.

    The four fixed Sakamoto coefficients (Asak=1.1500276, Bsak=0.02840,
    Csak=-0.3101349, Dsak=0.001222) parameterize the absorbance of
    seasalt at 35 salinity vs temperature (Sakamoto et al. 2009, Eq. 4).
    """
    n_data_packets = data_in.shape[0]

    # make sure that the dimensionalities of wllower and wlupper are consistent
    # regardless of whether or not they are specified in the argument list.
    if np.isscalar(wllower):
        wllower = np.tile(wllower, n_data_packets)

    if np.isscalar(wlupper):
        wlupper = np.tile(wlupper, n_data_packets)

    # coefficients to equation 4 of Sakamoto et al 2009 that give the
    # absorbance of seasalt at 35 salinity versus temperature
    Asak = 1.1500276
    Bsak = 0.02840
    Csak = -0.3101349
    Dsak = 0.001222

    NO3_conc = np.ones(n_data_packets)

    for i in range(0, n_data_packets):

        if frame_type[i] == 'SDB' or frame_type[i] == 'SDF' or frame_type[i] == "NDF":

            # change this to output nans instead.
            NO3_conc[i] = np.nan

        else:

            # Find wavelength bins that fall between the upper and lower
            # limits for spectra fit
            useindex = np.logical_and(wllower[i] <= wl[i, :], wl[i, :] <= wlupper[i])

            # subset data so that we only use wavelengths between wllower & wlupper
            WL = wl[i, useindex]
            ENO3 = eno3[i, useindex]
            ESWA = eswa[i, useindex]
            DI = np.array(di[i, useindex], dtype='float64')
            SW = np.array(data_in[i, useindex], dtype='float64')

            # correct each SW intensity for dark current
            SWcorr = SW - dark_value[i]

            # calculate absorbance
            Absorbance = np.log10(DI / SWcorr)

            # now estimate molar absorptivity of seasalt at in situ temperature
            # use Satlantic calibration and correct as in Sakamoto et al. 2009.
            SWA_Ext_at_T = (ESWA * ((Asak + Bsak * ctd_t[i]) / (Asak + Bsak * cal_temp[i]))
                            * np.exp(Dsak * (ctd_t[i] - cal_temp[i]) * (WL - 210.0)))

            # absorbance due to seasalt
            A_SWA = ctd_sp[i] * SWA_Ext_at_T
            # subtract seasalt absorbance from measured absorbance
            Acomp = np.array(Absorbance - A_SWA, ndmin=2).T

            # ENO3 plus a linear baseline
            subset_array_size = np.shape(ENO3)
            # for the constant in the linear baseline
            Ones = np.ones((subset_array_size[0],), dtype='float64') / 100
            M = np.vstack((ENO3, Ones, WL / 1000)).T

            # C has NO3, baseline constant, and slope (vs. WL)
            C = np.dot(np.linalg.pinv(M), Acomp)

            NO3_conc[i] = C[0, 0]

    return NO3_conc
