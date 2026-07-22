#!/usr/bin/env python
"""
Module containing Seafloor Properties data processing functions for the
Ocean Observatories Initiative. Covers the THSPH, TRHPH, PRESF, and PREST
instrument families.
"""

import numpy as np
from scipy.interpolate import RectBivariateSpline

from ion_functions.data.generic_functions import replace_fill_with_nan
# used by def sfl_trhph_chloride
from ion_functions.data.sfl_functions_surface import tdat, sdat, cdat

# .............................................................................
# THSPH data products: THSPHHC, THSPHHS, THSPHPH (4 PH products) ..............
# .............................................................................


def sfl_thsph_ph(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl,
                 arr_hgo, arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl):
    """
    Compute vent fluid pH (THSPHPH-PH_L2) using measured AgCl reference
    electrode and chloride from TRHPHCC_L2.

    Parameters
    ----------
    counts_ysz : array_like
        Raw counts from YSZ electrode (THSPHPH-YSZ_L0) [counts].
    counts_agcl : array_like
        Raw counts from AgCl reference electrode (THSPHPH-AGCL_L0) [counts].
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].
    e2l_ysz : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting YSZ electrode engineering values to lab
        calibrated values.
    e2l_agcl : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting AgCl electrode engineering values to lab
        calibrated values.
    arr_hgo : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for electrode material response to temperature.
    arr_agcl : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for AgCl electrode material response to temperature.
    arr_tac : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tac (=tbc0).
    arr_tbc1 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc1.
    arr_tbc2 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc2.
    arr_tbc3 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc3.
    chl : array_like
        Vent fluid chloride concentration from TRHPHCC_L2 [mmol kg^-1].

    Returns
    -------
    pH : ndarray
        Vent fluid pH (THSPHPH-PH_L2) [dimensionless].

    Notes
    -----
    Uses measured AgCl reference electrode and chloride from TRHPHCC_L2.
    Values with electrode potential outside [-0.7, 0.0] V or pH outside
    [3.0, 7.0] are set to NaN. See calculate_vent_pH for the shared core
    algorithm.
    """
    # calculate lab calibrated electrode response [V]
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)
    # AgCl reference electrode
    v_labcal_agcl = v_labcal(counts_agcl, e2l_agcl)

    # calculate chloride activity
    act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl)

    pH = calculate_vent_pH(v_labcal_ysz, v_labcal_agcl, temperature, arr_hgo, arr_agcl, act_chl)

    return pH


def sfl_thsph_ph_acl(counts_ysz, counts_agcl, temperature, e2l_ysz, e2l_agcl,
                     arr_hgo, arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3):
    """
    Compute vent fluid pH (THSPHPH-PH-ACL_L2) using measured AgCl reference
    electrode and assumed chloride concentration.

    Parameters
    ----------
    counts_ysz : array_like
        Raw counts from YSZ electrode (THSPHPH-YSZ_L0) [counts].
    counts_agcl : array_like
        Raw counts from AgCl reference electrode (THSPHPH-AGCL_L0) [counts].
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].
    e2l_ysz : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting YSZ electrode engineering values to lab
        calibrated values.
    e2l_agcl : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting AgCl electrode engineering values to lab
        calibrated values.
    arr_hgo : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for electrode material response to temperature.
    arr_agcl : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for AgCl electrode material response to temperature.
    arr_tac : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tac (=tbc0).
    arr_tbc1 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc1.
    arr_tbc2 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc2.
    arr_tbc3 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc3.

    Returns
    -------
    pH : ndarray
        Vent fluid pH (THSPHPH-PH-ACL_L2) [dimensionless].

    Notes
    -----
    Uses measured AgCl reference electrode and a default chloride
    concentration of 250.0 mmol/kg (set in chloride_activity). Values with
    electrode potential outside [-0.7, 0.0] V or pH outside [3.0, 7.0] are
    set to NaN. See calculate_vent_pH for the shared core algorithm.
    """
    # calculate lab calibrated electrode response [V]
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)
    # AgCl reference electrode
    v_labcal_agcl = v_labcal(counts_agcl, e2l_agcl)

    # chloride activity assuming the default value for chloride concentration
    # set in the chloride_activity subroutine
    act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)

    pH = calculate_vent_pH(v_labcal_ysz, v_labcal_agcl, temperature, arr_hgo, arr_agcl, act_chl)

    return pH


def sfl_thsph_ph_noref(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                       arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl):
    """
    Compute vent fluid pH (THSPHPH-PH-NOREF_L2) using a theoretical reference
    electrode potential and chloride from TRHPHCC_L2.

    Parameters
    ----------
    counts_ysz : array_like
        Raw counts from YSZ electrode (THSPHPH-YSZ_L0) [counts].
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].
    arr_agclref : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate the theoretical reference electrode
        potential from temperature.
    e2l_ysz : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting YSZ electrode engineering values to lab
        calibrated values.
    arr_hgo : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for electrode material response to temperature.
    arr_agcl : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for AgCl electrode material response to temperature.
    arr_tac : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tac (=tbc0).
    arr_tbc1 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc1.
    arr_tbc2 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc2.
    arr_tbc3 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc3.
    chl : array_like
        Vent fluid chloride concentration from TRHPHCC_L2 [mmol kg^-1].

    Returns
    -------
    pH : ndarray
        Vent fluid pH (THSPHPH-PH-NOREF_L2) [dimensionless].

    Notes
    -----
    Uses a theoretical reference electrode potential computed from vent
    temperature (arr_agclref) in place of a measured AgCl electrode signal,
    and chloride from TRHPHCC_L2. Values with electrode potential outside
    [-0.7, 0.0] V or pH outside [3.0, 7.0] are set to NaN. See
    calculate_vent_pH for the shared core algorithm.
    """
    # calculate lab calibrated electrode response [V]
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)

    # theoretical reference value calculated from vent temperature
    e_refcalc = eval_poly(temperature, arr_agclref)
    # calculate chloride activity
    act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3, chl)

    pH = calculate_vent_pH(v_labcal_ysz, e_refcalc, temperature, arr_hgo, arr_agcl, act_chl)

    return pH


def sfl_thsph_ph_noref_acl(counts_ysz, temperature, arr_agclref, e2l_ysz, arr_hgo,
                           arr_agcl, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3):
    """
    Compute vent fluid pH (THSPHPH-PH-NOREF-ACL_L2) using a theoretical
    reference electrode potential and assumed chloride concentration.

    Parameters
    ----------
    counts_ysz : array_like
        Raw counts from YSZ electrode (THSPHPH-YSZ_L0) [counts].
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].
    arr_agclref : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate the theoretical reference electrode
        potential from temperature.
    e2l_ysz : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting YSZ electrode engineering values to lab
        calibrated values.
    arr_hgo : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for electrode material response to temperature.
    arr_agcl : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for AgCl electrode material response to temperature.
    arr_tac : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tac (=tbc0).
    arr_tbc1 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc1.
    arr_tbc2 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc2.
    arr_tbc3 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc3.

    Returns
    -------
    pH : ndarray
        Vent fluid pH (THSPHPH-PH-NOREF-ACL_L2) [dimensionless].

    Notes
    -----
    Uses a theoretical reference electrode potential computed from vent
    temperature (arr_agclref) in place of a measured AgCl electrode signal,
    and a default chloride concentration of 250.0 mmol/kg (set in
    chloride_activity). Values with electrode potential outside [-0.7, 0.0] V
    or pH outside [3.0, 7.0] are set to NaN. See calculate_vent_pH for the
    shared core algorithm.
    """
    # calculate lab calibrated electrode response [V]
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)

    # theoretical reference value calculated from vent temperature
    e_refcalc = eval_poly(temperature, arr_agclref)
    # chloride activity assuming the default value for chloride concentration
    # set in the subroutine
    act_chl = chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3)

    pH = calculate_vent_pH(v_labcal_ysz, e_refcalc, temperature, arr_hgo, arr_agcl, act_chl)

    return pH


def calculate_vent_pH(e_ph, e_ref, temperature, arr_hgo, arr_agcl, act_chl):
    """
    Compute vent fluid pH for the THSPH instrument.

    Called by sfl_thsph_ph, sfl_thsph_ph_acl, sfl_thsph_ph_noref, and
    sfl_thsph_ph_noref_acl.

    Parameters
    ----------
    e_ph : array_like
        Lab-calibrated YSZ electrode potential [V].
    e_ref : array_like
        Reference electrode potential, either measured (AgCl) or
        theoretical [V].
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].
    arr_hgo : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for electrode material response to temperature.
    arr_agcl : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for AgCl electrode material response to temperature.
    act_chl : array_like
        Chloride activity computed by chloride_activity [dimensionless].

    Returns
    -------
    pH : ndarray
        Vent fluid pH [dimensionless]. Values with electrode potential
        outside [-0.7, 0.0] V or pH outside [3.0, 7.0] are set to NaN.

    Notes
    -----
    Out-of-range checks on both the electrode potential difference and the
    final pH are applied as specified in the unreleased DPS (1341-00190).
    """
    # fill value local to this function to avoid python warnings when nans are encountered
    # in boolean expressions. the masking will convert values derived from this local fill
    # back to nans.
    unphysical_pH_fill_value = -99999.0

    # calculate intermediate quantities that depend upon temperature
    e_nernst = nernst(temperature)
    e_hgo = eval_poly(temperature, arr_hgo)
    e_agcl = eval_poly(temperature, arr_agcl)

    # calculate pH potential
    e_phcalc = e_ph - e_ref
    # check for unphysical values as specified in the DPS.
    # logical indexing with boolean arrays is faster than integer indexing using np.where.
    # ok to apply mask at end of calculation.
    e_phcalc[np.isnan(e_phcalc)] = unphysical_pH_fill_value
    bad_eph_mask = np.logical_or(np.less(e_phcalc, -0.7), np.greater(e_phcalc, 0.0))

    # final data product calculation
    act_chl[act_chl <= 0.0] = np.nan  # trap out python warning
    pH = (e_phcalc - e_agcl + e_hgo) / e_nernst + np.log10(act_chl)

    # second check for unphysical values, as specified in the DPS
    pH[np.isnan(pH)] = unphysical_pH_fill_value
    bad_ph_mask = np.logical_or(np.less(pH, 3.0), np.greater(pH, 7.0))

    # set all out-of-range values to fill values
    pH[np.logical_or(bad_eph_mask, bad_ph_mask)] = np.nan

    return pH


def sfl_thsph_sulfide(counts_hs, counts_ysz, temperature, e2l_hs, e2l_ysz, arr_hgo,
                      arr_logkfh2g, arr_eh2sg, arr_yh2sg):
    """
    Compute vent fluid hydrogen sulfide concentration (THSPHHS_L2).

    Parameters
    ----------
    counts_hs : array_like
        Raw counts from sulfide electrode (THSPHHS_L0) [counts].
    counts_ysz : array_like
        Raw counts from YSZ electrode (THSPHPH-YSZ_L0) [counts].
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].
    e2l_hs : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting sulfide electrode engineering values to lab
        calibrated values.
    e2l_ysz : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting YSZ electrode engineering values to lab
        calibrated values.
    arr_hgo : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for electrode material response to temperature.
    arr_logkfh2g : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for equilibrium hydrogen fugacity as a function of
        temperature.
    arr_eh2sg : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for theoretical potential of gas phase H2S; pad unused
        high-degree terms with zeros: [0., 0., 0., 0., c1, c0].
    arr_yh2sg : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for the fugacity/concentration quotient yh2sg.

    Returns
    -------
    h2s : ndarray
        Hydrogen sulfide concentration at the vent (THSPHHS_L2) [mmol kg^-1].

    Notes
    -----
    The DPS document for THSPHHS (1341-00200) was never publicly released.
    The algorithm is documented from the code and code comments only.
    """
    # calculate lab calibrated electrode responses [V]
    v_labcal_hs = v_labcal(counts_hs, e2l_hs)
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)

    # calculate intermediate products that depend upon temperature
    e_nernst = nernst(temperature)
    e_hgo = eval_poly(temperature, arr_hgo)
    e_h2sg = eval_poly(temperature, arr_eh2sg)
    log_kfh2g = eval_poly(temperature, arr_logkfh2g)
    # y_h2sg depends on temperature because hydrogen fugacity depends on temperature
    y_h2sg = eval_poly(log_kfh2g, arr_yh2sg)

    # explicitly follow the DPS calculation for clarity:

    # measured potential of the sulfide electrode [V]
    e_h2s = v_labcal_ysz - v_labcal_hs

    # (common) log of measured hydrogen sulfide fugacity
    log_fh2sg = 2.0 * (e_h2s - e_hgo + e_h2sg) / e_nernst

    # final data product, hydrogen sulfide concentration, [mmol/kg]
    # in the DPS, this is 1000 * 10^( logfh2sg - log( yh2sg ) )
    h2s = 1000.0 * (10.0 ** (log_fh2sg)) / y_h2sg

    return h2s


def sfl_thsph_hydrogen(counts_h2, counts_ysz, temperature, e2l_h2, e2l_ysz, arr_hgo,
                       arr_logkfh2g):
    """
    Compute vent fluid hydrogen concentration (THSPHHC_L2).

    Parameters
    ----------
    counts_h2 : array_like
        Raw counts from hydrogen electrode (THSPHHC_L0) [counts].
    counts_ysz : array_like
        Raw counts from YSZ electrode (THSPHPH-YSZ_L0) [counts].
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].
    e2l_h2 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting hydrogen electrode engineering values to lab
        calibrated values.
    e2l_ysz : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients converting YSZ electrode engineering values to lab
        calibrated values.
    arr_hgo : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for electrode material response to temperature.
    arr_logkfh2g : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients for equilibrium hydrogen fugacity as a function of
        temperature.

    Returns
    -------
    h2 : ndarray
        Hydrogen concentration at the vent (THSPHHC_L2) [mmol kg^-1].

    Notes
    -----
    The DPS document for THSPHHC (1341-00210) was never publicly released.
    The algorithm is documented from the code and code comments only.
    """
    # calculate lab calibrated electrode responses [V]
    v_labcal_h2 = v_labcal(counts_h2, e2l_h2)
    v_labcal_ysz = v_labcal(counts_ysz, e2l_ysz)

    # calculate intermediate products that depend upon temperature
    e_nernst = nernst(temperature)
    e_hgo = eval_poly(temperature, arr_hgo)
    log_kfh2g = eval_poly(temperature, arr_logkfh2g)

    # explicitly follow the DPS calculation for clarity:

    # measured potential of the h2 electrode [V]
    e_h2 = v_labcal_ysz - v_labcal_h2

    # (common) log of measured hydrogen fugacity
    log_fh2g = 2.0 * (e_h2 - e_hgo) / e_nernst

    # final data product, hydrogen concentration, [mmol/kg]
    h2 = 1000.0 * (10.0 ** (log_fh2g - log_kfh2g))

    return h2


def chloride_activity(temperature, arr_tac, arr_tbc1, arr_tbc2, arr_tbc3,
                      chloride=250.0):
    """
    Compute chloride activity for THSPHPH_L2 data products.

    Parameters
    ----------
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].
    arr_tac : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tac (=tbc0).
    arr_tbc1 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc1.
    arr_tbc2 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc2.
    arr_tbc3 : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients to calculate tbc3.
    chloride : array_like, optional
        Vent fluid chloride concentration from TRHPHCC_L2 [mmol kg^-1].
        Defaults to 250.0 mmol/kg when not supplied.

    Returns
    -------
    act_chl : ndarray
        Chloride activity [dimensionless].
    """
    # find number of data packets to be processed;
    # this also works if temperature is not an np.array.
    nvalues = np.array([temperature]).shape[-1]

    # change units of chloride from mmol/kg to mol/kg
    chloride = chloride/1000.0

    # if chloride is not given in the argument list,
    # replicate its default value into a vector with
    # the same number of elements as temperature;
    # do so without using a conditional
    nreps = nvalues // np.array([chloride]).shape[-1]
    chloride = np.tile(chloride, nreps)

    # calculate the 4 coefficients needed to calculate the chloride activity from temperature
    tbc0 = eval_poly(temperature, arr_tac)
    tbc1 = eval_poly(temperature, arr_tbc1)
    tbc2 = eval_poly(temperature, arr_tbc2)
    tbc3 = eval_poly(temperature, arr_tbc3)

    # form these coeffs into a 2D array for the eval_poly routine.
    # need to pad the first two columns with zeros
    zeros = np.array([np.tile(0.0, nvalues)]).T
    arr_chloride_coeff = np.hstack((zeros, zeros, tbc3[:, np.newaxis], tbc2[:, np.newaxis],
                                    tbc1[:, np.newaxis], tbc0[:, np.newaxis]))

    # evaluate the activity
    act_chl = eval_poly(chloride, arr_chloride_coeff)

    return act_chl


def v_labcal(counts, array_e2l_coeff):
    """
    Convert raw THSPH electrode counts to lab-calibrated voltage.

    Used by all THSPH L2 data products (THSPHHC, THSPHHS, THSPHPH).

    Parameters
    ----------
    counts : array_like
        L0 decimal counts from one of the four THSPH electrodes
        (THSPHPH-YSZ_L0, THSPHPH-AGC_L0, THSPHHC_L0, THSPHHS_L0)
        [counts].
    array_e2l_coeff : array_like
        Calibration coefficient. 6-element array of 5th-degree polynomial
        coefficients (descending order) converting engineering values to lab
        calibrated values.

    Returns
    -------
    v_labcal_electrode : ndarray
        Lab-calibrated electrode voltage (V_actual) [V].

    Notes
    -----
    System fill values in counts are replaced with NaN before conversion,
    so fill value replacement applies to all THSPH L2 data products that
    call this function.
    """
    counts = replace_fill_with_nan(None, counts)

    # transform decimal counts to engineering values [volts]
    v_eng = (counts * 0.25 - 2048.0) / 1000.0

    # transform engineering values to lab calibrated values [volts]
    v_labcal_electrode = eval_poly(v_eng, array_e2l_coeff)

    # in the DPSs, these values are designated as "V_actual"
    return v_labcal_electrode


def nernst(temperature):
    """
    Compute the temperature-dependent term of the Nernst equation.

    Used by all THSPH L2 data products (THSPHHC, THSPHHS, THSPHPH).

    Parameters
    ----------
    temperature : array_like
        Temperature near sample inlet (THSPHTE-TH_L1) [deg_C].

    Returns
    -------
    e_nernst : ndarray
        Temperature-dependent Nernst factor [V].
    """
    # e_nernst = ln(10) * (gas constant) * (temperature, Kelvin)/(Faraday's constant)
    #          = 2.30259 * 8.31446 [J/mole/K] / 96485.3 [coulombs/mole] * (T + 273.15)
    return 1.9842e-4 * (temperature + 273.15)


# .............................................................................
# THSPH data products: THSPHTE -TH, -TL, -TCH, -TCL, -REF, -INT ...............
# .............................................................................


def sfl_thsph_temp_th(tc_rawdec_H, e2l_H, l2s_H, ts_rawdec_r, e2l_r, l2s_r, s2v_r):
    """
    Compute vent fluid temperature at the sample inlet (THSPHTE-TH_L1).

    Parameters
    ----------
    tc_rawdec_H : array_like
        H thermocouple decimal counts (THSPHTE-TCH_L0) [counts].
    e2l_H : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        H thermocouple engineering values to lab calibrated values.
    l2s_H : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        H thermocouple lab calibrated values to scientific values [deg_C].
    ts_rawdec_r : array_like
        Reference thermistor decimal counts (THSPHTE-REF_L0) [counts].
    e2l_r : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        reference thermistor engineering values to lab calibrated values.
    l2s_r : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        reference thermistor lab calibrated values to scientific values
        [deg_C].
    s2v_r : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        reference thermistor scientific values to thermocouple equivalent
        voltage [mV].

    Returns
    -------
    T_H : ndarray
        Final temperature at position H near sample inlet
        (THSPHTE-TH_L1) [deg_C].
    """
    # calculate intermediate product V_tc_actual_H (= V_tc_labcal_H)
    V_tc_actual_H = sfl_thsph_temp_labcal_h(tc_rawdec_H, e2l_H)

    # calculate intermediate products T_ts_r, then V_ts_r (June 2014 DPS)
    T_ts_r = sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r)
    V_ts_r = eval_poly(T_ts_r, s2v_r)

    # Correct thermocouple temperature to account for offset from cold junction as
    # measured by the reference thermistor
    T_H = eval_poly((V_tc_actual_H + V_ts_r), l2s_H)

    return T_H


def sfl_thsph_temp_tl(tc_rawdec_L, e2l_L, l2s_L, ts_rawdec_r, e2l_r, l2s_r, s2v_r):
    """
    Compute vent fluid temperature near the vent (THSPHTE-TL_L1).

    Parameters
    ----------
    tc_rawdec_L : array_like
        L thermocouple decimal counts (THSPHTE-TCL_L0) [counts].
    e2l_L : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        L thermocouple engineering values to lab calibrated values.
    l2s_L : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        L thermocouple lab calibrated values to scientific values [deg_C].
    ts_rawdec_r : array_like
        Reference thermistor decimal counts (THSPHTE-REF_L0) [counts].
    e2l_r : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        reference thermistor engineering values to lab calibrated values.
    l2s_r : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        reference thermistor lab calibrated values to scientific values
        [deg_C].
    s2v_r : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        reference thermistor scientific values to thermocouple equivalent
        voltage [mV].

    Returns
    -------
    T_L : ndarray
        Final temperature at position L near vent (THSPHTE-TL_L1) [deg_C].
    """
    # calculate intermediate product V_tc_actual_L (= V_tc_labcal_L)
    V_tc_actual_L = sfl_thsph_temp_labcal_l(tc_rawdec_L, e2l_L)

    # calculate intermediate products T_ts_r, then V_ts_r (June 2014 DPS)
    T_ts_r = sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r)
    V_ts_r = eval_poly(T_ts_r, s2v_r)

    # Correct thermocouple temperature to account for offset from cold junction as
    # measured by the reference thermistor
    T_L = eval_poly((V_tc_actual_L + V_ts_r), l2s_L)

    return T_L


def sfl_thsph_temp_tch(tc_rawdec_H, e2l_H, l2s_H):
    """
    Compute intermediate thermocouple temperature at position H
    (THSPHTE-TCH_L1).

    Parameters
    ----------
    tc_rawdec_H : array_like
        H thermocouple decimal counts (THSPHTE-TCH_L0) [counts].
    e2l_H : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        H thermocouple engineering values to lab calibrated values.
    l2s_H : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        H thermocouple lab calibrated values to scientific values [deg_C].

    Returns
    -------
    T_tc_H : ndarray
        Intermediate thermocouple temperature at position H
        (THSPHTE-TCH_L1) [deg_C].
    """
    # convert raw decimal output to lab calibrated values [mV]
    V_tc_actual_H = sfl_thsph_temp_labcal_h(tc_rawdec_H, e2l_H)

    # convert lab calibrated values to scientific values [degC]
    T_tc_H = eval_poly(V_tc_actual_H, l2s_H)

    return T_tc_H


def sfl_thsph_temp_tcl(tc_rawdec_L, e2l_L, l2s_L):
    """
    Compute intermediate thermocouple temperature at position L
    (THSPHTE-TCL_L1).

    Parameters
    ----------
    tc_rawdec_L : array_like
        L thermocouple decimal counts (THSPHTE-TCL_L0) [counts].
    e2l_L : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        L thermocouple engineering values to lab calibrated values.
    l2s_L : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        L thermocouple lab calibrated values to scientific values [deg_C].

    Returns
    -------
    T_tc_L : ndarray
        Intermediate thermocouple temperature at position L
        (THSPHTE-TCL_L1) [deg_C].
    """
    # convert raw decimal output to lab calibrated values [mV]
    V_tc_actual_L = sfl_thsph_temp_labcal_l(tc_rawdec_L, e2l_L)

    # convert lab calibrated values to scientific values [degC]
    T_tc_L = eval_poly(V_tc_actual_L, l2s_L)

    return T_tc_L


def sfl_thsph_temp_ref(ts_rawdec_r, e2l_r, l2s_r):
    """
    Compute reference thermistor temperature (THSPHTE-REF_L1).

    Parameters
    ----------
    ts_rawdec_r : array_like
        Reference thermistor decimal counts (THSPHTE-REF_L0) [counts].
    e2l_r : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        reference thermistor engineering values to lab calibrated values.
    l2s_r : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        reference thermistor lab calibrated values to scientific values
        [deg_C].

    Returns
    -------
    T_ts_r : ndarray
        Reference thermistor temperature (THSPHTE-REF_L1) [deg_C].
    """
    ts_rawdec_r = replace_fill_with_nan(None, ts_rawdec_r)

    # convert raw decimal output to engineering values [ohms]
    ts_rawdec_r_scaled = ts_rawdec_r * 0.125
    denom = 2048.0 - ts_rawdec_r_scaled
    denom[denom == 0.0] = np.nan
    R_ts_eng_r = 10000.0 * ts_rawdec_r_scaled / denom

    # convert engineering values to lab calibrated values [ohms]
    R_ts_actual_r = eval_poly(R_ts_eng_r, e2l_r)

    # convert lab calibrated values to scientific values [degC]
    R_ts_actual_r[R_ts_actual_r <= 0.0] = np.nan
    pval = eval_poly(np.log(R_ts_actual_r), l2s_r)
    T_ts_r = 1.0 / pval - 273.15

    return T_ts_r


def sfl_thsph_temp_int(ts_rawdec_b, e2l_b, l2s_b):
    """
    Compute internal board thermistor temperature (THSPHTE-INT_L1).

    Parameters
    ----------
    ts_rawdec_b : array_like
        Board thermistor decimal counts (THSPHTE-INT_L0) [counts].
    e2l_b : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        board thermistor engineering values to lab calibrated values.
    l2s_b : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        board thermistor lab calibrated values to scientific values [deg_C].

    Returns
    -------
    T_ts_b : ndarray
        Board thermistor temperature (THSPHTE-INT_L1) [deg_C].
    """
    ts_rawdec_b = replace_fill_with_nan(None, ts_rawdec_b)

    # convert raw decimal output to engineering values [ohms]
    ts_rawdec_b_scaled = ts_rawdec_b * 0.125
    denom = 2048.0 - ts_rawdec_b_scaled
    denom[denom == 0.0] = np.nan
    R_ts_eng_b = 10000.0 * ts_rawdec_b_scaled / denom

    # convert engineering values to lab calibrated values [ohms]
    R_ts_actual_b = eval_poly(R_ts_eng_b, e2l_b)

    # convert lab calibrated values to scientific values [degC]
    R_ts_actual_b[R_ts_actual_b <= 0.0] = np.nan
    pval = eval_poly(np.log(R_ts_actual_b), l2s_b)

    T_ts_b = 1.0 / pval - 273.15

    return T_ts_b


def sfl_thsph_temp_labcal_h(tc_rawdec_H, e2l_H):
    """
    Convert H thermocouple raw counts to lab-calibrated voltage.

    Called internally by sfl_thsph_temp_th and sfl_thsph_temp_tch.

    Parameters
    ----------
    tc_rawdec_H : array_like
        H thermocouple decimal counts (THSPHTE-TCH_L0) [counts].
    e2l_H : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        H thermocouple engineering values to lab calibrated values.

    Returns
    -------
    V_tc_labcal_H : ndarray
        Lab-calibrated H thermocouple voltage [mV].
    """
    tc_rawdec_H = replace_fill_with_nan(None, tc_rawdec_H)

    # convert raw decimal output to engineering values [mV]
    # leave constants as is for clarity
    V_tc_eng_H = (tc_rawdec_H * 0.25 - 1024.0) / 61.606

    # convert engineering values to lab calibrated values [mV]
    V_tc_labcal_H = eval_poly(V_tc_eng_H, e2l_H)

    return V_tc_labcal_H


def sfl_thsph_temp_labcal_l(tc_rawdec_L, e2l_L):
    """
    Convert L thermocouple raw counts to lab-calibrated voltage.

    Called internally by sfl_thsph_temp_tl and sfl_thsph_temp_tcl.

    Parameters
    ----------
    tc_rawdec_L : array_like
        L thermocouple decimal counts (THSPHTE-TCL_L0) [counts].
    e2l_L : array_like
        Calibration coefficient. Array of polynomial coefficients converting
        L thermocouple engineering values to lab calibrated values.

    Returns
    -------
    V_tc_labcal_L : ndarray
        Lab-calibrated L thermocouple voltage [mV].
    """
    tc_rawdec_L = replace_fill_with_nan(None, tc_rawdec_L)

    # convert raw decimal output to engineering values [mV]
    # leave constants as is for clarity
    V_tc_eng_L = (tc_rawdec_L * 0.25 - 1024.0) / 61.606

    # convert engineering values to lab calibrated values [mV]
    V_tc_labcal_L = eval_poly(V_tc_eng_L, e2l_L)

    return V_tc_labcal_L


# .............................................................................
# THSPH data products: eval_poly, used for all THSPH data products ............
# .............................................................................


def eval_poly(x, c):
    """
    Evaluate a 5th-degree polynomial using Horner's algorithm.

    Supports both scalar and vectorized (per-element calibration coefficient)
    evaluation for THSPH data products.

    Parameters
    ----------
    x : scalar or array_like
        Argument(s) at which to evaluate the polynomial.
    c : array_like
        Polynomial coefficients in descending degree order. If x is a
        scalar, c is a 1D array of length 6. If x is a vector of length N,
        c is a 2D array of shape (N, 6) where row j contains the
        coefficients for x[j].

    Returns
    -------
    val : ndarray
        Evaluated polynomial values:
        c[:,0]*x^5 + c[:,1]*x^4 + ... + c[:,4]*x + c[:,5].

    Notes
    -----
    The np.atleast_2d call allows both single and vectorized calls to work
    with the same implementation.
    """
    # the "c = np.atleast_2d(c)" statement is necessary so that both single and
    # "vectorized" (in the OOI CI sense) calls to the eval_poly subroutine work.
    c = np.atleast_2d(c)

    # Horner's algorithm
    val = c[:, 5] + x * (c[:, 4] + x * (c[:, 3] + x * (c[:, 2] + x * (c[:, 1] + x * c[:, 0]))))

    return val

# .............................................................................
# TRHPH data products .........................................................
# .............................................................................


def sfl_trhph_vfltemp(V_ts, V_tc, tc_slope, ts_slope,
                      c0=0.015, c1=0.0024, c2=7.00e-5, c3=-1.00e-6):
    """
    Compute vent fluid temperature from TRHPH (TRHPHTE_L1).

    Parameters
    ----------
    V_ts : array_like
        Thermistor voltage (TRHPHVS_L0) [V].
    V_tc : array_like
        Thermocouple voltage (TRHPHVC_L0) [V].
    tc_slope : array_like
        Calibration coefficient. Thermocouple slope.
    ts_slope : array_like
        Calibration coefficient. Thermistor slope.
    c0 : array_like, optional
        Calibration coefficient. 3rd-degree polynomial correction term.
        Defaults to 0.015.
    c1 : array_like, optional
        Calibration coefficient. 3rd-degree polynomial correction term.
        Defaults to 0.0024.
    c2 : array_like, optional
        Calibration coefficient. 3rd-degree polynomial correction term.
        Defaults to 7.00e-5.
    c3 : array_like, optional
        Calibration coefficient. 3rd-degree polynomial correction term.
        Defaults to -1.00e-6.

    Returns
    -------
    T : ndarray
        Vent fluid temperature (TRHPHTE_L1) [deg_C].

    Notes
    -----
    Three cases apply based on V_tc and thermistor temperature T_ts:
    when V_tc <= 0, T = T_ts; when V_tc > 0 and T_ts > 10 deg_C, a
    polynomial correction is applied; when V_tc > 0 and 0 < T_ts <= 10
    deg_C, slope-based corrections are used. The default polynomial
    coefficients (c0-c3) are fixed constants from DPS 1341-00150.
    """
    # Test if polynomial coefficients are scalars (set via defaults), set to
    # same size as other inputs if required. Assumes if 'a' is a default, they
    # all are.
    if np.isscalar(c0):
        c0 = np.tile(c0, (V_ts.shape))

    if np.isscalar(c1):
        c1 = np.tile(c1, (V_ts.shape))

    if np.isscalar(c2):
        c2 = np.tile(c2, (V_ts.shape))

    if np.isscalar(c3):
        c3 = np.tile(c3, (V_ts.shape))

    # raw thermistor temperature
    T_ts = 27.50133 - 17.2658 * V_ts + 15.83424 / V_ts

    # where V_tc is less than or equal to 0, T = T_ts, otherwise...
    T = T_ts

    # Adjust raw thermistor temperature when V_tc is greater than 0 and ...
    tFlag = (V_tc > 0) & (T_ts > 10)  # T_ts is greater than 10
    poly = (c3[tFlag] * T_ts[tFlag]**3 + c2[tFlag] * T_ts[tFlag]**2 +
            c1[tFlag] * T_ts[tFlag] + c0[tFlag])
    T[tFlag] = (V_tc[tFlag] + poly) * 244.97

    tFlag = (V_tc > 0) & ((T_ts > 0) & (T_ts <= 10))  # T_ts is greater than 0 and less than 10
    T[tFlag] = (V_tc[tFlag] + V_tc[tFlag] * 244.97 * tc_slope[tFlag] + T_ts[tFlag]
                * ts_slope[tFlag]) * 244.97

    return T


def sfl_trhph_vfl_thermistor_temp(V_ts):
    """
    Compute TRHPH thermistor reference temperature (TRHPHTE-T_TS-AUX).

    This auxiliary product is the thermistor temperature alone, without
    thermocouple correction. It is the same as T_ts computed internally
    by sfl_trhph_vfltemp and is useful as an instrument diagnostic.

    Parameters
    ----------
    V_ts : array_like
        Thermistor voltage (TRHPHVS_L0) [V].

    Returns
    -------
    T_ts : ndarray
        Thermistor reference temperature (TRHPHTE-T_TS-AUX) [deg_C].
    """
    # thermistor temperature
    T_ts = 27.50133 - 17.2658 * V_ts + 15.83424 / V_ts
    return T_ts


def sfl_trhph_vflorp(V, offset, gain):
    """
    Compute vent fluid oxidation-reduction potential from TRHPH
    (TRHPHEH_L1).

    Parameters
    ----------
    V : array_like
        ORP sensor voltage (TRHPHVO_L0) [V].
    offset : array_like
        Calibration coefficient. Electronic offset [mV].
    gain : array_like
        Calibration coefficient. Gain multiplier [dimensionless].

    Returns
    -------
    ORP : ndarray
        Oxidation-reduction potential (TRHPHEH_L1) [mV].

    Notes
    -----
    Because the reference electrode is not a Standard Hydrogen Electrode,
    ORP values from this instrument cannot be directly compared with
    standard in-situ ORP measurements made with a Pt-Ag/AgCl electrode
    pair. The primary use of this product is to quantify change in ORP
    with respect to time (DPS 1341-00170).
    """
    # convert sensor voltage V from volts to mV;
    # subtract offset; undo gain multiplier.
    ORP = np.round((V * 1000.0 - offset)/gain)

    return ORP


def sfl_trhph_chloride(V_R1, V_R2, V_R3, T):
    """
    Compute vent fluid chloride concentration from TRHPH (TRHPHCC_L2).

    Parameters
    ----------
    V_R1 : array_like
        Resistivity voltage 1 (TRHPHR1_L0) [V].
    V_R2 : array_like
        Resistivity voltage 2 (TRHPHR2_L0) [V].
    V_R3 : array_like
        Resistivity voltage 3 (TRHPHR3_L0) [V].
    T : array_like
        Vent fluid temperature from TRHPH (TRHPHTE_L1) [deg_C].

    Returns
    -------
    Cl : ndarray
        Vent fluid chloride concentration (TRHPHCC_L2) [mmol kg^-1].
        Returns NaN where T is outside the calibration surface bounds
        (103 to 382 deg_C).

    Notes
    -----
    The three resistivity voltage channels are scaled to different
    measurement ranges; the optimal channel is selected based on V_R2.
    Chloride is determined by interpolating a conductivity isotherm from
    the Larson et al. (2007) temperature-conductivity-chloride calibration
    surface (imported from sfl_functions_surface.py as tdat, sdat, cdat)
    at the observed temperature, then mapping the measured conductivity
    to chloride concentration using RectBivariateSpline and np.interp.
    """
    # load sfl_functions_surface.py This loads the 3-dimensional calibration
    # surface of temperature, chloride, and conductivity reproduced as numpy
    # arrays from Larson_2007surface.mat.
    #
    # imported at module top
    # from ion_functions.data.sfl_functions_surface import tdat, sdat, cdat

    # select the optimal L0 Resistivity voltage.
    V_R = V_R3 / 5.0

    vflag = np.where((V_R2 >= 0.75) & (V_R2 < 3.90))
    V_R[vflag] = V_R2[vflag]

    vflag = np.where(V_R2 >= 3.90)
    V_R[vflag] = V_R1[vflag] * 5.0

    # convert resistivity to conductivity
    C = 1. / V_R

    # initialize product array Cl [mmol/kg] values to nans
    Cl = np.zeros(len(C)) + np.nan
    # set up chloride ['S' in units of mol/kg] range
    Scurve = np.linspace(np.min(sdat), np.max(sdat), 100,
                         endpoint='True')
    # create bivariate spline for interpolation
    f = RectBivariateSpline(tdat, sdat, cdat.T, kx=1, ky=1, s=0)

    # Note that when T is out-of-range, the interpolation np.interp does not
    # always give nan values for Cl as is required. Since Cl has been
    # initialized to nan values, iterate only over good T values, which also
    # improves speed.
    for ii in np.where(np.logical_and(T >= min(tdat), T <= max(tdat)))[0]:
        # find conductivity Ccurve as f(T=constant, chloride).
        Ccurve = f(T[ii], Scurve)
        # now interpolate measured conductivity C into (Ccurve,Scurve) to get
        # Cl. the conditional statement is in the DPS and therefore retained.
        if (np.all(np.isfinite(Ccurve))):
            Cl[ii] = np.interp(C[ii], Ccurve[0], Scurve, left=np.nan, right=np.nan)

    # change units to mmol/kg; round to required # of sigfigs as specified in
    # the DPS
    Cl = np.round(Cl * 1000.)

    return Cl


# .............................................................................
# PRESF data products .........................................................
# .............................................................................
def sfl_sflpres_rtime(p_psia):
    """
    Compute real-time seafloor pressure from PRESF (SFLPRES-RTIME_L1).

    Parameters
    ----------
    p_psia : array_like
        Real-time pressure (SFLPRES-RTIME_L0) [psia].

    Returns
    -------
    rtime : ndarray
        Real-time seafloor pressure (SFLPRES-RTIME_L1), hydrostatic plus
        atmospheric [dbar].
    """
    rtime = p_psia * 0.689475728
    return rtime


def sfl_sflpres_tide(p_dec_tide, b, m, slope=1.0, offset=0.0):
    """
    Compute tidal seafloor pressure from post-recovery PRESF data
    (SFLPRES-TIDE_L1).

    Parameters
    ----------
    p_dec_tide : array_like
        Tidal pressure decimal number (SFLPRES-TIDE_L0) [counts].
    b : array_like
        Calibration coefficient. Pressure scaling parameter B.
    m : array_like
        Calibration coefficient. Pressure scaling parameter M.
    slope : array_like, optional
        Calibration coefficient. Slope correction factor. Defaults to 1.0.
    offset : array_like, optional
        Calibration coefficient. Offset correction factor. Defaults to 0.0.

    Returns
    -------
    tide : ndarray
        Tidal seafloor pressure (SFLPRES-TIDE_L1), hydrostatic plus
        atmospheric [dbar].
    """
    # replace type integer fill values with nans
    p_dec_tide = replace_fill_with_nan(None, p_dec_tide)

    psia = slope * ((p_dec_tide - b) / m) + offset
    tide = 0.689475728 * psia
    return tide


def sfl_sflpres_wave(ptcn, p_dec_wave, u0, y1, y2, y3, c1, c2, c3, d1, d2,
                     t1, t2, t3, t4, poff, slope=1.0, offset=0.0):
    """
    Compute wave burst seafloor pressure from post-recovery PRESF data
    (SFLPRES-WAVE_L1).

    Parameters
    ----------
    ptcn : array_like
        Pressure temperature compensation number [counts].
    p_dec_wave : array_like
        Wave burst pressure decimal numbers (SFLPRES-WAVE_L0) [counts].
        May be 1D (single burst) or 2D (multiple bursts).
    u0 : array_like
        Calibration coefficient.
    y1 : array_like
        Calibration coefficient.
    y2 : array_like
        Calibration coefficient.
    y3 : array_like
        Calibration coefficient.
    c1 : array_like
        Calibration coefficient.
    c2 : array_like
        Calibration coefficient.
    c3 : array_like
        Calibration coefficient.
    d1 : array_like
        Calibration coefficient.
    d2 : array_like
        Calibration coefficient.
    t1 : array_like
        Calibration coefficient.
    t2 : array_like
        Calibration coefficient.
    t3 : array_like
        Calibration coefficient.
    t4 : array_like
        Calibration coefficient.
    poff : array_like
        Calibration coefficient. Pressure offset.
    slope : array_like, optional
        Calibration coefficient. Slope correction factor. Defaults to 1.0.
    offset : array_like, optional
        Calibration coefficient. Offset correction factor. Defaults to 0.0.

    Returns
    -------
    wave : ndarray
        Wave burst seafloor pressure (SFLPRES-WAVE_L1), hydrostatic plus
        atmospheric [dbar].
    """
    # replace type integer fill values with nans
    p_dec_wave, ptcn = replace_fill_with_nan(None, p_dec_wave, ptcn)

    # if p_dec_wave is a 1D array make it into a row vector
    p_dec_wave = np.atleast_2d(p_dec_wave)
    n_time_points, n_values_in_burst = p_dec_wave.shape

    # compute the pressure temperature compensation frequency (PTCF) and
    # pressure frequency (PF) from raw inputs
    PTCF = ptcn / 256.0
    PF = p_dec_wave / 256.0

    # use calibration coefficients to compute scale factors.
    U = ((1.0 / PTCF) * 1000000) - u0
    C = c1 + (c2 * U) + (c3 * U**2)
    D = d1 + d2
    T0 = (t1 + t2 * U + t3 * U**2 + t4 * U**3) / 1000000
    # broadcast T0 to the shape of PF
    T0 = np.tile(T0, (n_values_in_burst, 1)).T
    W = 1.0 - (T0**2 * PF**2)
    # broadcast C, D, and poff to the shape of W
    C = np.tile(C, (n_values_in_burst, 1)).T
    D = np.tile(D, (n_values_in_burst, 1)).T
    poff = np.tile(poff, (n_values_in_burst, 1)).T
    # broadcast slope and offset to the shape of W if not a scalar
    if np.atleast_1d(slope).shape[0] != 1:
        slope = np.tile(slope, (n_values_in_burst, 1)).T
        offset = np.tile(offset, (n_values_in_burst, 1)).T

    # compute the wave pressure data in dbar
    psia = slope * ((C * W * (1.0 - D * W)) + poff) + offset
    wave = 0.689475728 * psia
    return wave


def sfl_sbe26plus_prestmp(t0):
    """
    Compute seawater temperature from SBE 26plus tide record (PRESTMP_L1).

    Parameters
    ----------
    t0 : array_like
        Temperature number from tide record [hex converted to decimal].

    Returns
    -------
    t : ndarray
        Seawater temperature (PRESTMP_L1) [deg_C].
    """

    t = t0 / 1000.0 - 10.0
    return t
