#!/usr/bin/env python
"""
Module containing CO2 instrument family data processing functions.

Covers two OOI instrument classes and four data products:

* PCO2W (Sunburst SAMI2-CO2) -- pCO2 in seawater at depth (PCO2WAT_L1)
* PCO2A (Pro-Oceanus CO2-Pro) -- pCO2 in air and surface seawater
  (PCO2ATM_L1, PCO2SSW_L1) and air-sea CO2 flux (CO2FLUX_L2) derived
  from the air and surface seawater measurements along with data from
  a co-located bulk meteorology (METBK) sensor suite.
"""

import numpy as np
from ion_functions.utils import fill_value


def pco2_abs434_ratio(light):
    """
    Extract the L0 absorbance ratio at 434 nm from a PCO2W light array.

    Parameters
    ----------
    light : array_like, shape (N, M)
        Array of raw light measurements from the SAMI2-CO2 instrument.
        Column index 6 contains the 434 nm ratio counts (Ratio434).

    Returns
    -------
    a434ratio : ndarray, shape (N,)
        Optical absorbance ratio at 434 nm (CO2ABS1_L0) [unitless].

    Notes
    -----
    Extracts column index 6 from the instrument light measurement array.
    This L0 value is used by pco2_blank and pco2_calc_pco2 to compute
    the blank-corrected absorbance at 434 nm.
    """
    light = np.atleast_2d(light)
    a434ratio = light[:, 6]
    return a434ratio


def pco2_abs620_ratio(light):
    """
    Extract the L0 absorbance ratio at 620 nm from a PCO2W light array.

    Parameters
    ----------
    light : array_like, shape (N, M)
        Array of raw light measurements from the SAMI2-CO2 instrument.
        Column index 7 contains the 620 nm ratio counts (Ratio620).

    Returns
    -------
    a620ratio : ndarray, shape (N,)
        Optical absorbance ratio at 620 nm (CO2ABS2_L0) [unitless].

    Notes
    -----
    Extracts column index 7 from the instrument light measurement array.
    This L0 value is used by pco2_blank and pco2_calc_pco2 to compute
    the blank-corrected absorbance at 620 nm.
    """
    light = np.atleast_2d(light)
    a620ratio = light[:, 7]
    return a620ratio


def pco2_blank(raw_blank):
    """
    Convert a raw PCO2W blank count to a normalized absorbance blank.

    Parameters
    ----------
    raw_blank : array_like
        Raw optical absorbance blank counts at 434 or 620 nm [counts].

    Returns
    -------
    blank : ndarray
        Normalized optical absorbance blank at 434 or 620 nm [unitless].

    Notes
    -----
    Divides raw blank counts by 16384 (2^14) to normalize to a
    dimensionless ratio consistent with the blank-correction step in
    pco2_calc_pco2. This normalization reflects the vendor-corrected
    algorithm as of 2018.
    """
    blank = raw_blank / 16384.
    return blank


def pco2_thermistor(traw, sami_bits=12):
    """
    Convert raw PCO2W thermistor counts to temperature in degrees C.

    Parameters
    ----------
    traw : array_like
        Raw thermistor counts (CO2THRM_L0) [counts].
    sami_bits : int, optional
        ADC resolution of the SAMI hardware. Use 12 for original
        hardware and 14 for upgraded hardware. Default is 12.

    Returns
    -------
    therm : ndarray
        Thermistor temperature [degC].

    Notes
    -----
    The conversion differs by hardware generation. For 12-bit hardware
    the full-scale count is 4096; for 14-bit hardware it is 16384. In
    both cases the thermistor resistance ratio is multiplied by 17400
    before taking the natural log. The log value is used in a three-term
    inverse-temperature polynomial to obtain temperature in Kelvin,
    which is converted to degrees C by subtracting 273.15.
    """
    traw = np.atleast_1d(traw)
    sami_bits = np.atleast_1d(sami_bits)

    # Conversion depends on whether the SAMI is 12-bit or 14-bit hardware
    if sami_bits[0] == 14:
        rt = np.log((traw / (16384.0 - traw)) * 17400.0)
    else:
        rt = np.log((traw / (4096. - traw)) * 17400.)
    inv_t = 0.0010183 + 0.000241 * rt + 0.00000015 * rt**3
    therm = (1. / inv_t) - 273.15
    return therm


def pco2_battery(braw, sami_bits):
    """
    Convert raw PCO2W battery counts to battery voltage in Volts.

    Parameters
    ----------
    braw : array_like
        Raw battery voltage counts [counts].
    sami_bits : int
        ADC resolution of the SAMI hardware. Use 12 for original
        hardware and 14 for upgraded hardware.

    Returns
    -------
    volts : ndarray
        Battery voltage [V].

    Notes
    -----
    The scaling differs by hardware generation. For 14-bit hardware
    the full-scale range is 3 V over 4000 counts. For 12-bit hardware
    the full-scale range is 15 V over 4096 counts.
    """
    braw = np.atleast_1d(braw)
    sami_bits = np.atleast_1d(sami_bits)

    if sami_bits[0] == 14:
        volts = braw * 3. / 4000.
    else:
        volts = braw * 15. / 4096.
    return volts


def pco2_pco2wat(mtype, light, therm, ea434, eb434, ea620, eb620,
                 calt, cala, calb, calc, a434blank, a620blank):
    """
    OOI wrapper for PCO2WAT_L1; returns fill value for blank records.

    Calls pco2_calc_pco2 for all records, then replaces results for
    blank measurement records (mtype == 5) with the fill value.

    Parameters
    ----------
    mtype : array_like
        Measurement type: 4 for pCO2 measurement, 5 for blank
        [unitless].
    light : array_like, shape (N, M)
        Array of raw light measurements from the SAMI2-CO2 instrument.
    therm : array_like
        PCO2W thermistor temperature (CO2THRM_L1) [degC].
    ea434 : array_like
        Calibration coefficient 1.
    eb434 : array_like
        Calibration coefficient 2.
    ea620 : array_like
        Calibration coefficient 3.
    eb620 : array_like
        Calibration coefficient 4.
    calt : array_like
        Calibration coefficient 5 (temperature reference).
    cala : array_like
        Calibration coefficient 6.
    calb : array_like
        Calibration coefficient 7.
    calc : array_like
        Calibration coefficient 8.
    a434blank : array_like
        Blank measurement at 434 nm (CO2ABS1_L0) [counts].
    a620blank : array_like
        Blank measurement at 620 nm (CO2ABS2_L0) [counts].

    Returns
    -------
    pco2 : ndarray
        Partial pressure of CO2 in seawater (PCO2WAT_L1) [uatm].
        Blank records (mtype == 5) are set to the system fill value.

    See Also
    --------
    pco2_calc_pco2 : Core algorithm; use directly for all-record access.
    """
    mtype = np.atleast_1d(mtype)
    light = np.atleast_2d(light)
    therm = np.atleast_1d(therm)
    ea434 = np.atleast_1d(ea434)
    eb434 = np.atleast_1d(eb434)
    ea620 = np.atleast_1d(ea620)
    eb620 = np.atleast_1d(eb620)
    calt = np.atleast_1d(calt)
    cala = np.atleast_1d(cala)
    calb = np.atleast_1d(calb)
    calc = np.atleast_1d(calc)
    a434blank = np.atleast_1d(a434blank)
    a620blank = np.atleast_1d(a620blank)

    pco2 = pco2_calc_pco2(light, therm, ea434, eb434, ea620, eb620,
                          calt, cala, calb, calc, a434blank, a620blank)

    # Reset blank measurement records to the fill value
    m = np.where(mtype == 5)[0]
    pco2[m] = fill_value

    return pco2


def pco2_calc_pco2(light, therm, ea434, eb434, ea620, eb620,
                   calt, cala, calb, calc, a434blank, a620blank):
    """
    Compute PCO2WAT_L1 from the Sunburst SAMI2-CO2 (PCO2W).

    Calculates the OOI Level 1 Partial Pressure of CO2 in Seawater
    (PCO2WAT_L1) from raw light measurements, thermistor temperature,
    blank corrections, and factory calibration coefficients.

    Parameters
    ----------
    light : array_like, shape (N, M)
        Array of raw light measurements from the SAMI2-CO2 instrument.
        Column 6 is Ratio434 and column 7 is Ratio620.
    therm : array_like
        PCO2W thermistor temperature (CO2THRM_L1) [degC].
    ea434 : array_like
        Calibration coefficient 1.
    eb434 : array_like
        Calibration coefficient 2.
    ea620 : array_like
        Calibration coefficient 3.
    eb620 : array_like
        Calibration coefficient 4.
    calt : array_like
        Calibration coefficient 5 (temperature reference).
    cala : array_like
        Calibration coefficient 6.
    calb : array_like
        Calibration coefficient 7.
    calc : array_like
        Calibration coefficient 8.
    a434blank : array_like
        Normalized blank at 434 nm (from pco2_blank) [unitless].
    a620blank : array_like
        Normalized blank at 620 nm (from pco2_blank) [unitless].

    Returns
    -------
    pco2 : ndarray
        Partial pressure of CO2 in seawater (PCO2WAT_L1) [uatm].
        Blank records (where Ratio434 == Ratio620 after blank
        correction) are set to the system fill value.

    Notes
    -----
    The e constants (e1=0.0043, e2=2.136, e3=0.2105) are fixed values
    provided by Sunburst Sensors based on lab determinations with the
    bromothymol blue (BTB) indicator solution. They are not recalculated
    per calibration cycle. The final pCO2 is computed from a quadratic
    calibration curve using coefficients CalA, CalB, and CalC (cala,
    calb, calc). All calibration coefficients are from factory
    calibration sheets.

    The Python implementation corrects errors present in the vendor-
    supplied Matlab code (DPS Appendix A, 2018). The corrected algorithm
    was verified against the vendor and is authoritative.
    """
    # Fixed equilibration constants provided by Sunburst Sensors
    e1 = 0.0043
    e2 = 2.136
    e3 = 0.2105

    # Extract 434 nm and 620 nm ratios from light array
    ratio434 = light[:, 6]
    ratio620 = light[:, 7]

    # Blank-correct the absorbance ratios
    ar434 = (ratio434 / a434blank)
    ar4620 = (ratio620 / a620blank)

    # Map blank records (ar434 == ar4620) to avoid log domain errors;
    # these will be set to fill_value at the end
    m = np.where(ar434 == ar4620)[0]
    ar434[m] = 0.99999
    ar4620[m] = 0.99999

    # Compute blank-corrected absorbances and their ratio
    a434 = -1 * np.log10(ar434)
    a620 = -1 * np.log10(ar4620)
    ratio = a620 / a434

    # Compute pCO2 via temperature-corrected SAMI response
    v1 = ratio - e1
    v2 = e2 - e3 * ratio
    rco21 = -1 * np.log10(v1 / v2)
    rco22 = (therm - calt) * 0.007 + rco21
    t_coeff = 0.0075778 + 0.0012389 * rco22 - 0.00048757 * rco22**2
    t_cor_rco2 = rco21 + t_coeff * (therm - calt)
    pco2 = 10.**((-1. * calb + (calb**2 - (4. * cala * (calc - t_cor_rco2)))**0.5) / (2. * cala))
    pco2[m] = fill_value

    return np.real(pco2)


def pco2_ppressure(xco2, p, std=1013.25):
    """
    Compute PCO2ATM_L1 or PCO2SSW_L1 from the Pro-Oceanus CO2-Pro.

    Converts the L0 CO2 mole fraction (XCO2ATM or XCO2SSW) and L0 gas
    stream pressure (PRESAIR) into the L1 partial pressure of CO2 in
    air (PCO2ATM_L1) or surface seawater (PCO2SSW_L1).

    Parameters
    ----------
    xco2 : array_like
        CO2 mole fraction in air or surface seawater (XCO2ATM_L0 or
        XCO2SSW_L0) [ppm].
    p : array_like
        Gas stream pressure (PRESAIR_L0) [mbar].
    std : float, optional
        Standard atmospheric pressure [mbar/atm]. Default is 1013.25.

    Returns
    -------
    ppres : ndarray
        Partial pressure of CO2 in air or surface seawater
        (PCO2ATM_L1 or PCO2SSW_L1) [uatm].

    Notes
    -----
    Because xco2 is in ppm (10^-6), the result of xco2 * p / std is
    directly in uatm without further scaling. The instrument computes
    xco2 internally with pressure, temperature, and humidity
    compensation applied by onboard firmware. All temperature, humidity,
    and pressure compensation is performed onboard; ion-functions
    applies only the final unit conversion. Instrument results are valid
    between 0 and 35 degC.
    """
    ppres = xco2 * p / std
    return ppres


def pco2_co2flux(pco2w, pco2a, u10, t, s):
    """
    Compute CO2FLUX_L2, the air-sea flux of CO2.

    Estimates the OOI Level 2 flux of CO2 from the ocean to the
    atmosphere using L1 partial pressures of CO2 in seawater and air,
    and L1/L2 bulk meteorology inputs.

    Parameters
    ----------
    pco2w : array_like
        Partial pressure of CO2 in seawater (PCO2SSW_L1) [uatm].
    pco2a : array_like
        Partial pressure of CO2 in air (PCO2ATM_L1) [uatm].
    u10 : array_like
        Normalized wind speed at 10 m height (WIND10M_L2) [m s^-1].
    t : array_like
        Sea surface temperature (TEMPSRF_L1) [degC].
    s : array_like
        Sea surface salinity (SALSURF_L2) [psu].

    Returns
    -------
    flux : ndarray
        Estimated flux of CO2 from the ocean to the atmosphere
        (CO2FLUX_L2) [mol m^-2 s^-1]. Positive values indicate flux
        from the ocean into the atmosphere.

    Notes
    -----
    Follows the bulk formula of Wanninkhof (1992) with the gas transfer
    velocity parameterization of Sweeney et al. (2007) and the CO2
    solubility of Weiss (1974). pco2w and pco2a are converted from uatm
    to atm before computing the flux. The Schmidt number polynomial uses
    Wanninkhof (1992) Table A1 coefficients. Gas transfer velocity k is
    in cm h^-1 from Sweeney et al. (2007) and converted to m s^-1. The
    volume-based solubility formulation of Weiss (1974) is used
    (mol atm^-1 m^-3). An inherent uncertainty of approximately 10% is
    expected in the flux estimate (Wanninkhof, 1992).
    """
    # Convert uatm to atm
    pco2a = pco2a / 1.0e6
    pco2w = pco2w / 1.0e6

    # Schmidt number (Wanninkhof, 1992, Table A1)
    Sc = 2073.1 - (125.62 * t) + (3.6276 * t**2) - (0.043219 * t**3)

    # Gas transfer velocity in cm h^-1 (Sweeney et al., 2007, Fig. 3 and Table 1),
    # converted to m s^-1
    k = 0.27 * u10**2 * np.sqrt(660.0 / Sc)
    k = k / (100.0 * 3600.0)

    # Absolute temperature
    T = t + 273.15

    # CO2 solubility, volume version, mol atm^-1 m^-3
    # (Weiss, 1974, Eqn. 12 and Table I)
    T100 = T / 100
    K0 = 1000 * np.exp(-58.0931 + (90.5069 * (100 / T)) + (22.2940 * np.log(T100)) +
                       s * (0.027766 - (0.025888 * T100) + (0.0050578 * T100**2)))

    # Air-sea CO2 flux (Wanninkhof, 1992, Eqn. A2)
    flux = k * K0 * (pco2w - pco2a)
    return flux
