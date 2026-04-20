#!/usr/bin/env python
"""
Module containing data processing functions for the Sea-Bird Scientific Deep
SeapHOx V2 pH instrument (PHSEN-G and PHSEN-H). Most functions with the
exception of ph_total were adapted from Sea-Bird Scientific's open-source
processing library: https://github.com/Sea-BirdScientific/seabirdscientific
"""

import numpy as np
from math import e
import gsw

# Various unit conversions
DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
OXYGEN_PHASE_TO_VOLTS = 39.457071
KELVIN_OFFSET_0C = 273.15
KELVIN_OFFSET_25C = 298.15
OXYGEN_MLPERL_TO_MGPERL = 1.42903
OXYGEN_MLPERL_TO_UMOLPERKG = 44660
# taken from https://blog.seabird.com/ufaqs/what-is-the-difference-in-temperature-expressions-between-ipts-68-and-its-90
ITS90_TO_IPTS68 = 1.00024
# micro moles of nitrate to milligrams of nitrogen per liter
UMNO3_TO_MGNL = 0.014007
# [J K^{-1} mol^{-1}] Gas constant from SBS application note 99
R = 8.3144621
# [Coulombs mol^{-1}] Faraday constant from SBS application note 99
F = 96485.365


def temperature_raw_conversion(temp_counts, a0, a1, a2, a3):
    """
    Convert raw SBE 37 temperature A/D counts to degrees C (ITS-90).

    Parameters
    ----------
    temp_counts : array_like
        Raw temperature A/D counts from the SBE 37 MicroCAT [counts].
    a0 : float
        Calibration coefficient 1 for the temperature sensor [unitless].
    a1 : float
        Calibration coefficient 2 for the temperature sensor [unitless].
    a2 : float
        Calibration coefficient 3 for the temperature sensor [unitless].
    a3 : float
        Calibration coefficient 4 for the temperature sensor [unitless].

    Returns
    -------
    temperature : ndarray
        Water temperature (ITS-90) [deg_C].

    Notes
    -----
    Calibration coefficients a0 through a3 are from factory calibration
    sheets supplied by Sea-Bird Scientific.
    """
    temperature_counts = temp_counts

    log_t = np.log(temperature_counts)
    temperature = (
        1 / (a0 + a1 * log_t + a2 * log_t ** 2 + a3 * log_t ** 3)
    ) - KELVIN_OFFSET_0C

    return temperature


def pressure_raw_conversion(pres_counts, compensation_voltage,
                            ptempa0, ptempa1, ptempa2,
                            ptca0, ptca1, ptca2,
                            ptcb0, ptcb1, ptcb2,
                            pa0, pa1, pa2):
    """
    Convert raw SBE 37 pressure A/D counts to pressure in dbar.

    Parameters
    ----------
    pres_counts : array_like
        Raw pressure A/D counts from the SBE 37 MicroCAT [counts].
    compensation_voltage : array_like
        Pressure temperature compensation voltage. Units are counts or
        volts depending on the instrument configuration [counts or V].
    ptempa0 : float
        Calibration coefficient 1 for the pressure temperature
        compensation [unitless].
    ptempa1 : float
        Calibration coefficient 2 for the pressure temperature
        compensation [unitless].
    ptempa2 : float
        Calibration coefficient 3 for the pressure temperature
        compensation [unitless].
    ptca0 : float
        Calibration coefficient 4 for the pressure sensor [unitless].
    ptca1 : float
        Calibration coefficient 5 for the pressure sensor [unitless].
    ptca2 : float
        Calibration coefficient 6 for the pressure sensor [unitless].
    ptcb0 : float
        Calibration coefficient 7 for the pressure sensor [unitless].
    ptcb1 : float
        Calibration coefficient 8 for the pressure sensor [unitless].
    ptcb2 : float
        Calibration coefficient 9 for the pressure sensor [unitless].
    pa0 : float
        Calibration coefficient 10 for the pressure sensor [unitless].
    pa1 : float
        Calibration coefficient 11 for the pressure sensor [unitless].
    pa2 : float
        Calibration coefficient 12 for the pressure sensor [unitless].

    Returns
    -------
    pressure : ndarray
        Sea pressure relative to one standard atmosphere (14.7 psia)
        [dbar].

    Notes
    -----
    Sea-level pressure (14.7 psia) is subtracted before unit conversion
    so the result is gauge pressure. Calibration coefficients are from
    factory calibration sheets supplied by Sea-Bird Scientific.
    """
    sea_level_pressure = 14.7

    t = (
        ptempa0
        + ptempa1 * compensation_voltage
        + ptempa2 * compensation_voltage ** 2
    )
    x = pres_counts - ptca0 - ptca1 * t - ptca2 * t ** 2
    n = x * ptcb0 / (ptcb0 + ptcb1 * t + ptcb2 * t ** 2)
    pressure = pa0 + pa1 * n + pa2 * n ** 2

    pressure -= sea_level_pressure
    pressure *= PSI_TO_DBAR

    return pressure


def conductivity_raw_conversion(cond_counts, temperature, pressure,
                                wbotc, g, h, i, j, ctcor, cpcor):
    """
    Convert raw SBE 37 conductivity A/D counts to conductivity in S/m.

    Parameters
    ----------
    cond_counts : array_like
        Raw conductivity A/D counts from the SBE 37 MicroCAT [counts].
    temperature : array_like
        Reference temperature from the SBE 37 MicroCAT [deg_C].
    pressure : array_like
        Reference pressure from the SBE 37 MicroCAT [dbar].
    wbotc : float
        Calibration coefficient 1 for the conductivity sensor
        [unitless].
    g : float
        Calibration coefficient 2 for the conductivity sensor
        [unitless].
    h : float
        Calibration coefficient 3 for the conductivity sensor
        [unitless].
    i : float
        Calibration coefficient 4 for the conductivity sensor
        [unitless].
    j : float
        Calibration coefficient 5 for the conductivity sensor
        [unitless].
    ctcor : float
        Calibration coefficient 6. Temperature correction factor for
        the conductivity sensor [unitless].
    cpcor : float
        Calibration coefficient 7. Pressure correction factor for the
        conductivity sensor [unitless].

    Returns
    -------
    conductivity : ndarray
        Conductivity [S/m].

    Notes
    -----
    Calibration coefficients are from factory calibration sheets
    supplied by Sea-Bird Scientific.
    """
    f = cond_counts * np.sqrt(1 + wbotc * temperature) / 1000.0
    numerator = g + h * f ** 2 + i * f ** 3 + j * f ** 4
    denominator = 1 + ctcor * temperature + cpcor * pressure

    return numerator / denominator


def internal_temperature(temp_counts):
    """
    Convert raw Deep SeapHOx V2 internal temperature counts to degrees C.

    Parameters
    ----------
    temp_counts : array_like
        Raw internal temperature counts from the Deep SeapHOx V2
        instrument [counts].

    Returns
    -------
    temperature : ndarray
        Internal instrument temperature [deg_C].

    Notes
    -----
    The conversion uses a 16-bit ADC with a fixed slope of 175.72 and
    offset of -46.85 derived from the SHT sensor calibration.
    """
    slope = 175.72
    offset = -46.85
    int_16bit = 2.0 ** 16
    temperature = temp_counts / int_16bit * slope + offset

    return temperature


def internal_humidity(humidity_counts, temperature):
    """
    Convert raw Deep SeapHOx V2 internal humidity counts to relative
    humidity in percent.

    Parameters
    ----------
    humidity_counts : array_like
        Raw relative humidity counts from the Deep SeapHOx V2
        instrument [counts].
    temperature : array_like
        Internal instrument temperature from internal_temperature
        [deg_C].

    Returns
    -------
    relative_humidity : ndarray
        Relative humidity, clipped to the range [0, 100] [percent].

    Notes
    -----
    Uncompensated relative humidity is first computed from the raw
    counts using a 16-bit ADC with slope 125 and offset -6. A
    temperature compensation of -0.15 percent per degree C deviation
    from 25 deg_C is then applied. Values outside [0, 119] percent
    before compensation are left unmodified; the final result is
    clipped to [0, 100] percent.
    """
    slope = 125
    offset = -6
    int_16bit = 2.0 ** 16
    max_humidity = 119
    temperature_coefficient = -0.15
    temperature_25c = 25

    # Uncompensated relative humidity
    relative_humidity = slope * humidity_counts / int_16bit + offset

    for n, humidity in enumerate(relative_humidity):
        # Theoretically, uncompensated relative humidity can be up to 119%
        if 0 <= humidity < max_humidity:
            relative_humidity[n] = humidity + temperature_coefficient * (
                temperature_25c - temperature[n]
            )

    np.clip(relative_humidity, a_min=0, a_max=100)

    return relative_humidity


def dissolved_oxygen(raw_oxygen_phase, thermistor, pressure, salinity,
                     c0, c1, c2, coeff_e, a0, a1, a2, b0, b1,
                     therm_ta0, therm_ta1, therm_ta2, therm_ta3,
                     lat, lon, thermistor_units="volts"):
    """
    Convert raw SBE 63 oxygen phase to dissolved oxygen in umol/kg.

    Parameters
    ----------
    raw_oxygen_phase : array_like
        Raw SBE 63 phase value [microseconds].
    thermistor : array_like
        Raw SBE 63 thermistor output used as the temperature reference.
        Units are volts when thermistor_units='volts' or degrees C when
        thermistor_units='C'.
    pressure : array_like
        Pressure from the co-located SBE 37 MicroCAT [dbar].
    salinity : array_like
        Practical salinity from the co-located SBE 37 MicroCAT
        [unitless].
    c0 : float
        Calibration coefficient 1 for the SBE 63 oxygen sensor
        [unitless].
    c1 : float
        Calibration coefficient 2 for the SBE 63 oxygen sensor
        [unitless].
    c2 : float
        Calibration coefficient 3 for the SBE 63 oxygen sensor
        [unitless].
    coeff_e : float
        Calibration coefficient 4. Pressure correction coefficient for
        the SBE 63 oxygen sensor [unitless].
    a0 : float
        Calibration coefficient 5 for the SBE 63 oxygen sensor
        [unitless].
    a1 : float
        Calibration coefficient 6 for the SBE 63 oxygen sensor
        [unitless].
    a2 : float
        Calibration coefficient 7 for the SBE 63 oxygen sensor
        [unitless].
    b0 : float
        Calibration coefficient 8 for the SBE 63 oxygen sensor
        [unitless].
    b1 : float
        Calibration coefficient 9 for the SBE 63 oxygen sensor
        [unitless].
    therm_ta0 : float
        Calibration coefficient 1 for the SBE 63 thermistor [unitless].
    therm_ta1 : float
        Calibration coefficient 2 for the SBE 63 thermistor [unitless].
    therm_ta2 : float
        Calibration coefficient 3 for the SBE 63 thermistor [unitless].
    therm_ta3 : float
        Calibration coefficient 4 for the SBE 63 thermistor [unitless].
    lat : float or array_like
        Latitude of the instrument [decimal degrees North].
    lon : float or array_like
        Longitude of the instrument [decimal degrees East].
    thermistor_units : str, optional
        Units of the thermistor input: 'volts' (default) or 'C'.

    Returns
    -------
    oxygen_umolkg : ndarray
        Dissolved oxygen concentration [umol/kg].

    Notes
    -----
    When thermistor_units is 'volts', the thermistor output is first
    converted to temperature using convert_sbe63_thermistor. The
    dissolved oxygen is computed from the phase value via a Stern-Volmer
    model with salinity and pressure corrections. The result is
    converted from ml/L to umol/kg using potential density computed
    from the GSW TEOS-10 library. Calibration coefficients are from
    factory calibration sheets supplied by Sea-Bird Scientific.
    """
    if thermistor_units == "volts":
        temperature = convert_sbe63_thermistor(
            thermistor, therm_ta0, therm_ta1, therm_ta2, therm_ta3
        )
    elif thermistor_units == "C":
        temperature = thermistor
    else:
        raise ValueError

    oxygen_volts = raw_oxygen_phase / OXYGEN_PHASE_TO_VOLTS  # from the manual

    ksv = c0 + c1 * temperature + c2 * temperature ** 2

    # The following correction coefficients are all constants
    sol_b0 = -6.24523e-3
    sol_b1 = -7.37614e-3
    sol_b2 = -1.0341e-2
    sol_b3 = -8.17083e-3
    sol_c0 = -4.88682e-7

    ts = np.log((KELVIN_OFFSET_25C - temperature) / (KELVIN_OFFSET_0C + temperature))
    s_corr_exp = (
        salinity * (sol_b0 + sol_b1 * ts + sol_b2 * ts ** 2 + sol_b3 * ts ** 3)
        + sol_c0 * salinity ** 2
    )
    s_corr = e ** s_corr_exp

    # temperature in Kelvin
    temperature_k = temperature + KELVIN_OFFSET_0C
    p_corr_exp = (coeff_e * pressure) / temperature_k
    p_corr = e ** p_corr_exp

    ox_val = (
        (((a0 + a1 * temperature + a2 * oxygen_volts ** 2)
          / (b0 + b1 * oxygen_volts) - 1.0) / ksv) * s_corr * p_corr
    )

    # Unit calculations to convert from ml/l to umol/kg
    absolute_salinity = gsw.SA_from_SP(salinity, pressure, lon, lat)
    conservative_temp = gsw.CT_from_t(salinity, temperature, pressure)
    pref = 0
    potential_density = gsw.rho(absolute_salinity, conservative_temp, pref)

    oxygen_umolkg = (ox_val * OXYGEN_MLPERL_TO_UMOLPERKG) / potential_density

    return oxygen_umolkg


def convert_sbe63_thermistor(instrument_output,
                             therm_ta0, therm_ta1, therm_ta2, therm_ta3):
    """
    Convert raw SBE 63 thermistor output to temperature in degrees C
    (ITS-90).

    Parameters
    ----------
    instrument_output : array_like
        Raw thermistor output from the SBE 63 optical dissolved oxygen
        sensor [V].
    therm_ta0 : float
        Calibration coefficient 1 for the SBE 63 thermistor [unitless].
    therm_ta1 : float
        Calibration coefficient 2 for the SBE 63 thermistor [unitless].
    therm_ta2 : float
        Calibration coefficient 3 for the SBE 63 thermistor [unitless].
    therm_ta3 : float
        Calibration coefficient 4 for the SBE 63 thermistor [unitless].

    Returns
    -------
    temperature : ndarray
        Temperature from the SBE 63 thermistor (ITS-90) [deg_C].

    Notes
    -----
    Calibration coefficients are from factory calibration sheets
    supplied by Sea-Bird Scientific.
    """
    log_raw = np.log((100000 * instrument_output) / (3.3 - instrument_output))
    temperature = (
        1 / (therm_ta0 + therm_ta1 * log_raw
             + therm_ta2 * log_raw ** 2 + therm_ta3 * log_raw ** 3)
        - KELVIN_OFFSET_0C
    )
    return temperature


def convert_ph_voltage_counts(ph_counts):
    """
    Convert raw ISFET pH voltage counts to volts for the Deep SeapHOx V2.

    Parameters
    ----------
    ph_counts : array_like
        Raw ISFET pH voltage counts from the Deep SeapHOx V2 instrument
        [counts].

    Returns
    -------
    ph_volts : ndarray
        ISFET external reference voltage (V_FET/REF) [Volts].

    Notes
    -----
    The conversion uses a 23-bit ADC with a 2.5 V reference and unity
    gain. The full-scale count is 8388608 (2^23).
    """
    adc_vref = 2.5
    gain = 1
    adc_23bit = 8388608.0
    ph_volts = adc_vref / gain * (ph_counts / adc_23bit - 1)
    return ph_volts


def ph_total(vrs_ext, degc, psu, dbar, k0, k2, f):
    """
    Compute the OOI L2 pH of seawater (PHWATER_L2) from the Sea-Bird
    Scientific Deep SeapHOx V2 (PHSEN-G and PHSEN-H).

    Parameters
    ----------
    vrs_ext : array_like
        External ISFET reference voltage (V_FET/REF) [Volts].
    degc : array_like
        In-situ temperature from the co-located CTD [deg_C].
    psu : array_like
        Practical salinity from the co-located CTD [unitless].
    dbar : array_like
        Pressure from the co-located CTD [dbar].
    k0 : float or array_like
        Calibration coefficient 1. Cell standard potential offset
        for the external reference [Volts].
    k2 : float or array_like
        Calibration coefficient 2. Temperature slope coefficient for
        the external reference [Volts deg_C^-1].
    f : array_like, shape (..., 6)
        Calibration coefficients 3-8. Coefficients f1 through f6 of
        the 6th-order pressure response polynomial f(P). Coefficient
        f0 is captured in k0 and is not used here [unitless].

    Returns
    -------
    p_h : ndarray
        pH of seawater on the total hydrogen ion scale (PHWATER_L2)
        [unitless].

    Notes
    -----
    The algorithm follows Sea-Bird Scientific Application Note 99
    (Johnson et al. 2016; Johnson et al. 2017). The ISFET external
    cell exhibits a Nernstian response to pH and is sensitive to
    chloride activity. The total pH is computed as:

        pH = (V_FET/REF - k0 - k2*t - f(P)) / S_nernst
             + log10(Cl_T) + 2*log10(gamma_HCl)_T&P
             - log10(1 + S_T/K_S,T&P)
             - log10((1000 - 1.005*S) / 1000)

    where S_nernst = R*T*ln(10)/F, Cl_T is total chloride, gamma_HCl
    is the HCl activity coefficient corrected for temperature and
    pressure, S_T is total sulfate, and K_S,T&P is the acid
    dissociation constant of HSO4- corrected for temperature and
    pressure. All intermediate quantities are computed from salinity,
    temperature, and pressure following Dickson et al. (2007),
    Khoo et al. (1977), Millero (1982, 1983), and Johnson et al.
    (2017).

    Calibration coefficients k0, k2, and f are from factory
    calibration sheets supplied by Sea-Bird Scientific.
    """
    f = np.atleast_2d(f)

    fp = (f[:, 0] * dbar + f[:, 1] * dbar ** 2 + f[:, 2] * dbar ** 3
          + f[:, 3] * dbar ** 4 + f[:, 4] * dbar ** 5 + f[:, 5] * dbar ** 6)

    bar = dbar * 0.10  # convert pressure from dbar to bar

    # Nernstian response of the pH electrode (slope of the response)
    r = 8.3144621  # J/(mol K) universal gas constant
    t = degc + 273.15  # temperature in Kelvin
    f = 9.6485365e4  # C/mol Faraday constant
    snerst = r * t * np.log(10) / f

    # total chloride in seawater (Dickson et al. 2007)
    cl_total = (0.99889 / 35.453) * (psu / 1.80655) * (1000 / (1000 - 1.005 * psu))

    # partial Molal volume of HCl (Millero 1983)
    vhcl = 17.85 + 0.1044 * degc - 0.0001316 * degc ** 2

    # Sample ionic strength (Dickson et al. 2007)
    i = (19.924 * psu) / (1000 - 1.005 * psu)

    # Debye-Huckel constant for activity of HCl (Khoo et al. 1977)
    adh = 0.0000034286 * degc ** 2 + 0.00067503 * degc + 0.49172143

    # log of HCl activity coefficient as a function of temperature
    # (Khoo et al. 1977)
    loghclt = ((-adh * np.sqrt(i)) / (1 + 1.394 * np.sqrt(i))) + (0.08885 - 0.000111 * degc) * i

    # log10 of HCl activity coefficient as a function of temperature
    # and pressure (Johnson et al. 2017)
    loghcltp = loghclt + (((vhcl * bar) / (np.log(10) * r * t * 10)) / 2)

    # total sulfate in seawater (Dickson et al. 2007)
    so4_total = (0.1400 / 96.062) * (psu / 1.80655)

    # acid dissociation constant of HSO4- (Dickson et al. 2007)
    ks = (1 - 0.001005 * psu) * np.exp(
        (-4276.1 / t) + 141.328 - 23.093 * np.log(t)
        + ((-13856 / t) + 324.57 - 47.986 * np.log(t)) * np.sqrt(i)
        + ((35474 / t) - 771.54 + 114.723 * np.log(t)) * i
        - (2698 / t) * i ** 1.5 + (1776 / t) * i ** 2)

    # partial Molal volume of HSO4- (Millero 1983)
    v_hso4 = -18.03 + 0.0466 * degc + 0.000316 * degc ** 2

    # compressibility of HSO4- (Millero 1983)
    kbar_s = (-4.53 + 0.09 * degc) / 1000

    # acid dissociation constant of HSO4- corrected for T and P
    # (Millero 1982)
    kstp = ks * np.exp((-v_hso4 * bar + 0.5 * kbar_s * bar ** 2) / (r * t * 10))

    # calculate total pH adjusted for pressure, temperature, and salinity
    p_h = (((vrs_ext - k0 - k2 * degc - fp) / snerst)
           + np.log10(cl_total) + 2 * loghcltp
           - np.log10(1 + (so4_total / kstp))
           - np.log10((1000 - 1.005 * psu) / 1000))

    return p_h
