#!/usr/bin/env python
"""
@package ion_functions.data.phsen_h_functions
@file ion_functions/data/phsen_h_functions.py
@author Samuel Dahlberg
@brief Module containing Cabled PHSEN H data-calculations. Instrument data is presented in raw units, and must be
calculated into L2 Data.
"""

import numpy as np
from math import e

# Various unit conversions
DBAR_TO_PSI = 1.450377
PSI_TO_DBAR = 0.6894759
OXYGEN_PHASE_TO_VOLTS = 39.457071
KELVIN_OFFSET_0C = 273.15
KELVIN_OFFSET_25C = 298.15
OXYGEN_MLPERL_TO_MGPERL = 1.42903
OXYGEN_MLPERL_TO_UMOLPERKG = 44660
# taken from https://blog.seabird.com/ufaqs/what-is-the-difference-in-temperature-expressions-between-ipts-68-and-its-90/
ITS90_TO_IPTS68 = 1.00024
# micro moles of nitrate to milligrams of nitrogen per liter
UMNO3_TO_MGNL = 0.014007
# [J K^{-1} mol^{-1}] Gas constant from SBS application note 99
R = 8.3144621
# [Coulombs mol^{-1}] Faraday constant from SBS application note 99
F = 96485.365


def temperature_raw_conversion(temp_counts, a0, a1, a2, a3):
    # The standard for the seaphox v2 is ITS90, C (degrees), and false on mv_r. These parameters will be removed since
    # this is only being used on seaphox v2
    # Source: https://github.com/Sea-BirdScientific/seabirdscientific/blob/808db391550ef574bc91240237b01cf46fb6e0e5/documentation/conversion-SeapHOxV2.ipynb


    temperature_counts = temp_counts

    log_t = np.log(temperature_counts)
    temperature = (
                          1 / (a0 + a1 * log_t + a2 * log_t ** 2 + a3 * log_t ** 3)
                  ) - KELVIN_OFFSET_0C

    return temperature


def pressure_raw_conversion(pres_counts, compensation_voltage, ptempa0, ptempa1, ptempa2, ptca0, ptca1, ptca2, ptcb0,
                            ptcb1, ptcb2, pa0, pa1, pa2):
    # The standard for the seaphox v2 is dbars. This parameter will be removed since this is only being used on seaphox v2
    # Source: https://github.com/Sea-BirdScientific/seabirdscientific/blob/808db391550ef574bc91240237b01cf46fb6e0e5/documentation/conversion-SeapHOxV2.ipynb

    sea_level_pressure = 14.7

    t = (
            ptempa0
            + ptempa1 * compensation_voltage
            + ptempa2 * compensation_voltage**2
    )
    x = pres_counts - ptca0 - ptca1 * t - ptca2 * t**2
    n = x * ptcb0 / (ptcb0 + ptcb1 * t + ptcb2 * t**2)
    pressure = pa0 + pa1 * n + pa2 * n**2

    pressure -= sea_level_pressure

    pressure *= PSI_TO_DBAR

    return pressure


def conductivity_raw_conversion(cond_counts, temperature, pressure, wbotc, g, h, i, j, ctcor, cpcor):
    f = cond_counts * np.sqrt(1 + wbotc * temperature) / 1000.0
    numerator = g + h * f ** 2 + i * f ** 3 + j * f ** 4
    denominator = 1 + ctcor * temperature + cpcor * pressure

    return numerator / denominator


def internal_temperature(temp_counts):
    """Converts the raw internal temperature counts to degrees Celcius

    :param temperature_counts: raw internal temperature counts
    :return: internal temperature in Celsius
    """
    slope = 175.72
    offset = -46.85
    int_16bit = 2.0 ** 16
    temperature = temp_counts / int_16bit * slope + offset

    return temperature


def internal_humidity(humidity_counts, temperature):
    """Convert relative humidity counts to percent

    :param humidity_counts: raw relative humidity counts
    :param temperature: converted internal temperature in Celcius
    :return: temperature compensated relative humidity in percent
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


def dissolved_oxygen(raw_oxygen_phase, thermistor, pressure, salinity, c0, c1, c2, coeff_e, a0, a1, a2, b0, b1,
                     therm_ta0, therm_ta1, therm_ta2, therm_ta3, thermistor_units="volts"):

    if thermistor_units == "volts":
        temperature = convert_sbe63_thermistor(thermistor, therm_ta0, therm_ta1, therm_ta2, therm_ta3)
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
            salinity * (sol_b0 + sol_b1 * ts + sol_b2 * ts ** 2 + sol_b3 * ts ** 3) + sol_c0 * salinity ** 2
    )
    s_corr = e ** s_corr_exp

    # temperature in Kelvin
    temperature_k = temperature + KELVIN_OFFSET_0C
    p_corr_exp = (coeff_e * pressure) / temperature_k
    p_corr = e ** p_corr_exp

    # fmt: off
    ox_val = (
            (((a0 + a1 * temperature + a2 * oxygen_volts ** 2)
              / (b0 + b1 * oxygen_volts) - 1.0) / ksv) * s_corr * p_corr
    )
    # fmt: on

    return ox_val


def convert_sbe63_thermistor(
        instrument_output,
        therm_ta0, therm_ta1, therm_ta2, therm_ta3
):
    """Converts a SBE63 thermistor raw output array to temperature in
    ITS-90 deg C.

    :param instrument_output: raw values from the thermistor
    :param coefs: calibration coefficients for the thermistor in the
        SBE63 sensor

    :return: converted thermistor temperature values in ITS-90 deg C
    """
    log_raw = np.log((100000 * instrument_output) / (3.3 - instrument_output))
    temperature = (
            1 / (therm_ta0 + therm_ta1 * log_raw + therm_ta2 * log_raw ** 2 + therm_ta3 * log_raw ** 3)
            - KELVIN_OFFSET_0C
    )
    return temperature


def ph_total(vrs_ext, degc, psu, dbar, k0, k2, f):
    """
    Calculate the total pH from the SeapHOx sensor. The total pH is calculated from the external voltage (vrs_ext),
    temperature (degC), salinity (psu), pressure (dbar), and the calibration coefficients (k0, k2, f).
    Source is Sea-Bird Scientific Application Note 99, "Calculating pH from ISFET pH Sensors".

    :param vrs_ext: external voltage from the FET sensor
    :param degc: temperature in degrees Celsius
    :param psu: salinity in practical salinity units
    :param dbar: pressure in decibars
    :param k0: calibration coefficient from vendor documentation
    :param k2: calibration coefficient from vendor documentation
    :param f: calibration coefficients (f0, f1, f2, f3, f4, f5) from the vendor documentation (as an array)
    :return: pH total
    """
    fp = f[0] * dbar + f[1] * dbar ** 2 + f[2] * dbar ** 3 + f[3] * dbar ** 4 + f[4] * dbar ** 5 + f[5] * dbar ** 6

    bar = dbar * 0.10  # convert pressure from dbar to bar

    # Nernstian response of the pH electrode (slope of the response)
    r = 8.3144621  # J/(mol K) universal gas constant
    t = degc + 273.15  # temperature in Kelvin
    f = 9.6485365e4  # C/mol Faraday constant
    snerst = r * t * np.log(10) / f

    # total chloride in seawater

    cl_total = (0.99889 / 35.453) * (psu / 1.80655) * (1000 / (1000 - 1.005 * psu))

    # partial Molal volume of HCl (calculated as Millero 1983)
    vhcl = 17.85 + 0.1044 * degc - 0.0001316 * degc ** 2

    # Sample ionic strength (calculated as Dickson et al. 2007)
    i = (19.924 * psu) / (1000 - 1.005 * psu)

    # Debye-Huckel constant for activity of HCl (calculated as Khoo et al. 1977)
    adh = 0.0000034286 * degc ** 2 + 0.00067503 * degc + 0.49172143

    # log of the activity coefficient of HCl as a function of temperature
    # (calculated as Khoo et al. 1977)
    loghclt = ((-adh * np.sqrt(i)) / (1 + 1.394 * np.sqrt(i))) + (0.08885 - 0.000111 * degc) * i

    # log10 of the activity coefficient of HCl as a function of temperature # and pressure (calculated as Johnson et al. 2017)
    loghcltp = loghclt + (((vhcl * bar) / (np.log(10) * r * t * 10)) / 2)

    # total sulfate in seawater (calculated as Dickson et al. 2007)
    so4_total = (0.1400 / 96.062) * (psu / 1.80655)

    # acid disassociation constant of HSO4- (calculated as Dickson et al.
    # 2007)
    ks = (1 - 0.001005 * psu) * np.exp(
        (-4276.1 / t) + 141.328 - 23.093 * np.log(t) + ((-13856 / t) + 324.57 - 47.986 * np.log(t)) * np.sqrt(i) + (
                (35474 / t) - 771.54 + 114.723 * np.log(t)) * i - (2698 / t) * i ** 1.5 + (1776 / t) * i ** 2)

    # partial Molal volume of HSO4- (calculated as Millero 1983)
    v_hso4 = -18.03 + 0.0466 * degc + 0.000316 * degc ** 2

    # compressibility of sulfate (calculated as Millero 1983)
    kbar_s = (-4.53 + 0.09 * degc) / 1000

    # acid dissociation constant of HSO4- as function of salinity,
    # temperature, and pressure (calculated as Millero 1982)
    kstp = ks * np.exp((-v_hso4 * bar + 0.5 * kbar_s * bar ** 2) / (r * t * 10))

    # calculate the pH total, adjusted for pressure, temperature and
    # salinity
    p_h = (((vrs_ext - k0 - k2 * degc - fp) / snerst) + np.log10(cl_total) + 2 * loghcltp - np.log10(
        1 + (so4_total / kstp)) - np.log10((1000 - 1.005 * psu) / 1000))

    return p_h
