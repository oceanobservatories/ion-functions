#!/usr/bin/env python
"""
@package ion_functions.data.pco2_functions
@file ion_functions/data/pco2_functions.py
@author Christopher Wingard
@brief Module containing CO2 instrument family related functions
"""

import numpy as np
from ion_functions.utils import fill_value


# wrapper functions to extract parameters from SAMI-II CO2 instruments (PCO2W)
# and process these extracted parameters to calculate pCO2
def pco2_abs434_ratio(light):
    """
    Description:

        Extract the absorbance ratio at 434 nm from the pCO2 instrument light
        measurements. This will extract the CO2ABS1_L0 data product from the
        instrument light measurement arrays.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.

    Usage:

        a434ratio = pco2_abs434_ratio(light)

            where

        a434ratio = optical absorbance Ratio at 434 nm (CO2ABS1_L0) [unitless]
        light = array of light measurements

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    light = np.atleast_2d(light)
    a434ratio = light[:, 6]
    return a434ratio


def pco2_abs620_ratio(light):
    """
    Description:

        Extract the absorbance ratio at 620 nm from the pCO2 instrument light
        measurements. This will extract the CO2ABS2_L0 data product from the
        instrument light measurement arrays.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.

    Usage:

        a620ratio = pco2_abs620_ratio(light)

            where

        a620ratio = optical absorbance Ratio at 620 nm (CO2ABS2_L0) [unitless]
        light = array of light measurements

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    light = np.atleast_2d(light)
    a620ratio = light[:, 7]
    return a620ratio


def pco2_blank(raw_blank):
    """
    Description:

        Calculates the absorbance blank at 434 or 620 nm from the SAMI2-pCO2
        instrument.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.
        2014-02-28: Christopher Wingard. Updated to except raw blank values
                    from a sparse array.
        2018-03-04: Christopher Wingard. Updated to correctly calculate the
                    blank based on new code from the vendor.

    Usage:

        blank = pco2_blank(raw_blank)

            where

        blank = optical absorbance blank at 434 or 620 nm [unitless]
        raw_blank = raw optical absorbance blank at 434 or 620 nm [counts]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    blank = raw_blank / 16384.
    return blank


def pco2_thermistor(traw, sami_bits=12):
    """
    Description:

        Convert the thermistor data from counts to degrees Centigrade from the
        pCO2 instrument.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.
        2023-01-12: Mark Steiner. Add sami_bits arg to handle hardware upgrades
        2023-08-15: Samuel Dahlberg. Renamed local variables to follow naming convention.
                    Replaced use of Numexpr with Numpy.

    Usage:

        therm = pco2_thermistor(traw, sami_bits)

            where

        therm = converted thermistor temperature [degC]
        traw = raw thermistor temperature (CO2THRM_L0) [counts]
        sami_bits = number of bits on the SAMI hardware

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    # reset inputs to arrays
    traw = np.atleast_1d(traw)
    sami_bits = np.atleast_1d(sami_bits)

    # convert raw thermistor readings from counts to degrees Centigrade
    # conversion depends on whether the SAMI is older 12 bit or newer 14 bit hardware
    if sami_bits[0] == 14:
        rt = np.log((traw / (16384.0 - traw)) * 17400.0)
    else:
        rt = np.log((traw / (4096. - traw)) * 17400.)
    inv_t = 0.0010183 + 0.000241 * rt + 0.00000015 * rt**3
    therm = (1. / inv_t) - 273.15
    return therm


def pco2_battery(braw, sami_bits):
    """
    Description:

        Convert the battery voltage data from counts to volts from the
        pCO2 instrument.

    Implemented by:

        2023-02-23: Mark Steiner. Initial code.

    Usage:

        volts = pco2_battery_voltage(vraw, sami_bits)

            where

        volts = converted battery voltage [degC]
        braw = raw battery voltage [counts]
        sami_bits = number of bits on the SAMI hardware
    """
    # reset inputs to arrays
    braw = np.atleast_1d(braw)
    sami_bits = np.atleast_1d(sami_bits)

    # convert raw battery readings from counts to Volts
    if sami_bits[0] == 14:
        volts = braw * 3. / 4000.
    else:
        volts = braw * 15. / 4096.
    return volts


def pco2_pco2wat(mtype, light, therm, ea434, eb434, ea620, eb620,
                 calt, cala, calb, calc, a434blank, a620blank):
    """
    Description:

        Function to calculate the L1 PCO2WAT core data from the pCO2 instrument
        if the measurement type is 4 (pCO2 measurement), otherwise it is a
        blank and return a fill value.

    Implemented by:

        2013-04-20: Christopher Wingard. Initial code.
        2014-02-19: Christopher Wingard. Updated comments.
        2014-03-19: Christopher Wingard. Optimized using feedback provided by
                    Chris Fortin.
        2017-04-04: Pete Cable. Updated algorithm to use thermistor/blank counts
                    as indicated in the DPS and the usage below.

    Usage:

        pco2 = pco2_pco2wat(mtype, light, therm, ea434, eb434, ea620, eb620,
                            calt, cala, calb, calc, a434blank, a620blank)

            where

        pco2 = measured pco2 in seawater (PCO2WAT_L1) [uatm]
        mtype = measurement type, where 4 == actual measurement and 5 == a
            blank measurement [unitless]
        light = array of light measurements
        therm = PCO2W thermistor temperature (CO2THRM_L0) [counts]
        ea434 = Reagent specific calibration coefficient
        eb434 = Reagent specific calibration coefficient
        ea620 = Reagent specific calibration coefficient
        eb620 = Reagent specific calibration coefficient
        calt = Instrument specific calibration coefficient for temperature
        cala = Instrument specific calibration coefficient for the pCO2 measurements
        calb = Instrument specific calibration coefficient for the pCO2 measurements
        calc = Instrument specific calibration coefficient for the pCO2 measurements
        a434blank = Blank measurements at 434 nm (CO2ABS1_L0) [counts]
        a620blank = Blank measurements to 620 nm (CO2ABS2_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    # reset inputs to arrays
    # measurements
    mtype = np.atleast_1d(mtype)
    light = np.atleast_2d(light)
    therm = np.atleast_1d(therm)
    # calibration coefficients
    ea434 = np.atleast_1d(ea434)
    eb434 = np.atleast_1d(eb434)
    ea620 = np.atleast_1d(ea620)
    eb620 = np.atleast_1d(eb620)
    calt = np.atleast_1d(calt)
    cala = np.atleast_1d(cala)
    calb = np.atleast_1d(calb)
    calc = np.atleast_1d(calc)
    # blank measurements
    a434blank = np.atleast_1d(a434blank)
    a620blank = np.atleast_1d(a620blank)

    # calculate the pco2 value
    pco2 = pco2_calc_pco2(light, therm, ea434, eb434, ea620, eb620,
                          calt, cala, calb, calc, a434blank, a620blank)

    # reset dark measurements to the fill value
    m = np.where(mtype == 5)[0]
    pco2[m] = fill_value

    return pco2


# L1a PCO2WAT calculation
def pco2_calc_pco2(light, therm, ea434, eb434, ea620, eb620,
                   calt, cala, calb, calc, a434blank, a620blank):
    """
    Description:

        OOI Level 1 Partial Pressure of CO2 (pCO2) in seawater core data
        product, which is calculated from the Sunburst SAMI-II CO2 instrument
        (PCO2W).

    Implemented by:

        20??-??-??: J. Newton (Sunburst Sensors, LLC). Original Matlab code.
        2013-04-20: Christopher Wingard. Initial python code.
        2014-02-19: Christopher Wingard. Updated comments.
        2014-03-19: Christopher Wingard. Optimized.
        2018-03-04: Christopher Wingard. Updated to correctly calculate pCO2 using
                    newly formulated code provided by the vendor. Original vendor code
                    incorrectly calculated the blank correction. Applies additional
                    corrections to calculations to avoid errors thrown when running a
                    blank measurement.
        2023-01-12: Mark Steiner. Arg therm in degrees C instead of counts
        2023-08-15: Samuel Dahlberg. Changed local variable names to follow naming convention.

    Usage:

        pco2 = pco2_pco2wat(light, therm, ea434, eb434, ea620, eb620,
                            calt, cala, calb, calc, a434blank, a620blank)

            where

        pco2 = measured pco2 in seawater (PCO2WAT_L1) [uatm]
        light = array of light measurements
        therm = PCO2W thermistor temperature (CO2THRM_L1) [degrees C]
        ea434 = Reagent specific calibration coefficient
        eb434 = Reagent specific calibration coefficient
        ea620 = Reagent specific calibration coefficient
        eb620 = Reagent specific calibration coefficient
        calt = Instrument specific calibration coefficient for temperature
        cala = Instrument specific calibration coefficient for the pCO2 measurements
        calb = Instrument specific calibration coefficient for the pCO2 measurements
        calc = Instrument specific calibration coefficient for the pCO2 measurements
        a434blank = Blank measurements at 434 nm (CO2ABS1_L0) [counts]
        a620blank = Blank measurements to 620 nm (CO2ABS2_L0) [counts]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00510.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)
    """
    # set constants -- original vendor formulation, reset below
    # ea434 = ea434 - 29.3 * calt
    # eb620 = eb620 - 70.6 * calt
    # e1 = ea620 / ea434
    # e2 = eb620 / ea434
    # e3 = eb434 / ea434

    # set the e constants, values provided by Sunburst
    e1 = 0.0043
    e2 = 2.136
    e3 = 0.2105

    # Extract variables from light array
    ratio434 = light[:, 6]     # 434nm Ratio
    ratio620 = light[:, 7]     # 620nm Ratio

    # correct the absorbance ratios using the blanks
    ar434 = (ratio434 / a434blank)
    ar4620 = (ratio620 / a620blank)

    # map out blank measurements and spoof the ratios to avoid throwing an error
    m = np.where(ar434 == ar4620)[0]
    ar434[m] = 0.99999
    ar4620[m] = 0.99999

    # Calculate the final absorbance ratio
    a434 = -1 * np.log10(ar434)  # 434 absorbance
    a620 = -1 * np.log10(ar4620)  # 620 absorbance
    ratio = a620 / a434          # Absorbance ratio

    # calculate pCO2
    v1 = ratio - e1
    v2 = e2 - e3 * ratio
    rco21 = -1 * np.log10(v1 / v2)
    rco22 = (therm - calt) * 0.008 + rco21
    t_coeff = 0.0075778 - 0.0012389 * rco22 - 0.00048757 * rco22**2
    t_cor_rco2 = rco21 + t_coeff * (therm - calt)
    pco2 = 10.**((-1. * calb + (calb**2 - (4. * cala * (calc - t_cor_rco2)))**0.5) / (2. * cala))
    pco2[m] = fill_value  # reset the blanks captured earlier to a fill value

    return np.real(pco2)


def pco2_ppressure(xco2, p, std=1013.25):
    """
    Description:

        OOI Level 1 Partial Pressure of CO2 in Air (PCO2ATM) or Surface
        Seawater (PCO2SSW) core date product is computed by using an
        equation that incorporates the Gas Stream Pressure (PRESAIR) and the
        CO2 Mole Fraction in Air or Surface Seawater (XCO2ATM or XCO2SSW,
        respectively). It is computed using data from the pCO2 air-sea (PCO2A)
        family of instruments.

    Implemented by:

        2014-10-27: Christopher Wingard. Initial python code.
        2023-08-15: Samuel Dahlberg. Removed use of Numexpr.

    Usage:

        ppres = pco2_ppressure(xco2, p, std)

            where

        ppres = partial pressure of CO2 in air or surface seawater [uatm]
                (PCO2ATM_L1 or PCO2SSW_L1)
        xco2 = CO2 mole fraction in air or surface seawater [ppm] (XCO2ATM_LO
               or XCO2SSW_L0)
        p = gas stream pressure [mbar] (PRESAIR_L0)
        std = standard atmospheric pressure set to default of 1013.25 [mbar/atm]

    References:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Air and Surface Seawater. Document Control Number 1341-00260.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00260_Data_Product_SPEC_PCO2ATM_PCO2SSW_OOI.pdf)
    """
    ppres = xco2 * p / std
    return ppres


def pco2_co2flux(pco2w, pco2a, u10, t, s):
    """
    Description:

        OOI Level 2 core date product CO2FLUX is an estimate of the CO2 flux
        from the ocean to the atmosphere. It is computed using data from the
        pCO2 air-sea (PCO2A) and bulk meteorology (METBK) families of
        instruments.

    Implemented by:

        2012-03-28: Mathias Lankhorst. Original Matlab code.
        2013-04-20: Christopher Wingard. Initial python code.

    Usage:

        flux = pco2_co2flux(pco2w, pco2a, u10, t, s)

            where

        flux = estimated flux of CO2 from the ocean to atmosphere [mol m-2 s-1]
               (CO2FLUX_L2)
        pco2w = partial pressure of CO2 in sea water [uatm] (PCO2SSW_L1)
        pco2a = partial pressure of CO2 in air [uatm] (PCO2ATM_L1)
        u10 = normalized wind speed at 10 m height from METBK [m s-1] (WIND10M_L2)
        t = sea surface temperature from METBK [deg_C] (TEMPSRF_L1)
        s = sea surface salinity from METBK [psu] (SALSURF_L2)

    References:

        OOI (2012). Data Product Specification for Flux of CO2 into the
            Atmosphere. Document Control Number 1341-00270.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00270_Data_Product_SPEC_CO2FLUX_OOI.pdf)
    """
    # convert micro-atm to atm
    pco2a = pco2a / 1.0e6
    pco2w = pco2w / 1.0e6

    # Compute Schmidt number (after Wanninkhof, 1992, Table A1)
    Sc = 2073.1 - (125.62 * t) + (3.6276 * t**2) - (0.043219 * t**3)

    # Compute gas transfer velocity (after Sweeney et al. 2007, Fig. 3 and Table 1)
    k = 0.27 * u10**2 * np.sqrt(660.0 / Sc)

    # convert cm h-1 to m s-1
    k = k / (100.0 * 3600.0)

    # Compute the absolute temperature
    T = t + 273.15

    # Compute solubility (after Weiss 1974, Eqn. 12 and Table I).
    # Note that there are two versions, one for units per volume and
    # one per mass. Here, the volume version is used.
    # mol atm-1 m-3
    T100 = T / 100
    K0 = 1000 * np.exp(-58.0931 + (90.5069 * (100/T)) + (22.2940 * np.log(T100)) +
                       s * (0.027766 - (0.025888 * T100) + (0.0050578 * T100**2)))

    # mol atm-1 kg-1
    #K0 = np.exp(-60.2409 + (93.4517 * (100/T)) + (23.3585 * np.log(T100)) +
    #            s * (0.023517 - (0.023656 * T100) + (0.0047036 * T100**2)))

    # Compute flux (after Wanninkhof, 1992, eqn. A2)
    flux = k * K0 * (pco2w - pco2a)
    return flux
