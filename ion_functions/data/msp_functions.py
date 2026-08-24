#!/usr/bin/env python
"""
ion_functions.data.msp_functions

Functions supporting the MASSP (Mass Spectrometer) instrument class. The
MASSP measures dissolved gas concentrations (DISSGAS, L1) in hydrothermal
vent and cold seep fluids from residual gas analyzer (RGA) mass spectral
scans, together with a set of auxiliary quality, timing, and pH-intensity
products, and a derived L2 total dissolved gas concentration product
(TOTLGAS) for hydrogen sulfide and carbon dioxide.

This module has never been used in an operational OOI data pipeline; see
docs/api/msp_functions.md for background and the DEPRECATIONS.md entry for
this module.
"""

# import main python modules
import numpy as np

from ion_functions import deprecated

@deprecated
def calc_l2_totlgas_smph2scon(port_timestamp_sampleint, L0_dissgas_sampleint,
                              gas_mode_sampleint, port_timestamp_sampleint_mcu,
                              ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass,
                              massp_rga_steps_per_amu, calibration_table,
                              l2_ph_calibration_table, sensor_depth, salinity):
    """
    Computes the L2 total dissolved hydrogen sulfide concentration
    in the sample water from the L1 DISSGAS-SMPH2SCON
    concentration, the equilibrated water pH from calc_l2_mswater_smpphval,
    and a temperature/salinity/pressure-dependent speciation
    correction.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 timing and intensity inputs; see calc_dissgas_smph2scon for
        the full parameter list, which this function passes
        through unchanged.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.
    l2_ph_calibration_table : ndarray
        Six-element L2 pH calibration coefficient array.
    sensor_depth : float
        In situ depth [m].
    salinity : float
        Estimated practical salinity [PSU] used in the speciation
        correction; the DPS-era code comment notes an assumed
        value of 35 where a measured value is unavailable.

    Returns
    -------
    totlgas_smph2scon : float
        TOTLGAS-SMPH2SCON, total dissolved hydrogen sulfide
        concentration [uM].

    Notes
    -----
    Speciation factor beta is 1 + (K1 / 10^-pH), with K1 an empirical
    """

    ph_temp_array = calc_l2_mswater_smpphval(port_timestamp_sampleint,
                                             L0_dissgas_sampleint,
                                             gas_mode_sampleint,
                                             port_timestamp_sampleint_mcu,
                                             ph_meter_sampleint_mcu,
                                             inlet_temp_sampleint_mcu,
                                             massp_rga_initial_mass,
                                             massp_rga_final_mass,
                                             massp_rga_steps_per_amu,
                                             calibration_table,
                                             l2_ph_calibration_table)

    t = ph_temp_array[1]
    ph = ph_temp_array[0]

    smph2scon = calc_dissgas_smph2scon(port_timestamp_sampleint,
                                       L0_dissgas_sampleint,
                                       gas_mode_sampleint,
                                       port_timestamp_sampleint_mcu,
                                       ph_meter_sampleint_mcu,
                                       inlet_temp_sampleint_mcu,
                                       massp_rga_initial_mass,
                                       massp_rga_final_mass,
                                       massp_rga_steps_per_amu,
                                       calibration_table, sensor_depth)

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #estimated salinity == 35
    PSU = salinity

    k1T = 10**(-19.83 - (930.8 / (t+273.15)) + (2.8 * np.log(t + 273.15)) -
              (np.sqrt(PSU) * (-0.2391 + 35.685 / (t+273.15))) - (PSU * (0.0109 - (0.3776
               / (t+273.15)))))

    r = ((11.07 + 0.009 * t + 0.000942 * t**2) * 0.0689475729 * pressure +
        (-6.869 * 10**(-6) + 1.2835 * 10**(-7) * t) * pressure**2) / ((t+273.15) * 83.131)

    k1 = np.exp(r) * k1T

    beta = 1 + (k1 / 10**-ph)

    totlgas_smph2scon = beta * smph2scon

    return totlgas_smph2scon

@deprecated
def calc_l2_totlgas_smpco2con(port_timestamp_sampleint, L0_dissgas_sampleint,
                              gas_mode_sampleint, port_timestamp_sampleint_mcu,
                              ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass,
                              massp_rga_steps_per_amu, calibration_table,
                              l2_ph_calibration_table, sensor_depth, salinity):
    """
    Computes the L2 total dissolved carbon dioxide concentration
    in the sample water from the L1 DISSGAS-SMPCO2CON
    concentration, the equilibrated water pH from calc_l2_mswater_smpphval,
    and a temperature/salinity/pressure-dependent speciation
    correction.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 timing and intensity inputs; see calc_dissgas_smpco2con for
        the full parameter list, which this function passes
        through unchanged.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.
    l2_ph_calibration_table : ndarray
        Six-element L2 pH calibration coefficient array.
    sensor_depth : float
        In situ depth [m].
    salinity : float
        Estimated practical salinity [PSU] used in the speciation
        correction; the DPS-era code comment notes an assumed
        value of 35 where a measured value is unavailable.

    Returns
    -------
    totlgas_smpco2con : float
        TOTLGAS-SMPCO2CON, total dissolved carbon dioxide
        concentration [uM].

    Notes
    -----
    Speciation factor alpha is 1 + (K1/10^-pH) + (K1*K2/10^-pH^2),
    """

    ph_temp_array = calc_l2_mswater_smpphval(port_timestamp_sampleint,
                                             L0_dissgas_sampleint,
                                             gas_mode_sampleint,
                                             port_timestamp_sampleint_mcu,
                                             ph_meter_sampleint_mcu,
                                             inlet_temp_sampleint_mcu,
                                             massp_rga_initial_mass,
                                             massp_rga_final_mass,
                                             massp_rga_steps_per_amu,
                                             calibration_table,
                                             l2_ph_calibration_table)

    t = ph_temp_array[1]
    ph = ph_temp_array[0]

    smpco2con = calc_dissgas_smpco2con(port_timestamp_sampleint,
                                       L0_dissgas_sampleint,
                                       gas_mode_sampleint,
                                       port_timestamp_sampleint_mcu,
                                       ph_meter_sampleint_mcu,
                                       inlet_temp_sampleint_mcu,
                                       massp_rga_initial_mass,
                                       massp_rga_final_mass,
                                       massp_rga_steps_per_amu,
                                       calibration_table, sensor_depth)

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #estimated salinity == 35
    PSU = salinity

    K1T = np.exp((2.83655 - 2307.1266 / (t + 273.15) - 1.5529413 * np.log(t + 273.15) -
                  (0.20760841 * 4.0484 / (t + 273.15)) * np.sqrt(PSU) + 0.0846834 *
                  PSU - 0.00654208 * np.sqrt(PSU**3) + np.log(1 - 0.001005 * PSU)))

    K2T = np.exp((-9.226508 - 3351.616 / (t + 273.15) - 0.2005743 * np.log(t + 273.15) -
                  (0.106901773 * 23.9722 / (t + 273.15)) * np.sqrt(PSU) + 0.1130822 *
                  PSU - 0.00846934 * np.sqrt(PSU**3) + np.log(1 - 0.001005 * PSU)))

    r1 = (pressure * (1.758163 - 0.008763 * t - pressure * ((7.32 * 10**-6) -
         (2.0845 * 10**-7 * t))) / ((t + 273.15) * 83.131))

    r2 = (pressure * (1.09075 + 0.00151 * t + pressure * ((2.69 * 10**-6) -
         (3.506 * 10**-7 * t))) / ((t + 273.15) * 83.131))

    K1 = np.exp(r1) * K1T

    K2 = np.exp(r2) * K2T

    alpha = 1 + (K1 / (10**-ph)) + ((K1 * K2) / (10**-ph)**2)

    totlgas_smpco2con = alpha * smpco2con

    return totlgas_smpco2con

@deprecated
def calc_timestamp_totlgas_smph2scon(port_timestamp_sampleint,
                                     L0_dissgas_sampleint,
                                     gas_mode_sampleint,
                                     port_timestamp_sampleint_mcu,
                                     ph_meter_sampleint_mcu,
                                     inlet_temp_sampleint_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    """
    OOI wrapper for TSTAMP-SMPH2SCON. Returns the Direct mode scan timestamp
    associated with the L2 total dissolved hydrogen sulfide concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the mean port
        timestamp over the Direct mode averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint,
                                        gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_timestamp_totlgas_smpco2con(port_timestamp_sampleint,
                                     L0_dissgas_sampleint,
                                     gas_mode_sampleint,
                                     port_timestamp_sampleint_mcu,
                                     ph_meter_sampleint_mcu,
                                     inlet_temp_sampleint_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    """
    OOI wrapper for TSTAMP-SMPCO2CON. Returns the Direct mode scan timestamp
    associated with the L2 total dissolved carbon dioxide concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the mean port
        timestamp over the Direct mode averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint,
                                        gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_l2_totlgas_bkgh2scon(port_timestamp_bkgndint, L0_dissgas_bkgndint,
                              gas_mode_bkgndint,
                              port_timestamp_bkgndint_mcu,
                              ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass,
                              massp_rga_steps_per_amu, calibration_table,
                              l2_ph_calibration_table, sensor_depth, salinity):
    """
    Computes the L2 total dissolved hydrogen sulfide concentration
    in the background water from the L1 DISSGAS-BKGH2SCON
    concentration, the equilibrated water pH from calc_l2_mswater_bkgphval,
    and a temperature/salinity/pressure-dependent speciation
    correction.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 timing and intensity inputs; see calc_dissgas_bkgh2scon for
        the full parameter list, which this function passes
        through unchanged.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.
    l2_ph_calibration_table : ndarray
        Six-element L2 pH calibration coefficient array.
    sensor_depth : float
        In situ depth [m].
    salinity : float
        Estimated practical salinity [PSU] used in the speciation
        correction; the DPS-era code comment notes an assumed
        value of 35 where a measured value is unavailable.

    Returns
    -------
    totlgas_bkgh2scon : float
        TOTLGAS-BKGH2SCON, total dissolved hydrogen sulfide
        concentration [uM].

    Notes
    -----
    Speciation factor beta is 1 + (K1 / 10^-pH), with K1 an empirical
    """

    ph_temp_array = calc_l2_mswater_bkgphval(port_timestamp_bkgndint,
                                             L0_dissgas_bkgndint,
                                             gas_mode_bkgndint,
                                             port_timestamp_bkgndint_mcu,
                                             ph_meter_bkgndint_mcu,
                                             inlet_temp_bkgndint_mcu,
                                             massp_rga_initial_mass,
                                             massp_rga_final_mass,
                                             massp_rga_steps_per_amu,
                                             calibration_table,
                                             l2_ph_calibration_table)

    t = ph_temp_array[1]
    ph = ph_temp_array[0]

    bkgh2scon = calc_dissgas_bkgh2scon(port_timestamp_bkgndint,
                                       L0_dissgas_bkgndint,
                                       gas_mode_bkgndint,
                                       port_timestamp_bkgndint_mcu,
                                       ph_meter_bkgndint_mcu,
                                       inlet_temp_bkgndint_mcu,
                                       massp_rga_initial_mass,
                                       massp_rga_final_mass,
                                       massp_rga_steps_per_amu,
                                       calibration_table, sensor_depth)

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #estimated salinity == 35
    PSU = salinity

    k1T = 10**(-19.83 - (930.8 / (t+273.15)) + (2.8 * np.log(t + 273.15)) -
              (np.sqrt(PSU) * (-0.2391 + 35.685 / (t+273.15))) - (PSU * (0.0109 - (0.3776
               / (t+273.15)))))

    r = ((11.07 + 0.009 * t + 0.000942 * t**2) * 0.0689475729 * pressure +
        (-6.869 * 10**(-6) + 1.2835 * 10**(-7) * t) * pressure**2) / ((t+273.15) * 83.131)

    k1 = np.exp(r) * k1T

    beta = 1 + (k1 / 10**-ph)

    totlgas_bkgh2scon = beta * bkgh2scon

    return totlgas_bkgh2scon

@deprecated
def calc_l2_totlgas_bkgco2con(port_timestamp_bkgndint, L0_dissgas_bkgndint,
                              gas_mode_bkgndint,
                              port_timestamp_bkgndint_mcu,
                              ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass,
                              massp_rga_steps_per_amu, calibration_table,
                              l2_ph_calibration_table, sensor_depth, salinity):
    """
    Computes the L2 total dissolved carbon dioxide concentration
    in the background water from the L1 DISSGAS-BKGCO2CON
    concentration, the equilibrated water pH from calc_l2_mswater_bkgphval,
    and a temperature/salinity/pressure-dependent speciation
    correction.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 timing and intensity inputs; see calc_dissgas_bkgco2con for
        the full parameter list, which this function passes
        through unchanged.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.
    l2_ph_calibration_table : ndarray
        Six-element L2 pH calibration coefficient array.
    sensor_depth : float
        In situ depth [m].
    salinity : float
        Estimated practical salinity [PSU] used in the speciation
        correction; the DPS-era code comment notes an assumed
        value of 35 where a measured value is unavailable.

    Returns
    -------
    totlgas_bkgco2con : float
        TOTLGAS-BKGCO2CON, total dissolved carbon dioxide
        concentration [uM].

    Notes
    -----
    Speciation factor alpha is 1 + (K1/10^-pH) + (K1*K2/10^-pH^2),
    """

    ph_temp_array = calc_l2_mswater_bkgphval(port_timestamp_bkgndint,
                                             L0_dissgas_bkgndint,
                                             gas_mode_bkgndint,
                                             port_timestamp_bkgndint_mcu,
                                             ph_meter_bkgndint_mcu,
                                             inlet_temp_bkgndint_mcu,
                                             massp_rga_initial_mass,
                                             massp_rga_final_mass,
                                             massp_rga_steps_per_amu,
                                             calibration_table,
                                             l2_ph_calibration_table)

    t = ph_temp_array[1]
    ph = ph_temp_array[0]

    bkgco2con = calc_dissgas_bkgco2con(port_timestamp_bkgndint,
                                       L0_dissgas_bkgndint,
                                       gas_mode_bkgndint,
                                       port_timestamp_bkgndint_mcu,
                                       ph_meter_bkgndint_mcu,
                                       inlet_temp_bkgndint_mcu,
                                       massp_rga_initial_mass,
                                       massp_rga_final_mass,
                                       massp_rga_steps_per_amu,
                                       calibration_table, sensor_depth)

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #estimated salinity == 35
    PSU = salinity

    K1T = np.exp((2.83655 - 2307.1266 / (t + 273.15) - 1.5529413 * np.log(t + 273.15) -
                  (0.20760841 * 4.0484 / (t + 273.15)) * np.sqrt(PSU) + 0.0846834 *
                  PSU - 0.00654208 * np.sqrt(PSU**3) + np.log(1 - 0.001005 * PSU)))

    K2T = np.exp((-9.226508 - 3351.616 / (t + 273.15) - 0.2005743 * np.log(t + 273.15) -
                  (0.106901773 * 23.9722 / (t + 273.15)) * np.sqrt(PSU) + 0.1130822 *
                  PSU - 0.00846934 * np.sqrt(PSU**3) + np.log(1 - 0.001005 * PSU)))

    r1 = (pressure * (1.758163 - 0.008763 * t - pressure * ((7.32 * 10**-6) -
         (2.0845 * 10**-7 * t))) / ((t + 273.15) * 83.131))

    r2 = (pressure * (1.09075 + 0.00151 * t + pressure * ((2.69 * 10**-6) -
         (3.506 * 10**-7 * t))) / ((t + 273.15) * 83.131))

    K1 = np.exp(r1) * K1T

    K2 = np.exp(r2) * K2T

    alpha = 1 + (K1 / (10**-ph)) + ((K1 * K2) / (10**-ph)**2)

    totlgas_bkgco2con = alpha * bkgco2con

    return totlgas_bkgco2con

@deprecated
def calc_timestamp_totlgas_bkgh2scon(port_timestamp_bkgndint,
                                     L0_dissgas_bkgndint,
                                     gas_mode_bkgndint,
                                     port_timestamp_bkgndint_mcu,
                                     ph_meter_bkgndint_mcu,
                                     inlet_temp_bkgndint_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    """
    OOI wrapper for TSTAMP-BKGH2SCON. Returns the Direct mode scan timestamp
    associated with the L2 total dissolved hydrogen sulfide concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the mean port
        timestamp over the Direct mode averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp

@deprecated
def calc_timestamp_totlgas_bkgco2con(port_timestamp_bkgndint,
                                     L0_dissgas_bkgndint,
                                     gas_mode_bkgndint,
                                     port_timestamp_bkgndint_mcu,
                                     ph_meter_bkgndint_mcu,
                                     inlet_temp_bkgndint_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    """
    OOI wrapper for TSTAMP-BKGCO2CON. Returns the Direct mode scan timestamp
    associated with the L2 total dissolved carbon dioxide concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the mean port
        timestamp over the Direct mode averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp


#Block of wrapper functions for calculating the pH intensity auxiliary data products and associated timestamps
@deprecated
def calc_l2_mswater_smpphval(port_timestamp_sampleint, L0_dissgas_sampleint,
                             gas_mode_sampleint, port_timestamp_sampleint_mcu,
                             ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass,
                             massp_rga_steps_per_amu, calibration_table,
                             l2_ph_calibration_table):
    """
    Converts the auxiliary products MSINLET-TEMP and MSINLET-
    SMPPHINT into the higher-level auxiliary product
    MSWATER-SMPPHVAL, the mass-spectrometer-equilibrated
    sample water pH
    value.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 timing and intensity inputs; see SamplePreProcess for the
        full parameter list, which this function passes through
        unchanged.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.
    l2_ph_calibration_table : ndarray
        Six-element L2 pH calibration coefficient array (A0, A1,
        A2, a1, a0, a2, in that storage order).

    Returns
    -------
    l2_smpphint : float
        MSWATER-SMPPHVAL, equilibrated water pH
        [dimensionless]; -9999999.0 if the computed value falls
        outside the valid pH 2-12 range.
    msinlet_temp : float
        MSINLET-TEMP [deg_C] used in the pH calculation, returned
        for reuse by calc_l2_totlgas_* functions.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #msinlet_temp is the 12th element of the preprocess_array array
    msinlet_temp = preprocess_array[11]
    #msinlet_smpphint is the 13th element of the preprocess_array array
    msinlet_smpphint = preprocess_array[12]

    A0 = l2_ph_calibration_table[0]
    A1 = l2_ph_calibration_table[1]
    A2 = l2_ph_calibration_table[2]
    a0 = l2_ph_calibration_table[4]
    a1 = l2_ph_calibration_table[3]
    a2 = l2_ph_calibration_table[5]

    pH = (A0 + (A1 * msinlet_temp) + (A2 * msinlet_temp**2)) * ((a2 * msinlet_smpphint**2) + (a1 * msinlet_smpphint) + a0 + 7)

    if pH < 2 or pH > 12:
        l2_msinlet_smpphint = -9999999.0
    else:
        l2_msinlet_smpphint = pH

    return l2_msinlet_smpphint, msinlet_temp

@deprecated
def calc_l2_mswater_bkgphval(port_timestamp_bkgndint, L0_dissgas_bkgndint,
                             gas_mode_bkgndint, port_timestamp_bkgndint_mcu,
                             ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass,
                             massp_rga_steps_per_amu, calibration_table,
                             l2_ph_calibration_table):
    """
    Converts the auxiliary products MSINLET-TEMP and MSINLET-
    BKGPHINT into the higher-level auxiliary product
    MSWATER-BKGPHVAL, the mass-spectrometer-equilibrated
    background water pH
    value.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 timing and intensity inputs; see BackgroundPreProcess for the
        full parameter list, which this function passes through
        unchanged.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.
    l2_ph_calibration_table : ndarray
        Six-element L2 pH calibration coefficient array (A0, A1,
        A2, a1, a0, a2, in that storage order).

    Returns
    -------
    l2_bkgphint : float
        MSWATER-BKGPHVAL, equilibrated water pH
        [dimensionless]; -9999999.0 if the computed value falls
        outside the valid pH 2-12 range.
    msinlet_temp : float
        MSINLET-TEMP [deg_C] used in the pH calculation, returned
        for reuse by calc_l2_totlgas_* functions.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #msinlet_temp is the 9th element of the preprocess_array array
    msinlet_temp = preprocess_array[8]
    #msinlet_bkgphint is the 10th element of the preprocess_array array
    msinlet_bkgphint = preprocess_array[9]

    A0 = l2_ph_calibration_table[0]
    A1 = l2_ph_calibration_table[1]
    A2 = l2_ph_calibration_table[2]
    a0 = l2_ph_calibration_table[4]
    a1 = l2_ph_calibration_table[3]
    a2 = l2_ph_calibration_table[5]

    pH = (A0 + (A1 * msinlet_temp) + (A2 * msinlet_temp**2)) * ((a2 * msinlet_bkgphint**2) + (a1 * msinlet_bkgphint) + a0 + 7)

    if pH < 2 or pH > 12:
        l2_msinlet_bkgphint = -9999999.0
    else:
        l2_msinlet_bkgphint = pH

    return l2_msinlet_bkgphint, msinlet_temp

@deprecated
def calc_msinlet_smpphint(port_timestamp_sampleint, L0_dissgas_sampleint,
                          gas_mode_sampleint, port_timestamp_sampleint_mcu,
                          ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass,
                          massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for MSINLET-SMPPHINT. Returns the Sample Water pH
    signal intensity [dimensionless] at the time of dissolved gas
    measurement, taken over the last minute of Direct mode.

    See Also
    --------
    SamplePreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint,
                                        gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #msinlet_smpphint is the 13th element of the preprocess_array array
    msinlet_smpphint = preprocess_array[12]

    return msinlet_smpphint

@deprecated
def calc_msinlet_smpphint_timestamp(port_timestamp_sampleint,
                                    L0_dissgas_sampleint, gas_mode_sampleint,
                                    port_timestamp_sampleint_mcu,
                                    ph_meter_sampleint_mcu,
                                    inlet_temp_sampleint_mcu,
                                    massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-SMPPHINT. Returns the Direct mode scan
    timestamp associated with MSINLET-SMPPHINT.

    See Also
    --------
    SamplePreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_msinlet_bkgphint(port_timestamp_bkgndint, L0_dissgas_bkgndint,
                          gas_mode_bkgndint, port_timestamp_bkgndint_mcu,
                          ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass,
                          massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for MSINLET-BKGPHINT. Returns the Background Water pH
    signal intensity [dimensionless] at the time of dissolved gas
    measurement, taken over the last minute of Nafion mode.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #msinlet_bkgphint is the 10th element of the preprocess_array array
    msinlet_bkgphint = preprocess_array[9]

    return msinlet_bkgphint

@deprecated
def calc_msinlet_bkgphint_timestamp(port_timestamp_bkgndint,
                                    L0_dissgas_bkgndint, gas_mode_bkgndint,
                                    port_timestamp_bkgndint_mcu,
                                    ph_meter_bkgndint_mcu,
                                    inlet_temp_bkgndint_mcu,
                                    massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-BKGPHINT. Returns the Nafion mode scan
    timestamp associated with MSINLET-BKGPHINT.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint,
                                            L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu,
                                            ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table,
                                            calibration_table)

    #the nafion mode timestamp is the 11th element of the preprocess_array array
    bkg_nafion_timestamp = preprocess_array[10]

    return bkg_nafion_timestamp

@deprecated
def calc_msinlet_cal1phint(port_timestamp_calint01, L0_dissgas_calint01,
                           gas_mode_calint01, port_timestamp_calint01_mcu,
                           ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass,
                           massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for MSINLET-CA1PHINT. Returns the Calibration Solution
    1 pH signal intensity [dimensionless] at the time of dissolved gas
    measurement, taken over the last minute of Nafion mode.

    See Also
    --------
    Cal1PreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01,
                                      L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu,
                                      ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table,
                                      calibration_table)

    #msinlet_cal1phint is the 5th element of the preprocess_array array
    msinlet_cal1phint = preprocess_array[4]

    return msinlet_cal1phint

@deprecated
def calc_msinlet_cal1phint_timestamp(port_timestamp_calint01,
                                     L0_dissgas_calint01,
                                     gas_mode_calint01,
                                     port_timestamp_calint01_mcu,
                                     ph_meter_calint01_mcu,
                                     inlet_temp_calint01_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    """
    OOI wrapper for TSTAMP-CA1PHINT. Returns the Nafion mode scan
    timestamp associated with MSINLET-CA1PHINT.

    See Also
    --------
    Cal1PreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01,
                                      L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu,
                                      ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table,
                                      calibration_table)

    #the nafion mode timestamp is the 6th element of the preprocess_array array
    cal1_nafion_timestamp = preprocess_array[5]

    return cal1_nafion_timestamp

@deprecated
def calc_msinlet_cal2phint(port_timestamp_calint02, L0_dissgas_calint02,
                           gas_mode_calint02, port_timestamp_calint02_mcu,
                           ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass,
                           massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for MSINLET-CA2PHINT. Returns the Calibration Solution
    2 pH signal intensity [dimensionless] at the time of dissolved gas
    measurement, taken over the last minute of Direct mode.

    See Also
    --------
    Cal2PreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #msinlet_cal2phint is the 5th element of the preprocess_array array
    msinlet_cal2phint = preprocess_array[4]

    return msinlet_cal2phint

@deprecated
def calc_msinlet_cal2phint_timestamp(port_timestamp_calint02,
                                     L0_dissgas_calint02,
                                     gas_mode_calint02,
                                     port_timestamp_calint02_mcu,
                                     ph_meter_calint02_mcu,
                                     inlet_temp_calint02_mcu,
                                     massp_rga_initial_mass,
                                     massp_rga_final_mass,
                                     massp_rga_steps_per_amu,
                                     calibration_table):
    """
    OOI wrapper for TSTAMP-CA2PHINT. Returns the Direct mode scan
    timestamp associated with MSINLET-CA2PHINT.

    See Also
    --------
    Cal2PreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02,
                                      L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu,
                                      ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table,
                                      calibration_table)

    #the direct mode timestamp is the 7th element of the preprocess_array array
    cal2_direct_timestamp = preprocess_array[6]

    return cal2_direct_timestamp


#Block of wrapper functions for calculating the nafion drier efficiency auxiliary data product and associated timestamp
@deprecated
def calc_smpnafeff(port_timestamp_sampleint, L0_dissgas_sampleint,
                   gas_mode_sampleint, port_timestamp_sampleint_mcu,
                   ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                   massp_rga_initial_mass, massp_rga_final_mass,
                   massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for NAFEFF. Returns the Nafion Drier Efficiency [%],
    the percentage of the mz 18 (water) signal seen in Nafion mode
    relative to Direct mode.

    See Also
    --------
    SamplePreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #smpnafeff is the 10th element of the preprocess_array array
    smpnafeff = preprocess_array[9]

    return smpnafeff

@deprecated
def calc_smpnafeff_timestamp(port_timestamp_sampleint, L0_dissgas_sampleint,
                             gas_mode_sampleint, port_timestamp_sampleint_mcu,
                             ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass,
                             massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-NAFEFF. Returns the Nafion mode scan
    timestamp associated with NAFEFF.

    See Also
    --------
    SamplePreProcess : Helper; computes this value directly.
    """

    mass_table = rga_status_process(massp_rga_initial_mass,
                                    massp_rga_final_mass,
                                    massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint,
                                        L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu,
                                        ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table,
                                        calibration_table)

    #the nafion mode timestamp is the 14th element of the preprocess_array array
    smp_nafion_timestamp = preprocess_array[13]

    return smp_nafion_timestamp


#Block of wrapper functions for calculating the L1 data products
@deprecated
def calc_dissgas_smpmethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-SMPMETHCON. Returns the in situ dissolved methane
    concentration [uM] in the sample water, computed from Nafion mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    SamplePreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 11th element of the preprocess_array array
    average_temperature = preprocess_array[10]

    smpmethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpmethcon = smpmethcon[0]

    return smpmethcon

@deprecated
def calc_dissgas_smpethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-SMPETHNCON. Returns the in situ dissolved ethane
    concentration [uM] in the sample water, computed from Direct mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    SamplePreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz30 is the 5th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[4]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 4
    last_column = 8

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpethcon = smpethcon[0]

    return smpethcon

@deprecated
def calc_dissgas_smph2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-SMPH2CON. Returns the in situ dissolved hydrogen
    concentration [uM] in the sample water, computed from Direct mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    SamplePreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz2 is the 3rd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[2]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 8
    last_column = 12

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smph2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smph2con = smph2con[0]

    return smph2con

@deprecated
def calc_dissgas_smparcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-SMPARCON. Returns the in situ dissolved argon
    concentration [uM] in the sample water, computed from Direct mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    SamplePreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz40 is the 8th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[7]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 12
    last_column = 16

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smparcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smparcon = smparcon[0]

    return smparcon

@deprecated
def calc_dissgas_smph2scon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-SMPH2SCON. Returns the in situ dissolved hydrogen
    sulfide concentration [uM] in the sample water, computed from Direct mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    SamplePreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz34 is the 7th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[6]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 16
    last_column = 20

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smph2scon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smph2scon = smph2scon[0]

    return smph2scon

@deprecated
def calc_dissgas_smpo2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-SMPO2CON. Returns the in situ dissolved oxygen
    concentration [uM] in the sample water, computed from Direct mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    SamplePreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz32 is the 6th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[5]
    #sample_mz34 is the 7th element of the preprocess_array array
    deconvolution_variable = preprocess_array[6]
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 20
    last_column = 24

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpo2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpo2con = smpo2con[0]

    return smpo2con

@deprecated
def calc_dissgas_smpco2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-SMPCO2CON. Returns the in situ dissolved carbon
    dioxide concentration [uM] in the sample water, computed from Direct mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    SamplePreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz44 is the 9th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[8]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpco2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpco2con = smpco2con[0]

    return smpco2con

@deprecated
def calc_dissgas_bkgmethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-BKGMETHCON. Returns the in situ dissolved methane
    concentration [uM] in the background water, computed from Nafion mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    BackgroundPreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz15 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 8th element of the preprocess_array array
    average_temperature = preprocess_array[7]

    bkgmethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgmethcon = bkgmethcon[0]

    return bkgmethcon

@deprecated
def calc_dissgas_bkgethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-BKGETHNCON. Returns the in situ dissolved ethane
    concentration [uM] in the background water, computed from Direct mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    BackgroundPreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz30 is the 3rd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[2]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 4
    last_column = 8

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgethcon = bkgethcon[0]

    return bkgethcon

@deprecated
def calc_dissgas_bkgh2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-BKGH2CON. Returns the in situ dissolved hydrogen
    concentration [uM] in the background water, computed from Direct mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    BackgroundPreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz2 is the 1st element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 8
    last_column = 12

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgh2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgh2con = bkgh2con[0]

    return bkgh2con

@deprecated
def calc_dissgas_bkgarcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-BKGARCON. Returns the in situ dissolved argon
    concentration [uM] in the background water, computed from Direct mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    BackgroundPreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz40 is the 6th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[5]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 12
    last_column = 16

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgarcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgarcon = bkgarcon[0]

    return bkgarcon

@deprecated
def calc_dissgas_bkgh2scon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-BKGH2SCON. Returns the in situ dissolved hydrogen
    sulfide concentration [uM] in the background water, computed from Direct
    mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    BackgroundPreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz34 is the 5th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[4]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 16
    last_column = 20

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgh2scon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgh2scon = bkgh2scon[0]

    return bkgh2scon

@deprecated
def calc_dissgas_bkgo2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-BKGO2CON. Returns the in situ dissolved oxygen
    concentration [uM] in the background water, computed from Direct mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    BackgroundPreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz32 is the 4th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[3]
    #sample_mz34 is the 5th element of the preprocess_array array
    deconvolution_variable = preprocess_array[4]
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 20
    last_column = 24

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgo2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgo2con = bkgo2con[0]

    return bkgo2con

@deprecated
def calc_dissgas_bkgco2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-BKGCO2CON. Returns the in situ dissolved carbon
    dioxide concentration [uM] in the background water, computed from Direct
    mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    BackgroundPreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz44 is the 7th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[6]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgco2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgco2con = bkgco2con[0]

    return bkgco2con

@deprecated
def calc_dissgas_cal1methcon(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                             port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                             calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-CA1METHCON. Returns the in situ dissolved methane
    concentration [uM] in Calibration Solution 1, computed from Nafion mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    Cal1PreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #cal1_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 3rd element of the preprocess_array array
    average_temperature = preprocess_array[2]

    cal1methcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                    first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    cal1methcon = cal1methcon[0]

    return cal1methcon

@deprecated
def calc_dissgas_cal1co2con(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                            port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-CA1CO2CON. Returns the in situ dissolved carbon
    dioxide concentration [uM] in Calibration Solution 1, computed from Direct
    mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    Cal1PreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #cal1_mz44 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 4th element of the preprocess_array array
    average_temperature = preprocess_array[3]

    cal1co2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    cal1co2con = cal1co2con[0]

    return cal1co2con

@deprecated
def calc_dissgas_cal2methcon(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                             port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                             calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-CA2METHCON. Returns the in situ dissolved methane
    concentration [uM] in Calibration Solution 2, computed from Nafion mode
    scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    Cal2PreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #Cal2_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 3rd element of the preprocess_array array
    average_temperature = preprocess_array[2]

    cal2methcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                    first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    cal2methcon = cal2methcon[0]

    return cal2methcon

@deprecated
def calc_dissgas_cal2co2con(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                            port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    """
    OOI wrapper for DISSGAS-CA2CO2CON. Returns the in situ dissolved carbon
    dioxide concentration [uM] in Calibration Solution 2, computed from Direct
    mode scans.

    See Also
    --------
    gas_concentration : Core algorithm; use directly for the
        corrected-intensity-to-concentration calculation and the
        associated CALRANG quality flag.
    Cal2PreProcess : Helper; produces the averaged mz
        intensity and mode-averaged temperature consumed here.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #Cal2_mz44 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 4th element of the preprocess_array array
    average_temperature = preprocess_array[3]

    cal2co2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    cal2co2con = cal2co2con[0]

    return cal2co2con


#Block of wrapper functions for calculating the timestamps of the L1 data products
@deprecated
def calc_timestamp_smpmethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                              port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-SMPMETH. Returns the Nafion mode scan
    timestamp associated with the dissolved methane concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the nafion mode timestamp is the 14th element of the preprocess_array array
    smp_nafion_timestamp = preprocess_array[13]

    return smp_nafion_timestamp

@deprecated
def calc_timestamp_smpethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                             port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-SMPETHN. Returns the Direct mode scan
    timestamp associated with the dissolved ethane concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_timestamp_smph2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-SMPH2. Returns the Direct mode scan
    timestamp associated with the dissolved hydrogen concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_timestamp_smparcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-SMPAR. Returns the Direct mode scan
    timestamp associated with the dissolved argon concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_timestamp_smph2scon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                             port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-SMPH2S. Returns the Direct mode scan
    timestamp associated with the dissolved hydrogen sulfide concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_timestamp_smpo2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-SMPO2. Returns the Direct mode scan
    timestamp associated with the dissolved oxygen concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_timestamp_smpco2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                             port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-SMPCO2. Returns the Direct mode scan
    timestamp associated with the dissolved carbon dioxide concentration
    in the sample water.

    See Also
    --------
    SamplePreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 15th element of the preprocess_array array
    smp_direct_timestamp = preprocess_array[14]

    return smp_direct_timestamp

@deprecated
def calc_timestamp_bkgmethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                              port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-BKGMETH. Returns the Nafion mode scan
    timestamp associated with the dissolved methane concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the nafion mode timestamp is the 11th element of the preprocess_array array
    bkg_nafion_timestamp = preprocess_array[10]

    return bkg_nafion_timestamp

@deprecated
def calc_timestamp_bkgethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                             port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-BKGETHN. Returns the Direct mode scan
    timestamp associated with the dissolved ethane concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp

@deprecated
def calc_timestamp_bkgh2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-BKGH2. Returns the Direct mode scan
    timestamp associated with the dissolved hydrogen concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp

@deprecated
def calc_timestamp_bkgarcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-BKGAR. Returns the Direct mode scan
    timestamp associated with the dissolved argon concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp

@deprecated
def calc_timestamp_bkgh2scon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                             port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-BKGH2S. Returns the Direct mode scan
    timestamp associated with the dissolved hydrogen sulfide concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp

@deprecated
def calc_timestamp_bkgo2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-BKGO2. Returns the Direct mode scan
    timestamp associated with the dissolved oxygen concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp

@deprecated
def calc_timestamp_bkgco2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                             port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-BKGCO2. Returns the Direct mode scan
    timestamp associated with the dissolved carbon dioxide concentration
    in the background water.

    See Also
    --------
    BackgroundPreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 12th element of the preprocess_array array
    bkg_direct_timestamp = preprocess_array[11]

    return bkg_direct_timestamp

@deprecated
def calc_timestamp_cal1methcon(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                               port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                               massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-CA1METH. Returns the Nafion mode scan
    timestamp associated with the dissolved methane concentration
    in Calibration Solution 1.

    See Also
    --------
    Cal1PreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #the nafion mode timestamp is the 6th element of the preprocess_array array
    cal1_nafion_timestamp = preprocess_array[5]

    return cal1_nafion_timestamp

@deprecated
def calc_timestamp_cal1co2con(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                              port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-CA1CO2. Returns the Direct mode scan
    timestamp associated with the dissolved carbon dioxide concentration
    in Calibration Solution 1.

    See Also
    --------
    Cal1PreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 7th element of the preprocess_array array
    cal1_direct_timestamp = preprocess_array[6]

    return cal1_direct_timestamp

@deprecated
def calc_timestamp_cal2methcon(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                               port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                               massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-CA2METH. Returns the Nafion mode scan
    timestamp associated with the dissolved methane concentration
    in Calibration Solution 2.

    See Also
    --------
    Cal2PreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #the nafion mode timestamp is the 6th element of the preprocess_array array
    cal2_nafion_timestamp = preprocess_array[5]

    return cal2_nafion_timestamp

@deprecated
def calc_timestamp_cal2co2con(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                              port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                              massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu, calibration_table):
    """
    OOI wrapper for TSTAMP-CA2CO2. Returns the Direct mode scan
    timestamp associated with the dissolved carbon dioxide concentration
    in Calibration Solution 2.

    See Also
    --------
    Cal2PreProcess : Helper; computes this timestamp as the
        mean port timestamp over the averaging window.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #the direct mode timestamp is the 7th element of the preprocess_array array
    cal2_direct_timestamp = preprocess_array[6]

    return cal2_direct_timestamp


#Block of wrapper functions for calculating the calibration ranges of the L1 data products
@deprecated
def calc_calrang_smpmethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                            port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-SMPMETH. Returns the MASSP Calibration
    Range quality flag for the dissolved methane concentration in
    the sample water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 11th element of the preprocess_array array
    average_temperature = preprocess_array[10]

    smpmethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpmethcon = smpmethcon[1]

    return smpmethcon

@deprecated
def calc_calrang_smpethcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-SMPETHN. Returns the MASSP Calibration
    Range quality flag for the dissolved ethane concentration in
    the sample water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz30 is the 5th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[4]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 4
    last_column = 8

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpethcon = smpethcon[1]

    return smpethcon

@deprecated
def calc_calrang_smph2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-SMPH2. Returns the MASSP Calibration
    Range quality flag for the dissolved hydrogen concentration in
    the sample water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz2 is the 3rd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[2]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 8
    last_column = 12

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smph2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smph2con = smph2con[1]

    return smph2con

@deprecated
def calc_calrang_smparcon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-SMPAR. Returns the MASSP Calibration
    Range quality flag for the dissolved argon concentration in
    the sample water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz40 is the 8th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[7]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 12
    last_column = 16

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smparcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smparcon = smparcon[1]

    return smparcon

@deprecated
def calc_calrang_smph2scon(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-SMPH2S. Returns the MASSP Calibration
    Range quality flag for the dissolved hydrogen sulfide concentration in
    the sample water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz34 is the 7th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[6]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 16
    last_column = 20

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smph2scon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smph2scon = smph2scon[1]

    return smph2scon

@deprecated
def calc_calrang_smpo2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                          port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-SMPO2. Returns the MASSP Calibration
    Range quality flag for the dissolved oxygen concentration in
    the sample water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz32 is the 6th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[5]
    #sample_mz34 is the 7th element of the preprocess_array array
    deconvolution_variable = preprocess_array[6]
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 20
    last_column = 24

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpo2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpo2con = smpo2con[1]

    return smpo2con

@deprecated
def calc_calrang_smpco2con(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                           port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-SMPCO2. Returns the MASSP Calibration
    Range quality flag for the dissolved carbon dioxide concentration in
    the sample water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                                        port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu,
                                        inlet_temp_sampleint_mcu, mass_table, calibration_table)

    #sample_mz44 is the 9th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[8]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 12th element of the preprocess_array array
    average_temperature = preprocess_array[11]

    smpco2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    smpco2con = smpco2con[1]

    return smpco2con

@deprecated
def calc_calrang_bkgmethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-BKGMETH. Returns the MASSP Calibration
    Range quality flag for the dissolved methane concentration in
    the background water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz15 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 8th element of the preprocess_array array
    average_temperature = preprocess_array[7]

    bkgmethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgmethcon = bkgmethcon[1]

    return bkgmethcon

@deprecated
def calc_calrang_bkgethcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-BKGETHN. Returns the MASSP Calibration
    Range quality flag for the dissolved ethane concentration in
    the background water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz30 is the 3rd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[2]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 4
    last_column = 8

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgethcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgethcon = bkgethcon[1]

    return bkgethcon

@deprecated
def calc_calrang_bkgh2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-BKGH2. Returns the MASSP Calibration
    Range quality flag for the dissolved hydrogen concentration in
    the background water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz2 is the 1st element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 8
    last_column = 12

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgh2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgh2con = bkgh2con[1]

    return bkgh2con

@deprecated
def calc_calrang_bkgarcon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-BKGAR. Returns the MASSP Calibration
    Range quality flag for the dissolved argon concentration in
    the background water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz40 is the 6th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[5]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 12
    last_column = 16

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgarcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgarcon = bkgarcon[1]

    return bkgarcon

@deprecated
def calc_calrang_bkgh2scon(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-BKGH2S. Returns the MASSP Calibration
    Range quality flag for the dissolved hydrogen sulfide concentration in
    the background water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz34 is the 5th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[4]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 16
    last_column = 20

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgh2scon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgh2scon = bkgh2scon[1]

    return bkgh2scon

@deprecated
def calc_calrang_bkgo2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                          port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                          massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                          calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-BKGO2. Returns the MASSP Calibration
    Range quality flag for the dissolved oxygen concentration in
    the background water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz32 is the 4th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[3]
    #sample_mz34 is the 5th element of the preprocess_array array
    deconvolution_variable = preprocess_array[4]
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 20
    last_column = 24

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgo2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                 first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgo2con = bkgo2con[1]

    return bkgo2con

@deprecated
def calc_calrang_bkgco2con(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                           port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                           massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                           calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-BKGCO2. Returns the MASSP Calibration
    Range quality flag for the dissolved carbon dioxide concentration in
    the background water. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                                            port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu,
                                            inlet_temp_bkgndint_mcu, mass_table, calibration_table)

    #sample_mz44 is the 7th element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[6]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 9th element of the preprocess_array array
    average_temperature = preprocess_array[8]

    bkgco2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    bkgco2con = bkgco2con[1]

    return bkgco2con

@deprecated
def calc_calrang_cal1methcon(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                             port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                             calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-CA1METH. Returns the MASSP Calibration
    Range quality flag for the dissolved methane concentration in
    Calibration Solution 1. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #cal1_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 3rd element of the preprocess_array array
    average_temperature = preprocess_array[2]

    ca1methcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    ca1methcon = ca1methcon[1]

    return ca1methcon

@deprecated
def calc_calrang_cal1co2con(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                            port_timestamp_calint01_mcu, ph_meter_calint01_mcu, inlet_temp_calint01_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-CA1CO2. Returns the MASSP Calibration
    Range quality flag for the dissolved carbon dioxide concentration in
    Calibration Solution 1. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                                      port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                                      inlet_temp_calint01_mcu, mass_table, calibration_table)

    #cal1_mz44 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 4th element of the preprocess_array array
    average_temperature = preprocess_array[3]

    ca1co2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    ca1co2con = ca1co2con[1]

    return ca1co2con

@deprecated
def calc_calrang_cal2methcon(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                             port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                             massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                             calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-CA2METH. Returns the MASSP Calibration
    Range quality flag for the dissolved methane concentration in
    Calibration Solution 2. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #Cal2_mz15 is the first element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[0]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 0
    last_column = 4

    #average inlet temperature (nafion mode) is the 3rd element of the preprocess_array array
    average_temperature = preprocess_array[2]

    ca2methcon = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                   first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    ca2methcon = ca2methcon[1]

    return ca2methcon

@deprecated
def calc_calrang_cal2co2con(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                            port_timestamp_calint02_mcu, ph_meter_calint02_mcu, inlet_temp_calint02_mcu,
                            massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu,
                            calibration_table, sensor_depth):
    """
    OOI wrapper for CALRANG-CA2CO2. Returns the MASSP Calibration
    Range quality flag for the dissolved carbon dioxide concentration in
    Calibration Solution 2. Values: -1 (intensity below calibration range),
    0 (intensity and temperature within range), 1 (intensity above
    range, temperature within range), 2 (intensity within range,
    temperature above range), 3 (intensity and temperature above
    range).

    See Also
    --------
    gas_concentration : Core algorithm; returns this flag alongside
        the concentration value.
    """

    mass_table = rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu)

    preprocess_array = Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                                      port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                                      inlet_temp_calint02_mcu, mass_table, calibration_table)

    #Cal2_mz44 is the 2nd element of the preprocess_array array
    intermediate_mass_ratio = preprocess_array[1]
    deconvolution_variable = 0
    #first and last column for this particular gas in the calibration table
    #This information is in table 1 of the DPS
    first_column = 24
    last_column = 28

    #average inlet temperature (direct mode) is the 4th element of the preprocess_array array
    average_temperature = preprocess_array[3]

    ca2co2con = gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                                  first_column, last_column, sensor_depth, average_temperature)

    #the first element of the array if the gas conc. the second is the calrange
    ca2co2con = ca2co2con[1]

    return ca2co2con


#Block of subfunctions called by the above wrapper functions that calculate the L1 and auxiliary data products
@deprecated
def gas_concentration(intermediate_mass_ratio, deconvolution_variable, calibration_table,
                      first_column, last_column, sensor_depth, average_temperature):
    """
    Converts a corrected mass spectral intensity into a final dissolved
    gas concentration and its MASSP Calibration Range (CALRANG) quality
    flag. This is the shared core algorithm called by every DISSGAS and
    CALRANG wrapper function, each of which supplies a gas-specific
    column range into the calibration table.

    Parameters
    ----------
    intermediate_mass_ratio : float
        Averaged mz intensity for the target gas, from average_mz.
    deconvolution_variable : float
        Secondary mz intensity used only for the oxygen deconvolution
        correction (mz 34); 0 for all other gases.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.
    first_column : int
        First column of the gas-specific range within calibration_table.
    last_column : int
        Last column (exclusive) of the gas-specific range within
        calibration_table.
    sensor_depth : float
        In situ depth [m]; converted internally to pressure [psi].
    average_temperature : float
        Mode-averaged MSINLET-TEMP [deg_C] for the gas and water source
        being processed.

    Returns
    -------
    final_conc : float
        Dissolved gas concentration, DISSGAS core data product [uM].
    calrange : int
        MASSP Calibration Range quality flag (-1, 0, 1, 2, or 3),
        CALRANG auxiliary data product.

    Notes
    -----
    Selects one or two calibration temperature columns bracketing
    average_temperature, evaluates a pressure-dependent polynomial and
    exponential fit at each, and interpolates between them when two
    columns are used.
    """

    #Converth depth (meters) to pressure (psi)
    pressure = (sensor_depth * 0.099204 + 1) * 14.695

    #extract the four columns of the cal table that I need for a particular gas.
    calibration_table = calibration_table[:, first_column:last_column]

    #Check to see if one of the 4 calibration temperatures == the averaged inlet temperature
    ind = np.where(calibration_table[0, :] == average_temperature)[0]
    if np.size(ind) == 1:
        ct1 = np.where(calibration_table[0, :] == average_temperature)[0][0]
        tempCalRange = 1
    #Check to see if the averaged inlet temperature is greater than the highest calibration temperatures
    elif average_temperature >= calibration_table[0, 3]:
        ct1 = 3
        tempCalRange = 2
    #Otherwise figure out which two columns in the calibration table are needed.
    else:
        #find first column
        ct1 = np.where(calibration_table[0, :] < average_temperature)[0][-1]
        #find second column
        ct2 = np.where(calibration_table[0, :] > average_temperature)[0][0]
        concT2_flag = 1
        tempCalRange = 0

    corrected_intensity = deconvolution_correction(intermediate_mass_ratio, deconvolution_variable, calibration_table)

    #Check to see if the corrected intensity falls within the calibration values
    #Minimum values, row 4 in the cal table
    if corrected_intensity < calibration_table[3, ct1] or corrected_intensity < calibration_table[3, ct2]:
        calrange = -1
    #Maximum values, row 5 in the cal table
    elif corrected_intensity > calibration_table[4, ct1] or corrected_intensity > calibration_table[4, ct2]:
        calrange = tempCalRange + 1
    else:
        calrange = tempCalRange

    #P0 is row 6 (with row 0 being the first row) in the cal table
    if corrected_intensity < calibration_table[5, ct1]:
        alpha = calibration_table[6, ct1] + (calibration_table[7, ct1] * pressure) + (calibration_table[8, ct1] * pressure**2) + (calibration_table[9, ct1] * pressure**3)
        beta = calibration_table[14, ct1] + (calibration_table[15, ct1] * pressure) + (calibration_table[16, ct1] * pressure**2) + (calibration_table[17, ct1] * pressure**3)
        delta = calibration_table[22, ct1] + (calibration_table[23, ct1] * pressure) + (calibration_table[24, ct1] * pressure**2) + (calibration_table[25, ct1] * pressure**3)
        gamma = calibration_table[28, ct1] + (calibration_table[29, ct1] * pressure) + (calibration_table[30, ct1] * pressure**2) + (calibration_table[31, ct1] * pressure**3)
        zeta = calibration_table[36, ct1] + (calibration_table[37, ct1] * pressure) + (calibration_table[38, ct1] * pressure**2) + (calibration_table[39, ct1] * pressure**3)
    elif corrected_intensity >= calibration_table[5, ct1]:
        alpha = calibration_table[10, ct1] + (calibration_table[11, ct1] * pressure) + (calibration_table[12, ct1] * pressure**2) + (calibration_table[13, ct1] * pressure**3)
        beta = calibration_table[18, ct1] + (calibration_table[19, ct1] * pressure) + (calibration_table[20, ct1] * pressure**2) + (calibration_table[21, ct1] * pressure**3)
        delta = calibration_table[26, ct1] * np.exp(calibration_table[26, ct1] * pressure)
        gamma = calibration_table[32, ct1] + (calibration_table[33, ct1] * pressure) + (calibration_table[34, ct1] * pressure**2) + (calibration_table[35, ct1] * pressure**3)
        zeta = calibration_table[40, ct1] + (calibration_table[41, ct1] * pressure) + (calibration_table[42, ct1] * pressure**2) + (calibration_table[43, ct1] * pressure**3)

    #Calculate concT1
    concT1 = (alpha * corrected_intensity**2) + (beta * corrected_intensity) + (delta * np.exp(zeta * corrected_intensity)) + gamma

    if concT2_flag == 1:
        if corrected_intensity < calibration_table[5, ct2]:
            alpha = calibration_table[6, ct2] + (calibration_table[7, ct2] * pressure) + (calibration_table[8, ct2] * pressure**2) + (calibration_table[9, ct2] * pressure**3)
            beta = calibration_table[14, ct2] + (calibration_table[15, ct2] * pressure) + (calibration_table[16, ct2] * pressure**2) + (calibration_table[17, ct2] * pressure**3)
            delta = calibration_table[22, ct2] + (calibration_table[23, ct2] * pressure) + (calibration_table[24, ct2] * pressure**2) + (calibration_table[25, ct2] * pressure**3)
            gamma = calibration_table[28, ct2] + (calibration_table[29, ct2] * pressure) + (calibration_table[30, ct2] * pressure**2) + (calibration_table[31, ct2] * pressure**3)
            zeta = calibration_table[36, ct2] + (calibration_table[37, ct2] * pressure) + (calibration_table[38, ct2] * pressure**2) + (calibration_table[39, ct2] * pressure**3)
        elif corrected_intensity >= calibration_table[5, ct2]:
            alpha = calibration_table[10, ct2] + (calibration_table[11, ct2] * pressure) + (calibration_table[12, ct2] * pressure**2) + (calibration_table[13, ct2] * pressure**3)
            beta = calibration_table[18, ct2] + (calibration_table[19, ct2] * pressure) + (calibration_table[20, ct2] * pressure**2) + (calibration_table[21, ct2] * pressure**3)
            delta = calibration_table[26, ct2] * np.exp(calibration_table[26, ct2] * pressure)
            gamma = calibration_table[32, ct2] + (calibration_table[33, ct2] * pressure) + (calibration_table[34, ct2] * pressure**2) + (calibration_table[35, ct2] * pressure**3)
            zeta = calibration_table[40, ct2] + (calibration_table[41, ct2] * pressure) + (calibration_table[42, ct2] * pressure**2) + (calibration_table[43, ct2] * pressure**3)
        #Calculate concT2
        concT2 = (alpha * corrected_intensity**2) + (beta * corrected_intensity) + (delta * np.exp(zeta * corrected_intensity)) + gamma
        #Calculate concT
        concT = concT1 + ((concT2 - concT1) * (average_temperature - calibration_table[0, ct1])) / (calibration_table[0, ct2] - calibration_table[0, ct1])
    else:
        #Calculate concT
        concT = concT1

    if calrange == -1:
        final_conc = 0
    else:
        final_conc = calibration_table[44, ct1] * (concT - calibration_table[45, ct1])

    return final_conc, calrange

def average_mz(mz, data_in, mass_table, window):
    """
    Averages a subset of RGA scans at a given mass-to-charge ratio to
    produce an intermediate mz intensity.

    Parameters
    ----------
    mz : int
        Target mass-to-charge ratio.
    data_in : ndarray
        Subset of L0 mass spectral scans (sample, background, or
        calibration fluid) for the time window being averaged.
    mass_table : ndarray
        RGA mass lookup table from rga_status_process.
    window : float
        Half-width of the mz +/- window searched in mass_table; taken
        from the calibration table's last row for the relevant gas.

    Returns
    -------
    intermediate_mass_ratio : float
        Mean, across scans, of the per-scan intensity at mz; negative
        values are set to zero before averaging.

    Notes
    -----
    The code takes the second-highest intensity within the mz window
    for each scan, not a median of three values as the original code
    comment states; see Additional Notes.
    """
    #find mz +/- window in the mass_table. The window value comes from the L1 Cal Table
    mz_ind = np.where((mass_table >= mz - window) & (mass_table <= mz + window))

    #subset the data_in array so that we are just dealing with the mz values
    #within the time period of interest
    temp_array = np.array(data_in[:, mz_ind])
    temp_array = np.squeeze(temp_array)

    #sort the array so that I can find the median of the three highest
    #values for each scan
    temp_array = np.sort(temp_array)

    #grab the median values
    median_array = temp_array[:, -2]
    #find and replace any negative values with zero
    median_ind = np.where(median_array < 0)
    median_array[median_ind] = 0
    #calculate the mean of the median values
    intermediate_mass_ratio = np.nanmean(median_array)

    return intermediate_mass_ratio

def deconvolution_correction(intermediate_mass_ratio, deconvolution_variable, calibration_table):
    """
    Applies the intensity deconvolution correction (DPS Equation 4) to
    an averaged mz intensity, using a second mz intensity and two
    calibration table coefficients.

    Parameters
    ----------
    intermediate_mass_ratio : float
        Averaged mz intensity for the target gas, from average_mz.
    deconvolution_variable : float
        Secondary mz intensity used in the correction; 0 where no
        deconvolution is required for the gas being processed.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table, already sliced to
        the gas-specific column range.

    Returns
    -------
    corrected_intensity : float
        Deconvolution-corrected intensity, consumed by
        gas_concentration.
    """
    #Equ 4 on page 13 of the DPS
    corrected_intensity = intermediate_mass_ratio - (calibration_table[2, 0] * deconvolution_variable) - calibration_table[1, 0]

    return corrected_intensity

def rga_status_process(massp_rga_initial_mass, massp_rga_final_mass, massp_rga_steps_per_amu):
    """
    Builds the RGA mass-to-charge lookup table from the residual gas
    analyzer scan configuration.

    Parameters
    ----------
    massp_rga_initial_mass : float
        RGA scan starting mass [amu].
    massp_rga_final_mass : float
        RGA scan ending mass [amu].
    massp_rga_steps_per_amu : float
        RGA steps per unit mass [amu^-1].

    Returns
    -------
    mass_table : ndarray
        Mass-to-charge value at each RGA scan step, rounded to one
        decimal place.
    """

    Tnb = int(((massp_rga_final_mass - massp_rga_initial_mass) * massp_rga_steps_per_amu) + 1)

    mass_table = np.ones(Tnb)
    mass_table[0] = massp_rga_initial_mass

    for x in range(1, Tnb):
        mass_table[x] = mass_table[x-1] + (1 / float(massp_rga_steps_per_amu))

    mass_table = np.around(mass_table, decimals=1)

    return mass_table

def SamplePreProcess(port_timestamp_sampleint, L0_dissgas_sampleint, gas_mode_sampleint,
                     port_timestamp_sampleint_mcu, ph_meter_sampleint_mcu, inlet_temp_sampleint_mcu,
                     mass_table, calibration_table):
    """
    Groups sample water scans into Direct and Nafion mode
    subsets and extracts the averaged mz intensities, mode-
    averaged temperatures, and pH intensity needed for the L1
    DISSGAS, CALRANG, TSTAMP, MSINLET, and (for the sample water
    only) NAFEFF auxiliary data products.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 port timestamps for the mass spectral data.
    L0_dissgas... : ndarray
        L0 mass spectral intensity data set [A].
    gas_mode... : ndarray
        GASMODE auxiliary data product for this water source.
    port_timestamp..._mcu : ndarray
        MCU timestamps for the pH and temperature sensors.
    ph_meter..._mcu : ndarray
        Raw pH signal intensity [dimensionless].
    inlet_temp..._mcu : ndarray
        Sample or background temperature [deg_C]; -127 and 85 are
        replaced with NaN as instrument not-sampling fill values.
    mass_table : ndarray
        RGA mass lookup table from rga_status_process.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.

    Returns
    -------
    tuple of float
    sample_mz15, sample_mz18naf, sample_mz2, sample_mz18, sample_mz30,
    sample_mz32, sample_mz34, sample_mz40, sample_mz44, nafeff, sample_Tnaf,
    sample_Tdir, msinlet_smpphint, nafion_mode_timestamp, and
    direct_mode_timestamp, in that order.
    """

    #replace bad data with nans
    inlet_temp_sampleint_mcu[inlet_temp_sampleint_mcu == -127] = np.nan
    inlet_temp_sampleint_mcu[inlet_temp_sampleint_mcu == 85] = np.nan

    #find gas_mode_sampleint == 0 (direct mode)
    ind_direct = np.where(gas_mode_sampleint == 0)[0]
    Tchange = port_timestamp_sampleint_mcu[ind_direct[0]]
    Tlast = port_timestamp_sampleint_mcu[ind_direct[-1]]

    #find gas_mode_sampleint == 1 (nafion mode)
    ind_nafion = np.where(gas_mode_sampleint == 1)[0]
    TlastScanNafion = port_timestamp_sampleint_mcu[ind_nafion[-1]]

    #ID timestamp closest to TlastScanNafion - 180
    idx = (np.abs(port_timestamp_sampleint-(TlastScanNafion - 180))).argmin()

    #subset the data collected in nafion mode
    nafion_samples_ind = np.where((port_timestamp_sampleint >= port_timestamp_sampleint[idx]) & (port_timestamp_sampleint <= TlastScanNafion))
    nafion_samples = np.squeeze(np.array(L0_dissgas_sampleint[nafion_samples_ind, :]))
    #DPS says to exclude the last scan at TlastScanNafion. This subsetted array then gets fed into the ave. routine.
    nafion_samples = nafion_samples[:-1, ]

    #calculate nafion mode timestamp
    nafion_mode_timestamp = np.squeeze(np.array(port_timestamp_sampleint[nafion_samples_ind]))
    #DPS says to exclude the last scan at TlastScanNafion.
    nafion_mode_timestamp = np.around(np.nanmean(nafion_mode_timestamp[:-1, ]))

    mass_charge_ratio = 15
    window = round(calibration_table[-1, 0], 1)
    sample_mz15 = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    mass_charge_ratio = 18
    #not sure that this window of 0.5 is OK but the 18mz window is not specified in the cal table
    window = round(calibration_table[-1, 8], 1)
    sample_mz18naf = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    #average MSINLET-TEMP for nafion time period
    nafion_samples_ind = np.squeeze(np.where((port_timestamp_sampleint_mcu >= port_timestamp_sampleint[idx]) & (port_timestamp_sampleint_mcu <= TlastScanNafion)))
    sample_Tnaf = np.nanmean(inlet_temp_sampleint_mcu[nafion_samples_ind[:-1]])

    #ID timestamp closest to Tlast - 180
    idx = (np.abs(port_timestamp_sampleint-(Tlast - 180))).argmin()

    #subset the data collected in direct mode
    direct_samples_ind = np.where((port_timestamp_sampleint >= port_timestamp_sampleint[idx]) & (port_timestamp_sampleint <= Tlast))
    direct_samples = np.squeeze(np.array(L0_dissgas_sampleint[direct_samples_ind, :]))

    #calculate direct mode timestamp
    direct_mode_timestamp = np.array(port_timestamp_sampleint[direct_samples_ind])
    direct_mode_timestamp = np.around(np.nanmean(np.squeeze(direct_mode_timestamp)))

    mass_charge_ratio = 2
    window = round(calibration_table[-1, 8], 1)
    sample_mz2 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 18
    #not sure that this window is true but the window is not specified in the cal table
    window = round(calibration_table[-1, 8], 1)
    sample_mz18 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 30
    window = round(calibration_table[-1, 4], 1)
    sample_mz30 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 32
    window = round(calibration_table[-1, 20], 1)
    sample_mz32 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 34
    window = round(calibration_table[-1, 16], 1)
    sample_mz34 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 40
    window = round(calibration_table[-1, 12], 1)
    sample_mz40 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 44
    window = round(calibration_table[-1, 24], 1)
    sample_mz44 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    #average MSINLET-TEMP for direct time period here, call it sample-Tdir
    direct_samples_ind = np.where((port_timestamp_sampleint_mcu >= port_timestamp_sampleint[idx]) & (port_timestamp_sampleint_mcu <= Tlast))
    sample_Tdir = np.nanmean(inlet_temp_sampleint_mcu[direct_samples_ind])

    #average ph_meter_value for time period Tlast-1min:Tlast, call it msinlet_smpphint
    direct_samples_ind = np.where((port_timestamp_sampleint_mcu >= Tlast - 60) & (port_timestamp_sampleint_mcu <= Tlast))
    msinlet_smpphint = np.absolute(np.nanmean(ph_meter_sampleint_mcu[direct_samples_ind]))

    #Calculate NAFEFF, which is an indicator of the drying efficiency of the nafion drier
    nafeff = int(100 * (sample_mz18naf / sample_mz18))

    return (sample_mz15, sample_mz18naf, sample_mz2, sample_mz18, sample_mz30,
            sample_mz32, sample_mz34, sample_mz40, sample_mz44, nafeff, sample_Tnaf,
            sample_Tdir, msinlet_smpphint, nafion_mode_timestamp, direct_mode_timestamp)

def BackgroundPreProcess(port_timestamp_bkgndint, L0_dissgas_bkgndint, gas_mode_bkgndint,
                         port_timestamp_bkgndint_mcu, ph_meter_bkgndint_mcu, inlet_temp_bkgndint_mcu,
                         mass_table, calibration_table):
    """
    Groups background water scans into Direct and Nafion mode
    subsets and extracts the averaged mz intensities, mode-
    averaged temperatures, and pH intensity needed for the L1
    DISSGAS, CALRANG, TSTAMP, MSINLET, and (for the sample water
    only) NAFEFF auxiliary data products.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 port timestamps for the mass spectral data.
    L0_dissgas... : ndarray
        L0 mass spectral intensity data set [A].
    gas_mode... : ndarray
        GASMODE auxiliary data product for this water source.
    port_timestamp..._mcu : ndarray
        MCU timestamps for the pH and temperature sensors.
    ph_meter..._mcu : ndarray
        Raw pH signal intensity [dimensionless].
    inlet_temp..._mcu : ndarray
        Sample or background temperature [deg_C]; -127 and 85 are
        replaced with NaN as instrument not-sampling fill values.
    mass_table : ndarray
        RGA mass lookup table from rga_status_process.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.

    Returns
    -------
    tuple of float
    bckgnd_mz2, bckgnd_mz15, bckgnd_mz30, bckgnd_mz32, bckgnd_mz34,
    bckgnd_mz40, bckgnd_mz44, bckgnd_Tnaf, bckgnd_Tdir, msinlet_bkgphint,
    nafion_mode_timestamp, and direct_mode_timestamp, in that order.
    """

    #replace bad data with nans
    inlet_temp_bkgndint_mcu[inlet_temp_bkgndint_mcu == -127] = np.nan
    inlet_temp_bkgndint_mcu[inlet_temp_bkgndint_mcu == 85] = np.nan

    #find gas_mode_bkgndint == 0 (direct mode)
    ind_direct = np.where(gas_mode_bkgndint == 0)[0]
    Tchange = port_timestamp_bkgndint_mcu[ind_direct[0]]
    Tlast = port_timestamp_bkgndint_mcu[ind_direct[-1]]

    #ID timestamp closest to Tlast - 180
    idx = (np.abs(port_timestamp_bkgndint-(Tlast - 180))).argmin()

    #subset the data collected in direct mode
    direct_samples_ind = np.where((port_timestamp_bkgndint >= port_timestamp_bkgndint[idx]) & (port_timestamp_bkgndint <= Tlast))
    direct_samples = np.squeeze(np.array(L0_dissgas_bkgndint[direct_samples_ind, :]))
    #DPS says to exclude the last scan at TlastScanDirect. This subsetted array then gets fed into the ave. routine.
    direct_samples = direct_samples[:-1, ]

    #calculate direct mode timestamp
    direct_mode_timestamp = np.squeeze(np.array(port_timestamp_bkgndint[direct_samples_ind]))
    #DPS says to exclude the last scan at TlastScanDirect.
    direct_mode_timestamp = np.around(np.nanmean(direct_mode_timestamp[:-1, ]))

    mass_charge_ratio = 2
    window = round(calibration_table[-1, 8], 1)
    bckgnd_mz2 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 30
    window = round(calibration_table[-1, 4], 1)
    bckgnd_mz30 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 32
    window = round(calibration_table[-1, 20], 1)
    bckgnd_mz32 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 34
    window = round(calibration_table[-1, 16], 1)
    bckgnd_mz34 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 40
    window = round(calibration_table[-1, 12], 1)
    bckgnd_mz40 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    mass_charge_ratio = 44
    window = round(calibration_table[-1, 24], 1)
    bckgnd_mz44 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    #average MSINLET-TEMP for direct time period here, call it bckgnd-Tdir
    direct_samples_ind = np.squeeze(np.where((port_timestamp_bkgndint_mcu >= port_timestamp_bkgndint[idx]) & (port_timestamp_bkgndint_mcu <= Tlast)))
    bckgnd_Tdir = np.nanmean(inlet_temp_bkgndint_mcu[direct_samples_ind[:-1]])

    #find gas_mode_bkgndint == 1 (nafion mode)
    ind_nafion = np.where(gas_mode_bkgndint == 1)[0]
    TlastScanNafion = port_timestamp_bkgndint_mcu[ind_nafion[-1]]

    #ID timestamp closest to TlastScanNafion - 180
    idx = (np.abs(port_timestamp_bkgndint-(TlastScanNafion - 180))).argmin()

    #subset the data collected in nafion mode
    nafion_samples_ind = np.where((port_timestamp_bkgndint >= port_timestamp_bkgndint[idx]) & (port_timestamp_bkgndint <= TlastScanNafion))
    nafion_samples = np.squeeze(np.array(L0_dissgas_bkgndint[nafion_samples_ind, :]))

    #calculate nafion mode timestamp
    nafion_mode_timestamp = np.array(port_timestamp_bkgndint[nafion_samples_ind])
    nafion_mode_timestamp = np.around(np.nanmean(np.squeeze(nafion_mode_timestamp)))

    mass_charge_ratio = 15
    window = round(calibration_table[-1, 0], 1)
    bckgnd_mz15 = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    #average MSINLET-TEMP for nafion time period here, call it bckgnd-Tnaf
    nafion_samples_ind = np.where((port_timestamp_bkgndint_mcu >= port_timestamp_bkgndint[idx]) & (port_timestamp_bkgndint_mcu <= TlastScanNafion))
    bckgnd_Tnaf = np.nanmean(inlet_temp_bkgndint_mcu[nafion_samples_ind])

    #average ph_meter_value for time period TlastScanNafion-1min:TlastScanNafion, call it msinlet_bkgphint
    nafion_samples_ind = np.where((port_timestamp_bkgndint_mcu >= TlastScanNafion - 60) & (port_timestamp_bkgndint_mcu <= TlastScanNafion))
    msinlet_bkgphint = np.absolute(np.nanmean(ph_meter_bkgndint_mcu[nafion_samples_ind]))

    return (bckgnd_mz2, bckgnd_mz15, bckgnd_mz30, bckgnd_mz32, bckgnd_mz34,
            bckgnd_mz40, bckgnd_mz44, bckgnd_Tnaf, bckgnd_Tdir, msinlet_bkgphint,
            nafion_mode_timestamp, direct_mode_timestamp)

def Cal1PreProcess(port_timestamp_calint01, L0_dissgas_calint01, gas_mode_calint01,
                   port_timestamp_calint01_mcu, ph_meter_calint01_mcu,
                   inlet_temp_calint01_mcu, mass_table, calibration_table):
    """
    Groups Calibration Solution 1 scans into Direct and Nafion mode
    subsets and extracts the averaged mz intensities, mode-
    averaged temperatures, and pH intensity needed for the L1
    DISSGAS, CALRANG, TSTAMP, MSINLET, and (for the sample water
    only) NAFEFF auxiliary data products.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 port timestamps for the mass spectral data.
    L0_dissgas... : ndarray
        L0 mass spectral intensity data set [A].
    gas_mode... : ndarray
        GASMODE auxiliary data product for this water source.
    port_timestamp..._mcu : ndarray
        MCU timestamps for the pH and temperature sensors.
    ph_meter..._mcu : ndarray
        Raw pH signal intensity [dimensionless].
    inlet_temp..._mcu : ndarray
        Sample or background temperature [deg_C]; -127 and 85 are
        replaced with NaN as instrument not-sampling fill values.
    mass_table : ndarray
        RGA mass lookup table from rga_status_process.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.

    Returns
    -------
    tuple of float
    cal1_mz15, cal1_mz44, cal1_Tnaf, cal1_Tdir, msinlet_cal1phint,
    nafion_mode_timestamp, and direct_mode_timestamp, in that order.
    """

    #replace bad data with nans
    inlet_temp_calint01_mcu[inlet_temp_calint01_mcu == -127] = np.nan
    inlet_temp_calint01_mcu[inlet_temp_calint01_mcu == 85] = np.nan

    #find gas_mode_calint01 == 0 (direct mode)
    ind_direct = np.where(gas_mode_calint01 == 0)[0]
    Tchange = port_timestamp_calint01_mcu[ind_direct[0]]
    Tlast = port_timestamp_calint01_mcu[ind_direct[-1]]

    #ID timestamp closest to Tlast - 60
    idx = (np.abs(port_timestamp_calint01-(Tlast - 60))).argmin()

    #subset the data collected in direct mode
    direct_samples_ind = np.where((port_timestamp_calint01 >= port_timestamp_calint01[idx]) & (port_timestamp_calint01 <= Tlast))
    direct_samples = np.squeeze(np.array(L0_dissgas_calint01[direct_samples_ind, :]))
    #DPS says to exclude the last scan at TlastScanDirect. This subsetted array then gets fed into the ave. routine.
    direct_samples = direct_samples[:-1, ]

    #calculate direct mode timestamp
    direct_mode_timestamp = np.squeeze(np.array(port_timestamp_calint01[direct_samples_ind]))
    #DPS says to exclude the last scan at TlastScanDirect.
    direct_mode_timestamp = np.around(np.nanmean(direct_mode_timestamp[:-1, ]))

    mass_charge_ratio = 44
    window = round(calibration_table[-1, 24], 1)
    cal1_mz44 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    #average MSINLET-TEMP for direct time period here, call it cal1-Tdir
    direct_samples_ind = np.squeeze(np.where((port_timestamp_calint01_mcu >= port_timestamp_calint01[idx]) & (port_timestamp_calint01_mcu <= Tlast)))
    cal1_Tdir = np.nanmean(inlet_temp_calint01_mcu[direct_samples_ind[:-1]])

    #find gas_mode_calint01 == 1 (nafion mode)
    ind_nafion = np.where(gas_mode_calint01 == 1)[0]
    TlastScanNafion = port_timestamp_calint01_mcu[ind_nafion[-1]]

    #ID timestamp closest to TlastScanNafion - 60
    idx = (np.abs(port_timestamp_calint01-(TlastScanNafion - 60))).argmin()

    #subset the data collected in nafion mode
    nafion_samples_ind = np.where((port_timestamp_calint01 >= port_timestamp_calint01[idx]) & (port_timestamp_calint01 <= TlastScanNafion))
    nafion_samples = np.squeeze(np.array(L0_dissgas_calint01[nafion_samples_ind, :]))

    #calculate nafion mode timestamp
    nafion_mode_timestamp = np.array(port_timestamp_calint01[nafion_samples_ind])
    nafion_mode_timestamp = np.around(np.nanmean(np.squeeze(nafion_mode_timestamp)))

    mass_charge_ratio = 15
    window = round(calibration_table[-1, 0], 1)
    cal1_mz15 = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    #average MSINLET-TEMP for nafion time period here, call it cal1-Tnaf
    nafion_samples_ind = np.where((port_timestamp_calint01_mcu >= port_timestamp_calint01[idx]) & (port_timestamp_calint01_mcu <= TlastScanNafion))
    cal1_Tnaf = np.nanmean(inlet_temp_calint01_mcu[nafion_samples_ind])

    #average ph_meter_value for time period TlastScanNafion-1min:TlastScanNafion, call it msinlet_cal1phint
    nafion_samples_ind = np.where((port_timestamp_calint01_mcu >= TlastScanNafion - 60) & (port_timestamp_calint01_mcu <= TlastScanNafion))
    msinlet_cal1phint = np.absolute(np.nanmean(ph_meter_calint01_mcu[nafion_samples_ind]))

    return (cal1_mz15, cal1_mz44, cal1_Tnaf, cal1_Tdir, msinlet_cal1phint,
            nafion_mode_timestamp, direct_mode_timestamp)

def Cal2PreProcess(port_timestamp_calint02, L0_dissgas_calint02, gas_mode_calint02,
                   port_timestamp_calint02_mcu, ph_meter_calint02_mcu,
                   inlet_temp_calint02_mcu, mass_table, calibration_table):
    """
    Groups Calibration Solution 2 scans into Direct and Nafion mode
    subsets and extracts the averaged mz intensities, mode-
    averaged temperatures, and pH intensity needed for the L1
    DISSGAS, CALRANG, TSTAMP, MSINLET, and (for the sample water
    only) NAFEFF auxiliary data products.

    Parameters
    ----------
    port_timestamp... : ndarray
        L0 port timestamps for the mass spectral data.
    L0_dissgas... : ndarray
        L0 mass spectral intensity data set [A].
    gas_mode... : ndarray
        GASMODE auxiliary data product for this water source.
    port_timestamp..._mcu : ndarray
        MCU timestamps for the pH and temperature sensors.
    ph_meter..._mcu : ndarray
        Raw pH signal intensity [dimensionless].
    inlet_temp..._mcu : ndarray
        Sample or background temperature [deg_C]; -127 and 85 are
        replaced with NaN as instrument not-sampling fill values.
    mass_table : ndarray
        RGA mass lookup table from rga_status_process.
    calibration_table : ndarray
        MASSP per-gas calibration coefficient table.

    Returns
    -------
    tuple of float
    cal2_mz15, cal2_mz44, cal2_Tnaf, cal2_Tdir, msinlet_cal2phint,
    nafion_mode_timestamp, and direct_mode_timestamp, in that order.
    """

    #replace bad data with nans
    inlet_temp_calint02_mcu[inlet_temp_calint02_mcu == -127] = np.nan
    inlet_temp_calint02_mcu[inlet_temp_calint02_mcu == 85] = np.nan

    #find gas_mode_calint02 == 0 (direct mode)
    ind_direct = np.where(gas_mode_calint02 == 0)[0]
    Tchange = port_timestamp_calint02_mcu[ind_direct[0]]
    Tlast = port_timestamp_calint02_mcu[ind_direct[-1]]

    #find gas_mode_calint02 == 1 (nafion mode)
    ind_nafion = np.where(gas_mode_calint02 == 1)[0]
    TlastScanNafion = port_timestamp_calint02_mcu[ind_nafion[-1]]

    #ID timestamp closest to TlastScanNafion - 60
    idx = (np.abs(port_timestamp_calint02-(TlastScanNafion - 60))).argmin()

    #subset the data collected in nafion mode
    nafion_samples_ind = np.where((port_timestamp_calint02 >= port_timestamp_calint02[idx]) & (port_timestamp_calint02 <= TlastScanNafion))
    nafion_samples = np.squeeze(np.array(L0_dissgas_calint02[nafion_samples_ind, :]))
    #DPS says to exclude the last scan at TlastScanNafion. This subsetted array then gets fed into the ave. routine.
    nafion_samples = nafion_samples[:-1, ]

    #calculate nafion mode timestamp
    nafion_mode_timestamp = np.squeeze(np.array(port_timestamp_calint02[nafion_samples_ind]))
    #DPS says to exclude the last scan at TlastScanNafion.
    nafion_mode_timestamp = np.around(np.nanmean(nafion_mode_timestamp[:-1, ]))

    mass_charge_ratio = 15
    window = round(calibration_table[-1, 0], 1)
    cal2_mz15 = average_mz(mass_charge_ratio, nafion_samples, mass_table, window)

    #average MSINLET-TEMP for nafion time period
    nafion_samples_ind = np.squeeze(np.where((port_timestamp_calint02_mcu >= port_timestamp_calint02[idx]) & (port_timestamp_calint02_mcu <= TlastScanNafion)))
    cal2_Tnaf = np.nanmean(inlet_temp_calint02_mcu[nafion_samples_ind[:-1]])

    #ID timestamp closest to Tlast - 60
    idx = (np.abs(port_timestamp_calint02-(Tlast - 60))).argmin()

    #subset the data collected in direct mode
    direct_samples_ind = np.where((port_timestamp_calint02 >= port_timestamp_calint02[idx]) & (port_timestamp_calint02 <= Tlast))
    direct_samples = np.squeeze(np.array(L0_dissgas_calint02[direct_samples_ind, :]))

    #calculate direct mode timestamp
    direct_mode_timestamp = np.array(port_timestamp_calint02[direct_samples_ind])
    direct_mode_timestamp = np.around(np.nanmean(np.squeeze(direct_mode_timestamp)))

    mass_charge_ratio = 44
    window = round(calibration_table[-1, 24], 1)
    cal2_mz44 = average_mz(mass_charge_ratio, direct_samples, mass_table, window)

    #average MSINLET-TEMP for direct time period here, call it cal2-Tdir
    direct_samples_ind = np.where((port_timestamp_calint02_mcu >= port_timestamp_calint02[idx]) & (port_timestamp_calint02_mcu <= Tlast))
    cal2_Tdir = np.nanmean(inlet_temp_calint02_mcu[direct_samples_ind])

    #average ph_meter_value for time period Tlast-1min:Tlast, call it msinlet_cal2phint
    direct_samples_ind = np.where((port_timestamp_calint02_mcu >= Tlast - 60) & (port_timestamp_calint02_mcu <= Tlast))
    msinlet_cal2phint = np.absolute(np.nanmean(ph_meter_calint02_mcu[direct_samples_ind]))
    #associate msinlet_smpphint with cal2_mz44 time stamp

    return (cal2_mz15, cal2_mz44, cal2_Tnaf, cal2_Tdir, msinlet_cal2phint,
            nafion_mode_timestamp, direct_mode_timestamp)

@deprecated
def GasModeDetermination(sample_valve1, sample_valve2, sample_valve3, sample_valve4):
    """
    Determines the auxiliary data product Gas Measurement Mode
    (GASMODE) from the four sample valve statuses.

    Parameters
    ----------
    sample_valve1 : ndarray of int
        Sample valve 1 status.
    sample_valve2 : ndarray of int
        Sample valve 2 status.
    sample_valve3 : ndarray of int
        Sample valve 3 status.
    sample_valve4 : ndarray of int
        Sample valve 4 status.

    Returns
    -------
    gasmode_array : ndarray of float
        GASMODE auxiliary data product: -1 for another operating mode,
        0 for Direct mode, 1 for Nafion mode; NaN where none of the
        valve combinations match.
    """

    data_array_size = np.shape(sample_valve1)
    gasmode_array = np.ones(data_array_size[0])
    gasmode_array[gasmode_array == 1] = np.nan

    ind = np.where(sample_valve4 == 1)
    gasmode_array[ind] = -1
    ind = np.where((sample_valve2 == 1) & (sample_valve1 == 0))
    gasmode_array[ind] = 0
    ind = np.where((sample_valve1 == 1) & (sample_valve2 == 0) & (sample_valve3 == 0))
    gasmode_array[ind] = 1

    return gasmode_array

@deprecated
def SmpModeDetermination(external_valve1_status, external_valve2_status,
                         external_valve3_status, external_valve4_status,
                         external_valve5_status):
    """
    Determines the auxiliary data product Sample Measurement Mode
    (SMPMODE) from the five external valve statuses.

    Parameters
    ----------
    external_valve1_status : ndarray of int
        External valve 1 status.
    external_valve2_status : ndarray of int
        External valve 2 status.
    external_valve3_status : ndarray of int
        External valve 3 status.
    external_valve4_status : ndarray of int
        External valve 4 status.
    external_valve5_status : ndarray of int
        External valve 5 status.

    Returns
    -------
    smpmode_array : ndarray of float
        SMPMODE auxiliary data product: 2, 1, -1, or -2 depending on
        the valve combination; NaN where none match.
    """

    data_array_size = np.shape(external_valve1_status)
    smpmode_array = np.ones(data_array_size[0])
    smpmode_array[smpmode_array == 1] = np.nan

    ind = np.where((external_valve1_status == 0) & (external_valve2_status == 0) &
                  (external_valve3_status == 0) & (external_valve4_status == 0) &
                  (external_valve5_status == 0))
    smpmode_array[ind] = 2

    ind = np.where((external_valve1_status == 1) & (external_valve2_status == 0) &
                  (external_valve3_status == 0) & (external_valve4_status == 0) &
                  (external_valve5_status == 1))
    smpmode_array[ind] = 1

    ind = np.where((external_valve1_status == 1) & (external_valve2_status == 1) &
                  (external_valve3_status == 0) & (external_valve4_status == 0) &
                  (external_valve5_status == 0))
    smpmode_array[ind] = -1

    ind = np.where((external_valve1_status == 1) & (external_valve2_status == 0) &
                  (external_valve3_status == 1) & (external_valve4_status == 0) &
                  (external_valve5_status == 0))
    smpmode_array[ind] = -2

    return smpmode_array
