#!/usr/bin/env python
"""
Module containing dissolved oxygen data processing functions for the
Ocean Observatories Initiative. Covers two instrument families:

    DOSTA -- Aanderaa Optode 4831, producing DOCONCS_L1 and DOXYGEN_L2.
    DOFST -- Sea-Bird Scientific SBE 43 and SBE 43F, producing DOCONCF_L2.

Processing configurations and data product levels are described in the
OOI Data Product Specifications for Oxygen Concentration from Stable
Instruments (1341-00520) and Fast Dissolved Oxygen (1341-00521).
"""
import numpy as np
import gsw

from ion_functions.data.generic_functions import replace_fill_with_nan

"""
DOSTA Processing Configurations:

    DOSTA configured for analog output of calphase and T_opt to CTD
    voltage channels:
        T_optode_degC = dosta_Topt_volt_to_degC(T_optode_volts)
        DOCONCS-DEG_L0 = dosta_phase_volt_to_degree(DOCONCS-VLT_L0)
        DOCONCS_L1 = do2_SVU(DOCONCS-DEG_L0, T_optode_degC, ...)
        DOXYGEN_L2 = do2_salinity_correction(DOCONCS_L1, ...)

    DOSTA configured for digital output of oxygen concentration to CTD
    RS-232:
        DOCONCS_L1 = o2_counts_to_uM(DOCONCS-CNT_L0)
        DOXYGEN_L2 = do2_salinity_correction(DOCONCS_L1, ...)

    DOSTA, autonomous operation, digital output of calphase and T_opt:
        DOCONCS_L1 = do2_SVU(DOCONCS-DEG_L0, T_optode_degC, ...)
        DOXYGEN_L2 = do2_salinity_correction(DOCONCS_L1, ...)


DOSTA DATA PRODUCTS:

    DOCONCS-CNT_L0 [counts]: oxygen concentration, uncorrected for
        salinity and pressure.
        (a) parsed from digital DOSTA output when routed through
            RS-232 CTD.
    DOCONCS-DEG_L0 [degrees]: CalPhase (calibrated phase). Two sources:
        (a) parsed from autonomous DOSTA digital output.
        (b) calculated by dosta_phase_volt_to_degree from
            DOCONCS-VLT_L0 when the DOSTA is connected to a CTD
            analog voltage channel.
    DOCONCS-VLT_L0 [volts]: CalPhase (calibrated phase).
        (a) parsed from analog DOSTA output when routed through CTD.

    DOCONCS_L1 [micro-mole/liter]: oxygen concentration, uncorrected
        for salinity and pressure.

        As of Aug 2015, SAF has not been updated to reflect the change
        of units from micro-mole/kg to micro-mole/liter.

        Two sources:
        (a) calculated by o2_counts_to_uM from DOCONCS-CNT_L0.
        (b) calculated by do2_SVU from DOCONCS-DEG_L0.

    DOXYGEN_L2 [micro-mole/kg]: oxygen concentration, corrected for
        salinity and pressure (and inherently temperature). Previously
        named DOCONCS_L2; renamed in DPS 1341-00520 version 1-02
        (2014-04-11).
        (a) calculated by do2_salinity_correction from DOCONCS_L1.


    Temperature in DOSTA data product calculations:

    The temperature input to the function do2_SVU should be the
    sensor's optode foil temperature [degC] measured by the optode's
    thermistor (Topt in the DPS). This variable is directly parsed from
    digital autonomous DOSTA data streams; if instead the DOSTA is
    connected to a CTD through an analog voltage channel, then
    Topt[degC] is calculated from Topt[V] by the function
    dosta_Topt_volt_to_degC.

    The temperature input to the function do2_salinity_correction is
    TEMPWAT from the co-located CTD.


DOFST DATA PRODUCTS:

    DOCONCF_L0 [counts]: represents either voltage_counts or frequency,
        depending on DOFST series.
    DOCONCF_L2 [micro-mole/kg]: oxygen concentration, corrected for
        temperature, salinity and pressure.
        (a) calculated by do2_dofst_volt and do2_dofst_frequency via
            the shared core function dofst_calc.
"""


def dosta_phase_volt_to_degree(phase_volt):
    """
    Convert DOCONCS-VLT_L0 analog voltage to DOCONCS-DEG_L0 calibrated
    phase in degrees for a DOSTA Aanderaa Optode connected to a SBE CTD
    0-5 V analog channel.

    Parameters
    ----------
    phase_volt : array_like
        Calibrated phase output from an Aanderaa optode converted to
        volts for analog transmission, DOCONCS-VLT_L0 [V].

    Returns
    -------
    phase_degree : ndarray
        Calibrated phase measured by the Aanderaa optode,
        DOCONCS-DEG_L0 [degrees].

    Notes
    -----
    The conversion coefficients (slope 12, offset 10) are universal for
    all Aanderaa optodes. Values were obtained from Shawn Sneddon at
    Xylem-Aanderaa.
    """
    # These coefficients to convert analog phase from volts to degrees
    # are universal for all Aanderaa optodes. Obtained from Shawn
    # Sneddon at Xylem-Aanderaa.
    phase_degree = 10.0 + 12.0 * phase_volt
    return phase_degree


def dosta_Topt_volt_to_degC(T_optode_volt):
    """
    Convert DOSTA foil temperature analog voltage to degrees C for an
    Aanderaa Optode connected to a SBE CTD 0-5 V analog channel.

    Parameters
    ----------
    T_optode_volt : array_like
        Optode foil temperature from an Aanderaa optode converted to
        volts for analog transmission [V].

    Returns
    -------
    t_optode_degc : ndarray
        Optode foil temperature measured by the Aanderaa optode
        thermistor [degC].

    Notes
    -----
    The conversion coefficients (slope 8, offset -5) are universal for
    all Aanderaa optodes. Values were obtained from Shawn Sneddon at
    Xylem-Aanderaa.

    The optode thermistor temperature is preferred for calculating
    oxygen concentration from DOSTA Aanderaa optodes because the
    permeability of the sensor foil to oxygen is sensitive to
    temperature and this measurement is situated directly at the foil.
    """
    # These coefficients to convert analog T_optode from volts to degC
    # are universal for all Aanderaa optodes. Obtained from Shawn
    # Sneddon at Xylem-Aanderaa.
    t_optode_degc = -5.0 + 8.0 * T_optode_volt
    return t_optode_degc


def o2_counts_to_uM(o2_counts):
    """
    Convert DOCONCS-CNT_L0 counts to DOCONCS_L1 dissolved oxygen
    concentration for a DOSTA Aanderaa Optode connected to a SBE 16+
    V2 CTD via RS-232.

    Parameters
    ----------
    o2_counts : array_like
        Raw oxygen concentration counts from the CTD,
        DOCONCS-CNT_L0 [counts].

    Returns
    -------
    DO : ndarray
        Dissolved oxygen concentration uncorrected for salinity and
        pressure, DOCONCS_L1 [micro-mole/L].

    Notes
    -----
    The DOCONCS_L1 data product has units of micromole/liter. The OOI
    Software Application Framework (SAF) incorrectly lists the units
    for this L1 product as micromole/kg.
    """
    # replace fill values with nan
    o2_counts = replace_fill_with_nan(None, o2_counts)

    DO = (o2_counts / 10000.0) - 10.0
    return DO


def do2_SVU(calphase, temp, csv, conc_coef=np.array([0.0, 1.0])):
    """
    Compute DOCONCS_L1 dissolved oxygen concentration from a DOSTA
    Aanderaa Optode using the Stern-Volmer-Uchida equation.

    Parameters
    ----------
    calphase : array_like
        Calibrated phase from the Aanderaa optode, DOCONCS-DEG_L0
        [degrees].
    temp : array_like
        Optode foil temperature from the optode thermistor [degC].
        See Notes regarding temperature source selection.
    csv : array_like, shape (..., 7)
        Stern-Volmer-Uchida calibration coefficients array. Elements
        are csv[0] through csv[6] corresponding to csv1 through csv7
        in OOI DPS 1341-00520 [various units].
    conc_coef : array_like, shape (..., 2), optional
        Secondary calibration coefficients [offset, slope] applied
        after the SVU equation. Default is [0.0, 1.0] (no correction).
        See Notes for details.

    Returns
    -------
    do : ndarray
        Dissolved oxygen concentration uncorrected for salinity and
        pressure, DOCONCS_L1 [micro-mole/L].

    Notes
    -----
    The Stern-Volmer-Uchida equation computes oxygen concentration as:

        Ksv = csv1 + csv2*T + csv3*T^2
        P0  = csv4 + csv5*T
        Pc  = csv6 + csv7*Pt
        O2  = [(P0/Pc) - 1] / Ksv

    where T is the optode foil temperature [degC] and Pt is the
    uncorrected (calibrated) phase [degrees].

    The secondary calibration (conc_coef) applies an offset and slope
    correction: DO = conc_coef[0] + conc_coef[1] * DO. Aanderaa uses
    two calibration procedures for the model 4831 optode. The primary
    multi-point calibration determines the SVU foil coefficients (csv).
    The secondary two-point calibration corrects the result against 0%
    and 100% oxygen data points to produce the conc_coef values. For
    new optodes or new SVU foil coefficient determinations, Aanderaa
    sets conc_coef to offset=0 and slope=1 by default. The conc_coef
    correction is applied automatically by the optode firmware but must
    be applied manually when oxygen is calculated externally, as in
    this function.

    The optode thermistor temperature should be used whenever possible.
    For OOI DOSTA model 4831 instruments the thermistor is situated
    directly at the sensor foil and the SVU calibration coefficients
    are derived in part to compensate for the change in oxygen
    permeability through the foil as a function of its temperature. On
    gliders, differences between CTD and optode temperature readings of
    1 degC can translate to approximately 5% differences in calculated
    oxygen concentration.

    The DOCONCS_L1 data product has units of micromole/liter. The OOI
    Software Application Framework (SAF) incorrectly lists the units
    for this L1 product as micromole/kg.
    """
    conc_coef = np.atleast_2d(conc_coef)
    # this will work for both old and new CI implementations of cal coeffs.
    csv = np.atleast_2d(csv)

    # Calculate DO using Stern-Volmer:
    ksv = csv[:, 0] + csv[:, 1]*temp + csv[:, 2]*(temp**2)
    p0 = csv[:, 3] + csv[:, 4]*temp
    pc = csv[:, 5] + csv[:, 6]*calphase
    do = ((p0/pc) - 1) / ksv

    # apply refurbishment calibration
    # conc_coef can be a 2D array of either 1 row or do.size rows.
    do = conc_coef[:, 0] + conc_coef[:, 1] * do
    return do


def do2_salinity_correction(DO, P, T, SP, lat, lon, sref=0, pref=0):
    """
    Correct DOCONCS_L1 for salinity and pressure to produce DOXYGEN_L2
    dissolved oxygen concentration.

    Parameters
    ----------
    DO : array_like
        Uncorrected dissolved oxygen concentration, DOCONCS_L1
        [micro-mole/L].
    P : array_like
        Sea pressure from the co-located CTD, PRESWAT_L1 [dbar].
        Interpolated to the DOSTA timestamp.
    T : array_like
        Water temperature from the co-located CTD, TEMPWAT_L1 [degC].
        Interpolated to the DOSTA timestamp.
    SP : array_like
        Practical salinity from the co-located CTD, PRACSAL_L2
        [unitless].
    lat : float or array_like
        Latitude of the instrument [decimal degrees North].
    lon : float or array_like
        Longitude of the instrument [decimal degrees East].
    sref : float, optional
        Reference salinity matching the Salinity setting configured
        in the Aanderaa optode firmware [unitless]. Default is 0.
    pref : float, optional
        Reference pressure for potential density calculation [dbar].
        Default is 0 dbar.

    Returns
    -------
    DOc : ndarray
        Dissolved oxygen concentration corrected for salinity and
        pressure, DOXYGEN_L2 [micro-mole/kg].

    Notes
    -----
    The correction proceeds in three steps.

    Step 1 -- Volume to mass unit conversion:
    Potential density (rho) referenced to pref is computed from
    absolute salinity and conservative temperature using the TEOS-10
    GSW library. The concentration is converted from per-volume to
    per-mass units:

        DO [umol/kg] = 1000 * DO [umol/L] / rho [kg/m^3]

    Step 2 -- Pressure correction (Uchida et al. 2008):

        DO = (1 + 0.032 * P / 1000) * DO

    Step 3 -- Salinity correction (Garcia and Gordon 1992, Table 1,
    combined fit):

        ts = ln[(298.15 - T) / (273.15 + T)]
        Bts = B0 + B1*ts + B2*ts^2 + B3*ts^3
        DO = exp[(SP - sref) * Bts + C0 * (SP^2 - sref^2)] * DO

    where B0=-6.24097e-3, B1=-6.93498e-3, B2=-6.90358e-3,
    B3=-4.29155e-3, C0=-3.11680e-7.

    The sref parameter corresponds to the preset Salinity setting in
    the Aanderaa optode configuration, typically 0 or 35.
    """
    # density calculation from GSW toolbox
    SA = gsw.SA_from_SP(SP, P, lon, lat)
    CT = gsw.CT_from_t(SA, T, P)
    pdens = gsw.rho(SA, CT, pref)  # potential referenced to p=0

    # Convert from volume to mass units:
    DO = 1000 * DO / pdens

    # Pressure correction:
    DO = (1 + (0.032 * P) / 1000) * DO

    # Salinity correction (Garcia and Gordon, 1992, combined fit):
    ts = np.log((298.15 - T) / (273.15 + T))
    B0 = -6.24097e-3
    B1 = -6.93498e-3
    B2 = -6.90358e-3
    B3 = -4.29155e-3
    C0 = -3.11680e-7
    Bts = B0 + B1*ts + B2*ts**2 + B3*ts**3
    DO = np.exp((SP - sref) * Bts + C0 * (SP**2 - sref**2)) * DO
    return DO


def do2_dofst_volt(voltage_counts, Voffset, Soc, A, B, C, E, P, T, SP, lat, lon):
    """
    OOI wrapper for DOCONCF_L2. Converts DOFST-A (SBE 43) voltage
    counts to dissolved oxygen concentration [micro-mole/kg].

    Parameters
    ----------
    voltage_counts : array_like
        Raw oxygen sensor voltage counts from the SBE 43,
        DOCONCF_L0 [counts].
    Voffset : float
        Calibration coefficient 1. Voltage offset [V].
    Soc : float
        Calibration coefficient 2. Oxygen signal slope [V^-1].
    A : float
        Calibration coefficient 3. Residual temperature correction
        factor A [unitless].
    B : float
        Calibration coefficient 4. Residual temperature correction
        factor B [unitless].
    C : float
        Calibration coefficient 5. Residual temperature correction
        factor C [unitless].
    E : float
        Calibration coefficient 6. Pressure correction factor
        [unitless].
    P : array_like
        Sea pressure from the co-located CTD, PRESWAT_L1 [dbar].
    T : array_like
        Water temperature from the co-located CTD, TEMPWAT_L1 [degC].
    SP : array_like
        Practical salinity from the co-located CTD, PRACSAL_L2
        [unitless].
    lat : float or array_like
        Latitude of the instrument [decimal degrees North].
    lon : float or array_like
        Longitude of the instrument [decimal degrees East].

    Returns
    -------
    do : ndarray
        Dissolved oxygen concentration corrected for temperature,
        salinity, and pressure, DOCONCF_L2 [micro-mole/kg].

    See Also
    --------
    dofst_calc : Core algorithm; use directly for multi-output access.
    """
    # replace fill values with nan
    voltage_counts = replace_fill_with_nan(None, voltage_counts)

    # convert voltage counts to volts
    volts = voltage_counts / 13107.

    do, do_int = dofst_calc(volts, Voffset, Soc, A, B, C, E, P, T, SP, lat, lon)
    return do


def do2_dofst_frequency(frequency, Foffset, Soc, A, B, C, E, P, T, SP, lat, lon):
    """
    OOI wrapper for DOCONCF_L2. Converts DOFST-K (SBE 43F) frequency
    to dissolved oxygen concentration [micro-mole/kg].

    Parameters
    ----------
    frequency : array_like
        Raw oxygen sensor frequency from the SBE 43F,
        DOCONCF_L0 [Hz].
    Foffset : float
        Calibration coefficient 1. Frequency offset [Hz].
    Soc : float
        Calibration coefficient 2. Oxygen signal slope [s].
    A : float
        Calibration coefficient 3. Residual temperature correction
        factor A [unitless].
    B : float
        Calibration coefficient 4. Residual temperature correction
        factor B [unitless].
    C : float
        Calibration coefficient 5. Residual temperature correction
        factor C [unitless].
    E : float
        Calibration coefficient 6. Pressure correction factor
        [unitless].
    P : array_like
        Sea pressure from the co-located CTD, PRESWAT_L1 [dbar].
    T : array_like
        Water temperature from the co-located CTD, TEMPWAT_L1 [degC].
    SP : array_like
        Practical salinity from the co-located CTD, PRACSAL_L2
        [unitless].
    lat : float or array_like
        Latitude of the instrument [decimal degrees North].
    lon : float or array_like
        Longitude of the instrument [decimal degrees East].

    Returns
    -------
    do : ndarray
        Dissolved oxygen concentration corrected for temperature,
        salinity, and pressure, DOCONCF_L2 [micro-mole/kg].

    See Also
    --------
    dofst_calc : Core algorithm; use directly for multi-output access.
    """
    # replace fill values with nan
    frequency = replace_fill_with_nan(None, frequency)

    do, do_int = dofst_calc(frequency, Foffset, Soc, A, B, C, E, P, T, SP, lat, lon)
    return do


# DOFST main sub-function
def dofst_calc(do_raw, offset, Soc, A, B, C, E, P, T, SP, lat, lon, freq=True):
    """
    Compute DOCONCF_L2 dissolved oxygen concentration from a DOFST
    (SBE 43 or SBE 43F) sensor using the Sea-Bird algorithm based on
    Owens and Millard (1985).

    Parameters
    ----------
    do_raw : array_like
        Oxygen sensor voltage [V] for the SBE 43, or frequency [Hz]
        for the SBE 43F.
    offset : float
        Voltage offset [V] for the SBE 43 (Voffset), or frequency
        offset [Hz] for the SBE 43F (Foffset). Calibration coefficient
        1 [V or Hz].
    Soc : float
        Oxygen signal slope. Calibration coefficient 2. Units are
        V^-1 for voltage input or s for frequency input.
    A : float
        Residual temperature correction factor A. Calibration
        coefficient 3 [unitless].
    B : float
        Residual temperature correction factor B. Calibration
        coefficient 4 [unitless].
    C : float
        Residual temperature correction factor C. Calibration
        coefficient 5 [unitless].
    E : float
        Pressure correction factor. Calibration coefficient 6
        [unitless].
    P : array_like
        Sea pressure from the co-located CTD, PRESWAT_L1 [dbar].
    T : array_like
        Water temperature from the co-located CTD, TEMPWAT_L1 [degC].
    SP : array_like
        Practical salinity from the co-located CTD, PRACSAL_L2
        [unitless].
    lat : float or array_like
        Latitude of the instrument [decimal degrees North].
    lon : float or array_like
        Longitude of the instrument [decimal degrees East].
    freq : bool, optional
        If True (default), do_raw is frequency [Hz] (SBE 43F path).
        If False, do_raw is voltage [V] (SBE 43 path). This flag is
        retained for compatibility but the conversion from counts to
        volts is performed by the calling wrapper before this function
        is invoked.

    Returns
    -------
    DO : ndarray
        Dissolved oxygen concentration corrected for temperature,
        salinity, and pressure, DOCONCF_L2 [micro-mole/kg].
    DO_int : ndarray
        Intermediate dissolved oxygen concentration before density
        conversion [mL/L].

    Notes
    -----
    The Sea-Bird algorithm (Owens and Millard 1985) with the derivative
    term removed (tau set to zero per DPS recommendation) is:

        Oxsol(T, S) = exp[A0 + A1*Ts + A2*Ts^2 + A3*Ts^3 + A4*Ts^4
                          + A5*Ts^5 + S*(B0 + B1*Ts + B2*Ts^2
                          + B3*Ts^3) + C0*S^2]

    where Ts = ln[(298.15 - T) / (273.15 + T)] and Oxsol is the
    oxygen saturation value (Garcia and Gordon 1992, Table 1, 1st
    column fit).

    The intermediate oxygen concentration in mL/L is:

        DO_int = Soc * (do_raw + offset) * Oxsol(T, S)
                 * (1 + A*T + B*T^2 + C*T^3) * exp(E*P / K)

    where K = T + 273.15. The final DOCONCF_L2 in micro-mole/kg is:

        DO = DO_int * 44660 / pot_rho_t

    where pot_rho_t is the potential density of seawater computed from
    absolute salinity, in-situ temperature, and sea pressure using the
    TEOS-10 GSW library (gsw.pot_rho_t_exact).

    The derivative term [tau(T,P) * dV/dt] is removed by setting
    tau = 0, as recommended by Sea-Bird and the OOI DPS, to prevent
    amplification of residual noise especially in deep water.
    """
    # Get potential density using the TEOS-10 toolbox
    SA = gsw.SA_from_SP(SP, P, lon, lat)
    pot_rho_t = gsw.pot_rho_t_exact(SA, T, P, 0)

    # Oxygen saturation value using Garcia and Gordon (1992) fit to
    # Benson and Krause data (Table 1, 1st column)
    #   empirical polynomial coefficients (not calibration coeffs)
    A0 = 2.00907
    A1 = 3.22014
    A2 = 4.0501
    A3 = 4.94457
    A4 = -0.256847
    A5 = 3.88767
    B0 = -0.00624523
    B1 = -0.00737614
    B2 = -0.010341
    B3 = -0.00817083
    C0 = -0.000000488682
    temp_K = T + 273.15  # temperature in Kelvin
    Ts = np.log((298.15 - T) / (temp_K))
    Oxsol = np.exp(
        A0 + A1*Ts + A2*Ts**2 + A3*Ts**3 + A4*Ts**4 + A5*Ts**5 +
        SP * (B0 + B1*Ts + B2*Ts**2 + B3*Ts**3) +
        C0*SP**2)

    if not freq:
        # convert voltage counts to volts
        do_raw = do_raw / 13107.

    # Intermediate step: Dissolved Oxygen concentration in [mL/L]
    DO_int = (Soc * (do_raw + offset) * Oxsol
              * (1.0 + A*T + B*T**2 + C*T**3)
              * np.exp((E * P)/temp_K))

    # Correct DO_int for Potential Density and convert to [micromole/Kg]
    DO = DO_int * 44660. / (pot_rho_t)
    return (DO, DO_int)
