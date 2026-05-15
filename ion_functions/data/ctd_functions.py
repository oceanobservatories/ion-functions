#!/usr/bin/env python
"""
CTD data processing functions for the Ocean Observatories Initiative.

Converts raw Sea-Bird Electronics CTD data into calibrated L1 engineering
products (TEMPWAT, PRESWAT, CONDWAT) and computes the L2 derived products
Practical Salinity (PRACSAL) and in-situ Density (DENSITY) using the TEOS-10
GSW library.
"""

# Import Numpy and the GSW library
import numpy as np
import gsw


def ctd_sbe16plus_tempwat(t0, a0, a1, a2, a3):
    """
    Compute TEMPWAT_L1 from SBE 16Plus raw temperature counts.

    Applies the Sea-Bird ITS-90 thermistor conversion for the SBE 16Plus
    V2 (OutputFormat 0). Covers all CTDBP instrument series and CTDPF
    series A and B.

    Parameters
    ----------
    t0 : array_like
        Raw temperature counts (TEMPWAT_L0) [counts].
    a0 : float
        Temperature calibration coefficient a0.
    a1 : float
        Temperature calibration coefficient a1.
    a2 : float
        Temperature calibration coefficient a2.
    a3 : float
        Temperature calibration coefficient a3.

    Returns
    -------
    t : ndarray
        Seawater temperature (TEMPWAT_L1) [degC, ITS-90].

    Notes
    -----
    Algorithm converts the 6-character hex integer t0 (already decoded
    to decimal counts by the CTD driver) to ITS-90 temperature via:

        MV = (t0 - 524288) / 1.6e7
        R  = (MV * 2.9e9 + 1.024e8) / (2.048e4 - MV * 2.0e5)
        T  = 1 / (a0 + a1*ln(R) + a2*ln(R)^2 + a3*ln(R)^3) - 273.15

    Calibration coefficients a0-a3 are from factory calibration sheets.
    """
    mv = (t0 - 524288) / 1.6e7
    r = (mv * 2.9e9 + 1.024e8) / (2.048e4 - mv * 2.0e5)
    t = 1 / (a0 + a1 * np.log(r) + a2 * np.log(r)**2 + a3 * np.log(r)**3) - 273.15
    return t


def ctd_sbe37im_tempwat_instrument_recovered(t0, a0, a1, a2, a3):
    """
    Compute TEMPWAT_L1 from SBE 37IM instrument-recovered counts.

    Applies the Sea-Bird ITS-90 thermistor conversion for data recovered
    directly from an SBE 37IM instrument (all series), where the raw
    value is stored as a decimal count rather than scaled engineering
    units.

    Parameters
    ----------
    t0 : array_like
        Raw temperature counts (TEMPWAT_L0) recovered directly from
        the CTD instrument [counts].
    a0 : float
        Temperature calibration coefficient a0.
    a1 : float
        Temperature calibration coefficient a1.
    a2 : float
        Temperature calibration coefficient a2.
    a3 : float
        Temperature calibration coefficient a3.

    Returns
    -------
    t : ndarray
        Seawater temperature (TEMPWAT_L1) [degC, ITS-90].

    Notes
    -----
    Unlike telemetered and recovered_host data (see
    ctd_sbe37im_tempwat), instrument-recovered data retains the raw
    count t0, which is used directly in the temperature equation:

        T = 1 / (a0 + a1*ln(t0) + a2*ln(t0)^2 + a3*ln(t0)^3) - 273.15

    As of June 2016 this processing path was not included in the
    TEMPWAT DPS (DCN 1341-00010). Calibration coefficients a0--a3 are
    from factory calibration sheets.
    """
    t = 1 / (a0 + a1 * np.log(t0) + a2 * np.log(t0)**2 + a3 * np.log(t0)**3) - 273.15
    return t


def ctd_sbe37im_tempwat(t0):
    """
    Compute TEMPWAT_L1 from SBE 37IM telemetered/recovered_host counts.

    Applies the linear scaling for SBE 37IM (all series) telemetered
    and recovered_host data streams (CTDMO instrument class). The
    instrument pre-scales raw counts to engineering units before
    transmission; no calibration coefficients are required.

    Parameters
    ----------
    t0 : array_like
        Raw temperature counts (TEMPWAT_L0) [counts].

    Returns
    -------
    t : ndarray
        Seawater temperature (TEMPWAT_L1) [degC, ITS-90].

    Notes
    -----
    The SBE 37IM encodes temperature as engineering units in
    hexadecimal (OutputFormat 0). The conversion is:

        T = t0 / 10000 - 10

    This function does not apply to instrument-recovered data; use
    ctd_sbe37im_tempwat_instrument_recovered for that stream.
    """
    t = t0 / 10000.0 - 10.0
    return t


def ctd_sbe52mp_tempwat(t0):
    """
    Compute TEMPWAT_L1 from SBE 52MP raw temperature counts.

    Applies the linear scaling for SBE 52MP instruments (CTDPF series
    C, K, and L).

    Parameters
    ----------
    t0 : array_like
        Raw temperature counts (TEMPWAT_L0) [counts].

    Returns
    -------
    t : ndarray
        Seawater temperature (TEMPWAT_L1) [degC, ITS-90].

    Notes
    -----
    The conversion is:

        T = t0 / 10000 - 5
    """
    t = t0 / 10000.0 - 5.0
    return t


def ctd_sbe16plus_preswat(p0, t0, ptempa0, ptempa1, ptempa2,
                          ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
                          pa0, pa1, pa2, offset=0):
    """
    Compute PRESWAT_L1 from SBE 16Plus strain-gauge pressure counts.

    Applies the Sea-Bird strain-gauge pressure conversion for SBE 16Plus
    instruments equipped with an internal strain-gauge pressure sensor
    (PType=1). Covers most CTDBP instrument series (exceptions: N and O)
    and CTDPF series A and B.

    Parameters
    ----------
    p0 : array_like
        Raw pressure counts (PRESWAT_L0) [counts].
    t0 : array_like
        Raw temperature counts from the pressure-sensor thermistor
        [counts].
    ptempa0 : float
        Strain-gauge pressure calibration coefficient PTEMPA0.
    ptempa1 : float
        Strain-gauge pressure calibration coefficient PTEMPA1.
    ptempa2 : float
        Strain-gauge pressure calibration coefficient PTEMPA2.
    ptca0 : float
        Strain-gauge pressure calibration coefficient PTCA0.
    ptca1 : float
        Strain-gauge pressure calibration coefficient PTCA1.
    ptca2 : float
        Strain-gauge pressure calibration coefficient PTCA2.
    ptcb0 : float
        Strain-gauge pressure calibration coefficient PTCB0.
    ptcb1 : float
        Strain-gauge pressure calibration coefficient PTCB1.
    ptcb2 : float
        Strain-gauge pressure calibration coefficient PTCB2.
    pa0 : float
        Strain-gauge pressure calibration coefficient PA0.
    pa1 : float
        Strain-gauge pressure calibration coefficient PA1.
    pa2 : float
        Strain-gauge pressure calibration coefficient PA2.
    offset : float, optional
        Druck sensor offset correction [dbar]. Default is 0.

    Returns
    -------
    p : ndarray
        Seawater pressure (PRESWAT_L1) [dbar].

    Notes
    -----
    Algorithm (from OOI DPS DCN 1341-00020, PType=1):

        t_v = t0 / 13107
        t   = PTEMPA0 + PTEMPA1*t_v + PTEMPA2*t_v^2
        x   = p0 - PTCA0 - PTCA1*t - PTCA2*t^2
        n   = x * PTCB0 / (PTCB0 + PTCB1*t + PTCB2*t^2)
        p_psi = PA0 + PA1*n + PA2*n^2
        P_L1  = p_psi * 0.689475729 - 10.1325 + offset

    All calibration coefficients are from factory calibration sheets.
    The optional offset corrects a known Druck sensor bias.
    """
    # compute calibration parameters
    tv = t0 / 13107.0
    t = ptempa0 + ptempa1 * tv + ptempa2 * tv**2
    x = p0 - ptca0 - ptca1 * t - ptca2 * t**2
    n = x * ptcb0 / (ptcb0 + ptcb1 * t + ptcb2 * t**2)

    # compute pressure in psi, rescale and compute in dbar and return
    p_psi = pa0 + pa1 * n + pa2 * n**2
    p_dbar = (p_psi * 0.689475729) - 10.1325
    return p_dbar + offset


def ctd_sbe16digi_preswat(p0, t0, C1, C2, C3, D1, D2, T1, T2, T3, T4, T5):
    """
    Compute PRESWAT_L1 from SBE 16Plus digiquartz pressure counts.

    Applies the Sea-Bird digiquartz pressure conversion for SBE 16Plus
    instruments equipped with a digiquartz pressure sensor (PType=3).
    Applies exclusively to CTDBP-N and CTDBP-O instruments.

    Parameters
    ----------
    p0 : array_like
        Raw pressure counts (PRESWAT_L0) [counts].
    t0 : array_like
        Raw temperature counts from the pressure-sensor thermistor
        [counts].
    C1 : float
        Digiquartz pressure calibration coefficient C1.
    C2 : float
        Digiquartz pressure calibration coefficient C2.
    C3 : float
        Digiquartz pressure calibration coefficient C3.
    D1 : float
        Digiquartz pressure calibration coefficient D1.
    D2 : float
        Digiquartz pressure calibration coefficient D2.
    T1 : float
        Digiquartz pressure calibration coefficient T1.
    T2 : float
        Digiquartz pressure calibration coefficient T2.
    T3 : float
        Digiquartz pressure calibration coefficient T3.
    T4 : float
        Digiquartz pressure calibration coefficient T4.
    T5 : float
        Digiquartz pressure calibration coefficient T5.

    Returns
    -------
    p : ndarray
        Seawater pressure (PRESWAT_L1) [dbar].

    Notes
    -----
    Algorithm (from OOI DPS DCN 1341-00020, PType=3):

        pf  = p0 / 256          (pressure frequency in Hz)
        t_v = t0 / 13107        (thermistor voltage)
        U   = 23.7*(t_v + 9.7917) - 273.15
        C   = C1 + C2*U + C3*U^2
        D   = D1 + D2*U
        T0  = T1 + T2*U + T3*U^2 + T4*U^3 + T5*U^4
        T   = (1/pf) * 1e6      (pressure period in microseconds)
        p_psi = C*(1 - T0^2/T^2)*(1 - D*(1 - T0^2/T^2))
        P_L1  = p_psi * 0.689475729 - 10.1325

    All calibration coefficients are from factory calibration sheets.
    """
    # Convert raw pressure input to frequency [Hz]
    pf = p0 / 256.0

    # Convert raw temperature input to voltage
    tv = t0 / 13107.0

    # Calculate U (thermistor temp):
    U = (23.7 * (tv + 9.7917)) - 273.15

    # Calculate calibration parameters
    C = C1 + C2 * U + C3 * U**2
    D = D1 + D2 * U
    T0 = T1 + T2 * U + T3 * U**2 + T4 * U**3 + T5 * U**4

    # Calculate T (pressure period, in microseconds):
    T = (1.0 / pf) * 1.0e6

    # compute pressure in psi, rescale and compute in dbar and return
    p_psi = C * (1.0 - T0**2 / T**2) * (1.0 - D * (1.0 - T0**2 / T**2))
    p_dbar = (p_psi * 0.689475729) - 10.1325
    return p_dbar


def ctd_sbe37im_preswat_instrument_recovered(p0, pt0, ptempa0, ptempa1, ptempa2,
                                             ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
                                             pa0, pa1, pa2):
    """
    Compute PRESWAT_L1 from SBE 37IM instrument-recovered pressure counts.

    Applies the Sea-Bird strain-gauge pressure conversion for data
    recovered directly from an SBE 37IM instrument (all series).

    Parameters
    ----------
    p0 : array_like
        Raw pressure counts (PRESWAT_L0) recovered directly from the
        CTD instrument [counts].
    pt0 : array_like
        Raw temperature counts from the pressure-sensor thermistor
        [counts].
    ptempa0 : float
        Strain-gauge pressure calibration coefficient PTEMPA0.
    ptempa1 : float
        Strain-gauge pressure calibration coefficient PTEMPA1.
    ptempa2 : float
        Strain-gauge pressure calibration coefficient PTEMPA2.
    ptca0 : float
        Strain-gauge pressure calibration coefficient PTCA0.
    ptca1 : float
        Strain-gauge pressure calibration coefficient PTCA1.
    ptca2 : float
        Strain-gauge pressure calibration coefficient PTCA2.
    ptcb0 : float
        Strain-gauge pressure calibration coefficient PTCB0.
    ptcb1 : float
        Strain-gauge pressure calibration coefficient PTCB1.
    ptcb2 : float
        Strain-gauge pressure calibration coefficient PTCB2.
    pa0 : float
        Strain-gauge pressure calibration coefficient PA0.
    pa1 : float
        Strain-gauge pressure calibration coefficient PA1.
    pa2 : float
        Strain-gauge pressure calibration coefficient PA2.

    Returns
    -------
    p : ndarray
        Seawater pressure (PRESWAT_L1) [dbar].

    Notes
    -----
    Unlike telemetered and recovered_host data (see ctd_sbe37im_preswat),
    instrument-recovered data retains raw strain-gauge counts. The
    algorithm matches the SBE 16Plus strain-gauge path but uses pt0
    directly (not scaled to voltage):

        t   = PTEMPA0 + PTEMPA1*pt0 + PTEMPA2*pt0^2
        x   = p0 - PTCA0 - PTCA1*t - PTCA2*t^2
        n   = x * PTCB0 / (PTCB0 + PTCB1*t + PTCB2*t^2)
        p_psi = PA0 + PA1*n + PA2*n^2
        P_L1  = p_psi * 0.689475729 - 10.1325

    As of June 2016 this processing path was not included in the
    PRESWAT DPS (DCN 1341-00020). Calibration coefficients are from
    factory calibration sheets.
    """
    # compute calibration parameters
    t = ptempa0 + ptempa1 * pt0 + ptempa2 * pt0**2
    x = p0 - ptca0 - ptca1 * t - ptca2 * t**2
    n = x * ptcb0 / (ptcb0 + ptcb1 * t + ptcb2 * t**2)

    # compute pressure in psi, rescale and compute in dbar and return
    p_psi = pa0 + pa1 * n + pa2 * n**2
    p_dbar = (p_psi * 0.689475729) - 10.1325
    return p_dbar


def ctd_sbe37im_preswat(p0, p_range_psia):
    """
    Compute PRESWAT_L1 from SBE 37IM telemetered/recovered_host counts.

    Applies the linear pressure scaling for SBE 37IM (all series)
    telemetered and recovered_host data streams (CTDMO instrument class).

    Parameters
    ----------
    p0 : array_like
        Raw pressure counts (PRESWAT_L0) [counts].
    p_range_psia : float
        Pressure range calibration coefficient [psia], a factory-set
        value stored in the instrument metadata.

    Returns
    -------
    p : ndarray
        Seawater pressure (PRESWAT_L1) [dbar].

    Notes
    -----
    Algorithm (from OOI DPS DCN 1341-00020, SBE 37IM):

        P_range_dbar = (p_range_psia - 14.7) * 0.6894757
        P_L1 = p0 * P_range_dbar / (0.85 * 65536) - 0.05 * P_range_dbar

    The pressure range is a factory-set calibration coefficient.
    This function does not apply to instrument-recovered data; use
    ctd_sbe37im_preswat_instrument_recovered for that stream.
    """
    # compute pressure range in units of dbar
    p_range_dbar = (p_range_psia - 14.7) * 0.6894757

    # compute pressure in dbar and return
    p_dbar = p0 * p_range_dbar / (0.85 * 65536.0) - 0.05 * p_range_dbar
    return p_dbar


def ctd_glider_preswat(pr_bar):
    """
    Compute PRESWAT_L1 from glider CTD pressure in bar.

    Converts pressure reported by a Seabird CTD installed on a glider
    (CTDGV instrument class) from bar to dbar.

    Parameters
    ----------
    pr_bar : array_like
        Seawater pressure from glider (PRESWAT_L0) [bar].

    Returns
    -------
    pr_dbar : ndarray
        Seawater pressure (PRESWAT_L1) [dbar].

    Notes
    -----
    Conversion:

        P_L1 [dbar] = pr_bar * 10
    """
    pr_dbar = pr_bar * 10.0
    return pr_dbar


def ctd_sbe52mp_preswat(p0):
    """
    Compute PRESWAT_L1 from SBE 52MP raw pressure counts.

    Applies the linear pressure scaling for SBE 52MP instruments
    (CTDPF series C, K, and L).

    Parameters
    ----------
    p0 : array_like
        Raw pressure counts (PRESWAT_L0) [counts].

    Returns
    -------
    p : ndarray
        Seawater pressure (PRESWAT_L1) [dbar].

    Notes
    -----
    Conversion:

        P_L1 = p0 / 100 - 10
    """
    p_dbar = p0 / 100.0 - 10.0
    return p_dbar


def ctd_sbe16plus_condwat(c0, t1, p1, g, h, i, j, cpcor, ctcor):
    """
    Compute CONDWAT_L1 from SBE 16Plus raw conductivity counts.

    Applies the Sea-Bird conductivity frequency conversion for SBE
    16Plus instruments. Covers all CTDBP instrument series and CTDPF
    series A and B.

    Parameters
    ----------
    c0 : array_like
        Raw conductivity counts (CONDWAT_L0) [counts].
    t1 : array_like
        Seawater temperature (TEMPWAT_L1) [degC].
    p1 : array_like
        Seawater pressure (PRESWAT_L1) [dbar].
    g : float
        Conductivity calibration coefficient g.
    h : float
        Conductivity calibration coefficient h.
    i : float
        Conductivity calibration coefficient i.
    j : float
        Conductivity calibration coefficient j.
    cpcor : float
        Conductivity calibration coefficient CPcor.
    ctcor : float
        Conductivity calibration coefficient CTcor.

    Returns
    -------
    c : ndarray
        Seawater conductivity (CONDWAT_L1) [S m^-1].

    Notes
    -----
    Algorithm (from OOI DPS DCN 1341-00030, SBE 16Plus):

        f    = (c0 / 256) / 1000    (conductivity frequency in kHz)
        C_L1 = (g + h*f^2 + i*f^3 + j*f^4) / (1 + CTcor*T + CPcor*P)

    where T is TEMPWAT_L1 (degC) and P is PRESWAT_L1 (dbar).
    All calibration coefficients are from factory calibration sheets.
    """
    # convert raw conductivity measurement to frequency
    f = (c0 / 256.0) / 1000.0

    # calculate conductivity [S m-1]
    c = (g + h * f**2 + i * f**3 + j * f**4) / (1 + ctcor * t1 + cpcor * p1)
    return c


def ctd_sbe37im_condwat_instrument_recovered(c0, t1, p1, g, h, i, j, cpcor, ctcor, wbotc):
    """
    Compute CONDWAT_L1 from SBE 37IM instrument-recovered counts.

    Applies the Sea-Bird conductivity frequency conversion for data
    recovered directly from an SBE 37IM instrument (all series).

    Parameters
    ----------
    c0 : array_like
        Raw conductivity counts (CONDWAT_L0) recovered directly from
        the CTD instrument [counts].
    t1 : array_like
        Seawater temperature (TEMPWAT_L1) [degC].
    p1 : array_like
        Seawater pressure (PRESWAT_L1) [dbar].
    g : float
        Conductivity calibration coefficient g.
    h : float
        Conductivity calibration coefficient h.
    i : float
        Conductivity calibration coefficient i.
    j : float
        Conductivity calibration coefficient j.
    cpcor : float
        Conductivity calibration coefficient CPcor.
    ctcor : float
        Conductivity calibration coefficient CTcor.
    wbotc : float
        Conductivity calibration coefficient wbotc.

    Returns
    -------
    c : ndarray
        Seawater conductivity (CONDWAT_L1) [S m^-1].

    Notes
    -----
    Algorithm applies a wbotc correction to the frequency before the
    standard polynomial evaluation:

        f    = (c0/256)/1000 * sqrt(1 + wbotc*T)
        C_L1 = (g + h*f^2 + i*f^3 + j*f^4) / (1 + CTcor*T + CPcor*P)

    where T is TEMPWAT_L1 (degC) and P is PRESWAT_L1 (dbar).
    As of June 2016 this processing path was not included in the
    CONDWAT DPS (DCN 1341-00030). Calibration coefficients are from
    factory calibration sheets.
    """
    # convert raw conductivity measurement to frequency
    f = (c0 / 256.0) / 1000.0 * np.sqrt(1.0 + wbotc * t1)

    # calculate conductivity [S m-1]
    c = (g + h * f**2 + i * f**3 + j * f**4) / (1 + ctcor * t1 + cpcor * p1)
    return c


def ctd_sbe37im_condwat(c0):
    """
    Compute CONDWAT_L1 from SBE 37IM telemetered/recovered_host counts.

    Applies the linear conductivity scaling for SBE 37IM (all series)
    telemetered and recovered_host data streams (CTDMO instrument class).

    Parameters
    ----------
    c0 : array_like
        Raw conductivity counts (CONDWAT_L0) [counts].

    Returns
    -------
    c : ndarray
        Seawater conductivity (CONDWAT_L1) [S m^-1].

    Notes
    -----
    The SBE 37IM encodes conductivity as engineering units in
    hexadecimal (OutputFormat 0). The conversion is:

        C_L1 = c0 / 100000 - 0.5

    This function does not apply to instrument-recovered data; use
    ctd_sbe37im_condwat_instrument_recovered for that stream.
    """
    c = c0 / 100000.0 - 0.5
    return c


def ctd_sbe52mp_condwat(c0):
    """
    Compute CONDWAT_L1 from SBE 52MP raw conductivity counts.

    Applies the linear conductivity scaling for SBE 52MP instruments
    (CTDPF series C, K, and L).

    Parameters
    ----------
    c0 : array_like
        Raw conductivity counts (CONDWAT_L0) [counts].

    Returns
    -------
    c : ndarray
        Seawater conductivity (CONDWAT_L1) [S m^-1].

    Notes
    -----
    Two-step conversion:

        C [mmho/cm] = c0 / 10000 - 0.5
        C [S m^-1]  = C [mmho/cm] * 0.1
    """
    c_mmho_cm = c0 / 10000.0 - 0.5
    c_S_m = 0.1 * c_mmho_cm
    return c_S_m


def ctd_pracsal(c, t, p):
    """
    Compute PRACSAL_L2 from L1 conductivity, temperature, and pressure.

    Calculates Practical Salinity on the PSS-78 scale using the TEOS-10
    GSW library function gsw.SP_from_C, which implements the UNESCO 1983
    PSS-78 algorithm with the Hill et al. (1986) extension for SP < 2.

    Parameters
    ----------
    c : array_like
        Seawater conductivity (CONDWAT_L1) [S m^-1].
    t : array_like
        Seawater temperature (TEMPWAT_L1) [degC, ITS-90].
    p : array_like
        Seawater pressure (PRESWAT_L1) [dbar].

    Returns
    -------
    SP : ndarray
        Practical salinity (PRACSAL_L2) [PSS-78, unitless].

    Notes
    -----
    Conductivity is converted from S m^-1 to mS cm^-1 (multiply by 10)
    before passing to gsw.SP_from_C(C, t, p), which expects mS cm^-1.
    """
    # Convert L1 Conductivity from S/m to mS/cm
    C10 = c * 10.0

    # Calculate the Practical Salinity (PSS-78) [unitless]
    SP = gsw.SP_from_C(C10, t, p)
    return SP


def ctd_density(SP, t, p, lat, lon):
    """
    Compute DENSITY_L2 from practical salinity, temperature, pressure,
    and position using TEOS-10.

    Calculates in-situ seawater density using the TEOS-10 GSW library
    via a three-step chain: Practical Salinity to Absolute Salinity,
    in-situ temperature to Conservative Temperature, then density from
    the computationally-efficient 48-term expression.

    Parameters
    ----------
    SP : array_like
        Practical salinity (PRACSAL_L2) [PSS-78, unitless].
    t : array_like
        Seawater temperature (TEMPWAT_L1) [degC, ITS-90].
    p : array_like
        Seawater pressure (PRESWAT_L1) [dbar].
    lat : array_like
        Latitude of observation [decimal degrees north].
    lon : array_like
        Longitude of observation [decimal degrees east].

    Returns
    -------
    rho : ndarray
        In-situ seawater density (DENSITY_L2) [kg m^-3].

    Notes
    -----
    TEOS-10 processing chain:

        SA  = gsw.SA_from_SP(SP, p, lon, lat)   (Absolute Salinity)
        CT  = gsw.CT_from_t(SA, t, p)           (Conservative Temp)
        rho = gsw.rho(SA, CT, p)                (density, 48-term)

    For moored instruments, latitude (lat) and longitude (lon)
    are from the mooring position metadata; for gliders, they are
    the vehicle's position at the time of each measurement.
    """
    # Calculate the density [kg m-3]
    sa = gsw.SA_from_SP(SP, p, lon, lat)  # absolute salinity
    ct = gsw.CT_from_t(sa, t, p)  # conservative temperature
    rho = gsw.rho(sa, ct, p)  # density
    return rho
