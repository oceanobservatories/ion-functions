#!/usr/bin/env python

# @package ion_functions.data.ctd_functions
# @file ion_functions/data/ctd_functions.py
# @author Christopher Wingard
# @brief Module containing CTD related data-calculations.

"""
# Outline
The CTD instruments measure temperature, conductivity, and pressure (depth), 
density, and practical salinity of the water column. The data products derived
from these measurements are:

  * \\(\\text{TEMPWAT_L1}\\): Water Temperature
  * \\(\\text{PRESWAT_L1}\\): Pressure (Depth)
  * \\(\\text{CONDWAT_L1}\\): Conductivity
  * \\(\\text{PRACSAL_L2}\\): Practical Salinity
  * \\(\\text{DENSITY_L2}\\): Density

General theoretical background for each is provided in the Theory section, while 
algorithm implementation specifics for each data product are provided in the 
corresponding function docstring.


# Theory

## Pressure
Absolute pressure is one of the variables directly measured by CTD instruments 
and is calculated from the raw hexadecimal data provided by the instrument. 
The \\(\\text{SBE 16plusV2}\\) model will output "raw frequencies and voltages in hexadecimal" 
as it is referred to in the Sea-Bird manuals. The \\(\\text{37IM}\\) model will output 
"engineering units in Hex". This requires that the L0 pressure data product, 
which is a hexadecimal string, be converted to decimal and scaled according to 
the CTD manual (different for each CTD make/model). Conversion and scaling 
(described herein) results in sea pressure in decibars (dbar). 
Note that the data product is named Pressure (Depth) because, though depth is 
not calculated, 1 dbar in pressure is approximately equal to one meter in depth.

It is worth noting that sea pressure is not exactly the same as the hydrostatic 
pressure of the water column. The CTD's internal pressure sensor measures 
absolute pressure: the pressure of the water column (hydrostatic pressure) 
plus the current atmospheric pressure at the sea surface. Sea pressure is 
calculated by the CTD's internal software by subtracting one standard 
atmosphere from the measured absolute pressure. Because the actual atmospheric 
pressure is not necessarily exactly one standard atmosphere, sea pressure is 
not necessarily identical to the hydrostatic pressure. 
Also, Sea-Bird uses a slightly different constant for one atmosphere: 
according to their manuals, the CTDs internal software converts absolute 
pressure to sea pressure by subtracting 
\\(14.7 \\text{psi} * 0.689476 \\text{dbar}/\\text{psi} = 10.1352972 \\text{dbar}\\). 
This is different from the TEOS-10 standard atmospheric pressure of 
\\(101325 \\text{Pa}\\) or \\(10.1325 \\text{dbar}\\). 
Sea pressure is the pressure measurement that is provided by the CTD (in hex 
and with some scaling) when the output format is set to deliver the data in 
"engineering units", as it is for the 37 IM CTDs. Therefore to back calculate 
absolute pressure in psia from dbar for the 37 IM data, use the 
Sea-Bird-specified conversion: 
$$
    \\text{P37IM} (\\text{psia}) = \\frac{\\text{P37IM} (\\text{dbar})}{0.689476} + 14.7
$$
The conversion between \\text{psia}\\ and \\text{dbar}\\ for the \\text{P16plusV2}\\ 
follows the TEOS-10 conventions using a "standard atmosphere":
$$
    \\text{P16plusV2} (\\text{dbar}) = 0.689475729 * \\text{P16plusV2} (\\text{psia}) - 10.1325
$$
Note that the raw data from the GPCTD make/model—the CTDs on board the gliders 
and autonomous underwater vehicles (AUVs)—are processed onboard the vehicles 
with proprietary software from the vehicle vendors. These data are presented 
already in decimal format in appropriate units (\\(\\text{°C}, \\text{Siemens/meter}, \\text{decibars}\\)), 
therefore processing raw hexadecimal data from the CTDGP is not included in the 
algorithm described in this document. 


## Conductivity

Water conductivity is one of the variables directly measured by CTD instruments and is calculated
from the raw hexadecimal data provided by the instrument. The \\(\\text{SBE 16plusV2}\\) model outputs
raw frequencies and voltages in hexadecimal. The \\(\\text{37IM}\\) model outputs engineering units in Hex.
This requires that the \\(\\text{L0}\\) conductivity data
product output from the CTD driver be converted from hexadecimal string to decimal
and scaled according to the CTD manual (different for each CTD make/model). Conversion and
scaling (described below) results in conductivity in \\(\\text{S/m}\\) based on factory calibration of
conductivity. Note that conductivity calibration scales with pressure and are handled in the
Conductivity Compensation DPS (DCN 1341-10030). See Sea-Bird App Note 10 (2008) for more
details. This is because conductivity cells measure the conductance of a specific geometry of
seawater within the conductivity cell and the ratio of the length of the cell to the cross sectional
area will be deformed slightly with changes to the ambient pressure surrounding the cell.

### Mathematical Theory
The CTD is received by OOI calibrated from Sea-Bird before its initial deployment. Post-initial
deployment calibration checks and/or recalibrations will be performed at Sea-Bird. The factory-calibrated
conductivity is converted to a hexadecimal string aboard the CTD instrument. 
This hexadecimal string, after it is separated from the rest of the data stream by the
CTD driver, is the L0 conductivity data product. The parsing instructions for the instrument
hexadecimal string can be found in Appendix A.
The L1 conductivity data product algorithm takes the L0 conductivity data product and converts it
into Siemens per meter \\(\\text{S/m}\\). Once the hexadecimal conductivity string is converted to decimal,
several simple calculations are necessary to produce the correct decimal floating representation
of the data in Siemens per meter \\(\\text{S/m}\\). The conversion differs by CTD make/model.
"""



import gsw

# Import Numpy and the GSW library
import numpy as np


def ctd_sbe16plus_tempwat(t0, a0, a1, a2, a3):
    """
    Calculates the OOI Level 1 Water Temperature data product using data from Sea-Bird Electronics CTD instruments.

    This function is for SBE 16Plus instruments and applies to CTDBP instruments (all series) and CTDPF instruments (series A and B).

    Parameters
    ----------
    t0 : float
        Raw temperature (TEMPWAT_L0) [counts].
    a0 : float
        Temperature calibration coefficient.
    a1 : float
        Temperature calibration coefficient.
    a2 : float
        Temperature calibration coefficient.
    a3 : float
        Temperature calibration coefficient.

    Returns
    -------
    float
        Sea water temperature (TEMPWAT_L1) [deg_C].

    Examples
    --------
    >>> ctd_sbe16plus_tempwat(123456, 1.0, 2.0, 3.0, 4.0)
    0.123456

    Historical References
    ---------------------
    [1] OOI (2012). Data Product Specification for Water Temperature. Document
       Control Number 1341-00010.
    """

    mv = (t0 - 524288) / 1.6e7
    r = (mv * 2.9e9 + 1.024e8)/(2.048e4 - mv * 2.0e5)
    t = 1 / (a0 + a1 * np.log(r) + a2 * np.log(r)**2 + a3 * np.log(r)**3) - 273.15
    return t


def ctd_sbe37im_tempwat_instrument_recovered(t0, a0, a1, a2, a3):
    """Calculates the OOI Level 1 Water Temperature data product from SBE 37IM 
    CTD instrument data. This function computes the sea water temperature using 
    raw temperature data recovered directly from the CTD instrument. 

    Parameters
    ----------
    t0 : array_like or float
        Raw temperature (TEMPWAT_L0) [counts] as recovered from the CTD itself.
    a0 : float
        Temperature calibration coefficient.
    a1 : float
        Temperature calibration coefficient.
    a2 : float
        Temperature calibration coefficient.
    a3 : float
        Temperature calibration coefficient.
    
    Returns
    -------
    t : ndarray or float
        Sea water temperature (TEMPWAT_L1) in degrees Celsius.

    Notes
    -----
    - Implemented for SBE 37IM instruments, all series.
    """


    t = 1 / (a0 + a1 * np.log(t0) + a2 * np.log(t0)**2 + a3 * np.log(t0)**3) - 273.15
    return t


def ctd_sbe37im_tempwat(t0):
    """
    Calculates the OOI Level 1 Water Temperature data product using data from Sea-Bird Electronics SBE 37IM CTD instruments.

    Parameters
    ----------
    t0 : float or array_like
        Raw temperature (TEMPWAT_L0) [counts].

    Returns
    -------
    float or ndarray
        Sea water temperature (TEMPWAT_L1) [deg_C].

    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Water Temperature. Document
        Control Number 1341-00010.
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00010_Data_Product_SPEC_TEMPWAT_OOI.pdf)
    """

    t = t0 / 10000.0 - 10.0
    return t


def ctd_sbe52mp_tempwat(t0):
    """
    Calculate sea water temperature from raw SBE 52MP CTD instrument counts
    for CTDPF instruments, C,K,L series.

    Parameters
    ----------
    t0 : array_like or float
        Raw temperature counts from the SBE 52MP CTD instrument (TEMPWAT_L0).

    Returns
    -------
    t : ndarray or float
        Sea water temperature in degrees Celsius (TEMPWAT_L1).

    Notes
    -----
    This function implements the OOI Level 1 Water Temperature data product
    calculation for SBE 52MP CTD instruments (CTDPF, C/K/L series).

    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Water Temperature.
    Document Control Number 1341-00010.
    """

    t = t0 / 10000.0 - 5.0
    return t


def ctd_sbe16plus_preswat(p0, t0, ptempa0, ptempa1, ptempa2,ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2,
        pa0, pa1, pa2, offset=0):
    """
    Calculates the OOI Level 1 Pressure (Depth) data product from raw sensor counts
    and calibration coefficients for Sea-Bird Electronics SBE 16Plus CTD instruments
    with a strain gauge pressure sensor.

    Parameters
    ----------
    p0 : float or array-like
        Raw pressure (PRESWAT_L0) [counts].
    t0 : float or array-like
        Raw temperature from pressure sensor thermistor [counts].
    ptempa0 : float
        Strain gauge pressure calibration coefficient.
    ptempa1 : float
        Strain gauge pressure calibration coefficient.
    ptempa2 : float
        Strain gauge pressure calibration coefficient.
    ptca0 : float
        Strain gauge pressure calibration coefficient.
    ptca1 : float
        Strain gauge pressure calibration coefficient.
    ptca2 : float
        Strain gauge pressure calibration coefficient.
    ptcb0 : float
        Strain gauge pressure calibration coefficient.
    ptcb1 : float
        Strain gauge pressure calibration coefficient.
    ptcb2 : float
        Strain gauge pressure calibration coefficient.
    pa0 : float
        Strain gauge pressure calibration coefficient.
    pa1 : float
        Strain gauge pressure calibration coefficient.
    pa2 : float
        Strain gauge pressure calibration coefficient.
    offset : float, optional
        Correction for Druck error [dbar]. Default is 0.

    Returns
    -------
    p : float or array-like
        Sea water pressure (PRESWAT_L1) [dbar].
    Notes
    -----
      - This data product is derived from SBE 16Plus instruments outfitted with
    a strain gauge pressure sensor. This is the default for most of the CTDBP
    instruments (the exceptions are series N and O) and for CTDPF instruments,
    series A and B.

    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Pressure (Depth). Document Control Number 1341-00020.
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
    Calculates the OOI Level 1 Pressure (Depth) data product for SBE 16Plus instruments
    outfitted with a digiquartz pressure sensor (CTDBP-N,O instruments only).

    Parameters
    ----------
    p0 : float or array_like
        Raw pressure (PRESWAT_L0) [counts].
    t0 : float or array_like
        Raw temperature from pressure sensor thermistor [counts].
    C1 : float
        Digiquartz pressure calibration coefficient.
    C2 : float
        Digiquartz pressure calibration coefficient.
    C3 : float
        Digiquartz pressure calibration coefficient.
    D1 : float
        Digiquartz pressure calibration coefficient.
    D2 : float
        Digiquartz pressure calibration coefficient.
    T1 : float
        Digiquartz pressure calibration coefficient.
    T2 : float
        Digiquartz pressure calibration coefficient.
    T3 : float
        Digiquartz pressure calibration coefficient.
    T4 : float
        Digiquartz pressure calibration coefficient.
    T5 : float
        Digiquartz pressure calibration coefficient.

    Returns
    -------
    float or ndarray
        Sea water pressure (PRESWAT_L1) [dbar].

    Notes
    -----
    - Applies to CTDBP-N,O instruments only.

    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Pressure (Depth). Document Control Number 1341-00020.
    OOI (2011). SeaBird 16Plus V2 User Manual. 1341-00020_PRESWAT Artifact.
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
    Calculates the OOI Level 1 Pressure (Depth) data product from SBE 37IM instrument data
    recovered directly from the CTD itself.

    Parameters
    ----------
    p0 : float or array_like
        Raw pressure (PRESWAT_L0) [counts] as recovered from the CTD itself.
    pt0 : float or array_like
        Raw temperature from pressure sensor thermistor [counts].
    ptempa0 : float
        Strain gauge pressure calibration coefficient.
    ptempa1 : float
        Strain gauge pressure calibration coefficient.
    ptempa2 : float
        Strain gauge pressure calibration coefficient.
    ptca0 : float
        Strain gauge pressure calibration coefficient.
    ptca1 : float
        Strain gauge pressure calibration coefficient.
    ptca2 : float
        Strain gauge pressure calibration coefficient.
    ptcb0 : float
        Strain gauge pressure calibration coefficient.
    ptcb1 : float
        Strain gauge pressure calibration coefficient.
    ptcb2 : float
        Strain gauge pressure calibration coefficient.
    pa0 : float
        Strain gauge pressure calibration coefficient.
    pa1 : float
        Strain gauge pressure calibration coefficient.
    pa2 : float
        Strain gauge pressure calibration coefficient.

    Returns
    -------
    p_dbar : float or ndarray
        Sea water pressure (PRESWAT_L1) [dbar].

    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Pressure (Depth). Document
        Control Number 1341-00020.
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
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
    Calculates the OOI Level 1 Pressure (Depth) data product using data from Sea-Bird Electronics SBE 37IM CTD instruments.

    Parameters
    ----------
    p0 : float or array_like
        Raw pressure (PRESWAT_L0) [counts].
    p_range_psia : float
        Pressure range calibration coefficient [psia].

    Returns
    -------
    float or ndarray
        Sea water pressure (PRESWAT_L1) [dbar].

    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Pressure (Depth). Document Control Number 1341-00020.
    (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
    1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
    """
    # compute pressure range in units of dbar
    p_range_dbar = (p_range_psia - 14.7) * 0.6894757

    # compute pressure in dbar and return
    p_dbar = p0 * p_range_dbar / (0.85 * 65536.0) - 0.05 * p_range_dbar
    return p_dbar


def ctd_glider_preswat(pr_bar):
    """
    Calculates the OOI Level 1 Pressure (Depth) data product for Seabird CTDs installed on gliders (CTDGV instruments).

    Parameters
    ----------
    pr_bar : float or array_like
        Sea water pressure value from glider [bar].

    Returns
    -------
    pr_dbar : float or ndarray
        Sea water pressure (PRESWAT_L1) [dbar].

    Historical References
    ---------------------
    OOI (2015). Data Product Specification for Coastal Glider Data Products (version 1-03).
        Document Control Number 1341-00020.
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
        1341-20001_Data_Product_SPEC_CSGLIDR_OOI.docx)
    """

    pr_dbar = pr_bar * 10.0
    return pr_dbar


def ctd_sbe52mp_preswat(p0):
    """
    Calculates the OOI Level 1 Pressure (Depth) data product using data from Sea-Bird Electronics SBE 52MP CTD instruments.

    Parameters
    ----------
    p0 : float or array_like
        Raw pressure (PRESWAT_L0) [counts].

    Returns
    -------
    p_dbar : float or ndarray
        Sea water pressure (PRESWAT_L1) [dbar].

    Notes
    -----
    - This data product is derived from SBE 52MP instruments and applies to CTDPF instruments, C/K/L series.

    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Pressure (Depth). Document Control Number 1341-00020.
    (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
    1341-00020_Data_Product_SPEC_PRESWAT_OOI.pdf)
    """

    p_dbar = p0 / 100.0 - 10.0
    return p_dbar


def ctd_sbe16plus_condwat(c0, t1, p1, g, h, i, j, cpcor, ctcor):
    """
    Calculates the OOI Level 1 Conductivity core data product using data from Sea-Bird Electronics SBE 16Plus CTD instruments.

    Parameters
    ----------
    c0 : float or array_like
        Raw conductivity (CONDWAT_L0) [counts].
    t1 : float or array_like
        Sea water temperature (TEMPWAT_L1) [deg_C].
    p1 : float or array_like
        Sea water pressure (PRESWAT_L1) [dbar].
    g : float
        Conductivity calibration coefficient.
    h : float
        Conductivity calibration coefficient.
    i : float
        Conductivity calibration coefficient.
    j : float
        Conductivity calibration coefficient.
    cpcor : float
        Conductivity calibration coefficient.
    ctcor : float
        Conductivity calibration coefficient.

    Returns
    -------
    c : float or ndarray
        Sea water conductivity (CONDWAT_L1) [S m-1].

    Notes
    -----
    - This data product is derived from SBE 16Plus instruments and applies to CTDBP instruments (all series) and CTDPF instruments (series A and B).

    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Conductivity. Document Control Number 1341-00030.
    (See: Company Home >> OOI >> Controlled >> 1000 System Level >>
    1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)
    """
    # convert raw conductivity measurement to frequency
    f = (c0 / 256.0) / 1000.0

    # calculate conductivity [S m-1]
    c = (g + h * f**2 + i * f**3 + j * f**4) / (1 + ctcor * t1 + cpcor * p1)
    return c


def ctd_sbe37im_condwat_instrument_recovered(c0, t1, p1, g, h, i, j, cpcor, ctcor, wbotc):
    """
    Calculates the OOI Level 1 Conductivity core data product using data from Sea-Bird Electronics SBE 37IM CTD instruments.

    Parameters
    ----------
    c0 : float or array_like
        Raw conductivity (CONDWAT_L0) [counts] as recovered from the CTD itself.
    t1 : float or array_like
        Sea water temperature (TEMPWAT_L1) [deg_C].
    p1 : float or array_like
        Sea water pressure (PRESWAT_L1) [dbar].
    g : float
        Conductivity calibration coefficient.
    h : float
        Conductivity calibration coefficient.
    i : float
        Conductivity calibration coefficient.
    j : float
        Conductivity calibration coefficient.
    cpcor : float
        Conductivity calibration coefficient.
    ctcor : float
        Conductivity calibration coefficient.
    wbotc : float
        Conductivity calibration coefficient.

    Returns
    -------
    c : float or ndarray
        Sea water conductivity (CONDWAT_L1) [S m-1].

    Notes
    -----
    - This data product is derived from SBE 37IM instruments, all series, and specifically processes data recovered directly from the CTD itself.

    Historical References
    ---------------------
    As of June 2016 the following DPS does not contain this specification.

    OOI (2012). Data Product Specification for Conductivity. Document Control Number 1341-00030.
    https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)
    """
    # convert raw conductivity measurement to frequency
    f = (c0 / 256.0) / 1000.0 * np.sqrt(1.0 + wbotc * t1)

    # calculate conductivity [S m-1]
    c = (g + h * f**2 + i * f**3 + j * f**4) / (1 + ctcor * t1 + cpcor * p1)
    return c


def ctd_sbe37im_condwat(c0):
    """
    Calculates the OOI Level 1 Conductivity core data product using data from Sea-Bird Electronics SBE 37IM CTD instruments.

    Parameters
    ----------
    c0 : float or array_like
        Raw conductivity (CONDWAT_L0) [counts].

    Returns
    -------
    c : float or ndarray
        Sea water conductivity (CONDWAT_L1) [S m-1].

    Notes
    -----
    - This data product is derived from SBE 37IM instruments and applies to CTDMO instruments, all series, specifically for telemetered and recovered_host (not recovered_instrument) data.
    Historical References
    ---------------------
    OOI (2012). Data Product Specification for Conductivity. Document Control Number 1341-00030.
    (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)
    """

    c = c0 / 100000.0 - 0.5
    return c


def ctd_sbe52mp_condwat(c0):
    """
    Calculates the OOI Level 1 Conductivity core data product using data from Sea-Bird Electronics SBE 52MP CTD instruments.

    Parameters
    ----------
    c0 : float or array_like
        Raw conductivity (CONDWAT_L0) [counts].

    Returns
    -------
    c_S_m : float or ndarray
        Sea water conductivity (CONDWAT_L1) [S m-1].

    Notes
    -----
    - This data product is derived from SBE 52MP instruments and applies to CTDPF instruments, C/K/L series.

    References
    ----------
    OOI (2012). Data Product Specification for Conductivity. Document Control Number 1341-00030.
    (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00030_Data_Product_SPEC_CONDWAT_OOI.pdf)
    """

    c_mmho_cm = c0 / 10000.0 - 0.5
    c_S_m = 0.1 * c_mmho_cm
    return c_S_m


def ctd_pracsal(c, t, p):
    """
    Calculates Practical Salinity (PSS-78) from conductivity, temperature, and pressure.
    This function computes the OOI Level 2 Practical Salinity core data product using the
    Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version 3.0, based on input from
    conductivity, temperature, and depth (CTD) instruments.

    Parameters
    ----------
    c : float or array-like
        Sea water conductivity (CONDWAT_L1) [S m-1].
    t : float or array-like
        Sea water temperature (TEMPWAT_L1) [deg_C].
    p : float or array-like
        Sea water pressure (PRESWAT_L1) [dbar].

    Returns
    -------
    SP : float or array-like
        Practical salinity, PSS-78, (PRACSAL_L2) unitless.

    Historical References
    ---------------------
      - OOI (2012). Data Product Specification for Salinity. Document Control Number 1341-00040.
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00040_Data_Product_SPEC_PRACSAL_OOI.pdf)
      - Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version 3.0
    """
    # Convert L1 Conductivity from S/m to mS/cm
    C10 = c * 10.0

    # Calculate the Practical Salinity (PSS-78) [unitless]
    SP = gsw.SP_from_C(C10, t, p)
    return SP


def ctd_density(SP, t, p, lat, lon):
    """
    Calculates the OOI Level 2 Density core data product for the conductivity, 
    temperature and depth (CTD) family of instruments.

    Parameters
    ----------
    SP : float or array_like
        Practical salinity PSS-78 (PRACSAL_L2) [unitless].
    t : float or array_like
        Sea water temperature (TEMPWAT_L1) [deg_C].
    p : float or array_like
        Sea water pressure (PRESWAT_L1) [dbar].
    lat : float or array_like
        Latitude where input data was collected [decimal degree].
    lon : float or array_like
        Longitude where input data was collected [decimal degree].

    Returns
    -------
    rho : float or ndarray
        Sea water density (DENSITY_L2) [kg m-3].

    Notes
    -----
      - Uses the Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version 3.0.

    Historical References
    ---------------------
    - OOI (2012). Data Product Specification for Density. Document Control Number 1341-00050.
        (See: Company Home >> OOI >> Controlled >> 1000 System Level >> 1341-00050_Data_Product_SPEC_DENSITY_OOI.pdf)
    - Thermodynamic Equations of Seawater - 2010 (TEOS-10) Version 3.0.
    """
    # Calculate the density [kg m-3]
    sa = gsw.SA_from_SP(SP, p, lon, lat)  # absolute salinity
    ct = gsw.CT_from_t(sa, t, p)  # conservative temperature
    rho = gsw.rho(sa, ct, p)  # density
    # rho = gsw.ctd_density(SP, t, p, lat, lon)

    return rho
