#!/usr/bin/env python
"""
Module containing Ocean Bottom Seismometer (OBS) instrument family data
processing functions for the Ocean Observatories Initiative. Converts
raw L0 data from Broadband (OBSBB) and Short Period (OBSSP) seismometers
into L1 ground motion data products.
"""
from ion_functions import deprecated


@deprecated
def obs_bb_ground_velocity(raw, gain=3.2, sensitivity=1500.):
    """
    Compute broadband ground velocity (GRNDVEL_L1) from raw OBS counts.

    Parameters
    ----------
    raw : array_like
        Raw time-series digitized in counts [counts] (GRNDVEL_L0).
    gain : float, optional
        Guralp DM24 fixed gain bit weight [uV/count]. Default is 3.2.
    sensitivity : float, optional
        Guralp CMG-1T sensor sensitivity [V/(m/s)]. Default is 1500.0.

    Returns
    -------
    grndvel : array_like
        Time-series broadband ground velocity [m/s] (GRNDVEL_L1).

    Notes
    -----
    Applied to OBSBB instruments deployed on the Regional Cabled Array.
    The gain is converted from uV/count to V/count before applying the
    calibration. Full algorithm details are in DPS DCN 1341-00090.

    History
    -------
    2014-07-09: Christopher Wingard. Initial Code
    2023-08-15: Samuel Dahlberg. Removed use of Numexpr library.
    2025-05-15: Christopher Wingard. Converted to NumPy docstring format;
        updated documentation.
    """
    # scale the gain and sensitivity ...
    gain = gain * 1.0e-6
    sense = 2. * sensitivity

    # ... and calculate the broadband ground velocity
    grndvel = raw * (gain / sense)
    return grndvel


@deprecated
def obs_bb_ground_acceleration(raw, gain=3.2, sensitivity=0.508):
    """
    Compute broadband ground acceleration (GRNDACC_L1) from raw OBS counts.

    Parameters
    ----------
    raw : array_like
        Raw time-series digitized in counts [counts] (GRNDACC_L0).
    gain : float, optional
        Guralp DM24 fixed gain bit weight [uV/count]. Default is 3.2.
    sensitivity : float, optional
        Guralp CMG-5T sensor sensitivity [V/(m/s^2)]. Default is 0.508.

    Returns
    -------
    grndacc : array_like
        Time-series broadband ground acceleration [m/s^2] (GRNDACC_L1).

    Notes
    -----
    Applied to OBSBB instruments deployed on the Regional Cabled Array.
    The gain is converted from uV/count to V/count before applying the
    calibration. Full algorithm details are in DPS DCN 1341-00100.

    History
    -------
    2014-07-09: Christopher Wingard. Initial Code
    2023-08-15: Samuel Dahlberg. Removed use of Numexpr library.
    2025-05-15: Christopher Wingard. Converted to NumPy docstring format;
        updated documentation.
    """
    # scale the gain and sensitivity ...
    gain = gain * 1.0e-6
    sense = 2. * sensitivity

    # ... and calculate the broadband ground acceleration
    grndacc = raw * (gain / sense)
    return grndacc


@deprecated
def obs_sp_ground_velocity(raw, gain=2.84, sensitivity=1200.):
    """
    Compute short period ground velocity (SGRDVEL_L1) from raw OBS counts.

    Parameters
    ----------
    raw : array_like
        Raw time-series digitized in counts [counts] (SGRDVEL_L0).
    gain : float, optional
        Guralp DM24 fixed gain bit weight [uV/count]. Default is 2.84.
    sensitivity : float, optional
        Guralp CMG-6T sensor sensitivity [V/(m/s)]. Default is 1200.0.

    Returns
    -------
    sgrdvel : array_like
        Time-series short period ground velocity [m/s] (SGRDVEL_L1).

    Notes
    -----
    Applied to OBSSP instruments deployed on the Regional Cabled Array
    at Hydrate Ridge and Axial Volcano. The gain is converted from
    uV/count to V/count before applying the calibration. Full algorithm
    details are in DPS DCN 1341-00110.

    History
    -------
    2014-07-09: Christopher Wingard. Initial Code
    2023-08-15: Samuel Dahlberg. Removed use of Numexpr library.
    2025-05-15: Christopher Wingard. Converted to NumPy docstring format;
        updated documentation.
    """
    # scale the gain and sensitivity ...
    gain = gain * 1.0e-6
    sense = 2. * sensitivity

    # ... and calculate the short period ground velocity
    sgrdvel = raw * (gain / sense)
    return sgrdvel
