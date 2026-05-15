#!/usr/bin/env python
"""
Module containing Hydrophone (HYD) instrument family data processing
functions for the Ocean Observatories Initiative. Converts raw L0 data
from Broadband Hydrophone (HYDBB) and Low Frequency Hydrophone (HYDLF)
instruments into L1 acoustic pressure wave data products.
"""

import numpy as np


def hyd_bb_acoustic_pwaves(wav, gain):
    """
    Compute the L1 Broadband Acoustic Pressure Waves product (HYDAPBB)
    from the Ocean Sonics icListen HF Broadband Hydrophone (HYDBB).

    Parameters
    ----------
    wav : array_like
        Raw L0 time-series from the HYDBB instrument encoded as
        normalized WAV floating-point values in the range [-1, 1]
        (HYDAPBB_L0). Shape (n_records, n_samples) or (n_samples,).
    gain : float or array_like
        External gain applied by the instrument signal conditioning
        chain [dBV]. Scalar or shape (n_records,).

    Returns
    -------
    tsv : ndarray
        L1 time-series voltage compensated for external gain
        (HYDAPBB_L1) [V]. Shape (n_records, n_samples).

    Notes
    -----
    The WAV format encodes samples as normalized values in [-1, 1].
    Multiplying by 3 recovers the physical voltage using the icListen
    HF full-scale range of +/-3 V. The gain in dBV is converted to a
    linear factor via 10^(gain/20) and divided out.

    The frequency-dependent OCVR (Open Circuit Voltage Response) of
    the hydrophone cannot be applied to a broadband time-series without
    prior frequency-domain filtering. HYDAPBB_L1 is therefore
    gain-compensated voltage, not calibrated pressure.
    """
    # shape inputs to correct dimensions
    wav = np.atleast_2d(wav)
    n_rec = wav.shape[0]

    if np.isscalar(gain) is True:
        gain = np.tile(gain, (n_rec, 1))
    else:
        gain = np.reshape(gain, (n_rec, 1))

    # Convert the gain from dB to a linear value
    gain = 10**(gain/20.)

    # convert the broadband acoustic pressure wave data to Volts
    volts = wav * 3.

    # and correct for the gain
    tsv = volts / gain
    return tsv


def hyd_lf_acoustic_pwaves(raw, gain=3.2):
    """
    Compute the L1 Low Frequency Acoustic Pressure Waves product
    (HYDAPLF) from the Low Frequency Hydrophone (HYDLF).

    Parameters
    ----------
    raw : array_like
        Raw L0 time-series digitized by the Guralp DM24 digitizer
        [counts] (HYDAPLF_L0).
    gain : float, optional
        Guralp DM24 fixed bit weight [uV/count]. Default is 3.2.

    Returns
    -------
    hydaplf : ndarray
        L1 time-series of low frequency acoustic pressure waves
        (HYDAPLF_L1) [V].

    Notes
    -----
    The bit weight converts raw counts to microvolts; the result is
    scaled by 1e-6 to return Volts. The frequency-dependent OCVR of
    the hydrophone cannot be applied to a broadband time-series without
    prior frequency-domain filtering. HYDAPLF_L1 is therefore
    gain-compensated voltage, not calibrated pressure.
    """
    # apply the gain correction to convert the signal from counts to V
    gain = gain * 1.0e-6
    hydaplf = raw * gain
    return hydaplf
