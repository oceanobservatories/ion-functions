#!/usr/bin/env python
"""
@package ion_functions.data.hyd_functions
@file ion_functions/data/hyd_functions.py
@author Christopher Wingard
@brief Module containing Hydrophone instrument family related functions
"""

import numpy as np


def hyd_bb_acoustic_pwaves(wav, gain):
    """
    Description:

        Calculates the OOI Level 1 (L1) Broadband Acoustic Pressure Waves core
        data product (HYDAPBB), using data from Broadband Hydrophone (HYDBB)
        instruments. The HYDBB instrument senses passive acoustic pressure
        waves from 5 Hz to 100 kHz, at 24-bit resolution.

    Implemented by:

        2014-05-16: Christopher Wingard. Initial Code
        2023-08-15: Samuel Dahlberg. Removed use of Numexpr library.
                    Changed local variable names to follow naming convention.

    Usage:

        tsv = hyd_bb_acoustic_pwaves(wav, gain)

            where

        tsv = time-series voltage compensated for external gain and wav format
            scaling [Volts] (HYDAPBB_L1)
        wav = raw time-series voltage [Volts] (HYDAPBB_L0)
        gain = external gain setting [dB]

    References:

        OOI (2013). Data Product Specification for Acoustic Pressure Waves.
            Document Control Number 1341-00820.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00820_Data_Product_SPEC_HYDAPBB_OOI.pdf)
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
    Description:

        Calculates the OOI Level 1 (L1) Low Frequency Acoustic Pressure Waves
        core data product (HYDAPLF), using data from the Low Frequency
        Hydrophone (HYDLF) instruments.

    Implemented by:

        2014-07-09: Christopher Wingard. Initial Code.
        2023-08-15: Samuel Dahlberg. Removed use of Numexpr library.

    Usage:

        hydaplf = hyd_lf_acoustic_pwaves(counts, gain)

            where

        hydaplf = time-series of low frequency acoustic pressure waves [V]
            (HYDAPLF_L1)
        raw = raw time-series digitizied in counts [counts] (HYDAPLF_L0)
        gain = Gurlap DM24 fixed gain bit weight [uV/count]

    References:

        OOI (2013). Data Product Specification for Low Frequency Acoustic
            Pressure Waves. Document Control Number 1341-00821.
            https://alfresco.oceanobservatories.org/ (See: Company Home >>
            OOI >> Controlled >> 1000 System Level >>
            1341-00821_Data_Product_SPEC_HYDAPLF_OOI.pdf)
    """
    # apply the gain correction to convert the signal from counts to V
    gain = gain * 1.0e-6
    hydaplf = raw * gain
    return hydaplf
