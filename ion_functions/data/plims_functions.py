#!/usr/bin/env python
"""
@package ion_functions.data.plims_functions
@file ion_functions/data/plims_functions.py
@author Joffrey Peters
@brief Module containing PLIMS/IFCB related data-calculations.
"""

"""
#...................................................................................
    The PLIMS/IFCB instrument outputs two data files for each sample: 
     * a header file (*.hdr) and
    * a data file from the ADCs (*.adc)
    There are additionally several values of interest to be calculated for each
    sample which combine values from the two files.
#...................................................................................
"""
import numpy as np


def plims_flow_rate(run_fast, run_fast_factor, sample_speed, sample_volume):
    run_fast = np.atleast_1d(run_fast).astype(bool)
    output = np.zeros(run_fast.shape)
    output[run_fast] = sample_volume[run_fast] / (sample_speed[run_fast] / run_fast_factor[run_fast])
    output[~run_fast] = sample_volume[~run_fast] / sample_speed[~run_fast]

    return output

