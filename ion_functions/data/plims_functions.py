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


def plims_flow_rate(run_fast, run_fast_factor):
    """
    Description:

        flow rate calculation for PLIMS/IFCB

    Implemented by:

        2025-02-04: Joffrey Peters. Initial code.
        2026-10-30: Samuel Dahlberg. Updated to remove sample speed and sample volume from calculation.
                    Redmine Ticket #16237

    Usage:

        output = plims_flow_rate(run_fast, run_fast_factor)

            where

        output = Standard injection rate 0.25 mL/min unless run fast with factor differing from 1

        run_fast = True when sampling speed is increased by the Run Fast Factor
        run_fast_factor = Factor multiplied by syringe base flow rate to obtain the fast sample sampling speed
    """

    run_fast = np.atleast_1d(run_fast).astype(bool)
    base_flow_rate = 0.25

    return base_flow_rate * np.where(run_fast, run_fast_factor, 1)

