#!/usr/bin/env python
"""
Module containing PLIMS (Plankton Imaging System) instrument family data
processing functions for the Ocean Observatories Initiative.
"""
import numpy as np


def plims_flow_rate(run_fast, run_fast_factor):
    """
    Compute the PLIMS sample flow rate (FLOWRT_L1).

    The Imaging FlowCytobot injects sample at a standard rate of
    0.25 mL/min. When the instrument is set to run fast, the standard
    rate is multiplied by the run fast factor.

    Parameters
    ----------
    run_fast : array_like
        Run fast setting (runSampleFast, from the IFCB header file);
        True or nonzero when the sampling speed is increased by the run
        fast factor.
    run_fast_factor : array_like
        Factor multiplied by the standard injection rate when run_fast
        is set (RunFastFactor, from the IFCB header file).

    Returns
    -------
    flow_rate : ndarray
        Sample flow rate (FLOWRT_L1) [mL min-1].
    """
    run_fast = np.atleast_1d(run_fast).astype(bool)
    base_flow_rate = 0.25

    return base_flow_rate * np.where(run_fast, run_fast_factor, 1)
