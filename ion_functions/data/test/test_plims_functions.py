#!/usr/bin/env python

"""
@package ion_functions.test.test_plims_functions
@file ion_functions/test/test_plims_functions.py
@author Joffrey Peters
@brief Unit tests for plims_functions module
"""

from ion_functions.data import plims_functions as plimsfunc
from ion_functions.test.base_test import BaseUnitTestCase
from nose.plugins.attrib import attr
import numpy as np

@attr('UNIT', group='func')
class TestPlimsFunctionsUnit(BaseUnitTestCase):
    # Test the flow_rate function
    def test_plims_flow_rate(self):
        """
        Test the plims_functions.flow_rate function
        """
        # Test the flow_rate function with run_fast = True
        run_fast = np.array([0, 0, 0, 1, 1])
        run_fast_factor = np.array([5.0, 5.0, 5.0, 5.0, 5.0])
        sample_speed = 20.0
        sample_volume = 5.0
        expected = np.array([0.25, 0.25, 0.25, 1.25, 1.25])
        result = plimsfunc.plims_flow_rate(run_fast, run_fast_factor, sample_speed, sample_volume)
        np.testing.assert_array_equal(result, expected)

        # Test the flow_rate function with run_fast = False
        run_fast = False
        run_fast_factor = 2.0
        sample_speed = 20.0
        sample_volume = 5.0
        expected = 0.25
        result = plimsfunc.plims_flow_rate(run_fast, run_fast_factor, sample_speed, sample_volume)
        self.assertEqual(result, expected)

