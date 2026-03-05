#!/usr/bin/env python

"""
@package ion_functions.test.test_plims_functions
@file ion_functions/test/test_plims_functions.py
@author Joffrey Peters
@brief Unit tests for plims_functions module
"""

from ion_functions.data import plims_functions as plimsfunc
from ion_functions.test.base_test import BaseUnitTestCase


@attr('UNIT', group='func')
class TestPlimsFunctionsUnit(BaseUnitTestCase):
    # Test the flow_rate function
    def test_plims_flow_rate(self):
        """
        Test the plims_functions.flow_rate function
        """
        # Test the flow_rate function with run_fast = True
        run_fast = True
        run_fast_factor = 2.0
        sample_speed = 0.5
        sample_volume = 1.0
        expected = 4.0
        result = plimsfunc.plims_flow_rate(run_fast, run_fast_factor, sample_speed, sample_volume)
        self.assertEqual(result, expected)

        # Test the flow_rate function with run_fast = False
        run_fast = False
        run_fast_factor = 2.0
        sample_speed = 0.5
        sample_volume = 1.0
        expected = 2.0
        result = plimsfunc.plims_flow_rate(run_fast, run_fast_factor, sample_speed, sample_volume)
        self.assertEqual(result, expected)

