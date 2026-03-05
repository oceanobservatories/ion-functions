#!/usr/bin/env python

"""
@package ion_functions.test.test_phsen_h_functions
@file ion_functions/test/test_phsen_h_functions.py
@author Samuel Dahlberg
@brief Unit tests for phsen_h_functions module
"""
import numpy as np
from nose.plugins.attrib import attr
from ion_functions.data import phsen_h_functions
from ion_functions.test.base_test import BaseUnitTestCase


@attr('UNIT', group='func')
class TestPhsenHFunctions(BaseUnitTestCase):

    def test_temperature_raw_conversion(self):
        """
        Tests the temperature_raw_conversion function with test units provided by SeaBird.
        """

        temp_counts = np.array([2864.51635696, 2864.61073548, 2864.79005751, 2864.89387722])
        a0 = 1.28015621e-003
        a1 = 2.58367774e-004
        a2 = -1.39527596e-008
        a3 = 1.39024630e-007

        expected = [20.4459, 20.4451, 20.4436, 20.4427]

        result = phsen_h_functions.temperature_raw_conversion(temp_counts, a0, a1, a2, a3)
        # print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=4)

    def test_pressure_raw_conversion(self):
        """
        Tests the pressure_raw_conversion function with test units provided by SeaBird.
        """

        pres_counts = np.array([533539, 533538, 533540, 533537])
        compensation_voltage = np.array([20625.0, 20626.0, 20626.0, 20626.0]) / 13107.0
        ptempa0 = -6.28624239e001
        ptempa1 = 5.41620689e001
        ptempa2 = -2.96026659e-001
        ptca0 = 5.24108391e005
        ptca1 = 5.47371611e000
        ptca2 = -1.53365246e-001
        ptcb0 = 2.59008750e001
        ptcb1 = 7.75000000e-004
        ptcb2 = 0.00000000e000
        pa0 = 6.51546669e-002
        pa1 = 1.54424116e-003
        pa2 = 6.13653149e-012

        expected = [-0.105, -0.106, -0.104, -0.107]

        result = phsen_h_functions.pressure_raw_conversion(pres_counts, compensation_voltage, ptempa0, ptempa1, ptempa2,
                                                           ptca0, ptca1, ptca2, ptcb0, ptcb1, ptcb2, pa0, pa1, pa2)
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=3)

    def test_conductivity_raw_conversion(self):
        """
        Tests the conductivity_raw_conversion function with test units provided by SeaBird.
        """

        cond_counts = np.array([6319.22, 6319.22, 6319.22, 6319.22])
        temperature = np.array([20.4459, 20.4451, 20.4436, 20.4427])
        pressure = np.array([-0.105, -0.106, -0.104, -0.107])
        wbotc = -2.319375e-007
        g = -9.996553e-001
        h = 1.213221e-001
        i = -1.856730e-004
        j = 2.812110e-005
        ctcor = 3.250000e-006
        cpcor = -9.570000e-008

        # Not the correct expected ouput, seabird did not provide proper test units for this code.
        expected = [3.84275, 2.715786, 2.715857, 2.715846]

        result = phsen_h_functions.conductivity_raw_conversion(cond_counts, temperature, pressure, wbotc, g, h, i, j,
                                                               ctcor, cpcor)
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=3)

    def test_internal_temperature(self):
        """
        Tests the internal_temperature function with test units provided by SeaBird.
        """

        temp_counts = np.array([25616, 25600])

        expected = np.array([21.8335, 21.7906])

        result = phsen_h_functions.internal_temperature(temp_counts)
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=4)

    def test_internal_humidity(self):
        """
        Tests the internal_humidity function with test units provided by SeaBird.
        """

        humidity_counts = np.array([24096, 24160])
        temperature = np.array([21.8335, 21.7906])

        expected = np.array([39.4845, 39.6001])

        result = phsen_h_functions.internal_humidity(humidity_counts, temperature)
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=4)

    def test_dissolved_oxygen(self):
        """
        Tests the dissolved_oxygen function with test units provided by SeaBird.
        """

        raw_oxygen_phase = np.array([31.06, 31.66, 32.59, 33.92, 34.82, 35.44, 35.44, 35.44])
        thermistor = np.array([30, 26, 20, 12, 6, 2, 2, 2])
        pressure = np.array([0, 0, 0, 0, 0, 0, 1000, 100])
        salinity = np.array([0, 0, 0, 0, 0, 0, 0, 35])
        c0 = 1.0355e-1
        c1 = 4.4295e-3
        c2 = 6.0011e-5
        coeff_e = 1.1e-2
        a0 = 1.0513
        a1 = -1.5e-3
        a2 = 4.1907e-1
        b0 = -2.5004e-1
        b1 = 1.6524
        therm_ta0 = 7.059180e-4
        therm_ta1 = 2.504670e-4
        therm_ta2 = 7.402389e-7
        therm_ta3 = 9.756123e-8
        lat = 44.65685
        lon = -124.09663

        # The expected outut is without the conversion to umol/kg. Test currently outputs false, and you have to
        # check the expected numbers with the calculations before umol/kg conversion.
        expected = np.array([0.706, 0.74, 0.799, 0.892, 1.005, 1.095, 1.1398, 0.8647])

        result = phsen_h_functions.dissolved_oxygen(raw_oxygen_phase, thermistor, pressure, salinity, c0, c1, c2,
                                                    coeff_e, a0, a1, a2, b0, b1, therm_ta0, therm_ta1, therm_ta2,
                                                    therm_ta3, lat, lon, thermistor_units="C")
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=4)

    def test_dissolved_oxygen_raw(self):
        """
        Tests the dissolved_oxygen function with test units provided by SeaBird.
        """

        raw_oxygen_phase = np.array([31.06, 31.66, 32.59, 33.92, 34.82, 35.44, 40.8, 25.64])
        thermistor = np.array([0.6, 0.5, 0.4, 0.35, 0.3, 0.25, .86, .92])
        pressure = np.array([0, 0, 0, 0, 0, 0, 79.23, 195.14])
        salinity = np.array([0, 0, 0, 0, 0, 0, 33.13, 33.85])

        c0 = 1.0355e-1
        c1 = 4.4295e-3
        c2 = 6.0011e-5
        coeff_e = 1.1e-2
        a0 = 1.0513
        a1 = -1.5e-3
        a2 = 4.1907e-1
        b0 = -2.5004e-1
        b1 = 1.6524
        therm_ta0 = 7.059180e-4
        therm_ta1 = 2.504670e-4
        therm_ta2 = 7.402389e-7
        therm_ta3 = 9.756123e-8
        lat = 44.65685
        lon = -124.09663

        # The expected outut is without the conversion to umol/kg. Test currently outputs false, and you have to
        # check the expected numbers with the calculations before umol/kg conversion.
        expected = np.array([0.93, 0.688, 0.459, 0.304, 0.206, 0.137, 0.00])

        result = phsen_h_functions.dissolved_oxygen(raw_oxygen_phase, thermistor, pressure, salinity, c0, c1, c2,
                                                    coeff_e, a0, a1, a2, b0, b1, therm_ta0, therm_ta1, therm_ta2,
                                                    therm_ta3, lat, lon)
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=2)

    def test_convert_sbe63_thermistor(self):
        """
        Tests the convert_sbe63_thermistor function with test units provided by SeaBird.
        """

        raw_temperature = np.array([1.12015, 1.12015, 1.12016, 1.12016, 0.99562, 0.82934, 0.64528])
        therm_ta0 = 7.059180e-4
        therm_ta1 = 2.504670e-4
        therm_ta2 = 7.402389e-7
        therm_ta3 = 9.756123e-8

        expected = np.array([2.0002, 2.0002, 1.9999, 1.9999, 5.9998, 12.0, 19.9998])

        result = phsen_h_functions.convert_sbe63_thermistor(raw_temperature, therm_ta0, therm_ta1, therm_ta2,
                                                            therm_ta3)
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=4)

    def test_convert_ph_voltage_counts(self):
        """
        Tests the convert_ph_voltage_counts function with test units provided by SeaBird.
        """

        internal_ph_counts = np.array([5105334, 5105384, 5105350, 5105505, 5105347])

        expected = np.array([-0.978492, -0.978477, -0.978487, -0.978441, -0.978488])

        result = phsen_h_functions.convert_ph_voltage_counts(internal_ph_counts)
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=4)

    def test_ph_total(self):
        """
        Tests the ph_total function.
        """

        vrs_ext = np.array([-0.885081, -0.885081, -0.885081])
        degc = np.array([23.4169, 23.4169, 23.4169])
        psu = np.array([34.812, 34.812, 34.812])
        dbar = np.array([100, 100, 100])
        k0 = -1.361736
        k2 = -1.07686e-3
        f = [[-8.31842e-6, -7.47152e-9, 1.91485e-11, -1.39273e-14, 4.48185e-18, -5.42588e-22]]

        expected = np.array([7.9394, 7.9394, 7.9394])

        result = phsen_h_functions.ph_total(vrs_ext, degc, psu, dbar, k0, k2, f)
        print(result)
        for x, y in zip(result, expected):
            self.assertAlmostEqual(x, y, places=4)
