"""
@package ion_functions.test.adcp_functions
@file ion_functions/test/test_adcp_functions.py
@author Christopher Wingard, Russell Desiderio, Craig Risien
@brief Unit tests for adcp_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np

from ion_functions.data import adcp_functions as af
from ion_functions.data.adcp_functions import ADCP_FILLVALUE
from ion_functions.data.generic_functions import SYSTEM_FILLVALUE


@attr('UNIT', group='func')
class TestADCPFunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        """
        Implemented by:
            2014-02-06: Christopher Wingard. Initial Code.
            2015-06-12: Russell Desiderio. Changed raw beam data to type int. This
                        change did not affect any previously written unit tests.
            2019-08-13: Christopher Wingard. Adds functionality to compute a 3-beam solution
                        and cleans up syntax used in the functions. Removes system fill value tests
                        and integer tests (uneccessary and redundant). Cleans up several other
                        tests and in general cleans and tightens up this body of unit tests.

        """
        # set test inputs -- values from DPS (reset to integers, matches format from instrument)
        self.b1 = np.array([[-30, -295, -514, -234, -188, 203, -325, 305, -204, -294]])
        self.b2 = np.array([[180, -132, 213, 309, 291, 49, 188, 373, 2, 172]])
        self.b3 = np.array([[-398, -436, -131, -473, -443, 188, -168, 291, -179, 8]])
        self.b4 = np.array([[-216, -605, -92, -58, 484, -5, 338, 175, -80, -549]])

        self.pg1 = np.array([[100, 100, 26,  24, 25, 25,  24, 26, 100, 100]])
        self.pg2 = np.array([[100, 100, 90, 80, 60, 50, 100,  24, 25, 26]])
        self.pg3 = np.array([[26, 100, 100, 90, 80, 70, 60, 50,  24, 25]])
        self.pg4 = np.array([[25, 26, 100,  24, 25, 26, 100, 100, 100,  24]])

        self.echo = np.array([[0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250]])

        self.sfactor = 0.45
        self.heading = 9841  # units are centidegrees
        self.pitch = 69  # units are centidegrees
        self.roll = -254  # units are centidegrees
        self.orient = 1
        self.lat = 50.0000
        self.lon = -145.0000
        self.depth = 0.0
        self.ntp = 3545769600.0    # May 12, 2012

        # set expected results -- velocity profiles in earth coordinates (values in DPS)
        self.u = np.array([[0.2175, -0.2814, -0.1002, 0.4831, 1.2380,
                            -0.2455, 0.6218, -0.1807, 0.0992, -0.9063]]) * 1000.
        self.v = np.array([[-0.3367, -0.1815, -1.0522, -0.8676, -0.8919,
                            0.2585, -0.8497, -0.0873, -0.3073, -0.5461]]) * 1000.
        self.w = np.array([[0.1401,  0.3977,  0.1870,  0.1637,  0.0091,
                            -0.1290,  0.0334, -0.3017, 0.1384, 0.1966]]) * 1000.
        self.e = np.array([[0.789762, 0.634704, -0.080630, 0.626434, 0.064090,
                            0.071326, -0.317352, 0.219148, 0.054787, 0.433129]]) * 1000.

        # set expected results -- magnetic variation correction applied
        self.u_cor = np.array([[0.1099, -0.3221, -0.4025, 0.2092, 0.9243,
                                -0.1595, 0.3471, -0.1983, 0.0053, -1.0261]])
        self.v_cor = np.array([[-0.3855, -0.0916, -0.9773, -0.9707, -1.2140,
                                0.3188, -0.9940, -0.0308, -0.3229, -0.2582]])

        # set the expected results -- echo intensity conversion from counts to dB
        self.dB = np.array([[0.00, 11.25, 22.50, 33.75, 45.00, 56.25, 67.50, 78.75, 90.00, 101.25, 112.50]])

    def test_adcp_beam(self):
        """
        Directly tests DPA functions adcp_beam_eastward, adcp_beam_northward,
        adcp_beam_vertical, and adcp_beam_error. Indirectly tests adcp_beam2ins,
        adcp_ins2earth and magnetic_correction functions. All three functions
        must return the correct output for final tests cases to work.

        Values based on those defined in DPS:

            OOI (2012). Data Product Specification for Velocity Profile and Echo
                Intensity. Document Control Number 1341-00750.
                https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
                >> Controlled >> 1000 System Level >>
                1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by:

            2013-04-10: Christopher Wingard. Initial code.
            2014-02-06: Christopher Wingard. Added tests to confirm arrays of
                        arrays can be processed (in other words, vectorized the
                        code).
            2015-06-23: Russell Desiderio. Revised documentation. Added unit test
                        for the function adcp_beam_error.
            2019-08-13: Christopher Wingard. Adds functionality to compute a 3-beam solution
                        and cleans up syntax used in the function.

        Notes:

            The original suite of tests within this function did not provide a
            test for adcp_beam_error. However, adcp_beam_error and vadcp_beam_error
            are identical functions, and vadcp_beam_error is implicitly tested in the
            test_vadcp_beam function when the 4th output argument of adcp_beam2inst
            is tested. Therefore values to directly test adcp_beam_error were
            then derived from the function itself and included as part of the unit
            test within this code (test_adcp_beam).
        """
        # adjust expected test values to include a 3-beam solution (computed in Matlab using modifications to functions
        # found in the scripts included with this repo).
        u_cor = np.array([[0.1099, -0.3221, -0.4025, np.nan, 0.9243, -0.1595, 0.5387, -0.0651, -0.0726, -0.4735]])
        v_cor = np.array([[-0.3854, -0.0916, -0.9773, np.nan, -1.2140, 0.3188, -0.5926, 0.2514, -0.2932, -0.5255]])
        w = np.array([[0.1401, 0.3977, 0.1870, np.nan, 0.0091, -0.1290, -0.0681, -0.2591, 0.1215, 0.0926]])
        e = np.array([[0.7898,  0.6347, -0.0806, np.nan, 0.0641, 0.07133, 0., 0., 0., 0.]])

        # single record case
        got_u_cor = af.adcp_beam_eastward(self.b1, self.b2, self.b3, self.b4,
                                          self.pg1, self.pg2, self.pg3, self.pg4,
                                          self.heading, self.pitch, self.roll, self.orient,
                                          self.lat, self.lon, self.ntp)
        got_v_cor = af.adcp_beam_northward(self.b1, self.b2, self.b3, self.b4,
                                           self.pg1, self.pg2, self.pg3, self.pg4,
                                           self.heading, self.pitch, self.roll, self.orient,
                                           self.lat, self.lon, self.ntp)
        got_w = af.adcp_beam_vertical(self.b1, self.b2, self.b3, self.b4,
                                      self.pg1, self.pg2, self.pg3, self.pg4,
                                      self.heading, self.pitch, self.roll, self.orient)
        got_e = af.adcp_beam_error(self.b1, self.b2, self.b3, self.b4, self.pg1, self.pg2, self.pg3, self.pg4,)

        # test results
        np.testing.assert_array_almost_equal(got_u_cor, u_cor, 4)
        np.testing.assert_array_almost_equal(got_v_cor, v_cor, 4)
        np.testing.assert_array_almost_equal(got_w, w, 4)
        np.testing.assert_array_almost_equal(got_e, e, 4)

        # reset the test inputs for multiple records
        b1 = np.tile(self.b1, (24, 1))
        b2 = np.tile(self.b2, (24, 1))
        b3 = np.tile(self.b3, (24, 1))
        b4 = np.tile(self.b4, (24, 1))
        pg1 = np.tile(self.pg1, (24, 1))
        pg2 = np.tile(self.pg2, (24, 1))
        pg3 = np.tile(self.pg3, (24, 1))
        pg4 = np.tile(self.pg4, (24, 1))
        heading = np.ones(24, dtype=np.int) * self.heading
        pitch = np.ones(24, dtype=np.int) * self.pitch
        roll = np.ones(24, dtype=np.int) * self.roll
        orient = np.ones(24, dtype=np.int) * self.orient
        lat = np.ones(24) * self.lat
        lon = np.ones(24) * self.lon
        ntp = np.ones(24) * self.ntp

        # reset outputs for multiple records
        u_cor = np.tile(u_cor, (24, 1))
        v_cor = np.tile(v_cor, (24, 1))
        w = np.tile(w, (24, 1))
        e = np.tile(e, (24, 1))

        # multiple record case
        got_u_cor = af.adcp_beam_eastward(b1, b2, b3, b4, pg1, pg2, pg3, pg4,
                                          heading, pitch, roll, orient,
                                          lat, lon, ntp)
        got_v_cor = af.adcp_beam_northward(b1, b2, b3, b4, pg1, pg2, pg3, pg4,
                                           heading, pitch, roll, orient,
                                           lat, lon, ntp)
        got_w = af.adcp_beam_vertical(b1, b2, b3, b4, pg1, pg2, pg3, pg4,
                                      heading, pitch, roll, orient)
        got_e = af.adcp_beam_error(b1, b2, b3, b4, pg1, pg2, pg3, pg4,)

        # test results
        np.testing.assert_array_almost_equal(got_u_cor, u_cor, 4)
        np.testing.assert_array_almost_equal(got_v_cor, v_cor, 4)
        np.testing.assert_array_almost_equal(got_w, w, 4)
        np.testing.assert_array_almost_equal(got_e, e, 4)


    def test_vadcp_b_beam(self):
        """
        Directly tests DPA functions vadcp_b_beam_eastward, vadcp_b_beam_northward,
        vadcp_b_beam_vertical_est, and vadcp_b_beam_vertical_true. Indirectly tests vadcp_b_beam2ins,
        vadcp_b_ins2earth and magnetic_correction functions. All three functions
        must return the correct output for final tests cases to work.

        Values based on those defined in DPS:

            OOI (2012). Data Product Specification for Velocity Profile and Echo
                Intensity. Document Control Number 1341-00750.
                https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
                >> Controlled >> 1000 System Level >>
                1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by:

            2013-04-10: Christopher Wingard. Initial code.
            2014-02-06: Christopher Wingard. Added tests to confirm arrays of
                        arrays can be processed (in other words, vectorized the
                        code).
            2015-06-23: Russell Desiderio. Revised documentation. Added unit test
                        for the function adcp_beam_error.
            2019-08-13: Christopher Wingard. Adds functionality to compute a 3-beam solution
                        and cleans up syntax used in the function.
            2025-02-19: Samuel Dahlberg. Duplicated original function (test_vadcp_beam) and
                        modified to work with the vadcp_b instrument.

        """

        # expected test result values computed from Nrotek provided Matlab code, then run through magnetic declination
        # and correction code in adcp_functions.
        u_cor = np.array([[-0.0254, -0.0413, -0.0047,  0.0872, -0.0342,  0.1109, -0.0482, 0.0283, -0.0804, -0.0345]])
        v_cor = np.array([[-0.0046,  0.0222,  0.0541,  0.0815, -0.0145,  0.0372,  0.1068, 0.0006,  0.0074,  0.1268]])
        w = np.array([[-0.0336, -0.0032, 0.0097, 0.0326, -0.0047, 0.0017, -0.0125, -0.0134, 0.0049, -0.0430]])
        w_true = np.array([[0.0531, -0.0316, -0.0504, 0.0444, 0.0160, 0.0058, -0.0687, 0.6178, 0.0285, 0.0094]])

        # Beam values, percent good values, transformation matrix, heading/pitch/roll, lat/lon, and date pulled from
        # instrument telemetered files
        b1 = np.array([[-0.0450, 0.0160, 0.0390, 0.0590, -0.0230, 0.0000, 0.0520, -0.0270, 0.0080, -0.0020]])
        b2 = np.array([[-0.0190, -0.0090, 0.0100, 0.0480, 0.0060, 0.0250, -0.0130, 0.0030, 0.0050, -0.0030]])
        b3 = np.array([[-0.0570, -0.0220, 0.0040, 0.0630, -0.0350, 0.0450, -0.0400, -0.0100, -0.0450, -0.0970]])
        b4 = np.array([[-0.0010, 0.0020, -0.0200, -0.0530, 0.0350, -0.0630, -0.0490, -0.0140, 0.0480, -0.0590]])
        b5 = np.array([[0.0530, -0.0320, -0.0510, 0.0440, 0.0160, 0.0060, -0.0700, 0.6180, 0.0280, 0.0080]])

        pg1 = np.array([[96.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 96.0]])
        pg2 = np.array([[100.0, 80, 100.0, 100.0, 100.0, 100.0, 99.0, 100.0, 95.0, 94.0]])
        pg3 = np.array([[97.0, 80, 88.0, 94.0, 91.0, 86.0, 83.0, 79.0, 79.0, 74.0]])
        pg4 = np.array([[85.0, 80, 74.0, 69.0, 57.0, 58.9, 61.0, 84.0, 95.0, 96.0]])
        pg5 = np.array([[63.0, 60.0, 62.0, 80.0, 87.0, 97.0, 98.0, 58.9, 84.0, 94.0]])

        tm = [1.1831, 0.0000, -1.1831, 0.0000, 0.0000, -1.1831, 0.0000, 1.1831, 0.5518, 0.0000, 0.5518, 0.0000,
              0.0000, 0.5518, 0.0000, 0.5518]

        heading = 300.1
        pitch = 0.6
        roll = -0.2

        lat = 44.5
        lon = -125.4

        ntp = 3948947793.3

        # single record case
        got_u_cor = af.vadcp_b_beam_eastward(b1, b2, b3, b4,
                                          pg1, pg2, pg3, pg4,
                                          heading, pitch, roll,
                                          lat, lon, ntp, tm)
        got_v_cor = af.vadcp_b_beam_northward(b1, b2, b3, b4,
                                           pg1, pg2, pg3, pg4,
                                           heading, pitch, roll,
                                           lat, lon, ntp, tm)
        got_w_est = af.vadcp_b_beam_vertical_est(b1, b2, b3, b4,
                                      pg1, pg2, pg3, pg4,
                                      heading, pitch, roll, tm)
        got_w_true = af.vadcp_b_beam_vertical_true(b1, b2, b3, b4, b5, pg1, pg2, pg3, pg4,
                                             pg5, heading, pitch, roll, tm)

        # test results
        np.testing.assert_array_almost_equal(got_u_cor, u_cor, 4)
        np.testing.assert_array_almost_equal(got_v_cor, v_cor, 4)
        np.testing.assert_array_almost_equal(got_w_est, w, 4)
        np.testing.assert_array_almost_equal(got_w_true, w_true, 4)

        # reset the test inputs for multiple records
        b1 = np.tile(b1, (24, 1))
        b2 = np.tile(b2, (24, 1))
        b3 = np.tile(b3, (24, 1))
        b4 = np.tile(b4, (24, 1))
        pg1 = np.tile(pg1, (24, 1))
        pg2 = np.tile(pg2, (24, 1))
        pg3 = np.tile(pg3, (24, 1))
        pg4 = np.tile(pg4, (24, 1))
        heading = np.ones(24, dtype=np.int) * heading
        pitch = np.ones(24, dtype=np.int) * pitch
        roll = np.ones(24, dtype=np.int) * roll
        lat = np.ones(24) * lat
        lon = np.ones(24) * lon
        ntp = np.ones(24) * ntp

        # reset outputs for multiple records
        u_cor = np.tile(u_cor, (24, 1))
        v_cor = np.tile(v_cor, (24, 1))
        w = np.tile(w, (24, 1))
        w_true = np.tile(w_true, (24, 1))

        # multiple record case
        got_u_cor = af.vadcp_b_beam_eastward(b1, b2, b3, b4,
                                             pg1, pg2, pg3, pg4,
                                          heading, pitch, roll,
                                          lat, lon, ntp, tm)
        got_v_cor = af.vadcp_b_beam_northward(b1, b2, b3, b4, pg1, pg2, pg3, pg4,
                                           heading, pitch, roll,
                                           lat, lon, ntp, tm)
        got_w_est = af.vadcp_b_beam_vertical_est(b1, b2, b3, b4, pg1, pg2, pg3, pg4,
                                      heading, pitch, roll, tm)
        got_w_true = af.vadcp_b_beam_vertical_true(b1, b2, b3, b4, b5, pg1, pg2, pg3, pg4,
                                             pg5, heading, pitch, roll, tm)

        # test results
        np.testing.assert_array_almost_equal(got_u_cor, u_cor, 4)
        np.testing.assert_array_almost_equal(got_v_cor, v_cor, 4)
        np.testing.assert_array_almost_equal(got_w_est, w, 4)
        np.testing.assert_array_almost_equal(got_w_true, w_true, 4)

    def test_adcp_earth(self):
        """
        Tests magnetic_correction function for ADCPs set to output data in the
        Earth Coordinate system. Also, tests simple scaling functions for the
        vertical and error velocity profiles.

        Values were not defined in DPS, were recreated using test values above.

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by:

            2014-02-06: Christopher Wingard. Initial code.
            2015-06-10: Russell Desiderio.
                            Changed adcp_ins2earth to require the units of the compass
                            data to be in centidegrees.
            2019-08-13: Christopher Wingard. Added functionality to test scaling functions
                        and incorporated instrument fill value (system fill values do
                        not occur).
        """
        # test the magnetic variation correction and simple scaling functions
        got_u_cor = af.adcp_earth_eastward(self.u, self.v, self.depth, self.lat, self.lon, self.ntp)
        got_v_cor = af.adcp_earth_northward(self.u, self.v, self.depth, self.lat, self.lon, self.ntp)
        got_w = af.adcp_earth_vertical(self.w)
        got_e = af.adcp_earth_error(self.e)

        np.testing.assert_array_almost_equal(got_u_cor, self.u_cor, 4)
        np.testing.assert_array_almost_equal(got_v_cor, self.v_cor, 4)
        np.testing.assert_array_almost_equal(got_w, self.w / 1000., 4)
        np.testing.assert_array_almost_equal(got_e, self.e / 1000., 4)

        # reset the test inputs for multiple records
        u = np.tile(self.u, (24, 1)).astype(int)
        u[:, (2, 5, 7)] = ADCP_FILLVALUE
        v = np.tile(self.v, (24, 1)).astype(int)
        v[:, (2, 5, 7)] = ADCP_FILLVALUE
        w = np.tile(self.w, (24, 1))
        e = np.tile(self.e, (24, 1))
        depth = np.ones(24) * self.depth
        lat = np.ones(24) * self.lat
        lon = np.ones(24) * self.lon
        ntp = np.ones(24) * self.ntp

        # reset expected results for multiple records
        u_cor = np.tile(self.u_cor, (24, 1))
        u_cor[:, (2, 5, 7)] = np.nan
        v_cor = np.tile(self.v_cor, (24, 1))
        v_cor[:, (2, 5, 7)] = np.nan

        # compute the results for multiple records
        got_u_cor = af.adcp_earth_eastward(u, v, depth, lat, lon, ntp)
        got_v_cor = af.adcp_earth_northward(u, v, depth, lat, lon, ntp)
        got_w = af.adcp_earth_vertical(w)
        got_e = af.adcp_earth_error(e)

        # test the magnetic variation correction (relaxing precision to account for integer/float round off errors)
        np.testing.assert_array_almost_equal(got_u_cor, u_cor, 3)
        np.testing.assert_array_almost_equal(got_v_cor, v_cor, 3)
        np.testing.assert_array_almost_equal(got_w, w / 1000., 4)
        np.testing.assert_array_almost_equal(got_e, e / 1000., 4)

    def test_adcp_backscatter(self):
        """
        Tests echo intensity scaling function (adcp_backscatter) for ADCPs
        in order to convert from echo intensity in counts to dB.

        Values were not defined in DPS, were created in Matlab using test values above.

        OOI (2012). Data Product Specification for Velocity Profile and Echo
            Intensity. Document Control Number 1341-00750.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00750_Data_Product_SPEC_VELPROF_OOI.pdf)

        Implemented by Christopher Wingard, 2014-02-06
                       Russell Desiderio, 2015-06-25. Added tests for fill values.
                       Christopher Wingard, 2019-08-13. Removed fill values tests. Not a valid scenario.
        """
        # the single record case
        got = af.adcp_backscatter(self.echo, self.sfactor)
        np.testing.assert_array_almost_equal(got, self.dB, 4)

        # the multi-record case -- inputs
        raw = np.tile(self.echo, (24, 1))
        sf = np.ones(24) * self.sfactor

        # the multi-record case -- outputs
        dB = np.tile(self.dB, (24, 1))
        got = af.adcp_backscatter(raw, sf)
        np.testing.assert_array_almost_equal(got, dB, 4)

    def test_vadcp_beam_vertical_true(self):
        """
        Tests vadcp_beam_vertical_true function for the specialized 5-beam ADCP. The other velocity
        estimates for the VADCP use the functions tested above in test_adcp_beam.

        Test values come from the function test_vadcp_beam, in this module.

        Implemented by:

            2014-06-22: Russell Desiderio. Initial code.
            2019-08-14: Christopher Wingard. Reset to only test the vadcp_beam_vertical_true function.
        """
        # inputs
        b1 = np.ones((10, 10), dtype=np.int) * -325
        b2 = np.ones((10, 10), dtype=np.int) * 188
        b3 = np.ones((10, 10), dtype=np.int) * 168
        b4 = np.ones((10, 10), dtype=np.int) * -338
        b5 = np.ones((10, 10), dtype=np.int) * -70
        pg1 = np.ones((10, 10), dtype=np.int) * 100
        pg2 = np.ones((10, 10), dtype=np.int) * 100
        pg3 = np.ones((10, 10), dtype=np.int) * 100
        pg4 = np.ones((10, 10), dtype=np.int) * 100
        pg5 = np.ones((10, 10), dtype=np.int) * 100

        # create some bad data points, which should yield NaNs for output
        pg1[:, (2, 3)] = 24
        pg2[:, (2, 3)] = 24
        pg5[:, 9] = 24

        # units of centidegrees
        heading = np.array([30, 30, 30, 30, 30,
                            32, 32, 32, 32, 32]) * 100
        pitch = np.array([0, 2, 3, 3, 1, 2, 2, 3, 3, 1]) * 100
        roll = np.array([0, 4, 3, 4, 3, 3, 4, 3, 4, 3]) * 100
        orient = np.ones(10, dtype=np.int)

        # expected outputs
        vlu_5bm_xpctd = np.array([[0.07000000, -0.00824854, -0.00804866, -0.02112871, 0.01775751,
                                   0.00485518, -0.00824854, -0.00804866, -0.02112871, 0.01775751]])
        vlu_5bm_xpctd = np.tile(vlu_5bm_xpctd.T, (1, 10))
        vlu_5bm_xpctd[:, (2, 3, 9)] = np.nan

        vlu_5bm_calc = af.vadcp_beam_vertical_true(b1, b2, b3, b4, b5, pg1, pg2, pg3, pg4, pg5,
                                                   heading, pitch, roll, orient)

        np.testing.assert_array_almost_equal(vlu_5bm_calc, vlu_5bm_xpctd, 6)

    def test_adcp_ins2earth_orientation(self):
        """
        Test the adcp worker routine adcp_inst2earth when the vertical orientation
        toggle switch codes for downward-looking (orient = 0).

        The instrument to earth coordinate transformation was coded in matlab using
        the January 2010 version of the Teledyne RD Instruments "ADCP Coordinate
        Transformation" manual. This matlab code was then used to generate the check
         values for the downward looking case.

        Implemented by:

            2014-07-02: Russell Desiderio. Initial code.
        """
        # input values: these are the output of adcp_beam2inst in test_vadcp_beam
        # (velocities are in instrument coordinates)
        u = np.array([[-749.95582864]])
        v = np.array([[-739.72251324]])
        w = np.array([[-81.67564404]])
        heading = np.array([3200])  # units of centidegrees
        pitch = np.array([300])  # units of centidegrees
        roll = np.array([400])  # units of centidegrees

        # check test: upwards looking case
        orient_1 = np.array([1])

        # expected outputs, earth coordinates, upwards case, from test_vadcp_beam which
        # agrees with the test matlab code values to (much) better than single precision.
        vle = np.array([[247.015599]])
        vln = np.array([[-1027.223026]])
        vlu = np.array([[-9.497397]])
        xpctd = np.hstack((vle, vln, vlu))

        # calculated upwards looking case
        uu, vv, ww = af.adcp_ins2earth(u, v, w, heading, pitch, roll, orient_1)
        calc = np.hstack((uu, vv, ww))

        # test results
        np.testing.assert_array_almost_equal(calc, xpctd, 6)

        # primary test: downwards looking case.
        orient_0 = np.array([0])

        # expected outputs, earth coordinates, downwards case, from matlab test code
        vle = np.array([[-1029.9328104]])
        vln = np.array([[-225.7064203]])
        vlu = np.array([[-67.7426771]])
        xpctd = np.hstack((vle, vln, vlu))

        # calculated downwards looking case
        uu, vv, ww = af.adcp_ins2earth(u, v, w, heading, pitch, roll, orient_0)
        calc = np.hstack((uu, vv, ww))

        # test results
        np.testing.assert_array_almost_equal(calc, xpctd, 6)

    def test_adcp_bin_depths_meters(self):
        """
        Test the adcp_bin_depths_meters function.

        Implemented by:
            Craig Risien, January 2015. Initial code.
            Russell Desiderio. 26-Jun-2015. Added time-vectorized unit test after modifying DPA.
                               30-Jun-2015. Added fill value unit test.
        """
        sfill = SYSTEM_FILLVALUE

        ### scalar time case (1) - adcp looking up
        # test inputs - note, CI will be sending these into the DPAs as ndarrays, not python scalars.
        adcp_orientation = 1
        bin_size = 400
        dist_first_bin = 900
        num_bins = 29
        sensor_depth = 450
        # expected outputs
        # note that the output should be a row vector, not a 1D array.
        xpctd_bins_up = np.array([[441., 437., 433., 429., 425., 421., 417., 413., 409., 405., 401., 397., 393., 389.,
                                  385., 381., 377., 373., 369., 365., 361., 357., 353., 349., 345., 341., 337., 333.,
                                  329.]])
         # calculate bin depths
        calc_bins_up = af.adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)

        # compare calculated results to expected results
        np.testing.assert_allclose(calc_bins_up, xpctd_bins_up, rtol=0.000001, atol=0.000001)

        ### scalar time case (2) - adcp looking down
        # test inputs
        adcp_orientation = np.array([0])
        bin_size = np.array([400])
        dist_first_bin = np.array([900])
        num_bins = np.array([29])
        sensor_depth = np.array([7])
        # expected outputs
        xpctd_bins_down = np.array([[16., 20., 24., 28., 32., 36., 40., 44., 48., 52., 56., 60., 64., 68., 72., 76., 80.,
                                    84., 88., 92., 96., 100., 104., 108., 112., 116., 120., 124., 128.]])
        # calculate bin depths
        calc_bins_down = af.adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)

        # compare calculated results to expected results
        np.testing.assert_allclose(calc_bins_down, xpctd_bins_down, rtol=0.000001, atol=0.000001)

        ### time-vectorized case; cat the above two scalar cases together.
        # inputs
        dist_first_bin = np.array([900, 900])
        bin_size = np.array([400, 400])
        num_bins = np.array([29, 29])
        sensor_depth = np.array([450, 7])
        adcp_orientation = np.array([1, 0])
        # expected
        xpctd_bins = np.vstack((xpctd_bins_up, xpctd_bins_down))
        # calculated
        calc_bins = af.adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)
        # compare calculated results to expected results
        np.testing.assert_allclose(calc_bins, xpctd_bins, rtol=0.000001, atol=0.000001)

        ### time-vectorized fill cases - test the action on a fill value in each of the 5 input data streams,
        # plus one instance of all good data.
        num_bins = np.array([29, 29, 29, sfill, 29, 29])  # NOTE: DPA uses only first num_bins value
        dist_first_bin = np.array([900, sfill, 900, 900, 900, 900])
        bin_size = np.array([400, 400, sfill, 400, 400, 400])
        sensor_depth = np.array([450, 7, 450, 7, 450, sfill])
        adcp_orientation = np.array([1, 0, 1, 0, sfill, 0])
        # 1st and 4th rows will have non-Nan data.
        xpctd_bins = np.tile(np.nan, (6, 29))
        xpctd_bins[0, :] = xpctd_bins_up
        xpctd_bins[3, :] = xpctd_bins_down
        # calculated
        calc_bins = af.adcp_bin_depths_meters(dist_first_bin, bin_size, num_bins, sensor_depth, adcp_orientation)
        # compare calculated results to expected results
        np.testing.assert_allclose(calc_bins, xpctd_bins, rtol=0.000001, atol=0.000001)

    def test_adcp_bin_depths_dapa(self):
        """
        Test the adcp_bin_depths_dapa function.

        Values based on z_from_p check values.

        Implemented by:
            Craig Risien, January 2015. Initial code.
            Russell Desiderio. 26-Jun-2015. Corrected pressure type and units and z_from_p usage.
                               30-Jun-2015. Added fill value unit test.
        """
        sfill = SYSTEM_FILLVALUE

        ### scalar time case (1) - adcp looking up
        # test inputs - note, CI will be sending these into the DPAs as ndarrays, not python scalars.
        # test inputs
        adcp_orientation = np.array([1])
        bin_size = np.array([400])
        dist_first_bin = np.array([900])
        latitude = np.array([4.0])
        num_bins = np.array([10])
        # input adcp pressure has units of decaPascals
        pressure = np.array([600000])
        # according to the z_from_p check value at 600db, this gives a depth of 595.8262 m
        # expected outputs
        # note that the output should be a row vector, not a 1D array.
        expected_bins_up = 0.8253480 + np.array([[586, 582, 578, 574, 570, 566, 562, 558, 554, 550]])
        # calculate bin depths
        calculated_bins = af.adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins_up, rtol=0.0, atol=0.000001)

        ### scalar time case (2) - adcp looking down
        # test inputs - should also work with python scalars, but this is not necessary
        adcp_orientation = 0
        bin_size = 400
        dist_first_bin = 900
        latitude = 4.0
        num_bins = 10
        pressure = 10000
        # expected depth from a pressure of 10000 decapascals is 9.94460074 m
        # expected outputs
        expected_bins_down = 0.9445834 + np.array([[18, 22, 26, 30, 34, 38, 42, 46, 50, 54]])
        # calculate bin depths
        calculated_bins = af.adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins_down, rtol=0.0, atol=0.000001)

        ### time-vectorized test - combine the 2 above
        adcp_orientation = np.array([1, 0])
        bin_size = np.array([400, 400])
        dist_first_bin = np.array([900, 900])
        latitude = np.array([4.0, 4.0])
        num_bins = np.array([10, 10])
        pressure = np.array([600000, 10000])
        #
        expected_bins = np.vstack((expected_bins_up, expected_bins_down))
        #
        calculated_bins = af.adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        #
        np.testing.assert_allclose(calculated_bins, expected_bins, rtol=0.0, atol=0.001)

        ### time-vectorized fill cases - test the action on a fill value in each of the 5 input data streams,
        # plus one instance of all good data.
        num_bins = np.array([10, 10, 10, sfill, 10, 10])  # NOTE: DPA uses only first num_bins value
        dist_first_bin = np.array([900, sfill, 900, 900, 900, 900])
        bin_size = np.array([400, 400, sfill, 400, 400, 400])
        pressure = np.array([600000, 10000, 600000, 10000, 600000, sfill])
        adcp_orientation = np.array([1, 0, 1, 0, sfill, 0])
        latitude = np.array([4.0, 4.0, 4.0, 4.0, 4.0, 4.0])
        # 1st and 4th rows will have non-Nan data.
        expected_bins = np.tile(np.nan, (6, 10))
        expected_bins[0, :] = expected_bins_up
        expected_bins[3, :] = expected_bins_down
        # calculated
        calculated_bins = af.adcp_bin_depths_dapa(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins, rtol=0.0, atol=0.001)

    def test_adcp_bin_depths_bar(self):
        """
        Test the adcp_bin_depths_bar function.

        Values based on z_from_p check values.

        Implemented by:
            Craig Risien, January 2015. Initial code.
            Russell Desiderio. 26-Jun-2015. Corrected pressure type and units and z_from_p usage.
                               30-Jun-2015. Added fill value unit test.
        """
        sfill = SYSTEM_FILLVALUE

        ### scalar time case (1) - adcp looking up
        # test inputs - note, CI will be sending these into the DPAs as ndarrays, not python scalars.
        # test inputs
        adcp_orientation = np.array([1])
        bin_size = np.array([400])
        dist_first_bin = np.array([900])
        latitude = np.array([4.0])
        num_bins = np.array([10])
        # water pressure of gliders has units of bar
        pressure = np.array([60])
        # according to the z_from_p check value at 600db, this gives a depth of 595.8262 m
        # expected outputs
        # note that the output should be a row vector, not a 1D array.
        expected_bins_up = 0.8253480 + np.array([[586, 582, 578, 574, 570, 566, 562, 558, 554, 550]])
        # calculate bin depths
        calculated_bins = af.adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins_up, rtol=0.0, atol=0.000001)

        ### scalar time case (2) - adcp looking down
        # test inputs - should also work with python scalars, but this is not necessary
        adcp_orientation = 0
        bin_size = 400
        dist_first_bin = 900
        latitude = 4.0
        num_bins = 10
        pressure = 1
        # expected depth from a pressure of 1 bar is 9.94460074 m
        # expected outputs
        expected_bins_down = 0.9445834 + np.array([[18, 22, 26, 30, 34, 38, 42, 46, 50, 54]])
        # calculate bin depths
        calculated_bins = af.adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins_down, rtol=0.0, atol=0.000001)

        ### time-vectorized test - combine the 2 above
        adcp_orientation = np.array([1, 0])
        bin_size = np.array([400, 400])
        dist_first_bin = np.array([900, 900])
        latitude = np.array([4.0, 4.0])
        num_bins = np.array([10, 10])
        pressure = np.array([60, 1])
        #
        expected_bins = np.vstack((expected_bins_up, expected_bins_down))
        #
        calculated_bins = af.adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        #
        np.testing.assert_allclose(calculated_bins, expected_bins, rtol=0.0, atol=0.001)

        ### time-vectorized fill cases - test the action on a fill value in each of the 5 input data streams,
        # plus one instance of all good data.
        num_bins = np.array([10, 10, 10, sfill, 10, 10])  # NOTE: DPA uses only first num_bins value
        dist_first_bin = np.array([900, sfill, 900, 900, 900, 900])
        bin_size = np.array([400, 400, sfill, 400, 400, 400])
        pressure = np.array([60, 1, 60, 1, 60, sfill])
        adcp_orientation = np.array([1, 0, 1, 0, sfill, 0])
        latitude = np.array([4.0, 4.0, 4.0, 4.0, 4.0, 4.0])
        # 1st and 4th rows will have non-Nan data.
        expected_bins = np.tile(np.nan, (6, 10))
        expected_bins[0, :] = expected_bins_up
        expected_bins[3, :] = expected_bins_down
        # calculated
        calculated_bins = af.adcp_bin_depths_bar(dist_first_bin, bin_size, num_bins, pressure, adcp_orientation, latitude)
        # compare calculated results to expected results
        np.testing.assert_allclose(calculated_bins, expected_bins, rtol=0.0, atol=0.001)

    def test_z_from_p(self):
        """
        Test the z_from_p function, which calculates depth in meters from pressure in decibars
        assuming a TEOS-10 standard ocean salinity and a TEOS-10 conservative temperature of 0 deg_C.

        Check values are taken from the TEOS-10 version 3.05 matlab documentation (see the function
        z_from_p in the adcp_functions.py module).

        There are no time-vectorized unit tests in this test function; the DPA that calls this
        function does have time-vectorized unit tests, however.

        Implemented by:

            2014-06-26: Russell Desiderio. Initial code.
            2015-07-01: Russell Desiderio. Updated check values to TEOS-10 ver. 3.05.
        """
        # test inputs
        p = np.array([10.0, 50.0, 125.0, 250.0, 600.0, 1000.0])
        lat = np.ones(6) * 4.0
        # outputs
        xpctd = np.array([-9.9445834469453,  -49.7180897012550, -124.2726219409978,
                          -248.4700576548589, -595.8253480356214, -992.0919060719987])
        calc = af.z_from_p(p, lat)
        # test relative accuracy
        np.testing.assert_allclose(calc, xpctd, rtol=0.00000001, atol=0.0)
        # test absolute accuracy
        np.testing.assert_allclose(calc, xpctd, rtol=0.0, atol=0.00000001)
