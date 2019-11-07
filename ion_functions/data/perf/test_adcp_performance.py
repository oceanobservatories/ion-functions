"""
@package ion_functions.data.perf.test_adcp_performance
@file ion_functions/data/perf/test_adcp_performance.py
@author Christopher Wingard
@brief Performance tests for adcp_functions module
"""

import numpy as np
from nose.plugins.attrib import attr

from ion_functions.data.perf.test_performance import PerformanceTestCase
from ion_functions.data import adcp_functions as af

# Note, the VADCP related data products use the same internal functions as the
# family of beam wrapper functions (e.g. adcp_beam_eastward). Thus, those
# functions won't be added to this test. The only way to really speed this
# process up any further is to set the wrapper functions to return all the data
# products for an instrument at once rather than singly. That way the functions
# only need to be run once rather than 4 times for each instrument (a factor of
# four improvement in performance).


@attr('PERF', group='func')
class TestADCPPerformance(PerformanceTestCase):

    def setUp(self):
        # set test inputs -- values from DPS
        # set test inputs -- values from DPS (reset to integers, matches format from instrument)
        self.b1 = np.array([[-30, -295, -514, -234, -188, 203, -325, 305, -204, -294]])
        self.b2 = np.array([[180, -132, 213, 309, 291, 49, 188, 373, 2, 172]])
        self.b3 = np.array([[-398, -436, -131, -473, -443, 188, -168, 291, -179, 8]])
        self.b4 = np.array([[-216, -605, -92, -58, 484, -5, 338, 175, -80, -549]])

        self.pg1 = np.array([[100, 100, 100,  24, 100, 100,  24, 100, 100, 100]])
        self.pg2 = np.array([[100, 100, 100, 100, 100, 100, 100,  24, 100, 100]])
        self.pg3 = np.array([[100, 100, 100, 100, 100, 100, 100, 100,  24, 100]])
        self.pg4 = np.array([[100, 100, 100,  24, 100, 100, 100, 100, 100,  24]])

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

        # set expected results -- velocity profiles in earth coordinates
        # (values in DPS)
        self.uu = np.array([0.2175, -0.2814, -0.1002, 0.4831, 1.2380,
                            -0.2455, 0.6218, -0.1807, 0.0992, -0.9063])
        self.vv = np.array([-0.3367, -0.1815, -1.0522, -0.8676, -0.8919,
                            0.2585, -0.8497, -0.0873, -0.3073, -0.5461])
        self.ww = np.array([0.1401,  0.3977,  0.1870,  0.1637,  0.0091,
                            -0.1290,  0.0334, -0.3017, 0.1384, 0.1966])

    def test_adcp_backscatter(self):
        stats = []

        echo = np.tile(self.echo, (10000, 1))
        sfactor = np.repeat(self.sfactor, 10000)
        self.profile(stats, af.adcp_backscatter, echo, sfactor)

    def test_adcp_beam_eastward(self):
        stats = []

        b1 = np.tile(self.b1, (10000, 1))
        b2 = np.tile(self.b2, (10000, 1))
        b3 = np.tile(self.b3, (10000, 1))
        b4 = np.tile(self.b4, (10000, 1))

        pg1 = np.tile(self.pg1, (10000, 1))
        pg2 = np.tile(self.pg2, (10000, 1))
        pg3 = np.tile(self.pg3, (10000, 1))
        pg4 = np.tile(self.pg4, (10000, 1))

        h = np.repeat(self.heading, 10000)
        p = np.repeat(self.pitch, 10000)
        r = np.repeat(self.roll, 10000)
        vf = np.repeat(self.orient, 10000)
        lat = np.repeat(self.lat, 10000)
        lon = np.repeat(self.lon, 10000)
        z = np.repeat(self.depth, 10000)
        dt = np.repeat(self.ntp, 10000)

        self.profile(stats, af.adcp_beam_eastward, b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt)

    def test_adcp_beam_northward(self):
        stats = []

        b1 = np.tile(self.b1, (10000, 1))
        b2 = np.tile(self.b2, (10000, 1))
        b3 = np.tile(self.b3, (10000, 1))
        b4 = np.tile(self.b4, (10000, 1))

        pg1 = np.tile(self.pg1, (10000, 1))
        pg2 = np.tile(self.pg2, (10000, 1))
        pg3 = np.tile(self.pg3, (10000, 1))
        pg4 = np.tile(self.pg4, (10000, 1))

        h = np.repeat(self.heading, 10000)
        p = np.repeat(self.pitch, 10000)
        r = np.repeat(self.roll, 10000)
        vf = np.repeat(self.orient, 10000)
        lat = np.repeat(self.lat, 10000)
        lon = np.repeat(self.lon, 10000)
        z = np.repeat(self.depth, 10000)
        dt = np.repeat(self.ntp, 10000)

        self.profile(stats, af.adcp_beam_northward, b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf, lat, lon, dt)

    def test_adcp_beam_vertical(self):
        stats = []

        b1 = np.tile(self.b1, (10000, 1))
        b2 = np.tile(self.b2, (10000, 1))
        b3 = np.tile(self.b3, (10000, 1))
        b4 = np.tile(self.b4, (10000, 1))

        pg1 = np.tile(self.pg1, (10000, 1))
        pg2 = np.tile(self.pg2, (10000, 1))
        pg3 = np.tile(self.pg3, (10000, 1))
        pg4 = np.tile(self.pg4, (10000, 1))

        h = np.repeat(self.heading, 10000)
        p = np.repeat(self.pitch, 10000)
        r = np.repeat(self.roll, 10000)
        vf = np.repeat(self.orient, 10000)

        self.profile(stats, af.adcp_beam_vertical, b1, b2, b3, b4, pg1, pg2, pg3, pg4, h, p, r, vf)

    def test_adcp_beam_error(self):
        stats = []

        b1 = np.tile(self.b1, (10000, 1))
        b2 = np.tile(self.b2, (10000, 1))
        b3 = np.tile(self.b3, (10000, 1))
        b4 = np.tile(self.b4, (10000, 1))

        pg1 = np.tile(self.pg1, (10000, 1))
        pg2 = np.tile(self.pg2, (10000, 1))
        pg3 = np.tile(self.pg3, (10000, 1))
        pg4 = np.tile(self.pg4, (10000, 1))

        self.profile(stats, af.adcp_beam_error, b1, b2, b3, b4, pg1, pg2, pg3, pg4)

    def test_adcp_earth_eastward(self):
        stats = []

        u = np.tile(self.uu, (10000, 1))
        v = np.tile(self.vv, (10000, 1))

        lat = np.repeat(self.lat, 10000)
        lon = np.repeat(self.lon, 10000)
        z = np.repeat(self.depth, 10000)
        dt = np.repeat(self.ntp, 10000)

        self.profile(stats, af.adcp_earth_eastward, u, v, z, lat, lon, dt)

    def test_adcp_earth_northward(self):
        stats = []

        u = np.tile(self.uu, (10000, 1))
        v = np.tile(self.vv, (10000, 1))

        lat = np.repeat(self.lat, 10000)
        lon = np.repeat(self.lon, 10000)
        z = np.repeat(self.depth, 10000)
        dt = np.repeat(self.ntp, 10000)

        self.profile(stats, af.adcp_earth_northward, u, v, z, lat, lon, dt)

    def test_adcp_earth_vertical(self):
        stats = []

        w = np.tile(self.ww, (10000, 1))

        self.profile(stats, af.adcp_earth_vertical, w)
        # adcp_earth_error is the same transform, so this test applies to both
