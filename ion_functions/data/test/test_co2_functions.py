#!/usr/bin/env python

"""
@package ion_functions.test.co2_functions
@file ion_functions/test/co2_functions.py
@author Christopher Wingard
@brief Unit tests for co2_functions module
"""

from nose.plugins.attrib import attr
from ion_functions.test.base_test import BaseUnitTestCase

import numpy as np
from ion_functions.data import co2_functions as co2func
from ion_functions.utils import fill_value


@attr('UNIT', group='func')
class Testpco2FunctionsUnit(BaseUnitTestCase):

    def setUp(self):
        ###### Test data for PCO2W ######
        raw_strings = np.array([
            '*BC2705D5A7C0E10082005A0CA9090E07CB08E82DCA4B1C0082005A0CA9090E07CD08EC0C3208C38A',
            '*BC2704D5A7E2B1007E005A0CB1022F07C40443099F226D007F005A0CAF022F07C404400C3F08BE2E',
            '*BC2704D5A7FEC90080005A0CAD028707CC03160B711800008300580CAC028607CC03160C4007389C',
            '*BC2704D5A8006F0083005B0CA1028D07DE02F70BA016AF0081005A0CA2028D07E202F70C4007A16D',
            '*BC2704D5A802150080005A0C98028D07E902DE0BAE15BF0081005A0C98028E07EB02DF0C40080C6B',
            '*BC2704D5A803BB007F00590C98028307EB02EA0B6B161B008100580C94028107EE02EB0C41085372',
            '*BC2704D5A80560007F005A0C9A027407E403050B26171F0083005C0C99027607E903060C4008862B',
            '*BC2704D5A80707007F005A0CA0027007DE03210AF9182A008000590CA0026D07DE03240C400895E4',
            '*BC2704D5A808AD0082005B0CA2026607D203470AC019A00080005C0CA6026607D403490C3F0899FF',
            '*BC2704D5A80A54007F00560CAB025407CC03860A701BCF0083005C0CAA025707D003840C3F089BE2',
            '*BC2704D5A80BF90080005A0CB3024907C603B60A2D1D99008100580CAF024707C603B90C3F089C58',
            '*BC2704D5A80E140080005B0CB3023807C0041009CA20C40082005A0CB3023807C004110C3F08A0D4',
            '*BC2704D5A8102F007F00580CB9022D07BA04350999222D0080005A0CB4022E07BC04370C3E08B067',
            '*BC2704D5A812E70081005B0CBE022107B00479094A24AD0080005A0CBC022007B8047B0C3E08CEE4',
        ])

        # reagent constants (instrument and reagent bag specific)
        self.calt = np.ones(14) * 4.6539
        self.cala = np.ones(14) * 0.0422
        self.calb = np.ones(14) * 0.6761
        self.calc = np.ones(14) * -1.5798

        # expected outputs
        self.therm = np.array([7.3151, 7.4258, 16.2306, 13.8108, 11.3900,
                               9.8019, 8.6674, 8.3345, 8.2458, 8.2014,
                               8.1792, 8.0905, 7.7359, 7.0716])
        self.pco2 = np.array([fill_value, 609.8626, 394.3221, 351.6737, 321.4986,
                              324.0670, 339.4317, 358.3335, 388.1735, 436.8431,
                              481.1713, 566.8172, 607.5355, 685.8555])

        # parse the data strings
        self.mtype = np.zeros(14)
        self.light = np.zeros((14, 14))
        self.traw = np.zeros(14)
        for i in range(14):
            # parse the raw strings into subelements, such as the driver would
            # provide.
            s = raw_strings[i]
            self.mtype[i] = int((s[5:7]), 16)
            self.traw[i] = int((s[75:79]), 16)
            strt = 15
            step = 4
            for j in range(14):
                self.light[i, j] = int((s[strt:strt+step]), 16)
                strt += step

            if self.mtype[i] == 5:
                a434blnk = self.light[i, 6]
                a620blnk = self.light[i, 7]

        self.a434blnk = np.ones(14) * a434blnk
        self.a620blnk = np.ones(14) * a620blnk

    def test_pco2_pco2wat(self):
        """
        Test pco2_pco2wat function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Seawater. Document Control Number 1341-00490.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00490_Data_Product_SPEC_PCO2WAT_OOI.pdf)

        Implemented by Christopher Wingard, April 2013

        Updated 2017-04-04 to pass raw thermistor and blank values into pco2_pco2wat
        as indicated in the function description and DPS.
        """

        # compute the thermistor temperature in deg_C, derive blanks and then
        # calculate pco2.

        ### bulk case ###
        tout = co2func.pco2_thermistor(self.traw)
        pco2out = co2func.pco2_pco2wat(self.mtype, self.light, self.traw,
                                       fill_value, fill_value, fill_value, fill_value,
                                       self.calt, self.cala, self.calb, self.calc,
                                       self.a434blnk, self.a620blnk)

        np.testing.assert_allclose(pco2out, self.pco2, rtol=1e-4, atol=1e-4)
        np.testing.assert_allclose(tout, self.therm, rtol=1e-4, atol=1e-4)

        ### single record case ###
        indx = 0
        for mtype in self.mtype:
            tout = co2func.pco2_thermistor(self.traw[indx])
            pco2out = co2func.pco2_pco2wat(mtype, self.light[indx, :], self.traw[indx],
                                           fill_value, fill_value, fill_value, fill_value,
                                           self.calt[indx], self.cala[indx],
                                           self.calb[indx], self.calc[indx],
                                           self.a434blnk[indx], self.a620blnk[indx])

            np.testing.assert_allclose(pco2out, self.pco2[indx], rtol=1e-4, atol=1e-4)
            np.testing.assert_allclose(tout, self.therm[indx], rtol=1e-4, atol=1e-4)

            indx += 1

    def test_pco2_ppressure(self):
        """
        Test pco2_ppressure function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Partial Pressure of CO2 in
            Air and Surface Seawater. Document Control Number 1341-00260.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00260_Data_Product_SPEC_PCO2ATM_PCO2SSW_OOI.pdf)

        Implemented by Christopher Wingard, October 2014
        """
        test_data = np.array([
            [674, 1000, 665.19],
            [619, 1000, 610.91],
            [822, 1000, 811.25],
            [973, 1000, 960.28],
            [941, 1000, 928.69],
            [863, 1000, 851.71],
            [854, 1000, 842.83],
            [833, 1000, 822.11],
            [826, 1000, 815.20],
            [814, 1000, 803.36],
            [797, 1000, 786.58],
            [782, 1000, 771.77],
            [768, 1000, 757.96],
            [754, 1000, 744.14],
            [740, 1000, 730.32]
        ])
        ppres = co2func.pco2_ppressure(test_data[:, 0], test_data[:, 1])
        np.testing.assert_allclose(ppres, test_data[:, 2], rtol=1e-2, atol=1e-2)

    def test_pco2_co2flux(self):
        """
        Test pco2_co2flux function.

        Values based on those described in DPS as available on Alfresco:

        OOI (2012). Data Product Specification for Flux of CO2 into the
            Atmosphere. Document Control Number 1341-00270.
            https://alfresco.oceanobservatories.org/ (See: Company Home >> OOI
            >> Controlled >> 1000 System Level >>
            1341-00270_Data_Product_SPEC_CO2FLUX_OOI.pdf)

        Implemented by Christopher Wingard, April 2013
        """
        test_data = np.array([
            [360, 390, 5, 0, 34, -2.063e-08],
            [360, 390, 5, 0, 35, -2.052e-08],
            [360, 390, 5, 10, 34, -1.942e-08],
            [360, 390, 5, 10, 35, -1.932e-08],
            [360, 390, 5, 20, 34, -1.869e-08],
            [360, 390, 5, 20, 35, -1.860e-08],
            [360, 390, 10, 0, 34, -8.250e-08],
            [360, 390, 10, 0, 35, -8.207e-08],
            [360, 390, 10, 10, 34, -7.767e-08],
            [360, 390, 10, 10, 35, -7.728e-08],
            [360, 390, 10, 20, 34, -7.475e-08],
            [360, 390, 10, 20, 35, -7.440e-08],
            [360, 390, 20, 0, 34, -3.300e-07],
            [360, 390, 20, 0, 35, -3.283e-07],
            [360, 390, 20, 10, 34, -3.107e-07],
            [360, 390, 20, 10, 35, -3.091e-07],
            [360, 390, 20, 20, 34, -2.990e-07],
            [360, 390, 20, 20, 35, -2.976e-07],
            [400, 390, 5, 0, 34, 6.875e-09],
            [400, 390, 5, 0, 35, 6.839e-09],
            [400, 390, 5, 10, 34, 6.472e-09],
            [400, 390, 5, 10, 35, 6.440e-09],
            [400, 390, 5, 20, 34, 6.229e-09],
            [400, 390, 5, 20, 35, 6.200e-09],
            [400, 390, 10, 0, 34, 2.750e-08],
            [400, 390, 10, 0, 35, 2.736e-08],
            [400, 390, 10, 10, 34, 2.589e-08],
            [400, 390, 10, 10, 35, 2.576e-08],
            [400, 390, 10, 20, 34, 2.492e-08],
            [400, 390, 10, 20, 35, 2.480e-08],
            [400, 390, 20, 0, 34, 1.100e-07],
            [400, 390, 20, 0, 35, 1.094e-07],
            [400, 390, 20, 10, 34, 1.036e-07],
            [400, 390, 20, 10, 35, 1.030e-07],
            [400, 390, 20, 20, 34, 9.966e-08],
            [400, 390, 20, 20, 35, 9.920e-08],
            [440, 390, 5, 0, 34, 3.438e-08],
            [440, 390, 5, 0, 35, 3.420e-08],
            [440, 390, 5, 10, 34, 3.236e-08],
            [440, 390, 5, 10, 35, 3.220e-08],
            [440, 390, 5, 20, 34, 3.114e-08],
            [440, 390, 5, 20, 35, 3.100e-08],
            [440, 390, 10, 0, 34, 1.375e-07],
            [440, 390, 10, 0, 35, 1.368e-07],
            [440, 390, 10, 10, 34, 1.294e-07],
            [440, 390, 10, 10, 35, 1.288e-07],
            [440, 390, 10, 20, 34, 1.246e-07],
            [440, 390, 10, 20, 35, 1.240e-07],
            [440, 390, 20, 0, 34, 5.500e-07],
            [440, 390, 20, 0, 35, 5.471e-07],
            [440, 390, 20, 10, 34, 5.178e-07],
            [440, 390, 20, 10, 35, 5.152e-07],
            [440, 390, 20, 20, 34, 4.983e-07],
            [440, 390, 20, 20, 35, 4.960e-07]
        ])

        # setup inputs and outputs
        pco2w = test_data[:, 0]
        pco2a = test_data[:, 1]
        u10 = test_data[:, 2]
        t = test_data[:, 3]
        s = test_data[:, 4]
        flux = test_data[:, 5]

        # compute the flux given the inputs
        out = co2func.pco2_co2flux(pco2w, pco2a, u10, t, s)

        # and compare the results
        self.assertTrue(np.allclose(out, flux, rtol=1e-9, atol=1e-9))
