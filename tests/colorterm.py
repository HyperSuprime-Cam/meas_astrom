#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Tests for Colorterm

Run with:
   python colorterm.py
or
   python
   >>> import colorterm; colorterm.run()
"""

import os, re
import unittest
import pickle

import lsst.utils.tests as utilsTests
from lsst.meas.photocal.colorterms import ColortermConfig, ColortermGroupConfig, ColortermLibraryConfig

# From the last page of http://www.naoj.org/staff/nakata/suprime/illustration/colorterm_report_ver3.pdf
# Transformation for griz band between SDSS and SC (estimated with GS83 SEDs)
mitll = ColortermLibraryConfig.fromValues(
    {"*": ColortermGroupConfig.fromValues({"g": ColortermConfig.fromValues("g", "r", -0.00569, -0.0427),
                                           "r": ColortermConfig.fromValues("r", "g", 0.00261, 0.0304),
                                           "i": ColortermConfig.fromValues("i", "r", 0.00586, 0.0827, -0.0118),
                                           "z": ColortermConfig.fromValues("z", "i", 0.000329, 0.0608, 0.0219),
                                       })})
hamamatsu = ColortermLibraryConfig.fromValues(
    {"*": ColortermGroupConfig.fromValues({"g": ColortermConfig.fromValues("g", "r", -0.00928, -0.0824),
                                           "r": ColortermConfig.fromValues("r", "i", -0.00282, -0.0498,
                                                                           -0.0149),
                                           "i": ColortermConfig.fromValues("i", "z", 0.00186, -0.140, -0.0196),
                                           "z": ColortermConfig.fromValues("z", "i", -4.03e-4, 0.0967, 0.0210),
                                       })})


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class ColortermTestCase(unittest.TestCase):
    """A test case for MaskedImage"""
    def setUp(self):
        self.sources = dict(g=0.0, r=0.0, true_g=-0.00928), dict(g=0.0, r=-1.0, true_g=-0.09168)
        self.colorterms = hamamatsu

    def tearDown(self):
        import lsst.meas.astrom.astrometry_net as an
        an.finalize()

    def testTransformSource(self):
        """Check if we can use colour terms"""

        ct = self.colorterms.selectColorTerm("g", version="foo")

        for s in self.sources:
            self.assertEqual(ct.transformSource(s), s["true_g"])

    def testTransformMags(self):
        """Check if we can use colour terms via transformMags"""

        filterName = "g"
        ct = self.colorterms.selectColorTerm("g", version="bar")

        for s in self.sources:
            self.assertEqual(ct.transformMags(s[ct.primary], s[ct.secondary]), s["true_g"])

    def testPickle(self):
        """Ensure color terms can be pickled"""
        colorterms = pickle.loads(pickle.dumps(self.colorterms))
        self.assertEqual(colorterms, self.colorterms)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(ColortermTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
