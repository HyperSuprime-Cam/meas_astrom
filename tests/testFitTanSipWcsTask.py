#!/usr/bin/env python

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
import math
import unittest

import numpy
try:
    import matplotlib
    matplotlib.use("Agg")
    import pylab
except ImportError:
    pass

import lsst.pipe.base
import lsst.utils.tests as tests
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
from lsst.meas.algorithms import LoadReferenceObjectsTask
from lsst.meas.base import SingleFrameMeasurementTask
from lsst.meas.astrom import FitTanSipWcsTask
from lsst.meas.astrom.sip import makeCreateWcsWithSip

class BaseTestCase(unittest.TestCase):
    """A test case for CreateWcsWithSip

    Use involves setting one class attribute:
    * MatchClass: match class, e.g. ReferenceMatch or SourceMatch

    This test is a bit messy because it exercises two templatings of makeCreateWcsWithSip,
    the underlying TAN-SIP WCS fitter, but only one of those is supported by FitTanSipWcsTask 
    """
    MatchClass = None

    def setUp(self):
        crval = afwCoord.IcrsCoord(afwGeom.PointD(44., 45.))
        crpix = afwGeom.PointD(0, 0)
        
        arcsecPerPixel = 1/3600.0
        CD11 = arcsecPerPixel
        CD12 = 0
        CD21 = 0
        CD22 = arcsecPerPixel
        
        self.tanWcs = afwImage.makeWcs(crval, crpix, CD11, CD12, CD21, CD22)

        S = 3000
        N = 5

        if self.MatchClass == afwTable.ReferenceMatch:
            refSchema = LoadReferenceObjectsTask.makeMinimalSchema(
                filterNameList = ["r"], addFluxSigma=True, addIsPhotometric=True)
            self.refCat = afwTable.SimpleCatalog(refSchema)
        elif self.MatchClass == afwTable.SourceMatch:
            refSchema = afwTable.SourceTable.makeMinimalSchema()
            self.refCat = afwTable.SourceCatalog(refSchema)
        else:
            raise RuntimeError("Unsupported MatchClass=%r" % (self.MatchClass,))
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        SingleFrameMeasurementTask(schema=srcSchema)
        self.srcCoordKey = srcSchema["coord"].asKey()
        self.srcCentroidKey = afwTable.Point2DKey(srcSchema["slot_Centroid"])
        self.srcCentroidKey_xSigma = srcSchema["slot_Centroid_xSigma"].asKey()
        self.srcCentroidKey_ySigma = srcSchema["slot_Centroid_ySigma"].asKey()
        self.sourceCat = afwTable.SourceCatalog(srcSchema)
        self.origSourceCat = afwTable.SourceCatalog(srcSchema) # undistorted copy
        self.matches = []


        for i in numpy.linspace(0., S, N):
            for j in numpy.linspace(0., S, N):
                src = self.sourceCat.addNew()
                origSrc = self.origSourceCat.addNew()
                refObj = self.refCat.addNew()

                src.set(self.srcCentroidKey, afwGeom.Point2D(i, j))
                src.set(self.srcCentroidKey_xSigma, 0.1)
                src.set(self.srcCentroidKey_ySigma, 0.1)
                origSrc.set(self.srcCentroidKey, afwGeom.Point2D(i, j))
                origSrc.set(self.srcCentroidKey_xSigma, 0.1)
                origSrc.set(self.srcCentroidKey_ySigma, 0.1)

                c = self.tanWcs.pixelToSky(afwGeom.Point2D(i, j))
                refObj.setCoord(c)
                
                if False:
                    print "x,y = (%.1f, %.1f) pixels -- RA,Dec = (%.3f, %.3f) deg" % \
                        (i, j, c.toFk5().getRa().asDegrees(), c.toFk5().getDec().asDegrees())

                self.matches.append(self.MatchClass(refObj, src, 0.0))

    def tearDown(self):
        del self.refCat
        del self.origSourceCat
        del self.sourceCat
        del self.matches
        del self.tanWcs

    def testTrivial(self):
        """Add no distortion"""
        self.doTest("testTrivial", lambda x, y: (x, y))

    def testOffset(self):
        """Add an offset"""
        self.doTest("testOffset", lambda x, y: (x + 5, y + 7))

    def testLinearX(self):
        """Scale x, offset y"""
        self.doTest("testLinearX", lambda x, y: (2*x, y + 7))

    def testLinearXY(self):
        """Scale x and y"""
        self.doTest("testLinearXY", lambda x, y: (2*x, 3*y))

    def testLinearYX(self):
        """Add an offset to each point; scale in y and x"""
        self.doTest("testLinearYX", lambda x, y: (x + 0.2*y, y + 0.3*x))

    def testQuadraticX(self):
        """Add quadratic distortion in x"""
        self.doTest("testQuadraticX", lambda x, y: (x + 1e-5*x**2, y), order=4, specifyBBox=True)

    def testRadial(self):
        """Add radial distortion"""
        radialTransform = afwGeom.RadialXYTransform([0, 1.01, 1e-8])
        def radialDistortion(x, y):
            x, y = radialTransform.forwardTransform(afwGeom.Point2D(x, y))
            return (x, y)
        self.doTest("testRadial", radialDistortion)

    def checkResults(self, fitRes, catsUpdated):
        """Check results

        @param[in] fitRes  a object with two fields:
            - wcs  fit TAN-SIP WCS, an lsst.afw.image.TanWcs
            - scatterOnSky  median on-sky scatter, an lsst.afw.geom.Angle
        @param[in] catsUpdated  if True then coord field of self.sourceCat and centroid fields of self.refCat
            have been updated
        """
        self.assertLess(fitRes.scatterOnSky.asArcseconds(), 0.001)
        tanSipWcs = fitRes.wcs
        maxAngSep = afwGeom.Angle(0)
        maxPixSep = 0
        refCoordKey = self.refCat.schema["coord"].asKey()
        if catsUpdated:
            refCentroidKey = afwTable.Point2DKey(self.refCat.schema["centroid"])
        for refObj, src, d in self.matches:
            srcPixPos = src.get(self.srcCentroidKey)
            refCoord = refObj.get(refCoordKey)
            if catsUpdated:
                refPixPos = refObj.get(refCentroidKey)
                srcCoord = src.get(self.srcCoordKey)
            else:
                refPixPos = tanSipWcs.skyToPixel(refCoord)
                srcCoord = tanSipWcs.pixelToSky(srcPixPos)

            angSep = refCoord.angularSeparation(srcCoord)
            maxAngSep = max(maxAngSep, angSep)
            self.assertLess(angSep, 0.001 * afwGeom.arcseconds)

            pixSep = math.hypot(*(srcPixPos - refPixPos))
            maxPixSep = max(maxPixSep, pixSep)
            self.assertLess(pixSep, 0.001)

        print("max angular separation = %0.4f arcsec" % (maxAngSep.asArcseconds(),))
        print("max pixel separation = %0.3f" % (maxPixSep,))

    def doTest(self, name, func, order=3, specifyBBox=False, doPlot=False):
        """Apply func(x, y) to each source in self.sourceCat, then fit and check the resulting WCS
        """
        bbox = afwGeom.Box2I()
        for refObj, src, d in self.matches:
            origPos = src.get(self.srcCentroidKey)
            x, y = func(*origPos)
            distortedPos = afwGeom.Point2D(*func(*origPos))
            src.set(self.srcCentroidKey, distortedPos)
            bbox.include(afwGeom.Point2I(afwGeom.Point2I(distortedPos)))
        
        if specifyBBox:
            sipObject = makeCreateWcsWithSip(self.matches, self.tanWcs, order, bbox)
        else:
            sipObject = makeCreateWcsWithSip(self.matches, self.tanWcs, order)
        tanSipWcs = sipObject.getNewWcs()
        fitRes = lsst.pipe.base.Struct(
            wcs = tanSipWcs,
            scatterOnSky = sipObject.getScatterOnSky(),
        )

        if doPlot:
            self.plotWcs(tanSipWcs, name=name)
        
        self.checkResults(fitRes, catsUpdated=False)

        if self.MatchClass == afwTable.ReferenceMatch:
            fitterConfig = FitTanSipWcsTask.ConfigClass()
            fitterConfig.order = order
            fitter = FitTanSipWcsTask(config=fitterConfig)
            if specifyBBox:
                fitRes = fitter.fitWcs(
                    matches = self.matches,
                    initWcs = self.tanWcs,
                    bbox = bbox,
                    refCat = self.refCat,
                    sourceCat = self.sourceCat,
                )
            else:
                fitRes = fitter.fitWcs(
                    matches = self.matches,
                    initWcs = self.tanWcs,
                    bbox = bbox,
                    refCat = self.refCat,
                    sourceCat = self.sourceCat,
                )

            self.checkResults(fitRes, catsUpdated=True)

    def plotWcs(self, tanSipWcs, name=""):
        fileNamePrefix = "testCreateWcsWithSip_%s_%s" % (self.MatchClass.__name__, name)
        pnum = 1

        xs,ys, xc,yc = [],[],[],[]
        rs,ds, rc,dc = [],[],[],[]
        for ref, src, d in self.matches:
            xs.append(src.getX())
            ys.append(src.getY())
            refPixPos = tanSipWcs.skyToPixel(ref.getCoord())
            xc.append(refPixPos[0])
            yc.append(refPixPos[1])
            rc.append(ref.getRa())
            dc.append(ref.getDec())
            srd = tanSipWcs.pixelToSky(src.getX(), src.getY()).toFk5()
            rs.append(srd.getRa())
            ds.append(srd.getDec())
        xs = numpy.array(xs)
        ys = numpy.array(ys)
        xc = numpy.array(xc)
        yc = numpy.array(yc)
            
        pylab.clf()
        pylab.plot(xs, ys, "r.")
        pylab.plot(xc, yc, "bx")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        pylab.savefig(fileName)
        print "Wrote", fileName
        pnum += 1

        pylab.clf()
        pylab.plot(xs, xc-xs, "b.")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        pylab.xlabel("x(source)")
        pylab.ylabel("x(ref - src)")
        pylab.savefig(fileName)
        print "Wrote", fileName
        pnum += 1

        pylab.clf()
        pylab.plot(rs, ds, "r.")
        pylab.plot(rc, dc, "bx")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        pylab.savefig(fileName)
        print "Wrote", fileName
        pnum += 1

        pylab.clf()
        for y in numpy.linspace(0, 4000, 5):
            x0,y0 = [],[]
            x1,y1 = [],[]
            for x in numpy.linspace(0., 4000., 401):
                x0.append(x)
                y0.append(y)
                rd = tanSipWcs.pixelToSky(x, y)
                xy = tanSipWcs.skyToPixel(rd)
                x1.append(xy[0])
                y1.append(xy[1])
            x0 = numpy.array(x0)
            x1 = numpy.array(x1)
            pylab.plot(x0, x1-x0, "b-")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        pylab.savefig(fileName)
        print "Wrote", fileName
        pnum += 1

def makeTestCase(_MatchClass):
    class CreateWcsWithSipTestCase(BaseTestCase):
        MatchClass = _MatchClass
    return CreateWcsWithSipTestCase

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(makeTestCase(afwTable.ReferenceMatch))
    suites += unittest.makeSuite(makeTestCase(afwTable.SourceMatch))
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
