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
import math, os, sys
import numpy as np

import lsst.meas.algorithms.utils as malgUtil
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConf
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9

from .colorterms import ColortermLibraryConfig

def checkSourceFlags(source, keys):
    """Return True if the given source has all good flags set and none of the bad flags set.

    @param[in] source    SourceRecord object to process.
    @param[in] keys      Struct of source catalog keys, as returned by PhotCalTask.getKeys()
    """
    for k in keys.goodFlags:
        if not source.get(k): return False
    for k in keys.badFlags:
        if source.get(k): return False
    return True

class PhotoCalConfig(pexConf.Config):

    magLimit = pexConf.Field(dtype=float, doc="Don't use objects fainter than this magnitude", default=22.0)
    outputField = pexConf.Field(
        dtype=str, optional=True, default="classification.photometric",
        doc="Name of the flag field that is set for sources used in photometric calibration"
        )
    fluxField = pexConf.Field(
        dtype=str, default="flux.sinc", optional=False,
        doc="Name of the source flux field to use.  The associated flag field\n"\
            "('<name>.flags') will be implicitly included in badFlags.\n"
        )
    applyColorTerms = pexConf.Field(
        dtype=bool, default=True,
        doc= "Apply photometric colour terms (if available) to reference stars",
        )
    goodFlags = pexConf.ListField(
        dtype=str, optional=False,
        default=[], 
        doc="List of source flag fields that must be set for a source to be used."
        )
    badFlags = pexConf.ListField(
        dtype=str, optional=False,
        default=["flags.pixel.edge", "flags.pixel.interpolated.any", "flags.pixel.saturated.any"], 
        doc="List of source flag fields that will cause a source to be rejected when they are set."
        )
    sigmaMax = pexConf.Field(dtype=float, default=0.25, optional=True,
                              doc="maximum sigma to use when clipping")
    nSigma = pexConf.Field(dtype=float, default=3.0, optional=False, doc="clip at nSigma")
    useMedian = pexConf.Field(dtype=bool, default=True,
                              doc="use median instead of mean to compute zeropoint")
    nIter = pexConf.Field(dtype=int, default=20, optional=False, doc="number of iterations")
    colorterms = pexConf.ConfigField(dtype=ColortermLibraryConfig, doc="Library of color terms")

class PhotoCalTask(pipeBase.Task):
    """Calculate the zero point of an exposure given a ReferenceMatchVector.
    """
    ConfigClass = PhotoCalConfig
    _DefaultName = "photoCal"

    def __init__(self, schema, **kwds):
        """Create the task, pulling input keys from the schema and adding a flag field
        'classification.photometric' that will be set for sources used to determine the
        photometric calibration.
        """
        pipeBase.Task.__init__(self, **kwds)
        if self.config.outputField is not None:
            self.output = schema.addField(self.config.outputField, type="Flag",
                                          doc="set if source was used in photometric calibration")
        else:
            self.output = None

    def getKeys(self, schema, prefix=""):
        """Return a struct containing the source catalog keys for fields used by PhotoCalTask."""
        flux = schema.find(prefix + self.config.fluxField).key
        fluxErr = schema.find(prefix + self.config.fluxField + ".err").key
        goodFlags = [schema.find(prefix + name).key for name in self.config.goodFlags]
        badFlags = [schema.find(prefix + self.config.fluxField + ".flags").key]
        badFlags.extend(schema.find(prefix + name).key for name in self.config.badFlags)
        return pipeBase.Struct(flux=flux, fluxErr=fluxErr, goodFlags=goodFlags, badFlags=badFlags)

    @pipeBase.timeMethod
    def selectMatches(self, matches, keys, frame=None):
        """Select reference/source matches according the criteria specified in the config.

        if frame is non-None, display information about trimmed objects on that ds9 frame:
            Bad:               red x
            Non-"photometric": blue +  (and a cyan o if a galaxy)
            Failed flux cut:   magenta *

        The return value is a ReferenceMatchVector that contains only the selected matches.
        If a schema was passed during task construction, a flag field will be set on sources 
        in the selected matches.

        An exception will be raised if there are no valid matches.

        @param[in] matches ReferenceMatchVector (not modified)
        @param[in] keys    Struct of source catalog keys, as returned by getKeys()
        @param[in] frame   ds9 frame number to use for debugging display
        @return Output ReferenceMatchVector
        """

        self.log.logdebug("Number of input matches: %d" % (len(matches)))
        if len(matches) == 0:
            raise ValueError("No input matches")

        # Only use stars for which the flags indicate the photometry is good.
        afterFlagCutInd = [i for i, m in enumerate(matches) if checkSourceFlags(m.second, keys)]
        afterFlagCut = [matches[i] for i in afterFlagCutInd]
        self.log.logdebug("Number of matches after source flag cuts: %d" % (len(afterFlagCut)))

        if len(afterFlagCut) != len(matches):
            if frame is not None:
                with ds9.Buffering():
                    for i, m in enumerate(matches):
                        if i not in afterFlagCutInd:
                            x, y = m.second.getCentroid()
                            ds9.dot("x", x,  y, size=4, frame=frame, ctype=ds9.RED)

            matches = afterFlagCut

        if len(matches) == 0:
            raise ValueError("All matches eliminated by source flags")

        refSchema = matches[0].first.schema
        try:
            refKey = refSchema.find("photometric").key
            try:
                stargalKey = refSchema.find("stargal").key
            except:
                stargalKey = None

            try:
                varKey = refSchema.find("var").key
            except:
                varKey = None
        except:
            self.log.warn("No 'photometric' flag key found in reference schema.")
            refKey = None

        if refKey is not None:
            afterRefCutInd = [i for i, m in enumerate(matches) if m.first.get(refKey)]
            afterRefCut = [matches[i] for i in afterRefCutInd]

            if len(afterRefCut) != len(matches):
                if frame is not None:
                    with ds9.Buffering():
                        for i, m in enumerate(matches):
                            if i not in afterRefCutInd:
                                x, y = m.second.getCentroid()
                                ds9.dot("+", x,  y, size=4, frame=frame, ctype=ds9.BLUE)

                                if stargalKey and not m.first.get(stargalKey):
                                    ds9.dot("o", x,  y, size=6, frame=frame, ctype=ds9.CYAN)
                                if varKey and m.first.get(varKey):
                                    ds9.dot("o", x,  y, size=6, frame=frame, ctype=ds9.MAGENTA)

                matches = afterRefCut

        self.log.logdebug("Number of matches after reference catalog cuts: %d" % (len(matches)))
        if len(matches) == 0:
            raise RuntimeError("No sources remain in match list after reference catalog cuts.")
        
        fluxKey = refSchema.find("flux").key
        if self.config.magLimit is not None:
            fluxLimit = 10.0**(-self.config.magLimit/2.5)

            afterMagCutInd = [i for i, m in enumerate(matches) if (m.first.get(fluxKey) > fluxLimit
                                                                   and m.second.get(keys.flux) > 0.0)]
        else:
            afterMagCutInd = [i for i, m in enumerate(matches) if m.second.get(keys.flux) > 0.0]

        afterMagCut = [matches[i] for i in afterMagCutInd]

        if len(afterMagCut) != len(matches):
            if frame is not None:
                with ds9.Buffering():
                    for i, m in enumerate(matches):
                        if i not in afterMagCutInd:
                            x, y = m.second.getCentroid()
                            ds9.dot("*", x,  y, size=4, frame=frame, ctype=ds9.MAGENTA)

            matches = afterMagCut
            
        self.log.logdebug("Number of matches after magnitude limit cuts: %d" % (len(matches)))

        if len(matches) == 0:
            raise RuntimeError("No sources remaining in match list after magnitude limit cuts.")

        if frame is not None:
            with ds9.Buffering():
                for m in matches:
                    x, y = m.second.getCentroid()
                    ds9.dot("o", x,  y, size=4, frame=frame, ctype=ds9.GREEN)

        result = afwTable.ReferenceMatchVector()
        for m in matches:
            if self.output is not None:
                m.second.set(self.output, True)
            result.append(m)
        return result

    @pipeBase.timeMethod
    def extractMagArrays(self, matches, filterName, keys, version=None):
        """Extract magnitude and magnitude error arrays from the given matches.

        @param[in] matches    ReferenceMatchVector object containing reference/source matches
        @param[in] filterName Name of filter being calibrated
        @param[in] keys       Struct of source catalog keys, as returned by getKeys()
        @param[in] version    astrometry_net_data version, or None (for color terms)

        @return Struct containing srcMag, refMag, srcMagErr, refMagErr, and errMag arrays.
        """
        srcFlux = np.array([m.second.get(keys.flux) for m in matches])
        srcFluxErr = np.array([m.second.get(keys.fluxErr) for m in matches])
        if not np.all(np.isfinite(srcFluxErr)):
            self.log.warn("Source catalog does not have flux uncertainties; using sqrt(flux).")
            srcFluxErr = np.sqrt(srcFlux)
        
        if not matches:
            raise RuntimeError("No reference stars are available")

        refSchema = matches[0].first.schema
        if self.config.applyColorTerms:
            ct = self.config.colorterms.selectColorTerm(filterName, version=version)
        else:
            ct = None

        if ct:                          # we have a colour term to worry about
            fluxNames = [ct.primary, ct.secondary]
            missingFluxes = [flux for flux in fluxNames if not flux in refSchema]
            if missingFluxes:
                self.log.warn("Source catalog does not have fluxes for %s; ignoring color terms" %
                              " ".join(missingFluxes))
                ct = None
                
        if not ct:
            if filterName in refSchema:
                fluxNames = [filterName]
            else:
                fluxNames = ["flux"]

        refFluxes = []
        refFluxErrors = []
        for flux in fluxNames:
            fluxKey = refSchema.find(flux).key
            refFlux = np.array([m.first.get(fluxKey) for m in matches])
            try:
                fluxErrKey = refSchema.find(flux + ".err").key
                refFluxErr = np.array([m.first.get(fluxErrKey) for m in matches])
            except KeyError:
                # Catalogue may not have flux uncertainties; HACK
                self.log.warn("Reference catalog does not have flux uncertainties for %s; using sqrt(flux)."
                              % flux)
                refFluxErr = np.sqrt(refFlux)

            refFluxes.append(refFlux)
            refFluxErrors.append(refFluxErr)

        srcMag = -2.5*np.log10(srcFlux)
        if ct:                          # we have a colour term to worry about
            refMag =  -2.5*np.log10(refFluxes[0]) # primary
            refMag2 = -2.5*np.log10(refFluxes[1]) # secondary

            refMag = ct.transformMags(refMag, refMag2)
            refFluxErr = ct.propagateFluxErrors(refFluxErrors[0], refFluxErrors[1])
        else:
            refMag = -2.5*np.log10(refFluxes[0])

        # Fitting with error bars in both axes is hard, so transfer all
        # the error to src, then convert to magnitude
        fluxErr = np.hypot(srcFluxErr, refFluxErr)
        magErr = fluxErr / (srcFlux * np.log(10))

        srcMagErr = srcFluxErr / (srcFlux * np.log(10))
        refMagErr = refFluxErr / (refFlux * np.log(10))

        return pipeBase.Struct(
            srcMag = srcMag,
            refMag = refMag,
            magErr = magErr,
            srcMagErr = srcMagErr,
            refMagErr = refMagErr
            )

    @pipeBase.timeMethod
    def run(self, exposure, matches, prefix=""):
        """Do photometric calibration - select matches to use and (possibly iteratively) compute
        the zero point.

        @param[in]  exposure   Exposure upon which the sources in the matches were detected.
        @param[in]  matches    Input ReferenceMatchVector (will not be modified).
        @param[in]  prefix     Field name prefix to be included in all flag field names.
        
        @return Struct of:
           calib ------- Calib object containing the zero point
           arrays ------ Magnitude arrays returned be extractMagArrays
           matches ----- Final ReferenceMatchVector, as returned by selectMathces.
        """
        global scatterPlot, fig
        import lsstDebug

        display = lsstDebug.Info(__name__).display
        displaySources = display and lsstDebug.Info(__name__).displaySources
        scatterPlot = display and lsstDebug.Info(__name__).scatterPlot

        if scatterPlot:
            from matplotlib import pyplot
            try:
                fig.clf()
            except:
                fig = pyplot.figure()

        if displaySources:
            frame = 1
            ds9.mtv(exposure, frame=frame, title="photocal")
        else:
            frame = None

        keys = self.getKeys(matches[0].second.schema, prefix=prefix)
        matches = self.selectMatches(matches, keys, frame=frame)
        arrays = self.extractMagArrays(matches, exposure.getFilter().getName(), keys)

        # Fit for zeropoint.  We can run the code more than once, so as to
        # give good stars that got clipped by a bad first guess a second
        # chance.
        # FIXME: these should be config values

        calib = afwImage.Calib()
        zp = None                           # initial guess
        r = self.getZeroPoint(arrays.srcMag, arrays.refMag, arrays.magErr, zp0=zp)
        zp = r.zp
        self.log.info("Magnitude zero point: %f +/- %f from %d stars" % (r.zp, r.sigma, r.ngood))

        flux0 = 10**(0.4*r.zp) # Flux of mag=0 star
        flux0err = 0.4*math.log(10)*flux0*r.sigma # Error in flux0
        
        calib.setFluxMag0(flux0, flux0err)

        return pipeBase.Struct(
            calib = calib,
            arrays = arrays,
            matches = matches,
            zp = r.zp,
            sigma = r.sigma,
            ngood = r.ngood,
            )

    def getZeroPoint(self, src, ref, srcErr=None, zp0=None):
        """Flux calibration code, returning (ZeroPoint, Distribution Width, Number of stars)

        We perform nIter iterations of a simple sigma-clipping algorithm with a a couple of twists:
        1.  We use the median/interquartile range to estimate the position to clip around, and the
        "sigma" to use.
        2.  We never allow sigma to go _above_ a critical value sigmaMax --- if we do, a sufficiently
        large estimate will prevent the clipping from ever taking effect.
        3.  Rather than start with the median we start with a crude mode.  This means that a set of magnitude
        residuals with a tight core and asymmetrical outliers will start in the core.  We use the width of
        this core to set our maximum sigma (see 2.)  
        """

        sigmaMax = self.config.sigmaMax

        dmag = ref - src

        i = np.argsort(dmag)
        dmag = dmag[i]
        
        if srcErr is not None:
            dmagErr = srcErr[i]
        else:
            dmagErr = np.ones(len(dmag))

        # need to remove nan elements to avoid errors in stats calculation with numpy
        ind_noNan = np.array([ i for i in range(len(dmag)) if (not np.isnan(dmag[i]) and not np.isnan(dmagErr[i])) ])
        dmag = dmag[ind_noNan]
        dmagErr = dmagErr[ind_noNan]

        IQ_TO_STDEV = 0.741301109252802;    # 1 sigma in units of interquartile (assume Gaussian)

        npt = len(dmag)
        ngood = npt
        if scatterPlot > 1:
            reply = None
        for i in range(self.config.nIter):
            if i > 0:
                npt = sum(good)

            center = None
            if i == 0:
                #
                # Start by finding the mode
                #
                nhist = 20
                try:
                    hist, edges = np.histogram(dmag, nhist, new=True)
                except TypeError:
                    hist, edges = np.histogram(dmag, nhist) # they removed new=True around numpy 1.5
                imode = np.arange(nhist)[np.where(hist == hist.max())]

                if imode[-1] - imode[0] + 1 == len(imode): # Multiple modes, but all contiguous
                    if zp0:
                        center = zp0
                    else:
                        center = 0.5*(edges[imode[0]] + edges[imode[-1] + 1])

                    peak = sum(hist[imode])/len(imode) # peak height

                    # Estimate FWHM of mode
                    j = imode[0]
                    while j >= 0 and hist[j] > 0.5*peak:
                        j -= 1
                    j = max(j, 0)
                    q1 = dmag[sum(hist[range(j)])]

                    j = imode[-1]
                    while j < nhist and hist[j] > 0.5*peak:
                        j += 1
                    j = min(j, nhist - 1)
                    j = min(sum(hist[range(j)]), npt - 1)
                    q3 = dmag[j]

                    if q1 == q3:
                        q1 = dmag[int(0.25*npt)]
                        q3 = dmag[int(0.75*npt)]

                    sig = (q3 - q1)/2.3 # estimate of standard deviation (based on FWHM; 2.358 for Gaussian)

                    if sigmaMax is None:
                        sigmaMax = 2*sig   # upper bound on st. dev. for clipping. multiplier is a heuristic

                    self.log.logdebug("Photo calibration histogram: center = %.2f, sig = %.2f" 
                                      % (center, sig))

                else:
                    if sigmaMax is None:
                        sigmaMax = dmag[-1] - dmag[0]

                    center = np.median(dmag)
                    q1 = dmag[int(0.25*npt)]
                    q3 = dmag[int(0.75*npt)]
                    sig = (q3 - q1)/2.3 # estimate of standard deviation (based on FWHM; 2.358 for Gaussian)

            if center is None:              # usually equivalent to (i > 0)
                gdmag = dmag[good]
                if self.config.useMedian:
                    center = np.median(gdmag)
                else:
                    gdmagErr = dmagErr[good]
                    center = np.average(gdmag, weights=gdmagErr)

                q3 = gdmag[min(int(0.75*npt + 0.5), npt - 1)]
                q1 = gdmag[min(int(0.25*npt + 0.5), npt - 1)]

                sig = IQ_TO_STDEV*(q3 - q1)     # estimate of standard deviation

            good = abs(dmag - center) < self.config.nSigma*min(sig, sigmaMax) # don't clip too softly

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            if scatterPlot:
                from matplotlib import pyplot
                try:
                    fig.clf()

                    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

                    axes.plot(ref[good], dmag[good] - center, "b+")
                    axes.errorbar(ref[good], dmag[good] - center, yerr=dmagErr[good], 
                                  linestyle='', color='b')

                    bad = np.logical_not(good)
                    if len(ref[bad]) > 0:
                        axes.plot(ref[bad], dmag[bad] - center, "r+")
                        axes.errorbar(ref[bad], dmag[bad] - center, yerr=dmagErr[bad], 
                                      linestyle='', color='r')

                    axes.plot((-100, 100), (0, 0), "g-")
                    for x in (-1, 1):
                        axes.plot((-100, 100), x*0.05*np.ones(2), "g--")

                    axes.set_ylim(-1.1, 1.1)
                    axes.set_xlim(24, 13)
                    axes.set_xlabel("Reference")
                    axes.set_ylabel("Reference - Instrumental")

                    fig.show()

                    if scatterPlot > 1:
                        while reply != "c":
                            try:
                                reply = raw_input("Next iteration? [ynhpc] ")
                            except EOFError:
                                reply = "n"

                            if reply == "h":
                                print >> sys.stderr, "Options: c[ontinue] h[elp] n[o] p[db] y[es]"
                                continue

                            if reply in ("", "c", "n", "p", "y"):
                                break
                            else:
                                print >> sys.stderr, "Unrecognised response: %s" % reply

                        if reply == "n":
                            break
                        elif reply == "p":
                            import pdb; pdb.set_trace()
                except Exception, e:
                    print >> sys.stderr, "Error plotting in PhotoCal.getZeroPoint: %s" % e

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            old_ngood = ngood
            ngood = sum(good)
            if ngood == 0:
                msg = "PhotoCal.getZeroPoint: no good stars remain"

                if i == 0:                  # failed the first time round -- probably all fell in one bin
                    center = np.average(dmag, weights=dmagErr)
                    msg += " on first iteration; using average of all calibration stars"

                self.log.log(self.log.WARN, msg)

                return pipeBase.Struct(
                    zp = center,
                    sigma = sig,
                    ngood = len(dmag)
                    )
            elif ngood == old_ngood:
                break

            if False:
                ref = ref[good]
                dmag = dmag[good]
                dmagErr = dmagErr[good]

        dmag = dmag[good]
        dmagErr = dmagErr[good]
        return pipeBase.Struct(
            zp = np.average(dmag, weights=dmagErr),
            sigma = np.std(dmag, ddof=1),
            ngood = len(dmag)
            )
