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


import os    
import pdb
import math
import numpy as np

import eups
import lsst.meas.astrom as measAst
import lsst.meas.astrom.sip as sip
import lsst.meas.algorithms.utils as malgUtil
from lsst.pex.logging import Log, Debug, LogRec, Prop

from lsst.meas.photocal.PhotometricMagnitude import PhotometricMagnitude

class fakelog(object):
    def log(self, x, msg):
        print msg

def calcPhotoCal(sourceMatch, log=None,
                 goodFlagValue=(malgUtil.getDetectionFlags()['BINNED1'] |
                                malgUtil.getDetectionFlags()['STAR'])):
    """Calculate photometric calibration, i.e the zero point magnitude"""

    if log is None:
        log = fakelog()
        
    if len(sourceMatch) == 0:
        raise ValueError("sourceMatch contains no elements")
        
    # Only use stars for which the flags indicate the photometry is good.
    allFlags = malgUtil.getDetectionFlags()

    for flag, val in allFlags.items():
        if val & goodFlagValue:
            log.log(Log.INFO, 'Flag value %s: 0x%x, %i' % (flag, val, val))
    log.log(Log.INFO, 'Good flag value: %i' % goodFlagValue)
    log.log(Log.INFO, "Number of sources: %d" % (len(sourceMatch)))

    cleanList = []
    for m in sourceMatch:
        if m.second.getFlagForDetection() & goodFlagValue == goodFlagValue:
            cleanList.append(m)
    sourceMatch = cleanList
    
    if len(sourceMatch) == 0:
        raise ValueError("flags indicate all elements of sourceMatch have bad photometry")
    
    #Convert fluxes to magnitudes
    out = getMagnitudes(sourceMatch)

    #Fit to get zeropoint
    lsf = robustFit(out["src"], out["cat"], order=2, plot=True)
    par = lsf.getParams()
    err = lsf.getErrors()
    
    #Sanity check output
    if not(par[1] - err[1] <= 1 and par[1] + err[1] >= 1):
        msg = "Slope of fitting function is not 1 (%g +- %g) " %(par[1], err[1])
        log.log(Log.WARN, msg)
    
    #Initialise and return a magnitude object
    medianInstMag = np.median(out["src"])
    medianFlux = np.power(10, -medianInstMag/2.5)
    mag = lsf.valueAt(medianInstMag)
    zp = mag - medianInstMag
    
    #return PhotometricMagnitude(zeroFlux=medianFlux, zeroMag=mag)
    return PhotometricMagnitude(zeroFlux=1.0, zeroMag=zp)


def getMagnitudes(sourceMatch):
    
    #Extract the fluxes as numpy arrays. Catalogues don't come with fluxes
    #so do an estimate using sqrt(counts)
    fluxCat = np.array(map(lambda x: x.first.getPsfFlux(), sourceMatch))
    fluxSrc = np.array(map(lambda x: x.second.getPsfFlux(), sourceMatch))
    fluxSrcErr = np.array(map(lambda x: x.second.getPsfFluxErr(), sourceMatch))
    
    
    #Remove objects where the source flux is bad
    idx = np.where(fluxSrc > 0)

    fluxSrc = fluxSrc[idx]
    fluxCat = fluxCat[idx]
    fluxSrcErr = fluxSrcErr[idx]
    
    #Remove the (unlikely) objects with bad catalogue fluxes
    idx = np.where(fluxCat > 0)
    fluxSrc = fluxSrc[idx]
    fluxCat = fluxCat[idx]
    fluxSrcErr = fluxSrcErr[idx]
    
    #Catalogue won't have flux uncertainties in general
    fluxCatErr = np.sqrt(fluxCat)
    
    #Convert to mags
    magSrc = -2.5*np.log10(fluxSrc)
    magCat = -2.5*np.log10(fluxCat)

    #Fitting with error bars in both axes is hard, so transfer all the error to 
    #src, then convert to magnitude
    fluxSrcErr = np.hypot(fluxSrcErr, fluxCatErr)
    magSrcErr = fluxSrcErr/fluxSrc/np.log(10)
    
    #I need to return three arrays, but am bound to get the order 
    #confused at some point, so use a dictionary instead
    out = dict()
    out["src"] = magSrc
    out["cat"] = magCat
    out["srcErr"] = magSrcErr
    
    return out
    


def robustFit(x, y, order=2, plot=False):
    """\brief Fit a polynomial to a dataset in a manner that is highly insensitive to outliers
    
    Proceedure is to bin the data into order+1 points. The x value of the bin is the mean x value
    of the points in the bin, and the y value of the bin is the *median* of the y values of the 
    points in the bin. This approach is very resistant to outliers affecting the fit.
    
    Input
    \param x       Array of ordinate values to fit
    \param y       Array of co-ordinate values
    \param order=2 Order of fit. Default (2) means to fit a straight line 
    """
    
    if len(x) == 0:
        raise ValueError("Input x array has zero length")
    
    if len(x) != len(y):
        raise ValueError("Input x and y arrays are of different length")
        
    if order <= 0:
        raise ValueError("Order must be >=1")
        
    if order > len(x)/3:
        #Hard to discriminate against outliers with only two points per bin
        raise ValueError("Order can be no greater than one third the number of data points")
        
    nBins = order+1
    idx = x.argsort()   #indices of the sorted array of x
    
    rx = chooseRobustCoord(x, idx, nBins)
    ry = chooseRobustCoord(y, idx, nBins)
    rs = np.ones(nBins)


    #if plot:
        #mpl.plot(x, y, 'ro')
        #mpl.plot(rx, ry, 'ks-')
        #lsf = sip.LeastSqFitter1dPoly(list(rx), list(ry), list(rs), order)
        #print lsf.getParams()
        #
        #fitx = range(-16, -9)
        #fity = map(lambda x: lsf.valueAt(x), fitx)
        #mpl.plot(fitx, fity, 'g-')
        #
        #mpl.show()
        
    return sip.LeastSqFitter1dPoly(list(rx), list(ry), list(rs), order)
    


def chooseRobustCoord(x, idx, nBins):
    """\brief Create nBins values of the ordinate based on the mean of groups of elements of x
    
    Inputs:
    \param x Ordinate to be binned
    \param idx Indices of x in sorted order, i.e x[idx[i]] <= x[idx[i+1]]
    \param nBins Number of desired bins
    """

    if len(x) == 0:
        raise ValueError("x array has no data")
        
    if len(x) != len(idx):
        raise ValueError("Length of x and idx don't agree")
        
    if nBins < 1:
        raise ValueError("nBins < 1")
        
    rSize = len(idx)/float(nBins)  #Note, a floating point number
    rx = np.zeros(nBins)
    
    for i in range(nBins):
        rng = range(int(rSize*i), int(rSize*(i+1)))
        rx[i] = np.median(x[idx[rng]])
    return rx
    

def clean(x, y, order=2, sigmaClip=3, maxIter=5):
    """\brief Remove outliers from the set of {(x,y)}
    
    Robust-fits a polynomial to y(x) and remove points that like more than sigmaClip times
    the scatter away from the line. Repeat until no more points are removed, or maxIter is reached
    The arrays x and y are modified.
    
    \param x ordinate array (numpy array)
    \param y coordinate array (numpy array)
    \param order  Order of polynomial to fit
    \param sigmaClip. Remove points more this number times the variance from the fit
    \param maxIter  Maximum number of iterations
    """
    
    if len(x) == 0:
        raise ValueError("Input x array has zero length")
    
    if len(x) != len(y):
        raise ValueError("Input x and y arrays are of different length")
        
    if order <= 0:
        raise ValueError("Order must be >=1")
        
    if order > len(x)/3:
        #Hard to discriminate against outliers with only two points per bin
        raise ValueError("Order can be no greater than one third the number of data points")

    if sigmaClip<0:
        raise ValueError("sigmaClip must be >=1")
        
    if maxIter<0:
        raise ValueError("maxIter must be >=1")
        
    i=0
    newSize = len(x)
    oldSize = newSize+1
    while newSize < oldSize and i<maxIter:
        lsf = robustFit(x,y,order)
        f = map(lambda x: lsf.valueAt(float(x)), x)
        
        
        sigma = (y-f).std()
        if sigma == 0:  #A perfect fit. Something odd with the data, but not our concern
            break
            
        deviance = np.fabs( (y - f) /sigma)
        idx = np.where(deviance < sigmaClip)

        x=x[idx]
        y=y[idx]
        
        oldSize=newSize
        newSize = len(x)
        
