# 
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

import fnmatch
import numpy as np

from lsst.pex.config import Config, Field, ConfigDictField


class ColortermNotFoundError(LookupError):
    """Exception class indicating we couldn't find a colorterm"""
    pass



class Colorterm(object):
    """A class to describe colour terms between photometric bands"""

    def __init__(self, primary, secondary, c0=0.0, c1=0.0, c2=0.0):
        """The transformed magnitude p' is given by
        p' = primary + c0 + c1*(primary - secondary) + c2*(primary - secondary)**2
        """
        self.primary = primary
        self.secondary = secondary
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2

    def __str__(self):
        return "Colorterm(%s, %s, %g, %g, %g)" % (self.primary, self.secondary, self.c0, self.c1, self.c2)

    def transformSource(self, source, reverse=False):
        """Transform the magnitudes in <source> and return it. The <source> must
support a get(band) (e.g. source.get("r")) method, as do afw::Source and dicts.
If reverse is True, return the inverse transformed magnitude
        """
        return self.transformMags(source.get(self.primary), source.get(self.secondary), reverse)

    def transformMags(self, primary, secondary, reverse=False):
        """Transform the magnitudes <primary> and <secondary> and return it.
If reverse is True, return the inverse transformed magnitude
        """
        if reverse:
            raise NotImplemented("reverse photometric transformations are not implemented")
        color = primary - secondary
        return primary + self.c0 + color*(self.c1 + color*self.c2)

    def propagateFluxErrors(self, primaryFluxErr, secondaryFluxErr, reverse=False):
        return np.hypot((1 + self.c1)*primaryFluxErr, self.c1*secondaryFluxErr)


class ColortermConfig(Config):
    """
    The transformed magnitude p' is given by
        p' = primary + c0 + c1*(primary - secondary) + c2*(primary - secondary)**2
    """
    primary = Field(dtype=str, doc="Primary band")
    secondary = Field(dtype=str, doc="Secondary band")
    c0 = Field(dtype=float, default=0.0, doc="Constant parameter")
    c1 = Field(dtype=float, default=0.0, doc="First-order parameter")
    c2 = Field(dtype=float, default=0.0, doc="Second-order parameter")

    @classmethod
    def fromValues(cls, primary, secondary, c0=0.0, c1=0.0, c2=0.0):
        self = cls()
        self.primary = primary
        self.secondary = secondary
        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        return self

    def getColorterm(self):
        """Return a Colorterm object based on this Config"""
        return Colorterm(self.primary, self.secondary, self.c0, self.c1, self.c2)


class ColortermGroupConfig(Config):
    """A group of color terms, e.g., for a particular camera used with a particular reference catalog

    Can't name the element "set" or "list", as these are python classes.
    """
    group = ConfigDictField(keytype=str, itemtype=ColortermConfig, default={},
                             doc="Mapping of band name to color terms")

    @classmethod
    def fromValues(cls, group):
        self = cls()
        for k,v in group.iteritems():
            self.group[k] = v
        return self


class ColortermLibraryConfig(Config):
    """A library of colorterm groups, for a particular camera with a variety of reference catalogs"""
    library = ConfigDictField(keytype=str, itemtype=ColortermGroupConfig, default={},
                              doc="Mapping of astrometry_net_data version (or glob) to group of color terms")

    @classmethod
    def fromValues(cls, library):
        self = cls()
        for k,v in library.iteritems():
            self.library[k] = v
        return self

    def selectColorTerm(self, filterName, version=None):
        """Select the appropriate Colorterm from the library

        We use the group of color terms in the library that matches the version.
        The version defaults to the astrometry_net_data version currently setup
        in eups.  If the version exactly matches an entry in the library, that
        group is used; otherwise if the version matches a single glob (shell syntax,
        e.g., "sdss-*" will match "sdss-dr8"), then that is used.  If there is no
        exact match and no unique match to the globs, we raise an exception.

        @param filterName: name of filter
        @param version: astrometry_net_data version name, or None to ask eups
        @return the appropriate Colorterm
        """
        if version is None:
            import eups
            version = eups.Eups().findSetupVersion("astrometry_net_data")[0]
        if version in self.library:
            group = self.library[version]
        else:
            matchList = [glob for glob in self.library if fnmatch.fnmatch(version, glob)]
            if len(matchList) > 1:
                raise RuntimeError("Multiple library globs match astrometry_net_data version %s: %s" %
                                   (version, matchList))
            if len(matchList) == 0:
                raise RuntimeError("No colorterm group matches astrometry_net_data version %s" % version)
            group = self.library[matchList.pop()].group
        if filterName not in group:
            # Check it's not an alias
            from lsst.afw.image import Filter
            filterName = Filter(Filter(filterName).getId()).getName()
            if filterName not in group:
                raise ColortermNotFoundError("No colorterm set for filter %s with "
                                             "astrometry_net_data version %s" %
                                             (filterName, version))
        return group[filterName].getColorterm()

    def getColorTerm(self, filterName, version=None):
        """Wrapper around selectColorTerm to always provide a color term

        A default color term (with zero-valued coefficients) is provided when one is not available.
        """
        try:
            ct = self.selectColorTerm(filterName, version=version)
        except ColortermNotFoundError as e:
            from lsst.pex.logging import getDefaultLog
            getDefaultLog().warn(str(e))
            ct = Colorterm(filterName, filterName)
        return ct

