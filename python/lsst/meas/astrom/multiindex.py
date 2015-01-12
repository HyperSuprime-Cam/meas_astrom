import os
import numpy
import pyfits
import itertools
from lsst.pex.logging import getDefaultLog
from .astrometry_net import multiindex_new, multiindex_add_index, INDEX_ONLY_LOAD_METADATA, healpixDistance
from .config import AstrometryNetDataConfig


def getIndexPath(fn):
    """Get absolute path for index file

    The ASTROMETRY_NET_DATA_DIR is prepended to a relative path.
    """
    if os.path.isabs(fn):
        return fn
    andir = os.getenv('ASTROMETRY_NET_DATA_DIR')
    if andir is not None:
        fn2 = os.path.join(andir, fn)
        if os.path.exists(fn2):
            return fn2
    if os.path.exists(fn):
        return os.path.abspath(fn)
    raise RuntimeError("Unable to determine path for index file")


class MultiIndexCache(object):
    """A wrapper for the multiindex_t, which only reads the data when it needs to

    The MultiIndexCache may be instantiated directly, or via the 'fromFilenameList'
    class method, which loads it from a list of filenames.
    """
    def __init__(self, filenameList, healpix, nside):
        """Constructor

        @param filenameList  List of filenames; first is the multiindex, then
                             follows the individual index files
        @param healpix       Healpix number
        @param nside         Healpix nside
        """
        self._filenameList = filenameList
        self._healpix = healpix
        self._nside = nside
        self._mi = None
        self._loaded = False
        self.log = getDefaultLog()

    @classmethod
    def fromFilenameList(cls, filenameList):
        """Construct from a list of filenames

        The list of filenames should contain the multiindex filename first,
        then the individual index filenames.  The healpix and nside are
        determined by reading the indices, so this is not very efficient.
        """
        self = cls(filenameList, 0, 0)
        self.read()
        healpix = set(self[i].healpix for i in range(len(self)))
        nside = set(self[i].hpnside for i in range(len(self)))
        assert len(healpix) == 1
        assert len(nside) == 1
        self._healpix = healpix.pop()
        self._nside = nside.pop()
        return self

    def read(self):
        """Read the indices"""
        if self._mi is not None:
            return
        fn = getIndexPath(self._filenameList[0])
        self._mi = multiindex_new(fn)
        if self._mi is None:
            raise RuntimeError('Failed to read stars from multiindex filename "%s"' % fn)
        for i, fn in enumerate(self._filenameList[1:]):
            if fn is None:
                self.log.logdebug('Unable to find index part of multiindex file %s' % fn)
                continue
            fn = getIndexPath(fn)
            self.log.logdebug('Reading index from multiindex file "%s"' % fn)
            if multiindex_add_index(self._mi, fn, INDEX_ONLY_LOAD_METADATA):
                raise RuntimeError('Failed to read index from multiindex filename "%s"' % fn)
            ind = self._mi[i]
            self.log.logdebug('  index %i, hp %i (nside %i), nstars %i, nquads %i' %
                              (ind.indexid, ind.healpix, ind.hpnside, ind.nstars, ind.nquads))

    def reload(self):
        """Reload the indices."""
        if self._loaded:
            return
        if self._mi is None:
            self.read()
        else:
            self._mi.reload()
        self._loaded = True

    def unload(self):
        """Unload the indices"""
        if not self._loaded:
            return
        self._mi.unload()
        self._loaded = False

    def isWithinRange(self, ra, dec, distance):
        """Is the index within range of the provided coordinates?

        ra, dec and distance are in degrees (astrometry.net standard).
        """
        return (self._nside == 0 or
                healpixDistance(self._healpix, self._nside, ra, dec).asDegrees() <= distance)

    def __getitem__(self, i):
        self.reload()
        return self._mi[i]

    def __len__(self):
        return len(self._filenameList) - 1 # The first is the multiindex; the rest are the indices

    def __iter__(self):
        self.reload()
        return iter(self._mi)


class AstrometryNetCatalog(object):
    """An interface to an astrometry.net catalog

    These should usually be constructed using the 'fromEnvironment'
    class method, which wraps the 'fromIndexFiles' and 'fromCache'
    alternative class methods.
    """

    def __init__(self, multiInds):
        """Constructor

        @param multiInds  List of multiindex objects
        """
        self._multiInds = multiInds

    @classmethod
    def fromIndexFiles(cls, andConfig):
        """Construct from the index files in an AstrometryNetDataConfig"""
        indexFiles = zip(andConfig.indexFiles, andConfig.indexFiles) + andConfig.multiIndexFiles
        multiInds = [MultiIndexCache.fromFilenameList(fnList) for fnList in indexFiles]
        return cls(multiInds)

    @classmethod
    def fromEnvironment(cls, allowCache=True):
        """Construct from the environment

        The ASTROMETRY_NET_DATA_DIR environment variable points to the
        root directory of the catalog.  If there is an "andCache.fits" file,
        and 'allowCache', we will construct from the cache.  Otherwise, we
        will construct from the index files in the "andConfig.py" file.
        """
        dirnm = os.environ.get('ASTROMETRY_NET_DATA_DIR')
        if dirnm is None:
            raise RuntimeError("astrometry_net_data is not setup")

        if allowCache:
            fn = os.path.join(dirnm, "andCache.fits")
            if os.path.exists(fn):
                return cls.fromCache(fn)

        andConfig = AstrometryNetDataConfig()
        fn = os.path.join(dirnm, 'andConfig.py')
        if not os.path.exists(fn):
            raise RuntimeError('astrometry_net_data config file \"%s\" required but not found' % fn)
        andConfig.load(fn)
        return cls.fromIndexFiles(andConfig)

    def writeCache(self, outName):
        """Write a cache file

        The cache file is a FITS file with all the required information to build the
        AstrometryNetCatalog quickly.  The first table extension contains a row for each multiindex,
        storing the healpix and nside values.  The second table extension contains a row
        for each filename in all the multiindexes.  The two may be JOINed through the
        'id' column.
        """
        numFilenames = sum(len(ind._filenameList) for ind in self._multiInds)
        maxLength = max(len(fn) for ind in self._multiInds for fn in ind._filenameList) + 1

        # First table
        first = pyfits.new_table([pyfits.Column(name="id", format="K"),
                                  pyfits.Column(name="healpix", format="K"),
                                  pyfits.Column(name="nside", format="K"),
                                  ], nrows=len(self._multiInds))
        first.data.field("id")[:] = numpy.arange(len(self._multiInds), dtype=int)
        first.data.field("healpix")[:] = numpy.array([ind._healpix for ind in self._multiInds])
        first.data.field("nside")[:] = numpy.array([ind._nside for ind in self._multiInds])

        # Second table
        second = pyfits.new_table([pyfits.Column(name="id", format="K"),
                                   pyfits.Column(name="filename", format="%dA" % (maxLength)),
                                   ], nrows=numFilenames)
        ident = second.data.field("id")
        filenames = second.data.field("filename")
        i = 0
        for j, ind in enumerate(self._multiInds):
            for fn in ind._filenameList:
                ident[i] = j
                filenames[i] = fn
                i += 1

        pyfits.HDUList([pyfits.PrimaryHDU(), first, second]).writeto(outName, clobber=True)

    @classmethod
    def fromCache(cls, filename):
        """Construct from a cache file

        Ingests the cache file written by the 'writeCache' method and
        uses that to quickly instantiate the AstrometryNetCatalog.
        """
        with pyfits.open(filename) as hduList:
            first = hduList[1].data
            second = hduList[2].data

            ident1 = first.field("id")
            ident2 = second.field("id")
            filenames = second.field("filename")

            multiInds = []
            for i, row1 in enumerate(first):
                where = numpy.where(ident2 == ident1[i]) # first JOIN second USING(id)
                multiInds.append(MultiIndexCache(filenames[where], row1.field("healpix"), row1.field("nside")))

        return cls(multiInds)

    def __getitem__(self, index):
        return self._multiInds[index]

    def __iter__(self):
        return iter(self._multiInds)

    def __len__(self):
        return len(self._multiInds)
