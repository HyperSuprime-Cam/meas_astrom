from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.spatial import cKDTree

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.geom as afwgeom
import lsst.afw.table as afwTable
from lsst.meas.algorithms.sourceSelector import sourceSelectorRegistry

from .pessimistic_pattern_matcher_b_3D import PessimisticPatternMatcherB

__all__ = ["MatchPessimisticBTask", "MatchPessimisticBConfig"]


class MatchTolerance(object):
    """ Stores match tolerances for use in AstrometryTask and later
    iterations of the matcher.

    Attributes
    ----------
    maxMatchDist : lsst.afw.geom.Angle
    autoMaxMatchDist : lsst.afw.geom.Angle
    maxShift : lsst.afw.geom.Angle
    lastMatchedPattern : int
    failedPatternList : list of ints
    """

    def __init__(self, maxMatchDist=None, autoMaxMatchDist=None,
                 maxShift=None, lastMatchedPattern=None,
                 failedPatternList=None):
        """ MatchPessimisticBTask relies on several state variables to be
        preserved over different iterations in the matchAndFitWcs loop of
        AstrometryTask.
        maxMatchDist is the current 2 sigma scatter from the
            previous matched wcs (if it exists. It is None if this is the first
            iteration.)
        autoMatxMatchDist is the the result of the automated match tolerance
            generation. maxShift is either None for the first iteration or is
            the magnitude of the previous iteration's wcs shift.
        lastMatchedPattern is an int reference to the position in the
            magnitude sorted source array where a successful pattern match was
            found.
        failedPatternList is a list of ints specifying indicies in the magnitude
            sourced source array to skip. These are skipped are pervious
            iterations that are likely false positives due to the code having
            to soften after a pattern is matched.
        """
        self.maxMatchDist = maxMatchDist
        self.autoMaxMatchDist = autoMaxMatchDist
        self.maxShift = maxShift
        self.lastMatchedPattern = lastMatchedPattern
        self.failedPatternList = failedPatternList


class MatchPessimisticBConfig(pexConfig.Config):
    """Configuration for MatchPessimisticBTask
    """
    numBrightStars = pexConfig.RangeField(
        doc="Number of bright stars to use. Sets the max number of patterns "
            "that can be tested.",
        dtype=int,
        default=200,
        min=2,
    )
    minMatchedPairs = pexConfig.RangeField(
        doc="Minimum number of matched pairs; see also minFracMatchedPairs.",
        dtype=int,
        default=30,
        min=2,
    )
    minFracMatchedPairs = pexConfig.RangeField(
        doc="Minimum number of matched pairs as a fraction of the smaller of "
            "the number of reference stars or the number of good sources; "
            "the actual minimum is the smaller of this value or "
            "minMatchedPairs.",
        dtype=float,
        default=0.3,
        min=0,
        max=1,
    )
    matcherIterations = pexConfig.RangeField(
        doc="Number of softening interations in matcher.",
        dtype=int,
        default=5,
        min=1,
    )
    maxOffsetPix = pexConfig.RangeField(
        doc="Maximum allowed shift of WCS, due to matching (pixel). "
            "MatchPessimisticBTask converts this to arseconds internally "
            "using the wcs pixel scale.",
        dtype=int,
        default=300,
        max=4000,
    )
    maxRotationDeg = pexConfig.RangeField(
        doc="Rotation angle allowed between sources and position reference "
            "objects (degrees).",
        dtype=float,
        default=1.0,
        max=6.0,
    )
    numPointsForShape = pexConfig.Field(
        doc="Number of points to define a shape for matching.",
        dtype=int,
        default=6,
    )
    numPointsForShapeAttempt = pexConfig.Field(
        doc="Number of points to try for creating a shape. This value should "
            "be greater than numPointsForShape.",
        dtype=int,
        default=7,
    )
    numPatternConsensus = pexConfig.Field(
        doc="Number of implied shift/rotations from patterns that must agree "
            "before it a given shift/rotation is accepted. This is only used "
            "after the first softening iteration fails and if both the "
            "number of reference and source objects is greater than "
            "numBrightStars.",
        dtype=int,
        default=3,
    )
    sourceSelector = sourceSelectorRegistry.makeField(
        doc="How to select sources for cross-matching. The default "
            "matcherSourceSelector removes objects with low S/N, bad "
            "saturated objects, edge objects, and interpolated objects.",
        default="matcher"
    )

    def setDefaults(self):
        sourceSelector = self.sourceSelector["matcher"]
        sourceSelector.setDefaults()

    # The paremeters below are not used in this matcher and are kept around
    # for compatability.
    allowedNonperpDeg = pexConfig.RangeField(
        doc="Allowed non-perpendicularity of x and y (degree). This field is "
            "not used in MatchPessimisticBTask and is kept for compatability "
            "purposes only.",
        dtype=float,
        default=3.0,
        max=45.0,
    )
    maxMatchDistArcSec = pexConfig.RangeField(
        doc="Maximum separation between reference objects and sources "
            "beyond which they will not be considered a match (arcsec). "
            "This field is not used in MatchPessimisticBTask. maxMatchDist "
            "is instead set automaticly based on the density of reference "
            "and source catlogs.",
        dtype=float,
        default=3,
        min=0,
    )
    maxDeterminant = pexConfig.Field(
        doc="maximum determinant of linear transformation matrix for a "
            "usable solution. This field is not used in MatchPessimisticB "
            "and is kept for compatability purposes only.",
        dtype=float,
        default=0.02,
    )

# The following block adds links to this task from the Task Documentation page.
# \addtogroup LSST_task_documentation
# \{
# \page measAstrom_MatchPessimisticBTask
# \ref MatchPessimisticBTask "MatchPessimisticBTask"
# Match sources to reference objects
# \}


class MatchPessimisticBTask(pipeBase.Task):
    """!Match sources to reference objects

    @anchor MatchPessimisticBTask_

    @section meas_astrom_MatchPessimisticB_Contents Contents

     - @ref meas_astrom_MatchPessimisticB_Purpose
     - @ref meas_astrom_MatchPessimisticB_Initialize
     - @ref meas_astrom_MatchPessimisticB_IO
     - @ref meas_astrom_MatchPessimisticB_Config
     - @ref meas_astrom_MatchPessimisticB_Example
     - @ref meas_astrom_MatchPessimisticB_Debug

    @section meas_astrom_MatchPessimisticB_Purpose  Description

    Match sources to reference objects. This is often done as a preliminary
    step to fitting an astrometric or photometric solution. For details about
    the matching algorithm see pessimistic_pattern_matcher_b_3D.py

    @section meas_astrom_MatchPessimisticB_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_MatchPessimisticB_IO       Invoking the Task

    @copydoc matchObjectsToSources

    @section meas_astrom_MatchPessimisticB_Config       Configuration
    parameters

    See @ref MatchPessimisticBConfig

    To modify the tests for good sources for matching, create a new
    sourceSelector class in meas_algorithms and use it in the config.

    @section meas_astrom_MatchPessimisticB_Example  A complete example of
    using MatchPessimisticBTask

    MatchPessimisticBTask is a subtask of AstrometryTask, which is called by
    PhotoCalTask. See \ref meas_photocal_photocal_Example.

    @section meas_astrom_MatchPessimisticB_Debug        Debug variables

    The @link lsst.pipe.base.cmdLineTask.CmdLineTask command line task@endlink
    interface supports a flag @c -d to import @b debug.py from your
    @c PYTHONPATH; see @ref baseDebug for more about @b debug.py files.

    The available variables in MatchPessimisticBTask are:
    <DL>
      <DT> @c verbose (bool)
      <DD> If True then the matcher prints debug messages to stdout
    </DL>

    To investigate the @ref meas_astrom_MatchPessimisticB_Debug, put something
    like
    @code{.py}
        import lsstDebug
        def DebugInfo(name):
            # N.b. lsstDebug.Info(name) would call us recursively
            debug = lsstDebug.getInfo(name)
            if name == "lsst.pipe.tasks.astrometry":
                debug.verbose = True

            return debug

        lsstDebug.Info = DebugInfo
    @endcode
    into your debug.py file and run this task with the @c --debug flag.
    """

    ConfigClass = MatchPessimisticBConfig
    _DefaultName = "matchObjectsToSources"

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("sourceSelector")

    def filterStars(self, refCat):
        """Extra filtering pass; subclass if desired
        """
        return refCat

    @pipeBase.timeMethod
    def matchObjectsToSources(self, refCat, sourceCat, wcs, refFluxField,
                              match_tolerance=None):
        """!Match sources to position reference stars

        @param[in] refCat  catalog of reference objects that overlap the
        exposure; reads fields for:
            - coord
            - the specified flux field
        @param[in] sourceCat  catalog of sources found on an exposure; reads
        fields for:
            - centroid
            - centroid flag
            - edge flag
            - saturated flag
            - aperture flux, if found, else PSF flux
        @param[in] wcs  estimated WCS
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] match_tolerance is a MatchTolerance class object or None.
            This this class is used to comunicate state between AstrometryTask
            and MatcherTask. AstrometryTask will also set the MatchTolerance
            class variable maxMatchDist based on the scatter AstrometryTask has
            found after fitting for the wcs.
        @return an lsst.pipe.base.Struct with fields:
        - matches  a list of matches, each instance of
            lsst.afw.table.ReferenceMatch
        - usableSourcCat  a catalog of sources potentially usable for
            matching.
        - match_tolerance a MatchTolerance object containing the resulting
            state variables from the match.
        """
        import lsstDebug
        debug = lsstDebug.Info(__name__)

        preNumObj = len(refCat)
        refCat = self.filterStars(refCat)
        numRefObj = len(refCat)

        if self.log:
            self.log.info(
                "filterStars purged %d reference stars, leaving %d stars" %
                (preNumObj - numRefObj, numRefObj))

        # If we get an empty tolerance struct create the variables we need for
        # this matcher.
        if match_tolerance is None:
            match_tolerance = MatchTolerance()

        # usableSourceCat: sources that are good but may be saturated
        numSources = len(sourceCat)
        selectedSources = self.sourceSelector.selectSources(sourceCat)
        goodSourceCat = selectedSources.sourceCat
        numUsableSources = len(goodSourceCat)
        self.log.info("Purged %d sources, leaving %d good sources" %
                      (numSources - numUsableSources, numUsableSources))

        if len(goodSourceCat) == 0:
            raise pipeBase.TaskError("No sources are good")

        # avoid accidentally using sourceCat; use goodSourceCat from now on
        del sourceCat

        minMatchedPairs = min(self.config.minMatchedPairs,
                              int(self.config.minFracMatchedPairs *
                                  min([len(refCat), len(goodSourceCat)])))

        # match usable (possibly saturated) sources and then purge saturated
        # sources from the match list
        doMatchReturn = self._doMatch(
            refCat=refCat,
            sourceCat=goodSourceCat,
            wcs=wcs,
            refFluxField=refFluxField,
            numUsableSources=numUsableSources,
            minMatchedPairs=minMatchedPairs,
            match_tolerance=match_tolerance,
            sourceFluxField=self.sourceSelector.fluxField,
            verbose=debug.verbose,
        )
        matches = doMatchReturn.matches
        match_tolerance = doMatchReturn.match_tolerance

        if len(matches) == 0:
            raise RuntimeError("Unable to match sources")

        self.log.info("Matched %d sources" % len(matches))
        if len(matches) < minMatchedPairs:
            self.log.warn("Number of matches is smaller than request")

        return pipeBase.Struct(
            matches=matches,
            usableSourceCat=goodSourceCat,
            match_tolerance=match_tolerance,
        )

    @pipeBase.timeMethod
    def _doMatch(self, refCat, sourceCat, wcs, refFluxField, numUsableSources,
                 minMatchedPairs, match_tolerance, sourceFluxField, verbose):
        """!Implementation of matching sources to position reference stars

        Unlike matchObjectsToSources, this method does not check if the sources
        are suitable.

        @param[in] refCat  catalog of position reference stars that overlap an
            exposure
        @param[in] sourceCat  catalog of sources found on the exposure
        @param[in] wcs  estimated WCS of exposure
        @param[in] refFluxField  field of refCat to use for flux
        @param[in] numUsableSources  number of usable sources (sources with
            known centroid that are not near the edge, but may be saturated)
        @param[in] minMatchedPairs  minimum number of matches
        @param[in] match_tolerance a MatchTolerance object containing
            variables specifying matcher tolerances and state from possible
            previous runs.
        @param[in] sourceInfo  SourceInfo for the sourceCat
        @param[in] verbose  true to print diagnostic information to std::cout

        @return a list of matches, an instance of
            lsst.afw.table.ReferenceMatch, a MatchTolerance object
        """

        # Set configurable defaults when we encounter None type or set
        # state based on pervious run of matchandfit.
        if match_tolerance.maxShift is None:
            maxShiftArcseconds = (self.config.maxOffsetPix *
                                  wcs.pixelScale().asArcseconds())
        else:
            maxShiftArcseconds = np.max(
                (match_tolerance.maxShift.asArcseconds(),
                 wcs.pixelScale().asArcseconds()))

        # Load the source and reference catalog as spherical points
        # in numpy array. We do this rather than relying on internal
        # lsst C objects for simplicity and because we require
        # objects contiguous in memory.
        ref_array = np.empty((len(refCat), 4))
        for ref_idx, refObj in enumerate(refCat):
            theta = np.pi/2 - refObj.getDec().asRadians()
            phi = refObj.getRa().asRadians()
            ref_array[ref_idx, 0] = np.sin(theta)*np.cos(phi)
            ref_array[ref_idx, 1] = np.sin(theta)*np.sin(phi)
            ref_array[ref_idx, 2] = np.cos(theta)
            ref_array[ref_idx, 3] = -2.5*np.log10(
                np.where(refObj[refFluxField] > 0,
                         refObj[refFluxField], 10**-32))

        src_array = np.empty((len(sourceCat), 4))
        for src_idx, srcObj in enumerate(sourceCat):
            coord = wcs.pixelToSky(srcObj.getCentroid())
            tmp_ra = coord.getLongitude().asRadians()
            tmp_dec = coord.getLatitude().asRadians()
            theta = np.pi/2 - tmp_dec
            phi = tmp_ra
            src_array[src_idx, 0] = np.sin(theta)*np.cos(phi)
            src_array[src_idx, 1] = np.sin(theta)*np.sin(phi)
            src_array[src_idx, 2] = np.cos(theta)
            src_array[src_idx, 3] = -2.5*np.log10(
                np.where(srcObj.getPsfFlux() > 0,
                         srcObj.getPsfFlux(), 10**-32))

        # If our tolerances are not set from a previous run, estiamte a
        # starting tolerance guess from statistics of patterns we can create
        # on both the source and reference catalog. We use the smaller of the
        # two.
        if match_tolerance.maxMatchDist is None:
            self.log.debug("Computing source statistics...")
            maxMatchDistArcSecSrc = self._get_pair_pattern_statistics(
                src_array)
            self.log.debug("Computing reference statistics...")
            maxMatchDistArcSecRef = self._get_pair_pattern_statistics(
                ref_array)
            maxMatchDistArcSec = np.min((maxMatchDistArcSecSrc,
                                         maxMatchDistArcSecRef))
            match_tolerance.autoMaxDist = afwgeom.Angle(maxMatchDistArcSec,
                                                        afwgeom.arcseconds)
        else:
            maxMatchDistArcSec = np.max(
                (wcs.pixelScale().asArcseconds(),
                 np.min((match_tolerance.maxMatchDist.asArcseconds(),
                         match_tolerance.autoMaxDist.asArcseconds()))))

        # Make sure the data we are considering is dense enough to require
        # the consensus mode of the matcher.
        numConsensus = self.config.numPatternConsensus
        if ref_array.shape[0] < self.config.numBrightStars or \
           src_array.shape[0] < self.config.numBrightStars:
            numConsensus = 1

        self.log.debug("Current tol maxDist: %.4f..." % maxMatchDistArcSec)
        self.log.debug("Current shift: %.4f" % maxShiftArcseconds)
        self.log.debug("Current rot: %.4f" % self.config.maxRotationDeg)

        # Create our matcher object.
        pyPPMb = PessimisticPatternMatcherB(ref_array[:, :3], self.log)

        # Start the ineration over our tolerances.
        for try_idx in xrange(self.config.matcherIterations):
            if try_idx == 0:
                # If we are on the first most stringent tolerance,
                # the matcher should behave like an optimistic pattern
                # matcher. Exiting at the first match.
                matcher_struct = pyPPMb.match(
                    source_catalog=src_array,
                    n_check=self.config.numPointsForShapeAttempt,
                    n_match=self.config.numPointsForShape,
                    n_agree=1,
                    max_n_patterns=self.config.numBrightStars,
                    max_shift=maxShiftArcseconds,
                    max_rotation=self.config.maxRotationDeg,
                    max_dist=maxMatchDistArcSec,
                    min_matches=minMatchedPairs,
                    pattern_skip_array=np.array(
                        match_tolerance.failedPatternList)
                )
            else:
                # Once we fail and start softening tolerences we switch to
                # pessimistic or consenus mode where 3 patterns must agree
                # on a rotation before exiting. We double the match dist
                # tolerance each round and add 1 to the pattern complexiy and
                # two to the number of candidate spokes to check.
                matcher_struct = pyPPMb.match(
                    source_catalog=src_array,
                    n_check=self.config.numPointsForShapeAttempt + 2 * try_idx,
                    n_match=self.config.numPointsForShape + try_idx,
                    n_agree=numConsensus,
                    max_n_patterns=self.config.numBrightStars,
                    max_shift=(self.config.maxOffsetPix *
                               wcs.pixelScale().asArcseconds()),
                    max_rotation=self.config.maxRotationDeg,
                    max_dist=maxMatchDistArcSec * 2. ** try_idx,
                    min_matches=minMatchedPairs,
                    pattern_skip_array=np.array(
                        match_tolerance.failedPatternList)
                )
            if matcher_struct is None and \
               match_tolerance.lastPattern is None and \
               try_idx == 0:
                # If we found a pattern on a previous run and can't
                # find an optimistic match with the harhest tolerances
                # the match we found previously was likely bad. We append
                # its index to a list so we do not use it in future match
                # iterations.
                match_tolerance.failedPatternList.append(
                    match_tolerance.lastPattern)
            elif matcher_struct is not None:
                # Match found, save a bit a state regarding this pattern
                # in the struct and exit.
                match_tolerance.maxShift = afwgeom.Angle(matcher_struct.shift,
                                                         afwgeom.arcseconds)
                match_tolerance.lastPattern = matcher_struct.pattern_idx
                break

        # A match has been found, return our list of matches and
        # return.
        matches = []
        if matcher_struct is not None:
            for match_ids, dist in zip(matcher_struct.matches,
                                       matcher_struct.distances):
                match = afwTable.ReferenceMatch()
                match.first = refCat[match_ids[1]]
                match.second = sourceCat[match_ids[0]]
                match.distance = match.first.getCoord().angularSeparation(
                    match.second.getCoord())
                matches.append(match)

        return pipeBase.Struct(
            matches=matches,
            match_tolerance=match_tolerance,
        )

    def _get_pair_pattern_statistics(self, cat_array):
        """ Convience functon for computing the tolerances for the
        matcher automatically by sorting the data as we would
        for matching and creating the patterns to match in order.
        The tolerance is then the average difference in length
        between the spokes of the pattern.
        """

        self.log.debug("Starting automated tolerance calculation...")

        # Create an empty array of all the patterns we possibly make
        # sorting from brightest to faintest.
        pattern_array = np.empty(
            (cat_array.shape[0] - self.config.numPointsForShape,
             self.config.numPointsForShape - 1))
        flux_args_array = np.argsort(cat_array[:, -1])

        # Sort our input array.
        tmp_sort_array = cat_array[flux_args_array]

        # Start making patterns.
        for start_idx in xrange(cat_array.shape[0] -
                                self.config.numPointsForShape):
            pattern_points = tmp_sort_array[start_idx:start_idx +
                                            self.config.numPointsForShape, :-1]
            pattern_delta = pattern_points[1:, :] - pattern_points[0, :]
            pattern_array[start_idx, :] = np.sqrt(
                pattern_delta[:, 0] ** 2 +
                pattern_delta[:, 1] ** 2 +
                pattern_delta[:, 2] ** 2)

            # When we store the length of each spoke in our pattern we
            # sort from shortest to longest so we have a defined space
            # to compare them in.
            pattern_array[start_idx, :] = pattern_array[
                start_idx, np.argsort(pattern_array[start_idx, :])]

        # Create a searchable tree object of the patterns and find
        # for any given pattern the closest pattern in the sorted
        # spoke length space.
        dist_tree = cKDTree(
            pattern_array[:, :(self.config.numPointsForShape - 1)])
        dist_nearest_array, ids = dist_tree.query(
            pattern_array[:, :(self.config.numPointsForShape - 1)], k=2)
        dist_nearest_array = dist_nearest_array[:, 1]
        dist_nearest_array.sort()

        # We use the two closest patterns to set our tolerance.
        dist_idx = 0
        dist_tol = (np.degrees(dist_nearest_array[dist_idx]) * 3600. /
                    (self.config.numPointsForShape - 1.))

        self.log.debug("Automated tolerance")
        self.log.debug("\tdistance/match tol: %.4f [arcsec]" % dist_tol)

        return dist_tol
