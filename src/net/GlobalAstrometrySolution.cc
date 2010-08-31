// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 

#include "lsst/afw/image/TanWcs.h"
#include "lsst/meas/astrom/net/GlobalAstrometrySolution.h"
#include <cstdio>
#include <iostream>

namespace lsst { 
namespace meas { 
namespace astrom { 
namespace net {


using namespace std;
namespace afwImg = lsst::afw::image;
namespace afwCoord = lsst::afw::coord;
namespace Except = lsst::pex::exceptions;
namespace Det = lsst::afw::detection;
namespace pexLog = lsst::pex::logging;


int const USE_ALL_STARS_FOR_SOLUTION = -1;

//Refuse to try to solve star lists with fewer than this many objects. Doing so
//increases the chances of returning a false match. The number twenty is
//recommended by astrometry.net as a minimum value.
static int defaultMinimumNumberOfObjectsToAccept = 20;


static vector<double> getTagAlongFromIndex(index_t* index, string fieldName, int *ids, int numIds);

//
//Constructors, Destructors
//
    GlobalAstrometrySolution::GlobalAstrometrySolution(const std::string policyPath, pexLog::Log mylog):
    _mylog(mylog),
    _indexList(NULL), 
    _solver(NULL), 
    _starxy(NULL), 
    _numBrightObjects(USE_ALL_STARS_FOR_SOLUTION),
    _minimumNumberOfObjectsToAccept(defaultMinimumNumberOfObjectsToAccept),
    _isSolved(false) {
    
    _solver = solver_new();

    setDefaultValues();
    
    lsst::pex::policy::Policy pol(policyPath);
    _equinox = pol.getDouble("equinox");
    _raDecSys = pol.getString("raDecSys");
    string pkgDir = lsst::utils::eups::productDir("ASTROMETRY_NET_DATA");

    //Add meta information about every index listed in the policy file
    std::vector<std::string> indexArray = pol.getStringArray("indexFile");

    _mylog.log(pexLog::Log::DEBUG, "Loading Astrometry.net index files (AKA astrometry_net_data)...");
    for (unsigned int i = 0; i < indexArray.size(); ++i) {
        index_t *meta = _loadIndexMeta(pkgDir + "/" + indexArray[i]);
        if (!meta)
            continue;
        // Check for duplicates.
        bool duplicate = false;
        for (unsigned int j=0; j<_indexList.size(); j++) {
            index_t* other = _indexList[j];
            //These three values [are supposed to] uniquely identify an index
            if (meta->indexid == other->indexid &&
                meta->healpix == other->healpix &&
                meta->hpnside == other->hpnside) {
                _mylog.format(pexLog::Log::WARN, "Index file \"%s\" is a duplicate (has same index id, healpix and healpix nside) as index file \"%s\"",
                              meta->indexname, other->indexname);
                duplicate = true;
                break;
            }
        }
        if (duplicate)
            continue;
        _indexList.push_back(meta);
    }
    _mylog.format(pexLog::Log::DEBUG, "Loaded %i Astrometry.net index files", _indexList.size());
}

void GlobalAstrometrySolution::loadIndices() {
    for (unsigned int i=0; i<_indexList.size(); i++)
        index_reload(_indexList[i]);
}

std::vector<const index_t*> GlobalAstrometrySolution::getIndexList() {
    std::vector<const index_t*> rtn;
    for (unsigned int i=0; i<_indexList.size(); i++)
        rtn.push_back(_indexList[i]);
    return rtn;
}


index_t *GlobalAstrometrySolution::_loadIndexMeta(std::string filename) {
    return index_load(filename.c_str(), INDEX_ONLY_LOAD_METADATA, NULL);
}


GlobalAstrometrySolution::~GlobalAstrometrySolution() {

    for (unsigned int i=0; i<_indexList.size(); i++) {
        index_free(_indexList[i]);
    }
    _indexList.clear();

    if ( _starxy != NULL) {
        starxy_free(_starxy);
        _starxy = NULL;
    }
    
    if (_solver != NULL){
        solver_free(_solver);
        _solver = NULL;
    }
}


///astrometry.net intialises the solver with some default values that guarantee failure in any
///attempted match. These values are more reasonable. 
void GlobalAstrometrySolution::setDefaultValues() {

    //Among other things, this sets the parity, positional uncertainty (_solver->verify_pix)                 <
    //matching accuracy (_solver->codetol                    
    solver_set_default_values(_solver);
    
    //Set image scale boundaries (in arcseconds per pixel) to non-zero and non-infinity.
    //These values still exceed anything you'll find in a real image
    setMinimumImageScale(1e-6);
    setMaximumImageScale(3600*360);  //2pi radians per pixel

    //How good must a match be to be considered good enough?  Log-odds.
    setMatchThreshold(log(1e12));

    // Reset counters and record of best match found so far.
    solver_cleanup_field(_solver);

    setParity(UNKNOWN_PARITY);
}

    
//
// Setup functions
//

void GlobalAstrometrySolution::setImageSize(int W, int H) {
	 solver_set_field_bounds(_solver, 0, W, 0, H);
}

///Set the image to be solved. The image is abstracted as a list of positions in pixel space
void GlobalAstrometrySolution::setStarlist(lsst::afw::detection::SourceSet vec ///<List of Sources
                                          ) {
    if (vec.empty()) {
        throw(LSST_EXCEPT(pexExcept::LengthErrorException, "Src list contains no objects"));
    }

    if (vec.size() < (unsigned int) _minimumNumberOfObjectsToAccept) {
        string msg = boost::str(boost::format("Source list should contain at least %i objects\n") % \
                 _minimumNumberOfObjectsToAccept);
        throw(LSST_EXCEPT(Except::LengthErrorException, msg));
    }

    // Step 1. Copy every valid element of the input vector into a tempory starlist structure
    // Valid means all of x, y and psfFlux are positive and finite.
    int const size = vec.size();
    starxy_t *tmpStarxy = starxy_new(size, true, false);   

    int i = 0;
    for (lsst::afw::detection::SourceSet::iterator ptr = vec.begin(); ptr != vec.end(); ++ptr) {
        double const x    = (*ptr)->getXAstrom();
        double const y    = (*ptr)->getYAstrom();
        double const flux = (*ptr)->getPsfFlux();

        //Only include objects where positions and fluxes are positive finite values
        if( isfinite(x)    && (x>=0) &&
            isfinite(y)    && (y>=0) &&
            isfinite(flux) && (flux>0)
          ) {
            starxy_set(tmpStarxy, i, x, y);
            starxy_set_flux(tmpStarxy, i, flux);
            ++i;
        }
    }

    if (i < (int) _minimumNumberOfObjectsToAccept) {
        string msg = boost::str(boost::format( \
            "Source list only has %i valid objects, needs %i\n") % i % _minimumNumberOfObjectsToAccept);
        msg += "Valid objects have positive, finite values for x, y and psfFlux";
        throw(LSST_EXCEPT(pexExcept::LengthErrorException, msg));
    }


    //Step 2. Copy these elements into a new starxy structure of the correct size. If we don't 
    //do this, starxy will report its size incorrectly, and we'll get confused later on.
    if (_starxy != NULL) {
        starxy_free(_starxy);
    }
    
    if(i == size) {
        //Every part of tmpStarxy is valid
        _starxy = tmpStarxy;
    } else {
        _starxy = starxy_new(i, true, false);
        
        for(int j=0; j<i; ++j) {
            starxy_set_x( _starxy, j, starxy_get_x( tmpStarxy, j));
            starxy_set_y( _starxy, j, starxy_get_y( tmpStarxy, j));
            starxy_set_flux( _starxy, j, starxy_get_flux( tmpStarxy, j));
        }
        
        starxy_free(tmpStarxy);
    }

    //Sort the array
    starxy_sort_by_flux(_starxy);
    _solverSetField();
}


///\brief Only find a solution using the brightest N objects
///Reducing the number of sources in the solution list can reduce the time taken to solve the image
///The truncated list is used by solve() to get the linear wcs, but all input sources are used when
///calculating the distortion terms of the SIP matrix.    
void GlobalAstrometrySolution::setNumBrightObjects(int N) {

    if (N <= 0) {
        throw(LSST_EXCEPT(pexExcept::RangeErrorException, "Illegal request. N must be greater than zero"));
    }

    _numBrightObjects = N;
    
    if (_starxy != NULL){
       _solverSetField();
    }
}


void GlobalAstrometrySolution::_solverSetField() {

    if ( ! _starxy) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    int starxySize = starxy_n(_starxy);
    if ( starxySize == 0){
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist has zero elements"));
    }

    int N = _numBrightObjects;  //Because I'm a lazy typist
    //The default value, -1, indicates that all objects should be used
    if (N == USE_ALL_STARS_FOR_SOLUTION) {
        N = starxySize;
    }

    starxy_t *shortlist = starxy_subset(_starxy, N);
    assert(shortlist);

    //Set the pointer in the solver to the new, smaller field
    starxy_free(solver_get_field(_solver));
    solver_free_field(_solver);
    solver_set_field(_solver, shortlist);
    solver_reset_field_size(_solver);

    //Find field boundaries and precompute kdtree
    solver_preprocess_field(_solver);
}   




///Set the plate scale of the image in arcsec per pixel
void GlobalAstrometrySolution::setImageScaleArcsecPerPixel(double imgScale ///< Plate scale of image
                                                          ) {
    //Note that the solver will fail if min==max, so we make them different by a small amount.    
    setMinimumImageScale(0.99*imgScale);
    setMaximumImageScale(1.01*imgScale);
}


///Set the verbosity level for astrometry.net. The higher the level the more information is returned.
///1 and 2 are typically good values to use. 4 will print so much to the screen that it slows execution
void GlobalAstrometrySolution::setLogLevel(int level) {
    if (level < 0 || level > 4) {
        throw(LSST_EXCEPT(pexExcept::DomainErrorException, "Logging level must be between 0 and 4"));
    }
    
    log_init((enum log_level) level);
}


///    
///How good does a match need to be to be accepted. Typical value is log(1e12) approximately 27
void GlobalAstrometrySolution::setMatchThreshold(double threshold) {
    solver_set_record_logodds(_solver, threshold);
}


///You can double the speed of a match if you know the parity, i.e whether the image is flipped or not.
///North up and East right (or some rotation thereof) is parity==NORMAL_PARITY, the opposite is
///parity==FLIPPED_PARITY. The default is UNKNOWN_PARITY
void GlobalAstrometrySolution::setParity(int parity){
    //Insist on legal values.
    if (solver_set_parity(_solver, parity)) {
        throw LSST_EXCEPT(pexExcept::DomainErrorException, "Illegal parity setting");
    }
}

//
// Solve functions
//

///Solve using a wcs object as an initial guess
bool GlobalAstrometrySolution::solve(const lsst::afw::image::Wcs::Ptr wcsPtr, 
                                     double imageScaleUncertaintyPercent) {

    //Rename the variable to something shorter to make the code easier to read
    double unc = imageScaleUncertaintyPercent/100.0;

    //This test is strictly unecessary as solverSetField throws the same exception
    if (! _starxy) {
        throw LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist hasn't been set yet");
    }

    double xc, yc;
    solver_get_field_center(_solver, &xc, &yc);

    //Get the central ra/dec and plate scale
    afwCoord::Coord::ConstPtr raDec = wcsPtr->pixelToSky(xc, yc);
    double const ra = raDec->getLongitude(afwCoord::DEGREES);
    double const dec = raDec->getLatitude(afwCoord::DEGREES);
    double const plateScaleArcsecPerPixel = sqrt(wcsPtr->pixArea(afwGeom::makePointD(xc, yc))) * 3600;
    setMinimumImageScale(plateScaleArcsecPerPixel*(1 - unc));
    setMaximumImageScale(plateScaleArcsecPerPixel*(1 + unc));
    
    _mylog.log(pexLog::Log::DEBUG,
               boost::format("Solving using initial guess at position of (%.7f %.7f)") % ra % dec);

    double lwr = plateScaleArcsecPerPixel*(1 - unc);
    double upr = plateScaleArcsecPerPixel*(1 + unc);
                     
    _mylog.log(pexLog::Log::DEBUG, boost::format("Exposure's WCS scale: %g arcsec/pix; setting scale range %.3f - %.3f arcsec/pixel") % plateScaleArcsecPerPixel % lwr % upr);

    if ( wcsPtr->isFlipped()) {
        setParity(FLIPPED_PARITY);
        _mylog.log(pexLog::Log::DEBUG, "Setting Flipped parity");        
    } else {
        setParity(NORMAL_PARITY);
        _mylog.log(pexLog::Log::DEBUG, "Setting Normal parity");        
    }
        
    return solve(ra, dec);
}


///Find a solution with an initial guess at the position    
bool GlobalAstrometrySolution::solve(const afw::image::PointD raDec   ///<Right ascension/declination
                                               ///in decimal degrees
                                          ) {
    return solve(raDec[0], raDec[1]);
}


    
///Find a solution with an initial guess at the position.
bool GlobalAstrometrySolution::solve(double ra,   ///<Right ascension in decimal degrees
                                     double dec   ///< Declination in decimal degrees
                                          )  {    
    string msg;

    // Tell the solver to only consider matches within the image size of the supposed RA,Dec.
    // The magic number 2.0 out front says to accept matches within 2 radii of the given *center* position.
    double maxRadius = 2.0 * arcsec2deg(solver_get_max_radius_arcsec(_solver));
    msg = boost::str(boost::format("Setting RA,Dec = (%g, %g), radius = %g deg") % ra % dec % maxRadius);
    _mylog.log(pexLog::Log::DEBUG, msg);
    solver_set_radec(_solver, ra, dec, maxRadius);

    _mylog.log(pexLog::Log::DEBUG, "Doing solve step");
    _callSolver(ra, dec);

    if (_isSolved){
        const char* indexname = solver_get_best_match_index_name(_solver);
        msg = boost::str(boost::format("Solved index is %s") % indexname);
    } else {
        msg = boost::str(boost::format("Failed to verify position (%.7f %.7f)\n") % ra % dec);
    }
    _mylog.log(pexLog::Log::DEBUG, msg);            

    return _isSolved;
}


///Find a solution blindly, with no initial guess. Go get a cup of tea, this function
///will take a while
bool GlobalAstrometrySolution::solve()  {    

    // Don't use any hints about the RA,Dec position.
    solver_clear_radec(_solver);

    _callSolver(NO_POSITION_SET, NO_POSITION_SET);

    if (_isSolved){
        _mylog.log(pexLog::Log::DEBUG, "Position Found");
        const char* indexname = solver_get_best_match_index_name(_solver);
        string msg = boost::str(boost::format("Solved index is %s") % indexname);
        _mylog.log(pexLog::Log::DEBUG, msg);            

        // Grab everything we need from the index file while it is still open!
        //index_t* index = _solver->best_match.index;
        // FIXME -- grab reference RA,Dec, fieldxy, tagalong data, etc now, or later.

    }
    else {
        _mylog.log(pexLog::Log::DEBUG, "Failed");
    }

    return _isSolved;
}

///Check that all the setup was done correctly, then set the solver to work
///By default, _callSolver does a blind solve, unless the optional ra and dec
///are given values. Sets _isSolved to true if a solution is found
bool GlobalAstrometrySolution::_callSolver(double ra, double dec) {
    //Throw exceptions if setup is incorrect
    if ( ! _starxy) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist hasn't been set yet"));
    }

    if (_indexList.size() == 0) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No index files loaded yet"));
    }
    
    if (_isSolved){
        string msg = "Solver indicated that a match has already been found. Do you need to reset?";
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
    }

    double lower = solver_get_pixscale_low(_solver);
    double upper = solver_get_pixscale_high(_solver);

    if (lower >= upper) {
        string msg = boost::str(boost::format("Minimum image scale (%g) must be strictly less than max scale (%g)") % lower % upper);
        throw(LSST_EXCEPT(pexExcept::DomainErrorException, msg));
    }

    //Calculate the best guess at image size
    double xSizePixels = solver_field_width(_solver);
    double ySizePixels = solver_field_height(_solver);
    double minSizePixels = min(xSizePixels, ySizePixels);
    double maxSizePixels = max(xSizePixels, ySizePixels);
    assert(maxSizePixels >= minSizePixels && minSizePixels > 0);

    //Set the range of sizes of quads to examine, in pixels.
    solver_set_quad_size_fraction(_solver, 0.1, 1.0);

    //Output some useful debugging info
    _mylog.format(pexLog::Log::DEBUG, "Image size %.0f x %.0f pixels", xSizePixels, ySizePixels);
    _mylog.format(pexLog::Log::DEBUG, "Searching plate scale range %.3f -- %.3f arcsec/pixel",
                  lower, upper);
    _mylog.format(pexLog::Log::DEBUG, "--> Image size %.3f x %.3f to %.3f x %.3f arcmin",
                  arcsec2arcmin(lower * xSizePixels), arcsec2arcmin(lower * ySizePixels), 
                  arcsec2arcmin(upper * xSizePixels), arcsec2arcmin(upper * ySizePixels));

    double qlo, qhi;
    solver_get_quad_size_range_arcsec(_solver, &qlo, &qhi);
    _mylog.format(pexLog::Log::DEBUG, "Using indices with quads in the range %.2f to %.2f arcmin",
                  arcsec2arcmin(qlo), arcsec2arcmin(qhi));

    _mylog.log(pexLog::Log::DEBUG, "Setting indices");
    _addSuitableIndicesToSolver(qlo, qhi, ra, dec);

    _mylog.log(pexLog::Log::DEBUG, "Doing solve step");

    //solver_print_to(_solver, stdout);

    solver_run(_solver);

    
    if (solver_did_solve(_solver)) {
        _isSolved = true;

        MatchObj* match = solver_get_best_match(_solver);
        _mylog.format(pexLog::Log::DEBUG, "Solved: %i matches, %i conflicts, %i unmatched; %i in index",
                      (int)match->nmatch, (int)match->nconflict, (int)match->ndistractor, (int)match->nindex);

        _mylog.log(pexLog::Log::DEBUG, "Calling tweak2() to tune up match...");
        _mylog.format(pexLog::Log::DEBUG, "Starting log-odds: %g", match->logodds);
        // Use "tweak2" to tune up this match, resulting in a better WCS and more catalog matches.
        // magic 1: only go to linear order (no SIP distortions).
        solver_tweak2(_solver, match, 1);

        _mylog.format(pexLog::Log::DEBUG, "After tweak2: %i matches, %i conflicts, %i unmatched",
                      (int)match->nmatch, (int)match->nconflict, (int)match->ndistractor);
        _mylog.format(pexLog::Log::DEBUG, "After tweak2: log-odds: %g", match->logodds);

    } else {
        _isSolved = false;
    }
    
    return _isSolved;
}



/// \brief Find indices that may contain a the correct solution, and add them to the solver.
/// 
/// Find indices that cover a suitable range of plate scales and optionally a suitable position
/// on the sky.
/// If these indices have been previously loaded from disk, add them to the solver, otherwise
/// load them from disk, then add them to the solver.
/// Because these files are typically large, caching them in memory can save a lot of 
/// initialisation time.
/// Because each index can take a long time to search, only using suitable ones speeds
/// the matching process.
/// 
/// If ra and dec are not supplied, the default values (see the header file) indicate to the 
/// function that position should not be used to determine the suitablity of an index
/// 
/// \param ra Optional right ascension of inital guess at solution position
/// \param dec Optional declination of inital guess at solution position
/// 
/// \return Number of indices loaded
int GlobalAstrometrySolution::_addSuitableIndicesToSolver(double quadSizeArcsecLwr, 
    double quadSizeArcsecUpr, double ra, double dec) {

    bool hasAtLeastOneIndexOfSuitableScale = false;
    int nMeta = _indexList.size();
    int nSuitable = 0;
    bool blind = (ra == NO_POSITION_SET) || (dec == NO_POSITION_SET);

    for (int i = 0; i<nMeta; ++i){
        index_t* index = _indexList[i];

        if (!index_overlaps_scale_range(index, quadSizeArcsecLwr, quadSizeArcsecUpr))
            continue;

        hasAtLeastOneIndexOfSuitableScale = true;

        //Is this either a blind solve, or does the index cover a suitable
        //patch of sky
        if (!(blind || index_is_within_range(index, ra, dec, arcsec2deg(quadSizeArcsecUpr))))
            continue;

        // Found a good one!

        // Load the data (not just metadata), if it hasn't been already...
        index_reload(index);

        string msg = boost::str(boost::format("Adding index %s") % index->indexname);
    	_mylog.log(pexLog::Log::DEBUG, msg);

        solver_add_index(_solver, index);
        nSuitable++;
    }
    
    if(nSuitable == 0) {
        string msg = "No suitable indices found for given input parameters:";
        
        if(hasAtLeastOneIndexOfSuitableScale) {
            msg += "Probably the ra/dec range isn't covered";
        }
        else {
            msg += "No indices of a suitable scale were found";
        }
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
    }
    
    return nSuitable;
}


//
//Return the solution
//              

///After solving, return a linear Wcs (i.e without distortion terms)
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getWcs()  {

    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    MatchObj* match = solver_get_best_match(_solver);

    lsst::afw::geom::PointD crpix = lsst::afw::geom::makePointD(match->wcstan.crpix[0],
                                                                match->wcstan.crpix[1]);   
    lsst::afw::geom::PointD crval = lsst::afw::geom::makePointD(match->wcstan.crval[0],
                                                                match->wcstan.crval[1]);
    
    int naxis = 2;   //This is hardcoded into the sip_t structure
    Eigen::Matrix2d CD;
    for (int i = 0; i<naxis; ++i) {
        for (int j = 0; j<naxis; ++j) {
            CD(i, j) = match->wcstan.cd[i][j];
        }
    }

    std::string const ctype1="RA---TAN";
    std::string const ctype2="DEC--TAN";
    return lsst::afw::image::Wcs::Ptr(new lsst::afw::image::Wcs(crval, crpix, CD,
                                                                ctype1, ctype2, _equinox, _raDecSys)); 
}


///
///After solving, return a full Wcs including SIP distortion matrics
lsst::afw::image::Wcs::Ptr GlobalAstrometrySolution::getDistortedWcs(int order)  {
    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    if (! _starxy){
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Starlist isn't set"));
    }

    //Generate an array of radec of positions in the field
    MatchObj* mo = solver_get_best_match(_solver);

    //Call the tweaking algorthim to generate the distortion coeffecients
    
    //jitter is a measure of how much we can expect the xy of stars to scatter from the expected
    //radec due to noise in our measurments.
    double jitterArcsec = tan_pixel_scale(&mo->wcstan) * solver_get_field_jitter(_solver);
    jitterArcsec = hypot(jitterArcsec, mo->index_jitter);
    int inverseOrder = order;
    int iterations = 5;        //blind.c:628 uses 5
    bool isWeighted = true;
    int skipShift = true;

    sip_t *sip = tweak_just_do_it(&mo->wcstan, _starxy, mo->refxyz,
                                  NULL, NULL, NULL, mo->nindex,
                                  jitterArcsec, order, inverseOrder,
                                  iterations, isWeighted, skipShift);

    //Check that tweaking worked.
    if (sip == NULL) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "Tweaking failed"));
    }

    lsst::afw::geom::PointD crpix = lsst::afw::geom::makePointD(sip->wcstan.crpix[0],
                                                                sip->wcstan.crpix[1]);
    lsst::afw::geom::PointD crval = lsst::afw::geom::makePointD(sip->wcstan.crval[0],
                                                                sip->wcstan.crval[1]);

    //Linear conversion matrix
    int naxis = 2;
    Eigen::Matrix2d CD;
    for (int i = 0; i<naxis; ++i) {
        for (int j = 0; j<naxis; ++j) {
            CD(i, j) = sip->wcstan.cd[i][j];
        }
    }


    //Forward distortion terms. In the SIP notation, these matrices are referred to
    //as A and B. I can find no documentation that insists that these matrices be
    //the same size, so I assume they aren't.
    int aSize = sip->a_order + 1;
    Eigen::MatrixXd sipA(aSize, aSize);
    for (int i = 0; i<aSize; ++i){
        for (int j = 0; j<aSize; ++j){
            sipA(i, j) = sip->a[i][j];
        }
    }

    //Repeat for B
    int bSize = sip->b_order + 1;
    Eigen::MatrixXd sipB(bSize, bSize);
    for (int i = 0; i<bSize; ++i){
        for (int j = 0; j<bSize; ++j){
            sipB(i, j) = sip->b[i][j];
        }
    }

    //Repeat for Ap, for the reverse transform
    int apSize = sip->ap_order + 1;
    Eigen::MatrixXd sipAp(apSize, apSize);
    for (int i = 0; i<apSize; ++i){
        for (int j = 0; j<apSize; ++j){
            sipAp(i, j) = sip->ap[i][j];
        }
    }

    //And finally, Bp, also part of the reverse transform
    int bpSize = sip->bp_order + 1;
    Eigen::MatrixXd sipBp(bpSize, bpSize);
    for (int i = 0; i<bpSize; ++i){
        for (int j = 0; j<bpSize; ++j){
            sipBp(i, j) = sip->bp[i][j];
        }
    }    
        
    lsst::afw::image::Wcs::Ptr
        wcsPtr(new lsst::afw::image::TanWcs(crval, crpix, CD, sipA, sipB, sipAp, sipBp, _equinox, _raDecSys));
    
    sip_free(sip);
    return wcsPtr;
}    



///\brief Return a list of the sources that match and the catalogue objects they matched to.
///
///For each object in the catalogue that matches an input source return a SourceMatch object
///containing a) the catalogue object, b) the source object, c) the distance between them in pixels
///For the two objects, X,YAstrom and Ra/Dec are set.
vector<Det::SourceMatch> GlobalAstrometrySolution::getMatchedSources(string filterName){
    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }
    
    string msg="";  //Used for exception messages
    MatchObj* match = solver_get_best_match(_solver);
    assert(match->nfield > 0);
    
    vector<Det::SourceMatch> sourceMatchSet;
    afwImg::Wcs::Ptr wcsPtr = this->getWcs();

    //Load magnitude information from catalogue
    vector<double> tagAlong = getTagAlongFromIndex(match->index, filterName, match->theta, match->nfield);

    const starxy_t* fieldxy = solver_get_field(_solver);
        
    for (int i=0; i<match->nfield; i++) {
        // "theta" is the mapping from image (aka field) stars to index (aka reference) stars.
        // negative means no match.
        if (match->theta[i] < 0)
            continue;
        
        //Matching input sources    
        lsst::afw::detection::Source::Ptr sPtr(new lsst::afw::detection::Source());
        double x = starxy_getx(fieldxy, i);
        double y = starxy_gety(fieldxy, i);
        double flux = starxy_get_flux(fieldxy, i);
        
        sPtr->setXAstrom(x);
        sPtr->setYAstrom(y);
        sPtr->setPsfFlux(flux);
        
        //Weight positions by confidence in the fact that they match
        double confidence = verify_logodds_to_weight(match->matchodds[i]); //Can be == 0
        sPtr->setXAstromErr( 1/(confidence + DBL_EPSILON));
        sPtr->setYAstromErr( 1/(confidence + DBL_EPSILON));
        
        afwCoord::Coord::Ptr cooPtr = wcsPtr->pixelToSky(x, y);
        sPtr->setRa(cooPtr->toFk5().getRa(afwCoord::DEGREES));
        sPtr->setDec(cooPtr->toFk5().getDec(afwCoord::DEGREES));
        
        //Matching catalogue objects
        lsst::afw::detection::Source::Ptr cPtr(new lsst::afw::detection::Source());
        //xyzarr2radecdeg() is defined in astrometry.net. It converts the position of the star
        //as a three dimensional unit vector to ra,dec.
        double ra, dec;
        xyzarr2radecdeg(match->refxyz + match->theta[i]*3, &ra, &dec);
        cPtr->setRa(ra);
        cPtr->setDec(dec);

        if( filterName != "") {
            double mag = tagAlong[i];
            cPtr->setPsfFlux( pow(10.0, -mag/2.5) );
        }
        
        
        lsst::afw::geom::PointD p = wcsPtr->skyToPixel(ra, dec);
        cPtr->setXAstrom(p[0]);
        cPtr->setYAstrom(p[1]);

        cPtr->setXAstromErr( 1/(confidence + DBL_EPSILON));
        cPtr->setYAstromErr( 1/(confidence + DBL_EPSILON));


        double dx = cPtr->getXAstrom() - sPtr->getXAstrom();
        double dy = cPtr->getYAstrom() - sPtr->getYAstrom();
        double dist = hypot(dx, dy);
        sourceMatchSet.push_back(Det::SourceMatch(cPtr, sPtr, dist));
    }


    return sourceMatchSet;
}


///Astrometry.net catalogues store additional data about objects in addition to their position (usually 
///magnitudes. This function returns the names of the strings used to describe those fields
///
///\note This function makes a potentially untrue assumption. It assumes that the meta data fields
///are uniform in all indices read from disk, so it only checks for metadata in one index file
vector<string> GlobalAstrometrySolution::getCatalogueMetadataFields() {

    vector<string> output;
    
    //Reload the index if necssary
    index_reload(_indexList[0]);
    if (! startree_has_tagalong(_indexList[0]->starkd) ) {
        _mylog.log(pexLog::Log::DEBUG, "No metadata found for index");        
        return output;
    }
    
    //Don't free this pointer, it points to memory that may be used later
    fitstable_t *table = startree_get_tagalong(_indexList[0]->starkd);
    assert(table != NULL);
    
    sl *nameList = fitstable_get_fits_column_names(table, NULL);
    
    int numNames = sl_size(nameList);
    for(int i=0; i< numNames; ++i) {
        string name(sl_pop(nameList));
        output.push_back(name);
    }
    
    sl_free2(nameList);
        
    return output;
}    

///Returns a sourceSet of objects that are nearby in an raDec sense to the best match solution
lsst::afw::detection::SourceSet 
GlobalAstrometrySolution::getCatalogueForSolvedField(string filterName, double margin) {
    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }
    MatchObj* match = solver_get_best_match(_solver);
    double* radec;
    int* starinds;
    int nstars;
    double scale;
    double r2;
    int i, W, H;
    int outi;
    Det::SourceSet out;

    // since _isSolved is True, we should have the match...
    assert(match);

    // arcsec/pix
    scale = tan_pixel_scale(&(match->wcstan));
    // add margin
    r2 = deg2distsq(match->radius_deg + arcsec2deg(scale * margin));

    // The index should still be open...
    assert(match->index);
    assert(match->index->starkd);
    
    startree_search_for(match->index->starkd, match->center, r2, NULL, &radec, &starinds, &nstars);
    if (nstars == 0)
        return out;

    W = match->wcstan.imagew;
    H = match->wcstan.imageh;
    outi = 0;
    for (i=0; i<nstars; i++) {
        double px, py;
        if (!tan_radec2pixelxy(&(match->wcstan), radec[2*i+0], radec[2*i+1], &px, &py))
            continue;
        // in bounds (+ margin) ?
        if (px < -margin || px > W+margin || py < -margin || py > H+margin)
            continue;
        // keep it
        Det::Source::Ptr src(new lsst::afw::detection::Source());
        src->setRa (radec[2*i+0]);
        src->setDec(radec[2*i+1]);
        out.push_back(src);
        // record the indices that are in-bounds so we can retrieve tag-along data below...
        starinds[outi] = starinds[i];
        outi++;
    }
    free(radec);

    vector<double> mag = getTagAlongFromIndex(match->index, filterName, starinds, outi);
    if (mag.size()) {
        for (unsigned int i=0; i<out.size(); i++) {
            // It seems crazy to convert mags back to fluxes...
            double flux = pow(10.0, -mag[i]/2.5);
            out[i]->setPsfFlux(flux);
            _mylog.format(pexLog::Log::DEBUG, "catalog obj %i: mag %.1f, flux %.3g", i, mag[i], out[i]->getPsfFlux());
        }
    } else {
        _mylog.format(pexLog::Log::DEBUG, "Tag-along field \"%s\" was not found in the Astrometry.net index (file \"%s\").  Catalog fluxes will not be set.",
                      filterName.c_str(), match->index->indexname);
    }
    free(starinds);

    return out;
}
    

///Returns a sourceSet of objects that are nearby in an raDec sense to the best match solution
lsst::afw::detection::SourceSet 
GlobalAstrometrySolution::getCatalogue(double radiusInArcsec, string filterName) {

    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    //Don't free this pointer. It points to memory that will still exist
    //when this function goes out of scope.
    MatchObj* match = solver_get_best_match(_solver);
    double *center = match->center;
    double ra, dec;
    xyzarr2radecdeg(center, &ra, &dec);
    
    return getCatalogue(ra, dec, radiusInArcsec, filterName);
}


///Returns a sourceSet of objects that are nearby in an raDec sense to the requested position. If
///filterName is not blank, we also extract out the magnitude information for that filter and 
///store (as a flux) in the returned SourceSet object. The value of filterName must match one of the
///strings returned by getCatalogueMetadataFields(). If you're not interested in fluxes, set
///filterName to ""
lsst::afw::detection::SourceSet GlobalAstrometrySolution::getCatalogue(double ra,
    double dec,
    double radiusInArcsec,
    string filterName) {

    double center[3];
    radecdeg2xyzarr(ra, dec, center);
    double radius2 = arcsec2distsq(radiusInArcsec);
    string msg;

    Det::SourceSet out;


    // FIXME -- this leads to many duplicate entries!!
    for (unsigned int i=0; i<_indexList.size(); i++) {
        index_t* index = _indexList[i];
        if (!index_is_within_range(index, ra, dec, arcsec2deg(radiusInArcsec)))
            continue;

        // Ensure the index is loaded...
        index_reload(index);

        //Find nearby stars
        double *radec = NULL;
        int *starinds = NULL;
        int nstars = 0;
        startree_search_for(index->starkd, center, radius2, NULL, &radec, &starinds, &nstars);

        if (nstars > 0) {
#if 0        
            double *mag = NULL;
            if (filterName != "") {
                // Grab tag-along data here. If it's not there, throw an exception
                if (! startree_has_tagalong(index->starkd) ) {
                    msg = boost::str(boost::format("Index file \"%s\" has no metadata") % index->indexname);
                    throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
                }
                mag = startree_get_data_column(index->starkd, filterName.c_str(), starinds, nstars);

                if (mag == NULL) {
                    msg = boost::str(boost::format("No meta data called %s found in index %s") %
                        filterName % index->indexname);
                    throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));            }
            }
#else
            vector<double> mag = getTagAlongFromIndex(index, filterName, starinds, nstars);
#endif            
            //Create a source for every position stored
            for (int j = 0; j<nstars; ++j) {
                Det::Source::Ptr ptr(new lsst::afw::detection::Source());

                ptr->setRa(radec[2*j]);
                ptr->setDec(radec[2*j + 1]);

                if(mag.size() > 0) {   //convert mag to flux
                    ptr->setPsfFlux( pow(10.0, -mag[j]/2.5) );
                }

                out.push_back(ptr);
            }
        }
        free(radec);
        free(starinds);
    }

    return out;

}


///A convenient interface to astrometry.net's startree_get_data_column. 
///Checks that the fieldName exists, and that something is returned.
///Returns an vector of length numIds
///
///\param index Astrometry.net index field to extract tagalong data from
///\param fieldName Name of tagalong column to extract. If this is empty (i.e ""), nothing is
///       extracted
///\param ids   Indices you want extracted. Get this from a match object, or startree_search_for
///\param numIds size of ids
static vector<double> getTagAlongFromIndex(index_t* index, string fieldName, int *ids, int numIds) {

    //Load magnitude information from catalogue
    double *tagAlong=NULL;
    string msg;

    if (fieldName == "")
        return vector<double>();
    
    // Grab tag-along data here. If it's not there, throw an exception
    if (!startree_has_tagalong(index->starkd)) {
        msg = boost::str(boost::format("Index file \"%s\" has no metadata") % index->indexname);
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
    }

    assert(index->starkd);
    tagAlong = startree_get_data_column(index->starkd, fieldName.c_str(), ids, numIds);

    if (tagAlong == NULL) {
        msg = boost::str(boost::format("No meta data called %s found in index %s") %
                         fieldName % index->indexname);
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, msg));
    }

    vector<double> out(tagAlong + 0, tagAlong + numIds);
    // the vector constructor makes a copy of the data, so free() it
    // http://www.sgi.com/tech/stl/Vector.html#1
    free(tagAlong);
    return out;
}


///Plate scale of solution in arcsec/pixel. Note this is different than getMin(Max)ImageScale()
///which return the intial guesses of platescale.
double GlobalAstrometrySolution::getSolvedImageScale(){
    if (! _isSolved) {
        throw(LSST_EXCEPT(pexExcept::RuntimeErrorException, "No solution found yet. Did you run solve()?"));
    }

    MatchObj* match = solver_get_best_match(_solver);
    return match->scale;
    //return tan_pixel_scale(&match->wcstan)
} 

/*
startree_t* GlobalAstrometrySolution::getSolvedStartree() {
    MatchObj* match = solver_get_best_match(_solver);
    if (!match)
        return NULL;
    if (!match->index)
        return NULL;
    index_reload(match->index);
    return match->index->starkd;
}
 */

MatchObj* GlobalAstrometrySolution::getMatchObject() {
    return solver_get_best_match(_solver);
}


///Reset the object so it's ready to match another field.
void GlobalAstrometrySolution::reset() {

    if (_solver != NULL) {
        solver_free(_solver);
        _solver = solver_new();
    }
    
    if (_starxy != NULL) {
        starxy_free(_starxy);
        _starxy = NULL;
    }

    _numBrightObjects = USE_ALL_STARS_FOR_SOLUTION;
    _isSolved = false;
        
    //I should probably be smarter than this and remember the actual values of
    //the settings instead of just resetting the defaults
    setDefaultValues();
}
                                                         
}}}}
