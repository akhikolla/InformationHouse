/*===========================================================================*\
 *
 *                         GEOS related (Helper-)Functions
 *
\*===========================================================================*/

#include <algorithm>    // std::max
#include <cstdarg>      // va_list, va_start, va_end
#include <cstdio>       // vsprintf
#include <cstring>      // strlen
#include <geos_c.h>     // GEOS functions
#include <math.h>       // sin, cos, ...
#include <random>       // random number generator
#include <Rcpp.h>       // Rcpp
#include <vector>       // vector

#include "geos_stuff.h"



/* GEOS Version
 *-----------------------------------------------------------------------------
 * adapted from sf (https://github.com/edzer/sfr/blob/master/src/geos.cpp)
 */
// [[Rcpp::export]]
std::string geos_version()
{
    return GEOS_VERSION;
    //return GEOS_CAPI_VERSION;
}


/* Start / End of GEOS Session
 *-----------------------------------------------------------------------------
 * adapted from sf (https://github.com/edzer/sfr/blob/master/src/geos.cpp)
 */
static void __errorHandler(const char *fmt, ...)
{
    char buf[BUFSIZ], *p;
    va_list(ap);
    va_start(ap, fmt);
    vsprintf(buf, fmt, ap);
    va_end(ap);
    p = buf + strlen(buf) - 1;
    if(strlen(buf) > 0 && *p == '\n') *p = '\0';

    Rcpp::Function error("stop");
    error(buf);

    return;
}

static void __warningHandler(const char *fmt, ...)
{
    char buf[BUFSIZ], *p;
    va_list(ap);
    va_start(ap, fmt);
    vsprintf(buf, fmt, ap);
    va_end(ap);
    p = buf + strlen(buf) - 1;
    if(strlen(buf) > 0 && *p == '\n') *p = '\0';

    Rcpp::Function warning("warning");
    warning(buf);

    return;
}

GEOSContextHandle_t geos_init(void)
{
    return initGEOS_r((GEOSMessageHandler) __warningHandler,
                      (GEOSMessageHandler) __errorHandler);
}

void geos_finish(GEOSContextHandle_t geosCtxtH)
{
    finishGEOS_r(geosCtxtH);
}


/* Calculate Distances / Ratios
 *-----------------------------------------------------------------------------
 * calculates distances and ratios of buffer(dist) inside study area
 *   for all patches of pattern closer than max_dist
 *
 * throws exceptions
 */
std::vector< std::vector<double> >
calc_distances(GEOSContextHandle_t geosCtxtH,
               const std::vector<GEOSGeometry*> pattern,
               const GEOSGeometry* area,
               double max_dist, bool verbose)
{
    double dist, buf_complete_len, buf_inside_len, ratio;
    GEOSGeometry *buf_complete, *buf_complete_Bd, *buf_inside_Bd;
    std::vector< std::vector<double> > dp(2);

    // roundness of buffers. default in postgis is 8, in sf it's 30.
    int quad_segs = 8;

    for(unsigned int i = 0; i < pattern.size(); i++)
    {
        for(unsigned int j = 0; j < pattern.size(); j++)
        {
            // no distance to self
            if (i == j)
                continue;

            // get distance
            if(!GEOSDistance_r(geosCtxtH, pattern[i], pattern[j], &dist))
                throw std::range_error("GEOSDistance failed.");

            // buffer and intersection is costly => done only if necessary
            if(dist <= max_dist)
            {

                // ratio of buffer boundary inside area
                buf_complete = GEOSBuffer_r(geosCtxtH, pattern[i], dist, quad_segs);
                if(!buf_complete)
                    throw std::range_error("GEOSBuffer failed.");
                buf_complete_Bd = GEOSBoundary_r(geosCtxtH, buf_complete);
                GEOSGeom_destroy_r(geosCtxtH, buf_complete);
                if(!buf_complete)
                    throw std::range_error("GEOSBoundary failed.");

                // no need for intersection if buffer completely within area
                if(GEOSContains_r(geosCtxtH, area, buf_complete_Bd))
                {
                    ratio = 1.0;

                    // clean up
                    GEOSGeom_destroy_r(geosCtxtH, buf_complete_Bd);
                }
                else
                {
                    buf_inside_Bd = GEOSIntersection_r(geosCtxtH, area, buf_complete_Bd);
                    if(!buf_inside_Bd)
                        throw std::range_error("GEOSIntersection failed.");

                    if(!GEOSLength_r(geosCtxtH, buf_complete_Bd, &buf_complete_len))
                        throw std::range_error("GEOSLength failed.");

                    if(!GEOSLength_r(geosCtxtH, buf_inside_Bd, &buf_inside_len))
                        throw std::range_error("GEOSLength failed.");


                    if(buf_complete_len <= 0 || buf_inside_len <= 0)
                        throw std::range_error("Length of buffer <= 0.");

                    ratio = buf_inside_len / buf_complete_len;

                    // clean up
                    GEOSGeom_destroy_r(geosCtxtH, buf_inside_Bd);
                    GEOSGeom_destroy_r(geosCtxtH, buf_complete_Bd);
                }

                // write results to vector of vectors
                dp[0].push_back(dist);
                dp[1].push_back(ratio);

                if(verbose)
                    Rcpp::Rcout << "[" << i << "->" << j << "]"
                                << "  dist: " << dist
                                << "  ratio: " << ratio
                                << std::endl;
            }
        }
    }
    return dp;
}


/* Randomize Pattern
 *-----------------------------------------------------------------------------
 * creates a new pattern by randomizing all patches of the given pattern
 *
 * Assumes patches in pattern are sorted from largest to smallest
 *
 * caller gains ownership of returned geometry
 * throws exceptions
 */
std::vector<GEOSGeometry*>
randomize_pattern(GEOSContextHandle_t geosCtxtH,
                  const std::vector<GEOSGeometry*> pattern,
                  const GEOSGeometry* area,
                  const unsigned int max_tries, const bool verbose)
{
    // static stuff: do only once ---------------------------------------------
    // create seed only once (either deterministic or using random_device)
    // static int seed = 20170101;
    static std::random_device rd;
    static int seed = rd();

    // get min/max of area, setup random number engine and desired distributions
    static std::vector<double> bbox = get_extent(geosCtxtH, area);
    static std::mt19937 mt(seed);
    static std::uniform_real_distribution<double> runif_angle(0.0, 2.0 * M_PI); // radians
    static std::uniform_real_distribution<double> runif_x(bbox[0], bbox[2]);
    static std::uniform_real_distribution<double> runif_y(bbox[1], bbox[3]);


    // random move all patches ------------------------------------------------
    GEOSGeometry* patch;
    std::vector<GEOSGeometry*> randomized;

    for(unsigned int i = 0; i < pattern.size(); i++)
    {
        // create random numbers and move => random_move
        patch = move_poly(geosCtxtH, pattern[i],
                          runif_angle(mt), runif_x(mt), runif_y(mt), verbose);

        unsigned int n_tries = 0;
        while( !location_okay(geosCtxtH, patch, area, randomized, false) )
        {
            // destroy old patch before creating a new one
            GEOSGeom_destroy_r(geosCtxtH, patch);

            patch = move_poly(geosCtxtH, pattern[i],
                              runif_angle(mt), runif_x(mt), runif_y(mt), verbose);

            // avoid infinite loop
            n_tries++;
            if(n_tries > max_tries)
            {
                Rcpp::Rcout << "Exceeded max_tries (" <<  max_tries << "). "
                            "Try again (stochastic process). "
                            "Maybe increase max_tries." << std::endl;
                throw std::range_error("Failed to randomize (exceeded max_tries).");
            }
        }
        randomized.push_back(patch);
    }

    return randomized;
}


/* Move Polygon
 *-----------------------------------------------------------------------------
 * rotates patch by angle (radians) and moves polygon so that the new centroid
 * has coordinates new_cent_x, new_cent_y
 *
 * caller gains ownership of returned geometry
 * throws exceptions
 */
GEOSGeometry*
move_poly(GEOSContextHandle_t geosCtxtH, const GEOSGeometry* geom,
          const double angle, const double new_cent_x, const double new_cent_y,
          const bool verbose)
{
    if(verbose)
    {
        Rcpp::Rcout << "angle: " << angle
                    << "\tcentX: " << new_cent_x
                    << "\tcentY: " << new_cent_y << std::endl;
    }

    // get centroid and calculate shift ---------------------------------------
    double cent_x, cent_y, shift_x, shift_y;
    GEOSGeometry* cent = GEOSGetCentroid_r(geosCtxtH, geom);
    if(!cent)
        throw std::range_error("GEOSGetCentroid failed.");
    if(GEOSGeomGetX_r(geosCtxtH, cent, &cent_x) == -1 ||
       GEOSGeomGetY_r(geosCtxtH, cent, &cent_y) == -1)
        throw std::range_error("GEOSGeomGetX / GetY failed.");

    shift_x = new_cent_x - cent_x;
    shift_y = new_cent_y - cent_y;

    // move every coordinate of geometry --------------------------------------
    unsigned int numCoords, numDims;
    double old_x, old_y, new_x, new_y;

    const GEOSGeometry* linearRing = GEOSGetExteriorRing_r(geosCtxtH, geom);
    if(!linearRing)
        throw std::range_error("GEOSGetExteriorRing failed.");

    const GEOSCoordSequence* coordSeq = GEOSGeom_getCoordSeq_r(geosCtxtH, linearRing);
    if(!coordSeq)
        throw std::range_error("GEOSGeom_getCoordSeq failed.");

    if(!GEOSCoordSeq_getSize_r(geosCtxtH, coordSeq, &numCoords))
        throw std::range_error("GEOSCoordSeq_getSize failed.");
    if(!GEOSCoordSeq_getDimensions_r(geosCtxtH, coordSeq, &numDims))
        throw std::range_error("GEOSCoordSeq_getDimensions failed.");

    GEOSCoordSequence* newCoordSeq = GEOSCoordSeq_create_r(geosCtxtH, numCoords, numDims);
    if(!newCoordSeq)
        throw std::range_error("GEOSCoordSeq_create failed.");

    for(unsigned int i = 0; i < numCoords; i++)
    {
        // read coordinates
        if(!GEOSCoordSeq_getX_r(geosCtxtH, coordSeq, i, &old_x) ||
           !GEOSCoordSeq_getY_r(geosCtxtH, coordSeq, i, &old_y))
            throw std::range_error("GEOSCoordSeq_getX / getY failed.");

        // rotate
        new_x = cent_x
                + cos(angle) * (old_x - cent_x)
                - sin(angle) * (old_y - cent_y);
        new_y = cent_y
                + sin(angle) * (old_x - cent_x)
                + cos(angle) * (old_y - cent_y);

        // shift
        new_x = new_x + shift_x;
        new_y = new_y + shift_y;

        // write coordinates
        if(!GEOSCoordSeq_setX_r(geosCtxtH, newCoordSeq, i, new_x) ||
           !GEOSCoordSeq_setY_r(geosCtxtH, newCoordSeq, i, new_y))
            throw std::range_error("GEOSCoordSeq_setX / setY failed.");
    }

    // create polygon w/o hole ------------------------------------------------
    GEOSGeometry *shell, *outGeom;
    shell   = GEOSGeom_createLinearRing_r(geosCtxtH, newCoordSeq);
    outGeom = GEOSGeom_createPolygon_r(geosCtxtH, shell, NULL, 0);
    if(!shell || !outGeom)
        throw std::range_error("GEOSGeom_createLinearRing / Polygon failed.");

    // clean-up ---------------------------------------------------------------
    // newCoordSeq, shell, outGeom are all part of return value
    // linearRing, coordSeq are pointer to internal storage => must not be destroyed
    GEOSGeom_destroy_r(geosCtxtH, cent);

    return outGeom;
}


/* Check Location of Polygon
 *-----------------------------------------------------------------------------
 * checks if location of moved patch is okay
 * - within area
 * - no overlaps with already placed patches in pattern
 *
 * throws exceptions
 */
bool
location_okay(GEOSContextHandle_t geosCtxtH, const GEOSGeometry* patch,
              const GEOSGeometry* area, const std::vector<GEOSGeometry*> pattern,
              const bool verbose)
{
    // geometry within study area?
    switch(GEOSContains_r(geosCtxtH, area, patch))
    {
        case 0:   // false
        {
            if(verbose)
                Rcpp::Rcout << "STOP (not in area)" << std::endl;
            return false;
        }
        case 1:   // true
            break;
        case 2:   // error code
        default:
            throw std::range_error("GEOSContains failed.");
    }

    // pattern empty?
    // assumes the entire vector contains NULLs if first position is NULL
    if(pattern.empty() || pattern[0] == NULL)
    {
        if(verbose)
            Rcpp::Rcout << "OKAY (pattern empty)" << std::endl;
        return true;
    }

    // check the current patch against all patches of the pattern
    for (unsigned int i = 0; i < pattern.size(); i++)
    {
        // assumes there are no valid geometries after the first NULL
        if(pattern[i] == NULL)
            break;

        // check for crosses and within
        switch(GEOSIntersects_r(geosCtxtH, patch, pattern[i]))
        {
            case 0:  // false
                break;
            case 1:  // true
            {
                if(verbose)
                    Rcpp::Rcout << "STOP (intersects)" << std::endl;
                return false;
            }
            case 2: // error code
            default:
                throw std::range_error("GEOSIntersects failed.");
        }
    }
    if(verbose)
        Rcpp::Rcout << "OKAY (pattern empty)" << std::endl;

    return true;
}


/* Get Extent of Geometry
 *-----------------------------------------------------------------------------
 * min and max in X- and Y-direction of a geometry
 *
 * works only for single geometry!
 *
 * Returns vector containing (xmin, ymin, xmax, ymax)
 *
 * throws exception
 */
std::vector<double> get_extent(GEOSContextHandle_t geosCtxtH,
                               const GEOSGeometry* geom, const bool verbose)
{
    double x, y;
    std::vector<double> bbox(4);

    GEOSGeometry* envelope = GEOSEnvelope_r(geosCtxtH, geom);
    if(!envelope)
        throw std::range_error("GEOSEnvelope failed.");

    const GEOSGeometry* envRing = GEOSGetExteriorRing_r(geosCtxtH, envelope);
    if(!envRing)
        throw std::range_error("GEOSGetExteriorRing failed.");

    const GEOSCoordSequence* envCoordSeq = GEOSGeom_getCoordSeq_r(geosCtxtH, envRing);
    if(!envCoordSeq)
        throw std::range_error("GEOSGeom_getCoordSeq failed.");

    unsigned int nCoords;
    if(!GEOSCoordSeq_getSize_r(geosCtxtH, envCoordSeq, &nCoords))
        throw std::range_error("GEOSCoordSeq_getSize failed.");

    for(unsigned int i=0; i < nCoords; i++)
    {
        if(!GEOSCoordSeq_getX_r(geosCtxtH, envCoordSeq, i, &x) ||
           !GEOSCoordSeq_getY_r(geosCtxtH, envCoordSeq, i, &y))
            throw std::range_error("GEOSCoordSeq_getX / getY failed.");

        if (i == 0)
        {
            bbox[0] = bbox[2] = x;
            bbox[1] = bbox[3] = y;
        }
        else
        {
            bbox[0] = std::min(bbox[0], x);
            bbox[1] = std::min(bbox[1], y);
            bbox[2] = std::max(bbox[2], x);
            bbox[3] = std::max(bbox[3], y);
        }
    }

    GEOSGeom_destroy_r(geosCtxtH, envelope);

    return bbox;
}


/*
 * adapted from sfr/src/geos.cpp @7621d66 (2017-01-31)
 * https://github.com/edzer/sfr/blob/7621d66fe881c8d238564deaa6ef7dd9675ae0e5/src/geos.cpp
 */
// Rcpp::LogicalVector CPL_geos_is_valid(Rcpp::List sfc) {
//     GEOSContextHandle_t hGEOSCtxt = CPL_geos_init();
//     std::vector<GEOSGeom> gmv = geometries_from_sfc(hGEOSCtxt, sfc);
//     Rcpp::LogicalVector out(gmv.size());
//     for (int i = 0; i < out.length(); i++) {
//         out[i] = chk_(GEOSisValid_r(hGEOSCtxt, gmv[i]));
//         GEOSGeom_destroy_r(hGEOSCtxt, gmv[i]);
//     }
//     CPL_geos_finish(hGEOSCtxt);
//     return out;
// }

// Rcpp::LogicalVector CPL_geos_is_simple(Rcpp::List sfc) {
//     GEOSContextHandle_t hGEOSCtxt = CPL_geos_init();
//     Rcpp::LogicalVector out(sfc.length());
//     std::vector<GEOSGeom> g = geometries_from_sfc(hGEOSCtxt, sfc);
//     for (size_t i = 0; i < g.size(); i++) {
//         out[i] = chk_(GEOSisSimple_r(hGEOSCtxt, g[i]));
//         GEOSGeom_destroy_r(hGEOSCtxt, g[i]);
//     }
//     CPL_geos_finish(hGEOSCtxt);
//     return out;
// }
