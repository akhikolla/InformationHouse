/*===========================================================================*\
 *
 *                      GEOS related (Helper-)Functions
 *
\*===========================================================================*/
#ifndef GEOS_STUFF_INCLUDED
#define GEOS_STUFF_INCLUDED

#include <vector>     // vector
#include <algorithm>  // sort, swap
#include <numeric>    // iota
#include <geos_c.h>   // GEOS functions



/* Start / End of GEOS Session
 *-----------------------------------------------------------------------------
 */
GEOSContextHandle_t geos_init(void);
void geos_finish(GEOSContextHandle_t geosCtxtH);


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
               double max_dist,
               bool verbose = false);


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
                  const unsigned int max_tries,
                  const bool verbose = false);


/* Move Polygon
 *-----------------------------------------------------------------------------
 * rotates patch by angle (radians) and moves polygon so that the new centroid
 * has coordinates new_cent_x, new_cent_y
 *
 * caller gains ownership of returned geometry
 * throws exceptions
 */
GEOSGeometry*
move_poly(GEOSContextHandle_t geosCtxtH,
          const GEOSGeometry* geom,
          const double angle, const double new_cent_x, const double new_cent_y,
          const bool verbose = false);


/* Check Location of Polygon
 *-----------------------------------------------------------------------------
 * checks if location of moved patch is okay
 * - within area
 * - no overlaps with already placed patches in pattern
 *
 * throws exceptions
 */
bool
location_okay(GEOSContextHandle_t geosCtxtH,
              const GEOSGeometry* patch,
              const GEOSGeometry* area,
              const std::vector<GEOSGeometry*> pattern,
              const bool verbose = false);


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
std::vector<double>
get_extent(GEOSContextHandle_t geosCtxtH,
           const GEOSGeometry* geom,
           const bool verbose = false);


/* Template Function for Sorting Indices
 *-----------------------------------------------------------------------------
 *
 * must live entirely in header file
 */

template <typename T>
std::vector<size_t> sort_indices(const std::vector<T> &vec, const bool asc=true)
{
    // initialize original index locations
    std::vector<size_t> idx(vec.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in vec
    if(asc){
        std::sort(idx.begin(), idx.end(),
                  [&vec](size_t a, size_t b)
                  {
                      return vec[a] < vec[b];
                  });
    }
    else
    {
        std::sort(idx.begin(), idx.end(),
                  [&vec](size_t a, size_t b)
                  {
                      return vec[a] > vec[b];
                  });
    }

    return idx;
}


/* Template Function for in Place Permutation
 *-----------------------------------------------------------------------------
 * http://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 * must live entirely in header file
 */
template <typename T>
void permutate_using_index(std::vector<T>& vec,
                           const std::vector<std::size_t>& idx)
{
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        if (done[i])
        {
            continue;
        }
        done[i] = true;
        std::size_t prev_j = i;
        std::size_t j = idx[i];
        while (i != j)
        {
            std::swap(vec[prev_j], vec[j]);
            done[j] = true;
            prev_j = j;
            j = idx[j];
        }
    }
}

#endif // GEOS_STUFF_INCLUDED
