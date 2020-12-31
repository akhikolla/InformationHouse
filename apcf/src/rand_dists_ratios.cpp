/*===========================================================================*\
 *
 *            Nullmodels and distance calculations for APCF
 *
\*===========================================================================*/

#include <vector>        // vector
#include <Rcpp.h>        // Rcpp
#include <geos_c.h>      // GEOS functions
#include <ogr_srs_api.h> // OSR functions

#include "import_export.h"
#include "geos_stuff.h"

// [[Rcpp::export]]
Rcpp::DataFrame
rand_dists_ratios(const char* fn_pattern,
                  const char* fn_area,
                  const double max_dist,
                  const unsigned int n_simulations = 99,
                  const unsigned int max_tries = 100000,
                  const bool save_randomized_patterns = false,
                  const char* save_basename = "./pattern",
                  const bool verbose = false)
{
    GEOSContextHandle_t geosCtxtH = geos_init();

    // import -----------------------------------------------------------------
    OGRSpatialReferenceH hSRS;
    GEOSGeometry* area = import_polys(geosCtxtH, fn_area, NULL)[0];
    std::vector<GEOSGeometry*> pattern = import_polys(geosCtxtH, fn_pattern, &hSRS);

    double area_size;
    if(!GEOSArea_r(geosCtxtH, area, &area_size))
        throw std::range_error("GEOSArea failed.");
    unsigned int n_patches = pattern.size();

    // check import data at a later time
    // check parameter

    // sort pattern largest to smallest ---------------------------------------
    std::vector<double> patch_areas(n_patches);
    for(unsigned int i = 0; i < n_patches; i++)
        GEOSArea_r(geosCtxtH, pattern[i], &patch_areas[i]);
    std::vector<size_t> index = sort_indices(patch_areas, false);
    std::vector<GEOSGeometry*> pattern_sorted(n_patches);
    for(unsigned int i = 0; i < n_patches; i++)
        pattern_sorted[i] = pattern[index[i]];

    // work horse loop --------------------------------------------------------
    // set up vectors for results
    std::vector<GEOSGeometry*> temp_pattern;
    std::vector< std::vector<double> > dp_v(2), sdp_v(3);

    if(verbose)
        Rcpp::Rcout << "dists: ";

    for(unsigned int i = 0; i <= n_simulations; i++)
    {
        // exit if user interrupts function within R
        Rcpp::checkUserInterrupt();

        if(i == 0)
        {
            // emperical pcf => no randomization
            if(verbose)
                Rcpp::Rcout << 'e';

            // do calc_distances of pattern here, to avoid deep copy of pattern
            // to temp_pattern which is needed to destroy temp_pattern later
            dp_v = calc_distances(geosCtxtH, pattern, area, max_dist, false);
        }
        else
        {
            if(verbose){
                if(i % 1000 == 0)
                    Rcpp::Rcout << 'M';
                else if(i % 500 == 0)
                    Rcpp::Rcout << 'D';
                else if(i % 100 == 0)
                    Rcpp::Rcout << 'C';
                else if(i % 10 == 0)
                    Rcpp::Rcout << 'x';
            }
            // randomize patches
            temp_pattern = randomize_pattern(geosCtxtH, pattern_sorted, area, max_tries, false);

            // save randomized patterns for debugging
            if(save_randomized_patterns)
            {
                export_polys(geosCtxtH, temp_pattern,
                             (save_basename + std::to_string(i) + ".shp").c_str(),
                             NULL, "ESRI Shapefile",  hSRS);
            }

            // calc distances and ratios
            dp_v = calc_distances(geosCtxtH, temp_pattern, area, max_dist, false);

            // make sure geometries are destroyed and temp_pattern is empty
            for (unsigned int k = 0; k < n_patches; k++)
                GEOSGeom_destroy_r(geosCtxtH, temp_pattern[k]);
            temp_pattern.clear();
        }
        // append results to collection
        sdp_v[0].insert(sdp_v[0].end(), dp_v[0].size(), i);
        sdp_v[1].insert(sdp_v[1].end(), dp_v[0].begin(), dp_v[0].end());
        sdp_v[2].insert(sdp_v[2].end(), dp_v[1].begin(), dp_v[1].end());
    }
    if(verbose)
        Rcpp::Rcout << std::endl;

    // clean up ---------------------------------------------------------------
    OSRDestroySpatialReference(hSRS);
    GEOSGeom_destroy_r(geosCtxtH, area);
    for (unsigned int i = 0; i < n_patches; i++)
        GEOSGeom_destroy_r(geosCtxtH, pattern[i]);

    geos_finish(geosCtxtH);

    // create data.frame ------------------------------------------------------
    Rcpp::DataFrame sdp_df = Rcpp::DataFrame::create(
                                 Rcpp::Named("sim")  = sdp_v[0],
                                 Rcpp::Named("dist") = sdp_v[1],
                                 Rcpp::Named("ratio") = sdp_v[2]
                             );
    sdp_df.attr("area")     = area_size;
    sdp_df.attr("n_obj")    = n_patches;
    sdp_df.attr("max_dist") = max_dist;
    sdp_df.attr("class")    = Rcpp::CharacterVector::create("dists", "data.frame");

    return sdp_df;
}
