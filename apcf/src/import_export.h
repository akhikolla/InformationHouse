/*===========================================================================*\
 *                   Reading and writing files via OGR
 *
\*===========================================================================*/

#ifndef IMPORT_EXPORT_INCLUDED
#define IMPORT_EXPORT_INCLUDED

#include <vector>        // vector
#include <gdal.h>        // GDAL &  OGR functions (includes ogr_api.h)
#include <ogr_srs_api.h> // OSR functions
#include <geos_c.h>      // GEOS functions


/* import polygons via OGR
 *
 * datasource: OGR data source (eg. path to Shapefile)
 * hSRS: pointer to SpatialRef. Can be NULL -> No SpatialRef returned
 *   ownership of hSRS passed on to caller.
 *   needs to be destroyed, eg. OSRDestroySpatialReference(hSRS);
 *
 * Returns a vector of GEOSGeometries
 *
 * throws exceptions
 */
std::vector<GEOSGeometry*>
import_polys(GEOSContextHandle_t geosCtxtH,
             const char *datasource,
             OGRSpatialReferenceH *hSRS);


/* export polygons via OGR
 *
 * dsn: OGR data sink (eg. path to Shapefile)
 * layer: layer name (usually just filename w/o extension. basename of dsn if NULL)
 * driver: output file format name (default is ESRI Shapefile)
 * hSRS: a SpatialRef object. Can be NULL -> No SpatialReference written
 *
 * throws exceptions
 */
void
export_polys(GEOSContextHandle_t geosCtxtH,
             const std::vector<GEOSGeometry*> vGeom,
             const char *dsn,
             const char *layer = NULL,
             const char *driver = "ESRI Shapefile",
             const OGRSpatialReferenceH hSRS = NULL);



/* import from R's SimpleFeature (package:sf)
 *
 * sfc: a simple feature collection. most generic simple feature type
 *
 * Returns a vector of GEOSGeometries
 *
 * throws exceptions ???
 */
// std::vector<GEOSGeom>
// geometries_from_sfc(GEOSContextHandle_t geosCtxtH,
//                     Rcpp::List sfc);



/* export to R's SimpleFeature (package:sf)
 *
 * geom: vector of GEOSGeometries
 *
 * Returns a Simple Feature Collection. disguised as a plain List
 *
 * throws exceptions ???
 */
// Rcpp::List
// sfc_from_geometry(GEOSContextHandle_t geosCtxtH,
//                   std::vector<GEOSGeom> geom);

#endif // IMPORT_EXPORT_INCLUDED
