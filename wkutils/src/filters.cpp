
#include "wk/wkt-writer.hpp"
#include "wk/wkt-reader.hpp"
#include "wk/wkb-writer.hpp"
#include "wk/wkb-reader.hpp"
#include "wk/rcpp-sexp-writer.hpp"
#include "wk/rcpp-sexp-reader.hpp"
#include "wk/filter.hpp"

#include <Rcpp.h>
#include "wk/rcpp-io.hpp"
using namespace Rcpp;


// --------- srid -------------

class  WKSetSridFilter: public WKMetaFilter {
public:
  WKSetSridFilter(WKGeometryHandler& handler, IntegerVector srid):
    WKMetaFilter(handler), srid(srid), featureSrid(NA_REAL) {}

  void nextFeatureStart(size_t featureId) {
    this->featureSrid = this->srid[featureId];
    WKMetaFilter::nextFeatureStart(featureId);
  }

  WKGeometryMeta newGeometryMeta(const WKGeometryMeta& meta, uint32_t partId) {
    WKGeometryMeta newMeta(meta);
    if (IntegerVector::is_na(this->featureSrid)) {
      newMeta.hasSRID = false;
    } else {
      newMeta.hasSRID = true;
      newMeta.srid = this->featureSrid;
    }

    return newMeta;
  }

private:
  IntegerVector srid;
  int featureSrid;
};


void set_srid_base(WKReader& reader, WKWriter& writer, IntegerVector srid) {
  WKSetSridFilter filter(writer, srid);
  reader.setHandler(&filter);

  while (reader.hasNextFeature()) {
    checkUserInterrupt();
    reader.iterateFeature();
  }
}

// [[Rcpp::export]]
CharacterVector cpp_wkt_set_srid(CharacterVector wkt, IntegerVector srid,
                                 int precision = 16, bool trim = true) {
  WKCharacterVectorProvider provider(wkt);
  WKTReader reader(provider);

  WKCharacterVectorExporter exporter(wkt.size());
  WKTWriter writer(exporter);
  exporter.setRoundingPrecision(precision);
  exporter.setTrim(trim);
  set_srid_base(reader, writer, srid);
  return exporter.output;
}

// [[Rcpp::export]]
List cpp_wkb_set_srid(List wkb, IntegerVector srid, int endian) {
  WKRawVectorListProvider provider(wkb);
  WKBReader reader(provider);

  WKRawVectorListExporter exporter(wkb.size());
  WKBWriter writer(exporter);
  writer.setEndian(endian);
  set_srid_base(reader, writer, srid);
  return exporter.output;
}

// [[Rcpp::export]]
List cpp_wksxp_set_srid(List wksxp, IntegerVector srid) {
  WKRcppSEXPProvider provider(wksxp);
  WKRcppSEXPReader reader(provider);

  WKRcppSEXPExporter exporter(wksxp.size());
  WKRcppSEXPWriter writer(exporter);
  set_srid_base(reader, writer, srid);
  return exporter.output;
}

// ----------- set z -------------


class WKSetZFilter: public WKMetaFilter {
public:
  WKSetZFilter(WKGeometryHandler& handler, NumericVector z):
  WKMetaFilter(handler), z(z), featureZ(NA_REAL) {}

  void nextFeatureStart(size_t featureId) {
    this->featureZ = this->z[featureId];
    WKMetaFilter::nextFeatureStart(featureId);
  }

  WKGeometryMeta newGeometryMeta(const WKGeometryMeta& meta, uint32_t partId) {
    WKGeometryMeta newMeta(meta);
    newMeta.hasZ = !NumericVector::is_na(this->featureZ);
    return newMeta;
  }

  void nextCoordinate(const WKGeometryMeta& meta, const WKCoord& coord, uint32_t coordId) {
    WKCoord newCoord(coord);
    newCoord.z = this->featureZ;
    newCoord.hasZ = !NumericVector::is_na(this->featureZ);
    WKMetaFilter::nextCoordinate(meta, newCoord, coordId);
  }

private:
  NumericVector z;
  double featureZ;
};


void set_z_base(WKReader& reader, WKWriter& writer, NumericVector z) {
  WKSetZFilter filter(writer, z);
  reader.setHandler(&filter);

  while (reader.hasNextFeature()) {
    checkUserInterrupt();
    reader.iterateFeature();
  }
}

// [[Rcpp::export]]
CharacterVector cpp_wkt_set_z(CharacterVector wkt, NumericVector z,
                              int precision = 16, bool trim = true) {
  WKCharacterVectorProvider provider(wkt);
  WKTReader reader(provider);

  WKCharacterVectorExporter exporter(wkt.size());
  WKTWriter writer(exporter);
  exporter.setRoundingPrecision(precision);
  exporter.setTrim(trim);
  set_z_base(reader, writer, z);
  return exporter.output;
}

// [[Rcpp::export]]
List cpp_wkb_set_z(List wkb, NumericVector z, int endian) {
  WKRawVectorListProvider provider(wkb);
  WKBReader reader(provider);

  WKRawVectorListExporter exporter(wkb.size());
  WKBWriter writer(exporter);
  writer.setEndian(endian);
  set_z_base(reader, writer, z);
  return exporter.output;
}

// [[Rcpp::export]]
List cpp_wksxp_set_z(List wksxp, NumericVector z) {
  WKRcppSEXPProvider provider(wksxp);
  WKRcppSEXPReader reader(provider);

  WKRcppSEXPExporter exporter(wksxp.size());
  WKRcppSEXPWriter writer(exporter);
  set_z_base(reader, writer, z);
  return exporter.output;
}

// ---------- transform -----------

class WKTransformFilter: public WKFilter {
public:
  // here, t is the 6-member affine transform in column-major format
  // (e.g., [1 0 0 1 tx ty])
  WKTransformFilter(WKGeometryHandler& handler, NumericVector t): WKFilter(handler),
    t1(t[0]), t2(t[1]), t3(t[2]), t4(t[3]), t5(t[4]), t6(t[5]) {}

  void nextCoordinate(const WKGeometryMeta& meta, const WKCoord& coord, uint32_t coordId) {
    WKCoord newCoord(coord);
    newCoord.x = t1 * coord.x + t3 * coord.y + t5;
    newCoord.y = t2 * coord.x + t4 * coord.y + t6;
    WKFilter::nextCoordinate(meta, newCoord, coordId);
  }

private:
  double t1, t2, t3, t4, t5, t6;
};

void transform_base(WKReader& reader, WKWriter& writer, NumericVector transform) {
  WKTransformFilter filter(writer, transform);
  reader.setHandler(&filter);

  while (reader.hasNextFeature()) {
    checkUserInterrupt();
    reader.iterateFeature();
  }
}

// [[Rcpp::export]]
CharacterVector cpp_wkt_transform(CharacterVector wkt, NumericVector transform,
                                  int precision = 16, bool trim = true) {
  WKCharacterVectorProvider provider(wkt);
  WKTReader reader(provider);

  WKCharacterVectorExporter exporter(wkt.size());
  WKTWriter writer(exporter);
  exporter.setRoundingPrecision(precision);
  exporter.setTrim(trim);
  transform_base(reader, writer, transform);
  return exporter.output;
}

// [[Rcpp::export]]
List cpp_wkb_transform(List wkb, NumericVector transform, int endian) {
  WKRawVectorListProvider provider(wkb);
  WKBReader reader(provider);

  WKRawVectorListExporter exporter(wkb.size());
  WKBWriter writer(exporter);
  writer.setEndian(endian);
 transform_base(reader, writer, transform);
  return exporter.output;
}

// [[Rcpp::export]]
List cpp_wksxp_transform(List wksxp, NumericVector transform) {
  WKRcppSEXPProvider provider(wksxp);
  WKRcppSEXPReader reader(provider);

  WKRcppSEXPExporter exporter(wksxp.size());
  WKRcppSEXPWriter writer(exporter);
  transform_base(reader, writer, transform);
  return exporter.output;
}
