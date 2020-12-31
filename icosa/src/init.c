#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void _locateTriangle_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _icosa_Aggregate_(SEXP, SEXP, SEXP);
extern SEXP _icosa_AllNeighboursTri_(SEXP, SEXP);
extern SEXP _icosa_AllShapeTri_(SEXP, SEXP);
extern SEXP _icosa_allTriangleCenters_(SEXP, SEXP, SEXP);
extern SEXP _icosa_ArcDist_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_ArcDistMany_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_ArcDistMat_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_centroidPoints_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_Collapse_(SEXP);
extern SEXP _icosa_CreateHexaSubfaces_(SEXP, SEXP, SEXP);
extern SEXP _icosa_DirectNeighboursTri_(SEXP);
extern SEXP _icosa_edgeListFromNeighbours_(SEXP);
extern SEXP _icosa_edgeMatTri_(SEXP, SEXP);
extern SEXP _icosa_edges_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_EdgesFromPoints_(SEXP, SEXP, SEXP);
extern SEXP _icosa_EdgesToFaces_(SEXP);
extern SEXP _icosa_EvenInterpolation_(SEXP, SEXP, SEXP);
extern SEXP _icosa_ExpandBoundariesToCols_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_expandEdges_(SEXP, SEXP, SEXP);
extern SEXP _icosa_ExpandEdgesByFacesTri_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_expandFacesToEdges_(SEXP);
extern SEXP _icosa_GetPatch_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_GreatCircle_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_HexaFaces_(SEXP);
extern SEXP _icosa_IcosahedronTesselation_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_InverseWeightByFaceCenter_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_NeighbourOfOneFace_(SEXP, SEXP);
extern SEXP _icosa_OccupiedCellUpSampling_(SEXP, SEXP);
extern SEXP _icosa_orderTriGrid_(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_Partition_(SEXP, SEXP, SEXP);
extern SEXP _icosa_pointLayerColorOrder_(SEXP);
extern SEXP _icosa_projectCloseToPoint_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_Refine2d_(SEXP, SEXP);
extern SEXP _icosa_RetrieveIndexMat_(SEXP);
extern SEXP _icosa_ShapeTri_(SEXP, SEXP, SEXP);
extern SEXP _icosa_SizeEstimate_(SEXP);
extern SEXP _icosa_SphericalTriangleCenter_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_SphericalTriangleSurface_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_spherTriSurfs(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_surfConvHullTri(SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_SplitArc_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_stl_sort(SEXP);
extern SEXP _icosa_SymmetricArcDistMat_(SEXP, SEXP, SEXP);
extern SEXP _icosa_TriangleTesselation_(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icosa_triMatTri_(SEXP, SEXP);
extern SEXP _icosa_whichMinVector_(SEXP);
extern SEXP _icosa_xxxxyyyyzzzz_(SEXP, SEXP);
extern SEXP _icosa_xyz1(SEXP);
extern SEXP _icosa_xyz1xyz1xyz1xyz1_(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"_locateTriangle_", (DL_FUNC) &_locateTriangle_, 11},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"_icosa_Aggregate_",                 (DL_FUNC) &_icosa_Aggregate_,                 3},
  {"_icosa_AllNeighboursTri_",          (DL_FUNC) &_icosa_AllNeighboursTri_,          2},
  {"_icosa_AllShapeTri_",               (DL_FUNC) &_icosa_AllShapeTri_,               2},
  {"_icosa_allTriangleCenters_",        (DL_FUNC) &_icosa_allTriangleCenters_,        3},
  {"_icosa_ArcDist_",                   (DL_FUNC) &_icosa_ArcDist_,                   4},
  {"_icosa_ArcDistMany_",               (DL_FUNC) &_icosa_ArcDistMany_,               4},
  {"_icosa_ArcDistMat_",                (DL_FUNC) &_icosa_ArcDistMat_,                4},
  {"_icosa_centroidPoints_",            (DL_FUNC) &_icosa_centroidPoints_,            4},
  {"_icosa_Collapse_",                  (DL_FUNC) &_icosa_Collapse_,                  1},
  {"_icosa_CreateHexaSubfaces_",        (DL_FUNC) &_icosa_CreateHexaSubfaces_,        3},
  {"_icosa_DirectNeighboursTri_",       (DL_FUNC) &_icosa_DirectNeighboursTri_,       1},
  {"_icosa_edgeListFromNeighbours_",    (DL_FUNC) &_icosa_edgeListFromNeighbours_,    1},
  {"_icosa_edgeMatTri_",                (DL_FUNC) &_icosa_edgeMatTri_,                2},
  {"_icosa_edges_",                     (DL_FUNC) &_icosa_edges_,                     4},
  {"_icosa_EdgesFromPoints_",           (DL_FUNC) &_icosa_EdgesFromPoints_,           3},
  {"_icosa_EdgesToFaces_",              (DL_FUNC) &_icosa_EdgesToFaces_,              1},
  {"_icosa_EvenInterpolation_",         (DL_FUNC) &_icosa_EvenInterpolation_,         3},
  {"_icosa_ExpandBoundariesToCols_",    (DL_FUNC) &_icosa_ExpandBoundariesToCols_,    5},
  {"_icosa_expandEdges_",               (DL_FUNC) &_icosa_expandEdges_,               3},
  {"_icosa_ExpandEdgesByFacesTri_",     (DL_FUNC) &_icosa_ExpandEdgesByFacesTri_,     4},
  {"_icosa_expandFacesToEdges_",        (DL_FUNC) &_icosa_expandFacesToEdges_,        1},
  {"_icosa_GetPatch_",                  (DL_FUNC) &_icosa_GetPatch_,                  4},
  {"_icosa_GreatCircle_",               (DL_FUNC) &_icosa_GreatCircle_,               5},
  {"_icosa_HexaFaces_",                 (DL_FUNC) &_icosa_HexaFaces_,                 1},
  {"_icosa_IcosahedronTesselation_",    (DL_FUNC) &_icosa_IcosahedronTesselation_,    4},
  {"_icosa_InverseWeightByFaceCenter_", (DL_FUNC) &_icosa_InverseWeightByFaceCenter_, 7},
  {"_icosa_NeighbourOfOneFace_",        (DL_FUNC) &_icosa_NeighbourOfOneFace_,        2},
  {"_icosa_OccupiedCellUpSampling_",    (DL_FUNC) &_icosa_OccupiedCellUpSampling_,    2},
  {"_icosa_orderTriGrid_",              (DL_FUNC) &_icosa_orderTriGrid_,              6},
  {"_icosa_Partition_",                 (DL_FUNC) &_icosa_Partition_,                 3},
  {"_icosa_pointLayerColorOrder_",      (DL_FUNC) &_icosa_pointLayerColorOrder_,      1},
  {"_icosa_projectCloseToPoint_",       (DL_FUNC) &_icosa_projectCloseToPoint_,       4},
  {"_icosa_Refine2d_",                  (DL_FUNC) &_icosa_Refine2d_,                  2},
  {"_icosa_RetrieveIndexMat_",          (DL_FUNC) &_icosa_RetrieveIndexMat_,          1},
  {"_icosa_ShapeTri_",                  (DL_FUNC) &_icosa_ShapeTri_,                  3},
  {"_icosa_SizeEstimate_",              (DL_FUNC) &_icosa_SizeEstimate_,              1},
  {"_icosa_SphericalTriangleCenter_",   (DL_FUNC) &_icosa_SphericalTriangleCenter_,   4},
  {"_icosa_SphericalTriangleSurface_",  (DL_FUNC) &_icosa_SphericalTriangleSurface_,  5},
  {"_icosa_spherTriSurfs",              (DL_FUNC) &_icosa_spherTriSurfs,              4},
  {"_icosa_surfConvHullTri",            (DL_FUNC) &_icosa_surfConvHullTri,            4},
  {"_icosa_SplitArc_",                  (DL_FUNC) &_icosa_SplitArc_,                  5},
  {"_icosa_stl_sort",                   (DL_FUNC) &_icosa_stl_sort,                   1},
  {"_icosa_SymmetricArcDistMat_",       (DL_FUNC) &_icosa_SymmetricArcDistMat_,       3},
  {"_icosa_TriangleTesselation_",       (DL_FUNC) &_icosa_TriangleTesselation_,       5},
  {"_icosa_triMatTri_",                 (DL_FUNC) &_icosa_triMatTri_,                 2},
  {"_icosa_whichMinVector_",            (DL_FUNC) &_icosa_whichMinVector_,            1},
  {"_icosa_xxxxyyyyzzzz_",              (DL_FUNC) &_icosa_xxxxyyyyzzzz_,              2},
  {"_icosa_xyz1",                       (DL_FUNC) &_icosa_xyz1,                       1},
  {"_icosa_xyz1xyz1xyz1xyz1_",          (DL_FUNC) &_icosa_xyz1xyz1xyz1xyz1_,          2},
  {NULL, NULL, 0}
};

void R_init_icosa(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
