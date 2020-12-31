
#' Modify well-known geometries
#'
#' @inheritParams wk::wkb_translate_wkt
#' @param srid An integer spatial reference identifier with a user-defined meaning.
#'   Use `NA` to unset this value.
#' @param z A Z value that will be assigned to every coordinate in each feature.
#'   Use `NA` to unset this value.
#' @param trans A 3x3 transformation matrix that will be applied to all coordinates
#'   in the input.
#'
#' @return An unclassed well-known vector with the same type
#'   as the input.
#' @export
#'
#' @examples
#' wkt_set_srid("POINT (30 10)", 1234)
#' wkt_set_z("POINT (30 10)", 1234)
#' wkt_transform(
#'   "POINT (0 0)",
#'   # translation +12 +13
#'   matrix(c(1, 0, 0, 0, 1, 0, 12, 13, 1), ncol = 3)
#' )
#'
wkt_set_srid <- function(wkt, srid, precision = 16, trim = TRUE)  {
  recycled <- vctrs::vec_recycle_common(wkt, srid)
  cpp_wkt_set_srid(recycled[[1]], recycled[[2]], precision, trim)
}

#' @rdname wkt_set_srid
#' @export
wkb_set_srid <- function(wkb, srid) {
  recycled <- vctrs::vec_recycle_common(wkb, srid)
  cpp_wkb_set_srid(recycled[[1]], recycled[[2]], wk_platform_endian())
}

#' @rdname wkt_set_srid
#' @export
wksxp_set_srid <- function(wksxp, srid) {
  recycled <- vctrs::vec_recycle_common(wksxp, srid)
  cpp_wksxp_set_srid(recycled[[1]], recycled[[2]])
}

#' @rdname wkt_set_srid
#' @export
wkt_set_z <- function(wkt, z, precision = 16, trim = TRUE)  {
  recycled <- vctrs::vec_recycle_common(wkt, z)
  cpp_wkt_set_z(recycled[[1]], recycled[[2]], precision, trim)
}

#' @rdname wkt_set_srid
#' @export
wkb_set_z <- function(wkb, z) {
  recycled <- vctrs::vec_recycle_common(wkb, z)
  cpp_wkb_set_z(recycled[[1]], recycled[[2]], wk_platform_endian())
}

#' @rdname wkt_set_srid
#' @export
wksxp_set_z <- function(wksxp, z) {
  recycled <- vctrs::vec_recycle_common(wksxp, z)
  cpp_wksxp_set_z(recycled[[1]], recycled[[2]])
}

#' @rdname wkt_set_srid
#' @export
wkt_transform <- function(wkt, trans, precision = 16, trim = TRUE)  {
  cpp_wkt_transform(wkt, as_trans_matrix(trans)[c(1, 2), ], precision, trim)
}

#' @rdname wkt_set_srid
#' @export
wkb_transform <- function(wkb, trans)  {
  cpp_wkb_transform(wkb, as_trans_matrix(trans)[c(1, 2), ], endian = wk_platform_endian())
}

#' @rdname wkt_set_srid
#' @export
wksxp_transform <- function(wksxp, trans)  {
  cpp_wksxp_transform(wksxp, as_trans_matrix(trans)[c(1, 2), ])
}

as_trans_matrix <- function(trans) {
  trans <- as.matrix(trans)
  stopifnot(ncol(trans) == 3, nrow(trans) == 3)
  trans
}
