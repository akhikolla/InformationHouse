
#' Flatten nested geometry structures
#'
#' @inheritParams wk::wkb_translate_wkt
#' @param keep_empty If `TRUE`, a GEOMETRYCOLLECTION EMPTY is left as-is
#'   rather than collapsing to length 0.
#' @param keep_multi If `TRUE`, MULTI* geometries are not expanded to sub-features.
#' @param max_depth The maximum recursive GEOMETRYCOLLECTION depth to unnest.
#'
#' @return An unclassed vector with attribute `lengths`, which is an integer vector
#'   with the same length as the input denoting the length to which each
#'   feature was expanded.
#' @export
#'
#' @examples
#' wkt_unnest("GEOMETRYCOLLECTION (POINT (1 2), POINT (3 4))")
#' wkt_unnest("GEOMETRYCOLLECTION EMPTY")
#' wkt_unnest("GEOMETRYCOLLECTION EMPTY", keep_empty = TRUE)
#'
wkt_unnest <- function(wkt, keep_empty = FALSE, keep_multi = TRUE, max_depth = 1) {
  cpp_wkt_unnest(wkt, keep_empty, keep_multi, max_depth)
}

#' @rdname wkt_unnest
#' @export
wkb_unnest <- function(wkb, keep_empty = FALSE, keep_multi = TRUE, max_depth = 1) {
  cpp_wkb_unnest(wkb, keep_empty, keep_multi, max_depth, endian = wk_platform_endian())
}

#' @rdname wkt_unnest
#' @export
wksxp_unnest <- function(wksxp, keep_empty = FALSE, keep_multi = TRUE, max_depth = 1) {
  cpp_wksxp_unnest(wksxp, keep_empty, keep_multi, max_depth)
}
