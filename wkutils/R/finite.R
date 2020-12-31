
#' Test well-known geometries for missing and non-finite coordinates
#'
#' Note that EMTPY geometries are considered finite and non-missing.
#' Use the `size` column of [wkt_meta()] to test for empty geometries.
#'
#' @inheritParams wk::wkb_translate_wkt
#'
#' @return A logical vector with the same length as the input.
#' @export
#'
#' @examples
#' wkt_has_missing("POINT (0 1)")
#' wkt_has_missing("POINT (nan nan)")
#' wkt_has_missing("POINT (inf inf)")
#'
#' wkt_is_finite("POINT (0 1)")
#' wkt_is_finite("POINT (nan nan)")
#' wkt_is_finite("POINT (inf inf)")
#'
wkt_has_missing <- function(wkt) {
  cpp_wkt_has_missing(wkt)
}

#' @rdname wkt_has_missing
#' @export
wkb_has_missing <- function(wkb) {
  cpp_wkb_has_missing(wkb)
}

#' @rdname wkt_has_missing
#' @export
wksxp_has_missing <- function(wksxp) {
  cpp_wksxp_has_missing(wksxp)
}

#' @rdname wkt_has_missing
#' @export
wkt_is_finite <- function(wkt) {
  !cpp_wkt_has_non_finite(wkt)
}

#' @rdname wkt_has_missing
#' @export
wkb_is_finite <- function(wkb) {
  !cpp_wkb_has_non_finite(wkb)
}

#' @rdname wkt_has_missing
#' @export
wksxp_is_finite <- function(wksxp) {
  !cpp_wksxp_has_non_finite(wksxp)
}
