
#' Generate grid geometries from well-known geometries
#'
#' Using [wkt_meta()] and [wkt_coords()], these functions create graphical objects
#' using the grid package. Vectors that contain geometries of a single dimension
#' are efficiently packed into a [grid::pointsGrob()], [grid::polylineGrob()],
#' or [grid::pathGrob()]. Vectors with mixed types and nested collections are encoded
#' less efficiently using a [grid::gTree()].
#'
#' @inheritParams wk::wkb_translate_wkt
#' @param ... Graphical parameters passed to [grid::gpar()]. These are recycled along
#'   the input. Dynamic dots (e.g., `!!!`) are supported.
#' @param rule Use "winding" if polygon rings are correctly encoded with a winding
#'   direction.
#' @param default.units Coordinate units, which may be defined by the viewport (see
#'   [grid::unit()]). Defaults to native.
#' @param name,vp Passed to [grid::pointsGrob()], [grid::polylineGrob()],
#'   [grid::pathGrob()], or [grid::gTree()] depending on the types of geometries
#'   in the input.
#'
#' @return A [graphical object][grid::grob]
#' @export
#'
#' @examples
#' grid::grid.newpage()
#' grid::grid.draw(wkt_grob("POINT (0.5 0.5)", pch = 16, default.units = "npc"))
#'
wkt_grob <- function(wkt, ..., rule = "evenodd", default.units = "native", name = NULL, vp = NULL) {
  grob_wk_possibly_nested(
    wkt,
    ...,
    unnest_fun = wkt_unnest,
    meta_fun = wkt_meta,
    coords_fun = wkt_coords,
    rule = rule,
    default.units = default.units,
    name = name,
    vp = vp
  )
}

#' @rdname wkt_grob
#' @export
wkb_grob <- function(wkt, ..., rule = "evenodd", default.units = "native", name = NULL, vp = NULL) {
  grob_wk_possibly_nested(
    wkt,
    ...,
    unnest_fun = wkb_unnest,
    meta_fun = wkb_meta,
    coords_fun = wkb_coords,
    rule = rule,
    default.units = default.units,
    name = name,
    vp = vp
  )
}

#' @rdname wkt_grob
#' @export
wksxp_grob <- function(wkt, ..., rule = "evenodd", default.units = "native", name = NULL, vp = NULL) {
  grob_wk_possibly_nested(
    wkt,
    ...,
    unnest_fun = wksxp_unnest,
    meta_fun = wksxp_meta,
    coords_fun = wksxp_coords,
    rule = rule,
    default.units = default.units,
    name = name,
    vp = vp
  )
}

grob_wk_possibly_nested <- function(x, ..., unnest_fun, meta_fun, coords_fun,
                                    rule, default.units, name, vp) {
  meta <- meta_fun(x)
  gpar_values_all <- vctrs::vec_recycle_common(..., .size = length(x))

  # if there are any collections, unnest everything
  if (any((meta$size > 0) & (meta$type_id == 7))) {
    unnested <- unnest_fun(x, keep_empty = FALSE, keep_multi = TRUE, max_depth = 10)

    lengths <- attr(unnested, "lengths")
    run_length_enc <- structure(
      list(lengths = lengths, values = seq_along(lengths)),
      class = "rle"
    )
    unnested_gpar <- lapply(gpar_values_all, "[", inverse.rle(run_length_enc))

    return(
      wkt_grob(
        unnested, !!!unnested_gpar,
        rule = rule,
        default.units = default.units,
        name = name,
        vp = vp
      )
    )
  }

  coords <- coords_fun(x)

  # grid doesn't do zero-length input, so if there are no coordinates, return
  # an empty grob
  if (nrow(coords) == 0) {
    return(grid::gTree(name = name, vp = vp, children = grid::gList()))
  }

  grob_wk_base(
    meta, coords, gpar_values_all,
    rule = rule,
    default.units = default.units,
    name = name,
    vp = vp
  )
}

grob_wk_base <- function(meta, coords, gpar_values_all, rule, default.units, name = NULL, vp = NULL) {
  # use meta to try to create the most efficient grob possible
  non_empty_types <- meta$type_id[meta$size > 0]

  # non-empty IDs are used by mixed types and polygons
  non_empty_features <- unique(meta$feature_id[meta$size > 0])

  if (all(non_empty_types %in% c(1, 4))) {
    # Need to get a tiny bit creative here, because pointsGrob
    # doesn't have a pathId argument like pathGrob. The key here is
    # to select gpar_values_all so that its rows match the
    # number of actual part_id values for each feature
    part_lengths <- rle(coords$part_id)
    part_row_start <- c(0, cumsum(part_lengths$lengths)) + 1
    part_feature_ids <- coords$feature_id[part_row_start[-length(part_row_start)]]
    gpar_values_all <- lapply(gpar_values_all, "[", part_feature_ids)

    pch <- gpar_values_all$pch %||% 1
    size <- gpar_values_all$size %||% grid::unit(1, "char")
    gpar_values_all$pch <- NULL
    gpar_values_all$size <- NULL

    grid::pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = pch,
      size = size,
      gp = do.call(grid::gpar, gpar_values_all),
      default.units = default.units,
      name = name,
      vp = vp
    )
  } else if (all(non_empty_types %in% c(2, 5))) {
    # Need to get a tiny bit creative here, because polylineGrob
    # doesn't have a pathId argument like pathGrob. The key here is
    # to select gpar_values_all so that its rows match the
    # number of actual part_id values for each feature
    part_lengths <- rle(coords$part_id)
    part_row_start <- c(0, cumsum(part_lengths$lengths)) + 1
    part_feature_ids <- coords$feature_id[part_row_start[-length(part_row_start)]]
    gpar_values_all <- lapply(gpar_values_all, "[", part_feature_ids)

    grid::polylineGrob(
      x = coords$x,
      y = coords$y,
      id.lengths = part_lengths$lengths,
      gp = do.call(grid::gpar, gpar_values_all),
      default.units = default.units,
      name = name,
      vp = vp
    )
  } else if (all(non_empty_types %in% c(3, 6))) {
    non_empty_gpar_values <- lapply(gpar_values_all, "[", non_empty_features)
    grid::pathGrob(
      x = coords$x,
      y = coords$y,
      id = coords$ring_id,
      pathId = coords$feature_id,
      gp = do.call(grid::gpar, non_empty_gpar_values),
      default.units = default.units,
      name = name,
      vp = vp,
      rule = rule
    )
  } else if (all(non_empty_types != 7)) {
    # Mixed input, but no collections (collections should be handled by unnest())
    # not very efficient, but better than failing
    grobs <- lapply(non_empty_features, function(feature_id) {
      grob_wk_base(
        meta[feature_id, ],
        coords[coords$feature_id == feature_id, ],
        gpar_values_all,
        rule = rule,
        default.units = default.units
        # name and vp are passed to the gTree
      )
    })

    grid::gTree(children = do.call(grid::gList, grobs), name = name, vp = vp)
  } else {
    stop("Can't use grob_wk_base() on a GEOMETRYCOLLECTION") # nocov
  }
}
