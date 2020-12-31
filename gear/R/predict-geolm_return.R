#' Return output of predict.geolm* in proper format
#'
#' @param kdtf A data.frame with pred, mspe, rmspe, and possibly sim
#' @param newcoords A matrix with the prediction coordinates
#' @param return_type The type of object to return
#' @param coordnames A vector of coordinate names
#' @noRd
return_predict_geolm = function(kdtf, newcoords, return_type, coordnames) {
  # return results
  if (requireNamespace("sp", quietly = TRUE) &&
      return_type == "SpatialPointsDataFrame") {
    return_SpatialPointsDataFrame(kdtf, newcoords, coordnames)
  } else if (!requireNamespace("sp", quietly = TRUE) &&
             return_type == "SpatialPointsDataFrame") {
    warning("sp package must be installed to enable this functionality. Returning a geardf object.")
    geardf(kdtf, newcoords)
  } else if (requireNamespace("sf", quietly = TRUE) &&
             return_type == "sf") {
    return_sf(kdtf, newcoords, coordnames)
  } else if (!requireNamespace("sf", quietly = TRUE) &&
             return_type == "sf") {
    warning("sf package must be installed to enable this functionality. Returning a geardf object.")
    geardf(kdtf, newcoords)
  } else if (return_type == "gearPredict") {
    geardf(kdtf, newcoords)
  } else if (return_type == "geardf") {
    geardf(kdtf, newcoords)
  } else {
    stop("Unexpected return_type case")
  }
}

# return_gearPredict = function(kdtf, newcoords, coordnames) {
#   kdtf = cbind(kdtf, newcoords)
#   attr(kdtf, "coordnames") = coordnames
#   class(kdtf) = c("gearPredict", "data.frame")
#   return(kdtf)
# }

return_SpatialPointsDataFrame = function(kdtf, newcoords, coordnames) {
  sp::coordinates(kdtf) = newcoords
  sp::coordnames(kdtf) = coordnames
  return(kdtf)
}

return_sf = function(kdtf, newcoords, coordnames) {
    kdtf = cbind(kdtf, newcoords)
    return(sf::st_as_sf(kdtf, coords = coordnames))
}

