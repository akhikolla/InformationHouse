# needed by roxygen to add Rcpp stuff to NAMESPACE
#' @useDynLib apcf, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL


.onAttach <- function(libname, pkgname){
  packageStartupMessage(paste0("Linking to GEOS ", geos_version(),
                               " and GDAL ", gdal_version()))
}

apcf_extSoftVersion <- function() {
  structure(c(geos_version(), gdal_version()),
            names = c("GEOS", "GDAL"))
}

.onUnload <- function (libpath) {
  library.dynam.unload("apcf", libpath)
}
