#' Climatic departure
#'
#' This function quantifies the amount of change between historical and future
#' climate conditions inside a species' habitat.
#'
#' @param x Raster* object, typically a brick or stack of historical climate
#'   raster layers or a brick of absolute differences (see Details)
#' @param y  Raster* object, future climate values with the same layers as \code{x}
#' @param s.dat SpatialPolygons*, sf, or cnfa object detailing species presence
#' @param field field of \code{s.dat} that specifies presence. This is
#'   equivalent to the \code{field} argument of \code{raster::rasterize}. Options
#'   are 'first', 'last' (default), and 'count'
#' @param fun function or character. Determines what values to assign to cells
#'   with multiple spatial features, similar to the \code{fun} argument in
#'   \code{\link[raster]{rasterize}}
#' @param center logical. If \code{TRUE} then the values of \code{x} and
#'   \code{y} will be centered on the means of the historical
#'   climate data
#' @param scale logical. If \code{TRUE} then the values of \code{x} and
#'   \code{y} will be scaled by the sds of the historical
#'   climate data
#' @param filename character. Optional filename to save the Raster* output to
#'   file. If this is not provided, a temporary file will be created for large \code{x}
#' @param progress logical. If \code{TRUE}, messages and progress bar will be
#'   printed
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Optional number of CPU cores to utilize for parallel processing
#' @param ... Additional arguments for \code{\link[raster]{clusterR}}
#'
#' @examples
#' dep1 <- departure(x = climdat.hist, y = climdat.fut, s.dat = ABPR, field = "CODE")
#'
#' # using difRaster as an initial step
#' # for multi-species comparison
#'
#' gld <- GLdeparture(x = climdat.hist, y = climdat.fut)
#' dep2 <- departure(x = gld, s.dat = ABPR, field = "CODE")
#'
#'# same results either way
#' all.equal(dep1@df, dep2@df)
#'
#' @return Returns an S4 object of class \code{departure} with the following slots:
#' \describe{
#'   \item{call}{Original function call}
#'   \item{df}{Departure factor. Vector of length p that describes the amount of
#'    departure between future and historical conditions for each climate variable}
#'   \item{departure}{Magnitude of the departure factor}
#'   \item{g.cov}{p x p historical global covariance matrix}
#'   \item{ras}{RasterBrick of climate departures, with p layers}
#'   \item{weights}{Raster layer of weights used for departure calculation}
#' }
#'
#' @details
#'  For comparisons of multiple species in the same study area, it is much more
#'  efficient to first construct a Raster* object of absolute differences between
#'  the historical and future values, so that the differences do not need to be
#'  recalculated for each species. This can be achieved with by passing \code{x}
#'  and \code{y} to the \code{difRaster} function, and then passing the
#'  results to the \code{departure} function.
#'
#'  When only one Raster* object is supplied, it is assumed that \code{x} is
#'  a Raster* object containing the absolute differences of a historical and
#'  future dataset.
#'
#' @references
#' Rinnan, D. Scott and Lawler, Joshua. Climate-niche factor analysis: a spatial
#' approach to quantifying species vulnerability to climate change. Ecography (2019):
#' \href{https://doi.org/10.1111/ecog.03937}{doi:10.1111/ecog.03937}.
#'
#' @include CENFA.R cnfa-class.R GLcenfa-class.R
#'
#' @export
#' @importFrom stats sd
#' @importFrom parallel detectCores
# @importFrom raster overlay crop

setGeneric("departure", function(x, y, s.dat, ...) {
  standardGeneric("departure")
})

#' @rdname departure
setMethod("departure",
          signature(x = "GLdeparture", y = "missing", s.dat = "cnfa"),
          function(x, s.dat, filename = '', ...){

            call <- sys.call(sys.parent())
            call <- match.call(departure, call)

            s.dat.ras <- s.dat@weights
            ras <- x@global_difras
            ext <- extent(ras)
            ext.s <- extent(s.dat.ras)

            if (!identicalCRS(ras, s.dat.ras)) stop("climate and species projections do not match")
            if (is.null(intersect(ext, ext.s))) stop("climate and species data do not overlap")
            if (raster::union(ext, ext.s) != ext) stop("extent of species data not contained within extent of climate data")

            filename <- trim(filename)
            if (!canProcessInMemory(ras)) {
              if (filename == '') filename <- rasterTmpFile()
            }

            Rg <- x@cov
            nm <- names(ras)
            x.dif <- crop(ras, s.dat.ras)
            names(x.dif) <- nm
            x.dif <- mask(x.dif, s.dat.ras, filename = filename, ...)
            w <- s.dat.ras / (cellStats(s.dat.ras, sum, na.rm = T) - 1)
            x.dif.w <- overlay(x = x.dif, y = w, fun = function(x,y) {return(x*y)})
            d <- cellStats(x.dif.w, sum)
            names(d) <- nm
            D <- tryCatch(sqrt(as.numeric(t(d) %*% d)),
                          error = function(e){
                            message("Warning: global covariance matrix not invertible. Overall departure will not be computed.")
                            return(as.numeric(NA))})

            depart <- methods::new("departure", call = call, df = d, departure = D, g.cov = Rg, ras = x.dif, weights = s.dat.ras)
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x = "GLdeparture", y = "missing", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", filename = '', ...){

            call <- sys.call(sys.parent())
            call <- match.call(departure, call)

            ras <- x@global_difras
            ext <- extent(ras)
            ext.s <- extent(s.dat)

            if (!inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if (is.null(intersect(ext, ext.s))) stop("climate and species data do not overlap")
            if (raster::union(ext, ext.s) != ext) stop("extent of species data not contained within extent of climate data")

            filename <- trim(filename)
            if (!canProcessInMemory(ras)) {
              if (filename == '') filename <- rasterTmpFile()
            }

            Rg <- x@cov
            s.dat.ras <- rasterize(s.dat, ras, field = field, fun = fun)

            nm <- names(x)
            x.dif <- crop(ras, s.dat.ras)
            names(x.dif) <- nm
            x.dif <- mask(x.dif, s.dat.ras, filename = filename, ...)
            w <- s.dat.ras / (cellStats(s.dat.ras, sum, na.rm = T) - 1)
            x.dif.w <- overlay(x = x.dif, y = w, fun = function(x,y) {return(x*y)})
            d <- cellStats(x.dif.w, sum)
            names(d) <- nm
            D <- tryCatch(sqrt(as.numeric(t(d) %*% d)),
                          error = function(e){
                            message("Warning: global covariance matrix not invertible. Overall departure will not be computed.")
                            return(as.numeric(NA))})

            depart <- methods::new("departure", call = call, df = d, departure = D, g.cov = Rg, ras = x.dif, weights = s.dat.ras)
            return(depart)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x = "Raster", y = "Raster", s.dat = "cnfa"),
          function(x, y, s.dat, center = TRUE, scale = TRUE, filename = '', progress = FALSE, parallel = FALSE, n = 1, ...) {

            call <- sys.call(sys.parent())
            call <- match.call(departure, call)

            GLdep <- GLdeparture(x, y, center = center, scale = scale, progress = progress, parallel = parallel, n = n)
            dep <- departure(x = GLdep, s.dat = s.dat, filename = filename, ...)
            dep@call <- call
            return(dep)
          }
)

#' @rdname departure
setMethod("departure",
          signature(x = "Raster", y = "Raster", s.dat = "Spatial"),
          function(x, y, s.dat, center = TRUE, scale = TRUE, filename = '', progress = FALSE, parallel = FALSE, n = 1, ...) {

            call <- sys.call(sys.parent())
            call <- match.call(departure, call)

            GLdep <- GLdeparture(x, y, center = center, scale = scale, progress = progress, parallel = parallel, n = n)
            dep <- departure(x = GLdep, s.dat = s.dat, filename = filename, ...)
            dep@call <- call
            return(dep)
          }
)
