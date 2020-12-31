#' Ecological-niche factor analysis
#'
#' Performs ecological-niche factor analysis using environmental raster data and
#' species presence data.
#'
#' @aliases print.enfa show.enfa
#'
#' @param x Raster* object, typically a brick or stack of ecological raster layers,
#'  or a \code{GLcenfa} object
#' @param s.dat RasterLayer, SpatialPolygons*, or SpatialPoints* object indicating
#'   species presence or abundance
#' @param field field of \code{s.dat} that specifies presence or abundance. This
#'   is equivalent to the \code{field} argument in the \code{raster} package
#' @param fun function or character. Determines what values to assign to cells
#'   with multiple spatial features, similar to the \code{fun} argument in
#'   \code{\link[raster]{rasterize}}.  Options are 'first', 'last' (default),
#'   and 'count' (see Details)
#' @param scale logical. If \code{TRUE} then the values of the Raster* object will
#'   be centered and scaled. Depending on the resolution of the climate data and
#'   the extent of the study area, this can be quite time consuming. If running
#'   this function for multiple species, it is recommended that the climate data
#'   be scaled beforehand using the \code{\link{GLcenfa}} function
#' @param filename character. Optional filename to save the Raster* output to
#'   file. If this is not provided, a temporary file will be created for large \code{x}
#' @param progress logical. If \code{TRUE}, messages and progress bar will be
#'   printed
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized for the
#'   calculation of the covariance matrices
#' @param n numeric. Number of CPU cores to utilize for parallel processing
#' @param cl optional cluster object
#' @param keep.open logical. If \code{TRUE} and \code{parallel = TRUE}, the
#'   cluster object will not be closed after the function has finished
#' @param ... Additional arguments for \code{\link[raster]{writeRaster}}
#'
#' @details
#' The \code{cnfa} function is not to be confused with the \code{\link{enfa}}
#' function. \code{enfa} performs ENFA as described by Hirzel et al. (2002) and
#' Basille et al. (2008), and is offered as an alternative to the \code{enfa}
#' function in the \code{adehabitatHS} package. \code{CENFA::enfa} will give
#' different results than \code{adehabitatHS::enfa} for versions of \code{adehabitatHS}
#' 0.3.13 or earlier, however, for two primary reasons.
#'
#' First, \code{CENFA::enfa} corrects a minor mistake in the calculation of
#' the species covariance matrix. This correction influences the values of the
#' coefficients of specialization in each ecological variable, which will lead to
#' a different interpretation of the degree of specialization. Second, we define
#' the overall marginality \eqn{M} as the norm of the marginality factor \code{mf},
#' rather than the square of the norm of \code{mf}.
#'
#' The default \code{fun = 'last'} gives equal weight to each occupied cell.
#' If multiple species observations occur in the same cell, the cell will only
#' be counted once. \code{fun = 'count'} will weight the cells by the number
#' of observations.
#'
#' If there is too much correlation between the layers of \code{x}, the global
#' covariance matrix will be singular, and the overall marginality will not be
#' meaningful. In this case, a warning is issued and \code{marginality} is
#' returned as \code{NA}.
#'
#' @examples
#' mod1 <- enfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#'
#' # using GLcenfa as an initial step
#' # for multi-species comparison
#'
#' glc <- GLcenfa(x = climdat.hist)
#' mod2 <- enfa(x = glc, s.dat = ABPR, field = "CODE")
#' all.equal(m.factor(mod1), m.factor(mod2))
#'
#' @return Returns an S4 object of class \code{enfa} with the following components:
#' \describe{
#'   \item{call}{Original function call}
#'   \item{mf}{Marginality factor. Vector that describes the location of the
#'    species Hutchinsonian niche relative to the global niche}
#'   \item{marginality}{Magnitude of the marginality factor}
#'   \item{sf}{Specialization factor. Vector of eigenvalues of specialization}
#'   \item{specialization}{Square root of the mean of the specialization factor}
#'   \item{sf.prop}{Vector representing the proportion of specialization found in each
#'   ENFA factor}
#'   \item{co}{A matrix describing the amount of marginality and specialization
#'    on each ENFA factor}
#'   \item{ras}{RasterBrick of transformed climate values, with p layers}
#'   \item{weights}{Raster layer of weights used for ENFA calculation}
#' }
#'
#' @references
#' Basille, Mathieu, et al. Assessing habitat selection using multivariate
#' statistics: Some refinements of the ecological-niche factor analysis. Ecological
#' Modelling 211.1 (2008): 233-240.
#'
#' Hirzel, Alexandre H., et al. Ecological-niche factor analysis: how to compute
#' habitat-suitability maps without absence data?. Ecology 83.7 (2002): 2027-2036.
#'
#' @seealso \code{\link{GLcenfa}}, \code{\link{cnfa}}
#'
# @importFrom raster rasterize
#' @importFrom sp identicalCRS
#' @export

setGeneric("enfa", function(x, s.dat, ...){
  standardGeneric("enfa")})

#' @rdname enfa
setMethod("enfa",
          signature(x = "GLcenfa", s.dat = "Raster"),
          function(x, s.dat, filename = "", progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

            call <- sys.call(sys.parent())
            call <- match.call(enfa, call)

            if(file.exists(filename)) {
              params <- list(...)
              if(length(params) > 0) {
                if(is.null(params$overwrite)) {
                  stop(paste0(filename, " exists. use 'overwrite=TRUE' if you want to overwrite it"))
                } else if(!(params$overwrite)) {
                  stop(paste0(filename, " exists. use 'overwrite=TRUE' if you want to overwrite it"))
                }
              } else stop(paste0(filename, " exists. use 'overwrite=TRUE' if you want to overwrite it"))
            }

            if (nlayers(s.dat) > 1) stop('"s.dat" should be a single RasterLayer')
            if (!identicalCRS(raster(x), s.dat)) stop("climate and species projections do not match")
            ras <- raster(x)
            ext <- extent(ras)
            ext.s <- extent(s.dat)
            if (is.null(intersect(ext, ext.s))) stop("climate and species data do not overlap")
            if (raster::union(ext, ext.s) != ext) stop("extent of species data not contained within extent of climate data")

            x.crop <- crop(ras, ext.s)

            filename <- trim(filename)
            if (!canProcessInMemory(x.crop) && filename == '') {
              filename <- rasterTmpFile()
            }

            if (canProcessInMemory(x.crop)){
              pres <- which(!is.na(values(s.dat)) & !is.na(values(max(x.crop))))
              S <- values(x.crop)[pres,]
              nS <- nrow(S)
              Rg <- x@cov
              p <- values(s.dat)[pres]
              p.sum <- sum(p)
              mar <- apply(S, 2, function(x) sum(x * p)) / p.sum
              if (progress) cat("Calculating species covariance matrix...\n")
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm)/ (p.sum - 1)
            } else {
              x.mask <- mask(x.crop, s.dat)
              Rg <- x@cov
              p.sum <- cellStats(s.dat, sum)
              DpS <- x.mask * s.dat
              mar <- cellStats(DpS, sum) / p.sum
              if (progress) cat("Calculating species covariance matrix...\n")
              if (parallel) {
                if (!keep.open) on.exit(closeAllConnections())
                if (missing(cl) && n > 1) cl <- snow::makeCluster(getOption("cl.cores", n))
              }
              Sm <- parScale(x.mask, center = mar, scale = F, parallel = parallel, n = n, progress = F, keep.open = keep.open, cl = cl)
              Rs <- parCov(x = Sm, w = s.dat, parallel = parallel, n = n, progress = progress, keep.open = keep.open, cl = cl)
            }

            cZ <- nlayers(ras)
            m <- sqrt(as.numeric(t(mar) %*% mar))
            if (max(Im(eigen(Rs)$values)) > 1e-05) stop("complex eigenvalues. Try removing correlated variables.")
            eigRs <- lapply(eigen(Rs), Re)
            keep <- (eigRs$values > 1e-09)
            Rs12 <- eigRs$vectors[, keep] %*% diag(eigRs$values[keep]^(-0.5)) %*% t(eigRs$vectors[, keep])
            #Rs12 <- eigRs$vectors %*% diag(eigRs$values^(-0.5)) %*% t(eigRs$vectors)
            W <- Rs12 %*% Rg %*% Rs12
            z <- Rs12 %*% mar
            y <- z/sqrt(sum(z^2))
            H <- (diag(cZ) - y %*% t(y)) %*% W %*% (diag(cZ) - y %*% t(y))
            sf <- eigen(H)$values[-cZ]
            s.p <- (t(mar) %*% Rg %*% mar) / (t(mar) %*% Rs %*% mar)
            s <- c(s.p, sf)
            spec <- sqrt(mean(s))
            s.p <- abs(s)/sum(abs(s))
            v <- Re(eigen(H)$vectors)
            U <- matrix(nrow = cZ, ncol = cZ)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            U[, -1] <- sweep(u, 2, norw, "/")
            U[, 1] <- mar
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if (canProcessInMemory(x.crop)) {
              s.ras <- brick(x.crop)
              if (progress) cat("Creating factor rasters...")
              values(s.ras)[pres, ] <- S %*% U
              names(s.ras) <- nm
            } else {
              if (progress) cat("Creating factor rasters...")
              if(parallel) {
                f1 <- function(x) x %*% U
                s.ras <- clusterR(x.mask, fun = .calc, args = list(fun = f1, forceapply = T, names = nm), cl = cl, filename = filename, ...)
              } else {
                s.ras <- .calc(x.mask, fun = f1, forceapply = T, filename = filename, names = nm, ...)
              }
            }
            colnames(U) <- names(s.p) <- names(s) <- nm
            rownames(U) <- names(mar) <- names(ras)

            mod <- methods::new("enfa", call = call, mf = mar, marginality = m, sf = s,
                                 specialization = spec, sf.prop = s.p, co = U, cov = Rs, ras = s.ras, weights = s.dat)
            return(mod)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(x = "GLcenfa", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", filename = "", progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

            call <- sys.call(sys.parent())
            call <- match.call(enfa, call)

            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if (!identicalCRS(raster(x), s.dat)) stop("climate and species projections do not match")
            ras <- raster(x)
            ext <- extent(ras)
            ext.s <- extent(s.dat)
            if (is.null(intersect(ext, ext.s))) stop("climate and species data do not overlap")
            if (raster::union(ext, ext.s) != ext) stop("extent of species data not contained within extent of climate data")

            x.crop <- crop(ras, ext.s)
            s.dat.ras <- rasterize(s.dat, x.crop, field = field, fun = fun)

            mod <- enfa(x = x, s.dat = s.dat.ras, filename = filename, progress = progress, parallel = parallel, n = n, ...)
            mod@call <- call
            return(mod)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(x = "Raster", s.dat = "Raster"),
          function(x, s.dat, scale = TRUE, filename = "", progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

            call <- sys.call(sys.parent())
            call <- match.call(enfa, call)

            if (nlayers(s.dat) > 1) stop('"s.dat" should be a single RasterLayer')
            if (!identicalCRS(x, s.dat)) stop("projections do not match")
            if (is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if (raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            if (scale) {
              if (progress) cat("Scaling raster data...\n")
              x <- GLcenfa(x = x, center = T, scale = T, progress = progress, parallel = parallel, n = n)
            } else x <- GLcenfa(x = x, center = F, scale = F, progress = progress, parallel = parallel, n = n)

            mod <- enfa(x = x, s.dat = s.dat, filename = filename, progress = progress, parallel = parallel, n = n, ...)
            mod@call <- call
            return(mod)
          }
)

#' @rdname enfa
setMethod("enfa",
          signature(x = "Raster", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

            call <- sys.call(sys.parent())
            call <- match.call(enfa, call)

            if (!inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if (!identicalCRS(x, s.dat)) stop("projections do not match")
            if (is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if (raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)
            if (scale) {
              if (progress) cat("Scaling raster data...\n")
              x <- GLcenfa(x = x, center = T, scale = T, progress = progress, parallel = parallel, n = n)
            } else x <- GLcenfa(x = x, center = F, scale = F, progress = progress, parallel = parallel, n = n)

            mod <- enfa(x = x, s.dat = s.dat.ras, filename = filename, progress = progress, parallel = parallel, n = n, ...)
            mod@call <- call
            return(mod)
          }
)

# setMethod("enfa",
#           signature(x = "Raster", s.dat = "sf"),
#           function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", progress = FALSE, parallel = FALSE, n = 1, ...){
#
#             call <- sys.call(sys.parent())
#
#             if (!identicalCRS(x, s.dat)) stop("projections do not match")
#             if (is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
#             if (raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")
#
#             s.dat <- as(s.dat, "Spatial")
#             if (!inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('geometry of "s.dat" should be of class "sfc_POLYGON", "sfc_MULTIPOLYGON", "sfc_POINT", or "sfc_MULTIPOINT"')
#             s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)
#             if (scale) {
#               if (progress) cat("Scaling raster data...\n")
#               x <- GLcenfa(x = x, center = T, scale = T, progress = progress, parallel = parallel, n = n)
#             } else x <- GLcenfa(x = x, center = F, scale = F, progress = progress, parallel = parallel, n = n)
#
#             enfa(x = x, s.dat = s.dat.ras, filename = filename, progress = progress, parallel = parallel, n = n, ...)
#           }
# )
