#' Climate-niche factor analysis
#'
#' Performs climate-niche factor analysis using climate raster data and species
#' presence data.
#'
#' @aliases print.cnfa show.cnfa
#'
#' @param x Raster* object, typically a brick or stack with p climate
#'   raster layers, or a \code{GLcenfa} object
#' @param s.dat RasterLayer, SpatialPolygons*, or SpatialPoints* object indicating
#'   species presence or abundance
#' @param field field of \code{s.dat} that specifies presence or abundance. This
#'   is equivalent to the \code{field} argument in \code{\link[raster]{rasterize}}
#' @param fun function or character. Determines what values to assign to cells
#'   with multiple spatial features, similar to the \code{fun} argument in
#'   \code{\link[raster]{rasterize}}.  Options are 'first', 'last' (default),
#'   and 'count' (see Details)
#' @param scale logical. If \code{TRUE} then the values of \code{x} will get
#'   centered and scaled. Depending on the resolution of the climate data and
#'   the extent of the study area, this can be quite time consuming. If running
#'   this function for multiple species, it is recommended that the data be
#'   scaled beforehand using the \code{\link{GLcenfa}} function
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
#' @examples
#' mod1 <- cnfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#'
#' # using GLcenfa as an initial step
#' # for multi-species comparison
#'
#' glc <- GLcenfa(x = climdat.hist)
#' mod2 <- cnfa(x = glc, s.dat = ABPR, field = "CODE")
#'
#'# same results either way
#' all.equal(m.factor(mod1), m.factor(mod2))
#' all.equal(s.factor(mod1), s.factor(mod2))
#'
#' @return Returns an S4 object of class \code{cnfa} with the following components:
#' \describe{
#'   \item{call}{Original function call}
#'   \item{mf}{Marginality factor. Vector of length p that describes the location
#'   of the species Hutchinsonian niche relative to the global niche}
#'   \item{marginality}{Magnitude of the marginality factor}
#'   \item{sf}{Sensitivity factor. Vector of length p that describes the amount of
#'    sensitivity for each climate variable}
#'   \item{sensitivity}{Square root of the mean of the sensitivity factor}
#'   \item{eig}{Named vector of eigenvalues of specialization for each CNFA factor}
#'   \item{co}{A p x p matrix describing the amount of marginality and specialization
#'    in each CNFA factor.}
#'   \item{cov}{p x p species covariance matrix}
#'   \item{g.cov}{p x p global covariance matrix}
#'   \item{ras}{RasterBrick of transformed climate values, with p layers}
#'   \item{weights}{Raster layer of weights used for CNFA calculation}
#' }
#'
#' @details
#' The \code{cnfa} function is not to be confused with the
#' \code{\link{enfa}} function. \code{enfa} performs ENFA as described by Hirzel
#' et al. (2002) and Basille et al. (2008), and is offered as an alternative to
#' the \code{enfa} function in the \code{adehabitatHS} package. There are
#' several key differences between ENFA and CNFA.
#'
#' Whereas ENFA returns a \strong{specialization factor} that describes
#' the specialization in each \strong{ENFA factor}, CNFA returns a
#' \strong{sensitivity factor} \code{sf} that describes the sensitivity in each
#' \strong{environmental variable}. This makes the sensitivity factor more
#' directly comparable to the marginality factor \code{mf}, because their
#' dimensions are identical. Sensitivity is calculated by a weighted sum
#' of the amount of specialization found in each CNFA factor, \emph{including}
#' the marginality factor. As such, the sensitivity factor offers a more complete
#' measure of specialization than ENFA's specialization factor, which does
#' not calculate the amount of specialization found in the marginality factor.
#' As such, CNFA's overall sensitivity (found in the slot \code{sensitivity})
#' offers a more complete measure of niche specialization than ENFA's overall
#' specialization (found in the slot \code{specialization}).
#'
#' The default \code{fun = 'last'} gives equal weight to each occupied cell.
#' If multiple species observations occur in the same cell, the cell will only
#' be counted once. \code{fun = 'count'} will weight the cells by the number
#' of observations.
#'
#' If there is too much correlation between the layers of \code{x}, the global
#' covariance matrix will be singular, and the overall marginality and overall
#' sensitivity will not be meaningful. In this case, a warning is issued,
#' and \code{marginality} and \code{sensitivity} are both returned as \code{NA}.
#'
#' @references
#' Rinnan, D. Scott and Lawler, Joshua. Climate-niche factor analysis: a spatial
#' approach to quantifying species vulnerability to climate change. Ecography (2019):
#' \href{https://doi.org/10.1111/ecog.03937}{doi:10.1111/ecog.03937}.
#'
#' Basille, Mathieu, et al. Assessing habitat selection using multivariate
#' statistics: Some refinements of the ecological-niche factor analysis. Ecological
#' Modelling 211.1 (2008): 233-240.
#'
#' Hirzel, Alexandre H., et al. Ecological-niche factor analysis: how to compute
#' habitat-suitability maps without absence data?. Ecology 83.7 (2002): 2027-2036.
#'
#' @seealso \code{\link{GLcenfa}}, \code{\link{enfa}}
#'
#' @export
#' @name cnfa
#'
#' @importFrom stats cov
#' @importFrom methods as is
# @importFrom raster rasterTmpFile intersect extent values mask crop values<-

setGeneric("cnfa", function(x, s.dat, ...){
  standardGeneric("cnfa")})

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "GLcenfa", s.dat = "Raster"),
          function(x, s.dat, filename = "", progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

            call <- sys.call(sys.parent())
            call <- match.call(cnfa, call)

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
            nS <- length(which(values(s.dat) > 0))
            if(nS == 1) stop("CNFA is not meaningful for single observations")
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

            if (canProcessInMemory(x.crop) && !parallel) {
              pres <- which(!is.na(values(s.dat)) & !is.na(values(max(x.crop))))
              S <- values(x.crop)[pres, ]
              nS <- nrow(S)
              Rg <- x@cov
              p <- values(s.dat)[pres]
              p.sum <- sum(p)
              mar <- apply(S, 2, function(x) sum(x * p)) / p.sum
              if (progress) cat("Calculating species covariance matrix...\n")
              Sm <- sweep(S, 2, mar)
              DpSm <- apply(Sm, 2, function(x) x * p)
              Rs <- crossprod(Sm, DpSm) * 1 / (p.sum - 1)
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
            s <- Re(eigen(H)$values)[-cZ]
            s.p <- (t(mar) %*% Rg %*% mar) / (t(mar) %*% Rs %*% mar)
            s <- c(s.p, s)
            #s.p <- abs(s) / sum(abs(s))
            v <- Re(eigen(H)$vectors)
            U <- matrix(nrow = cZ, ncol = cZ)
            u <- as.matrix((Rs12 %*% v)[, 1:(cZ-1)])
            norw <- sqrt(diag(t(u) %*% u))
            U[, -1] <- sweep(u, 2, norw, "/")
            #co[, 1] <- mar / m
            U[, 1] <- mar
            V <- sweep(abs(U), 2, colSums(abs(U)), "/")
            sf <- as.numeric(V %*% s)
            #sens <- sqrt(as.numeric(t(sf) %*% sf)/cZ)
            sens <- sqrt(mean(sf))
            nm <- c("Marg", paste0("Spec", (1:(cZ-1))))
            if (canProcessInMemory(x.crop) & !parallel){
              s.ras <- brick(x.crop)
              if (progress) cat("Creating factor rasters...")
              values(s.ras)[pres, ] <- S %*% U
              names(s.ras) <- nm
            } else {
              if (progress) cat("Creating factor rasters...")
              f1 <- function(x) x %*% U
              if(parallel) {
                s.ras <- clusterR(x.mask, fun = .calc, args = list(fun = f1, forceapply = T, names = nm), cl = cl, filename = filename, ...)
              } else {
                s.ras <- .calc(x.mask, fun = f1, forceapply = T, filename = filename, names = nm, ...)
              }
            }
            colnames(U) <- names(s) <- nm
            rownames(U) <- names(sf) <- names(mar) <- names(ras)
            if (!keep.open || missing(cl)) snow::stopCluster(cl)

            mod <- methods::new("cnfa", call = call, mf = mar, marginality = m, sf = sf,
                                 sensitivity = sens, eig = s, co = U, cov = Rs, g.cov = Rg, ras = s.ras, weights = s.dat)
            return(mod)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "GLcenfa", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", filename = "", progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

            call <- sys.call(sys.parent())
            call <- match.call(cnfa, call)

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

            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if (!identicalCRS(raster(x), s.dat)) stop("climate and species projections do not match")
            ras <- raster(x)
            ext <- extent(ras)
            ext.s <- extent(s.dat)
            if (is.null(intersect(ext, ext.s))) stop("climate and species data do not overlap")
            if (raster::union(ext, ext.s) != ext) stop("extent of species data not contained within extent of climate data")

            x.crop <- crop(ras, ext.s)
            s.dat.ras <- rasterize(s.dat, x.crop, field = field, fun = fun)

            mod <- cnfa(x = x, s.dat = s.dat.ras, filename = filename, progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open, ...)
            mod@call <- call
            return(mod)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "Raster", s.dat = "Raster"),
          function(x, s.dat, scale = TRUE, filename = "", progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = keep.open, ...){

            call <- sys.call(sys.parent())
            call <- match.call(cnfa, call)

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
            if(!identicalCRS(x, s.dat)) stop("projections do not match")
            if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if(raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            if (scale) {
              if (progress) cat("Scaling raster data...\n")
              #x <- parScale(x, center = T, scale = T, parallel = parallel, n = n, progress = progress)
              x <- GLcenfa(x = x, center = T, scale = T, progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open)
            } else x <- GLcenfa(x = x, center = F, scale = F, progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open)

            mod <- cnfa(x = x, s.dat = s.dat, filename = filename, progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open, ...)
            mod@call <- call
            return(mod)
          }
)

#' @rdname cnfa
setMethod("cnfa",
          signature(x = "Raster", s.dat = "Spatial"),
          function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

            call <- sys.call(sys.parent())
            call <- match.call(cnfa, call)

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

            if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('"s.dat" should be a "SpatialPolygons*" or "SpatialPoints*" object')
            if(!identicalCRS(x, s.dat)) stop("projections do not match")
            if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
            if(raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")

            s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)
            if (scale) {
              x <- GLcenfa(x = x, center = T, scale = T, progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open)
            } else x <- GLcenfa(x = x, center = F, scale = F, progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open)

            mod <- cnfa(x = x, s.dat = s.dat.ras, filename = filename, progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open, ...)
            mod@call <- call
            return(mod)
          }
)

# setMethod("cnfa",
#           signature(x = "Raster", s.dat = "sf"),
#           function(x, s.dat, field, fun = "last", scale = TRUE, filename = "", progress = FALSE, parallel = FALSE, n = 1, ...){
#             if (!requireNamespace("sf")) {
#               warning('cannot do this because sf is not available')
#             }
#
#             call <- sys.call(sys.parent())
#
#             if(!identicalCRS(x, s.dat)) stop("projections do not match")
#             if(is.null(intersect(extent(x), extent(s.dat)))) stop("climate and species data do not overlap")
#             if(raster::union(extent(x), extent(s.dat)) != extent(x)) stop("extent of species data not contained within extent of climate data")
#
#             s.dat <- as(s.dat, "Spatial")
#             if (! inherits(s.dat, c('SpatialPolygons', 'SpatialPoints'))) stop('geometry of "s.dat" should be of class "sfc_POLYGON", "sfc_MULTIPOLYGON", "sfc_POINT", or "sfc_MULTIPOINT"')
#
#             s.dat.ras <- rasterize(s.dat, raster(x), field = field, fun = fun)
#             if (scale) {
#               if (progress) cat("Scaling raster data...\n")
#               x <- GLcenfa(x = x, center = T, scale = T, progress = progress, parallel = parallel, n = n)
#             } else x <- GLcenfa(x = x, center = F, scale = F, progress = progress, parallel = parallel, n = n)
#
#             cnfa(x = x, s.dat = s.dat.ras, filename = filename, progress = progress, parallel = parallel, n = n, ...)
#           }
# )
