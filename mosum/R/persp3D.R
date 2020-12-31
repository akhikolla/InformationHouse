#' 3D Visualisation of multiscale MOSUM statistics
#' 
#' 3D Visualisation of multiscale MOSUM statistics.
#' @param x a numeric input data vector
#' @param threshold string indicating which threshold should be used for normalisation of
#' MOSUM statistics computed with different bandwidths.
#' By default, it is chosen from the asymptotic distribution at the given significance level \code{alpha}.
#' Alternatively it is possible to parse a user-defined numerical value with \code{threshold.custom}; see also Details.
#' @param mosum.args a named list containing further arguments
#' to be parsed to the respective \code{mosum} function calls, see \link[mosum]{mosum};
#' the bandwidths are chosen by the function and should not be given as an argument in \code{mosum.args}
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff \code{threshold = "critical.value"}
#' @param threshold.function function object of form \code{function(G)}, to compute a
#' threshold of significance for different bandwidths G; use iff \code{threshold='custom'}
#' @param pal.name a string containing the name of the ColorBrewer palette to be used; 
#' sequential palettes are recommended.
#' See \code{RColorBrewer::brewer.pal.info} for details
#' @param expand expansion factor applied to the z coordinates
#' @param theta azimuthal angle defining the viewing direction
#' @param phi colatitude angle defining the viewing direction
#' @param xlab,ylab,zlab,ticktype graphical parameters
#' @param NAcol coloring parameter
#' @param ... further arguments to be passed to function call of \link[plot3D]{persp3D}
#' @return see \link[plot3D]{persp3D}
#' @details The visualisation is based on \link[plot3D]{persp3D}.
#' MOSUM statistics computed with different bandwidths are rescaled
#' for making them visually comparable.
#' Rescaling is done either by dividing by their respective critical value at the significance level \code{alpha}
#' (iff \code{threshold = "critical.value"}) or by a custom value given by \code{threshold.function}
#' (iff \code{threshold = "custom"}).
#' By default, \code{clim} argument of \link[plot3D]{persp3D} is given so that the three lightest 
#' (for sequential palettes) hues indicate insignificance of the corresponding MOSUM statistics,
#' while darker hues indicate the presence of significant changes.
#' @examples
#' \dontrun{
#' # If you run the example be aware that this may take some time
#' print("example may take some time to run")
#' 
#' x <- testData(model = "blocks", seed = 1234)$x
#' persp3D.multiscaleMosum(x, mosum.args = list(boundary.extension = FALSE))
#' }
#' @importFrom plot3D persp3D
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @export
persp3D.multiscaleMosum <- function(x, mosum.args=list(), 
                                    threshold = c('critical.value', 'custom')[1],
                                    alpha=0.1, threshold.function = NULL,
                                    pal.name='YlOrRd',
                                    expand=.2, theta=120, phi=20, xlab='G', 
                                    ylab='time', zlab='MOSUM', ticktype='detailed', 
                                    NAcol='#800000FF', ...) {
  stopifnot(is.null(mosum.args$G))
  stopifnot(is.null(mosum.args$x))
  stopifnot(is.null(mosum.args$G.right))
  stopifnot(is.element(pal.name, rownames(RColorBrewer::brewer.pal.info)))
  col.number <- brewer.pal.info$maxcolors[is.element(rownames(brewer.pal.info), pal.name)]
  mycols <- brewer.pal(col.number, pal.name)
  clim <- c(0, col.number/3)
  
  n <- length(x)
  G <- 8:floor(3*sqrt(n))
  
  if (class(G) == 'integer' || class(G) == 'numeric') {
    grid <- multiscale.grid(G, method='concatenate')
  } else stop('Expecting a vector of numbers')

  # Collect all statistics for visualization
  m <- array(list(), nrow(grid$grid))
  for (i in seq_len(nrow(grid$grid))) {
    G <- grid$grid[[i,1]]
    argList <- mosum.args
    argList$x <- x
    argList$G <- G
    if(threshold == "custom"){
      argList$threshold <- "custom"
      argList$threshold.custom <- threshold.function(G, n, alpha)
    }
    m[[i]] <- do.call(mosum, argList)
  }
  
  zz <- t(sapply(m, function(z) z$stat / z$threshold.value))
  xx <- grid$grid[,1]
  yy <- seq_len(length(x)) 
  persp3D(x=xx, y=yy, z=zz, expand=expand, 
          theta=theta, phi=phi, xlab=xlab, ylab=ylab, 
          zlab=zlab, ticktype=ticktype, 
          clim=clim, NAcol=NAcol, col=mycols, ...)
}