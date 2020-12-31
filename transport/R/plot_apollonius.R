#' Plot Apollonius Diagram
#'
#' Plots the Apollonius diagram, a.k.a. (additively) weighted Voronoi diagram, based
#' on a matrix of points (centers) in 2d and their weights.
#'
#' @param points A two-column matrix containing the 2d points.
#' @param weights A vector of weights for the points.
#' @param show_points Logical. Should the points be displayed in the plot?
#'        Defaults to TRUE.
#' @param show_weights Logical. Should the weights be displayed in the plot?
#'        Defaults to TRUE.
#' @param add_to_weights A value added to the weights to make the plot
#'        more informative.
#' @param add Logical. Should the plot be added to the current device?
#'        Defaults to FALSE.
#' @param col The colour for the cell boundaries.
#' @param lwd The line width for the cell boundaries.
#' @param ... Further parameters to the base plot if \code{add} is \code{FALSE}.
#'
#' @details
#' For points \eqn{x_1, \ldots, x_n}{x_1, ..., x_n} with weights \eqn{w_1, \ldots, w_n}
#' The $i$-th cell of the Apollonius diagram contains all the points x that satisfy
#' \deqn{\|x-x_i\|-w_i < \|x-x_j\|-w_j}{||x-x_i|| - w_i < ||x-x_j|| - w_j} 
#' for all  \eqn{j \neq i}{j != i}. Its boundaries are hyperbola segments.
#' 
#' If \code{show_weights} is \code{TRUE}, grey circles of radii \code{weights + add_to_weights}
#' are plotted around the points. Negative radii are set to zero.
#'
#' @note This function requires the Computational Geometry Algorithms Library (CGAL),
#'       available at \url{https://www.cgal.org}. Adapt the file src/Makevars according
#'       to the instructions given there and re-install from source.
#'
#' @author Valentin Hartmann \email{valentin.hartmann@epfl.ch} (most of the code)\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de} (R-port)
#'
#' @references Menelaos Karavelas and Mariette Yvinec. 2D Apollonius Graphs 
#'             (Delaunay Graphs of Disks). In CGAL User and Reference Manual.
#'             CGAL Editorial Board, 4.12 edition, 2018
#'
#' @examples
#' \dontrun{
#' w <- c(0.731, 0.0372, 0.618, 0.113, 0.395, 0.222, 0.124, 0.101, 0.328, 0)
#' points <- matrix(runif(20), 10, 2)
#' plot_apollonius(points, w, add_to_weights = -0.1)}
#'
#' @export
# The assumption is that the points live on [0,S]x[0,T]
# And are not concentrated on a minuscule part of it for precision reasons
#
plot_apollonius <- function(points, weights, show_points = TRUE, show_weights = TRUE,
                            add_to_weights = 0, add = FALSE, col=4, lwd=1.5, ...) {
stopifnot(dim(points)[1] == length(weights) && length(weights) > 0)
if (length(weights) == 1) {
  warning("Nothing to plot.")
  return(invisible())
}
  
nocgal <- !as.logical(.Call('_transport_cgal_present', PACKAGE = 'transport'))
if (nocgal) {
  stop("Computation of semidiscrete optimal transport for p=1 requires CGAL (the Computational
Geometry Algorithms Library). Install it from https://www.cgal.org/download.html and 
adapt the file src/Makevars of this package according to the instructions given there.
Then re-install 'transport' from source as usual.")
}
  
xmin <- min(points[,1])
ymin <- min(points[,2])
sites <- cbind(points[,1]-xmin, points[,2]-ymin, weights)
# minus xmin and ymin to bring to [0,S]x[0,T]
intersections <- create_diagram(sites)  
  
resolution <- 1000
EPS = 1e-6;
PSEUDOINF = EPS^-1;

# dist1 <- function(x,y) {as.numeric(dist(rbind(x,y)))}
dist1 <- function(x,y) {sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)}

if (!add) {
  plot(points[,1],points[,2], type="n", pch=16, asp=1, xaxs="i", yaxs="i", xlab="", ylab="", ...)
}
  
for (i in 1:dim(intersections)[1]) {
# pp and qq are the points between which the current hyperbola segment lies.
  pp <- intersections[i, 1:2]
  qq <- intersections[i, 4:5]
  wp <- intersections[i, 3]
  wq <- intersections[i, 6]
  end1 <- intersections[i, 7:8]
  end2 <- intersections[i, 9:10]

  if (wp == wq) {
    # We have a straight line instead of a hyperbola.
    segments(end1[1],end1[2],end2[1],end2[2], col=col, lwd=lwd)
  } else {
    dist_qp <- dist1(pp,qq)
    diff_qp_1 <- qq[1] - pp[1]
    if (diff_qp_1 == 0) {
      sign_qp_1 <- -1
    } else {
      sign_qp_1 <- sign(diff_qp_1)
    }

    cos_qp <- abs(diff_qp_1) / dist_qp
    sin_qp <- sign_qp_1 * (qq[2] - pp[2]) / dist_qp

    a <- abs(wp-wq) / 2
    b <- sqrt(dist_qp^2/4 - a^2)
    M <- (1/2) * (qq + pp)
    A <- matrix(0, 2, 2)
    A[1,1] <- a*cos_qp + b*sin_qp
    A[2,1] <- a*sin_qp - b*cos_qp
    A[1,2] <- a*cos_qp - b*sin_qp
    A[2,2] <- a*sin_qp + b*cos_qp
    A <- (1/2) * A
  
    A_inv <- matrix(0, 2, 2)
    A_inv[1,1] <- cos_qp/a + sin_qp/b
    A_inv[2,1] <- cos_qp/a - sin_qp/b
    A_inv[1,2] <- sin_qp/a - cos_qp/b
    A_inv[2,2] <- sin_qp/a + cos_qp/b

    # Does the respective endpoint exist or does the hyperbola go to
    # infinity? (i.e. our endpoint lies on the boundary of the realArea)
    infinit_1 = is.nan(end1[1])
    infinit_2 = is.nan(end2[1])

    # transform the endpoints to the graph of x -> 1/x or x -> -1/x
    if (pp[1] >= qq[1]) {
      if (wp < wq) {
        if (infinit_1) {
          end1[1] <- EPS
        } else {
          end1 <- A_inv %*% (end1 - M)
        }
        if (infinit_2) {
          end2[1] <- PSEUDOINF
        } else {
          end2 <- A_inv %*% (end2 - M)
        }
      } else {
        if (infinit_1) {
          end1[1] = -PSEUDOINF
        } else {
          end1 = A_inv %*% (end1 - M)
        }
        if (infinit_2) {
          end2[1] = -EPS
        } else {
          end2 <- A_inv %*% (end2 - M)
        }
      }
    } else {
      if (wp > wq) {
        if (infinit_1) {
          end1[1] <- PSEUDOINF
        } else {
          end1 <- A_inv %*% (end1 - M)
        }
        if (infinit_2) {
          end2[1] <- EPS
        } else {
          end2 <- A_inv %*% (end2 - M)
        }
      } else {
        if (infinit_1) {
          end1[1] = -EPS
        } else {
          end1 = A_inv %*% (end1 - M)
        }
        if (infinit_2) {
          end2[1] = -PSEUDOINF
        } else {
          end2 <- A_inv %*% (end2 - M)
        }
      }
    }

    x <- seq(min(end1[1], end2[1]), max(end1[1], end2[1]), length.out=resolution);
    y <- 1/x

    hyperbola <- rbind(x,y)
    hyperbola <- matrix(M,2,length(x)) + A %*% hyperbola

    # all the visible plotting:
    if (show_points) {
      points(points[,1],points[,2], pch=20)
    } 
    
    if (show_weights) {
      phi <- seq(0,2*pi,length.out=300)
      r <- pmax(0, add_to_weights + weights)
      for (i in 1:dim(points)[1]) {
        polygon(points[i,1]+r[i]*cos(phi), points[i,2]+r[i]*sin(phi), border=grey(0.6))
      }
    }
    
    hyperbola[1,] <- hyperbola[1,]+xmin
    hyperbola[2,] <- hyperbola[2,]+ymin
    lines(hyperbola[1,], hyperbola[2,], col=col, lwd=lwd)
  }
}

invisible()
}
