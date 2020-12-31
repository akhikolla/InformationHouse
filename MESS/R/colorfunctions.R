#' Shade an RGB color
#'
#' Shades an RBG color
#'
#' This function shades an RGB color and returns the shaded RGB color (with alpha channel added)
#'
#' @param col  a vector of RGB color(s)
#' @param shade numeric value between 0 and 1. Zero means no change and 1 results in black
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Ekstrom, CT (2011) \emph{The R Primer}.
#' @keywords iplot
#' @examples
#'
#' newcol <- col.shade("blue")
#'
#' @importFrom grDevices col2rgb rgb
#' @export
col.shade <- function(col, shade=.5) {

    if(missing(col))
        stop("a vector of colours is missing")

    if (shade<0 | shade>1)
        stop("shade must be between 0 and 1")

    mat <- t(col2rgb(col, alpha=TRUE) * c(rep(1-shade, 3), 1))
    rgb(mat, alpha=mat[,4], maxColorValue=255)
}


#' Tint an RGB color
#'
#' Tints an RBG color
#'
#' This function tints an RGB color and returns the tinted RGB color (with alpha channel added)
#'
#' @param col  a vector of RGB color(s)
#' @param tint numeric value between 0 and 1. Zero results in white and 1 means no change
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Ekstrom, CT (2011) \emph{The R Primer}.
#' @keywords iplot
#' @examples
#'
#' newcol <- col.tint("blue")
#'
#' @export
col.tint <- function(col, tint=.5) {

    if(missing(col))
        stop("a vector of colours is missing")

    if (tint<0 | tint>1)
        stop("shade must be between 0 and 1")

    mat <- t(col2rgb(col, alpha=TRUE)  +  c(rep(1-tint, 3), 0)*(255-col2rgb(col, alpha=TRUE)))
    rgb(mat, alpha=mat[,4], maxColorValue=255)
}

#' Add and set alpha channel for RGB color
#'
#' Add and set alpha channel
#'
#' This function adds and set an alpha channel to a RGB color
#'
#' @param col a vector of RGB color(s)
#' @param alpha numeric value between 0 and 1. Zero results fully transparent and 1 means full opacity
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Ekstrom, CT (2011) \emph{The R Primer}.
#' @keywords iplot
#' @examples
#'
#' newcol <- col.alpha("blue", .5)
#'
#' @export
col.alpha <- function(col, alpha=1) {

    if(missing(col))
        stop("a vector of colours is missing")

    if (alpha<0 | alpha>1)
        stop("alpha must be between 0 and 1")

    apply(col2rgb(col, alpha=TRUE)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))

}


