#' Panel plot of histogram and density curve
#'
#' Prints the histogram and corresponding density curve
#'
#' This function prints a combined histogram and density curve for use with
#' the pairs function
#' @param x a numeric vector of x values
#' @param col.bar the color of the bars
#' @param ... options passed to hist
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Ekstrom, CT (2011) \emph{The R Primer}.
#' @keywords iplot
#' @examples
#'
#' pairs(~ Ozone + Temp + Wind + Solar.R, data=airquality,
#'       lower.panel=panel.smooth, diag.panel=panel.hist,
#'       upper.panel=panel.r2)
#'
#' @export
panel.hist <- function(x, col.bar="gray", ...)  {
    funArgs <- list(...)
    funArgs$col <- NULL
    ## Set user coordinates of plotting region
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    ## Do not start new plot
    par(new=TRUE)
    myhist <- function(..., col) hist(..., col=col.bar)
    myhist(x, prob=TRUE, axes=FALSE, xlab="", ylab="", main=NULL)
    lines(density(x, na.rm=TRUE))    ## Add density curve
}


#' Panel plot of R2 values for pairs
#'
#' Prints the R2 with text size depending on the size of R2
#'
#' This function is a slight modification of the panel.cor
#' function defined on the pairs help page. It calculated and
#' prints the squared correlation, R2, with text size depending
#' on the proportion of explained variation.
#'
#' @param x a numeric vector of x values
#' @param y a numeric vector of y values
#' @param digits a numeric value giving the number of digits to present
#' @param cex.cor scaling fator for the size of text
#' @param ... extra options (not used at the moment)
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @references Ekstrom, CT (2011) \emph{The R Primer}.
#' @keywords iplot
#' @examples
#'
#' pairs(~ Ozone + Temp + Wind + Solar.R, data=airquality,
#'       lower.panel=panel.smooth, upper.panel=panel.r2)
#'
#' @export
panel.r2 <- function(x, y, digits=2, cex.cor, ...) {
    ## Set user coordinates of plotting region
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use="complete.obs")**2  # Compute R^2
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * (r/2 + .5 ))
}
