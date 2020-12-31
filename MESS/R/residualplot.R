#' @rdname residualplot
#' @export
residualplot.default <- function(x, y=NULL, candy=TRUE, bandwidth = 0.3, xlab="Fitted values", ylab="Std.res.", col.sd="blue", col.alpha=0.3, ylim=NA, ...) {

    if (is.null(y))
        stop("y must be specified")

    if (!candy) {
        plot(x, y, xlab = xlab, ylab=ylab, ...)
    } else {

        ## Set the colors
        if (col.alpha == 0)
            col.trans <- col.sd
        else col.trans <- sapply(col.sd, FUN = function(x) do.call(rgb,
                                                                   as.list(c(col2rgb(x)/255, col.alpha))))

        ## Compute the range of the interval
        uniqx <- sort(unique(x))
        vary <- length(uniqx)
        window <- bandwidth * (max(x) - min(x))/2
        for (i in 1:length(uniqx)) {
            vary[i] <- 1.96 * sd(y[abs(x - uniqx[i]) <= window])
        }
        vary[is.na(vary)] <- 0


        if (any(is.na(ylim))) {
            ylim <- c(min(y, -vary), max(y, vary))
        }
        
        plot(x, y, xlab = xlab, ylab = ylab, pch = 1 + 15 * (abs(y) > 1.96), ylim=ylim, ...)
               
        if (length(uniqx) > 3) {
            lines(smooth.spline(x, y, df = 3), lty = 2, lwd = 2.5,
                  col = "black")
        }

        polygon(c(uniqx, rev(uniqx)), c(vary, -(rev(vary))),
                col = col.trans, border = NA)
    }
    return(invisible(NULL))
}


#' @rdname residualplot
#' @export
residualplot.lm <- function(x, y, candy=TRUE, bandwidth = 0.3, xlab="Fitted values", ylab="Stud.res.", col.sd="blue", col.alpha=0.3,...) {
  y <- rstudent(x)
  x <- predict(x)
  residualplot(x, y, candy, bandwidth, xlab, ylab, col.sd, col.alpha, ...)
}


#' @rdname residualplot
#' @export
residualplot.glm <- function(x, y, candy=TRUE, bandwidth = 0.4, xlab="Fitted values", ylab="Std. dev. res.", col.sd="blue", col.alpha=0.3,...) {
    
#    family <- family(x)$family
#    if (!(family %in% c("binomial", "poisson", "gaussian")))
#        stop(paste0("Residualplot for family ", family, " in glm is not implemented yet"))

#    y <- rstudent(x, type="response")
    y <- rstandard(x) # Deviance residuals
    x <- predict(x, type="link")
    residualplot(x, y, candy, bandwidth, xlab, ylab, col.sd, col.alpha, ...)
}





#' Plots a standardized residual
#'
#' Plots a standardized residual plot from an lm or glm object and provides additional
#' graphics to help evaluate the variance homogeneity and mean.
#'
#' The y axis shows the studentized residuals (for lm objects) or
#' standardized deviance residuals (for glm objects). The x axis shows the linear predictor, i.e., the
#' predicted values for lm objects.
#' 
#' The blue area is a smoothed estimate of 1.96*SD of the standardized
#' residuals in a window around the predicted value. The blue area should
#' largely be rectangular if the standardized residuals have more or less the
#' same variance.
#'
#' The dashed line shows the smoothed mean of the standardized residuals and
#' should generally follow the horizontal line through (0,0).
#'
#' Solid circles correspond to standardized residuals outside the range from [-1.96; 1.96] while open circles are inside that interval. Roughly 5% of the observations should be outside the interval and the points should be evenly distributed.
#'
#' @aliases residualplot residualplot.lm residualplot.default
#' @param x lm object or a numeric vector
#' @param y numeric vector for the y axis values
#' @param candy logical. Should a lowess curve and local standard deviation of
#' the residual be added to the plot. Defaults to \code{TRUE}
#' @param bandwidth The width of the window used to calculate the local
#' smoothed version of the mean and the variance. Value should be between 0 and
#' 1 and determines the percentage of the window width used
#' @param xlab x axis label
#' @param ylab y axis label
#' @param col.sd color for the background residual deviation
#' @param col.alpha number between 0 and 1 determining the transprency of the
#' standard deviation plotting color
#' @param ylim pair of observations that set the minimum and maximum of the y axis. If set to NA (the default) then the limits are computed from the data. 
#' @param ... Other arguments passed to the plot function
#' @return Produces a standardized residual plot
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @seealso \code{\link{rstandard}}, \code{\link{predict}}
#' @keywords hplot
#' @examples
#'
#' # Linear regression example
#' data(trees)
#' model <- lm(Volume ~ Girth + Height, data=trees)
#' residualplot(model)
#' model2 <- lm(Volume ~ Girth + I(Girth^2) + Height, data=trees)
#' residualplot(model2)
#'
#' @export
residualplot <- function(x, y=NULL, candy=TRUE, bandwidth = 0.3,
	                 xlab="Fitted values", ylab="Std.res.",
                         col.sd="blue", col.alpha=0.3, ylim=NA, ...) {
  UseMethod("residualplot")
}







#' Plots a standardized residual
#'
#' Plots a standardized residual plot from an lm or glm object and provides additional
#' graphics to help evaluate the variance homogeneity and mean.
#'
#' The y axis shows the studentized residuals (for lm objects) or
#' standardized deviance residuals (for glm objects). The x axis shows the linear predictor, i.e., the
#' predicted values for lm objects.
#' 
#' The blue area is a smoothed estimate of 1.96*SD of the standardized
#' residuals in a window around the predicted value. The blue area should
#' largely be rectangular if the standardized residuals have more or less the
#' same variance.
#'
#' The dashed line shows the smoothed mean of the standardized residuals and
#' should generally follow the horizontal line through (0,0).
#'
#' Solid circles correspond to standardized residuals outside the range from [-1.96; 1.96] while open circles are inside that interval. Roughly 5% of the observations should be outside the interval and the points should be evenly distributed.
#'
#' @aliases residual_plot residual_plot.lm residual_plot.default
#' @param x lm object or a numeric vector
#' @param y numeric vector for the y axis values
#' @param candy logical. Should a lowess curve and local standard deviation of
#' the residual be added to the plot. Defaults to \code{TRUE}
#' @param bandwidth The width of the window used to calculate the local
#' smoothed version of the mean and the variance. Value should be between 0 and
#' 1 and determines the percentage of the window width used
#' @param xlab x axis label
#' @param ylab y axis label
#' @param col.sd color for the background residual deviation
#' @param alpha number between 0 and 1 determining the transprency of the
#' standard deviation plotting color
#' @param ylim pair of observations that set the minimum and maximum of the y axis. If set to NA (the default) then the limits are computed from the data. 
#' @param ... Other arguments passed to the plot function
#' @return Produces a standardized residual plot
#' @author Claus Ekstrom <claus@@rprimer.dk>
#' @seealso \code{\link{rstandard}}, \code{\link{predict}}
#' @keywords hplot
#' @examples
#'
#' # Linear regression example
#' data(trees)
#' model <- lm(Volume ~ Girth + Height, data=trees)
#' residual_plot(model)
#' model2 <- lm(Volume ~ Girth + I(Girth^2) + Height, data=trees)
#' residual_plot(model2)
#'
#' # Add extra information about points by adding geom_text to the object produced
#'
#' m <- lm(mpg ~ hp + factor(vs), data=mtcars)
#' residual_plot(m) + ggplot2::geom_point(ggplot2::aes(color=factor(cyl)), data=mtcars) 
#'
#' @import ggplot2 ggformula
#' @export
residual_plot <- function(x, y=NULL, candy=TRUE, bandwidth = 0.3,
	                 xlab="Fitted values", ylab="Std.res.",
                         col.sd="blue", alpha=0.1, ylim=NA, ...) {
  UseMethod("residual_plot")
}


#' @rdname residual_plot
#' @export
residual_plot.default <- function(x, y=NULL, candy=TRUE, bandwidth = 0.3, 
     xlab="Fitted values", ylab="Std.res.", col.sd="blue", alpha=0.1, ylim=NA, ...) {

    if (is.null(y))
        stop("y must be specified")

    outside <- NULL
    ymin <- NULL
    ymax <- NULL
    dummyDF <- data.frame(x, y, outside=abs(y)>1.96)
  
    # Make the default plot 
    p <- ggplot(data=dummyDF, aes(x=x, y=y)) + xlab(xlab) + ylab(ylab) 


if (candy) {

    p <- p + geom_point(alpha=alpha) + ggformula::geom_spline(df=3)

    ## Compute the range of the interval
    uniqx <- sort(unique(x))
    vary <- length(uniqx)
    window <- bandwidth * (max(x) - min(x))/2

    vary <- sapply(1:length(uniqx), function(i) { 1.96 * sd(y[abs(x - uniqx[i]) <= window]) })
    vary[is.na(vary)] <- 0

    if (any(is.na(ylim))) {
       ylim <- c(min(y, -vary), max(y, vary))
    }
 
    DF2 <- data.frame(x=uniqx, y=uniqx, ymin=-vary, ymax=vary)

    p <- p + geom_ribbon(data=DF2, aes(x=x, ymin=ymin, ymax=ymax), fill="blue", alpha=alpha) 

    # Add outlier colouring if any

    if (any(dummyDF$outside)) {
        p <- p + geom_point(data=dummyDF[dummyDF$outside,]) 
    }

}

    return(p + geom_point(alpha=alpha))
}


#' @rdname residual_plot
#' @export
residual_plot.lm <- function(x, y, candy=TRUE, bandwidth = 0.3, xlab="Fitted values", ylab="Stud.res.", col.sd="blue", alpha=0.1,...) {
  y <- rstudent(x)
  x <- predict(x)
  residual_plot(x, y, candy, bandwidth, xlab, ylab, col.sd, alpha=alpha, ...)
}


#' @rdname residual_plot
#' @export
residual_plot.glm <- function(x, y, candy=TRUE, bandwidth = 0.4, xlab="Fitted values", ylab="Std. dev. res.", col.sd="blue", alpha=0.1,...) {
    
    y <- rstandard(x) # Deviance residuals
    x <- predict(x, type="link")
    residual_plot(x, y, candy, bandwidth, xlab, ylab, col.sd, alpha=alpha, ...)
}

