\name{plotVIMP}
\alias{plotVIMP}
\title{
Variable Importance (VIMP) plot
}
\description{
Barplot displaying variable importance for the main effect.
}
\usage{
plotVIMP(vimp_Object,
         xvar.names = NULL,
         cex.xlab = NULL,
         ymaxlim = 0,
         yminlim = 0,
         main = "Variable Importance (\%)",
         col = grey(0.8),
         cex.lab = 1.5,
         ylbl = NULL,
         legend_placement = NULL,
         plot.it = TRUE,
         path_saveplot = NULL,
         Verbose = TRUE)
}
\arguments{
  \item{vimp_Object}{
List with number of elements equal to the number of response variables.
}
\item{xvar.names}{
Names of the covariates. If NULL, names will be pulled from \code{vimp_Object}.
}
  \item{cex.xlab}{
Magnification of the names of the covariates for the barplot.
}
  \item{ymaxlim}{
By default, we use the range of the vimp values for the barplot limit on the y-axis. If one wants to extend the limit, add the amount with which the limit will extend above the x-axis.
}
  \item{yminlim}{
Similar to \code{ymaxlim}, this will add the amount with which the limit will extend below the x-axis.
}
  \item{main}{
Main title for the plot.
}
  \item{col}{
Color of the plot.
}
  \item{cex.lab}{
Magnification of the x and y lables.
}
  \item{ylbl}{
Label for the y-axis.
}
\item{legend_placement}{
Do you want name of the covariates on top of the each barplot? If so, use
default setting; else set value on the negative direction of y-axis
which arrange covariate name beneath the barplot.
}
\item{plot.it}{Should the VIMP plot be displayed?}
\item{path_saveplot}{
Provide the location where plot should be saved. By default the plot will be
saved at temporary folder.
}
\item{Verbose}{
Display the path where the plot is saved?
}

}
\details{
Barplot displaying VIMP for each response. Barplot will be save as
pdf file in the working directory.
}
\author{
  Amol Pande and Hemant Ishwaran
}
\examples{
\donttest{
##-----------------------------------------------------------------
## VIMP plot for multivariate longitudinal response
##-----------------------------------------------------------------

# Simulate data involves 3 response and 4 covariates

dta <- simLong(n = 100, ntest = 100 ,N = 5, rho =.80, model = 1, q_x = 0, 
                                  q_y = 0,type = "corCompSym")
dtaL <- dta$dtaL
trn <- dta$trn
# Boosting call: Raw values of covariates, B-spline for time, 
# no shrinkage, no estimate of rho and phi

boost.grow <- BoostMLR(x = dtaL$features[trn,], tm = dtaL$time[trn], 
                      id = dtaL$id[trn], y = dtaL$y[trn,], M = 100, VarFlag = FALSE)

boost.pred <- predictBoostMLR(Object = boost.grow, x = dtaL$features[-trn,], 
                               tm = dtaL$time[-trn], id = dtaL$id[-trn], 
                               y = dtaL$y[-trn,], importance = TRUE)

# Plot VIMP
plotVIMP(vimp_Object = boost.pred$vimp,ymaxlim = 20,plot.it = FALSE)
}
}
\keyword{plot}
