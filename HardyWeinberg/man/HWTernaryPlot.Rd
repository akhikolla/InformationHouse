\name{HWTernaryPlot}
\alias{HWTernaryPlot}
\title{ Ternary plot with the Hardy-Weinberg acceptance region }
\description{
  \code{HWTernaryPlot} is a routine that draws a ternary plot for three-way genotypic compositions (AA,AB,BB), and represents
  the acceptance region for different tests for Hardy-Weinberg equilibrium (HWE) in the plot. This allows for graphical
  testing of a large set of markers (e.g. SNPs) for HWE. The (non) significance of the test
  for HWE can be inferred from the position of the marker in the ternary plot. Different statistical tests for HWE
  can be done graphically with this routine: the ordinary chisquare test, the chisquare test with continuity
  correction and the Haldane's exact test.
}
\usage{
HWTernaryPlot(X, n = NA, addmarkers = TRUE, newframe = TRUE, hwcurve = TRUE, 
vbounds = TRUE, mafbounds = FALSE, mafvalue = 0.05, axis = 0, region = 1, 
vertexlab = colnames(X), alpha = 0.05, vertex.cex = 1, pch = 19, cc = 0.5, 
markercol = "black", markerbgcol = "black", cex = 0.75, axislab = "", 
verbose = FALSE, markerlab = NULL, markerpos = NULL, mcex = 1, connect =
FALSE, curvecols = rep("black",5), signifcolour = TRUE, curtyp =
"solid", ssf = "max", pvaluetype = "selome", ...)
}
\arguments{
  \item{X}{ a matrix of \code{n} genotypic compositions or counts. If it
    is a matrix of compositions, X should have (\code{n} rows that sum
    1, and 3 columns, with the relative frequencies of AA, AB and BB
    respectively. Argument \code{n} should be supplied as well. If X is
    a matrix of raw genotypic counts, it should have 3 columns with the
    absolute counts of AA, AB and BB respectively. Argument \code{n} may
    be supplied and will be used for painting acceptance regions. If not
    supplied \code{n} is computed from the data in \code{X}. }
  \item{n}{ the samples size (for a complete composition with no missing data). }
  \item{addmarkers}{ represent markers by dots in the triangle (\code{addmarkers=TRUE}) or not \cr 
                    (\code{addmarkers=FALSE}). }
  \item{newframe}{ allows for plotting additional markers in an already existing ternary plot. Overplotting
                   is achieved by setting \code{newframe} to \code{FALSE}. Setting \code{newframe = TRUE} (default) will 
                   create a new ternary plot. }
  \item{hwcurve}{ draw the HW parabola in the plot (\code{hwcurve=TRUE)} or not (\code{hwcurve=FALSE}). }
  \item{vbounds}{ indicate the area corresponding to expected counts > 5 (\code{vbounds=TRUE}) or not 
                  (\code{vbounds=FALSE}).}
  \item{mafbounds}{ indicate the area corresponding to MAF < \code{mafvalue}. }
  \item{mafvalue}{ a critical value for the minor allele frequency (MAF).}
  \item{axis}{ draw a vertex axis \cr
   0 = no axis is drawn \cr
   1 = draw the AA axis \cr
   2 = draw the AB axis \cr
   3 = draw the BB axis \cr}
  \item{region}{ the type of acceptance region to be delimited in the triangle \cr 
   0 = no acceptance region is drawn \cr
   1 = draw the acceptance region corresponding to a Chi-square test \cr 
   2 = draw the acceptance region corresponding to a Chi-square test with continuity correction \cr
   3 = draw the acceptance region corresponding to a Chi-square test with continuity correction for D > 0 \cr
   4 = draw the acceptance region corresponding to a Chi-square test with continuity correction for D < 0 \cr 
   5 = draw the acceptance regions for all preceding tests simultaneously \cr 
   6 = draw the acceptance region corresponding to a Chi-square test with continuity correction with the upper
       limit for D > 0 and the lower limit for D < 0 \cr
   7 = draw the acceptance region corresponding to a two-sided exact test \cr
}
  \item{vertexlab}{ labels for the three vertices of the triangle }
  \item{alpha}{ significance level (0.05 by default) }
  \item{vertex.cex}{ character expansion factor for the labels of the vertices of the triangle. }
  \item{pch}{ the plotting character used to represent the markers. }
  \item{cc}{ value for the continuity correction parameter (0.5 by default). }
  \item{markercol}{ vector with colours for the marker points in the triangle. }
  \item{markerbgcol}{ vector with background colours for the marker points in the triangle. }
  \item{cex}{ expansion factor for the marker points in the triangle. }
  \item{axislab}{ a label to be put under the horizontal axis. }
  \item{verbose}{ print information on the numerically found cut-points between curves of the acceptance region and 
                  the edges of the triangle. }
  \item{markerlab}{ labels for the markers in the triangle. }
  \item{markerpos}{ positions for the marker labels in the triangle
    (1,2,3 or 4).}
  \item{mcex}{ character expansion factor for the labels of the markers in the ternary plot. }
  \item{connect}{ connect the represented markers by a line in the ternary plot. }
  \item{curvecols}{ a vector with four colour specifications for the different curves that can be used
        to delimit the HW acceptance region. E.g. \code{curvecols=c("red",
        "green","blue","black","purple")} will paint
        the Hardy-Weinberg curve red, the limits of the acceptance region for an ordinary chi-square test
        for HWE green, the limits of the acceptance region for a chi-square test with continuity correction
        when D > 0 blue and the limits of the acceptance region for a chi-square test with continuity 
        correction when D < 0 black, and the limits of the exact acceptance region purple. }
  \item{signifcolour}{ colour the marker points automatically according to the result of a signifance test 
        (green markers non-siginficant, red markers significant). 
        \code{signifcolour} only takes effect if \code{region} is set to
	1, 2 or 7.}
  \item{curtyp}{ style of the drawn curves (\code{"dashed","solid","dotted",...})}
  \item{ssf}{ sample size function (\code{"max","min","mean","median",...}). Indicates how the sample size for
              drawing acceptance regions is determined from the matrix
	      of counts.}
  \item{pvaluetype}{ method to compute p-values in an exact test
    (\code{"dost"} or \code{"selome"})}	    
  \item{\dots}{ other arguments passed on to the plot function (e.g. \code{main} for a main title). }
}
\value{
  \item{minp }{minimum allele frequency above which testing for HWE is appropriate (expected counts exceeding 5).}
  \item{maxp }{maximum allele frequency below which testing for HWE is appropriate.}
  \item{inrange }{number of markers in the appropriate range.}
  \item{percinrange }{percentage of markers in the appropriate.}
  \item{nsignif }{number of significant markers (only if \code{region} equals 1,2 or 7.)}
}
\details{
   \code{HWTernaryPlot} automatically colours significant markers in
   red, and non-significant markers in green if \code{region} is set to
   1, 2 or 7.
   }
\references{
Graffelman, J. and Morales, J. (2008) Graphical Tests for Hardy-Weinberg Equilibrium 
Based on the Ternary Plot. \emph{Human Heredity} 65(2):77-84. 

Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The HardyWeinberg Package.
\emph{Journal of Statistical Software} 64(3): 1-23. \url{https://www.jstatsoft.org/v64/i03/}.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}} }
\examples{

n <- 100 # sample size
m <- 100 # number of markers

X <- HWData(n,m)
HWTernaryPlot(X,100,region=1,hwcurve=TRUE,vbounds=FALSE,vertex.cex=2)
}

\keyword{aplot}

