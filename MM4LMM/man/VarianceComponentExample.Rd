\name{VarianceComponentExample}
\alias{VarianceComponentExample}
\docType{data}
\title{
Variance Component Example
%%   ~~ data name/kind ... ~~
}
\description{
This corresponds to (a sample of) the dataset presented in Giraud et al. (2017).
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("VarianceComponentExample")}
\format{
  The format is:
List of 3
\describe{
    \item{\code{Data}}{'data.frame': 432 obs. of 5 variables
		\describe{
		    \item{\code{Trial}}{a factor with 2 levels}
		    \item{\code{CodeHybrid}}{a factor with 177 levels}
		    \item{\code{CodeDent}}{a factor with 116 levels}
		    \item{\code{CodeFlint}}{a factor with 122 levels}
		    \item{\code{Trait}}{a numeric vector}
		}
  }
\item{\code{KinshipD}}{num [1:116,1:116]}
\item{\code{KinshipF}}{num [1:122,1:122]}
}
}
\details{
The list contains three elements:
\itemize{
\item Data: a data frame containing the information about hybrids (trials, hybrid names, dent parental lines, flint parental lines and phenotypes)
\item KinshipD: a matrix of simple relatedness coefficients between dent lines
\item KinshipF: a matrix of simple relatedness coefficients between flint lines
}
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{\url{http://www.genetics.org/content/207/3/1167.figures-only}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Giraud, H., Bauland, C., Falque, M., Madur, D., Combes, V., Jamin, P., ... & Blanchard, P. (2017). Reciprocal Genetics: Identifying QTL for General and Specific Combining Abilities in Hybrids Between Multiparental Populations from Two Maize (Zea mays L.) Heterotic Groups. Genetics, 207(3), 1167-1180.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(VarianceComponentExample)
names(VarianceComponentExample)
## maybe str(VarianceComponentExample) ; plot(VarianceComponentExample) ...
}
\keyword{datasets}
