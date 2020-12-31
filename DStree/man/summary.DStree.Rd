\name{summary.DStree}
\alias{summary.DStree}
\title{
  Summarize a Fitted DStree Object 
}
\description{
  Returns a detailed listing of a fitted \code{DStree} object.
}
\usage{
\method{summary}{DStree}(object, cp = 0, digits = getOption("digits"), file, \dots)
}
\arguments{
\item{object}{
  fitted model object of class \code{"DStree"}.  This is assumed to be the result
  of some function that produces an object with the same named
  components as that returned by the \code{DStree} function.
}
\item{digits}{
    Number of significant digits to be used in the result.
}
\item{cp}{
  trim nodes with a complexity of less than \code{cp} from the listing.
}
\item{file}{
  write the output to a given file name. (Full listings of a tree are
  often quite long).
}
\item{\dots}{
  arguments to be passed to or from other methods.
}}
\details{
  The function prints the call, the table shown by \code{\link{printcp}}, the
  variable importance (summing to 100) and details for each node.
}
\seealso{
  \code{\link{summary}}, \code{\link{DStree.object}}, \code{\link{printcp}}.
}
\examples{
data(UnempDur)
tree <- DStree(spell~ui+age+tenure+logwage,status="censor1",data=UnempDur)
summary(tree)

}
\keyword{tree}
