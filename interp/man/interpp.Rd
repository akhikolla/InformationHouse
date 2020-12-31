\name{interpp}
\alias{interpp}
\title{
  Pointwise interpolate irregular gridded data
}
\description{
  This function implements bivariate interpolation onto a set of points
  for irregularly spaced input data.

  This function is meant for backward compatibility to package
  \code{akima}, please use \code{\link{interp}} with its \code{output}
  argument set to \code{"points"} now.
}
\usage{
interpp(x, y = NULL, z, xo, yo = NULL, linear = TRUE,
  extrap = FALSE, duplicate = "error", dupfun = NULL,
  deltri = "shull")
}
\arguments{

  \item{x}{
    vector of x-coordinates of data points or a
    \code{SpatialPointsDataFrame} object.
    Missing values are not accepted.
  }
  \item{y}{
    vector of y-coordinates of data points.
    Missing values are not accepted.

    If left as NULL indicates that \code{x} should be a
    \code{SpatialPointsDataFrame} and \code{z} names the variable of
    interest in this dataframe.
  }
  \item{z}{
    vector of z-coordinates of data points or a character variable
    naming the variable of interest in the
    \code{SpatialPointsDataFrame} \code{x}.

    Missing values are not accepted.

    \code{x}, \code{y}, and \code{z} must be the same length
    (execpt if \code{x} is a \code{SpatialPointsDataFrame}) and may
    contain no fewer than four points. The points of \code{x} and \code{y}
    cannot be collinear, i.e, they cannot fall on the same line (two vectors
    \code{x} and \code{y} such that \code{y = ax + b} for some \code{a},
    \code{b} will not be accepted).
  }
  \item{xo}{
    vector of x-coordinates of points at which to evaluate the interpolating
    function. If \code{x} is a \code{SpatialPointsDataFrame} this has
    also to be a \code{SpatialPointsDataFrame}.
  }
  \item{yo}{
    vector of y-coordinates of points at which to evaluate the interpolating
    function.

    If operating on \code{SpatialPointsDataFrame}s this is left as \code{NULL}
  }
  \item{linear}{logical -- indicating wether linear or spline
    interpolation should be used.
  }
  \item{extrap}{
    logical flag: should extrapolation be used outside of the
    convex hull determined by the data points? Not possible for linear
    interpolation.}
  \item{duplicate}{
    indicates how to handle duplicate data points. Possible values are
    \code{"error"} - produces an error message, \code{"strip"} - remove
    duplicate z values, \code{"mean"},\code{"median"},\code{"user"}  -
    calculate  mean , median or user defined function of duplicate z
    values.
  }
  \item{dupfun}{
    this function is applied to duplicate points if \code{duplicate="user"}
  }
  \item{deltri}{
    triangulation method used, this argument will later be moved
    into a control set together with others related to the spline
    interpolation!
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  a list with 3 components:
  \item{x,y}{
    If \code{output="grid"}:
    vectors of \eqn{x}- and \eqn{y}-coordinates of output grid, the same
    as the input
    argument \code{xo}, or \code{yo}, if present.  Otherwise, their
    default, a vector 40 points evenly spaced over the range of the
    input \code{x} and \code{y}.

    If \code{output="points"}: vectors of \eqn{x}- and \eqn{y}-coordinates of
    output points as given by \code{xo} and \code{yo}.
    }
  \item{z}{
    If \code{output="grid"}:
    matrix of fitted \eqn{z}-values.  The value \code{z[i,j]} is computed
    at the point \eqn{(xo[i], yo[j])}. \code{z} has
    dimensions \code{length(xo)} times \code{length(yo)}.

    If \code{output="points"}: a vector with the calculated z values for
    the output points as given by \code{xo} and \code{yo}.

    If the input was a \code{SpatialPointsDataFrame} a
    \code{SpatialPixelssDataFrame} is returned for \code{output="grid"}
    and a  \code{SpatialPointsDataFrame}  for \code{output="points"}.
    }
}
\references{
  Moebius, A. F. (1827) Der barymetrische Calcul. Verlag v. Johann
  Ambrosius Barth, Leipzig,
  https://books.google.at/books?id=eFPluv_UqFEC&hl=de&pg=PR1#v=onepage&q&f=false

  Franke, R., (1979). A critical comparison of some methods for interpolation
  of scattered data. Tech. Rep. NPS-53-79-003, Dept. of Mathematics, Naval
  Postgraduate School, Monterey, Calif.


}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\note{
  This is only a call wrapper meant for backward compatibility, see
  \code{\link{interp}} for more details!
}


\seealso{
  \code{\link{interp}}
}
\examples{
### Use all datasets from Franke, 1979:
### calculate z at shifted original locations.
data(franke)
for(i in 1:5)
    for(j in 1:3){
        FR <- franke.data(i,j,franke)
        IL <- with(FR, interpp(x,y,z,x+0.1,y+0.1,linear=TRUE))
        str(IL)
    }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

