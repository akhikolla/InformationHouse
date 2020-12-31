\name{interp}
\alias{interp}
\title{
  Interpolation function
}
\description{
  This function currently implements piecewise linear interpolation
  (=barycentric interpolation).
}
\usage{
interp(x, y = NULL, z, xo = seq(min(x), max(x), length = nx),
       yo = seq(min(y), max(y), length = ny),
       linear = (method == "linear"), extrap = FALSE,
       duplicate = "error", dupfun = NULL,
       nx = 40, ny = 40, input="points", output = "grid",
       method = "linear", deltri = "shull")
}
\arguments{
  \item{x}{
    vector of \eqn{x}-coordinates of data points or a
    \code{SpatialPointsDataFrame} object.
    Missing values are not accepted.
}
  \item{y}{
    vector of \eqn{y}-coordinates of data points.
    Missing values are not accepted.

    If left as NULL indicates that \code{x} should be a
    \code{SpatialPointsDataFrame} and \code{z} names the variable of
    interest in this dataframe.
}
  \item{z}{
    vector of \eqn{z}-values at data points or a character variable
    naming the variable of interest in the
    \code{SpatialPointsDataFrame} \code{x}.

    Missing values are not accepted.

    \code{x}, \code{y}, and \code{z} must be the same length
    (execpt if \code{x} is a \code{SpatialPointsDataFrame}) and may
    contain no fewer than four points. The points of \code{x} and
    \code{y} should not be collinear, i.e, they should not fall on the same line
    (two vectors \code{x} and \code{y} such that \code{y = ax + b} for
    some \code{a}, \code{b} will not produce menaningful results).


    \code{interp} is
    meant for cases in which you have \eqn{x}, \eqn{y} values
    scattered over a plane and a \eqn{z} value for each.  If, instead,
    you are trying to evaluate a mathematical function, or get a
    graphical interpretation of relationships that can be described by a
    polynomial, try \code{\link{outer}}.
  }
  \item{xo}{
    If \code{output="grid"} (default): sequence of \eqn{x} locations for
    rectangular output grid, defaults to \code{nx} points between
    \code{min(x)} and \code{max(x)}.

    If \code{output="points"}: vector of \eqn{x} locations for output points.
  }
  \item{yo}{
    If \code{output="grid"} (default): sequence of \eqn{y} locations for
    rectangular output grid, defaults to \code{ny} points between
    \code{min(y)} and \code{max(y)}.

    If \code{output="points"}: vector of \eqn{y} locations for output
    points. In this case it has to be same length as \code{xo}.
  }
  \item{input}{
    text, possible values are \code{"grid"} (not yet implemented) and
    \code{"points"} (default).

    This is used to distinguish between regular and irregular gridded data.
  }
  \item{output}{
    text, possible values are \code{"grid"} (=default) and
    \code{"points"}.

    If \code{"grid"} is choosen then \code{xo} and \code{yo} are
    interpreted as vectors spanning a rectangular grid of points
    \eqn{(xo[i],yo[j])}, \eqn{i=1,...,nx}, \eqn{j=1,...,ny}. This
    default behaviour matches how \code{akima::interp} works.

    In the case of \code{"points"}  \code{xo} and \code{yo} have to be
    of same length and are taken as possibly irregular spaced output
    points \eqn{(xo[i],yo[i])}, \eqn{i=1,...,no} with
    \code{no=length(xo)}. \code{nx} and \code{ny} are ignored in this
    case. This case is meant as replacement for the pointwise
    interpolation done by \code{akima::interpp}. If the input \code{x}
    is a \code{SpatialPointsDataFrame} and \code{output="points"} then
    \code{xo} has to be a \code{SpatialPointsDataFrame}, \code{yo} will
    be ignored.
  }
  \item{linear}{
    logical, only for backward compatibility with \code{akima::interp},
    indicates if piecewise linear interpolation or Akima splines should be
    used. Warning: in this release only \code{linear=TRUE} is
    implemented!

    Please use the new \code{method} argument instead!
  }
  \item{method}{
    text, possible methods are (currently only, more is under
    developement) \code{"linear"} (piecewise linear interpolation within
    the triangles of the Delauney triangulation, also
    referred to as barycentric interpolation based on barycentric 
    coordinates).

    This replaces the old \code{linear} argument of \code{akima::interp}.
  }
  \item{extrap}{
    logical, indicates if extrapolation outside the convex hull is
    intended, will not work for piecewise linear interpolation!
  }
  \item{duplicate}{
    character string indicating how to handle duplicate
    data points. Possible values are
    \describe{
      \item{\code{"error"}}{produces an error message,}
      \item{\code{"strip"}}{remove duplicate z values,}
      \item{\code{"mean"},\code{"median"},\code{"user"}}{calculate
	mean , median or user defined function (\code{dupfun}) of duplicate
	\eqn{z} values.}
    }
  }
  \item{dupfun}{
    a function, applied to duplicate points if
    \code{duplicate= "user"}.}
  \item{nx}{
    dimension of output grid in x direction
  }
  \item{ny}{
    dimension of output grid in y direction
  }
  \item{deltri}{
    triangulation method used, this argument will later be moved
    into a control set together with others related to the spline
    interpolation! Possible values are \code{"shull"} (default, sweep 
    hull algorithm) and \code{"deldir"} (uses package\code{deldir}).
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
    \code{SpatialPixelsDataFrame} is returned for \code{output="grid"}
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
%\note{
%%  ~~further notes~~
%}

\seealso{
  \code{\link{interpp}}
}
\examples{
### Use all datasets from Franke, 1979:
data(franke)
for(i in 1:5)
    for(j in 1:3){
        FR <- franke.data(i,j,franke)
        IL <- with(FR, interp(x,y,z,method="linear"))
        image(IL)
        contour(IL,add=TRUE)
        with(FR,points(x,y))
    }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dplot }
\keyword{ math }
