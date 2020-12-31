\name{locpoly}
\alias{locpoly}
\title{
  Local polynomial fit.
}
\description{
  This function performs a local polynomial fit of up to order 3 to
  bivariate data. It returns estimated values of the regression function
  as well as estimated partial derivatives up to order 3.
}
\usage{
locpoly(x, y, z, xo = seq(min(x), max(x), length = nx), yo = seq(min(y),
 max(y), length = ny), nx = 40, ny = 40, input = "points", output = "grid",
 h = 0, kernel = "uniform", solver = "QR", degree = 3, pd = "")
}
\arguments{
  \item{x}{
    vector of \eqn{x}-coordinates of data points.

    Missing values are not accepted.
}
  \item{y}{
    vector of \eqn{y}-coordinates of data points.

    Missing values are not accepted.

}
  \item{z}{
    vector of \eqn{z}-values at data points.

    Missing values are not accepted.

    \code{x}, \code{y}, and \code{z} must be the same length
    
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
    of same lenght and are taken as possibly irregular spaced output
    points \eqn{(xo[i],yo[i])}, \eqn{i=1,...,no} with
    \code{no=length(xo)}. \code{nx} and \code{ny} are ignored in this
    case. 
  }
  \item{nx}{
    dimension of output grid in x direction
  }
  \item{ny}{
    dimension of output grid in y direction
  }

  \item{h}{
    bandwidth parameter, between 0 and 1.
    If a scalar is given it is interpreted as ratio applied to the dataset size 
    to determine a local search neighbourhood, if set to 0 a minimum useful 
    search neighbourhood is choosen (e.g. 10 points for a cubic trend function 
    to determine all 10 parameters).

    If a vector of length 2 is given both components are interpreted as
    ratio of the \eqn{x}- and \eqn{y}-range and taken as global
    bandwidth. 
  }
  \item{kernel}{
    Text value, implemented kernels are \code{uniform} (default),
    \code{epanechnikov} and \code{gaussian}.
}
\item{solver}{
  Text value, determines used solver in fastLM algorithm used by this 
  code 

  Possible values are \code{LLt}, \code{QR} (default), \code{SVD},
  \code{Eigen} and
  \code{CPivQR} (compare \code{\link[RcppEigen]{fastLm}}).
  %% FIXME: translate their integer codes to our string values!
}
  \item{degree}{
    Integer value, degree of polynomial trend, maximum allowed value is 3.
}
  \item{pd}{
    Text value, determines which partial derivative should be returned,
    possible values are  \code{""} (default, the polynomial itself),
    \code{"x"}, \code{"y"}, \code{"xx"}, \code{"xy"}, \code{"yy"}, 
    \code{"xxx"}, \code{"xxy"}, \code{"xyy"}, \code{"yyy"} or \code{"all"}.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
If \code{pd="all"}:
\item{x }{\eqn{x} coordinates}
\item{y }{\eqn{y} coordinates}
\item{z }{estimates of \eqn{z}}
\item{zx }{estimates of \eqn{dz/dx}}
\item{zy }{estimates of \eqn{dz/dy}}
\item{zxx }{estimates of \eqn{d^2z/dx^2}}
\item{zxy }{estimates of \eqn{d^2z/dxdy}}
\item{zyy }{estimates of \eqn{d^2z/dy^2}}
\item{zxxx }{estimates of \eqn{d^3z/dx^3}}
\item{zxxy }{estimates of \eqn{d^3z/dx^2dy}}
\item{zxyy }{estimates of \eqn{d^3z/dxdy^2}}
\item{zyyy }{estimates of \eqn{d^3z/dy^3}}

If \code{pd!="all"} only the elements \code{x}, \code{y} and the desired
derivative will be returned, e.g. \code{zxy} for \code{pd="xy"}. 
}
\references{
  Douglas Bates, Dirk Eddelbuettel (2013). Fast and Elegant Numerical
  Linear Algebra Using the RcppEigen Package. Journal of Statistical
  Software, 52(5), 1-24. URL http://www.jstatsoft.org/v52/i05/.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\note{
  Function \code{\link[KernSmooth]{locpoly}} of package
  \code{KernSmooth} performs a similar task for univariate data.
}

\seealso{
  \code{\link[KernSmooth]{locpoly}}, \code{\link[RcppEigen]{fastLm}}
}
\examples{

## choose a kernel
knl <- "gaussian"

## choose global and local bandwidth 
bwg <- 0.25 # *100% of x- y-range
bwl <- 0.1  # *100% of data set

## a bivariate polynomial of degree 5:
f <- function(x,y) 0.1+ 0.2*x-0.3*y+0.1*x*y+0.3*x^2*y-0.5*y^2*x+y^3*x^2+0.1*y^5

## degree of model
dg=3 

## part 1:
## regular gridded data:
ng<- 21 # x/y size of a square data grid

## build and fill the grid with the theoretical values:

xg<-seq(0,1,length=ng)
yg<-seq(0,1,length=ng)

# xg and yg as matrix matching fg
nx <- length(xg)
ny <- length(yg)
xx <- t(matrix(rep(xg,ny),nx,ny))
yy <- matrix(rep(yg,nx),ny,nx)

fg   <- outer(xg,yg,f)

## local polynomial estimate
## global bw:
ttg <- system.time(pdg <- locpoly(xg,yg,fg,
  input="grid", pd="all", h=c(bwg,bwg), solver="QR", degree=dg, kernel=knl))
## time used:
ttg

## local bw:
ttl <- system.time(pdl <- locpoly(xg,yg,fg,
  input="grid", pd="all", h=bwl, solver="QR", degree=dg, kernel=knl))
## time used:
ttl

image(pdg$x,pdg$y,pdg$z)
contour(pdl$x,pdl$y,pdl$zx,add=TRUE,lty="dotted")
contour(pdl$x,pdl$y,pdl$zy,add=TRUE,lty="dashed")
points(xx,yy,pch=".")


## part 2:
## irregular data,
## results will not be as good as with the regular 21*21=231 points.

nd<- 41 # size of data set

## random irregular data
oldseed <- set.seed(42)
x<-runif(ng)
y<-runif(ng)
set.seed(oldseed)

z <- f(x,y)

## global bw:
ttg <- system.time(pdg <- interp::locpoly(x,y,z, xg,yg, pd="all",
  h=c(bwg,bwg), solver="QR", degree=dg,kernel=knl))

ttg

## local bw:
ttl <- system.time(pdl <- interp::locpoly(x,y,z, xg,yg, pd="all",
  h=bwl, solver="QR", degree=dg,kernel=knl))

ttl

image(pdg$x,pdg$y,pdg$z)
contour(pdl$x,pdl$y,pdl$zx,add=TRUE,lty="dotted")
contour(pdl$x,pdl$y,pdl$zy,add=TRUE,lty="dashed")
points(x,y,pch=".")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ regression }% __ONLY ONE__ keyword per line
