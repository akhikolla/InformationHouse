\name{outer.convhull}
\title{Version of outer which operates only in a convex hull}
\usage{outer.convhull(cx,cy,px,py,FUN,duplicate="remove",...)
}
\alias{outer.convhull}
\arguments{
  \item{cx}{x cordinates of grid}
  \item{cy}{y cordinates of grid}
  \item{px}{vector of x coordinates of points}
  \item{py}{vector of y coordinates of points}
  \item{FUN}{function to be evaluated over the grid}
  \item{duplicate}{indicates what to do with duplicate
    \eqn{(px_i,py_i)} points, default \code{"remove"}.}
  \item{...}{additional arguments for \code{FUN}}
}
\description{This version of \code{outer} evaluates \code{FUN}
  only on that part of the grid \eqn{cx} times \eqn{cy} that is enclosed within
  the convex hull of the points \eqn{(px,py)}.
  
  This can be useful for spatial estimation if no extrapolation is wanted. 
}
\value{Matrix with values of \code{FUN} (\code{NA}s if outside the
  convex hull).
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}
\seealso{\code{\link{in.convex.hull}}
}
\examples{
x<-runif(20)
y<-runif(20)
z<-runif(20)
z.lm<-lm(z~x+y)
f.pred<-function(x,y)
  {predict(z.lm,data.frame(x=as.vector(x),y=as.vector(y)))}
xg<-seq(0,1,0.05)
yg<-seq(0,1,0.05)
image(xg,yg,outer.convhull(xg,yg,x,y,f.pred))
points(x,y)
}
\keyword{spatial}

