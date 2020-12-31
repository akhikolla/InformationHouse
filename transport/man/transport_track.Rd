\name{transport_track}
\alias{transport_track}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a Dynamic Visualization of a Transference Plan Between Two pgrids
}
\description{
  Given two objects \code{source} and \code{target} of class \code{pgrid} and a transference plan, typically
  the result of a call to \code{\link{transport}}, create an animation of the dynamic transference plan
  (a.k.a. displacement interpolation)
}
\usage{
  transport_track(source, target, tplan, K = 50, scmult = 1, smooth = FALSE,
     H = matrix(c(1,0,0,1),2,2), create.file = c("none","gif_im"), 
     file.name = "Rtransport.gif", fps = 20, cut = FALSE,
     col=grey((0:1000)/1000),width=800,height=800)
}
\arguments{
  \item{source, target}{
    objects of class \code{pgrid}.
}
  \item{tplan}{
    a transference plan between \code{source} and \code{target}, typically an optimal transference plan
    obtained by a call to \code{\link{transport}}.
  }
    \item{K}{
    the number of intermediate frames to be produced between \code{source} and \code{target}. 
  }
  \item{scmult}{
  	the factor by which the number of pixels in each dimension is multiplied to obtain a smoother rendering 
  	of the dynamic transference plan.
  }
  \item{smooth}{
    logical. Whether a kernel smoothing or a linear binning procedure is used to generate the images. Defaults to \code{FALSE}.
  }
  \item{H}{
  	the bandwith matrix used to perform the two dimensional kernel density estimation or the linear binning respectively.
  }
  \item{create.file}{
    the file type to be created or \code{"none"} to return only an array of intermediate mass distributions.
  }
  \item{file.name}{
    the path for the output file. Ignored if \code{create.file} is \code{"none"}.
  }
  \item{fps}{
  	the number of frames per second in the generated gif. The default is 20 frames per second.
  }
  \item{cut}{
  	logical. Whether the boundary pixels are cut off. Currently the only way to deal with the edge effect (see Details).
  }
  \item{col}{
    the vector of RGB colours which is used to generate the gif, if create.file is not "none". See the 
    documentation of image for more details.
  }
  \item{width}{
    interger specifying the width of the images used to generate the output gif, if create.file is not "none".
  }
  \item{height}{
    interger specifying the width of the images used to generate the output gif, if create.file is not "none".
  }
}
\details{
  The intermediate frames are produced by the interpolation formula \eqn{[(1-t)\mathrm{pr}_1 + t\mathrm{pr}_2]_{\#} \pi}{[(1-t)*pr_1 + t*pr_2]_# \pi}, 
  where \eqn{\pi} is the transference plan, \eqn{\mathrm{pr}_1}{pr_1} and \eqn{\mathrm{pr}_2}{pr_2} are the first and second coordinate projections of \eqn{\mathbf{R}^2 \times \mathbf{R}^2}{R^2 x R^2} onto \eqn{\mathbf{R}^2}{R^2}, and \eqn{t \in \{0,1/(K+1), \ldots, K/(K+1),1\}}{t in {0, 1/(K+1), \ldots, K/(K+1), 1}}. If \eqn{\pi} is an optimal transference plan this yields the displacement interpolation, at least if we assume 
  as underlying cost function the Euclidean metric to the \eqn{p}-th power, where \eqn{p = 1,2}.
	
  The kernel smoothing procedure gives usually nicer animations, but takes several orders of magnitudes longer. 

  There are currently visible edge effects in both the kernel smoothing and the linear binning procedure that lead to darker pixels at the boundary of the image. The cut parameter may be used to remove the boundary pixels completely and thus produce a smaller output. The edge will be dealt with more adequatly in future versions.
  
  Conversion to an animated gif is performed by a system call to the convert tool of ImageMagick. The latter may have to be 
  installed first.
}
\section{Warning}{Running this function with \code{smooth=TRUE} and even moderate \code{K} can take a long time!}
\value{
  An array containing the various interpolation images.

  Unless \code{create.file="none"}, the function is mainly used for its side effect (saving a file to the specified path).
  So the array is returned invisibly.
}

%\references{
%  [Displacement interpolation?].
%}
\author{
  Florian Heinemann \email{florian.heinemann@uni-goettingen.de} \cr
  (slightly modified by Dominic Schuhmacher \email{dschuhm1@uni-goettingen.de})
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{Function \link{transport}} for computing optimal transference plans.
}
\examples{
tplan <- transport(random32a,random32b)
series <- transport_track(random32a, random32b, tplan, scmult=3, create.file="none")
dev.new(width=16,height=8)
oldpar <- par(mfrow=c(5,10), mai=rep(0.01,4))
for (i in 1:50) { 
  image(series[,,i], col=grey(seq(0,1,0.005)), asp=1, axes=FALSE,zlim=c(min(series),max(series)))
}
par(oldpar)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
