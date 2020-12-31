\name{transport}
\alias{transport}
\alias{transport.default}
\alias{transport.pgrid}
\alias{transport.pp}
\alias{transport.wpp}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Find Optimal Transport Plan Between Two Objects
}
\description{
  Given two objects \code{a} and \code{b} that specify distributions of mass and an object that specifies (a way to compute) costs,
  find the transport plan for going from \code{a} to \code{b} that minimizes the total cost.
}
\usage{
transport(a, b, ...)
\method{transport}{default}(a, b, costm, method = c("networkflow", "shortsimplex", "revsimplex",
"primaldual"), fullreturn=FALSE, control = list(), threads=1, ...)
\method{transport}{pgrid}(a, b, p = NULL, method = c("auto", "networkflow", "revsimplex", "shortsimplex",
  "shielding", "aha", "primaldual"), fullreturn=FALSE,
  control = list(), threads=1,...)
\method{transport}{pp}(a, b, p = 1, method = c("auction", "auctionbf", "networkflow", "shortsimplex",
  "revsimplex", "primaldual"),fullreturn=FALSE, 
  control = list(),threads=1, ...)
\method{transport}{wpp}(a, b, p = 1, method = c("networkflow", "revsimplex", "shortsimplex",
"primaldual"), fullreturn=FALSE ,control = list(),threads=1, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a, b}{
    two objects that describe mass distributions, between which the optimal transport map is to be computed. For the default
    method these are vectors of non-negative values. For the other three methods these are objects of the respective classes.
    It is also possible to have \code{a} of class \code{pgrid} and \code{b} of class \code{wpp}.
}
  \item{costm}{
    for the default method a \code{length(a)} by \code{length(b)} matrix specifying the cost of transporting single
    units of mass between the corresponding source and destination points. 
}
  \item{p}{
    for the three specialized methods the power \eqn{\geq 1}{>=1} to which the Euclidean distance between points is taken in order
    to compute costs.
}
  \item{method}{
    the name of the algorithm to use. See details below.
}
  \item{fullreturn}{
    A boolean specifying whether the output of the function should also include the dual solution, the optimal transport cost between a and b and the transport plan in matrix form should be returned as well.
}
  \item{control}{
  	a named list of parameters for the chosen method or the result of a call to \code{\link{trcontrol}}. Any parameters
  	that are not set by the control argument will get reasonable (sometimes problem specific) defaults.
  }
  \item{threads}{ An Integer specifying the number of threads used in parallel computing. Currently only
availible for the method "networkflow".}
  \item{...}{
  	currently without effect.
  }
}
\details{
  There is a number of algorithms that are currently implemented and more will be added in future versions of the package. 
  The following is a brief description of each key word used. Much more details can be found in the cited references 
  and in a forthcoming package vignette.

  \code{aha}: The Aurenhammer--Hoffmann--Aronov (1998) method with the multiscale approach presented in Mérigot (2011). The original theory was limited to \eqn{p=2}. We refer by \code{aha} also to the extension of the same idea for \eqn{p=1} as presented in Hartmann and Schuhmacher (2017) and for more general \eqn{p} (currently not implemented).

  \code{auction}: The auction algorithm by Bertsekas (1988) with epsilon-scaling, see Bertsekas (1992).

  \code{auctionbf}: A refined auction algorithm that combines forward and revers auction, see Bertsekas (1992).
  
  \code{networkflow}: The fast implementation of the network simplex algorithm by Nicolas Bonneel based on the LEMON Library (see citations below).
  
  \code{primaldual}: The primal-dual algorithm as described in Luenberger (2003, Section 5.9).

  \code{revsimplex}: The revised simplex algorithm as described in Luenberger and Ye (2008, Section 6.4) with various speed improvements, including a multiscale approach.

  \code{shielding}: The shielding (or shortcut) method, as described in Schmitzer (2016).
  
  \code{shortsimplex}: The shortlist method based an a revised simplex algorithm, as described in Gottschlich and Schuhmacher (2014).

  The order of the \emph{default} key words specified for the argument \code{method} gives a rough idea of the relative efficiency of the algorithms for the corresponding class of objects. For a given \code{a} and \code{b} the actual computation times may deviate significantly from this order.
  For class \code{pgrid} the default method is \code{"auto"}, which resolves to \code{"revsimplex"} if \code{p} is not 2 or the problem is very small, and to \code{"shielding"} otherwise.

  The following table gives information about the applicability of the various algorithms (or sometimes rather
  their current implementations).
  
  \tabular{lllll}{
  	                 \tab default \tab pgrid \tab  pp  \tab  wpp   \cr
  	aha (p=1 or 2!)  \tab   -     \tab   +   \tab   -   \tab   @   \cr
    auction          \tab   -     \tab   -   \tab   +   \tab   -   \cr
    auctionbf        \tab   -     \tab   -   \tab   +   \tab   -   \cr
    networkflow      \tab   +     \tab   +   \tab   +   \tab   +   \cr  
    primaldual       \tab   *     \tab   *   \tab   *   \tab   +   \cr  
    revsimplex       \tab   +     \tab   +   \tab   *   \tab   +   \cr      
    shielding (p=2!) \tab   -     \tab   +   \tab   -   \tab   -   \cr
    shortsimplex     \tab   +     \tab   +   \tab   *   \tab   +   \cr  
  } 
  where: + recommended, * applicable (may be slow), - no implementation planned or combination does not make sense; @ indicates that the aha algorithm is available in the special combination where \code{a} is a \code{pgrid} object and \code{b} is a \code{wpp} object (and \code{p} is 2). For more details on this combination see the function \code{\link{semidiscrete}}.   

  Each algorithm has certain parameters supplied by the \code{control} argument. The following table gives an overview of parameter names and 
  their applicability.
  
  \tabular{llll}{
  	                               \tab \code{start} \tab multiscale \tab individual parameters    \cr

    \code{aha} (\eqn{p=2}!)        \tab   -     \tab   +   \tab   \code{factr}, \code{maxit}  \cr  	                        
  	\code{auction}                 \tab   -     \tab   -   \tab   \code{lasteps}, \code{epsfac}   \cr
  	\code{auctionbf}               \tab   -     \tab   -   \tab   \code{lasteps}, \code{epsfac}   \cr
  	\code{networkflow}             \tab   -     \tab   -   \tab      \cr
  	\code{primaldual}              \tab   -     \tab   -   \tab      \cr
  	\code{revsimplex}              \tab   +     \tab   +   \tab      \cr  
  	\code{shielding} (\eqn{p=2}!)  \tab   -     \tab   +   \tab      \cr  
    \code{shortsimplex}            \tab   -     \tab   -   \tab   \code{slength}, \code{kfound}, \code{psearched}   \cr  
  } 
 
  \code{start} specifies the algorithm for computing a starting solution (if needed). Currently the Modified Row Minimum Rule
  (\code{start="modrowmin"}), the North-West Corner Rule (\code{start="nwcorner"}) and the method by Russell (1969) (\code{start="russell"})
  are implemented. When \code{start="auto"} (the default) the ModRowMin Rule is chosen. However,
  for \code{transport.pgrid} and \code{p} larger than 1, there are two cases where an automatic multiscale procedure is also performed, i.e. the optimal transport is first computed on coarser grids and information from these solutions is then used for the finer girds.
  This happens for 
  \code{method = "revsimplex"}, where a single coarsening at factor \code{scmult=2} is performed, and for \code{method = "shielding"}, where a number of coarsenings adapted to the dimensions of the array is performed.   
  
  For \code{p=1} and \code{method="revsimplex"}, as well as \code{p=2} and \code{method="aha"} there are multiscale versions of
  the corresponding algorithms that allows for finer control via the parameters
  \code{nscales}, \code{scmult} and \code{returncoarse}. The default value of \code{nscales=1} suppresses
  the multiscale version. For larger problems it is advisable to use the multiscale version, which currently is only implemented for
  square pgrids in two dimensions. The algorithm proceeds then by coarsening the pgrid \code{nscales-1} times, summarizing
  each time \code{scmult^2} pixels into one larger pixels, and then solving the various transport problems starting from the coarsest and
  using each previous problem to compute a starting solution to the next finer problem. If \code{returncoarse} is \code{TRUE}, the coarser
  problems and their solutions are returned as well (\code{revsimplex} only).
  
  \code{factr}, \code{maxit} are the corresponding components of the \code{control} argument in the \code{\link{optim}} L-BFGS-B method.
  
  \code{lasteps}, \code{epsfac} are parameters used for epsilon scaling in the auction algortihm. The algorithm starts with a \dQuote{transaction cost} per bid of \code{epsfac^k * lasteps} for some reasonable \code{k} generating finer and finer approximate solutions as the \code{k} counts down to zero. Note that in order for the procedure to make sense, \code{epsfac} should be larger than one (typically two- to three-digit) and in order for the final solution to be exact \code{lasteps} should be smaller than \code{1/n}, where \code{n} is the total number of points in either of the point patterns.
  \code{slength}, \code{kfound}, \code{psearched} are the shortlist length, the number of pivot candidates needed, and the percentage of
  shortlists searched, respectively.
}
\value{
A data frame with columns \code{from}, \code{to} and \code{mass} that specifies from which element of \code{a} to which element of \code{b} what amount of mass is sent in the optimal transport plan. For class \code{pgrid} elements are specified as vector indices in terms of the usual column major enumeration of the matrices \code{a$mass} and \code{b$mass}. There are \code{\link[=plot.pgrid]{plot}} methods for the classes \code{pgrid} and \code{pp}, which can plot this solution.

If \code{returncoarse} is \code{TRUE} for the \code{revsimplex} method, a list with components \code{sol} and \code{prob} giving the solutions and problems on the various scales considered. The solution on the finest scale (i.e. the output we obtain when setting \code{returncoarse} to \code{FALSE}) is in \code{sol[[1]]}.

If \code{a} is of class \code{pgrid} and \code{b} of class \code{wpp} (and \code{p=2}), an object of class \code{power_diagram} as described in the help for the function \code{\link{semidiscrete}}. The \code{\link[=plot.pgrid]{plot}} method for class \code{pgrid} can plot this solution.
}
\references{
F. Aurenhammer, F. Hoffmann and B. Aronov (1998). Minkowski-type theorems and least-squares clustering. Algorithmica 20(1), 61--76.

D. P. Bertsekas (1988). The auction algorithm: a distributed relaxation method for the assignment problem. Annals of Operations Research 14(1), 105--123.

D. P. Bertsekas (1992). Auction algorithms for network flow problems: a tutorial introduction. Computational Optimization and Applications 1, 7--66. 

N. Bonneel (2018). Fast Network Simplex for Optimal Transport. Github repository, \href{https://github.com/nbonneel/network_simplex}{nbonneel/network_simplex}.

N. Bonneel, M. van de Panne, S. Paris and W. Heidrich (2011). Displacement interpolation using Lagrangian mass transport. ACM Transactions on Graphics (SIGGRAPH ASIA 2011) 30(6).

Egervary Research Group on Combinatorial Optimization, EGRES (2014). LEMON Graph Library v1.3.1.
\href{https://lemon.cs.elte.hu/trac/lemon}{lemon.cs.elte.hu/trac/lemon}.

C. Gottschlich and D. Schuhmacher (2014). The shortlist method for fast computation of the earth mover's distance and finding optimal solutions to transportation problems. PLOS ONE 9(10), e110214. \href{http://dx.doi.org/10.1371/journal.pone.0110214}{doi:10.1371/journal.pone.0110214}

V. Hartmann and D. Schuhmacher (2017). Semi-discrete optimal transport --- the case p=1. Preprint \href{https://arxiv.org/abs/1706.07650}{arXiv:1706.07650}

D.G. Luenberger (2003). Linear and nonlinear programming, 2nd ed. Kluwer.

D.G. Luenberger and Y. Ye (2008). Linear and nonlinear programming, 3rd ed. Springer.

Q. Mérigot (2011). A multiscale approach to optimal transport. Computer Graphics Forum 30(5), 1583--1592. \href{http://dx.doi.org/10.1111/j.1467-8659.2011.02032.x}{doi:10.1111/j.1467-8659.2011.02032.x}

B. Schmitzer (2016). A sparse multiscale algorithm for dense optimal transport. J. Math. Imaging Vision 56(2), 238--259. \href{https://arxiv.org/abs/1510.05466}{https://arxiv.org/abs/1510.05466}
}
\author{
  Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
  
  Björn Bähre \email{bjobae@gmail.com} (code for \code{\link{aha}}-method for \code{p=2})
  
  Nicolas Bonneel \email{nicolas.bonneel@liris.cnrs.fr} \cr (adaption of LEMON code for fast \code{networkflow} method)
  
  Carsten Gottschlich \email{gottschlich@math.uni-goettingen.de} \cr (original java code for \code{shortlist} and \code{revsimplex} methods)
  
  Valentin Hartmann \email{valentin.hartmann@epfl.ch} (code for \code{\link[=semidiscrete1]{aha}} method for \code{p=1})
  
  Florian Heinemann \email{florian.heinemann@uni-goettingen.de} \cr (integration of \code{networkflow} method)
  
  Bernhard Schmitzer \email{schmitzer@uni-muenster.de} (code for \code{\link{shielding}}-method)
}
\section{Use of CPLEX}{The combination of the shielding-method with the CPLEX numerical solver outperforms the other algorithms by an order of magnitude for large problems (only applicable for \code{p=2} and objects of class \code{"pgrid"}). If a local installation of CPLEX is available, the transport package can be linked against it during installation. See the file src/Makevars in the source package for instructions.
}
\section{Use of CGAL}{The combination of the aha-method with \code{p=1} requires the use of CGAL (the Computational Geometry Algorithms Library) for dealing with Apollonius diagrams. If you require this functionality, install it from https://www.cgal.org/download.html and adapt the file src/Makevars of this package according to the instructions given in that file. Then re-install 'transport' from source as usual.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[=plot.pgrid]{plot}}, \code{\link{wasserstein}}
}
\examples{
#
# example for the default method
#
a <- c(100, 200, 80, 150, 50, 140, 170, 30, 10, 70)
b <- c(60, 120, 150, 110, 40, 90, 160, 120, 70, 80)
set.seed(24)
costm <- matrix(sample(1:20, 100, replace=TRUE), 10, 10)  
res <- transport(a,b,costm)

# pretty-print solution in matrix form for very small problems:
transp <- matrix(0,10,10)
transp[cbind(res$from,res$to)] <- res$mass
rownames(transp) <- paste(ifelse(nchar(a)==2," ",""),a,sep="")
colnames(transp) <- paste(ifelse(nchar(b)==2," ",""),b,sep="")
print(transp)	

	
#
# example for class 'pgrid'
#
dev.new(width=9, height=4.5)
par(mfrow=c(1,2), mai=rep(0.1,4))
image(random32a$mass, col = grey(0:200/200), axes=FALSE)
image(random32b$mass, col = grey(0:200/200), axes=FALSE)
res <- transport(random32a,random32b)
dev.new()
par(mai=rep(0,4))
plot(random32a,random32b,res,lwd=1)


#
# example for class 'pp'
#
set.seed(27)
x <- pp(matrix(runif(400),200,2))
y <- pp(matrix(runif(400),200,2))
res <- transport(x,y)
dev.new()
par(mai=rep(0.02,4))
plot(x,y,res)


#
# example for class 'wpp'
#
set.seed(30)
m <- 30
n <- 60
massx <- rexp(m)
massx <- massx/sum(massx)
massy <- rexp(n)
massy <- massy/sum(massy)
x <- wpp(matrix(runif(2*m),m,2),massx)
y <- wpp(matrix(runif(2*n),n,2),massy)
res <- transport(x,y,method="revsimplex")
plot(x,y,res)


#
# example for semidiscrete transport between class
# 'pgrid' and class 'wpp' (p=2)
#
set.seed(33)
n <- 100
massb <- rexp(n)
massb <- massb/sum(massb)*1e5
b <- wpp(matrix(runif(2*n),n,2),massb)
res <- transport(random32a,b,p=2)
plot(random32a,b,res)


#
# example for semidiscrete transport between class
# 'pgrid' and class 'wpp' (p=1)
#
if (transport:::cgal_present()) {
  set.seed(33)
  n <- 30
  massb <- rexp(n)
  massb <- massb/sum(massb)*1e5
  b <- wpp(matrix(runif(2*n),n,2),massb)
  res <- transport(random32a,b,p=1)
  plot(random32a,b,res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aha}
\keyword{auction}
\keyword{auctionbf}
\keyword{networkflow}
\keyword{primaldual}
\keyword{revsimplex}
\keyword{shielding}
\keyword{shortsimplex}
\keyword{main function}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
