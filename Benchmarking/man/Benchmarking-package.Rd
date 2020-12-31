% $Id: Benchmarking-package.Rd 232 2020-08-07 08:47:59Z lao $
\name{Benchmarking-package}
\alias{Benchmarking-package}
\alias{Benchmarking}
\docType{package}

\title{ Data Envelopment Analyses (DEA) and Stochastic Frontier
Analyses (SFA) -- Model Estimations and Efficiency Measuring }

\description{

The Benchmarking package contains methods to estimate technologies and
measure efficiencies using DEA and SFA. Data Envelopment Analysis (DEA)
are supported under different technology assumptions (fdh, vrs, drs, crs,
irs, add), and using different efficiency measures (input based, output
based, hyperbolic graph, additive, super, directional). Peers are
available, partial price information can be included, and optimal cost,
revenue and profit can be calculated. Evaluation of mergers are also 
supported. Comparative methods for estimating stochastic frontier function
(SFA) efficiencies and for convex nonparametric least squares here for
convex functions (StoNED) are also included. The methods can solve not
only standard models, but also many other model variants, and they can
be modified to solve new models.

The package also support simple plots of DEA technologies with two
goods; either as a transformation curve (2 outputs), an isoquant (2
inputs), or a production function (1 input and 1 output). When more
inputs and outputs are available they are aggregated using weights
(prices, relative prices).

The package complements the book, Bogetoft and Otto,
\emph{Benchmarking with DEA, SFA, and R}, Springer-Verlag 2011, but
can of course also be used as a stand-alone package.

}

\details{

\tabular{ll}{
Package: \tab Benchmarking\cr
Type: \tab Package\cr
Version: \tab 0.29 ($Revision: 232 $)\cr
Date: \tab $Date: 2020-08-07 10:47:59 +0200 (fr, 07 aug 2020) $\cr
License: \tab Copyright\cr
}

\tabular{ll}{
\code{dea} \tab DEA input or output efficience measures, peers,
     lambdas and slacks \cr

\code{dea.dual} \tab Dual weights (prices), including restrictions on
      weights \cr

\code{dea.direct} \tab Directional efficiency \cr

\code{sdea} \tab Super efficiency. \cr

\code{dea.add} \tab Additive efficiency; sum of slacks in DEA
     technology. \cr

\code{mea} \tab Multidirectional efficiency analysis or potential
     improvements. \cr

\code{eff} \tab Efficiency from an object returned from any of the dea
     or sfa functions. \cr

\code{slack} \tab Slacks in DEA models \cr

\code{excess} \tab Calculates excess input or output compared to DEA
     frontier. \cr

\code{peers} \tab  get the peers for each firm. \cr

\code{dea.boot} \tab Bootstrap DEA models \cr

\code{cost.opt} \tab Optimal input for given output and prices. \cr

\code{revenue.opt} \tab Optimal output for given input and prices. \cr

\code{profit.opt} \tab Optimal input and output for given input and
     output prices.  \cr

\code{dea.plot} \tab Graphs of DEA technologies under alternative
     technology assumptions. \cr

\code{dea.plot.frontier} \tab Specialized for 1 input and 1
     output. \cr

\code{dea.plot.isoquant} \tab Specialized for 2 inputs. \cr

\code{dea.plot.transform} \tab Specialized for 2 outputs. \cr

\code{eladder} \tab Efficiency ladder for a single firm.  \cr

\code{eladder.plot} \tab Plot efficiency ladder for a single firm.  \cr

\code{make.merge} \tab Make an aggregation matrix to perform mergers.  \cr

\code{dea.merge} \tab Decompose efficiency from a merger of firms \cr

\code{sfa} \tab Stochastic frontier analysis, production, distance,
     and cost function (SFA) \cr
     
\code{stoned} \tab Convex nonparametric least squares here for convex 
    function function \cr     

\code{outlierC.ap, outlier.ap} \tab Detection of outliers  \cr

\code{eff.dens} \tab Estimate and plot kernel density of efficiencies \cr

\code{critValue} \tab Critical values calculated from bootstrap DEA
     models.  \cr

\code{typeIerror} \tab Probability of a type I error for a test in
      bootstrap DEA models.  \cr

}
}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\references{

Bogetoft and Otto; \emph{Benchmarking with DEA, SFA, and R}; Springer
2011

Paul W. Wilson (2008), \dQuote{FEAR 1.0: A Software Package for
Frontier Efficiency Analysis with R,} \emph{Socio-Economic Planning
Sciences} 42, 247--254

}


\note{ 

The interface for the methods are very much like the interface to the
methods in the package \pkg{FEAR} (Wilson 2008). One change is that
the data now are transposed to reflect how data is usually available
in applications, i.e. we have firms on rows, and inputs and output in
the columns. Also, the argument for the options \code{RTS} and
\code{ORIENTATION} can be given as memotechnical strings, and there
are more options to control output.

The input and output matrices can contain negative numbers, and the
methods can thereby manage restricted or fixed input or output.

The return is not just the efficiency, but also slacks, dual values
(shadow prices), peers, and lambdas (weights). 

}



\examples{
# Plot of different technologies
x <- matrix(c(100,200,300,500),ncol=1,dimnames=list(LETTERS[1:4],"x"))
y <- matrix(c(75,100,300,400),ncol=1,dimnames=list(LETTERS[1:4],"y"))
dea.plot(x,y,RTS="vrs",ORIENTATION="in-out",txt=rownames(x))
dea.plot(x,y,RTS="drs",ORIENTATION="in-out",add=TRUE,lty="dashed",lwd=2)
dea.plot(x,y,RTS="crs",ORIENTATION="in-out",add=TRUE,lty="dotted")
                      
dea.plot(x,y,RTS="fdh",ORIENTATION="in-out",txt=rownames(x),main="fdh")
dea.plot(x,y,RTS="irs",ORIENTATION="in-out",txt=TRUE,main="irs")
dea.plot(x,y,RTS="irs2",ORIENTATION="in-out",txt=rownames(x),main="irs2")
dea.plot(x,y,RTS="add",ORIENTATION="in-out",txt=rownames(x),main="add")

#  A quick frontier with 1 input and 1 output
dea.plot(x,y, main="Basic plot of frontier")

# Calculating efficiency
dea(x,y, RTS="vrs", ORIENTATION="in")
e <- dea(x,y, RTS="vrs", ORIENTATION="in")
e
eff(e)
peers(e)
peers(e, NAMES=TRUE)
print(peers(e, NAMES=TRUE), quote=FALSE)
lambda(e)
summary(e)


# Calculating super efficiency
esuper <- sdea(x,y, RTS="vrs", ORIENTATION="in")
esuper
print(peers(esuper,NAMES=TRUE),quote=FALSE)
# Technology for super efficiency for firm number 3/C 
# Note that drop=FALSE is necessary for XREF and YREF to be matrices
# when one of the dimensions is or is reduced to 1.
e3 <- dea(x,y, XREF=x[-3,,drop=FALSE], YREF=y[-3,,drop=FALSE])
dea.plot(x[-3],y[-3],RTS="vrs",ORIENTATION="in-out",txt=LETTERS[c(1,2,4)])
points(x[3],y[3],cex=2)
text(x[3],y[3],LETTERS[3],adj=c(-.75,.75))
e3 <- dea(x,y, XREF=x[-3,,drop=FALSE], YREF=y[-3,,drop=FALSE])
eff(e3)
peers(e3)
print(peers(e3,NAMES=TRUE),quote=FALSE)
lambda(e3)
e3$lambda

# Taking care of slacks
x <- matrix(c(100,200,300,500,100,600),ncol=1,
        dimnames=list(LETTERS[1:6],"x"))
y <- matrix(c(75,100,300,400,50,400),ncol=1,
        dimnames=list(LETTERS[1:6],"y"))

# Phase one, calculate efficiency
e <- dea(x,y)
print(e)
peers(e)
lambda(e)
# Phase two, calculate slacks (maximize sum of slacks)
sl <- slack(x,y,e)
data.frame(sl$sx,sl$sy)
peers(sl)
lambda(sl)
sl$lambda
summary(sl)

# The two phases in one function call
e2 <- dea(x,y,SLACK=TRUE)
print(e2)
data.frame(eff(e2),e2$slack,e2$sx,e2$sy,lambda(e2))
peers(e2)
lambda(e2)
e2$lambda
}


\keyword{package}
\keyword{models}
\keyword{efficiency}

