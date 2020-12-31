% $Id: dea.boot.Rd 228 2020-07-03 15:06:03Z lao $
\name{dea.boot}
\alias{dea.boot}
\alias{boot.fear}


\title{ Bootstrap DEA models }

\description{ The function \code{dea.boot} bootstrap DEA models and
returns bootstrap of Farrell efficiencies. This function is slower than the 
\code{boot.sw89} from the package \pkg{FEAR}. 
The faster function \code{boot.fear} is a wrapper for \code{boot.sw89} from 
the package \pkg{FEAR} returning results directly as Farrell measures.
}


\usage{
dea.boot(X, Y, NREP = 200, EFF = NULL, RTS = "vrs", ORIENTATION="in", 
         alpha = 0.05, XREF = NULL, YREF = NULL, FRONT.IDX=NULL, 
         EREF = NULL, DIRECT = NULL, TRANSPOSE = FALSE, 
         SHEPHARD.INPUT = TRUE, LP, CONTROL=NULL)

boot.fear(X, Y, NREP = 200, EFF = NULL, RTS = "vrs", ORIENTATION = "in", 
         alpha = 0.05, XREF = NULL, YREF = NULL, EREF = NULL)
}

\arguments{

  \item{X}{Inputs of firms to be evaluated, a K x m matrix of
      observations of K firms with m inputs (firm x input) }

  \item{Y}{Outputs of firms to be evaluated, a K x n matrix of
      observations of K firms with n outputs (firm x input). }

  \item{NREP}{Number of bootstrap replicats}

  \item{EFF}{Efficiencies for (X,Y) relative to the technology
      generated from (XREF,YREF).}

  \item{RTS}{The returns to scale assumptions as in \code{\link{dea}},
      only works for "vrs", "drs", and "crs"; more to come.}

  \item{ORIENTATION}{Input efficiency "in" (1), output
      efficiency "out" (2), and graph efficiency "graph" (3). }

  \item{alpha}{One minus the size of the confidence interval for the
      bias corrected efficiencies}

  \item{XREF}{Inputs of the firms determining the technology, defaults
      to \code{X}. }

  \item{YREF}{Outputs of the firms determining the technology,
      defaults to \code{Y}. }
      
  \item{FRONT.IDX}{Index for firms determining the technology.}

  \item{EREF}{Efficiencies for the firms in XREF, YREF.}

  \item{DIRECT}{Does not yet work and is therefore not used.}

  \item{TRANSPOSE}{Input and output matrices are K x m and K x n for
      the default value \code{TRANSPOSE=FALSE}; this is standard in R
      for statistical models.  When \code{TRANSPOSE=TRUE} data
      matrices are m x K and n x K.}
      
  \item{SHEPHARD.INPUT}{The bootstrap of the Farrell input
    efficiencies is done as a Shephard input distance function,
    the inverse Farrell input efficiency. The option is only
    relevant for input and graph directions. }

  \item{LP}{Only for debugging purposes.}
  
  \item{CONTROL}{Possible controls to \pkg{lpSolveAPI}, see the
      documentation for that package. For examples of use see the
      function \code{\link{dea}}.}

}


\details{ The details are lightly explained in Bogetoft and Otto
   (2011) Chap. 6, and with more mathematical details in Dario and
   Simar (2007) Sect. 3.4 and in Simar and Wilson (1998).

   The bootstrap at the moment does not work for any kind of
   directional efficiency.  

   The returned confidence intervals are for the bias corrected
   efficiencies; to get confidence intervals for the uncorrected
   efficiencies add the biases to both upper and lower values for the
   intervals.
   
   Under the default option \code{SHEPHARD.INPUT=TRUE} bias and bias
   corrected efficiencies are calculated for Shephard input distance
   function and then transformed to Farrell input efficiencies to
   avoid possible negative biased corrected input efficiencies. If
   this is not wanted use the option \code{SHEPHARD.INPUT=FALSE}. This
   option is only relevant for input and graph oriented directions.}



\value{ The returned values from both functions are as follows:

  \item{eff}{Efficiencies}

  \item{eff.bc}{Bias-corrected efficiencies}

  \item{bias}{An array of bootstrap bias estimates for the K firms}

  \item{conf.int}{\code{K x 2} matrix with confidence interval for the
      estimated efficiencies}

  \item{var}{An array of bootstrap variance estimates for the K firms}

  \item{boot}{ The replica bootstrap estimates of the Farrell
      efficiencies, a \code{K x NREP} matrix }

}


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
   DEA, SFA, and R}; Springer 2011.

   Cinzia Dario and L. Simar; \emph{Advanced Robust and Nonparametric
   Methods in Efficiency Analysis}.  Methodology and Applications;
   Springer 2007.

   Leopold Simar and Paul .W. Wilson (1998), \dQuote{Sensitivity
   analysis of efficiency scores: How to bootstrap in nonparametric
   frontier models}, \emph{Management Science} 44, 49--61.

   Paul W. Wilson (2008), \dQuote{FEAR 1.0: A Software Package for
   Frontier Efficiency Analysis with R,} \emph{Socio-Economic Planning
   Sciences} 42, 247--254

}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\note{ 
   %\emph{The function \code{boot.fear} only works if the package
   %\pkg{FEAR} is installed}; it does not need to be loaded.  The
   %\code{boot.sw89} from the package \pkg{FEAR} is designed for
   %Shephard distance functions, and \code{boot.fear} is a wrapper to
   %recalculate the results suitable for Farrell efficiencies.

   The function \code{dea.boot} does not depend on the FEAR package
   and can therefore be used on computers where the package FEAR is
   not available. This, however, comes with a time penalty as it takes
   around 4 times longer to run compared to using FEAR directly% or
   %indirectly by the wrapper \code{boot.fear}
   %.

   The returned bootstrap estimates from \code{FEAR::boot.sw98} %\code{boot.fear} 
   of efficiencies are sorted for each firm individually. 
   % as this is done    by \code{FEAR::boot.sw98}. 
   Unfortunately, this means that the component of replicas is not the efficiencies for
   the same bootstrap replica, but could easily be from differencts
   bootstrap replicas. This also means that this function can
   \emph{not} be used to bootstrap tests for statistical hypotheses
   where the statistics involves summing of firms efficiencies.}

\seealso{ The documentation for \code{boot.sw98} in the package
   \pkg{FEAR}. }


\examples{
x <- matrix(c(100,200,300,500,100,200,600),ncol=1)
y <- matrix(c( 75,100,300,400, 25, 50,400),ncol=1)

e <- dea(x,y)
eff(e)

dea.plot.frontier(x,y,txt=TRUE)

#  To bootstrap for real, NREP should be at least 2000. Run the
#  following lines a couple of times with nrep=100 and see how the
#  bootstrap frontier changes from one run to the next. Try the same
#  with NREP=2000 even though is does take a longer time to run,
#  especially for dea.boot.
nrep <- 5
# nrep <- 2000

# if ( "FEAR" \%in\% .packages(TRUE) )  {
##  The following only works if the package FEAR is installed; it does
##  not have to be loaded.
#  b <- boot.fear(x,y, NREP=nrep)
# } else {
  b <- dea.boot(x,y, NREP=nrep)
# }

#  bias corrected frontier
dea.plot.frontier(b$eff.bc*x, y, add=TRUE, lty="dashed")
#  outer 95% confidence interval frontier for uncorrected frontier
dea.plot.frontier((b$conf.int[,1]+b$bias)*x, y, add=TRUE, lty="dotted")


## Test of hypothesis in DEA model
# Null hypothesis is that technology is CRS and the alternative is VRS
# Bogetoft and Otto (2011) pages 183--185.
ec <- dea(x,y, RTS="crs")
Ec <- eff(ec)
ev <- dea(x,y, RTS="vrs")
Ev <- eff(ev)
# The test statistic; equation (6.1)
S <- sum(Ec)/sum(Ev)

# To calculate CRS and VRS efficiencies in the same bootstrap replicas
# we reset the random number generator before each call of the
# function dea.boot.

# To get the an initial value for the random number generating process
# we save its state (seed)
save.seed <- sample.int(1e9,1)

# The bootstrap and calculate CRS and VRS under the assumption that
# the true technology is CRS (the null hypothesis) and such that the
# results correponds to the case where CRS and VRS are calculated for
# the same reference set of firms; to make this happen we set the
# random number generator to the same state before the calls.
set.seed(save.seed)
bc <- dea.boot(x,y, nrep,, RTS="crs")
set.seed(save.seed)
bv <- dea.boot(x,y, nrep,, RTS="vrs", XREF=x,YREF=y, EREF=ec$eff)

# Calculate the statistic for each bootstrap replica
bs <- colSums(bc$boot)/colSums(bv$boot)
# The critical value for the test (default size \code{alpha} of test is 5%)
critValue(bs, alpha=.1)
S
# Accept the hypothesis at 10% level?
critValue(bs, alpha=.1) <= S

# The probability of observering a smaller value of S when the
# hypothesis is true; the p--value.
typeIerror(S, bs)
# Accept the hypothesis at size level 10%?
typeIerror(S, bs) >= .10
}


\keyword{ efficiencies }
\keyword{ bootstrap }
