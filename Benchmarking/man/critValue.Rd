% $Id: critValue.Rd 129 2014-06-16 13:31:01Z b002961 $
\name{critValue}
\alias{critValue}

\title{Critical values from bootstrapped DEA models}

\description{Calculates critical value for test using bootstrap output
in DEA models%; an after method for the method boot.sw98 in the FEAR
%package
} 

\usage{ critValue(s, alpha=0.05) }

\arguments{

\item{s}{Vector with calculated values of the statistic for each of
the \code{NREP} bootstraps; NREP is from \code{boot.sw98} }

\item{alpha}{The size of the test}

}

\details{Needs bootstrapped values of the test statistic}

\value{Returns the critical value}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\seealso{\code{boot.sw98} in \pkg{FEAR}, Paul W. Wilson (2008),
\dQuote{FEAR 1.0: A Software Package for Frontier Efficiency Analysis
with R,} \emph{Socio-Economic Planning Sciences} 42, 247--254 }


\examples{
# The critical value for two-sided test in normal distribution found
# by simulation.
x <- rnorm(1000000)
critValue(x,.975)
}

% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
