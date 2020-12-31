% $Id: typeIerror.Rd 129 2014-06-16 13:31:01Z b002961 $
\name{typeIerror}
\alias{typeIerror}

\title{Probability of type I error for test in a bootstrap DEA model}

\description{Calculates the probability of a type I error for a test
   in bootstrapped DEA models. %; gives an after-method for the method
   %\code{boot.sw98} in the \pkg{FEAR} package
   }

\usage{typeIerror(shat,s) }

\arguments{

\item{shat}{The value of the statistic for which the probability of a
     type I error is to be calculated}

\item{s}{Vector with calculated values of the statistic for each of
     the \code{NREP} bootstraps; NREP is from \code{dea.boot} }

}

\details{Needs bootstrapped values of the test statistic}

\value{Returns the probability of a type I error}

\seealso{\code{boot.sw98} in \pkg{FEAR}, Paul W. Wilson (2008),
   \dQuote{FEAR 1.0: A Software Package for Frontier Efficiency
   Analysis with R,} \emph{Socio-Economic Planning Sciences} 42,
   247--254 }



\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\examples{
# Probability of getting something larger than 1.96 in 10000 random
# standard normal variates.
x <- rnorm(10000)
typeIerror(1.96,x)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
