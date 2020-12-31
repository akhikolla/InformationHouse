% $Id: sfa.Rd 231 2020-07-26 20:56:26Z lao $
\name{sfa}
\alias{sfa}
\alias{sfa.cost}
\alias{lambda.sfa}
\alias{logLik.sfa}
\alias{print.sfa}
\alias{residuals.sfa}
\alias{sigma2.sfa}
\alias{sigma2u.sfa}
\alias{sigma2v.sfa}
\alias{summary.sfa}
\alias{coef.sfa}
\alias{te.sfa}
% \alias{te1.sfa}
% \alias{te2.sfa}
\alias{teJ.sfa}
\alias{teBC.sfa}
\alias{teMode.sfa}
\alias{te.add.sfa}

% \alias{efficiencies.sfa}
% \alias{eff.sfa}


\title{Stochastic frontier estimation}

\description{
Estimate a stochastic frontier production or cost function using a maximum 
likelihood method.
}

\usage{
sfa(x, y, beta0 = NULL, lambda0 = 1, resfun = ebeta, 
    TRANSPOSE = FALSE, DEBUG=FALSE,
    control=list(), hessian=2)

sfa.cost(W, Y, COST, beta0 = NULL, lambda0 = 1, resfun = ebeta, 
    TRANSPOSE = FALSE, DEBUG=FALSE,
    control=list(), hessian=2)

% \method{efficiencies}{sfa}(object)
% \method{eff}{sfa}(object)
te.sfa(object)
teBC.sfa(object)
teMode.sfa(object)
teJ.sfa(object)

te.add.sfa(object)

sigma2u.sfa(object)
sigma2v.sfa(object)
sigma2.sfa(object)

lambda.sfa(object)

}


\arguments{

  \item{x}{Input as a K x m matrix of observations on m inputs from K
    firms; (firm x input); MUST be a matrix. No constant for the
    intercept should be included in x as it is added by default.}

  \item{y}{Output; K times 1 matrix (one output)  }

  \item{Y}{Output; K times n matrix for m outputs; only to be used in cost function estimation.  }
  
  \item{W}{Input prices as a K x m matrix.}
  
  \item{COST}{Cost as a K array for the K firms}

  \item{beta0}{Optional initial parameter values}

%  \item{sigma0}{Optional initial estimate of total variance}

  \item{lambda0}{Optional initial ratio of variances}

  \item{resfun}{Function to calculate the residuals, default is a
  linear model with an intercept. Must be called as
  \code{resfun(x,y,parm)} where \code{parm=c(beta,lambda)} or
  \code{parm=c(beta)}, and return the residuals as an array of length
  corresponding to the length of output \code{y}.}

%   \item{t.value}{Calculates the t-values and write the results with
%   the parameters; also calculates standard errors of the parameters}
% 
%   \item{std.err}{Calculates the standard errors of the parameter
%   estimates and write the results along the parameters}
% 
%   \item{vcov}{Calculates the variance-covarians matrix of the parameters}

  \item{TRANSPOSE}{If TRUE, data is transposed, i.e. input is now m x K matrix}

  \item{DEBUG}{Set to TRUE to get various debugging information written
  on the console}

  \item{control}{List of control parameters to \code{ucminf}}

  \item{hessian}{How the Hessian is delivered, see the ucminf
      documentation}


  \item{object}{Object of class \sQuote{sfa} as output from the
      function \code{sfa}} }


\details{ The optimization is done by the R method \code{ucminf} from
   the package with the same name. The efficiency terms are assumed to
   be half--normal distributed.
   
   Changing the maximum steplength, the trust region,  might be important,
    and this can be done by the option 'control = list(stepmax=0.1)'. The
    default value is 0.1 and that value is suitable for parameters
    around 1; for smaller parameters a lower value should be used.
    Notice that the step length is updated by the optimizing program and
    thus must be set for every call of the function sfa if it is to be
    set.

The generic functions \code{print.sfa}, \code{summary.sfa},
   \code{fitted.sfa}, \code{residuals.sfa}, \code{logLik.sfa}, and
   \code{coef.sfa} all work as expected.

The methods \code{te.sfa}, \code{teMode.sfa} etc. calculates the
   efficiency corresponding to different methods

}


\value{

The values returned from \code{sfa} is the same as for \code{ucminf},
i.e. a list with components plus some especially relevant for sfa:

\item{par}{The best set of parameters found \code{c(beta,lambda)}.}

\item{value}{The value of minus log-likelihood function corresponding to 'par'.}

\item{beta}{The parameters for the function}

\item{sigma2}{The estimate of the total variance}

\item{lambda}{The estimate of lambda}

\item{N}{The number of observations}

\item{df}{The degrees of freedom for the model}

\item{residuals}{The residuals as a K times 1 matrix/vector, 
can also be obtained by \cr \code{residuals(sfa-object)}}

\item{fitted.values}{Fitted values}

\item{vcov}{The variance-covarians matrix for all estimated parameters
 incl. lambda}

\item{convergence}{An integer code. '0' indicates successful
          convergence.  Some of the error codes taken from
          \code{ucminf} are

'1' Stopped by small gradient (grtol). 

'2' Stopped by small step (xtol). 

'3' Stopped by function evaluation limit (maxeval). 

'4' Stopped by zero step from line search 

More codes are found in \code{\link{ucminf}}
}

\item{message}{A character string giving any additional information
    returned by the optimizer, or 'NULL'.}


\item{o}{The object returned by \code{ucminf}, for further information
on this see \code{ucminf}}

}



\references{

Bogetoft and Otto; \emph{Benchmarking with DEA, SFA, and R}, Springer
2011; chapters 7 and 8.

}


\seealso{ See the method \code{ucminf} for the possible optimization
   methods and further options to use in the option \code{control}.

   The method \code{sfa} in the package \pkg{frontier} gives another
   way to estimate stochastic production functions.

}


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\note{ Calculation of technical efficiencies for each unit can be done
     by the method te.sfa as shown in the examples.

\code{te.sfa(sfaObject), teBC.sfa(sfaObject)}: Efficiencies estimated
     by minimizing the mean square error; Eq. (7.21) in Bogetoft and
     Otto (2011, 219) and Battese and Coelli (1988, 392)

\code{teMode.sfa(sfaObject), te1.sfa(sfaObject)}: Efficiencies
     estimates using the conditional mode approach; Bogetoft and Otto
     (2011, 219), Jondrow et al. (1982, 235).

\code{teJ.sfa(sfaObject), te2.sfa(sfaObject)}: Efficiencies estimates
     using the conditional mean approach Jondrow et al. (1982, 235).

\code{te.add.sfa(sfaObject)} Efficiency in the additive model,
      Bogetoft and Otto (2011, 219)

The variance pf the distribution of efficiency can be calculated by
      \code{sigma2u.sfa(sfaObject)}, the variance of the random
      error by \code{sigma2v.sfa(sfaObject)}, and the total variance
      (sum of variances of efficiency and random noise) by
      \code{sigma2.sfa}.

The ratio of variances of the efficiency and the random noise can be
      found from the method \code{lambda.sfa}

The generic method \code{summary} prints the parameters, standard
      errors, t-values, and a few more statistics from the
      optimization.

}



\examples{
# Example from the book by Coelli et al.
# d <- read.csv("c:/0work/rpack/front41Data.csv", header = TRUE, sep = ",")
# x <- cbind(log(d$capital), log(d$labour))
# y <- matrix(log(d$output))

n <- 50
x1 <- 1:50 + rnorm(n, 0, 10)
x2 <- 100 + rnorm(n, 0, 10)
x <- cbind(x1, x2)
y <- 0.5 + 1.5*x1 + 2*x2 + rnorm(n, 0, 1) - pmax(0, rnorm(n, 0, 1))
sfa(x,y)
summary(sfa(x,y))


# Estimate efficiency for each unit
o <- sfa(x,y)
eff(o)

te <- te.sfa(o)
teM <- teMode.sfa(o)
teJ <- teJ.sfa(o)
cbind(eff(o),te,Mode=eff(o, type="Mode"),teM,teJ)[1:10,]


sigma2.sfa(o)       # Estimated varians
lambda.sfa(o)       # Estimated lambda
}

\keyword{sfa}
\keyword{efficiency}%
