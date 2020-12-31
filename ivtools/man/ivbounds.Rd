\name{ivbounds}

\alias{ivbounds}

\title{
Bounds for counterfactual outcome probabilities in instrumental variables scenarios 
}

\description{
\code{ivbounds} computes non-parametric bounds for counterfactual outcome probabilities
in instrumental variables scenarios. Let \eqn{Y}, \eqn{X}, and \eqn{Z} 
be the outcome, exposure, and instrument, respectively. \eqn{Y} and \eqn{X} must be binary, 
whereas \eqn{Z} can be either binary or ternary. 
Ternary instruments are common in, for instance, Mendelian randomization. 
Let \eqn{p(Y_x=1)} be the counterfactual probability of the outcome, had all 
subjects been exposed to level \eqn{x}. \code{ivbounds} computes bounds for the 
counterfactuals probabilities \eqn{p(Y_1=1)} and \eqn{p(Y_0=1)}. Below, we define
\eqn{p_{yx.z}=p(Y=y,X=x|Z=x)}. 
}  

\usage{
ivbounds(data, Z, X, Y, monotonicity=FALSE, weights)
}

\arguments{
  \item{data}{
either a data frame containing the variables in the model, or a named vector  
\code{(p00.0,...,p11.1)} when \eqn{Z} is binary, or 
a named vector \code{(p00.0,...,p11.2)} when \eqn{Z} is ternary.
}
  \item{Z}{
a string containing the name of the instrument \eqn{Z} in \code{data} if
\code{data} is a data frame. In this case \eqn{Z} has to be coded as (0,1) when
binary, and coded as (0,1,2) when ternary. \code{Z} is not specified if 
\code{data} is a vector of probabilities.  
}
  \item{X}{
a string containing the name of the exposure \eqn{X} in \code{data} if
\code{data} is a data frame. In this case \eqn{X} has to be coded as (0,1). 
\code{X} is not specified if \code{data} is a vector of probabilities. 
}
  \item{Y}{
a string containing the name of the outcome \eqn{Y} in \code{data} if
\code{data} is a data frame. In this case \eqn{Y} has to be coded as (0,1). 
\code{Y} is not specified if \code{data} is a vector of probabilities. 
}
\item{monotonicity}{
logical. It is sometimes realistic to make the monotonicity assumption 
\eqn{z \geq z' \Rightarrow X_z \geq X_{z'}}. Should the bounds be computed 
under this assumption?
 
}

  \item{weights}{
an optional vector of `prior weights' to be used in the fitting process. 
Should be NULL or a numeric vector. Only applicable if \code{data} is 
a data frame.  
}

}

\details{
\code{ivbounds} uses linear programming techniques to bound the counterfactual probabilities
\eqn{p(Y_1=1)} and \eqn{p(Y_0=1)}. Bounds for a causal effect, defined as a contrast between these,
are obtained by plugging in the bounds for \eqn{p(Y_1=1)} and \eqn{p(Y_0=1)} into the 
contrast. For instance, bounds for the causal risk difference \eqn{p(Y_1=1)-p(Y_0=1)}
are obtained as \eqn{[min\{p(Y_1=1)\}-max\{p(Y_0=1)\},max\{p(Y_1=1)\}-min\{p(Y_0=1)\}]}.
In addition to the bounds, \code{ivbounds} evaluates the IV inequality
\deqn{\max\limits_{x}\sum_{y}\max\limits_{z}p_{yx.z}\leq 1.} 
}

\value{
An object of class \code{"ivbounds"} is a list containing 
\item{call}{
  the matched call. 
  }
\item{p0}{
  a named vector with elements \code{"min"} and \code{"max"}, containing the 
  evaluated lower and upper bounds for \eqn{p(Y_0=1)}, respectively. 
  }
\item{p1}{
  a named vector with elements \code{"min"} and \code{"max"}, containing the 
  evaluated lower and upper bounds for \eqn{p(Y_1=1)}, respectively. 
  }
\item{p0.symbolic}{
  a named vector with elements \code{"min"} and \code{"max"}, containing the 
  lower and upper bounds for \eqn{p(Y_0=1)}, respectively, on a symbolic form
  (i.e. as strings). 
  }
\item{p1.symbolic}{
  a named vector with elements \code{"min"} and \code{"max"}, containing the 
  lower and upper bounds for \eqn{p(Y_1=1)}, respectively, on a symbolic form
  (i.e. as strings). 
  }
  \item{IVinequality}{
  logical. Does the IV inequality hold?  
  }
\item{conditions}{
  a character vector containing the violated condiations, if \code{IVinequality=FALSE}.  
  }
}

\references{
Balke, A. and Pearl, J. (1997). \emph{Bounds on treatment effects from studies with imperfect compliance}.
\emph{Journal of the American Statistical Association} \bold{92}(439), 1171-1176. 

Sjolander A., Martinussen T. (2019). Instrumental variable estimation with the R package ivtools.
\emph{Epidemiologic Methods} \bold{8}(1), 1-20.

}

\author{
Arvid Sjolander. 
}

\examples{

##Vitamin A example from Balke and Pearl (1997).
n000 <- 74
n001 <- 34
n010 <- 0
n011 <- 12
n100 <- 11514
n101 <- 2385
n110 <- 0
n111 <- 9663
n0 <- n000+n010+n100+n110
n1 <- n001+n011+n101+n111

#with data frame...
data <- data.frame(Y=c(0,0,0,0,1,1,1,1), X=c(0,0,1,1,0,0,1,1), 
  Z=c(0,1,0,1,0,1,0,1))
n <- c(n000, n001, n010, n011, n100, n101, n110, n111)
b <- ivbounds(data=data, Z="Z", X="X", Y="Y", weights=n)
summary(b)

#...or with vector of probabilities
p <- n/rep(c(n0, n1), 4)
names(p) <- c("p00.0", "p00.1", "p01.0", "p01.1", 
  "p10.0", "p10.1", "p11.0", "p11.1") 
b <- ivbounds(data=p)
summary(b)



}
