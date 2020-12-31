\name{ah}

\alias{ah}

\title{
Fitting semiparametric additive hazards regression models. 
}

\description{
\code{ah} is a wrapper around the \code{ahaz} function in the \pkg{ahaz}
package, with a more user-friendly and standard interface. Refer to the 
manual for \code{ahaz} for details. 
}     

\usage{
ah(formula, data, weights, robust=FALSE)
}

\arguments{
\item{formula}{
an object of class "formula": a symbolic description of the model to be fitted. 
}
  \item{data}{
a data frame containing the variables in the model. 
}
  \item{weights}{
an optional vector of prior weights to be used in the fitting process. 
}
  \item{robust}{
robust calculation of variance; see manual for \code{ahaz}.
}
}

\details{ 
 See manual for \code{ahaz}.
 
}

\value{
An object of class \code{"ah"} is a list containing the same elements as 
an object of class \code{"ahaz"}, plus 
\item{call}{
  the matched call. 
  }
\item{formula}{
   the formula argument.
  }
\item{coefficients}{
  a named vector of estimated coefficients. 
  }
\item{vcov}{
  the variance-covariance matrix for the estimated coefficients.  
  }
\item{incl}{
  the \code{ahaz} function does not allow for missing data. Thus, if \code{data}
  contains rows with missing data on any of the variables in the model, then 
  these rows are excluded before calling \code{ahaz}. \code{incl} is a vector 
  containing the rownames of those rows that are included in the analysis, 
  that is, the rows with no missing data on any of the variables in the model.          
  }


}

\references{
Lin D.Y., Ying Z. (1994). Semiparametric analysis of the additive risk model. 
\emph{Biometrika} \bold{81}(1), 61-71.
}

\author{
Arvid Sjolander. 
}

\note{
  The \code{ahaz} function does not allow for ties. Thus, before calling 
  \code{ah} any ties have to be manually broken.  
}

\examples{

require(ahaz)

##This example is adapted from the example given for the ahaz function 
##in the ahaz package 

data(sorlie)

# Break ties
set.seed(10101)
sorlie$time <- sorlie$time+runif(nrow(sorlie))*1e-2

# Fit additive hazards model 
fit <- ah(formula=Surv(time, status)~X13+X14+X15+X16+X17+X18+X19+X20+X21+X22, 
  data=sorlie)
summary(fit)


}
