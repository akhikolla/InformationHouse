\name{logLikelihood.Guilds.Conditional}
\alias{logLikelihood.Guilds.Conditional}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood of the Guilds sampling formula, conditional on guild size
}
\description{
This function calculates the likelihood of the guilds model, conditional on guild size; provided abundance data and parameter values.
}
\usage{
logLikelihood.Guilds.Conditional(parameters, model, sadx, sady, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
 \arguments{
  \item{parameters}{\code{parameters} corresponds to a vector of parameter values depending on the provided model: \cr
		- model: "D0" \code{parameters} = c(theta, alpha) \cr
		- model: "D1" \code{parameters} = c(theta, alpha X, alpha Y) \cr
	}
  \item{model}{ The chosen model to calculate the likelihood for, please note that the vector of parameters should contain the corresponding parameters in the right order. The user can pick one of these models:\cr
		-  "D0" \cr
		-  "D1" \cr
	}
  \item{sadx}{ The Species Abundance Distribution of guild X}
  \item{sady}{ The Species Abundance Distribution of guild Y} 
  \item{verbose}{TRUE/FALSE flag, indicates whether intermediate output is shown on screen}
}
\value{returns the LogLikelihood}

\author{
 Thijs Janzen
}

\examples{
exampleData <- generate.Guilds.Cond(theta = 200, alpha_x = 0.005, alpha_y = 0.001,
                                    JX = 1000, JY = 2000)
parametervals <- c(200,0.005,0.001) #theta = 200, alpha X = 0.005, alpha Y = 0.001
LL = logLikelihood.Guilds.Conditional(parametervals, model="D1",
				      exampleData$guildX, exampleData$guildY, verbose=TRUE)
}
