\name{GUILDS-package}
\alias{GUILDS}
\docType{package}
\title{Package implementing the Guilds sampling formula for the Neutral Theory of Biodiversity
}
\description{ The GUILDS package contains a number of sampling formula's being the Etienne Sampling Formula (Etienne 2005), the GUILDS sampling formula (Janzen et al. 2014) and the GUILDS sampling formula conditioned on guild Size (Janzen et al. 2015). Furthermore it contains functions to generate data given the guilds model, with or without conditioning on guild size. C++ Code to obtain Sterling numbers of the first kind was adopted from the Tetame program by Jabot et al. (2008).   
\cr\cr
Updates \cr\cr
Version 1.3   : GUILDS is now on GitHub: https://github.com/thijsjanzen/GUILDS\cr
Version 1.3   : Wrote code tests to check code integrity, code coverage is >95\% \cr
Version 1.3   : Modified maximum likelihood functions to take into account theta_x = theta_y = theta / 2 \cr
Version 1.3   : Added a plotting function to plot Preston style plots
Version 1.2.1 : Updated the User manual\cr
Version 1.2   : fixed memory leak issues by adding extra vector access checks \cr
Version 1.2   : fixed memory leak issues by introducing vectors in KDA code\cr
Version 1.2   : renamed logLik to avoid shadowing of the function logLik in the package stats\cr
Version 1.1   : removed malloc header from KDA code \cr
}
\details{
\tabular{ll}{
Package: \tab GUILDS\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
License: \tab GPL 2.0 \cr
}

}
\author{
Thijs Janzen

Maintainer: Thijs Janzen <thijsjanzen@gmail.com>
}
\references{
Janzen, T., Haegeman B., Etienne, R.S. (2015) A sampling formula for communities with multiple dispersal syndromes. Journal of Theoretical Biology 374: 94-106

Etienne, R.S. (2005). A new sampling formula for neutral biodiversity. Ecology Letters, 8(3), 253-260.

Jabot, F., Etienne, R.S., & Chave, J. (2008). Reconciling neutral community models and environmental filtering: theory and an empirical test. Oikos 117: 1308-1320
}

\keyword{GUILDS}
\keyword{Etienne Sampling Formula}
\keyword{Neutral Theory}