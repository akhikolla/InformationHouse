\name{generate.Guilds.Cond}
\alias{generate.Guilds.Cond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Artificial data under the GUILDS model, conditioned on Guild size

}
\description{
 Using this function it is possible to generate a community dataset consisting of two separate abundance vectors for each guild, where the data generated adhere to the Guilds model. Data generated is conditioned on guild size.
}
\usage{
generate.Guilds.Cond(theta, alpha_x, alpha_y, JX, JY)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{ Fundamental Biodiversity Number theta

}
  \item{alpha_x}{ Dispersal Ability of Guild X

}
  \item{alpha_y}{Dispersal Ability of Guild Y
}

  \item{JX}{Total number of individuals in Guild X

}
  \item{JY}{Total number of individuals in Guild Y
}
}

\value{
    \item{guildX}{Vector containing the unlabeled abundances of species in Guild X}
    \item{guildY}{Vector containing the unlabeled abundances of species in Guild Y}
}

\author{
 Thijs Janzen
}



\examples{
generate.Guilds.Cond(theta=200,alpha_x = 0.005, alpha_y = 0.001,JX=15000,JY=5000);
}
\keyword{ Guilds }
\keyword{ Neutral Theory }% __ONLY ONE__ keyword per line
