\name{surv}
\alias{surv}
\title{Print fitted Survival and Hazard Probabilites of a DStree Object}
\usage{
surv(object)
}
\arguments{
  \item{object}{fitted model object of class
  \code{"DStree"}. This is assumed to be the result of some
  function that produces an object with the same named
  components as that returened by the \code{DStree}
  function.}
}
\value{
Two matrices containing \code{Survival} and \code{Hazard} probabilities,
where each row denotes fitted probabilities
in the respective terminal leaves. The first column refers to the
probability of the first observed timepoint, the second column
to the second timepoint, etc.
}
\description{
This function prints the fitted survival and hazard
probabilities of every terminal node of a \code{DStree}
object.
}
\examples{
# Grow tree
tree <- DStree(spell~ui+age,status="censor1",data=UnempDur)

# Print fitted probabilites
surv(tree)
}

