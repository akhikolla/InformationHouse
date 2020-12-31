\name{snip}
\alias{snip}
\title{Cost-complexity Pruning of a DStree Object}
\usage{
snip(tree, cp)
}
\arguments{
  \item{tree}{a fitted object of class \code{"DStree"}}

  \item{cp}{complexity parameter to which the \code{DStree}
  object will be trimmed}
}
\value{
A new \code{DStree} object that is trimmed to the value
\code{cp}.
}
\description{
Determines a nested sequence of subtrees of the supplied
\code{DStree} object by recursively snipping off the
least important splits, based on the complexity parameter
(\code{cp}).
}
\examples{
data(cost)
d.cost <- dis.cost(cost)

fit <- DStree(time~prevStroke+age+sex+alcohol+smoke,status="status",data=d.cost)

sfit<-snip(fit,cp=0.02)
plot(sfit) #plot smaller DStree object
}

