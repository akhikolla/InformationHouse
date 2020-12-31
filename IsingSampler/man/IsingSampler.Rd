\name{IsingSampler}
\alias{IsingSampler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample states from the Ising model
}
\description{
This function samples states from the Ising model using one of three methods. See details.
}
\usage{
IsingSampler(n, graph, thresholds, beta = 1, nIter = 100, responses = c(0L, 1L), 
    method = c("MH", "CFTP", "direct"), CFTPretry = 10, constrain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of states to draw
}
  \item{graph}{
Square matrix indicating the weights of the network. Must be symmetrical with 0 as diagonal.
}
  \item{thresholds}{
Vector indicating the thresholds, also known as the external field.
}
  \item{beta}{
Scalar indicating the inverse temperature.
}
  \item{nIter}{
Number of iterations in the Metropolis and exact sampling methods.
}
  \item{responses}{
Response options. Typically set to \code{c(-1L, 1L)} or \code{c(0L, 1L)} (default). Must be integers!
}
  \item{method}{
The sampling method to use. Must be \code{"MH"}, \code{"CFTP"} or \code{"direct"}. See details.
}
 \item{CFTPretry}{
 The amount of times a sample from CFTP may be retried. If after 100 couplings from the past the chain still results in \code{NA} values the chain is reset with different random numbers. Be aware that data that requies a lot of CFTP resets might not resemble exact samples anymore.
 }
 \item{constrain}{A (number of samples) by (number of nodes) matrix with samples that need be constrained; \code{NA} indicates that the sample is unconstrained. Defaults to a matrix of \code{NA}s.}
}
\details{
This function uses one of three sampling methods. \code{"MH"} can be used to sample using a Metropolis-Hastings algorithm. The chain is initiated with random values from the response options, then for each iteration for each node a node is set to the second response option with the probability of that node being in the second response option given all other nodes and parameters. Typically, 100 of such iterations should suffice for the chain to converge.

The second method, \code{"CFTP"} enhances the Metropolis-Hastings algorithm with Coupling from the Past (CFTP; Murray, 2007) to draw exact samples from the distribution. This is slower than the default Metropolis-Hastings but guarantees exact samples. However, it does depend on the graph structure and the number of nodes if these exact samples can be obtained in feasable time.

The third option, \code{"direct"}, simply computes for every possibly state the probability and draws samples directly from the distribution of states by using these probabilities. This also guarantees exact samples, but quickly becomes intractible (roughly above 10 nodes).
}
\value{
A matrix containing samples of states.
}
\references{
Murray, I. (2007). Advances in Markov chain Monte Carlo methods.}
\author{
Sacha Epskamp (mail@sachaepskamp.com)
}

\seealso{
\code{\link{IsingSampler-package}} for examples
}
\examples{
## See IsingSampler-package help page
}
