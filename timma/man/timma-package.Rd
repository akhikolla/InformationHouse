\name{timma-package}
\alias{timma-package}
\alias{TIMMA}
\docType{package}
\title{
Target Inhibition inference using Maximization and Minimization Averaging
}
\description{
Due to the exponentially increasing number of potential drug and target combinations, it is meaningful to select the most promising combinations based on computational models. The TIMMA model was proposed to utilize drug-target interaction data and drug sensitivity data to infer the effects of drug combinations. This R package TIMMA is the implementation of the TIMMA model. It consists of the following components: (a) model selection using the sffs algorithm; (b) model construction using the maximization and minimization averaging rules; (c) ranking of drug combinations according to their synergy scores and a target inhibition network.
}
\details{
\tabular{ll}{
Package: \tab TIMMA\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-10-07\cr
License: \tab Artistic License 2.0\cr
}

}
\author{
Liye He \email{liye.he@helsinki.fi}
}
\references{
Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T.
Target inhibition networks: predicting selective combinations of druggable targets to block cancer
survival pathways. PLOS Computational Biology 2013; 9: e1003226.
}
\keyword{ drug discovery }
\examples{
\dontrun{
data(tyner_interaction_binary)
data(tyner_sensitivity)
median_sensitivity<-tyner_sensitivity[, 1]
results<-timma(tyner_interaction_binary, median_sensitivity)
}
}
