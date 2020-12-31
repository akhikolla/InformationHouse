\name{GeneUni}
\docType{data}
\alias{GeneUni}
\title{Example of genetic data for univariate error-prone covariates}
\description{
  This data set gives an example data for the illustration of usage of \code{\link{augSIMEX}} function. The data set is adapted from the outbred Carworth Farms White mice data (Parker et. al., 2016). The dataset contains main data and validation data. We are interested to study the association between the response genotype of rs38916331 and the tibia length, which is subject to mismeasurement. 
  
  Tibia_T: the precise measurement of tibia length
  Tibia_M: the error-prone measurement of tibia length
  Batch_T: the precisely labeled batch effect
  Batch_O: the observed batch effect
  Weight:  the body weights of the mice
}
\usage{data(GeneUni)}
\references{
Parker C C, Gopalakrishnan S, Carbonetto P, et al. Genome-wide association study of behavioral, physiological and gene expression traits in outbred CFW mice[J]. Nature genetics, 2016, 48(8): 919.
}
\format{A list of two data frames. Main data (rs38916331, Weight, Tibia_M, Batch_O) and validation data (Tibia_T, Batch_T, Weight, Tibia_M, Batch_O). }
\keyword{datasets}