\name{GeneRepeat}
\docType{data}
\alias{GeneRepeat}
\title{Example data for univariate error-prone covariates in repeated measurements case}
\description{
 This dataset gives an example data for the illustration of usage of \code{\link{augSIMEX}} function. The data set is adapted from the outbred Carworth Farms White mice data (Parker et. al., 2016). The dataset contains main data and validation data. We are interested to study the association between the response genotype of rs223979909 and the locomotor response to methamphetamine injections, which is subject to mismeasurement. 
  
  Tibia_5 - Tibia_30: the repeated measurement of tibia length
  Tibia_M: the error-prone measurement of tibia length
  Batch_T: the precisely labeled batch effect
  Batch_O: the observed batch effect
  Weight:  the body weights of the mice
}
\usage{data(GeneRepeat)}
\references{
Parker C C, Gopalakrishnan S, Carbonetto P, et al. Genome-wide association study of behavioral, physiological and gene expression traits in outbred CFW mice[J]. Nature genetics, 2016, 48(8): 919.
}
\format{A list of two data frames. Main data (rs38916331, Totdist5,  Totdist10, Totdist15, Totdist20, Totdist25, Totdist30, Weight, Batch_O) and validation data (Batch_T, Batch_O).}
\keyword{datasets}

