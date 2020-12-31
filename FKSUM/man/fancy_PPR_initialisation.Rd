\name{fancy_PPR_initialisation}
\alias{fancy_PPR_initialisation}
\title{Initialisation for PPR based on Ridge LM after GAM type smoothing}
\description{
  Computes an initial projection for PPR by replacing each covariate by the smoothed response against it. Then fit a
  standard ridge estimator to the transformed data. Not intended for general use.
}
\keyword{internal}
