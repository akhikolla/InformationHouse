\name{DD_weight}
\alias{DD_weight}

\title{One iteration to run Newton Raphson to get Ds-optimal weights

}
\description{This is a sub-part of running Newton Raphon method to search Ds-optimal weights for the given design points. This is function probvides the obtained Ds-optimal weights after one iteration of Newton Raphson method.

}
\usage{
DD_weight(W, T, X, d, I4, I5, order)
}

\arguments{
  \item{W}{A numeric vector. The first K-1 weights for a given design

}
  \item{T}{A numeric vector. Model parameter values

}
  \item{X}{A numeric vector. K dose levels for a given design

}
  \item{d}{Numeric. Step adjustment

}
  \item{I4}{A numeric matrix. Information matrix for the 4PL model

}
  \item{I5}{A numeric matrix. Information matrix for the 5PL model

}
  \item{order}{Numeric. The number of model paraemters

}
}
