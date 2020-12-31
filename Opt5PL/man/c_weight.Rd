\name{c_weight}
\alias{c_weight}

\title{One iteration to run Newton Raphson to get c-optimal weights

}
\description{This is a sub-part of running Newton Raphon method to search c-optimal weights for the given design points. This is function probvides the obtained c-optimal weights after one iteration of Newton Raphson method.
}
\usage{
c_weight(W, T, X, d, p, order, UB, I)
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
  \item{p}{Numeric. Define EDp

}
  \item{order}{Numeric. The number of model paraemters

}
  \item{UB}{Numeric. The upper bound of the design points

}
  \item{I}{A numeric matrix. Identity matrix.

}
}
