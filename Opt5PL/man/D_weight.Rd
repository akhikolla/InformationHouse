\name{D_weight}
\alias{D_weight}

\title{One iteration to run Newton Raphson to get D-optimal weights

}
\description{This is a sub-part of running Newton Raphon method to search D-optimal weights for the given design points. This is function probvides the obtained D-optimal weights after one iteration of Newton Raphson method.

}
\usage{
D_weight(W, T, X, d, q)
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
  \item{q}{A numeric vector. User select weights defined in the function RDOPT

}
}
