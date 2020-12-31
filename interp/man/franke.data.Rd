\name{franke.data}
\alias{franke.data}
\alias{franke.fn}
\alias{franke}
\title{
Test datasets from Franke for interpolation of scattered data
}
\description{
\code{franke.data} generates the test datasets from Franke, 1979, see references.
}
\usage{
franke.data(fn = 1, ds = 1, data)
franke.fn(x, y, fn = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{
    function number, from 1 to 5.
  }
  \item{x}{'x' value}
  \item{y}{'y' value}
  \item{ds}{
    data set number, from 1 to 3. Dataset 1 consists of 100 points,
    dataset 2 of 33 points and dataset 3 of 25 points scattered in the
    square \eqn{[0,1]\times[0,1]}{[0,1]x[0,1]}. (and partially slightly
    outside).
  }
  \item{data}{
    A list of dataframes with 'x' and 'y' to choose from, dataset
    \code{franke} should be used here.
  }
}
\details{
  These datasets are mentioned in Akima, (1996) as a testbed for the
  irregular scattered data interpolator.

  Franke used the five functions:

  \deqn{0.75e^{-\frac{(9x-2)^2+(9y-2)^2}{4}}+
    0.75e^{-\frac{(9x+1)^2}{49}-\frac{9y+1}{10}}+
    0.5e^{-\frac{(9x-7)^2+(9y-3)^2}{4}}-
      0.2e^{-((9x-4)^2-(9y-7)^2)}
    }{0.75*exp(-((9*x-2)^2+(9*y-2)^2)/4)+
      0.75*exp(-((9*x+1)^2)/49-(9*y+1)/10)+
      0.5*exp(-((9*x-7)^2+(9*y-3)^2)/4)-
      0.2*exp(-(9*x-4)^2-(9*y-7)^2)}


  \deqn{\frac{\mbox{tanh}(9y-9x)+1}{9}}{(tanh(9*y-9*x)+1)/9}

  \deqn{\frac{1.25+\cos(5.4y)}{6(1+(3x-1)^2)}}{(1.25+cos(5.4*y))/(6*(1+(3*x-1)^2))}

  \deqn{e^{-\frac{81((x-0.5)^2+\frac{(y-0.5)^2}{16})}{3}}}{exp(-81*((x-0.5)^2+(y-0.5)^2)/16)/3}

  \deqn{e^{-\frac{81((x-0.5)^2+\frac{(y-0.5)^2}{4})}{3}}}{exp(-81*((x-0.5)^2+(y-0.5)^2)/4)/3}

  \deqn{\frac{\sqrt{64-81((x-0.5)^2+(y-0.5)^2)}}{9}-0.5}{sqrt(64-81*((x-0.5)^2+(y-0.5)^2))/9-0.5}

  and evaluated them on different more or less dense grids over \eqn{[0,1]\times[0,1]}{[0,1]x[0,1]}.
}
\value{
A data frame with components
  \item{x }{'x' coordinate}
  \item{y }{'y' coordinate}
  \item{z }{'z' value}
}
\note{
  The datasets have to be generated via \code{franke.data} before
  use, the dataset \code{franke} only contains a list of 3 dataframes of
  'x' and 'y' coordinates for the above mentioned irregular grids.
  Do not forget to load the \code{franke} dataset first.

  The 'x' and 'y' values have been taken from Akima (1996).
}
\references{
FRANKE, R., (1979). A critical comparison of some methods for interpolation
of scattered data. Tech. Rep. NPS-53-79-003, Dept. of Mathematics, Naval
Postgraduate School, Monterey, Calif.

Akima, H. (1996). Algorithm 761: scattered-data surface fitting that has
the accuracy of a cubic polynomial.
ACM Transactions on Mathematical Software \bold{22}, 362--371.
}
\author{
  Albrecht Gebhardt <albrecht.gebhardt@aau.at>,
  Roger Bivand <roger.bivand@nhh.no>
}

\seealso{
\code{\link{interp}}
}
\examples{
## generate Frankes data set for function 2 and dataset 3:
data(franke)
F23 <- franke.data(2,3,franke)
str(F23)
}
\keyword{ datagen }

