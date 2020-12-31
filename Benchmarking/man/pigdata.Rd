% $Id: pigdata.Rd 104 2011-02-06 11:56:03Z Lars $
\name{pigdata}
\alias{pigdata}
\docType{data}

\title{ Data: Multi-output pig producers }

\description{ Input and output data for 248 pig producers that
also produces crop, i.e. a multi--output data set.  }

\usage{data(pigdata)}

\format{
  A data frame with 248 observations on the following 16 variables.
  \describe{
    \item{\code{firm}}{Serial number for pig producer}
    \item{\code{x1}}{Input fertilizer}
    \item{\code{x2}}{Input feedstuf}
    \item{\code{x3}}{Input land}
    \item{\code{x4}}{Input labour}
    \item{\code{x5}}{Input machinery}
    \item{\code{x6}}{Input other capital}
    \item{\code{y2}}{Output crop}
    \item{\code{y4}}{Output pig}
    \item{\code{w1}}{Price of fertilizer}
    \item{\code{w2}}{Price of feedstuf}
    \item{\code{w3}}{Price of land}
    \item{\code{w4}}{Price of labour}
    \item{\code{w5}}{Price of michenery}
    \item{\code{w6}}{Price of other capital}
    \item{\code{p2}}{Price of pig}
    \item{\code{p4}}{Price of crop}
    \item{\code{cost}}{Total cost, w1*x1+...+w6*x6.}
    \item{\code{rev}}{Total revenue, p2*y2+p4*y4.}
  }
}

\details{ In raising pigs, most farmers also produce crops to feed the
pigs. Labor and capital are used not just directly for pig-raising but
also on the field.  }

\note{ Data as \code{.csv} are loaded by the command \code{data} using
\code{ read.table(..., header = TRUE, sep = ";")} such that this file
is a semicolon separated file and not a comma separated file.  }

\source{ Farmers accounting data converted to index. }

\examples{
data(pigdata)
## maybe str(pigdata) ; plot(pigdata) ...
}
\keyword{datasets}
