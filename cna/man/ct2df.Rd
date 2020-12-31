\name{ct2df}
\alias{ct2df}

\title{
Transform a configuration table into a data frame
}

\description{
Transform a configuration table into a data frame. This is the converse function of \code{\link{configTable}}.
}

\usage{
ct2df(ct, tt)
}

\arguments{
  \item{ct}{A \code{configTable}.}
  \item{tt}{Argument \code{tt} is deprecated in \code{ct2df()}; use \code{ct} instead.}
}

\details{
Rows in the \code{configTable} corresponding to several cases are rendered as multiple rows in the resulting data frame.
}

\value{
A data frame.
}

\seealso{\code{\link{configTable}}
}

\examples{
ct.educate <- configTable(d.educate[1:2])
ct.educate
ct2df(ct.educate)

dat1 <- some(configTable(allCombs(c(2, 2, 2, 2, 2)) - 1), n = 200, replace = TRUE)
dat2 <- selectCases("(A*b + a*B <-> C)*(C*d + c*D <-> E)", dat1)
dat2 
ct2df(dat2)

dat3 <- data.frame(
  A = c(1,1,1,1,1,1,0,0,0,0,0),
  B = c(1,1,1,0,0,0,1,1,1,0,0),
  C = c(1,1,1,1,1,1,1,1,1,0,0),
  D = c(1,0,0,1,0,0,1,1,0,1,0),
  E = c(1,1,0,1,1,0,1,0,1,1,0)
  )
ct.dat3 <- configTable(dat3, frequency = c(4,3,5,7,4,6,10,2,4,3,12))
ct2df(ct.dat3)
}