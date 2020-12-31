% $Id: make.merge.Rd 118 2011-07-22 23:28:53Z Lars $
\name{make.merge}
\alias{make.merge}

\title{ Make an aggregation matrix to perform mergers }

\description{ Make an aggregation matrix to perform mergers of
 firms. The matrix can be post multiplied (matrix multiplication) to
 input and output matrices to make merged input and output matrices.
 }

\usage{
make.merge(grp, nFirm = NULL, X = NULL, names = NULL)
}


\arguments{

 \item{grp}{ Either a list of length Kg for Kg firms after mergers;
       each component of the list is a (named) list with the firm numbers or
       names going into this merger.

       Or a factor of length K with Kg levels where where each level
       determines a merger; to exclude firms for mergers set the
       factor value to NA. }

 \item{nFirm}{ Number of firms before the mergers }

 \item{X}{A matrix of inputs or outputs where the rows corresponds to
       the number of original (starting) firms }

 \item{names}{A list with names of all firms, only needed if the
       mergers are given as a list of names, i.e. \code{grp} is a list
       of names.}

}


\details{ 

Either \code{nFirm} or \code{X} must be present; if both are present
then \code{nFirm} must be equal to the number of rows in \code{X}, the
number of firms.

When \code{X} is an input matrix of dimension \code{K x m}, \code{K}
firms and \code{m} inputs, and \code{M <- make.merge(gr,K)} then
\code{M \%*\% X} is the input matrix for the merged firms.

}


\value{ Returns an aggregation matrix of dimension Kg times K where
rows corresponds to new merged firms and columns are 1 for firms to be
included and 0 for firms to be excluded in the the given merger as
defined by the row. }


\note{ The argument \code{TRANSPOSE} has not been implemented for this
function. If you need transposed matrices you must transpose the
merger matrix yourself. If you define mergers via factors there is no
need to transpose in the arguments; just do not use \code{X} in the
arguments.  }



\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com}}


\seealso{ \code{\link{dea.merge}} }


\examples{
# To merge firms 1,3, and 5; and 2 and 4 of 7 firms into 2 new firms
# the aggregation matrix is M; not all firms are involved in a merger.
M <- make.merge(list(c(1,3,5), c(2,4)),7)
print(M)

# Merge 1 and 2, and 4 and 5, and leave 3 alone, total of 5 firms.
# Using a list
M1 <- make.merge(list(c(1,2), c(4,5)), nFirm=5)
print(M1)

# Using a factor
fgr <- factor(c("en","en",NA,"to","to"))
M2 <- make.merge(fgr)
print(M2)

# Name of mergers
M3 <- make.merge(list(AB=c("A","B"), DE=c("D","E")), names=LETTERS[1:5])
print(M3)

# No name of mergers
M4 <- make.merge(list(c("A","B"), c("D","E")), names=LETTERS[1:5])
print(M4)
}

