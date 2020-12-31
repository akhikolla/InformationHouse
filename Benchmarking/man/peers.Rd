% $Id: peers.Rd 205 2019-07-10 14:22:51Z lao $
\name{peers}
\alias{peers}
\alias{get.number.peers}
\alias{get.which.peers}
\alias{get.peers.lambda}


\title{
Find peer firms and units
}


\description{ The function \code{peers} finds for each firm its peers,
   \code{get.number.peers} finds for each peer the number of times
   this peer apears as a peer, and \code{get.which.peers} determines
   for one or more peers the firms they apear as peers for.   Also
   include a function \code{get.peers.lambda} to calculate for fimrs the 
   importances (lambdas) of peers.
   }


\usage{
  peers(object, NAMES = FALSE, N=1:dim(object$lambda)[1], LAMBDA=0)
  get.number.peers(object, NAMES = FALSE, N=1:dim(object$lambda)[2], LAMBDA=0)
  get.which.peers(object, N = 1:dim(object$lambda)[2], LAMBDA=0)
  get.peers.lambda(object,  N=1:dim(object$lambda)[1], LAMBDA=0)
}


\arguments{

  \item{object}{An object of class Farrell as returned by the functions
      \code{\link{dea}}, \code{\link{dea.direct}} et al.}

  \item{NAMES}{If true then names for the peers are returned if names
      are available otherwise the unit index numbers are used. If \code{NAMES}
		is a list of names with length equal to the number of units then it is 
		used as names for peers. }

  \item{N}{The firm(s) or peer(s) for which to get the results.}
  
  \item{LAMBDA}{Minimum weight for extracted peers, i.e. the extracted peers have 
  lambda values larger than \code{LAMBDA}. }

}


\details{ The returned values are index of the firms and can be used
   by itself, but can also by used as an index for a variable with
   names of the firms. 

   The \code{peers} returns a matrix with numbers for the peers for
   each firm; for firms with efficiency 1 the peers are just the firm
   itself.  If there is slack in the evaluation of a firm with
   efficiency 1, this can be found with a call to \code{\link{slack}},
   either directly or by the argument \code{SLACK} when a function
   \code{\link{dea}} was called to generate the Farrell \code{object}. 

   The \code{get.number.peers} returns the number of firms that a peer
   serves as a peer for.
   
   The \code{get.peers.lambda} returns a list of firms with the peers 
   and corresponding value af lambda.
}


\value{ The return values are firm numbers. If the argument
   \code{NAMES=TRUE} is used in the function \code{peers} the return
   is a list of names of the peers if names for the firms are avilable
   as rownames.}


\references{ Peter Bogetoft and Lars Otto; \emph{Benchmarking with
   DEA, SFA, and R}; Springer 2011. Sect. 4.6 page 93 }


\author{ Peter Bogetoft and Lars Otto \email{larsot23@gmail.com} }


\note{ Peers are defines as firms where the coresponding lambdas are
  positive.

  Note that peers might change between a Farrell object return from
  \code{dea} with \code{SLACK=FALSE} and a call with \code{SLACK=TRUE}
  or a following call to the function \code{slack} because a peer on
  the frontier with slack might by the call to \code{dea} be seen as a
  peer for itself whereas this will not happen when slacks are
  calculated.  }



\seealso{
\code{\link{dea}}
}


\examples{
x <- matrix(c(100,200,300,500,100,200,600),ncol=1)
y <- matrix(c(75,100,300,400,25,50,400),ncol=1)

e <- dea(x,y)
peers(e)
get.number.peers(e)

# Who are the firms that firm 1 and 4 is peers for
get.which.peers(e, c(1,4))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
