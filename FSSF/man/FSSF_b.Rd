\name{fssf_b}
\alias{fssf_b}
\title{Generate fully-sequential maximin designs inside a unit hypercube.}
\description{
  Produces a random fully-sequential design (a nested sequence of designs with points added one at a time) inside a unit hypercube such that the design points are as far away from each other as possible (Shang and Apley, 2019). The "b" stands for "backward", since the algorithm begins with the largest design and removes points one at a time. \code{fssf_b} is much faster than \code{fssf_f} when the design size is large.
}

\usage{
	fssf_b(d, nMax, N=-1, eps=0.5)
}

\arguments{
	\item{d}{
		The dimension of the design space.
	}
	\item{nMax}{
		The largest design size required by the user.
	}
	\item{N}{
		Size of the candidate set used to generate the design points. -1 corresponds to the default setting, and the candidate set size will be calculated as  \eqn{1000 \times d + 2 \times nMax}{1000d+2nMax}. Using large \code{N} will make the design more space-filling, but will slow down the program.
	}
	\item{eps}{
		The error bound for approximate nearest neighbor searching. Default value is 0.5 and \code{eps} must be greater or equal to 0. Using large \code{eps} will make the program run faster, but will cause worse space-filling performances. 
	}
}

\value{
	A \eqn{nMax \times d}{nMax x d} matrix with the \eqn{i^{th}}{i-th} row corresponding to the \eqn{i^{th}}{i-th} design point.
}

\details{
	The \code{fssf_b} function calls portions of the Approximate Nearest Neighbor Library, version 1.1.2,  written by David M. Mount and Sunil Arya to do the neareast neighbor search. Some changes have been made to the original ANN library to suit the needs of the \code{fssf_b} function. More information about the ANN library can be found in the ANN Programming Manual at \url{http://www.cs.umd.edu/~mount/ANN}
}


\references{
Shang, B. and Apley, D.W. (2019), "Large-Scale Fully-Sequential Space-filling Algorithms for Computer Experiments", Journal of Quality Technology (in press). doi:10.1080/00224065.2019.1705207.
		
Kennard, R.W. and Stone, L.A. (1969). "Computer aided design of experiments". Technometrics 11.1, pp. 137-148.
}

\author{Boyang Shang \email{boyangshang2015@u.northwestern.edu} \cr Daniel W. Apley \email{apley@northwestern.edu}}


\examples{
##Generate a design using the fssf_b function.
Design <- fssf_b(d=2, nMax = 320)
plot(Design[,1], Design[,2])
}

