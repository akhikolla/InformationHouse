\name{fssf_f}
\alias{fssf_f}
\title{Generate fully-sequential maximin designs inside a unit hypercube.}
\description{
  Produces a random fully-sequential design (a nested sequence of designs with points added one at a time) inside a unit hypercube such that the design points are as far away from each other as possible (Shang and Apley, 2019). The "f" stands for "forward", since the algorithm begins with the smallest design and adds points one at a time.
}

\usage{
	fssf_f(d, nMax, N=-1, ScaleVector = NULL, Dinit = NULL)
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
	\item{ScaleVector}{Array of the lengthscale parameters of different inputs. Default is NULL, which corresponds to the \code{ScaleVector} being a unit vector of length \code{d}. When \code{ScaleVector} is not NULL, for instance, \code{ScaleVector} is \eqn{(\theta_1, \cdots, \theta_d)}, the distance between point \eqn{(x_1, \cdots, x_d)} and point \eqn{(y_1, \cdots, y_d)} will be computed as \eqn{\sum_{j=1}^{d} {\Big( \frac{y_j - x_j}{\theta_j} \Big)}^2}}  
	\item{Dinit}{Numerical Matrix with \eqn{n_{init}} rows and \code{d} columns, where \eqn{n_{init}} is a user-specified parameter. This is an optional initial design with size \eqn{n_{init}} provided by the user. Default is NULL, which corresponds to no initial design. If \code{Dinit} is not NULL, then the algorithm will select \code{nMax} additional design points taking into account of this initial design.}
}

\value{
	A \eqn{nMax \times d}{nMax x d} matrix with the \eqn{i^{th}}{i-th} row corresponding to the \eqn{i^{th}}{i-th} design point.
}

\details{
	The \code{fssf_f} uses a similar idea as proposed by Kennard and Stone(1969); modifications have been made to improve the space-filling performance as well as efficiency.
}


\references{
Shang, B. and Apley, D.W. (2019), "Large-Scale Fully-Sequential Space-filling Algorithms for Computer Experiments", Journal of Quality Technology (in press). doi:10.1080/00224065.2019.1705207. 
		
Kennard, R.W. and Stone, L.A. (1969). "Computer aided design of experiments". Technometrics 11.1, pp. 137-148.
}

\author{Boyang Shang \email{boyangshang2015@u.northwestern.edu} \cr Daniel W. Apley \email{apley@northwestern.edu}}


\examples{
##Generate a design using the fssf_f function with no scaling vector and no initial design.
Design <- fssf_f(d=2, nMax = 320)
plot(Design[,1], Design[,2])

##Generate a design using the fssf_f function with scaling vector and no initial design.
d = 2
n = 100
ScaleVector = c(1.0, 20.0)*0.5
Design = fssf_f(d = d, nMax = n, ScaleVector = ScaleVector)
plot(Design[,1], Design[,2])

##Generate a design using the fssf_f function with a scaling vector and with an initial design
d = 2
n = 100
Dinit =  fssf_f(d=2, nMax = 40)
ScaleVector = c(1.0, 20.0)*0.5
Design = fssf_f(d = d, nMax = n, ScaleVector = ScaleVector, Dinit = Dinit)
plot(Design[,1], Design[,2])
points(Dinit[,1], Dinit[,2], col="red")
}

