\docType{package}
\name{FSSF-package}
\alias{FSSF}
\alias{FSSF-package}
\title{Large-Scale Fully-Sequential Space-filling Algorithms for Computer Experiments}
\description{
  Provides large-scale computer experimental design tools to generate fully-sequential (a nested sequence of designs with points added one at a time) space-filling designs inside a unit hypercube (Shang and Apley, 2019).
}

\details{
  Fully-sequential (i.e., with design points added one at a time) space-filling designs are useful for global surrogate modeling of expensive computer experiments when the number of design points required to achieve a suitable accuracy is unknown in advance. We provide three fully-sequential space-filling (FSSF) design algorithms that are conceptually simple and computationally efficient and that achieve much better space-filling properties than alternative methods such as Sobol sequences and more complex batch-sequential methods based on sliced or nested optimal Latin hypercube designs (LHDs). 
  
  Brief descriptions of the main functions are provided below:
  
  \code{fssf_f()} generates maximin designs using a similar idea with Kennard and Stone (1969). The "f" stands for "forward", since the algorithm begins with the smallest design and adds points one at a time.
  
  \code{fssf-b()} generates large-size maximin designs efficiently. The "b" stands for "backward", since the algorithm begins with the largest design and removes points one at a time.
  
  \code{fssf-fr()} uses a heuristic to generate what is closer to a minimax design. The "fr" stands for "forward reflected", since the algorithm uses a modification of the criterion used in the fssf-f algorithm. 
}

\references{
Shang, B. and Apley, D.W. (2019), "Large-Scale Fully-Sequential Space-filling Algorithms for Computer Experiments", Journal of Quality Technology (in press). doi:10.1080/00224065.2019.1705207. 
		
Kennard, R.W. and Stone, L.A. (1969). "Computer aided design of experiments". Technometrics 11.1, pp. 137-148.
}

\author{Boyang Shang \email{boyangshang2015@u.northwestern.edu} \cr Daniel W. Apley \email{apley@northwestern.edu}}


\examples{
#
# See the examples in the help pages for the main functions mentioned above.
#
}

