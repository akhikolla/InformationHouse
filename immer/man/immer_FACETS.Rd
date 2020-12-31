%% File Name: immer_FACETS.Rd
%% File Version: 0.36

\name{immer_FACETS}
\alias{immer_FACETS}
\alias{summary.immer_FACETS}


%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Wrapper to FACDOS
}

\description{
This Function is a wrapper do the DOS version of FACETS (Linacre, 1999).
}

\usage{
immer_FACETS(title=NULL, convergence=NULL, totalscore=NULL, facets=NULL,
  noncenter=NULL,  arrange=NULL,entered_in_data=NULL, models=NULL,
  inter_rater=NULL, pt_biserial=NULL, faire_score=NULL, unexpected=NULL,
  usort=NULL, positive=NULL, labels=NULL, fileinput=NULL, data=NULL,
  path.dosbox=NULL, path.facets="", model.name=NULL, facetsEXE=NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{title}{
  title of the analysis
  }
  \item{convergence}{
  convergence criteria
  }
  \item{totalscore}{
  show the total score with each observation
  }
  \item{facets}{
  number of specified facets
  }
  \item{noncenter}{
  specified the non centered facet here
  }
  \item{arrange}{
  control the ordering in each table/output
  }
  \item{entered_in_data}{
  optional specification for facets
  }
  \item{models}{
  model to be used in the analysis
  }
  \item{inter_rater}{
  Specify rater facet number for the agreement report among raters
  }
  \item{pt_biserial}{
  correlation between the raw-score for each element
  }
  \item{faire_score}{
  intended for communicating the measures as adjusted ratings
  }
  \item{unexpected}{
  size of smallest standardized residual
  }
  \item{usort}{
  order in which the unexpected observation are listed
  }
  \item{positive}{
  specifies which facet is positively oriented
  }
  \item{labels}{
  name of each facet, followed by a list of elements
  }
  \item{fileinput}{
  optional argument, if your data are stored within a separate file
  }
  \item{data}{
  Input of the data in long-format
  }
  \item{path.dosbox}{
Path to the installed DOSBox. If \code{NULL}: the function assumed that you have
purchased FACETS and would like to use this version (currently only for Windows-User)
}
  \item{path.facets}{
Path to FACDOS or FACETS if the path.dosbox is "NULL"
}
  \item{model.name}{
Name of the configuration file for FACETS
}
  \item{facetsEXE}{
optional argument to specifie specific FACETS.exe
}
}

\details{
Within the function \code{immer_FACETS} it is either possible to pass existing
FACETS input files or to specify the Input within the function.
To run the estimation in FACETS it is necessary to provide both the path to the
DosBox and FACDOS (it is recommended to use the function \code{immer_install} for the
installation process). After the estimation process is finished the Exports are
in the Facets folder.
}

%\value{
%...need to be defined....
%}


\references{
Linacre, J. M. (1999). \emph{FACETS} (Version 3.17)[Computer software]. Chicago: MESA.
}

%% Browne, W. J., & Draper, D. (2006). A comparison of Bayesian and
%% likelihood-based methods for fitting multilevel models.
%% \emph{Bayesian Analysis}, \bold{1(3)}, 473-514.


%\author{
%Jan Steinfeld <j.steinfeld@bifie.at>
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Install FACDOS and DOSBox \code{\link{immer_install}}.
}


\examples{
\dontrun{
################################
# 1. Example on Windows
################################
# define data generating parameters
set.seed(1997)
N <- 500  # number of persons
I <- 4    # number of items
R <- 3    # number of raters
K <- 3    # maximum score
sigma <- 2  # standard deviation
theta <- rnorm( N, sd=sigma )  # abilities
# item intercepts
b <- outer( seq( -1.5, 1.5, len=I), seq( -2, 2, len=K), "+" )
# item loadings
a <- rep(1,I)
# rater severity parameters
phi <- matrix( c(-.3, -.2, .5), nrow=I, ncol=R, byrow=TRUE )
phi <- phi + rnorm( phi, sd=.3 )
phi <- phi - rowMeans(phi)
# rater variability parameters
psi <- matrix( c(.1, .4, .8), nrow=I, ncol=R, byrow=TRUE )
# simulate HRM data
data <- immer::immer_hrm_simulate( theta, a, b, phi=phi, psi=psi )

# prepare data for FACETS
data2FACETS <- function(data){
  tmp <- match(c("pid","rater"),colnames(data))
  items <- grep("I",colnames(data))
cbind(data[, match(c("pid","rater"),
colnames(data))],gr=paste0("1-",length(items)),data[,items])
}
facets_in <- data2FACETS(data)

# Example of FACETS
mod1.a <- immer::immer_FACETS(
  title="Example 1 with simulated data",
  convergence=NULL,
  totalscore="YES",
  facets=3,
  noncenter=NULL,
  arrange="m,N",
  entered_in_data="2,1,1",
  models="?$,?$,?$,R4",
  inter_rater=NULL,
  pt_biserial=NULL,
  faire_score="Zero",
  unexpected=2,
  usort=NULL,
  positive=1,
  labels=c("1,Persons","1-500","2,Rater","1-3","3,Item","1-4"),
  fileinput=NULL,
  data=facets_in,
  path.dosbox=NULL,
  path.facets="C:\\\\Facets",
  model.name="Example.SD",
  facetsEXE=NULL
  )

################################
# 2. Example on Windows using existing input-files of FACETS
################################
data(data.immer06)

mod1b <- immer::immer_FACETS(
  fileinput=data.immer06,
  path.facets="C:\\\\Facets",
  model.name="Example.SD",
  facetsEXE=NULL
  )
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{Hierarchical rater model}

% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
