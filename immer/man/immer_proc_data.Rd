%% File Name: immer_proc_data.Rd
%% File Version: 0.26

\name{immer_proc_data}
\alias{immer_proc_data}
\alias{immer_create_design_matrix_formula}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Processing Datasets and Creating Design Matrices for Rating Data
}

\description{
The function \code{immer_proc_data} processes datasets containing rating data
into a dataset into a long format of pseudoitems (item \eqn{\times} raters).

The function \code{immer_create_design_matrix_formula} creates a design matrix
for a processed dataset and a provided formula.
}

\usage{
immer_proc_data(dat, pid=NULL, rater=NULL, weights=NULL, maxK=NULL)

immer_create_design_matrix_formula( itemtable, formulaA )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Datasets with integer item responses
}
  \item{pid}{
Vector with person identifiers
}
  \item{rater}{
Vector with rater identifiers
}
  \item{weights}{
Vector with sampling weights
}
\item{maxK}{Optional vector with maximum category per item}
\item{itemtable}{Processed item table. The table must include the column
    \code{item} (an integer item identifier) and \code{maxK} (maximum number of
    categories per item). Optional columns are \code{rater}
    (an integer rater identifier),\code{item_name} and \code{rater_name}.
}
\item{formulaA}{An \R formula. The facets \code{item}, \code{step} and
\code{rater} are treated as numeric. However, numeric transformation can
be applied for the step parameter by using the arguments \code{item_num},
\code{step_num} or \code{rater_num} in \code{formulaA}.
}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\value{
The output of \code{immer_proc_data} is a list with several entries (selection)

\item{dat2}{Dataset containing pseudoitems}
\item{dat2.resp}{Dataset containing response indicators for pseudoitems}
\item{dat2.NA}{Dataset containing pseudoitems and missing responses coded as \code{NA}}
\item{dat}{Original dataset}
\item{person.index}{Person identifiers}
\item{rater.index}{Rater identifiers}
\item{VV}{Number of items}
\item{N}{Number of persons}
\item{RR}{Number of raters}
\item{dat2.ind.resp}{Array containing indicators of pseudoitems and categories}
\item{ND}{Number of person-rater interactions}
\item{itemtable}{Information about processed data}

The output of \code{immer_create_design_matrix_formula} is a list with several
entries (selection)

\code{A}{design matrix}
\code{itemtable2}{Processed item table}
}


%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}


\examples{
#############################################################################
# EXAMPLE 1: Processing rating data
#############################################################################

data(data.immer01a, package="immer")
dat <- data.immer01a

res <- immer::immer_proc_data( dat=dat[,paste0("k",1:5)], pid=dat$idstud,
             rater=dat$rater)
str(res, max.level=1)

\dontrun{
#############################################################################
# EXAMPLE 2: Creating several design matrices for rating data
#############################################################################

data(data.ratings1, package="sirt")
dat <- data.ratings1
resp <- dat[,-c(1,2)]
#- redefine the second and third item such that the maximum category score is 2
for (vv in c(2,3)){
    resp[ resp[,vv] >=2,vv ] <- 2
}

#--- process data
res0 <- immer::immer_proc_data( dat=resp, pid=dat$idstud,  rater=dat$rater)

#--- rating scale model
des1 <- immer::immer_create_design_matrix_formula( itemtable=res0$itemtable,
                formulaA=~ item + step )
des1$des

#--- partial scale model
des2 <- immer::immer_create_design_matrix_formula( itemtable=res0$itemtable,
                formulaA=~ item + item:step )
des2$des

#--- multi-facets Rasch model
des3 <- immer::immer_create_design_matrix_formula( itemtable=res0$itemtable,
                formulaA=~ item + item:step + rater )
des3$des

#--- polytomous model with quadratic step effects
des4 <- immer::immer_create_design_matrix_formula( itemtable=res0$itemtable,
                formulaA=~ item + item:I(step_num^2) )
des4$des
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
