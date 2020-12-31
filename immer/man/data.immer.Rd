%% File Name: data.immer.Rd
%% File Version: 0.38

\name{data.immer}

\alias{data.immer01a}
\alias{data.immer01b}
\alias{data.immer02}
\alias{data.immer03}
\alias{data.immer04a}
\alias{data.immer04b}
\alias{data.immer05}
\alias{data.immer06}
\alias{data.immer07}
\alias{data.immer08}
\alias{data.immer09}
\alias{data.immer10}
\alias{data.immer11}
\alias{data.immer12}

\docType{data}


\title{
Some Example Datasets for the \pkg{immer} Package
}


\description{
Some example rating datasets for the \pkg{immer} package.
}


\usage{
data(data.immer01a)
data(data.immer01b)
data(data.immer02)
data(data.immer03)
data(data.immer04a)
data(data.immer04b)
data(data.immer05)
data(data.immer06)
data(data.immer07)
data(data.immer08)
data(data.immer09)
data(data.immer10)
data(data.immer11)
data(data.immer12)
}


\format{

\itemize{

\item The format of the dataset \code{data.immer01a} is:

\code{'data.frame':   23904 obs. of  8 variables:} \cr
\code{ $ idstud: int  10001 10001 10003 10003 10003 10004 10004 10005 10005 10006 ...} \cr
\code{ $ type  : Factor w/ 2 levels "E","I": 1 2 1 1 2 1 2 1 2 1 ...} \cr
\code{ $ rater : Factor w/ 57 levels "R101","R102",..: 1 36 33 20 21 57 36 9 31 21 ...} \cr
\code{ $ k1    : int  2 1 0 0 0 2 2 1 2 0 ...} \cr
\code{ $ k2    : int  1 1 0 0 0 1 1 1 2 0 ...} \cr
\code{ $ k3    : int  1 1 0 0 0 1 1 1 2 1 ...} \cr
\code{ $ k4    : int  2 2 1 0 0 1 1 1 2 1 ...} \cr
\code{ $ k5    : int  1 2 0 0 0 2 1 2 3 2 ...} \cr

\item The format of the dataset \code{data.immer01b} is:

\code{'data.frame':   4244 obs. of  8 variables:} \cr
\code{ $ idstud: int  10001 10003 10005 10007 10009 10016 10018 10022 10024 10029 ...} \cr
\code{ $ type  : Factor w/ 1 level "E": 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ rater : Factor w/ 20 levels "R101","R102",..: 1 20 9 5 14 19 20 6 10 10 ...} \cr
\code{ $ k1    : int  2 0 1 2 2 2 3 1 3 2 ...} \cr
\code{ $ k2    : int  1 0 1 2 2 1 3 2 2 1 ...} \cr
\code{ $ k3    : int  1 0 1 1 3 2 2 1 3 1 ...} \cr
\code{ $ k4    : int  2 0 1 2 3 2 2 2 3 2 ...} \cr
\code{ $ k5    : int  1 0 2 1 3 1 2 3 3 1 ...}

This dataset is a subset of \code{data.immer01a}. \cr



\item The format of the dataset \code{data.immer02} is:

\code{'data.frame':   6105 obs. of  6 variables:} \cr
\code{ $ idstud: int  10002 10004 10005 10006 10007 10008 10009 10010 10013 10014 ...} \cr
\code{ $ rater : Factor w/ 44 levels "DR101","DR102",..: 43 15 12 21 9 3 35 24 11 17 ...} \cr
\code{ $ a1    : int  3 1 2 1 0 2 1 2 1 1 ...} \cr
\code{ $ a2    : int  3 0 3 1 0 3 0 2 2 1 ...} \cr
\code{ $ a3    : int  1 2 0 1 2 3 2 2 1 1 ...} \cr
\code{ $ a4    : int  2 1 2 1 1 3 1 2 2 1 ...} \cr




\item The format of the dataset \code{data.immer03} is:

\code{'data.frame':   6466 obs. of  6 variables:} \cr
\code{ $ idstud: int  10001 10002 10003 10004 10005 10006 10007 10009 10010 10012 ...} \cr
\code{ $ rater : Factor w/ 44 levels "R101","R102",..: 18 10 8 25 19 31 16 22 29 6 ...} \cr
\code{ $ b1    : int  1 2 1 3 3 2 3 2 2 1 ...} \cr
\code{ $ b2    : int  2 1 0 3 3 1 1 2 2 1 ...} \cr
\code{ $ b3    : int  2 3 1 2 3 1 2 2 2 2 ...} \cr
\code{ $ b4    : int  1 2 0 2 2 2 3 2 3 1 ...} \cr


\item The format of the dataset \code{data.immer04a} is:

\code{'data.frame':   25578 obs. of  7 variables:} \cr
\code{ $ idstud: int  10001 10001 10001 10002 10002 10002 10003 10003 10004 10004 ...} \cr
\code{ $ task  : Factor w/ 4 levels "l1","l2","s1",..: 1 4 4 1 1 3 1 3 2 2 ...} \cr
\code{ $ rater : Factor w/ 43 levels "R101","R102",..: 14 31 25 39 35 19 43 27 12 4 ...} \cr
\code{ $ TA    : int  5 2 4 0 0 0 2 6 5 3 ...} \cr
\code{ $ CC    : int  4 1 3 1 0 0 2 6 4 3 ...} \cr
\code{ $ GR    : int  4 1 2 1 0 0 1 7 5 2 ...} \cr
\code{ $ VOC   : int  4 2 3 1 0 0 1 6 5 3 ...} \cr



\item The format of the dataset \code{data.immer04b} is:

\code{'data.frame':   2975 obs. of  7 variables:} \cr
\code{ $ idstud: int  10002 10004 10010 10013 10015 10016 10024 10025 10027 10033 ...} \cr
\code{ $ task  : Factor w/ 1 level "s1": 1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ rater : Factor w/ 20 levels "R101","R102",..: 19 1 5 16 13 13 8 10 19 5 ...} \cr
\code{ $ TA    : int  0 3 5 5 3 2 3 6 4 5 ...} \cr
\code{ $ CC    : int  0 3 4 5 4 1 4 7 3 3 ...} \cr
\code{ $ GR    : int  0 3 3 6 5 2 3 6 3 2 ...} \cr
\code{ $ VOC   : int  0 2 4 6 5 2 3 6 3 2 ...}

This dataset is a subset of \code{data.immer04a}. \cr



\item The format of the dataset \code{data.immer05} is:

\code{'data.frame':   21398 obs. of  9 variables:} \cr
\code{ $ idstud     : int  10001 10001 10002 10002 10003 10003 10004 10004 10005 10005 ...} \cr
\code{ $ type       : Factor w/ 2 levels "l","s": 2 1 2 1 2 1 2 1 2 1 ...} \cr
\code{ $ task       : Factor w/ 6 levels "l1","l4","l5",..: 5 2 6 3 5 1 5 1 5 2 ...} \cr
\code{ $ rater      : Factor w/ 41 levels "ER101","ER102",..: 1 40 38 23 37 33 2 33 21 27 ...} \cr
\code{ $ idstud_task: Factor w/ 19484 levels "10001l4","10001s3",..: 2 1 4 3 6 5 8 7 10 9 ...} \cr
\code{ $ TA         : int  3 4 6 6 4 2 0 3 1 3 ...} \cr
\code{ $ CC         : int  5 4 5 5 3 3 0 2 5 3 ...} \cr
\code{ $ GR         : int  4 4 5 6 5 3 0 4 5 4 ...} \cr
\code{ $ VO         : int  6 4 6 6 4 3 0 3 4 3 ...} \cr

\item The dataset \code{data.immer06} is a string containing
an input syntax for the FACETS program.

\item The format of the dataset \code{data.immer07} is:

\code{'data.frame':   1500 obs. of  6 variables:} \cr
\code{ $ pid  : int  1 1 1 2 2 2 3 3 3 4 ...} \cr
\code{ $ rater: chr  "R1" "R2" "R3" "R1" ...} \cr
\code{ $ I1   : num  1 1 2 1 1 1 0 1 1 2 ...} \cr
\code{ $ I2   : num  0 1 1 2 1 2 1 1 2 1 ...} \cr
\code{ $ I3   : num  1 1 2 0 0 1 1 0 2 1 ...} \cr
\code{ $ I4   : num  0 0 1 0 0 1 0 1 2 0 ...} \cr



\item The format of the dataset \code{data.immer08} (example in
Schuster & Smith, 2006) is

\code{'data.frame':   16 obs. of  3 variables:} \cr
\code{ $ Facility: int  1 1 1 1 2 2 2 2 3 3 ...} \cr
\code{ $ Research: int  1 2 3 4 1 2 3 4 1 2 ...} \cr
\code{ $ weights : int  40 6 4 15 4 25 1 5 4 2 ...} \cr


\item The dataset \code{data.immer09} contains reviewer ratings for
conference papers (Kuhlisch et al., 2016):

\code{'data.frame':   128 obs. of  3 variables:} \cr
\code{ $ idpaper   : int  1 1 1 2 2 3 3 3 4 4 ...} \cr
\code{ $ idreviewer: int  11 15 20 1 10 11 15 20 13 16 ...} \cr
\code{ $ score     : num  7 7 7 7 7 7 7 7 7 7 ...} \cr

\item The dataset \code{data.immer10} contains standard setting
ratings of 13 raters on 61 items (including item identifier \code{item}
and item difficulty \code{itemdiff})

\code{'data.frame':   61 obs. of  15 variables:} \cr
\code{ $ item    : chr  "I01" "I02" "I03" "I04" ...} \cr
\code{ $ itemdiff: num  380 388 397 400 416 425 427 434 446 459 ...} \cr
\code{ $ R01     : int  1 3 2 2 1 3 2 2 3 1 ...} \cr
\code{ $ R02     : int  1 1 1 1 1 2 1 2 2 1 ...} \cr
\code{ $ R03     : int  1 1 1 1 1 1 2 2 3 1 ...} \cr
\code{ $ R04     : int  1 2 1 3 2 2 2 2 3 2 ...} \cr
\code{ $ R05     : int  1 1 2 1 1 1 2 2 3 2 ...} \cr
\code{ $ R06     : int  1 2 1 1 1 2 2 2 3 2 ...} \cr
\code{ $ R07     : int  1 2 1 2 1 1 2 1 3 1 ...} \cr
\code{ $ R08     : int  2 2 1 2 1 1 2 2 3 2 ...} \cr
\code{ $ R09     : int  2 1 1 2 1 2 1 2 3 1 ...} \cr
\code{ $ R10     : int  2 2 2 2 1 2 2 3 3 2 ...} \cr
\code{ $ R11     : int  2 2 1 2 1 2 2 2 3 2 ...} \cr
\code{ $ R12     : int  2 2 1 3 1 2 2 2 3 2 ...} \cr
\code{ $ R13     : int  1 1 1 1 1 1 1 1 2 1 ...} \cr

\item The dataset \code{data.immer11} contains ratings of 148 cases (screening
mammogram samples) diagnoses by 110 raters (Zhang & Petersen, xxxx). The codes of the
polytomous rating are normal (code 0), benign (code 1), probably benign (code 2),
possibly malignant (code 3), and probably malignant (code 4). The dataset was extracted
from an image plot in Figure 2 by using the processing function
\code{png::readPNG}. The format of the dataset is

\code{'data.frame':   148 obs. of  110 variables:} \cr
\code{ $ R001: num  2 1 3 2 1 2 0 0 0 2 ...} \cr
\code{ $ R002: num  1 3 4 4 0 4 0 0 3 0 ...} \cr
\code{ $ R003: num  0 0 0 4 0 2 3 0 0 0 ...} \cr
\code{ $ R004: num  1 2 1 4 2 2 2 0 4 4 ...} \cr
\code{ [... ]} \cr


\item The dataset \code{data.immer12} contains ratings of the 2002 olympic pairs
figure skating competition. This dataset has been used in Lincare (2009).  The items
are \code{ST} (short program, technical merit),
\code{SA} (short program, artistic impression),
\code{FT} (free program, technical merit), and
\code{FA} (free program, artistic impression). The format of the dataset is

\code{'data.frame':   180 obs. of  7 variables:} \cr
\code{ $ idpair: int  1 1 1 1 1 1 1 1 1 2 ...} \cr
\code{ $ pair  : chr  "BB-Svk" "BB-Svk" "BB-Svk" "BB-Svk" ...} \cr
\code{ $ judge : chr  "RUS" "CHI" "USA" "FRA" ...} \cr
\code{ $ ST    : int  58 57 57 56 55 55 50 51 51 47 ...} \cr
\code{ $ SA    : int  58 57 57 56 55 55 50 51 51 47 ...} \cr
\code{ $ FT    : int  58 57 57 56 55 55 50 51 51 47 ...} \cr
\code{ $ FA    : int  58 57 57 56 55 55 50 51 51 47 ...} \cr

}
}





%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}


%\seealso{
%Some free datasets can be obtained from \cr
%}

%\source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
%}



\references{
Kuhlisch, W., Roos, M., Rothe, J., Rudolph, J., Scheuermann, B., &
Stoyan, D. (2016). A statistical approach to calibrating the scores of
biased reviewers of scientific papers.
\emph{Metrika, 79}, 37-57.

Linacre, J. M. (2009). Local independence and residual covariance: A study of
Olympic figure skating ratings. \emph{Journal of Applied Measurement, 10}(2), 157-169.

Schuster, C., & Smith, D. A. (2006). Estimating with a latent
class model the reliability of nominal judgments upon which
two raters agree. \emph{Educational and Psychological Measurement, 66}(5), 739-747.

Zhang, S., & Petersen, J. H. (XXXX). Quantifying rater variation for ordinal data using
a rating scale model. \emph{Statistics in Medicine, XX}(xx), xxx-xxx.
}


%\examples{
%# xx
% }


\keyword{datasets}
