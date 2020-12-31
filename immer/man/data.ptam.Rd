%% File Name: data.ptam.Rd
%% File Version: 0.19

\name{data.ptam}

\alias{data.ptam1}
\alias{data.ptam2}
\alias{data.ptam3}
\alias{data.ptam4}
\alias{data.ptam4long}
\alias{data.ptam4wide}

\docType{data}


\title{
Example Datasets for Robitzsch and Steinfeld (2018)
}


\description{
Example datasets for Robitzsch and Steinfeld (2018).
}


\usage{
data(data.ptam1)
data(data.ptam2)
data(data.ptam3)
data(data.ptam4)
data(data.ptam4long)
data(data.ptam4wide)
}


\format{

\itemize{

\item The dataset \code{data.ptam1} is a subset of the dataset from Example 3
of the ConQuest manual and
contains 9395 ratings for 6877 students and 9 raters
on 2 items (\code{OP} and \code{TF}). The format is

\code{'data.frame':   9395 obs. of  4 variables:} \cr
\code{ $ pid  : int  1508 1564 1565 1566 1567 1568 1569 1629 1630 1631 ...} \cr
\code{ $ rater: num  174 124 124 124 124 124 124 114 114 114 ...} \cr
\code{ $ OP   : int  2 1 2 1 1 1 2 2 2 3 ...} \cr
\code{ $ TF   : int  3 1 2 2 1 1 2 2 2 3 ...} \cr


\item The dataset \code{data.ptam2} contains 1043 ratings for
262 students and 17 raters on 19 items (\code{A1}, ..., \code{D9}). The format is

\code{'data.frame':   1043 obs. of  21 variables:} \cr
\code{ $ idstud : int  1001 1001 1001 1001 1002 1002 1002 1002 1003 1003 ...} \cr
\code{ $ idrater: int  101 108 212 215 104 108 209 211 103 104 ...} \cr
\code{ $ A1     : int  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ A2     : int  1 1 1 1 0 0 0 1 1 1 ...} \cr
\code{ $ A3     : int  1 1 1 1 1 1 0 1 0 0 ...} \cr
\code{ [...]} \cr
\code{ $ D9     : int  2 2 2 2 2 2 2 2 1 0 ...} \cr


\item The dataset \code{data.ptam3} contains 523 ratings for
262 students and 8 raters on 23 items (\code{A1}, ..., \code{J0}). The format is

\code{'data.frame':   523 obs. of  25 variables:} \cr
\code{ $ idstud : int  1001 1001 1002 1002 1003 1003 1004 1004 1005 1005 ...} \cr
\code{ $ idrater: int  101 108 104 108 103 104 102 104 102 108 ...} \cr
\code{ $ A1     : int  1 1 1 1 1 1 1 1 1 1 ...} \cr
\code{ $ A2     : int  1 1 0 0 1 1 NA 0 1 1 ...} \cr
\code{ $ A3     : int  1 1 1 1 0 0 0 0 0 0 ...} \cr
\code{ [...]} \cr
\code{ $ J0     : int  2 3 3 2 0 0 2 2 0 1 ...} \cr

\item The dataset \code{data.ptam4} contains 592 ratings for
209 students and 10 raters on 3 items (\code{crit2}, \code{crit3} and
\code{crit4}). The format is

\code{'data.frame':   592 obs. of  5 variables:} \cr
\code{ $ idstud: num  10005 10009 10010 10010 10014 ...} \cr
\code{ $ rater : num  802 802 844 802 837 824 820 803 816 844 ...} \cr
\code{ $ crit2 : int  3 2 0 2 1 0 2 1 1 0 ...} \cr
\code{ $ crit3 : int  3 2 1 2 2 2 2 2 2 2 ...} \cr
\code{ $ crit4 : int  2 1 2 1 2 2 2 2 2 2 ...} \cr


\item The dataset \code{data.ptam4long} is the dataset \code{data.ptam4} which has been
converted into a long format for analysis with mixed effects models in the
\pkg{lme4} package. The format is

\code{'data.frame':   1776 obs. of  17 variables:} \cr
\code{ $ idstud : num  10005 10005 10005 10009 10009 ...} \cr
\code{ $ rater  : num  802 802 802 802 802 802 844 802 844 802 ...} \cr
\code{ $ item   : Factor w/ 3 levels "crit2","crit3",..: 1 2 3 1 2 3 1 1 2 2 ...} \cr
\code{ $ value  : int  3 3 2 2 2 1 0 2 1 2 ...} \cr
\code{ $ I_crit2: num  1 0 0 1 0 0 1 1 0 0 ...} \cr
\code{ $ I_crit3: num  0 1 0 0 1 0 0 0 1 1 ...} \cr
\code{ $ I_crit4: num  0 0 1 0 0 1 0 0 0 0 ...} \cr
\code{ $ R_802  : num  1 1 1 1 1 1 0 1 0 1 ...} \cr
\code{ $ R_803  : num  0 0 0 0 0 0 0 0 0 0 ...} \cr
\code{ [...] } \cr
\code{ $ R_844  : num  0 0 0 0 0 0 1 0 1 0 ...} \cr

\item The dataset \code{data.ptam4wide} contains multiple ratings of
40 students from the dataset \code{data.ptam4} from the item \code{crit2}.
Each column corresponds to one rater. The format is

\code{'data.frame':   40 obs. of  11 variables:} \cr
\code{ $ pid : chr  "10014" "10085" "10097" "10186" ...} \cr
\code{ $ R802: int  2 3 2 2 2 1 1 2 2 2 ...} \cr
\code{ $ R803: int  1 1 3 1 2 0 0 0 1 0 ...} \cr
\code{ $ R810: int  1 2 2 2 1 0 1 1 2 1 ...} \cr
\code{ $ R816: int  1 2 3 2 2 0 1 1 2 1 ...} \cr
\code{ $ R820: int  2 2 2 2 1 1 1 1 1 1 ...} \cr
\code{ $ R824: int  0 3 2 3 2 0 0 1 2 1 ...} \cr
\code{ $ R831: int  1 2 2 2 1 0 0 0 1 1 ...} \cr
\code{ $ R835: int  0 1 2 2 1 1 0 0 2 1 ...} \cr
\code{ $ R837: int  1 2 3 2 2 0 1 1 2 2 ...} \cr
\code{ $ R844: int  0 2 3 2 2 0 0 0 1 3 ...} \cr


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
Robitzsch, A., & Steinfeld, J. (2018). Item response models for human ratings: Overview,
estimation methods, and implementation in R.
\emph{Psychological Test and Assessment Modeling, 60}(1), 101-139.
}


%\examples{
%# xx
% }


\keyword{datasets}
