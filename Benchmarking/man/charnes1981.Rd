% $Id: charnes1981.Rd 96 2010-11-27 22:47:25Z Lars $
\name{charnes1981}
\alias{charnes1981}
\docType{data}

\title{ Data: Charnes et al. (1981): Program follow through }


\description{

 The data set is from an US federally sponsored program for providing
remedial assistance to disadvantaged primary school students. The
firms are 70 school sites, and data are from entire sites. The
variables consists of results from three different kind of tests, a
reading score, \code{y1}, a math score, \code{y2}, and a self--esteem
score, \code{y3}, which are considered outputs in the model, and five
different variables considered to be inputs, the education level of
the mother, \code{x1}, the highest occupation of a family member,
\code{x2}, parental visits to school, \code{x3}, time spent with
children in school-related topics, \code{x4}, and the number of
teachers at the site, \code{x5}. 

}


\usage{data(charnes1981) }

\format{
  A data frame with 70 school sites with the following variables.
  \describe{
    \item{\code{firm}}{school site number}
    \item{\code{x1}}{education level of the mother}
    \item{\code{x2}}{highest occupation of a family member}
    \item{\code{x3}}{parental visits to school}
    \item{\code{x4}}{time spent with children in school-related topics}
    \item{\code{x5}}{the number of teachers at the site}
    \item{\code{y1}}{reading score}
    \item{\code{y2}}{math score}
    \item{\code{y3}}{self--esteem score}
    \item{\code{pft}}{=1 if in program (program follow through) and 
                      =0 if not in program}
    \item{\code{name}}{Site name}
  }
}


\details{ The command \code{data(charnes1981)} will create a data
frame named \code{charnes1981} with the above data.

Beside input and output varianles there is further information in the
data set, that the first 50 school sites followed the program and that
the last 20 are the results for sites not following the program. This
is showed by the variable \code{pft}.  }


\note{ Data as \code{.csv} are loaded by the command \code{data} using
\code{ read.table(..., header=TRUE, sep=";")} such that this file
is a semicolon separated file and not a comma separated file.

Therefore, to read the file from a script the command must be
\code{read.csv("charnes1981.csv", sep=";")} or
\code{read.csv2("charnes1981.csv")}.

Thus the data can be read either as \code{charnes1981 <-}
\cr\code{read.csv2(paste(.Library, "Benchmarking/data",
"charnes1981.csv",  sep ="/"))} \cr or as \code{data(charnes1981)} if
the package \pkg{Benchmarking} is loaded. In both cases the data will
be in the data frame \code{charnes1981}.

}


\source{ Charnes, Cooper, and Rhodes, \dQuote{Evaluating Program and
Managerial Efficiency: An Application of Data Envelopment Analysis to
Program Follow Through}, \emph{Management Science}, volume 27, number
6, June 1981, pages 668--697.  }


\examples{
data(charnes1981)
x <- with(charnes1981, cbind(x1,x2,x3,x4,x5))
y <- with(charnes1981, cbind(y1,y2,y3))

# Farrell inpout efficiency; vrs technology
e <- dea(x,y)
# The number of times each peer is a peer
np <- get.number.peers(e) 
# Peers that are peers for more than 20 schools, and the number of
# times they are peers
np[which(np[,2]>20),]

# Plot first input against first output and emphasize the peers that
# are peers for more than 20 schools in the model with five inputs and
# three outputs
inp <- np[which(np[,2]>20),1]
dea.plot(x[,1],y[,1])
points(x[inp,1], y[inp,1], pch=16, col="red")
}

\keyword{datasets}
