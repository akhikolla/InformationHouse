\name{tankPerformance}
\alias{tankPerformance}
\title{
A function to calculate difference performance from simulated tank behaviour
}
\description{
A function to calculate difference performance from simulated tank behaviour
}
\usage{tankPerformance(data=NULL,
                      roofArea=50,   
                      nPeople=1,
                      tankVol=3000,
		      firstFlush=1,
                      write.file=TRUE,
                      fnam="tankperformance.csv")
}
\arguments{
  \item{data}{A dataframe of observed climate data in the form \emph{Year Month Day P Temp}.}
  \item{roofArea}{ roof area in m2}
  \item{nPeople}{number of people using water}
  \item{tankVol}{tank volume in L}
  \item{firstFlush}{first flush depth over roof in mm}
  \item{write.file}{logical. write output tank timeseries to file T/F?}
  \item{fnam}{string indicating name of file}
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
%\examples{
%data(modelTag)
%## maybe str(modelTag) ; plot(modelTag) ...
%}