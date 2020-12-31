\name{Laboratory_Data}
\docType{data}
\alias{Laboratory_Data}
\title{
Laboratory Data
}
\description{
The laboratory data is based on 459 patients who were listed
for heart transplant and were put on mechanical circulatory
system through device implantation from December 1991 to
July 2009 at Cleveland Clinic. These patients had periodic
measurements of their bilirubin and creatinine levels.
Data from 459 patients includes 18285 measurements
of bilirubin and creatinine with
an average of 39 measurements per patient.
}
\format{
  Laboratory data has 4 parts:
  \enumerate{
    \item  A total of 41 x-variables.
    \item  Time points (time).
    \item  Patient identifier (id).
    \item  Longitudinal responses (tbili_po and creat_po).
  }
}
\references{
Rajeswaran J., Blackstone E.H. and Bernard J.
Evolution of association between renal and liver function
while awaiting for the heart transplant: An application using
bivariate multiphase nonlinear mixed effect model.
\emph{Statistical methods in medical research} 27(7):2216--2230, 2018.
}
\examples{data(Laboratory_Data, package = "BoostMLR")}
\keyword{datasets}
