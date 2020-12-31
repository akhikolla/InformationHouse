\name{modCalibrator}
\alias{modCalibrator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
modCalibrator
}
\description{
Calibrates weather generator models specified using modelTag.
}
\usage{
modCalibrator(obs = NULL,
              modelTag = NULL,
              window = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{A dataframe of observed climate data in the form \emph{Year Month Day P Temp}.}
  \item{modelTag}{A character vector of which stochastic models to use to create each climate variable. Supported tags are shown in under details below.}
  \item{window}{moving average window to calibrate daily gamma parameters for the modelTag \code{"P-har-WGEN"}.}
}
\details{
modelTag provides the main function with requested models. modelTag is vector of any of the following supported models:
  \itemize{
    \item \code{"Simple-ann"}{a simple annual scaling}
    \item \code{"P-ann-wgen"}{a four parameter annual rainfall model}
    \item \code{"P-seas-wgen"}{a 16 parameter seasonal rainfall model}
    \item \code{"P-har6-wgen"}{a harmonic rainfall model with 6 periods}
    \item \code{"P-har12-wgen"}{a harmonic rainfall model}
    \item \code{"P-har12-wgen-FS"}{a harmonic rainfall model where seasonality is fixed (phase angles must be specified via modelInfoMod=list("P-har12-wgen-FS"=fixedPars=c(x,x,x,x))}
    \item \code{"P-har26-wgen"}{a harmonic rainfall model}
    \item \code{"P-2har26-wgen"}{a double harmonic rainfall model}
    \item \code{"Temp-har26-wgen"}{a harmonic temperature model not conditional on rainfall}
    \item \code{"Temp-har26-wgen-wd"}{a harmonic temperature model dependent on wet or dry day}
    \item \code{"Temp-har26-wgen-wdsd"}{a harmonic temperature model where standard deviation parameters are dependent on wet or dry day}
    \item \code{"PET-har12-wgen"}{a harmonic potential evapotranspiration model}
    \item \code{"PET-har26-wgen"}{a harmonic potential evapotranspiration model}
    \item \code{"PET-har26-wgen-wd"}{a harmonic potential evapotranspiration model dependent on wet or dry day}
    \item \code{"Radn-har26-wgen"}{a harmonic solar radiation model (MJ/m2)}
  }
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
data(tankDat)                               #Load tank data (tank_obs)
modelTag=c("P-ann-wgen","Temp-har26-wgen")  #Select a rainfall and a temperature generator
out<- modCalibrator(obs = tank_obs,         #Calibrate models
                    modelTag = modelTag)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
%\keyword{  }% __ONLY ONE__ keyword per line
