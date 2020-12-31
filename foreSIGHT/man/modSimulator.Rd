\name{modSimulator}
\alias{modSimulator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
modSimulator
}
\description{
Simulates using weather generator models specified using modelTag.
}
\usage{
modSimulator(datStart=NULL,
             datFinish=NULL,
             modelTag=NULL,      
             parS=NULL,
             seed=NULL,
             file=NULL,
             IOmode="suppress"
             )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datStart}{A date string in an accepted date format e.g. \emph{"01-10-1990"}.}
  \item{datFinish}{A date string in an accepted date format e.g. \emph{"01-10-1990"}. Must occur after datStart.}
  \item{modelTag}{A character vector of which stochastic models to use to create each climate variable. Supported tags are shown in details below.}
  \item{parS}{A list (names must match supplied modelTags) containing numeric vectors of model parameters.} 
  \item{seed}{Numeric. Seed value supplied to weather generator.}
  \item{file}{Character. Specifies filename for simulation output.}
  \item{IOmode}{A string that specifies the input-output mode for the time series = "verbose", "dev" or "suppress".} 
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
\dontrun{
data(tankDat); obs=tank_obs                     #Get observed data
modelTag=c("P-har12-wgen","Temp-har26-wgen")    #Select models
pars=modCalibrator(obs=obs,modelTag=modelTag)   #Calibrate models
sim=modSimulator(datStart="1970-01-01",         #Simulate!
                 datFinish="1999-12-31",
                 modelTag=modelTag,
                 parS=pars,
                 seed=123,
                 file=paste0("tester.csv"),
                 IOmode="verbose")
plot(sim$P[1:365])                             #Plot first year of rainfall
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}% use one of  RShowDoc("KEYWORDS")
%\keyword{  }% __ONLY ONE__ keyword per line
