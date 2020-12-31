%% File Name: immer_install.Rd
%% File Version: 0.12

\name{immer_install}
\alias{immer_install}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Support for the installation of the DOS-version from FACETS
}

\description{
This function supports the installation process of the DOS-version from FACETS
and also the necessary DOSBox in Windows, Linux (Ubuntu) and OS X
}

\usage{
immer_install(DosBox_path=NULL, Facets_path=NULL )


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DosBox_path}{
optional argument for the specification of the path where the DosBox should be saved}
  \item{Facets_path}{
optional argument for the specification of the path where FACETS should be saved
}

}

\details{
This function provides assistance for the installation process of the FACDOS
(DOS version of FACETS) and the required DosBox. Currently supported operating
systems are: Windows, Mac OS X and Ubuntu (Linux).
}

%\value{
%...need to be defined....
%}


\references{
Linacre, J. M. (1999). \emph{FACETS} (Version 3.17) [Computer software]. Chicago: MESA.

Veenstra, P., Froessman, T., Wohlers, U. (2015): \emph{DOSBox} (Version 0.74)
[Computer Software]. Arizona: Scottsdale.
}



%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Install FACDOS and DOSBox \code{\link{immer_FACETS}}.
}

\examples{
\dontrun{
    immer::immer_install( DosBox_path=NULL, Facets_path=NULL )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{install FACDOS}
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
