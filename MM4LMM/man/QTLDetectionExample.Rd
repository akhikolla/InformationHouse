\name{QTLDetectionExample}
\alias{QTLDetectionExample}
\docType{data}
\title{
QTL Detection Example
%%   ~~ data name/kind ... ~~
}
\description{
This corresponds to (a sample of) the dataset presented in Rincent et al. (2014).
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("QTLDetectionExample")}
\format{
  The format is:
List of 3
\describe{
    \item{\code{Phenotype}}{Named num [1:259]}
\item{\code{Genotype}}{int [1:259,1:10]}
\item{\code{Kinship}}{num [1:259,1:259]}
}


}
\details{
The list contains three elements:
\itemize{
\item Phenotype: a numeric vector containing phenotypes of individuals
\item Genotype: a matrix containing the genotypes of indviduals over 10 biallelic markers
\item Kinship: a matrix of simple relatedness coefficients between individuals
}
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{\url{https://link.springer.com/article/10.1007\%2Fs00122-014-2379-7}
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Rincent, R., Nicolas, S., Bouchet, S., Altmann, T., Brunel, D., Revilla, P., ... & Schipprack, W. (2014). Dent and Flint maize diversity panels reveal important genetic potential for increasing biomass production. Theoretical and applied genetics, 127(11), 2313-2331.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(QTLDetectionExample)
names(QTLDetectionExample)
## maybe str(QTLDetectionExample) ; plot(QTLDetectionExample) ...
}
\keyword{datasets}
