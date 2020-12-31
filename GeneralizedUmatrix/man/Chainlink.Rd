\name{Chainlink}
\alias{Chainlink}

\docType{data}
\title{
Chainlink is part of the Fundamental Clustering Problem Suit (FCPS) [Thrun/Ultsch, 2020].
}
\description{
linear not separable dataset of two interwined chains.
}
\usage{data("Chainlink")}

\details{
Size 1000, Dimensions 3, stored in \code{Chainlink$Data}

Teo clusters, stored in \code{Chainlink$Cls}

Published in [Ultsch et al.,1994] in German and [Ultsch 1995] in English.
}

\references{
[Thrun/Ultsch, 2020]	Thrun, M. C., & Ultsch, A.: Clustering Benchmark Datasets Exploiting the Fundamental Clustering Problems, Data in Brief,Vol. 30(C), pp. 105501, DOI 10.1016/j.dib.2020.105501 , 2020. 

[Ultsch 1995]	Ultsch, A.: Self organizing neural networks perform different from statistical k-means clustering,
Proc. Society for Information and Classification (GFKL), Vol. 1995, Basel 8th-10th March, 1995.

[Ultsch et al.,1994]	Ultsch, A., Guimaraes, G., Korus, D., & Li, H.: Knowledge extraction from artificial neural
networks and applications, Parallele Datenverarbeitung mit dem Transputer, pp. 148-16Chainlink, Springer, 1994.
}
\examples{
data(Chainlink)
str(Chainlink)

\dontrun{
require(DataVisualizations)
DataVisualizations::Plot3D(Chainlink$Data,Chainlink$Cls)
}
}
\keyword{datasets}
\keyword{Chainlink}
\keyword{FCPS}