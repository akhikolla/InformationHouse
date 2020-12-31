\name{iUmapIsland}
\alias{iUmapIsland}
\title{iUmapIsland}

\description{The toroid Umatrix is usually drawn 4 times, so that connected areas on borders can be seen as a whole. An island is a manual cutout of such a tiled visualization, that is selected such that all connected areas stay intact. This 'shiny' tool allows the user to do this manually.}

\arguments{
  \item{Umatrix}{Matrix of Umatrix Heights}
  \item{BestMatches}{Array with positions of BestMatches}
  \item{Cls}{Classification of the BestMatches}
}
\value{Boolean Matrix that represents the island within the tiled Umatrix}

\examples{\dontrun{
data("Hepta")
e = esomTrain(Hepta$Data, Key = 1:nrow(Hepta$Data))
Imx = iUmapIsland(e$Umatrix, e$BestMatches)
plotMatrix(e$Umatrix, e$BestMatches, Imx = Imx$Imx)
}}

\references{
Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016.
}
