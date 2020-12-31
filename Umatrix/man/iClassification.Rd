\name{iClassification}
\alias{iClassification}
\title{GUI for manual classification}

\description{This tool is a 'shiny' GUI that visualizes a given Umatrix and allows the user to select areas and mark them as clusters.}

\arguments{
  \item{Umatrix}{Matrix of Umatrix Heights}
  \item{BestMatches}{Array with positions of Bestmatches}
  \item{Cls}{Classification of the Bestmatches}
  \item{Imx}{Matrix of an island that will be cut out of the Umatrix}
  \item{Toroid}{Are BestMatches placed on a toroid grid? TRUE by default}
}
\value{A vector containing the selected class ids. The order is corresponding to the given Bestmatches }

\examples{\dontrun{
data("Hepta")
e = esomTrain(Hepta$Data, Key = 1:nrow(Hepta$Data))
cls = iClassification(e$Umatrix, e$BestMatches)
}}

\references{
Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016.
}
