\name{interactiveClustering}
\alias{interactiveClustering}
\title{GUI for interactive cluster analysis}

\description{This tool is an interactive shiny tool that visualizes a given generalized Umatrix and allows the user to select areas and mark them as clusters to improve a projection based clustering.}

\arguments{
  \item{Umatrix}{[1:Lines,1:Columns] Matrix of Umatrix Heights}
  \item{Bestmaches}{[1:n,1:2]Array with positions of Bestmatches}
  \item{Cls}{[1:n]Classification of the Bestmatches}
  \item{Imx}{[1:4*Lines,1:4*Columns]Matrix of an island that will be cut out of the umatrix, use package Umatrix for generation.}
  \item{Toroid}{Are Bestmatches placed on a toroid grid? TRUE by default.}
}
\value{Cls[1:n]: A vector containing the selected class ids. The order is corresponding to the given Bestmatches }

\note{ If you want to verifiy your clustering result externally, you can use \code{Heatmap} or \code{SilhouettePlot} of the CRAN package \code{DataVisualizations}.
}
\details{
Clicking on "Quit" returns the Cls vector to the workspace.
}
\author{Florian Lerch, Michael Thrun}
\references{
[Thrun/Ultsch, 2017] Thrun, M.C., Ultsch, A.: Projection based Clustering, Conf. Int. Federation of Classification Societies (IFCS),DOI:10.13140/RG.2.2.13124.53124, Tokyo, 2017.

[Thrun, 2018]  Thrun, M. C.: Projection Based Clustering through Self-Organization and Swarm Intelligence, doctoral dissertation 2017, Springer, Heidelberg, ISBN: 978-3-658-20539-3, \url{https://doi.org/10.1007/978-3-658-20540-9}, 2018. 
}
\examples{
data('Hepta')
#2d projection
# Visualizuation of GeneralizedUmatrix
\donttest{
projectionpoints=NeRV(Hepta$Data)
#Computation of Generalized Umatrix
library(GeneralizedUmatrix)
visualization=GeneralizedUmatrix(Data = Hepta$Data,projectionpoints)
}
## Semi-Automatic Clustering done interactivly in a  shiny gui
\dontrun{
Cls = interactiveClustering(visualization$Umatrix, visualization$Bestmatches)
##Plotting
plotTopographicMap(visualization$Umatrix,visualization$Bestmatches,Cls)
}
}