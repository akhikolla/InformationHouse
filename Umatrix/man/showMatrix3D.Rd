\name{showMatrix3D}
\alias{showMatrix3D}
\title{showMatrix3D}

\arguments{
 \item{Matrix}{Matrix to be plotted}
\item{BestMatches}{Positions of BestMatches to be plotted onto the matrix}
  \item{Cls}{Class identifier for the BestMatch at the given point}
 \item{Imx}{a mask (island) that will be used to cut out the Umatrix}
 \item{Toroid}{Should the Matrix be drawn 4 times (in a toroid view)}
 \item{HeightScale}{Optional. Scaling Factor for Mountain Height}
 \item{BmSize}{Size of drawn BestMatches}
 \item{RemoveOcean}{Remove as much area sourrounding an island as possible}
 \item{ColorStyle}{Either "Umatrix" or "Pmatrix" respectevily for their colors}
 \item{ShowAxis}{Draw an axis arround the drawn matrix}
 \item{SmoothSlope}{Try to increase the island size, to get smooth slopes around the island}
 \item{ClsColors}{Vector of colors that will be used for classes}
 \item{FileName}{Name for a stl file to write the Matrix to}
}

\description{Visualizes the matrix(Umatrix/Pmatrix) in an interactive window in 3D.}

\details{
The heightScale is set at the proportion of the 1 percent quantile against the 99 percent quantile of the Matrix values.
}

\examples{\dontrun{
data("Hepta")
e = esomTrain(Hepta$Data, Key = 1:nrow(Hepta$Data))
showMatrix3D(e$Umatrix)
}}

\references{
Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016.
}
