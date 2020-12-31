\name{iUstarmatrix}
\alias{iUstarmatrix}
\title{iUstarmatrix}
\description{Calculates the Ustarmatrix by combining a Umatrix with a Pmatrix.}
\arguments{
  \item{Weights}{Weights that were trained by the ESOM algorithm}
  \item{Lines}{Height of the used grid}
  \item{Columns}{Width of the used grid}
  \item{Data}{Matrix of Data that was used to train the ESOM. One datapoint per row}
  \item{Imx}{Island mask that will be cut out from displayed Umatrix}
  \item{Cls}{Classification of the Bestmatches}
  \item{Toroid}{Are weights placed on a toroid grid?}
}
\value{
  \item{Ustarmatrix}{matrix with height values of the Ustarmatrix}
}

\references{
Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016.
}
