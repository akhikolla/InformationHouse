\name{iEsomTrain}
\alias{iEsomTrain}
\title{iEsomTrain}
\description{Trains the ESOM and shows the Umatrix.}
\arguments{
  \item{Data}{Matrix of Data that will be used to learn. One DataPoint per row}
  \item{BestMatches}{Array with positions of Bestmatches}
  \item{Cls}{Classification of the Bestmatches as a  vector}
  \item{Key}{Numeric vector of keys matching the Bestmatches}
  \item{Toroid}{Are BestMatches placed on a toroid grid? TRUE by default}
}
\value{
List with
  \item{Umatrix}{matrix with height values of the umatrix}
  \item{BestMatches}{matrix containing the bestmatches}
  \item{Lines}{number of lines of the chosen ESOM}
  \item{Columns}{number of columns of the chosen ESOM}
  \item{Epochs}{number of epochs of the chosen ESOM}
  \item{Weights}{List of weights}
  \item{Toroid}{True if a toroid grid was used}
  \item{EsomDetails}{Further details describing the chosen ESOM parameters}
}

\references{
Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016.
}
