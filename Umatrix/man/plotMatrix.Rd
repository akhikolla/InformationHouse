\name{plotMatrix}
\alias{plotMatrix}
\title{plotMatrix}
\arguments{
\item{Matrix}{Umatrix or Pmatrix to be plotted}
\item{BestMatches}{Positions of BestmMtches to be plotted onto the Umatrix}
\item{Cls}{Class identifier for the BestMatches}
\item{ClsColors}{Vector of colors that will be used to colorize the different classes}
\item{ColorStyle}{If "Umatrix" the colors of a Umatrix (Blue -> Green -> Brown -> White) will be used; If "Pmatrix" the colors of a Pmatrix (White -> Yellow -> Red) will be used}
\item{Toroid}{Should the Umatrix be drawn 4times?}
\item{BmSize}{Integer between 0.1 and 5, magnification factor of the drawn BestMatch circles}
\item{DrawLegend}{If TRUE, a color legend will be drawn next to the plot}
\item{FixedRatio}{If TRUE, the plot will be drawn with a fixed ratio of x and y axis}
\item{CutoutPol}{Only draws the area within given polygon}
\item{Nrlevels}{Number of height levels that will be used within the Umatrix}
\item{TransparentContours}{Use half transparent contours. Looks better but is slow}
\item{Imx}{Mask to cut out an island. Every value should be either 1 (stays in) or 0 (gets cut out)}
\item{Clean}{If TRUE axis, margins, ... surrounding the Umatrix image will be removed}
\item{RemoveOcean}{If TRUE, the surrounding blue area around an island will be reduced as much as possible (while still maintaining a rectangular form)}
\item{TransparentOcean}{If TRUE, the surrounding blue area around an island will be transparent}
\item{Title}{A title that will be drawn above the plot}
\item{BestMatchesLabels}{Vector of strings corresponding to the order of BestMatches which will be drawn on the plot as labels}
\item{BestMatchesShape}{Numeric value of Shape that will be used. Responds to the usual shapes of ggplot}
\item{MarkDuplicatedBestMatches}{If TRUE, BestMatches that are shown more than once within an island, will be marked}
\item{YellowCircle}{If TRUE, a yellow circle is drawn around Bestmatches to distinct them better from background}
}

\value{A 'ggplot' of a Matrix}
\details{The heightScale (nrlevels) is set at the proportion of the 1 percent quantile against the 99 percent quantile of the matrix values.}

\description{
 Draws a plot based of given Umatrix or Pmatrix.}

\references{
Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016.

Ultsch, A.: Maps for the visualization of high-dimensional data spaces, Proc. Workshop on Self organizing Maps (WSOM), pp. 225-230, Kyushu, Japan, 2003.

Siemon, H.P., Ultsch,A.: Kohonen Networks on Transputers: Implementation and Animation, in: Proceedings Intern. Neural Networks, Kluwer Academic Press, Paris, pp. 643-646, 1990.
}
\examples{
data("Hepta")
e = esomTrain(Hepta$Data, Key = 1:nrow(Hepta$Data))
plotMatrix(e$Umatrix,e$BestMatches)
}
