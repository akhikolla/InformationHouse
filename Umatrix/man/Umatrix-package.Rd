\name{Umatrix-package}
\alias{Umatrix-package}
\docType{package}
\title{
Umatrix-package
}
\description{
The \bold{ESOM}(emergent self organizing map) is an improvement of the regular \bold{SOM}(self organizing map) which allows for toroid grids of neurons and is intended to be used in combination with the \bold{Umatrix}. The set of neurons is referred to as \bold{weights} within this package, as they represent the values within the high dimensional space. The neuron with smallest distance to a datapoint is called a \bold{Bestmatch} and can be considered as projection of said datapoint. As the Umatrix is usually toroid, it is drawn four consecutive times to remove border effects. An island, or Imx, is a filter mask, which cuts out a subset of the Umatrix, which shows every point only a single time while avoiding border effects cutting through potential clusters.
Finally the Pmatrix shows the density structures within the grid, by a set radius. It can be combined with the Umatrix resulting in the UStarMatrix, which is therefore a combination of density based structures as well as clearly divided ones.
}

\references{
Ultsch, A.: Data mining and knowledge discovery with emergent self-organizing feature maps for multivariate time series, In Oja, E. & Kaski, S. (Eds.), Kohonen maps, (1 ed., pp. 33-46), Elsevier, 1999.

Ultsch, A.: Maps for the visualization of high-dimensional data spaces, Proc. Workshop on Self organizing Maps (WSOM), pp. 225-230, Kyushu, Japan, 2003.

Ultsch, A.: U* C: Self-organized Clustering with Emergent Feature Maps, Lernen, Wissensentdeckung und Adaptivitaet (LWA), pp. 240-244, Saarbruecken, Germany, 2005.

Loetsch, J., Ultsch, A.: Exploiting the Structures of the U-Matrix, in Villmann, T., Schleif, F.-M., Kaden, M. & Lange, M. (eds.), Proc. Advances in Self-Organizing Maps and Learning Vector Quantization, pp. 249-257, Springer International Publishing, Mittweida, Germany, 2014.

Ultsch, A., Behnisch, M., Loetsch, J.: ESOM Visualizations for Quality Assessment in Clustering, In Merenyi, E., Mendenhall, J. M. & O'Driscoll, P. (Eds.), Advances in Self-Organizing Maps and Learning Vector Quantization: Proceedings of the 11th International Workshop WSOM 2016, pp. 39-48, Houston, Texas, USA, January 6-8, 2016, (10.1007/978-3-319-28518-4_3), Cham, Springer International Publishing, 2016.

Thrun, M. C., Lerch, F., Loetsch, J., Ultsch, A.: Visualization and 3D Printing of Multivariate Data of Biomarkers, in Skala, V. (Ed.), International Conference in Central Europe on Computer Graphics, Visualization and Computer Vision,Plzen, 2016.
}

\keyword{Umatrix}
\keyword{U-Matrix}
\keyword{ESOM}
\keyword{SOM}
\keyword{Projection}
\keyword{Visualization}
\keyword{Highdimensional}
\keyword{High-dimensional}
\keyword{High-Dimensional Data}
