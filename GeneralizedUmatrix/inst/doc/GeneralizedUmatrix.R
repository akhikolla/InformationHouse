## ----setup, include=FALSE-----------------------------------------------------
library(rgl)
setupKnitr()
knitr::opts_chunk$set(echo = TRUE,
                      warning = FALSE,
                      dpi=100,
                      fig.align='center',
                      message = FALSE
                      )
library(GeneralizedUmatrix)

## ----testgl, webgl=TRUE, fig.width=6,fig.height = 6,fig.keep="none"-----------
data(Chainlink)
Data=Chainlink$Data
Cls=Chainlink$Cls
require(DataVisualizations)
DataVisualizations::Plot3D(
  Data,
  Cls,
  type = 's',
  radius = 0.1,
  box = F,
  aspect = T,
  top = T
)
rgl::grid3d(c("x", "y", "z"))

## ---- fig.width=1,fig.height = 1,fig.show='hide'------------------------------
InputDistances = as.matrix(dist(Data))
model = cmdscale(
d = InputDistances,
k = 2,
eig = TRUE,
add = FALSE,
x.ret = FALSE
)
ProjectedPoints = as.matrix(model$points)

## ---- fig.width=6,fig.height = 6,fig.keep='first'-----------------------------
plot(ProjectedPoints, col = Cls)

## ----results=FALSE,fig.show='hide',fig.keep="none"----------------------------
genUmatrix = GeneralizedUmatrix(Data, ProjectedPoints)

## ---- fig.width=7,fig.height = 7,fig.keep='high',fig.show='asis', webgl=TRUE----
plotTopographicMap(genUmatrix$Umatrix, 
                   genUmatrix$Bestmatches, 
                   NoLevels = 10
)

