## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(rgdal)
library(raster)
library(knitr)
library(rgl)
library(icosa)
knit_hooks$set(rgl = hook_rgl)

## ----packages-----------------------------------------------------------------
library(rgdal)
library(raster)
library(rgl)
library(icosa)

## ----first, rgl=TRUE,dev='png',dpi=300, fig.width=6, fig.height=4-------------
# create a trigrid class object
tri <- trigrid()

# the show() method displays basic information of the package
tri

# plot the object in 3d
plot3d(tri, guides=F)

## ----crs, rgl=TRUE,dev='png',dpi=300------------------------------------------
tri@proj4string

## ----trigrid, rgl=TRUE,dev='png',dpi=300--------------------------------------
# create a trigrid class object
gLow <- trigrid(tessellation=c(4,4))

# plot the object in 3d
plot3d(gLow, guides=F)

## ----hexagrid1, rgl=TRUE,dev='png',dpi=300------------------------------------
# create a hexagrid object
hLow <- hexagrid()

# plot it in 3d
plot3d(hLow, guides=F)

## ----hexagrid, rgl=TRUE,dev='png',dpi=300-------------------------------------
# create a hexagrid object
hLow <- hexagrid(c(4,4))

# plot it in 3d
plot3d(hLow)

## ----tri2, rgl=TRUE,dev='png',dpi=300-----------------------------------------
# the beginning of the vertices matrix
head(gLow@vertices)

# the beginning of the faces matrix
head(gLow@faces)

## ----tri3, rgl=TRUE,dev='png',dpi=300-----------------------------------------
# grid radius
gLow@r

# grid center
gLow@center

## ----tri5, rgl=TRUE,dev='png',dpi=300-----------------------------------------
head(gLow@faceCenters)

## ----tri7, rgl=TRUE,dev='png',dpi=300-----------------------------------------
head(gLow@edges)

## ----tri6, rgl=TRUE,dev='png',dpi=300-----------------------------------------
plot3d(gLow)
gridlabs3d(gLow, type="v", col="blue", cex=0.6)

## ----tri8, rgl=TRUE,dev='png',dpi=300-----------------------------------------
gLow2 <- rotate(gLow) # random rotation
plot3d(gLow2)
guides3d(col="green")

## ----first2, rgl=TRUE,dev='png',dpi=300, fig.width=6, fig.height=4------------
plot3d(tri, guides=F, arcs=T, sphere=6300)

## ----tri10, plot=TRUE, ,dev='png',dpi=300-------------------------------------
hLow <- hexagrid(c(4,4), sp=TRUE)
# After this procedure finishes, a regular 2d plotting function can be invoked:
plot(hLow)

## ----tri23, plot=TRUE,dev='png',dpi=300,echo=TRUE, results='hide'-------------
# use the rgdal package
library(rgdal)

# file path
file <- system.file("extdata", "land_polygons_z1.shx", package = "icosa")

# read in the shape file
wo <- readOGR(file, "land_polygons_z1")

## ----examplePaleo, plot=TRUE, ,dev='png',dpi=300------------------------------
# transform the land data to long-lat coordinates
wo <- spTransform(wo, gLow@proj4string)

#triangular grid
gLow<-newsp(gLow)

# # load in a map
# plot the grid (default longitude/latitude)
plot(gLow, border="gray", lty=1)
 
# the reconstruction
lines(wo, lwd=2, col="blue")

## ----tri12, plot=TRUE, ,dev='png',dpi=300-------------------------------------
# a very low resolution hexagrid
hVeryLow<-hexagrid(c(4))
# add 2d component
hVeryLow<-newsp(hVeryLow)
# the Robinson projection
robin <- CRS("+proj=robin")
# plot with labels
plot(hVeryLow, projargs=robin)
gridlabs(hVeryLow, type="f", cex=0.6,projargs=robin)

## ----tri13, echo=TRUE---------------------------------------------------------
pos(hLow, c("P2", "F12", NA))

## ----tri14, echo=TRUE---------------------------------------------------------
fl1<-facelayer(gLow) # the argument is the grid object to which the layer is linked
fl1
str(fl1)

## ----tri15, echo=TRUE---------------------------------------------------------
length(fl1)

## ----tri16, echo=TRUE---------------------------------------------------------
values(fl1) <-1:length(fl1)
values(fl1)[1:10]

## ----tri18, rgl=TRUE,dev='png',dpi=300----------------------------------------
a <-facelayer(gLow)
values(a) <- sample(c(T,F), length(a), replace=T)
# plot the grid first
plot3d(gLow, guides=F)

# invoke lower level plotting for the facelayer 
# (draws on previously plotted rgl environemnts)
faces3d(a, col="green")

## ----tri19b, rgl=TRUE,dev='png',dpi=300---------------------------------------
# new layer
# grid frame
plot3d(gLow)
# the heatmap
faces3d(a, col=c("green", "brown")) 

## ----categ1, rgl=TRUE,dev='png',dpi=300---------------------------------------
# new layer
catLayer<-facelayer(hLow)

# assign random information
catLayer@values<-sample(c("one","two","three"),length(catLayer), replace=T)

plot(catLayer)

## ----tri20, rgl=TRUE,dev='png',dpi=300----------------------------------------
# generate 5000 random coordinates on a sphere of default radius
pointdat <- rpsphere(5000)

# and locate them on the grid 'gLow'
cells<-locate(gLow, pointdat)

# the return of this function is vector of cell names
head(cells)

## ----tri20b, rgl=TRUE,dev='png',dpi=300---------------------------------------
tCell <- table(cells)
fl <- facelayer(gLow,0)
# [] invokes a method that save the values to places that 
# correspond to the names attribute of tCell
fl[] <-tCell #
# heat map of the point densities
plot3d(fl)
 

## ----tri22, rgl=TRUE,dev='png',dpi=300----------------------------------------
# run function only on the first 300
fl<-occupied(hLow, pointdat[1:300,])

# the plot function can also be applied to the facelayer object
plot(fl, col="blue")

# show the points as well
points(CarToPol(pointdat[1:300,]), col="red", pch=3, cex=0.7)

## ----tri22b, rgl=TRUE,dev='png',dpi=300---------------------------------------
# the plot function can also be applied to the facelayer object
plot(fl, col="blue")

points(CarToPol(pointdat[1:300,]), col="red", pch=3, cex=0.7)
lines(hLow, col="gray")

## ----vic1, rgl=TRUE,dev='png',dpi=300-----------------------------------------
# calculate a very coarse resolution grid
gVeryLow<-trigrid(8, sp=T)
# names of faces that are neighbours to face F125
facenames<-vicinity(gVeryLow, "F125")
# plot a portion of the grid
plot(gVeryLow, xlim=c(0,180), ylim=c(0,90))
# plot the original and the neighbouring faces
plot(gVeryLow@sp[facenames], col="red", add=T)
# the names of all the cells
gridlabs(gVeryLow, type="f", cex=0.5)

## ----ggraphAttach, results='hide'---------------------------------------------
# attach igraph
library(igraph)

## ----ggraph4,dev='png',dpi=300, echo=TRUE, results='hide'---------------------
faces<-paste("F", 1:10, sep="")
subGraph <- induced_subgraph(gVeryLow@graph,faces) 
plot(subGraph)

## ----ggraph5,dev='png',dpi=300------------------------------------------------
lowGraph<-gLow[1:12]@graph

## ----ggraph4b, dev='png',dpi=300----------------------------------------------
# look up the polygons
landFaces<-occupied(hLow, wo)

# create a new grid from a facelayer
landGraph<-gridgraph(landFaces)
plot(landFaces, col="brown")

## ----ggraph6, dev='png',dpi=300-----------------------------------------------
# shortest path in igraph
path <- shortest_paths(landGraph, from="F432", to="F1073", output="vpath")
# the names of the cells in order
cells<-path$vpath[[1]]$name
# plot the map
plot(landFaces, col="brown", xlim=c(0,90), ylim=c(0,90))
# make a subset of the grid - which corresponds to the path
routeGrid<-hLow[cells]
# plot the path
plot(routeGrid, col="red", add=T)

## ----ggraph6b, dev='png',dpi=300----------------------------------------------
# plot the map
plot(landFaces, col="brown", xlim=c(0,90), ylim=c(0,90))
# create a random walk from source cell with a given no. of steps
randomWalk <- random_walk(landGraph, steps=100, start="F432")
# the names of the cells visited by the random walker
cells<-randomWalk$name
# the source cell
plot(hLow["F432"], col="green",add=T)
# the centers of these faces
centers<-CarToPol(hLow@faceCenters[cells,], norad=T)
# draw the lines of the random walk
for(i in 2:nrow(centers)){
	segments(x0=centers[i-1,1], y0=centers[i-1,2], x1=centers[i,1], y1=centers[i,2], lwd=2)
}

