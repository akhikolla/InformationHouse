#' 3d plotting of an icosahedral grid, its subset or a data layer
#' 
#'  
#' @usage plot3d(x,...)
#' @rdname plot3d
#' @export plot3d
"plot3d"

if(requireNamespace("rgl", quietly = TRUE)){
	# manually copy S3 generic
	plot3d <- rgl::plot3d
}else{
	# manually create S3 generic
	plot3d <-function(x,...){
		if(!requireNamespace("rgl", quietly = TRUE)) stop("Install the 'rgl' package and reload 'icosa' to use this function.")
		UseMethod("plot3d")
	}
}


#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}) Object to be plotted.
#' 
#' @param type (\code{character}) Value specifying the part of the grid to be plotted by the call of the function. 
#' \code{"v"} plots the grid vertex points. 
#' \code{"e"} draws the grid edges.
#' \code{"f"} draws the grid faces.
#' \code{"c"} draws the face centers of the grid.
#' 
#' @param sphere (\code{numeric}) Defaults to \code{NULL}, adding a central white sphere to the plot. Assigning a \code{numeric} value will draw a new sphere with the given radius,
#'		\code{FALSE} does not plot the sphere. 
#' @param guides (\code{logical}) Value indicating whether the guidelines of the polar coordinate system shall be plotted.
#' @param add (\code{logical}) Value indicating whether a new plot shall be drawn, or the currently plotted information should be added to the active \code{rgl} device.
#' 
#' @param ... Further graphical parameters passed to (see \code{\link[rgl]{plot3d}}).
#' 
#' @return The function does not return any value.
#'
#'
#' @examples
#' # create a hexagonal grid
#'     g <- hexagrid(c(2,2))
#' # plot the grid in 3d space
#'     plot3d(g, col="blue")
#' # make a subset to select faces
#'    subG <- subset(g, c("F5", "F2"))
#' # plot the subset defined above
#'     plot3d(subG, type="f", col=c("orange"), add=TRUE, lwd=1)
#' @rdname plot3d
#' @exportS3Method rgl::plot3d trigrid
#' @exportS3Method plot3d trigrid
#' @export plot3d.trigrid
plot3d.trigrid <- function(x, type=c("l"),sphere=NULL,  add=FALSE, guides=TRUE, ...){
	
	#create new plot?
	if(add==FALSE)
	{
		#checking plotting
		rgl::plot3d(x@vertices, type="n", box=FALSE, axes=FALSE, xlab="", ylab="", zlab="")
		
		#default sphere plotting
		if(is.null(sphere)){
			#get the radius
			fc<-apply(x@vertices[x@faces[1,],],2,mean)-x@center
			rad<-sqrt(fc[1]^2+fc[2]^2+fc[3]^2)-15
			blankSphere(x@center[1],x@center[2], x@center[3], radius = rad, color ="white", ng=200, box=FALSE, axes=FALSE)
		}else{
			if(sphere){
				blankSphere(x@center[1],x@center[2], x@center[3], radius = sphere, color ="white", ng=200, box=FALSE, axes=FALSE)
			}
		}
	}
	
	
	if(type=="p")
	{
		#single point
		if(length(x@vertices)==3)
		{
			rgl::points3d(x=x@vertices[1],y=x@vertices[2],z=x@vertices[3],...)
		}else{
			rgl::points3d(x@vertices, ...)
		}
	}
		
	if(type=="l")
	{
		lines3d(x, ...)
	}
	
	if(type=="f")
	{
		faces3d(x,...)
	}
	
	if(type=="n"){
	
	}
	# guides
	if(guides){
		guides3d(col="green", origin=x@center, radius=x@r, lwd=2)
	}
	
	
}

#' 3d plotting of an icosahedral grid or its subset
#' @rdname plot3d
#' @param color (\code{character}) Only for the hexagrid plotting: value/values passed to the \code{\link{faces3d}} function instead of \code{col}.
#' @exportS3Method rgl::plot3d hexagrid
#' @exportS3Method plot3d hexagrid
#' @export plot3d.hexagrid
plot3d.hexagrid <- function(x, type=c("l"),sphere=NULL, color="gray70", add=FALSE, guides=TRUE, ...){
	#create new plot?
	if(add==FALSE)
	{
		#empty plotting plotting
		rgl::plot3d(x@vertices, type="n", box=FALSE, axes=FALSE, xlab="", ylab="", zlab="")
		
		#default sphere plotting
		if(is.null(sphere)){
			fVect<-x@faces[2,!is.na(x@faces[2,])]
			#get the radius
			fc<-apply(x@vertices[fVect,],2,mean)-x@center
			rad<-sqrt(fc[1]^2+fc[2]^2+fc[3]^2)-10
			blankSphere(x@center[1],x@center[2], x@center[3], radius = rad, color ="white", ng=200, box=FALSE, axes=FALSE)
		}else{
			if(sphere){
				blankSphere(x@center[1],x@center[2], x@center[3], radius = sphere, color ="white", ng=200, box=FALSE, axes=FALSE)
			}
		}
	}
	
	
	if(type=="p")
	{
		#single point
		if(length(x@vertices)==3)
		{
			rgl::points3d(x=x@vertices[1],y=x@vertices[2],z=x@vertices[3], col=color,...)
		}else{
			rgl::points3d(x@vertices, col=color, ...)
		}
	}
		
	if(type=="l")
	{
		lines3d(x, ...)
	}
	
	if(type=="f")
	{
		faces3d(x,...)
	}
	
	if(type=="c")
	{
		#single point
		if(length(x@faceCenters)==3)
		{
			rgl::points3d(x=x@faceCenters[1],y=x@faceCenters[2],z=x@faceCenters[3], col=color, ...)
		}else{
			rgl::points3d(x@faceCenters, col=color, ...)
		}
	}
	
	if(type=="t")
	{
		rgl::text3d(x@faceCenters, texts=rownames(x@faceCenters),col=color, ...)
	}
	
	if(type=="n"){
	
	}
	if(guides){
		guides3d(col="green", origin=x@center, radius=x@r, lwd=2)
	}
	
	
}

#' Methods of 3d line plotting
#' 
#' This is a generic function used to plot the edge lines of either a \code{trigrid} or a \code{hexagrid} object, a \code{facelayer}, or \code{Spatial} objects in 3d space. The method is also implemented for 
#' the object classes defined by the package 'sp'.
#' 
#' The function is built on the openGL renderer of the R package \code{rgl}, which needs to be installed for the function to run. Although the function is works without attaching rgl, note that if you want to attach both \code{icosa} and \code{rgl},the \code{rgl} package has to be loaded ifrst otherwise the function will not be usable.
#'  
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}}, \code{\link{facelayer}} or \code{sp}) Object to be plotted.
#' @param arcs \code{logical} Value setting whether great circle arcs or segments shall be drawn betwenn the points of the grid.
#' 
#' @param ... Further graphical parameters passed to (see \code{\link[rgl]{plot3d}}).
#' 
#' @return The function does not return any value.
#'
#' @examples
#' # create a hexagonal grid
#'   g <- hexagrid(c(2,2))
#' # plot the grid in 3d space
#'   lines3d(g, col="blue")
#' @rdname lines3d
"lines3d"
if(requireNamespace("rgl", quietly = TRUE)){
	setGeneric("lines3d", def=rgl::lines3d, package="rgl")
}else{
	setGeneric(
		name="lines3d",
		def=function(x,y=NULL,z=NULL, ...){
			if(!requireNamespace("rgl", quietly = TRUE)) stop("Install the 'rgl' package and reload 'icosa' to use this function.")
			standardGeneric("lines3d")
		}
	)
}


#' lines3d method for the trigrid class
#' @rdname lines3d
#' @exportMethod lines3d
setMethod(
	"lines3d",
	signature="trigrid",
	definition=function(x, arcs=FALSE, ...){
		v<-x@skeleton$v
		e<-x@skeleton$e[x@skeleton$aE,]
		
		#create edgeMat with a simple Rccp function
		edgeMat<-.Call(Cpp_icosa_edgeMatTri_, v=v, e=e)
		
		# get the list of additional arguments
		newArgs<-list(...)
		
		
		if(prod(x@tessellation)<16 & arcs){
			res<-10
			edgeMat<-.Call(Cpp_icosa_expandEdges_, edgeMat, x@center, res)
		}
		edgeMatExp<-edgeMat
		rgl::segments3d(x=edgeMatExp[,1],y=edgeMatExp[,2],z=edgeMatExp[,3], ...)
	}
)


#' Methods of 3D face plotting.
#' 
#' This function is used to plot the faces of either a \code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}} object in 3D space. 
#' 
#' The function is built on the openGL renderer of the R package \code{rgl}.
#'  
#' @name faces3d
#' @param x The \code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}} object to be plotted.
#' 
#' @param ... Further graphical parameters passed to (see \code{\link[rgl]{plot3d}}) and the \code{\link{heatMapLegend}} function.
#' 
#' @return The function does not return any value.
#'
#' @examples
#' # create a hexagonal grid
#'     g <- hexagrid(c(2,2))
#' # plot the grid in 3d space
#'     faces3d(g)
#' @exportMethod faces3d
#' @rdname faces3d
setGeneric(
	name="faces3d",
	package="icosa",
	def=function(x,...){
		if(!requireNamespace("rgl", quietly = TRUE)) stop("Install the 'rgl' package and reload 'icosa' to use this function.")
		standardGeneric("faces3d")
		
	}
)

#' Methods of 3d face plotting.
#' 
#'  
#' @rdname faces3d
setMethod(
	"faces3d",
	signature="trigrid",
	definition=function(x, ...){
		
		v<-x@skeleton$v
		f<-x@skeleton$f[as.logical(x@skeleton$aF),1:3]
		
		#create edgeMat with a simple Rccp function
		triMat<- .Call(Cpp_icosa_triMatTri_, v, f)
		
		rgl::triangles3d(x=triMat[,1],y=triMat[,2],z=triMat[,3],...)
			
	}
)

#' Methods of 3d face plotting.
#' 
#' @rdname faces3d
setMethod(
	"faces3d",
	signature="hexagrid",
	definition=function(x,...){
		v<-x@skeleton$plotV
		f<-x@skeleton$f[as.logical(x@skeleton$aSF),1:3]
		
	#	f2<-f[order(f[,1]),]
	#	
	#	f2<-unique(f2)
		
		#create edgeMat with a simple Rccp function
		triMat<- .Call(Cpp_icosa_triMatTri_, v, f)
		
		rgl::triangles3d(x=triMat[,1],y=triMat[,2],z=triMat[,3],...)
	
	
	}
)


#' Display the names of the grid elements in 3d plots.
#' 
#' This function will display the names of vertices, faces and edges on 3d plots.
#' 
#' @name gridlabs3d
#'  
#' @param gridObj (\code{\link{trigrid}}, \code{\link{hexagrid}}) An icosahedral grid.
#' 
#' @param type (\code{character}) Vector containing either \code{"f"}, \code{"e"} or \code{"v"}, rendering the names
#' of either the faces, edges or vertives respectively.
#'
#' @param ... Additional arguments passed to \code{\link[rgl:texts]{text3d}} function of the \code{rgl} package.
#' 
#' @return The function does not return any value.
#'
#' @examples
#' # create a hexagonal grid
#' g <- hexagrid(c(2,2))
#' # plot the grid in 3d space
#' lines3d(g, guides=FALSE)
#' # labels
#' gridlabs3d(g)
#' @exportMethod gridlabs3d
#' @rdname gridlabs3d
setGeneric(
	name="gridlabs3d",
	package="icosagrid",
	def=function(gridObj,...){
		if(!requireNamespace("rgl", quietly = TRUE)) stop("Install the 'rgl' package and reload 'icosa' to use this function.")
		standardGeneric("gridlabs3d")
		
	}
)

#' @rdname gridlabs3d
setMethod(
	"gridlabs3d",
	signature="trigrid",
	definition=function(gridObj,type="f",...){
		
		# vertex names
		if("v"%in%type){
			rgl::text3d(texts=rownames(gridObj@vertices), gridObj@vertices*1.005,...)
		}
		
		
		if("f"%in%type){
			rgl::text3d(texts=rownames(gridObj@faceCenters), gridObj@faceCenters*1.005,...)
		}
		
		if("e"%in%type){
			#the coordinates
			coords<-apply(gridObj@edges,1,function(x){
				apply(gridObj@vertices[x,], 2, mean)
			})
			
			rgl::text3d(t(coords)*1.005, texts=rownames(gridObj@edges),...)
		}
	
	
	}
)
	
#' @rdname gridlabs3d
setMethod(
	"gridlabs3d",
	signature="hexagrid",
	definition=function(gridObj,type="f",...){
		
		# vertex names
		if("v"%in%type){
			rgl::text3d(texts=rownames(gridObj@vertices), gridObj@vertices*1.005,...)
		}
		
		
		if("f"%in%type){
			rgl::text3d(texts=rownames(gridObj@faceCenters), gridObj@faceCenters*1.005,...)
		}
		
		if("e"%in%type){
			#the coordinates
			coords<-apply(gridObj@edges,1,function(x){
				apply(gridObj@vertices[x,], 2, mean)
			})
			
			rgl::text3d(t(coords)*1.005, texts=rownames(gridObj@edges),...)
		}
	
	
	}
)

