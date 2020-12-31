setMethod(	
	f="show",
	signature="obj3d",
	definition= function(object){
		cat(paste("A/An ", class(object), " object with ", 
			object@length[1], " vertices, ",
			object@length[2], " edges and ",
			object@length[3], " faces.\n",
			sep=""))
		if(class(object)%in% c("trigrid", "hexagrid")){
			cat(paste("The mean grid edge length is ", 
			round(object@edgeLength[1],2), " km or ",
			round(object@edgeLength[2],2), " degrees.\n", sep=""))
		}
		cat("Use plot3d() to see a 3d render.\n")
	}
)


#' The number of faces in a \code{trigrid} or \code{hexagrid} class object. 
#'
#' The length of the object is interpreted as the number of faces it contains.
#'
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}) The object.
#' @return An integer value.
#' @rdname length
#' @exportMethod length
setMethod(	
	f="length",
	signature="trigrid",
	definition= function(x){
		return(x@length[3])
	}
)


#' The vertices of an icosahedral grid object
#'
#' Shorthand function to return the vertices slot of an icosahedral grid or a grid linked to a facelayer. 
#' @name vertices
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}) The icosahedral grid, or linked data object.
#' @param ... Additional arguments passed to class-specific methods.
#' @param output (\code{character}) The coordinate system of output.
#' @rdname vertices
#' @exportMethod vertices
#' @examples
#' a <- trigrid(1)
#' vertices(a)
setGeneric(
	name="vertices",
	def=function(x,...){
		standardGeneric("vertices")
	}
)

#' @rdname vertices
setMethod(	
	f="vertices",
	signature=c("trigrid"),
	definition= function(x, output="polar"){
		if(output=="polar"){
			return(CarToPol(x@vertices, origin=x@center, norad=TRUE))
		
		}else{
			return(x@vertices)
		}
	}
)

#' @rdname vertices
setMethod(	
	f="vertices",
	signature=c("facelayer"),
	definition= function(x, output="polar"){
		actGrid <- get(x@grid)
		if(output=="polar"){
			return(CarToPol(actGrid@vertices, origin=actGrid@center, norad=TRUE))
		
		}else{
			return(actGrid@vertices)
		}
	}
)

#' The faces of a 3d object
#'
#' Shorthand function to get the faces slot of an icosahedral grid or a grid linked to a \code{\link{facelayer}}. 
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}) The grid or facelayer object.
#' @name faces
#' @rdname faces
#' @return The faces of the grid as a \code{character} matrix.
#' @exportMethod faces
setGeneric(
	name="faces",
	def=function(x){
		standardGeneric("faces")
	}
)

#' @rdname faces
setMethod(	
	f="faces",
	signature="trigrid",
	definition= function(x){
		return(x@faces)
	}
)

#' @rdname faces
setMethod(	
	f="faces",
	signature="gridlayer",
	definition= function(x){
		actGrid<-get(x@grid)
		return(actGrid@faces)
	}
)

#	#' @rdname faces
#	#' @exportMethod faces
#	setMethod(	
#		f="faces",
#		signature="facelayer",
#		definition= function(x){
#			actGrid<-get(x@grid)
#			return(actGrid@faces)
#		}
#	)


#' The edges of a 3d object
#'
#' Shorthand function to get the edges slot of an icosahedral grid or a grid linked to a facelayer. 
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}) The grid or linked data object.

#' @name edges
#' @rdname edges
#' @exportMethod edges
setGeneric(
	name="edges",
	def=function(x){
		standardGeneric("edges")
	}
)

#' @rdname edges
setMethod(	
	f="edges",
	signature="obj3d",
	definition= function(x){
		return(x@edges)
	}
)

#' @exportMethod edges
#' @return The edges of the grid, as a \code{character} matrix. 
#' @rdname edges
setMethod(	
	f="edges",
	signature="facelayer",
	definition= function(x){
		actGrid<-get(x@grid)
		return(actGrid@edges)
	}
)

#' The face centers of an icosahedral grid object
#'
#' Shorthand function to return the \code{@faceCenters} slot of an icosahedral grid or a grid linked to a facelayer. 
#' @name centers
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}). The grid or linked data layer object.
#' @param ... Arguments passed to the class specific methods.
#' @rdname centers
#' @return The coordinates of the face centers as a \code{numeric} matrix.
#' @examples
#' a <- trigrid()
#' centers(a)
#' @exportMethod centers
setGeneric(
	name="centers",
	def=function(x,...){
		standardGeneric("centers")
	}
)
#' The face centers of a trigrid or hexagrid class object
#'
#' @param output (\code{character}) The coordinate system of the output. Either \code{"polar"} or \code{"cartesian"}.
#' @rdname centers
setMethod(	
	f="centers",
	signature="trigrid",
	definition= function(x, output="polar"){
		if(output=="polar"){
			return(CarToPol(x@faceCenters, origin=x@center, norad=TRUE))
		
		}else{
			return(x@faceCenters)
		}
	}
)

#' The face centers of a trigrid or hexagrid class object that is linked to a facelayer
#'
#' @rdname centers
#' @exportMethod centers
setMethod(	
	f="centers",
	signature="facelayer",
	definition= function(x, output="polar"){
		actGrid <- get(x@grid)
		if(output=="polar"){
			return(CarToPol(actGrid@faceCenters, origin=actGrid@center, norad=TRUE))
		
		}else{
			return(actGrid@faceCenters)
		}
	}
)

	
#' Extracting and setting the grid orientation
#' 
#' @name orientation
#' 
#' @param x (\code{\link{trigrid}} or \code{\link{hexagrid}}): Input grid. 
#' @param display (\code{character}) The output unit. In case it is set to \code{"deg"} the output will be in degrees, in case it is \code{"rad"}, then radians.
#' @exportMethod orientation
#' 
#' @rdname orientation
setGeneric(
	name="orientation",
	package="icosa",
	def=function(x,...){
		standardGeneric("orientation")
	}
)

#' @rdname orientation
setMethod(
	"orientation",
	signature="trigrid",
	definition=function(x,display="deg",...){
		if(display=="rad"){
			names(x@orientation)<-c("x (rad)", "y (rad)", "z (rad)")
			return(x@orientation)
		}
		if(display=="deg"){
			names(x@orientation)<-c("x (deg)", "y (deg)", "z (deg)")
			return(x@orientation/pi*180)
			
		}
	
	}
	
)




#' @name orientation<-
#' 
#' @param value (\code{numeric}) The vector of rotation. Passed as the \code{angles} argument of \code{\link[icosa]{rotate}}.
#' @param ... Values passed on to the \code{\link[icosa]{rotate}} function.
#' 
#' @return In case the function returns does, it returns the orientation angles of the grid (as \code{numeric}).
#' @exportMethod orientation<-
#' 
#' @rdname orientation
setGeneric(
	name="orientation<-",
	def=function(x,value){
		standardGeneric("orientation<-")
	}
	
)


#' @rdname orientation
setReplaceMethod(
	"orientation",
	signature="trigrid",
	definition=function(x, value){
		x<-rotate(x, angles=value)
		return(x)
	
	}
)




#' Lengths of grid edges
#' 
#' This function will return the length of all edges in the specified grid object.
#' 
#' @name edgelength
#' @param gridObj (\code{\link{trigrid}} or \code{{hexagrid}}) A grid object. 
#' 
#' @param ... Arguments passed to the class specific methods.
#' 
#' @examples
#' g <- trigrid(3)
#' edges <- edgelength(g, output="deg")
#' edges
#' 
#' @return A named \code{numeric} vector.
#' 	
#' @exportMethod edgelength
#' @rdname edgelength
setGeneric(
	name="edgelength",
	def=function(gridObj,...){
		standardGeneric("edgelength")
	}

)

#' Lengths of grid edges
#' 
#' 
#' @param output (\code{character}) The type of the output. \code{"distance"} will give back the distance
#'	in the metric that was fed to the function in the coordinates or the radius.
#'	\code{"deg"} will output the the distance in degrees, \code{"rad"} will do
#'	so in radians.
#' 
#' @rdname edgelength
setMethod(
	"edgelength", 
	signature="trigrid",
	definition=function(gridObj, output="distance"){
		if(!output%in%c("distance", "rad", "deg")) stop("Invalid distance method.")
		v<-gridObj@skeleton$v
		e<-gridObj@skeleton$e
		
		if(output=="distance"){
			met<-1
		}else{
			met<-0
		}
		sizes<-  .Call(Cpp_icosa_edges_, v, e, origin=as.integer(gridObj@center), method=as.logical(met))
		
		names(sizes)<-rownames(gridObj@edges)
		
		if(output=="deg") sizes<-sizes/pi*180
		
		return(sizes)
		
	}
)


#' Areas of grid cell surfaces
#' 
#' This function will return the areas of all cells in the specified grid object.
#' 
#' @name surfacearea
#' @param gridObj (\code{\link{trigrid}} or \code{\link{hexagrid}}) Object. 
#' 
#' 
#' @examples
#' g <- trigrid(3)
#' surfaces <- surfacearea(g)
#' surfaces
#' 
#' @return A named \code{numeric} vector, in the metric that was given to the function in the coordinates or the radius. \code{"deg"} will output the the distance in degrees, \code{"rad"} will do so in radians.
#' 	
#' @rdname surfacearea
#' @exportMethod surfacearea
setGeneric(
	name="surfacearea",
	package="icosa",
	def=function(gridObj){
		standardGeneric("surfacearea")
		
	}
)

#' @rdname surfacearea
setMethod(
	"surfacearea", 
	signature="trigrid", 
	def=function(gridObj){
		# get the highest resolution faces
		newF <- gridObj@skeleton$f[as.logical(gridObj@skeleton$aF),1:3]
		v <- gridObj@skeleton$v
		
		# call the surface calculation function
		surfInner <-  .Call(Cpp_icosa_spherTriSurfs,
			v=v, 
			f=newF, 
			origin=gridObj@center, 
			pi=pi
		)
		
		# reorganize the faces: outer representation
		ord<-gridObj@skeleton$aF[as.logical(gridObj@skeleton$aF)]
		
		surfOuter<-surfInner
		surfOuter[ord]<- surfInner
		
		names(surfOuter) <- rownames(gridObj@faces)
		
		return(surfOuter)
	}
)

#' @rdname surfacearea
setMethod(
	"surfacearea", 
	signature="hexagrid", 
	def=function(gridObj){
		# get the highest resolution faces
		newF <- gridObj@skeleton$f[as.logical(gridObj@skeleton$aSF),1:3]
		v <- gridObj@skeleton$v
		
		# call the surface calculation function
		surfInner <-  .Call(Cpp_icosa_spherTriSurfs, 
			v=v, 
			f=newF, 
			origin=gridObj@center, 
			pi=pi
		)
		
		# the subfaces belong to these face IDs in the outer representation
		aS<-gridObj@skeleton$aSF[as.logical(gridObj@skeleton$aSF)]
		
		# calculate the sums of all subface areas in a face, and order them
		doubleSurf<-tapply(INDEX=aS, X=surfInner, sum)
		
		# each subface occurs two times in the f matrix, divide area by 2
		singleSurf <- doubleSurf/2
		
		# augment the names attributes
		names(singleSurf)<- paste("F", names(singleSurf), sep="")
		
		return(singleSurf)
	}
)


#' Shape distortions of the triangular faces and subfaces
#' 
#' This function will return a value that is proportional to the irregularity of a triangonal face or subface. The ratio of the lengths of the shortest and the longest edges.
#' 
#' The value is exactly \code{1} for an equilateral triangle, and becomes \code{0} as one of the edges approach \code{0}.
#'
#' @name trishape
#' @param gridObj (\code{\link{trigrid}}, \code{\link{hexagrid}}) Object. 
#' 
#' @examples
#' g <- trigrid(3)
#' shape <- trishape(g)
#' 
#' 
#' @return A named \code{numeric} vector, one value for every face of the grid.
#' 	
#' @rdname trishape
#' @exportMethod trishape
setGeneric(
	name="trishape",
	def=function(gridObj){
		standardGeneric("trishape")
	}
)

#' @rdname trishape
setMethod(
	"trishape",
	signature="trigrid",
	definition=function(gridObj){
		# center back to origin if not there already
		if(gridObj@center[1]!=0 | gridObj@center[2]!=0 | gridObj@center[3]!=0){
			gridObj<-translate(gridObj,-gridObj@center)
		}
		
		v <- gridObj@skeleton$v
		f <- gridObj@skeleton$f[as.logical(gridObj@skeleton$aF),]
		
		# the shapes in the inner order
		innerShape <-.Call(Cpp_icosa_AllShapeTri_, v, f)
		outerShape <- innerShape[gridObj@skeleton$uiF]

		return(outerShape)
	}
)

#' @rdname trishape
setMethod(
	"trishape",
	signature="hexagrid",
	definition=function(gridObj){
		# center back to origin if not there already
		if(gridObj@center[1]!=0 | gridObj@center[2]!=0 | gridObj@center[3]!=0){
			gridObj<-translate(gridObj,-gridObj@center)
		}
		
		v <- gridObj@skeleton$v
		f <- gridObj@skeleton$f[as.logical(gridObj@skeleton$aSF),]
		
		# the shapes in the inner order
		innerShape <-.Call(Cpp_icosa_AllShapeTri_, v, f)
		outerShape<- tapply(
			X=innerShape, 
			INDEX=gridObj@skeleton$aSF[as.logical(gridObj@skeleton$aSF)],
			FUN=mean)	
		outerShape<-as.numeric(outerShape)
		names(outerShape)<-paste("F", names(outerShape), sep="")
		
		return(outerShape)
		
	
	}
)

		