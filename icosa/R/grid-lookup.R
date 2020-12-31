## wrapper function around the OccupiedFaces generic, to get the occupied grid cells. return out a facelayer object
#export!

#' Faces occupied by the specified object
#'
#' This function will return a \code{\link{facelayer}} class object showing which faces are occupied by the input object.
#'
#' This is a wrapper function on the \code{OccupiedFaces} methods that are specific to grid class and input data. 
#'
#' @param gridObj (\code{\link{trigrid}} or \code{\link{hexagrid}}) An icoshedral grid.
#' 
#' @param data (\code{matrix}, \code{data.frame} or \code{Spatial}) The queried data.
#'
#' @param ... Arguments passed to the class specific methods
#'
#' @return The function Returns a \code{\link{facelayer}}-class object. 
#'
#' @examples
#'	# create a grid
#'	g <- trigrid(8, sp=TRUE)
#'
#'	# create random points
#'	randPoints <- rpsphere(100,output="polar")
#'
#'	# the facelayer occupied by these points
#'	randomLayer <- occupied(g, randPoints)
#'	plot(randomLayer)
#'	points(randPoints, col="blue", pch="+")
#'	
#'
#' @export	
occupied  <- function(gridObj, data,...){
	
	# do spatial transformation if a CRS is present
	if(methods::.hasSlot(data, "proj4string")){
		# and only if it is not NA
		if(!is.na(data@proj4string)){
			# need rgdal
			if(requireNamespace("rgdal", quietly = TRUE)){
				data<-sp::spTransform(data, gridObj@proj4string)
			} else{
				stop("The rgdal package is required to appropriately project this object. ")
			}
		}
	}
	
	#exectue the appropriate searching procedure
	boolVec<-OccupiedFaces(gridObj, data,...)
	
	# construct a facelayer
	endObj<-facelayer(gridObj)
	
	#outer ordering for the hexagrid
	translNum<-which(boolVec)
		
	# replace
	endObj@values<-rep(FALSE, length(endObj))
	endObj@values[translNum]<-TRUE

	# get the name of the grid - 
	endObj@grid<-deparse(substitute(gridObj))
	return(endObj)

}



#OccupiedFaces
#returns a boolean for all the faces (ordered as the internal representation in skeleton$f)
#use this when you do not need to know which face contains which point
setGeneric(
		name="OccupiedFaces",
		def=function(gridObj,data,...){
			standardGeneric("OccupiedFaces")
		}
	
	)
	

# occupied method for the trigrid v6.0
setMethod(
	"OccupiedFaces", 
	signature=c("trigrid", "matrix"),
	definition=function(gridObj, data){
		#locate the cells
		occCells<-locate(gridObj, data, randomborder=FALSE, output="ui")
		
		# the logical vector indicating the face
		boolVec<-rep(FALSE, nrow(gridObj@faces))
		
		boolVec[rownames(gridObj@faces)%in%occCells] <- TRUE
		
		return(boolVec)
	}
)
	
# for spatial points
setMethod(
	"OccupiedFaces",
	signature=c("trigrid", "SpatialPoints"),
	definition=function(gridObj, data){
		borders<-NA
		# basic method for matrices
		OccupiedFaces(gridObj, data@coords)
	}
)
	
# for polygon occupation development
# v2.0 - using igraph
# 2017.02.22.
setMethod(
	"OccupiedFaces", 
	signature=c("trigrid", "Polygon"),
	definition=function(gridObj, data){
		#if no @graph found
		if(suppressWarnings(is.na(gridObj@graph)[1])){
			stop("Slot @graph is empty. Use newgraph() to add an igraph respresentation. ")
		}
		#get the number of faces occupied by the line
			lin<-PolToCar(data@coords, origin=gridObj@center, radius=gridObj@r)
			lin2<- .Call(Cpp_icosa_EvenInterpolation_, lin, gridObj@center, gridObj@edgeLength[2]/180*pi/15)
			
			lineCells<-unique(locate(gridObj,lin2))
		
		# get all the faces
			allFaces<-rownames(gridObj@faces)
			subFaces<-allFaces[!allFaces%in%lineCells]
			
			subGraph<-igraph::induced_subgraph(gridObj@graph, subFaces)
			clusters <- igraph::membership(igraph::clusters(subGraph))
		
		
		#sample the middle part
		middleSample<-sp::spsample(data, type="regular", n=25)
		middleCells<-unique(locate(gridObj, middleSample@coords))
			
		# the group ID of this unit
			clusterIDs<-clusters[names(clusters)%in%middleCells]
			
		# the inner part of the faces
		innerFaces <- names(clusters)[clusters%in%clusterIDs]
		fLayer <- rep(F, nrow(gridObj@faces))
		fLayer[rownames(gridObj@faces)%in%c(lineCells, innerFaces)]<-TRUE
		
		return(fLayer)
	}
)

#for Polygons 
setMethod(
	"OccupiedFaces",
	signature=c("trigrid", "Polygons"),
	definition=function(gridObj, data, n=10000,...){
		borders<-NA
		#faces on the line
		coordLine<-lines3d(data, plot=FALSE)
		coordLine<-coordLine[!is.na(coordLine[,1]),]
		#look these up
		lineFaces<-OccupiedFaces(gridObj, coordLine)
		#get all the sampling points
		all<-sp::spsample(data, type="regular", n=n)

		inFaces<-OccupiedFaces(gridObj, all@coords)
		fl <- inFaces | lineFaces
		return(fl)
	}
)

#for SpatialLines 
setMethod(
	"OccupiedFaces",
	signature=c("trigrid", "SpatialLines"),
	definition=function(gridObj, data, f=5){
		borders<-NA
		# increase resolution
		data<-linIntCoords(data, res=f)
		
		# get the coordinates
		coords<-lines3d(data, plot=TRUE)
		
		#get rid of NAs
		coords<-coords[!is.na(coords[,1]),]
		
		# the faces occupied by the line
		occupiedByLine<-OccupiedFaces(gridObj, coords)
		
		return(occupiedByLine)
	}
)
		
#for SpatialPolygons 
# v. 3.0
setMethod(
	"OccupiedFaces",
	signature=c("trigrid", "SpatialPolygons"),
	definition=function(gridObj, data){
		if(!requireNamespace("raster", quietly = TRUE)) stop("Install the 'raster' package to run this function.")
				
		borders<-NA
		
		#faces on the line
		coordLine<-lines3d(data, plot=FALSE)
		coordLine<-coordLine[!is.na(coordLine[,1]),]
		#look these up
		lineFaces<-OccupiedFaces(gridObj, coordLine)
		
		# create a raster from the SpatialPolygons
		r <-raster::raster()
		
		# set the resolution to that of the grid
		raster::res(r)<-min(edgelength(gridObj, output="deg"))/4
		
		#rasterize it
		data<-raster::rasterize(data,r)
				
		# use the OccupiedFaces method of the raster
		inFaces<-OccupiedFaces(gridObj, data)
		fl <- inFaces | lineFaces
		return(fl)
	}
)

#for SpatialPolygonsDataFrame
setMethod(
	"OccupiedFaces",
	signature=c("trigrid", "SpatialPolygonsDataFrame"),
	definition=function(gridObj, data){
		borders<-NA
		temp<-methods::as(data,"SpatialPolygons")
		fl <- OccupiedFaces(gridObj, temp)
		return(fl)
	}
)


# for spatial points
setMethod(
	"OccupiedFaces",
	signature=c("trigrid", "RasterLayer"),
	definition=function(gridObj, data){
		if(!requireNamespace("raster", quietly = TRUE)) stop("Install the 'raster' package to run this function.")
		
		borders<-NA
		resGrid<-mean(edgelength(gridObj,"deg"))
		# if the default resolution of the raster is too coarse for the trigrid
		if(resGrid<(4*raster::res(data)[1]) | resGrid<(4*raster::res(data)[2])){
			#upscale
			r<-data
			raster::res(r)<-resGrid/4
			data<-raster::resample(data, r, "ngb")
		}
		
		xmin<-data@extent@xmin
		xmax<-data@extent@xmax
		ymin<-data@extent@ymin
		ymax<-data@extent@ymax
		
		xres<-raster::res(data)[1]
		yres<-raster::res(data)[2]
		
		xs<-seq(xmin+xres/2, xmax-xres/2,xres)
		ys<-seq(ymax-yres/2, ymin+yres/2,-yres)
		
		x<-rep(xs, length(ys))
		y<-rep(ys, each=length(xs))
		mat<-cbind(x,y)
		
		cells<-locate(gridObj, mat)
		
		occup<-tapply(X=raster::values(data), INDEX=cells, function(x){sum(!is.na(x))})
		occupiedCells<-names(occup)[occup>0]
		
		fl<-rep(FALSE, length(gridObj))
		fl[rownames(gridObj@faces)%in%occupiedCells]<-T
		
		return(fl)
	}
)



#' Basic lookup function of coordinates on an icosahedral grid
#'
#' @name locate
#' @return The function returns the cell names (as \code{character}) where the input coordinates fall.
#'
#' @param x (\code{trigrid}, \code{hexagrid}) Icosahedral grid object.
#' @param y (\code{matrix}, \code{data.frame}, \code{numeric} or \code{Spatial}) Coordinates of individual points. Can be either a two-dimensional 
#' matrix of long-lat coordinates, a three-dimensional matrix of XYZ coordinates, 
#' or a set of points with class \code{\link[sp]{SpatialPoints}} or \code{\link[sp:SpatialPoints]{SpatialPointsDataFrame}}.
#'
#' @param randomborder (\code{logical}) Defaults to \code{FALSE}. If \code{TRUE}, then the points
#' falling on vertices and edges will be randomly assigned, otherwise they will be kept as \code{NA}s.
#'
#' @param output (\code{character}) Either \code{"ui"} or \code{"skeleton"}. \code{"ui"} returns the face 
#' 	names used in the user interface, while \code{"skeleton"} returns their 
#' 	indices used in back-end procedures.
#' @param ... Arguments passed to class specific methods.
#' @examples
#'	# create a grid 
#'	g <- trigrid(4)
#'	# some random points
#'	randomPoints<-rpsphere(4, output="polar")
#'	# cells
#'	locate(g, randomPoints)
#' @rdname locate
#' @exportMethod locate
setGeneric(
	name="locate",
	def=function(x,y,...){
		standardGeneric("locate")
	}

)


# locate method for the trigrid v6.0
# this version uses my own c++ function for point in tetrahedron testing
#' @rdname locate
setMethod(
	"locate", 
	signature=c(x="trigrid", y="matrix"),
	definition=function(x, y, randomborder=FALSE, output="ui"){

	#the tetrahedron algorithm does not find vertices
	if(!is.logical(randomborder)){
		stop("Invalid randomborder argument.")
	}
	
	if(!output%in%c("ui", "skeleton")){
		stop("Invalid value for output argument.")
	}
	
	#data argument
	# which formatting?
	if(ncol(y)==2){
	
		# transform the two columns
		y<-PolToCar(y, origin=x@center, radius=x@r)
	}
	
	# does the data include NAs?
	boolResultNoNA<-!is.na(y[,1]) & !is.na(y[,2]) & !is.na(y[,3]) 
	y<-y[boolResultNoNA,, drop=FALSE]
	
	
	#project the coordinates out from the origin
	#access the skeleton of the grid
		v<-x@skeleton$v*1.5
		f<-x@skeleton$f[,1:3]
		origin<-x@center
	
		d<-x@div
		
	#organize vertices to the linear coordinate + add 1s for the determinants
		#written with C++ for speed
		vtsBig<- .Call(Cpp_icosa_xyz1xyz1xyz1xyz1_, v, f)
		
		
	#check whether the point is one of the vertices!!!!! here
		
	#the queried data in a similar linear format x*n,y*n, z*n
		qrs<-.Call(Cpp_icosa_xyz1, y)
		
		nQrs<-as.integer(nrow(y))
				
			
		# allocate some memory to the results vector 
		queryIndex<-rep(-9, nrow(y))
		faceIndex<-rep(-9, nrow(y))
		foundMiddle<-rep(0, nrow(y)*12)
		faceContainer<-rep(0, max(d)+1)
		offset<-rep(0,length(d)+1)
		tempF<-rep(0, max(d)+1)
			
	#invoke the C function
	#written with pass by reference object manipulation to evade speed loss with copying
	Output = .C(Cpp_locateTriangle_,
		allVertices=as.double(vtsBig),
		divs=as.integer(d),
		nDivs=as.integer(length(d)),
		queries=as.double(qrs), 
		nQrs=as.integer(nQrs),
		queryIndex=as.integer(queryIndex),
		faceIndex=as.integer(faceIndex),
		offset=as.integer(offset),
		faceContainer=as.integer(faceContainer),
		foundMiddle=as.integer(foundMiddle),
		tempF=as.integer(tempF)
	)
	# and 1 to the 0 indexing
	fi<-Output$faceIndex+1
	qi<-Output$queryIndex+1
	
	#1. in case no values are passed to the function
	if(nQrs==0){
		return(NULL)	
	}else{
		# clean up results: indicate the points the program was unable to assign
		# delete empty entries
		fi<-fi[qi>0]
		qi<-qi[qi>0]
		
		# in case of a duplicate - just get rid of the first
		fi<-fi[!duplicated(qi)]
		qi<-qi[!duplicated(qi)]
		
		# create new container for the face indices
		newFi<-rep(NA, nQrs)
		newFi[qi] <- fi
		fi <- newFi
		qi <- 1:nQrs
	}
	
	# 2. do different stuff depending on the borders argument
	if(sum(is.na(fi))>0 & randomborder){
		# this is the more difficult case
		dubiousIndex<-which(is.na(fi))
		
		# the coordinates of these points
		weirdPoints<-y[dubiousIndex,, drop=FALSE]
		
		# repeat locate on randomly generated close points
		addFi<-apply(weirdPoints, 1, approximateFace, n=20, d=2e-8, gridObj=x, onlyOne=FALSE, output="skeleton")
		
		# add these points to the rest
		fi[dubiousIndex]<-addFi
		
	}
	
	
	if(output=="ui"){
		# translate the inner C representation to the UI
		fiUI<-x@skeleton$aF[x@skeleton$offsetF+fi]
		
		options(scipen=999)
		fiUI[!is.na(fiUI)]<-paste("F", fiUI[!is.na(fiUI)], sep="")
		options(scipen=0)
		resVec<-rep(NA, length(boolResultNoNA))
		resVec[boolResultNoNA]<-fiUI
		
		return(resVec)
	}
	if(output=="skeleton"){
		resVec<-rep(NA, length(boolResultNoNA))
		resVec[boolResultNoNA]<-fi
		
		return(resVec)
	}
}
)

# locate-method of trigrid-numeric
#' @rdname locate
setMethod(
	"locate",
	signature=c(x="trigrid", y="numeric"),
	function(x,y,...){
		# if
		if(length(y)!=2 & length(y)!=3) stop("Please provide a matrix, or vector with 2 or 3 values.")
		y <- matrix(y, nrow=1)
		locate(x, y, ...)
	}
)

# locate-method of trigrid - data.frame
#' @rdname locate
setMethod(
	"locate",
	signature=c(x="trigrid", y="data.frame"),
	function(x,y,...){
		# if
		if(ncol(y)!=2 & ncol(y)!=3) stop("Please provide a data.frame or matrix with 2 or 3 columns")
		for(i in 1:ncol(y)){
			if(!is.numeric(y[,i])) stop("One of the columns of 'x' is not numeric.")
		}
		# if all the checks are passed, pass the data frame as a matrix
		newY <- as.matrix(y)
		locate(x, newY, ...)
	}
)

# locate method of trigrid - SpatialPoints
#' @rdname locate
setMethod(
	"locate",
	signature=c(x="trigrid", y="SpatialPoints"),
	function(x,y,...){
		# if it has a proj4
		if(methods::.hasSlot(y, "proj4string")){
			# and it's not NA
			if(!is.na(y@proj4string)){
				# need rgdal
				if(requireNamespace("rgdal", quietly = TRUE)){
					y<-sp::spTransform(y, x@proj4string)@coords
				} else{
					stop("The 'rgdal' package is required to appropriately project this object. ")
				}
			}else{
				y <- y@coords 
			}
		}else{
			y <- y@coords 
		}

		locate(x, y, ...)
	}
)

# trigrid-SPDF method
#' @rdname locate
setMethod(
	"locate",
	signature=c(x="trigrid", y="SpatialPointsDataFrame"),
	function(x,y,...){
		# if it has a proj4
		if(methods::.hasSlot(y, "proj4string")){
			# and it's not NA
			if(!is.na(y@proj4string)){
				# need rgdal
				if(requireNamespace("rgdal", quietly = TRUE)){
					y<-sp::spTransform(y, x@proj4string)@coords
				} else{
					stop("The 'rgdal' package is required to appropriately project this object. ")
				}
			}else{
				y <- y@coords 
			}
		}else{
			y <- y@coords 
		}

		locate(x, y, ...)
	}
)


# locate() method for the hexagrid v6.0 - written for matrix
# This is only the y="matrix" method, inheritance will take care of the rest. 
# This version uses my own c++ function for point in tetrahedron testing
#' @param forceNA (\code{logical}) Suppressing the recursive lookup of points falling on subface boundaries.
#' @rdname locate
#' @exportMethod locate
setMethod(
	"locate", 
	signature=c(x="hexagrid",y="matrix"),
	definition=function(x, y, output="ui", randomborder=FALSE, forceNA=FALSE){

	#the tetrahedron algorithm does not find vertices
	if(!is.logical(randomborder)){
		stop("Invalid randomborder argument.")
	}
	
	if(!output%in%c("ui", "skeleton")){
		stop("Invalid value for output argument.")
	}
	
	#data argument
	# which formatting?
	if(ncol(y)==2){
		# transform the two columns
		y<-PolToCar(y, origin=x@center, radius=x@r)
	}
	
	# does the data include NAs?
	boolResultNoNA<-!is.na(y[,1]) & !is.na(y[,2]) & !is.na(y[,3]) 
	y<-y[boolResultNoNA,, drop=FALSE]
	
	#project the coordinates out from the origin
	#access the skeleton of the grid
		v<-x@skeleton$v*1.5
		f<-x@skeleton$f[,1:3]
		origin<-x@center
	
		d<-x@div
		d<-c(d,6)
		
	#organize vertices to the linear coordinate + add 1s for the determinants
		#written with C++ for speed
		vtsBig<- .Call(Cpp_icosa_xyz1xyz1xyz1xyz1_, v, f)
		
		
	#check whether the point is one of the vertices!!!!! here
		
	#the queried data in a similar linear format x*n,y*n, z*n
		qrs<-.Call(Cpp_icosa_xyz1, y)
		
		nQrs<-as.integer(nrow(y))
				
			
		# allocate some memory to the results vector 
		queryIndex<-rep(-9, nrow(y)*6)
		faceIndex<-rep(0, nrow(y)*6)
			
	#invoke the C function
	#written with direct C object manipulation to evade speed loss with copying
	Output = .C(Cpp_locateTriangle_,
		allVertices=as.double(vtsBig),
		divs=as.integer(d),
		nDivs=as.integer(length(d)),
		queries=as.double(qrs), 
		nQrs=as.integer(nQrs),
		queryIndex=as.integer(queryIndex),
		faceIndex=as.integer(faceIndex),
		offset=as.integer(rep(0,length(d)+1)),
		faceContainer=as.integer(rep(0, max(d)+1)),
		foundMiddle=as.integer(rep(0,12)),
		tempF=as.integer(rep(0, max(d)+1))
	)
	# and 1 to the 0 indexing
	fi<-Output$faceIndex+1
	qi<-Output$queryIndex+1
	
	#1. in case no values are passed to the function
	if(nQrs==0){
		return(NULL)	
	}else{
		# clean up results: indicate the points the program was unable to assign
		# delete empty entries
		fi<-fi[qi>0]
		qi<-qi[qi>0]
		
		# in case of a duplicate:
		if(sum(duplicated(qi))>0){
			# delete both and reinvestigate
			tqi<-table(qi)
			duplicateBullshit<-as.numeric(names(tqi[tqi>1]))
			fi<-fi[!qi%in%duplicateBullshit]
			qi<-qi[!qi%in%duplicateBullshit]
		}
	
		# create new container for the face indices
		newFi<-rep(NA, nQrs)
		newFi[qi] <- fi
		fi <- newFi
		qi <- 1:nQrs
	}
	
	if(output=="ui"){

		# translate the inner C representation to the UI
		fiUI<-x@skeleton$aSF[Output$offset[length(d)]+fi]
		
		# add the labels
		#temporarily supress scientific notation
		options(scipen=999)
		fiUI[!is.na(fiUI)]<-paste("F", fiUI[!is.na(fiUI)], sep="")
		options(scipen=0)
		
		# this section needs to be here, otherwise it won't recognize the same faces separated to different subfaces
		# 2. do different stuff depending on the borders argument
		# stop it for base case of recursion
		if(!forceNA){
			if(sum(is.na(fiUI))>0){
				# this is the more difficult case
				dubiousIndex<-which(is.na(fiUI))
				
				# the coordinates of these points
				weirdPoints<-y[dubiousIndex,, drop=FALSE]
				
				# repeat locate on randomly generated close points
				addFiUI<-apply(weirdPoints, 1, approximateFace, n=20, gridObj=x, d=2e-10, onlyOne=!randomborder, output="ui")
				addFiUI
				
				# add these points to the rest
				fiUI[dubiousIndex]<-addFiUI
				
			}
		}
		resVec<-rep(NA, length(boolResultNoNA))
		resVec[boolResultNoNA]<-fiUI
		
		return(resVec)
	}
	if(output=="skeleton"){
		fiInner<-x@skeleton$f[Output$offset[length(d)]+fi,1]
		
		# stop it for base case of recursion
		if(!forceNA){
			if(sum(is.na(fiInner))>0){
				# this is the more difficult case
				dubiousIndex<-which(is.na(fiInner))
				
				# the coordinates of these points
				weirdPoints<-y[dubiousIndex,, drop=FALSE]
			
				# repeat locate on randomly generated close points
				addFiInner<-apply(weirdPoints, 1, approximateFace, n=20, gridObj=x, d=2e-10, onlyOne=!randomborder, output="skeleton")
				
				# add these points to the rest
				fiInner[dubiousIndex]<-addFiInner
				
			}
		}
	
		resVec<-rep(NA, length(boolResultNoNA))
		resVec[boolResultNoNA]<-fiInner
		
		return(resVec)
	}
	}
)




#' Position of face centers and vertices on a grid
#' 
#' This function will retrieve the position of a vertex or a face on a \code{\link{hexagrid}} or \code{\link{trigrid}} object.
#' 
#' Vertex and face names can be mixed in a single \code{names} argument.
#' 
#' @param gridObj a (\code{\link{hexagrid}} or \code{\link{trigrid}}) Icosahedral grid object.
#' 
#' @param names (\code{character}) Vector of the names that are to be looked up.
#' 
#' @param output (\code{character}) The coordinate system in which the names are to be shown: use \code{"polar"} for longitude-latitude and \code{"cartesian"} for XYZ output.
#' 
#' @return A \code{numeric} matrix.
#' 
#' @examples
#' g <- trigrid(c(4,4))
#' pos(g, c("F2", "P6", "dummyname"))
#' 
#' 
#' @export
pos<-function(gridObj, names, output="polar"){
	
	if(!class(gridObj)%in%c("trigrid", "hexagrid")) stop("Invalid gridObj argument.")
	if(!output%in%c("polar", "cartesian")) stop("Invalid output argument.")
	
	names<-as.character(names)
	
	#facecenters
	fBool<-names%in%rownames(gridObj@faces)
	fcs<-gridObj@faceCenters[names[fBool],]

	#vertices
	vBool<-names%in%rownames(gridObj@vertices)
	vs<-gridObj@vertices[names[vBool],]
	
	#result
	res<-matrix(NA, nrow=length(names), ncol=3)
	
	res[fBool,] <- fcs
	res[vBool,] <- vs
	
	if(output=="cartesian"){
		rownames(res)<-names
		colnames(res)<-c("x","y", "z")
	}
	
	if(output=="polar"){
		res<-CarToPol(res, norad=TRUE, origin=gridObj@center)
		rownames(res)<-names
		colnames(res)<-c("long", "lat")
	}
	return(res)
}


# small, fast utility function for the lookup of vertices (used in lookup!)
whichVertices<-function(vertices, data){
	# result if no vertex found
	vertIndex<-NULL
	
	#the only total check 
	ctrl1<-data[,1]%in%vertices[,1]
	
	#first coordinate match
	if(sum(ctrl1)>0){
		datSub<-data[ctrl1,,drop=FALSE]
		indSub<-which(ctrl1)
		ctrl2<-datSub[,2]%in%vertices[,2]
		
		#second coordinate match too
		if(sum(ctrl2)>0){
			indSub2<-indSub[ctrl2]
			datSub2<-datSub[ctrl2,,drop=FALSE ]
			ctrl3<-datSub2[,3]%in%vertices[,3]
			
			#third matches as well: vertex!
			if(sum(ctrl3)>0){
				vertIndex<-indSub2[ctrl3]
			}
		}
	
	}

	return(vertIndex)
}




#only one for the hexagrid: subface boundaries should give back an entry even if the randomborder is FALSE
approximateFace<-function(coords, n, d, gridObj, onlyOne=FALSE, output="skeleton"){
	# test variables
#	d<-10e-6
#	n<-15
	# matrix of random coordinates around this point
	randMat<-cbind(coords[1]+stats::rnorm(n,0,d), coords[2]+stats::rnorm(n,0,d), coords[3]+stats::rnorm(n,0,d))
	
	if(class(gridObj)=="hexagrid"){
		temp<-suppressWarnings(locate(gridObj, randMat, output=output, randomborder=FALSE, forceNA=TRUE))
	}else{
		temp<-suppressWarnings(locate(gridObj, randMat, output=output, randomborder=FALSE))
	}
	temp<-unique(temp[!is.na(temp)])
	
	if(onlyOne){
		if(length(temp)==1){
			return(temp)
		}else{
			return(NA)
		}
	}else{
		return(sample(temp,1))
	}
}
		
