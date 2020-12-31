#' Calculation of point coordinates along an arc
#' 
#' This function calculates points along an arc between two points and a circle center.
#' 
#' The function always returns the smaller arc, with angle alpha < pi.
#' 
#' @param p1 (\code{numeric}) Vector, XYZ or longitude-latitude coordinates of the first point along the arc.
#' 
#' @param p2 (\code{numeric}) Vector, XYZ or longitude-latitude coordinates of the last point along the arc.
#' 
#' @param origin (\code{numeric}) vector, The center of the circle in XYZ coordinates (default is \code{c(0,0,0)}).
#' 
#' @param breaks (\code{integer}) The number of points inserted between \code{p1} and \code{p2}. Has to be positive.
#' 
#' @param onlyNew (\code{logical}) Should \code{p1} and \code{p2} be omitted from the result?
#' 
#' @param output (\code{character}) The coordinate system of the output points. Can either be \code{"polar"} for
#' 	longitude-latitude or \code{"cartesian"} for XYZ data.
#'
#' @param radius (\code{numeric}) Single value, the radius of the circle in case the input points have only polar coordinates.
#'	Unused when XYZ coordinates are entered. Defaults to the authalic radius of Earth ca. 6371.007km.
#' 
#' @return Either an XYZ or a long-lat numeric matrix.
#' 
#' @examples
#'	# empty plot
#'	plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90,90))
#'	# then endpoints of the arc
#'	point1<-c(-45,-70)
#'	point2<-c(130,65)
#'	points(arcpoints(point1, point2, breaks=70, output="polar"))
#'
#' @export arcpoints
arcpoints<-function(p1,p2,breaks=2,origin=c(0,0,0), onlyNew=FALSE, output="cartesian", radius=authRadius){
	if(!output%in%c("cartesian", "polar")) stop("Invalid \'output\' argument.")
	
	if(!is.logical(onlyNew)| length(onlyNew)>1) "Invalid \'onlyNew\' argument."

	if(!is.numeric(breaks) | length(breaks)!=1 | breaks%%1!=0) stop("Invalid \'breaks\' argument")
	
	if(!is.numeric(p1) | !is.numeric(p2) | !is.numeric(origin)) stop("Invalid coordinate input.")
	
	if(!sum(is.finite(p1))%in%c(2,3) | 
		!sum(is.finite(p2))%in%c(2,3) | 
		!sum(is.finite(origin))%in%c(2,3)) stop("Invalid coordinate input.")
	
	if(length(p1)==2){
		p1<-PolToCar(matrix(p1, ncol=2, nrow=1),radius=radius, origin=origin)
		p1<-p1[1,]
	}
	
	if(length(p2)==2){
		p2<-PolToCar(matrix(p2, ncol=2, nrow=1),radius=radius, origin=origin)
		p2<-p2[1,]
	}
		
	#invoke Rcpp function
	temp<-.Call(Cpp_icosa_SplitArc_, p1, p2, origin, breaks, onlyNew)
	
	colnames(temp)<-c("x", "y", "z")
	
	if(output=="polar"){
		temp<- CarToPol(temp, norad=TRUE, origin=origin )
	}
	
	return(temp)
}

#' Calculation of distances along arcs
#' 
#' This function calculates the shortest arc distance between two points.
#' 
#' @param p1 (\code{numeric}) Vector, XYZ or longitude-latitude coordinates of the first point along the arc.
#' 
#' @param p2 (\code{numeric}) Vector, XYZ or longitude-latitude coordinates of the last point along the arc.
#' 
#' @param origin (\code{numeric}) Vector, the center of the circle in XYZ coordinates (default is \code{c(0,0,0)}).
#' 
#' @param output (\code{character}) The type of the output value. \code{"distance"} will give the distance
#'	in the metric that was fed to the function for the coordinates or the radius.
#'	\code{"deg"} will output the the distance in degrees, \code{"rad"} will do
#'	so in radians.
#'
#' @param radius (\code{numeric}) The radius of the circle in case the input points have polar coordinates only.
#'	Unused when XYZ coordinates are entered. Defaults to the authalic radius of Earth ca. 6371.007km.
#'
#' @return A single \code{numeric} value.
#'
#' @examples 
#'	# coordinates of two points
#'	point1<- c(0,0)
#'	point2<- c(180,0)
#'	arcdist(point1,point2,"distance")
#' @export arcdist
arcdist <- function(p1, p2, output="distance", origin=c(0,0,0), radius=authRadius) {
    if(!is.numeric(p1) | !is.numeric(p2) | !is.numeric(origin)) stop("Invalid coordinate input.")
	
	if(!output%in%c("distance", "deg", "rad")) stop("Invalid \'output\' argument.")
	
	if(!sum(is.finite(p1))%in%c(2,3) | 
		!sum(is.finite(p2))%in%c(2,3) | 
		!sum(is.finite(origin))%in%c(2,3)) stop("Invalid coordinate input.")
	
	if(length(p1)==2){
		p1<-PolToCar(matrix(p1, ncol=2, nrow=1),radius=radius, origin=origin)
		p1<-p1[1,]
	}
	
	if(length(p2)==2){
		p2<-PolToCar(matrix(p2, ncol=2, nrow=1),radius=radius, origin=origin)
		p2<-p2[1,]
	}
		
	if(output=="distance") method<-T
	if(output%in%c("deg", "rad")) method<-F
	
	result<-.Call(Cpp_icosa_ArcDist_, p1, p2, origin, method)
	
	#acos() out of domain error
		if(is.nan(result) & output=="distance"){
			v<-p1-origin
			result <- pi*sqrt(v[1]^2+v[2]^2+v[3]^2)
		} 
		
		if(is.nan(result) & output%in%c("deg","rad")){
			v<-p1-origin
			result <- pi
		} 
	
	if(output=="deg"){
		result<-result*180/pi
	}

	return(result)
	
}

#' Calculation of distance matrices along arcs
#' 
#' This function calculates the shortest arc distance matrix between two sets of points.
#' 
#' This function will create all possible shortest arc distances between points in the two sets,
#' 	but not between the points within the sets. The function is useful for great circle distance calculations.
#' 	For a symmetrical distance matrix leave the \code{points2} argument empty.
#' 
#' @param points1 (\code{numeric}) Matrix, XYZ or longitude-latitude coordinates of the first set of points.
#' 
#' @param points2 (\code{numeric}) Matrix, XYZ or longitude-latitude coordinates of the second set of points. 
#'	Leave this empty if you want all the arc distances between a set of points	
#' 
#' @param origin (\code{numeric}) Vector, the center of the circle in XYZ coordinates (default is \code{c(0,0,0)}).
#' 
#' @param output (\code{character}) The type of the output value. \code{"distance"} will give back the distance
#' 	in the metric that was fed to the function in the coordinates or the radius.
#' 	\code{"deg"} will output the the distance in degrees, \code{"rad"} will do
#' 	so in radians.
#' 
#' @param radius (\code{numeric}) The radius of the circle in case the input points have polar coordinates only.
#' 	Unused when XYZ coordinates are entered. Defaults to the authalic radius of Earth ca. 6371.007km.
#' 
#' @return A single \code{numeric} value.
#' 
#' @examples
#' g <- trigrid(c(4))
#' res <- arcdistmat(g@vertices)
#' 
#' rand<-rpsphere(500)
#' res2 <- arcdistmat(g@vertices, rand)
#'
#'	@export
arcdistmat<-function(points1, points2=NULL, origin=c(0,0,0), output="distance", radius=authRadius){
	# output argument
	if(!output%in%c("distance", "deg", "rad")) stop("Invalid \'output\' argument.")
	if(output=="distance") method<-T
	if(output%in%c("deg", "rad")) method<-F
	
	present<-T
	if(is.null(points2)){
		points2<-points1
		present<-F
	}
	
	if(!is.numeric(points1) | !is.numeric(points2) | !is.numeric(origin)) stop("Invalid coordinate input.")
	
	
	if(!ncol(points1)%in%c(2,3) | 
		!ncol(points2)%in%c(2,3) | 
		!sum(is.finite(origin))%in%c(2,3)) stop("Invalid coordinate input.")
	
	if(sum(is.na(points1))>0 | sum(is.na(points2))>0) stop("The coordinates include NAs")
	if(ncol(points1)==2){
		points1<-PolToCar(matrix(points1, ncol=2),radius=radius, origin=origin)
		points1<-points1[1,]
	}
	
	
	
	if(ncol(points2)==2){
		points2<-PolToCar(matrix(points2, ncol=2),radius=radius, origin=origin)
		points2<-points2[1,]
	}
		
	if(present){
		distMat<- .Call(Cpp_icosa_ArcDistMat_, points1, points2, origin, method)
	}else{
		distMat<- .Call(Cpp_icosa_SymmetricArcDistMat_, points1, origin, method)
	}
	
	rownames(distMat)<-rownames(points1)
	colnames(distMat)<-rownames(points2)
	
	#acos() out of domain error
		if(sum(is.nan(distMat))>0 & output=="distance"){
			v<-points1[1,]-origin
			distMat[is.nan(distMat)] <- pi*sqrt(v[1]^2+v[2]^2+v[3]^2)
		} 
		
		if(sum(is.nan(distMat))>0 & output%in%c("deg","rad")){
			v<-points1[1,]-origin
			distMat[is.nan(distMat)] <- pi
		} 

		
	if(output=="deg"){
		distMat<-distMat*180/pi
	}
	return(distMat)
}



#function to rotate a single point
rotateOnePoint<-function(coords, angles,origin)
{
	#coords<-c(0,1,0)
	#angles<-c(pi/2,pi/2,pi)
	
	#the rotation matrix
	rotMat<-function(theta)
	{
		mat<-matrix(NA,ncol=2,nrow=2)
		mat[1,1]<-cos(theta)
		mat[1,2]<--sin(theta)
		mat[2,1]<-sin(theta)
		mat[2,2]<-cos(theta)
		return(mat)
	}
	
	#location vector
	locVec<-coords-origin
			
	#first rotation
	#around x
	xMat<-matrix(locVec[2:3],ncol=1, nrow=2)
	xMat<-rotMat(angles[1])%*%xMat
	locVec<-c(locVec[1],as.numeric(xMat))
	
	#second rotation
	yMat<-matrix(locVec[c(1,3)],ncol=1, nrow=2)
	yMat<-rotMat(angles[2])%*%yMat
	locVec<-c(yMat[1], locVec[2],yMat[2])
	
	#third rotation
	zMat<-matrix(locVec[c(1,2)],ncol=1, nrow=2)
	zMat<-rotMat(angles[3])%*%zMat
	locVec<-c(zMat[1:2], locVec[3])
	
	return(locVec)

}
	

# function to create random points on the sphere
#' Random point generation on the surface of a sphere
#' 
#' This function will create a predefined number of points randomly distributed
#' on the surface of a sphere with a given radius.
#' 
#' The function uses a three dimension normal distribution to generate points, 
#' which are then projected to the surface of the sphere.
#' 
#' @param n (\code{numeric}) The number of random points to be created.
#' 
#' @param radius (\code{numeric}) The radius of the sphere
#' 
#' @param origin (\code{numeric}) The center of the sphere (XYZ coordinates).
#' 
#' @param output (\code{character}) The coordinate system of the new points. Can either be 
#'	\code{"cartesian"} for XYZ coordiates or \code{"polar"} for spherical, 
#'	longitude-latitudes coordinates.
#'
#' @return A 3-column (XYZ) or a 2-column (long-lat) \code{numeric} matrix.
#' 
#' @examples
#'  randomPoints <- rpsphere(2000, output="polar")
#' # observe latitudinal pattern
#'  plot(randomPoints, xlim=c(-180, 180), ylim=c(-90, 90))
#' 
#' @export
rpsphere <- function(n=1, output="cartesian", radius=authRadius, origin=c(0,0,0)){
	if(!is.numeric(radius) |
	length(radius)!=1) stop("Invalid input for argument \'radius\'.")
	
	if(!is.numeric(n) |
	length(n)!=1) stop("Invalid input for argument \'n\'.")

	if(!is.numeric(origin) |
	length(origin)!=3) stop("Invalid input for argument \'origin\'.")
	
	
	if(!output%in%c("polar", "cartesian")) 
	
	stop("Invalid input for argument \'output\'.")
		#random variables
		x1 <- stats::rnorm(n, 0, 1)
		y1 <- stats::rnorm(n, 0, 1)
		z1 <- stats::rnorm(n, 0, 1)
		origPoints<-cbind(x1,y1,z1)
	
	# location vectors
	vectors<-origPoints
	vectors[,1]<-origPoints[,1]
	vectors[,2]<-origPoints[,2]
	vectors[,3]<-origPoints[,3]
	
	# distances from the origin
	dists<-sqrt(vectors[,1]^2+vectors[,2]^2+vectors[,3]^2)
	
	#project the point to the sphere
	newPoints<-origPoints
	newPoints[,1]<-origPoints[,1]*radius/dists
	newPoints[,2]<-origPoints[,2]*radius/dists
	newPoints[,3]<-origPoints[,3]*radius/dists
		
	# column names
	colnames(newPoints)<-c("x", "y", "z")
	
	#outputs
	if(output=="polar"){
		result<-CarToPol(newPoints, norad=TRUE, origin=origin)
		colnames(result) <- c("long", "lat")
	}else{
		result<-newPoints
		result[,1] <- result[,1]+origin[1]
		result[,2] <- result[,2]+origin[2]
		result[,3] <- result[,3]+origin[3]
	}
	
	return(result)
}





#' Surface centroid point of a spherical point cloud
#' 
#' This function the projected place of the centroid from a pointset on the sphere.
#' 
#' The function implements great circle calculations to infer on the place of the centroid, which makes it resource demanding. This is necessary
#'	to avoid a particual error that frequently occurrs with other methods for centroid calculation, namely that the place of the centroid is right,
#' 	but on the opposite hemisphere.
#' 
#' @param x (\code{matrix} or \code{data.frame}) Numeric data, XYZ or longitude-latitude coordinates of the set of points.
#' 
#' @param output (\code{character}) The coordinate system of the output points. Can either be \code{"polar"} for
#' 	longitude-latitude or \code{"cartesian"} for XYZ data.
#'
#' @param center (\code{numeric}) The center of the sphere in XYZ coordinates (default is 0,0,0).
#' 
#' @param radius (\code{numeric}) The radius of the circle in case the input points have only polar coordinates.
#'	Unused when XYZ coordinates are entered. Defaults to the authalic radius of Earth ca. 6371.007km.
#'
#' @param ... Arguments passed to the \code{matrix}-method.
#' @return Either an XYZ or a long-lat \code{numeric} vector.
#' 
#' @examples
#'	# generate some random points
#'	allData <- rpsphere(1000)
#'	# select only a subset
#'	points<-allData[allData[,2]>1500,]
#' # transform to 2d
#'  points2 <- CarToPol(points, norad=TRUE)
#'	# the spherical centroid
#'	sc <- surfacecentroid(points2, output="polar")
#'	sc
#'	
#'	#3d plot
#'	plot(points2, xlim=c(-180, 180), ylim=c(-90, 90))
#'	points(sc[1], sc[2], col="red", cex=5)
#'
#' @exportMethod surfacecentroid
#' @rdname surfacecentroid
setGeneric(
	"surfacecentroid", 
	function(x,...) standardGeneric("surfacecentroid")
)

#' Matrix-method of surfacecentroid()
#' @rdname surfacecentroid
setMethod(
	"surfacecentroid",
	signature=c(x="matrix"), 
	function(x, output="polar", center=c(0,0,0), radius=authRadius){
		if(nrow(x)<2) return(x)
		#data argument
		# which formatting?
		if(ncol(x)==2){
			# transform the two columns
			x<-PolToCar(x, origin=center, radius=radius)
		}
		if (ncol(x)==3){
			radVec<-x[1,]-center
			rad<-sqrt(radVec[1]^2+radVec[2]^2+radVec[3]^2)
		}

		#the 3d centroid of the point cloud
			centroid3d<-apply(x, 2, mean, na.rm=TRUE)
			if(output=="cartesian"){
				radVec<-(centroid3d-center)
				retCentroid<-centroid3d/(sqrt(radVec[1]^2+radVec[2]^2+radVec[3]^2))*rad
				return(retCentroid)
			}
		
		#transform back to spherical coordinates
		if(output=="polar"){
			#the longitude problem!!!
			xSign<-sign(centroid3d[1])
		
			theta<-atan(centroid3d[2]/centroid3d[1])
			phi<-atan(sqrt(centroid3d[1]^2+centroid3d[2]^2)/centroid3d[3])
		
		#transform spherical coordinates to long/lat
			theta<-theta/pi*180
			phi<-phi/pi*180
		
			#convert to lat-long
			if(phi>=0) lat<-90-phi
			if(phi<0) lat<--90-phi
			
			if(xSign<0 & theta<=0) long<-180+theta
			if(xSign<0 & theta>0) long<--180+theta
			if(xSign>=0) long<-theta
			
			#return value
			centroidLongLat<-c(long, lat)
			names(centroidLongLat)<-c("long", "lat")
		
		
		#return value
			return(centroidLongLat)
		}		
	}
)

#' df-method of surfacecentroid()
#' @rdname surfacecentroid
setMethod(
	"surfacecentroid", 
	signature=c(x="data.frame"),
	function(x,...){
		newX <- as.matrix(x)
		surfacecentroid(newX,...)
	}
)

#' SP-method of surfacecentroid()
#' @rdname surfacecentroid
setMethod(
	"surfacecentroid",
	signature=c(x="SpatialPoints"),
	function(x,...){
		# if it has a proj4
		if(methods::.hasSlot(x, "proj4string")){
			# and it's not NA
			if(!is.na(x@proj4string)){
				# need rgdal
				if(requireNamespace("rgdal", quietly = TRUE)){
					x<-sp::spTransform(x, sp::CRS("+proj=longlat"))@coords
				} else{
					stop("The 'rgdal' package is required to appropriately project this object. ")
				}
			}else{
				x <- x@coords 
			}
		}else{
			x <- x@coords 
		}

		locate(x, ...)
	}
)




#' Spherical convex hull. 
#' 
#' This function calculates a possible implementation of the spherical convex hull.
#' 
#' With the method \code{centroidprojection} the function calls the \code{\link{surfacecentroid}} 
#'	function to get the a reference point from the shape. Then all the points are 'projected' 
#'	close to this point using the great circles linking them to the reference point.
#'	Each such great circle will be devided to an equal number of points and the closest
#'	 will replace the original point coordinates in the convex hull algorithm implemented in \code{\link[grDevices]{chull}}. 
#' 
#' @param data  (\code{numeric}) Matrix, XYZ or longitude-latitude coordinates of the set of points.
#' 
#' @param center (\code{numeric}) Vector, The center of the sphere in XYZ coordinates (default is 0,0,0).
#' @param radius (\code{numeric}) Single value, indicating the radius of the sphere. Defaults to the R2 radius of Earth (6371.007km).
#' @param param (\code{numeric}) Single positive integer, indicates the number of divisions in the centroid projection method. The higher the number, the closer the replacement points are to the centroid.
#' 
#' @param strict (\code{logical}) Strictly convex output is required.
#' @return The indices of the data points forming the convex hull as a (\code{numeric}) vector.
#' 
#' @examples
#'	# generate some random points
#'	allData <- rpsphere(1000)
#'	# select only a subset
#'	points<-allData[allData[,1]>3000,]
#'	chullsphere(points)
#'	
#'
#' @export chullsphere
chullsphere<-function(data, center=c(0,0,0), radius=authRadius, param=200, strict=TRUE)
{
	if(ncol(data)==2){
		# transform the two columns
		data<-PolToCar(data, radius, origin=center)
	}
	if (ncol(data)==3){
		radVec<-data[1,]-center
		rad<-sqrt(radVec[1]^2+radVec[2]^2+radVec[3]^2)
	}
	
	#calculate the group centroid
		centroid<-surfacecentroid(data, output="cartesian", center=center, radius)

	#shrink the 2d surface proportionally! to the vicinity of the reference point of latLong 2d space(roughly planar area)
	projectedPoints<-.Call(Cpp_icosa_projectCloseToPoint_, data, centroid, center, param)
	
	#omit NA's (including the case where the centroid is among the points!)
	boolMiss<-is.na(projectedPoints[,1])& is.na(projectedPoints[,2])
	if(any(boolMiss)){
		for(i in which(boolMiss)){
			projectedPoints[i,] <- centroid
		}
	}
	
	projP<- CarToPol(projectedPoints,norad=TRUE)
	convHull <- grDevices::chull(projP)
	
	# do another quick check based on the dihedral angles, as some points are seen as convex from the centroid
	if(strict){
		# do iteratively!
		more <- TRUE
		while(more){
			# test results
			testIndex <- rep(FALSE, length(convHull))
	
			# for all hull points
			for(i in 1:length(convHull)){
				if(i==1){
					one <- length(convHull)
				}else{
					one <- i-1
				}
				if(i==length(convHull)){
					three <- 1
				}else{
					three<-i+1
				}
				# the first dihedral angle
				xa <- arcdist(centroid, data[convHull[one],], center, output="rad")
				ab <- arcdist(data[convHull[one],], data[convHull[i],], center, output="rad")
				bx <- arcdist(data[convHull[i],], centroid, center, output="rad")
		
				firstAng = acos((cos(xa)-(cos(ab)*cos(bx)))/(sin(ab)*sin(bx)));
				
				# the second dihedral angle
				xc <- arcdist(centroid, data[convHull[three],], center, output="rad")
				cb <- arcdist(data[convHull[three],], data[convHull[i],], center, output="rad")
		
				secondAng = acos((cos(xc)-(cos(cb)*cos(bx)))/(sin(cb)*sin(bx)));
				
				if((secondAng+firstAng)<=pi) testIndex[i] <- TRUE
			}
	
			# result
			if(sum(testIndex) < length(convHull) & sum(testIndex)>2){
				more<-TRUE
			}else{
				more <-FALSE
			}
			
			convHull <- convHull[testIndex]
		}
	}

	return(convHull)

}


#	#' Surface area of a spherical convex-hull defined by a set of points
#	#' 
#	#' The function returns the area covered by the spherical convex hull of a pointset in square kilometers
#	#' 
#	#' The function assumes a round Earth, but SpatialPoints with \code{CRS} entries will be projected to a sphere before the calulations. 
#	#' 
#	#' @param data Coordinates of individual points. Can be either a two-dimensional 
#	#' matrix of long-lat coordinates, a three-dimensional matrix of XYZ coordinates, 
#	#' or a set of points with class 'SpatialPoints'.
#	#' @param origin Numeric vector of length 3, defining the center of the sphere. Defaults to c(0,0,0).
#	#' @param radius Numeric value, the radius of the circle in case the input points have only polar coordinates.
#	#' @examples
#	#' # simple example with a hexagrid 
#	#' a<-hexagrid(c(4), sp=T)
#	#' b<-a[c(lomax=40, lomin=-40, lamax=30, lamin=-30)]
#	#' data <- centers(b)
#	#' surfacechullsphere(data)
#	#' 
#	#' @export
#	surfacechullsphere <-function(data, origin=c(0,0,0), radius=authRadius){
#		
#		# 1. make sure you have cartesian coordinates
#		
#		# for the SpatialPoints
#		if(class(data)=="SpatialPoints"){
#			# if it has a proj4
#			if(methods::.hasSlot(data, "proj4string")){
#				# and it's not NA
#				if(!is.na(data@proj4string)){
#					# need rgdal
#					if(requireNamespace("rgdal", quietly = TRUE)){
#						data<-sp::spTransform(data, gridObj@proj4string)@coords
#					} else{
#						stop("The rgdal package is required to appropriately project this object. ")
#					}
#				}
#			}
#		}
#		
#		#data argument
#		# which formatting?
#		if(ncol(data)==2){
#			# transform the two columns
#			data<-PolToCar(data, radius, origin)
#		}
#		
#	
#		# 2. calculate surface centroid
#		surfcent <- surfacecentroid(data, output="cartesian", center=origin, radius)
#	
#		# 3. calculate spherical convex hulls
#		indices <- rev(chullsphere(data, center=surfcent, radius))
#	
#		# 4. calculate the area
#		if(length(indices)>2){
#			surfarea <-.Call(Cpp_icosa_surfConvHullTri,
#				data[indices,],
#				surfcent,
#				origin,
#				pi)
#		}else{
#			surfarea <- 0
#		}
#		return(surfarea)
#	}
