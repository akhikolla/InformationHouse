blankSphere <- function (x, y=NULL, z=NULL, res=100, radius = authRadius, color="white", add=FALSE, ...) {
	if(!requireNamespace("rgl", quietly = TRUE)) stop("Install the 'rgl' package and reload 'icosa' to use this function.")
	
	# lat and long values
	lat <- matrix(seq(90, -90, len = res)*pi/180, res, res, byrow = TRUE)
	long <- matrix(seq(-180, 180, len = res)*pi/180, res, res)
	
	# the center of the sphere
	center <-c(x,y,z)
	
	add2 <- add
	x <- center[1] + radius*cos(lat)*cos(long)
	y <- center[2] + radius*cos(lat)*sin(long)
	z <- center[3] + radius*sin(lat)
	rgl::persp3d(x, y, z, specular="white", add=add2, color=color, ...)
	
}



#' Guides for 3d spherical plotting.
#' 
#' This function plots 3d guidelines for navigation on the surface of the sphere,
#' 	 includings the rotational axis and a polar coordinate system.
#' 
#' The function is built on the openGL renderer of the R package \code{rgl}.
#'  
#' @param axis (\code{numeric}) Draws the -90(lat. deg. ) +90 (lat. deg.) axis. The plotted radius will be \code{axis} times the authalic radius, ca. 6371km.
#' 
#' @param polgrid (\code{numeric}) with the length of \code{2}, where the first argument specifies
#' the size of the longitudinal and the second the latitudinal divisions (degrees). Setting this argument to \code{NULL} will turn this feature off.
#'
#' @param res (\code{numeric}) Graphical resolution of the curves:
#' the distance in degrees between the points of the rendered guides. 
#'	@param textPG (\code{logical}) Flag indicating whether the coordinate values should be added to the 3d render.
#' 
#' @param origin (\code{numeric}) Vector of length=3. Indicates the center of the guiding sphere.
#' @param radius (\code{numeric}) Values indicating the radius of the guiding sphere. Defaults to the R2 radius of Earth (6371.007km). 
#' @param drad (\code{numeric}) Value, indicates the position of coordinate 3d text relative to the guiding sphere radius.
#' @param ... Additional arguments passed to \code{\link[rgl:3dobjects]{segments3d}}, \code{\link[rgl:3dobjects]{lines3d}} and \code{\link[rgl:texts]{text3d}}.
#' @return The function does not return any value.
#'
#' @examples
#' # create a hexagonal grid
#'   g <- hexagrid(c(2,2))
#' # plot the grid in 3d space
#'	 plot3d(g, guides=FALSE)
#' # plot the rotational axis in blue
#'   guides3d(axis=2, polgrid=NULL, col="blue")
#' # plot the polar grid at 10 degree resolution
#'   guides3d(axis=NULL, polgrid=c(10,10), col="red")
#' # plot some coordinates
#'   guides3d(axis=NULL, polgrid=c(30,30), textPG=TRUE, col="orange", cex=1.4)
#' @export
guides3d <- function(axis=1.5, polgrid=c(30,30), textPG=FALSE, res=1,  origin=c(0,0,0), radius=authRadius, drad=1.1, ...){
	if(!requireNamespace("rgl", quietly = TRUE)) stop("Install the 'rgl' package and reload 'icosa' to use this function.")
	if(!is.null(axis)){
		rgl::segments3d(x=c(0,0)+origin[1], y=c(0,0)+origin[2], z=c(-axis*radius,axis*radius)+origin[3],...)
		rgl::segments3d(x=c(200,0)+origin[1], y=c(0,0)+origin[2], z= c(axis*radius-500,axis*radius)+origin[3],...)
		rgl::segments3d(x=c(-200,0)+origin[1], y=c(0,0)+origin[2], z= c(axis*radius-500,axis*radius)+origin[3],...)
	}
	if(!is.null(polgrid[1])){

		if(360%%polgrid[1]!=0) 
			stop(paste("360 is not divisble by ", polgrid[1],sep=""))
		if(180%%polgrid[2]!=0) 
			stop(paste("180 is not divisble by ", polgrid[2],sep=""))
		
		#meridians
		#division
		usedLongs<-seq(-180,180,polgrid[1])
		a<-usedLongs
		#resolution
		b<-c(seq(-90,90,res), NA)
		
		lngs<-rep(a,each=length(b))
		lngs[1:length(a)*length(b)]<-NA
		
		lats<-rep(b, length(a))
		merid<-cbind(lngs, lats)
		merid3d<-PolToCar(merid,origin=origin, radius=radius)
		
		#lat circles
		#division
		usedLats<-seq(-90,90,polgrid[2])
		a<-usedLats
		#resolution
		b<-c(seq(-180,180,res),NA)
		
		
		lats<-rep(a,each=length(b))
		lats[1:length(a)*length(b)]<-NA
		
		lngs<-rep(b, length(a))
		latC<-cbind(lngs, lats)
		
		latC3d<-PolToCar(latC,origin=origin, radius=radius)
		
		rgl::lines3d(merid3d,...)
		rgl::lines3d(latC3d,...)
		
		if(textPG){
			# the latitudes
			latTab<-cbind(rep(0,length(usedLats)), usedLats)
			coordLat<-PolToCar(latTab, origin=origin, radius=radius)*drad
			rgl::text3d(coordLat, text=usedLats, ...)
		
			# the longitudes
			# you need only 180, no -180
			usedLongs<- usedLongs[-length(usedLongs)]
			longTab<-cbind(usedLongs,rep(0,length(usedLongs)))
			coordLong<-PolToCar(longTab, origin=origin, radius=radius)*drad
			rgl::text3d(coordLong, text=usedLongs, ...)
		
		}
		
	}
	
}
