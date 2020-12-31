#' Resampling a facelayer to a different resolution
#' 
#' The function applies different resampling algorithms. Currently there are only two implemented methods, one for upscaling and one for downscaling. The downscaling method "average" will tabluate all face centers from the high resolution grid that fall on a coarse resolution cell and average them. The upscaling method "ebaa" (edge breakpoint area approximation) will estimate the areas covered by the high resolution cells using the number of edge breakpoints.
#' 
#' @param x (\code{\link[raster:raster]{RasterLayer}}, \code{\link{facelayer}}) Object to resample.  
#' 
#' @param y  (\code{\link{hexagrid}} or \code{\link{trigrid}}) Object describing the target structure.
#' @param res (\code{numeric}) Value indicating the precision of area estimation during the upscaling (\code{facelayer}-method). In case the \code{"ebaa"} method is chosen, the variable indicate the number of breaking points on an edge.
#' 
#' @param method (\code{character}) The name of the algorithm used for resampling. 
#' 
#' @return A named \code{numeric} vector.
#' 
#' @examples
#' # create a grid
#' g <- trigrid(c(4,4))
#' # create a data layer
#' fl <- facelayer(g)
#' fl@values<-rnorm(length(fl))
#' # target structure
#' h <- trigrid(4)
#' # resampling
#' res <- resample(fl, h)
#' fl2<-facelayer(h)
#' fl2@values[] <- res
#'
#' @exportMethod resample
#' @rdname resample
setMethod(
	"resample",
	signature=c("facelayer", "trigrid"),
	definition=function(x, y,method=NULL,res=5){
	#x:layer, y new grid
		oldGrid<-get(x@grid)
		values<-x@values
		#only if values is a numeric vector!
		if(!is.numeric(x@values)){
			
		}
		
		# do some checking immediately: upscaling or downscaling
		# downscaling
		if(length(oldGrid)>length(y)){
			downscale <- TRUE
			if(is.null(method)){
				method <- "average"
			}else{
				# perform a check of the potential resampling methods
				if(!method%in%c("average")){
					stop("Invalid downscaling method. ")
				}
			}
			
		} else {
			downscale <- FALSE
			if(is.null(method)){
				method <- "ebaa"
			}else{
				# perform a check of the potential resampling methods
				if(!method%in%c("ebaa")){
					stop("Invalid upscaling method.")
				}
			}
			
		}
		
		# do some checking
		#the two grids should have the same centers and radius
		
#		# reorder the values to the internal order of oldGrid
#		intValues<-values
#		if(class(oldGrid)=="trigrid"){
#			intValues[oldGrid@skeleton$uiF]<-values
#		}
#		
#		if(class(oldGrid)=="hexagrid"){
#			intValues<-values[oldGrid@skeleton$aF]
#		}
		
		# easier case:
		# given that downsampling occurs
		if(downscale){
			if(method=="average"){
				#oldgrids facecenters
				loc<-locate(y, oldGrid@faceCenters)
				newVals<-tapply(INDEX=loc, X=values, function(x){mean(x, na.rm=TRUE)})
				newVals<-newVals[rownames(y@faces)]
				nVNames<-names(newVals)
				newVals<-as.numeric(newVals)
				names(newVals)<-nVNames
			}
		
		}
		
		#given that upsampling occurs
		# edge breakpoint area approximation 
		if(!downscale){
			if(method=="ebaa"){
				
				# basic algorithm
				# 1. break down every single edge of the finer, new grid using SplitArc() and produce a table with xyz coordinates
				if(class(y)=="trigrid"){
					newGridF<-y@skeleton$f[y@skeleton$f[,4]==(length(y@tessellation)-1),]
					
					# c++ function will look up produce these points
					newPoints<-.Call(Cpp_icosa_ExpandEdgesByFacesTri_, y@skeleton$v, newGridF, y@center, res)
					newPoints[,4] <- newPoints[,4]+1
					# outer order
					fNames<-names(sort(y@skeleton$uiF))
					newNames <- fNames[newPoints[,4]]
				}
				if(class(y)=="hexagrid"){
					newGridF<-y@skeleton$f[y@skeleton$f[,4]==-6,]
				
					# c++ function will look up produce these points
					newPoints<-.Call(Cpp_icosa_ExpandEdgesByFacesTri_, y@skeleton$v, newGridF, y@center, res)
					newPoints[,4] <- newPoints[,4]+1
					
					# create a proper UIF table
					tempUIF <- as.numeric(t(y@skeleton$uiF))
					names(tempUIF)<-rep(rownames(y@skeleton$uiF), each=12)
					
					# get rid of NAs
					tempUIF<-tempUIF[!is.na(tempUIF)]
					
					# outer order
					fNames<-names(sort(tempUIF))
					newNames <- fNames[newPoints[,4]]
				
				}
				
				
				# 2. look up these new coordinates in the old grid, it is not a problem if they are not found (locate produces NAs)
					oldCells<-locate(oldGrid, newPoints[,1:3])
					
					# 3. get the values from the facelayer that corresponds to these face designations
					names(x@values) <- x@names
					oldVals<-x@values[oldCells]
					
					# 4. average them out with tapply()
					newVals<-tapply(INDEX=newNames, X=oldVals, mean, na.rm=TRUE)
				
					# 5. you have an estimation based on the approximate coverages
					nVNames<-names(newVals)
					newVals<-as.numeric(newVals)
					names(newVals)<-nVNames
				
			}
		}
		
		return(newVals)
	}
	
)

