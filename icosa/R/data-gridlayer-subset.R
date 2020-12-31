#' Subsetting a gridlayer object.
#'
#' The function extracts subsets of the \code{gridlayer} depending on different criteria.
#'
#' The following methods are incorporated into the function: If \code{i} argument is a vector of integers, they will be interpreted as indices. If the \code{numeric} \code{i} contains either the lamin, lamax, lomin or lomax names, the subsetting will be done using the latitude-longitude coordinates outlined by these 4 values. Logical subsetting and subsetting by face names are also possible.
#'
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}) The object to be subsetted.
#' @param i (\code{logical}, \code{numeric} or \code{character}) The subscript vector, specifying the faces that are used for subsetting. As in \code{\link[base]{subset}}.

#'
#' @rdname subset
#' @exportMethod subset
setMethod(
	"subset",
	signature="gridlayer",
	definition=function(x, i){
		if(is.numeric(i)){
			#add checking for lat/long subsetting
			# lat-long mode of subsetting
			potConds<-c("lamin", "lamax", "lomin", "lomax")
			if(sum(names(i)%in%potConds)>0){
				#if it contains an unitelligible names
				if(sum(!names(i)%in%potConds)>0) 
					warning("Some subscript condition names were not recognized.")
				
				
				#in case you want something at the dateline
				normal <- T
				if(sum(c("lomax", "lomin")%in%names(i))==2){
					if(i["lomin"]>i["lomax"]){
						normal<- F
					}
				}
				
				#get the facecenters
				actGrid<-get(x@grid)
				pol <- CarToPol(actGrid@faceCenters, norad=TRUE, origin=actGrid@center)
				
				boolSelect<-rep(T, nrow(pol))
				
				#longitude
				if(normal){
					#minimum longitude condition
					if("lomin"%in%names(i)){
						boolSelect <- boolSelect & pol[,1]>=i["lomin"]
					}
					
					#maximum longitude condition
					if("lomax"%in%names(i)){
						boolSelect <- boolSelect & pol[,1]<=i["lomax"]
					}
				}else{
					#minimum longitude condition
					if("lomin"%in%names(i)){
						boolSelect <- boolSelect & pol[,1]>=i["lomin"]
					}
					
					#maximum longitude condition
					if("lomax"%in%names(i)){
						boolSelect <- boolSelect | pol[,1]<=i["lomax"]
					}
				
				}
				
				#minimum latitude condition
				if("lamin"%in%names(i)){
					boolSelect <- boolSelect & pol[,2]>=i["lamin"]
				}
				
				#minimum latitude condition
				if("lamax"%in%names(i)){
					boolSelect <- boolSelect & pol[,2]<=i["lamax"]
				}
				
				i<-rownames(actGrid@faceCenters)[boolSelect]
				# control will pass over to the subsetting by facenames

			}else{
			
			# index subsetting
				y<-x
				y@names<-y@names[i]
				y@values<-y@values[i]
				y@length<-length(y@values)
			}
		}
		if(is.logical(i)){
			if(length(i)==(length(x@names))){
				i<-x@names[i]
			}else{
				stop("Length of logical subscript does not match the facelayer.")
			}
		
		}
		
		if(is.character(i)){
			if(sum(i%in%x@names)==length(i)){
				y<-x
				y@names<-i
				y@values<-y@values[x@names%in%i]
				y@length<-length(y@values)
			}
		
		}
		
		return(y)
		
	}
)	


#subsetting for layers
#' Extraction from a gridlayer using indices
#' 
#' Shorthand to the \code{\link[icosa]{subset}} function.
#' 
#' @param x (\code{\link{facelayer}}) The object to be subsetted.
#' @param i (\code{logical}, \code{numeric} nor \code{extent}) The subscript vector, or extent, specifying the faces that are used for subsetting. As in \code{\link[base]{subset}}.
#' @exportMethod "["
#' @return The extraction methods return \code{\link{facelayer}}-class objects.
#' @rdname extract-methods
setMethod(
	"[",
	signature=c("gridlayer","ANY", "missing"),
	definition=function(x,i){
		subset(x, i)
	
	}
)
 
#' @exportMethod "["
#' @rdname extract-methods
setMethod(
	"[",
	signature=c("gridlayer","Extent", "missing"),
	definition=function(x,i){
		#check the extent object
		
		actGrid <- get(x@grid)
		pol <- CarToPol(actGrid@faceCenters, origin=actGrid@center)
		
		boolLong<-pol[,1]>=i@xmin & pol[,1]<=i@xmax
		boolLat<-pol[,2]>=i@ymin & pol[,2]<=i@ymax
		
		nm<-rownames(pol)[boolLong & boolLat]
	
	
		subset(x, nm)
	
	}
)

# this method produces a warning without the aliases!!!


#' Replacement of elements in a gridlayer object.
#' 
#' Function to replace specific elements in a gridlayer object 
#' 
#' All these methods are implementing direct replacement in the \code{@values} slot of a layer, depending on criteria used for subsetting. 
#'
#' @param value The replacement values.
#'
#' @docType methods
#' @aliases [<-,gridlayer-method
#' @exportMethod "[<-"
#' @rdname extract-methods
setReplaceMethod(
	"[",
	signature=c("gridlayer"),
#	definition=function(x,i,j,..., value){
	definition=function(x,i,value){
		y<-x
		#named vector replacement
		if(length(names(value))>0 & missing(i)){
			if(sum(names(value)%in%y@names)==length(value)){
				u<-y@values
				names(u)<-y@names
				u[names(value)]<-value
				y@values<-u
			}
		}else{
		#numeric
			
			if(is.numeric(i)){
				#add checking for lat/long subsetting
				# lat-long mode of subsetting
				potConds<-c("lamin", "lamax", "lomin", "lomax")
				if(sum(names(i)%in%potConds)>0){
					#if it contains an unitelligible names
					if(sum(!names(i)%in%potConds)>0) 
						warning("Some subscript condition names were not recognized.")
					
					
					#in case you want something at the dateline
					normal <- T
					if(sum(c("lomax", "lomin")%in%names(i))==2){
						if(i["lomin"]>i["lomax"]){
							normal<- F
						}
					}
					
					#get the facecenters
					actGrid<-get(x@grid)
					pol <- CarToPol(actGrid@faceCenters, norad=TRUE, origin=actGrid@center)
					
					boolSelect<-rep(T, nrow(pol))
					
					#longitude
					if(normal){
						#minimum longitude condition
						if("lomin"%in%names(i)){
							boolSelect <- boolSelect & pol[,1]>=i["lomin"]
						}
						
						#maximum longitude condition
						if("lomax"%in%names(i)){
							boolSelect <- boolSelect & pol[,1]<=i["lomax"]
						}
					}else{
						#minimum longitude condition
						if("lomin"%in%names(i)){
							boolSelect <- boolSelect & pol[,1]>=i["lomin"]
						}
						
						#maximum longitude condition
						if("lomax"%in%names(i)){
							boolSelect <- boolSelect | pol[,1]<=i["lomax"]
						}
					
					}
					
					#minimum latitude condition
					if("lamin"%in%names(i)){
						boolSelect <- boolSelect & pol[,2]>=i["lamin"]
					}
					
					#minimum latitude condition
					if("lamax"%in%names(i)){
						boolSelect <- boolSelect & pol[,2]<=i["lamax"]
					}
					
					i<-rownames(actGrid@faceCenters)[boolSelect]
					# control will pass over to the subsetting by facenames
	
				}else{
					# index subsetting
					actGrid<-get(x@grid)
	
					subGrid<-subset(actGrid,i)
					i<-rownames(subGrid@faces)
				}
			}
		
			# pass on from the numeric too!
			if(is.character(i)){
				if(sum(i%in%y@names)==length(i)){
					u<-y@values
					names(u)<-y@names
					u[i]<-value
					y@values<-u
				}else{
					stop("Invalid character subscript.")
				}
			}
			if(is.logical(i)){
				y@values[i]<-value
			}
		}
		
		return(y)
	
	}
)

