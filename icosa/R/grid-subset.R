#' Subsetting an icosahedral grid or data layers organized with them
#' 
#' This is a generic function used to access data from either a triangular or hexagonal grid using the names of the faces, integers or logical vectors. 
#' 
#' The function returns subsets of the grid pertaining to the specified faces that can be used for additional operations (e.g. plotting). 
#' The subscript vector can be either a logical, character or numeric one. The character vector should contain the names of faces, the logical subscript should have 
#' the same length as the number of faces in the order in which the faces are present in the \code{faces} slot. 
#' The numeric vector can either refer to indices to the rownames of faces in the faces slot, or
#' to surfaces bounded by longitude/latitude data. In the latter case, the the vector should contain an element with a names of at least one of the \code{"lomax"}, \code{"lamax"}, 
#' \code{"lomin"} or \code{"lamin"} strings (lo for longitude, la: latitude, min: minimum, max: maximum). In case a subset around the dateline is needed a larger longitude to a smaller longitude value is needed (e.g. between \code{150}° to \code{-150}°). 
#' 
#' 
#' 
#' @examples
#' #create a triangular grid
#' g <- trigrid(c(2,2))
#' 
#' #make a subset pertaining to the faces
#' subG1 <- subset(g, c("F1", "F33"))
#'     
#' #additional way of subsetting
#' subG2 <- g[1:15] # selects faces F1 through F15
#' logicalSub<-sample(c(TRUE,FALSE), nrow(g@faces), replace=TRUE)
#' subG3 <- g[logicalSub]
#' 
#' #plot the subset in 3d space
#' plot3d(subG3)
#' 
#' # previously mentioned case around the dateline
#' gDateLine<-g[c(lomax=-150, lomin=150)]
#' plot3d(gDateLine)
#' @rdname subset
#' @return Subset of the input grid. The class of the original object is retained, the \code{@skeleton} slot contains all previous information.
#' @exportMethod subset
"subset"


#' Subset method of the trigrid class
#' @rdname subset
setMethod(
	"subset",
	signature="trigrid",
	definition=function(x, i){
	
	#	i<-c("F295", "F300")
	#	x<-grid
		
		#checking
		if(is.numeric(i)){
			#add checking for lat/long subsetting
			# lat-long mode of subsetting
			potConds<-c("lamin", "lamax", "lomin", "lomax")
			if(sum(names(i)%in%potConds)>0){
				#if it contains an unitelligible names
				if(sum(!names(i)%in%potConds)>0) 
					warning("Some subcript condition names were not recognized.")
				
				
				#in case you want something at the dateline
				normal <- T
				if(sum(c("lomax", "lomin")%in%names(i))==2){
					if(i["lomin"]>i["lomax"]){
						normal<- F
					}
				}
				
				#get the facecenters
				pol <- CarToPol(x@faceCenters, norad=TRUE, origin=x@center)
				
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
				
				i<-rownames(x@faceCenters)[boolSelect]
				# control will pass over to the subsetting by facenames

			}else{
			
			# index subsetting
				i<-paste("F", i, sep="")
			}
		}		
		
		if(is.logical(i)){
			if(length(i)==nrow(x@faces)){
				i<-rownames(x@faces)[i]
			}else{
				stop("Logical subscript has wrong length.")
			}
		}
		
		#check whether there are NAs in the i variables
		if(sum(is.na(i))){
			warning("The specified face names contain NA values.")
			#omit
			i<-i[!is.na(i)]
		}
		
		if(sum(!i%in%rownames(x@faces))>0) stop("Invalid face names entered.")
		
	
		#triGrid subset
		#1. subset of faces
			subsetFaces<-x@faces[i, ,drop=FALSE]
			
			# skeleton - faces - use skeleton$uiF to deactivate faces in skeleton$aF
			#false vector
			tempFacesLog<-rep(F, length(x@skeleton$aF))
			#what should not be removed
			tempFacesLog[x@skeleton$offsetF+x@skeleton$uiF[which(names(x@skeleton$uiF)%in%i)]] <- TRUE
			#remove everthing but that
			x@skeleton$aF[!tempFacesLog]<- FALSE
		
		#2. points
			pointnames<-unique(as.character(subsetFaces))
			subsetVertices<-x@vertices[pointnames,, drop=FALSE]
			
			# skeleton - use skeleton$uiV to deactivate faces in skeleton$aV
			tempVerticesLog<-rep(F,length(x@skeleton$aV))
			tempVerticesLog[x@skeleton$uiV[which(names(x@skeleton$uiV)%in%pointnames)]] <- TRUE
			
			x@skeleton$aV[!tempVerticesLog] <- FALSE
		
		#3. subset of edges
			#the original logical of the edges
			edgeTemp<-x@skeleton$aE
			
			#logical for the points - kept or not?
			logE<-matrix(as.logical(x@skeleton$aV)[as.numeric(x@skeleton$e+1)], ncol=2)
			
			#where both points are needed in the subset, keep both!
			aE<-apply(logE,1, sum)==2
			x@skeleton$aE<-aE
			
			#use that but only where it was subsetted previous - for UI
			subsetEdges<-x@edges[aE[edgeTemp],, drop=FALSE]
			
			
		
		#4. subset of faceCenters
			subsetFaceCenters<-x@faceCenters[i,, drop=FALSE]
			
		
		#6. copy over the original object
		#	y<-x
			
			#and change it accordingly
			x@vertices=subsetVertices
			x@faces=subsetFaces
			x@edges=subsetEdges
			x@faceCenters=subsetFaceCenters
			
			x@length<- c(
			"vertices"=nrow(x@vertices),
			"edges"=nrow(x@edges),
			"faces"=nrow(x@faces))
		
		#7. if the @sp slot contains data, subset it too
		if(suppressWarnings(!is.na(x@sp))){
			x@sp <- x@sp[rownames(x@faces)]
		}
		
		if(suppressWarnings(!is.na(x@graph))[1]){
			x@graph <- igraph::induced_subgraph(x@graph,rownames(x@faces))
		}
		
		
		return(x)
		
	}
)

#' Subset method of the hexagrid class
#' @rdname subset
setMethod(
	"subset",
	signature="hexagrid",
	definition=function(x, i){
	
		#checking
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
				pol <- CarToPol(x@faceCenters, norad=TRUE, origin=x@center)
				
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
				
				i<-rownames(x@faceCenters)[boolSelect]
				# control will pass over to the subsetting by facenames

			}else{
			
			# index subsetting
				i<-paste("F", i, sep="")
			}
		}		
		
		if(is.logical(i)){
			if(length(i)==nrow(x@faces)){
				i<-rownames(x@faces)[i]
			}else{
				stop("Logical subscript has wrong length.")
			}
		}
		
	#	i<-c("F295", "F300")
	#	x<-grid
	
		#check whether therer are NAs in the i variables
		if(sum(is.na(i))){
			warning("The specified face names contain NA values")
			#omit
			i<-i[!is.na(i)]
		}
		
		
		if(sum(!i%in%rownames(x@faces))>0) stop("Invalid face names entered.")
		
		# order the i vector to avoid any potential errors
			i<-sort(i)
		
		
		# hexaGrid subset
		#1. subset of faces
			#ui representation
			subsetFaces<-x@faces[i, ,drop=FALSE]
			
			# indexing of $f
			subFaceIndex<-as.numeric(x@skeleton$uiF[i,])
			aSF<-x@skeleton$aSF
			aSF[!1:length(x@skeleton$aSF)%in%subFaceIndex]<-0
		
			aF<-x@skeleton$aF
			aF[!paste("F",aF,sep="")%in%i]<-0
			
		#2. points 
		#2a.(vertices)
			pointnames<-unique(as.character(subsetFaces))
			subsetVertices<-x@vertices[rownames(x@vertices)%in%pointnames,]
		
			#indexing in $v
			aV<-x@skeleton$aV
			aV[!1:length(aV)%in%x@skeleton$uiV[rownames(subsetVertices)]]<-0
		
		
		#3b. subset of edges
			edgeTemp<-x@skeleton$aE
			
			logE<-matrix(as.logical(aV)[as.numeric(x@skeleton$e+1)], ncol=2)
			aE<-apply(logE,1, sum)==2
			
			subsetEdges<-x@edges[aE[edgeTemp],]
		
		#4. subset of faceCenters
			subsetFaceCenters<-x@faceCenters[i,, drop=FALSE]
			
		#5. copy over the original object
			y<-x
			
			#and change it accordingly
			y@vertices=subsetVertices
			y@faces=subsetFaces
			y@edges=subsetEdges
			y@faceCenters=subsetFaceCenters
			
			y@skeleton$aV<-aV
			y@skeleton$aSF<-aSF
			y@skeleton$aE<-aE
			y@skeleton$aF<-aF
			
			y@length<- c(
			"vertices"=nrow(y@vertices),
			"edges"=nrow(y@edges),
			"faces"=nrow(y@faces))
			#7. if the @sp slot contains data, subset it too
		
		if(suppressWarnings(!is.na(x@sp))){
			y@sp <- y@sp[rownames(y@faces)]
		}
		
		if(suppressWarnings(!is.na(x@graph))[1]){
			y@graph <- igraph::induced_subgraph(y@graph,rownames(y@faces))
		}
			return(y)
		
	}
)


#' Extract subset faces of a trigrid or hexagrid object using.
#' @rdname subset
setMethod(
	"[",
	signature="trigrid",

	definition=function(x,i){
		subset(x, i)
	
	}
)

