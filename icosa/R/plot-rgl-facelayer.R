#' 3d plotting of a facelayer of an icosahedral grid or its subset
#'
#' The function is built on the openGL renderer of the R package \code{rgl}. The default plotting window size is \code{800x800} pixels. In case you want to override this, please
#' use the function with \code{defaultPar3d=FALSE} after running \code{\link[rgl]{par3d}}\code{(windowRect=<>)}. 
#'  
#' @param defaultPar3d (\code{logical}) Flag indicating whether the default settings for \code{\link[rgl]{par3d}} are to be used \code{(windowRect = c(50, 60, 800, 800), zoom=0.8)}.
#' @param frame (\code{logical}) If set to TRUE the grid line structure will be plotted.
#' 
#' @exportS3Method rgl::plot3d facelayer
#' @exportS3Method plot3d facelayer
#' @rdname plot3d
#' @export plot3d.facelayer
plot3d.facelayer <- function(x,type="f",frame=TRUE, guides=TRUE, defaultPar3d=TRUE, ...){
	
		# default par3d options
		if(defaultPar3d){
			rgl::par3d(windowRect = c(50, 60, 800, 800), zoom=0.8)
		}
			
		actGrid  <- get(x@grid)
		checkLinkedGrid(actGrid, x)
			
		#do not allow arguments to pass through!
		if(frame==TRUE){
			plot3d(actGrid, guides=guides, col="gray50")
		}else{
			plot3d(actGrid, type="n", guides=guides)
			
			#boundaries
			if(type=="l"){
				lines3d(x,...)
		
			}
		}
		if(type=="f"){
			faces3d(x,specular="black",...)
		}
		
		#add additional types of plotting to this method
		#no plotting
		if(type=="n"){
		}
		
	
	}

#' Methods of 3d face plotting.
#' 
#' @param col (\code{character}) Graphical parameter indicating the colours of the faces. A single value is accepted for \code{logical} values. Multiple colors will be passed to \code{\link[grDevices:colorRamp]{colorRampPalette}}, to create palettes for heat maps in case of \code{numeric} values. The default plotting method in this case is the reversed \code{\link[grDevices:palettes]{heat.colors}}. In case of categorical data, random colors will be chosen.
#' @param breaks (\code{numeric}) Vector stating the breakpoints between the plotted levels. The argument is passed to the \code{\link[base]{cut}} function. 
#' @param inclusive (\code{logical}): If there are values beyond the limits of breaks, should these be represented in the plot (\code{TRUE}) or left out completely \code{FALSE}?
#' @param legend (\code{logical}) Should the heatmap legend be plotted?
#'
#' @exportMethod faces3d
#' @rdname faces3d
#' @examples
#' h <- hexagrid(8)
#' b <- facelayer(h)
#' values(b)<- rnorm(length(b))
setMethod(	
	f="faces3d",
	signature="facelayer",
	definition= function(x,col="heat",breaks=NULL, inclusive=TRUE, legend=TRUE,  ...){
		# extract the grid that needs to be plotted:
		actGrid  <- get(x@grid)
	#	checkLinkedGrid(actGrid, x)
		discrete<- FALSE
		#check whether the  grid is actually updated
		if(sum(x@names%in%rownames(actGrid@faces))!=length(x)) 
		stop("The facenames in thelinked grid does not match the facelayer object.")
		
		#when the valuues are logical
		#FALSEs do not plot; NAs do not plot, TRUEs plot
		
		# defend 'breaks'
		if(!is.null(breaks)){
			if(!is.numeric(breaks)) stop("The 'breaks' argument has to be numeric.")
			if(length(breaks)<3) stop("You need to provide at least three values to the 'breaks' argument.")
		}


		# if the grid is numerical and it has only one value, make it logical
		if(class(x@values)%in%c("integer","double", "numeric")){
			if(length(unique(x@values[!is.na(x@values)]))==1){
				x@values<-as.logical(x@values)
			}
			
		}
		if(is.logical(x@values)){
			#just add NAs where the values are 0
			x@values[x@values==FALSE]<-NA
		}
		
		#if the number of values does not match the grid face no
		boolPresent1<-rep(T,nrow(actGrid@faces))
		if(length(x)!=nrow(actGrid@faces)){
			boolPresent1<-rownames(actGrid@faces)%in%x@names
			actGrid<-subset(actGrid, rownames(actGrid@faces)[boolPresent1])
		}
		
		# in case there are NAs, do a subsetting before going on
		# rgl does not understand col=NA as omission of plotting
		if(sum(is.na(x@values))>0){
			# select only the faces that are available
			boolPresent<-!is.na(x@values)
			#1. the values
			x@values<-x@values[boolPresent]
			#2. the names too
			x@names<- x@names[boolPresent]
			#3. number
			x@length <- sum(boolPresent)
			
		#	#do a pseudo subsetting!
		#	tempFacesLog<-rep(F, length(actGrid@skeleton$aF))
		#	#what should not be removed
		#	tempFacesLog[actGrid@skeleton$offsetF+actGrid@skeleton$uiF[which(names(actGrid@skeleton$uiF)%in%x@names)]] <- TRUE
		#	#remove everthing but that
		#	actGrid@skeleton$aF[!tempFacesLog]<- FALSE
		
		
			actGrid<-subset(actGrid, x@names) # the real subsetting
			
		}
		#when the values are logical
		if(class(x@values)=="logical"){
			#set default color value
			faces3d(actGrid,col=col,...)
		}
		
		# when  numerical values are added to the facelayer object, do a heatmap!
		if(class(x@values)%in%c("integer","double", "numeric")){
			
			
			# calculate the breaking vector
			if(is.null(breaks)){
				minimum <- min(x@values)
				maximum <- max(x@values)
				steps <- length(x)+1
				
				# the vector used to cut the plottted variable
				useBreaks <- seq(minimum, maximum,length.out=steps)
			}else{
				minimum <- min(breaks)
				maximum <- max(breaks)
				useBreaks <- breaks
			}

			# still need to include limitations
			bMax <- FALSE
			bMin <- FALSE
			if(inclusive){
				# values that are beyond the minimum boundary set by breaks
				beyondMax <- which(x@values>maximum)
				if(length(beyondMax)>1){
					x@values[beyondMax] <- maximum
					bMax <- TRUE
				}
				# values that are beyond the minimum boundary set by breaks
				beyondMin <- which(x@values<minimum)
				if(length(beyondMin)>1){
					x@values[beyondMin] <- minimum
					bMin <- TRUE
				}
			}


			#do a heatmap!
			#create a ramp, with a given number of colours
			#the color vector will control the heatmap
			if(length(col)==1){
				# predefined
				if(col=="heat"){
#					col<-c("red","orange","yellow", "white")
					cols <- rev(grDevices::heat.colors(length(useBreaks)-1))
				
					cols<-substring(cols, 1,7)
				}else{
					
					if(length(col)==1){
						stop("You specified only one color.")
					}
				
				}
			} else{
			#do a heatmap!
				ramp<-grDevices::colorRampPalette(col, bias=2, space="Lab")
				# produce as many colours as there are values
				cols <- ramp(length(useBreaks)-1)
			}

			# do the cutting
			alreadyCut <- base::cut(x@values, breaks=useBreaks, include.lowest=TRUE)

			# transfer the factor to indices
			trans2 <- as.numeric(alreadyCut)

			# this is the ui sequence	
			faceColors<-cols[trans2]
			
			if(class(actGrid)=="trigrid"){
				
				#in the inner sequence
				#create a source vector as if it was complete
					faceColors2<-rep(NA, length(actGrid@skeleton$uiF))
					names(faceColors2)<-paste("F", 1:length(faceColors2), sep="")
					faceColors2[names(x)]<-faceColors
				
				#order them
					faceColors3<-rep(NA, length(faceColors2))
					faceColors3[actGrid@skeleton$uiF]<-faceColors2
				
				#and get rid of the NAs
				faceColors3<-faceColors3[!is.na(faceColors3)]
			
				
			}
			if(class(actGrid)=="hexagrid"){
				
				tu <- as.numeric(t(actGrid@skeleton$uiF[names(x),]))
				
				empty<-rep(NA, nrow(actGrid@skeleton$f))
				
				fc<-rep(faceColors, each=12)
				
				empty[tu[!is.na(tu)]] <-fc[!is.na(tu)]
				
				noNA<-empty[as.logical(actGrid@skeleton$aSF)]
			
			###	
				# get the subfaces where there is information
			#	f<-as.data.frame(actGrid@skeleton$f[as.logical(actGrid@skeleton$aSF),1:3])
				
				#which outer faces do the subfaces belong?
				aas<- actGrid@skeleton$aSF[as.logical(actGrid@skeleton$aSF)]
				
				#create a vector fro all the total colors (as if the grid was full)
				totCol<-rep(NA,nrow(actGrid@skeleton$uiF))
				names(totCol) <- paste("F", 1:length(totCol), sep="")
				
				# insert the information
				totCol[names(x)]<-faceColors
				
				#reorder the colors to the subfaces
				faceColors3<-totCol[aas]
				
				
			#	
			#	temp<-cbind(f, newCol, stringsAsFactors=FALSE)
			#
			#	temp2<-temp[order(temp[,1]),]
			#	temp2<-unique(temp2)
				
				
				
			}
			faces3d(actGrid,col=rep(faceColors3, each=3),...)
			
			
		# numeric heatmap!			
			# increase the resolution when you plot the legend
			currentset<-rgl::par3d("windowRect")
			currentset2<-currentset
			currentset2[3]<-currentset[3]*1.5
			currentset2[4]<-currentset[4]*1.5
			
			
			# what should be passed to the heatmaplegend
			if(!discrete){
				tickLabs <- useBreaks
			}else{
				tickLabs <-  (useBreaks+useBreaks[2:(length(useBreaks)+1)])/2
				tickLabs <- tickLabs[!is.na(tickLabs)]
			}
				

			# double the resolution
			rgl::par3d(windowRect=currentset2)
			# plot the background
			rgl::bgplot3d(
				# turn off the graphical parameters warning bullshit
				suppressWarnings(
					if(legend) heatMapLegend(cols,vals=tickLabs,...)
				)
			)
				
			rgl::par3d(windowRect=currentset)

	
		}
		
		#when all the values are colors
		#plot faces as 
		if(class(x@values)=="character" & sum(x@values%in%grDevices::colors())==x@length){
			faces3d(actGrid, col=x@values, plot="faces",...)
			
		}
		
		# when the values are text | they are not colors
		if(class(x@values)=="character" & !sum(x@values%in%grDevices::colors())==x@length){
			# state the labels in 3d on the face (using the centers of the faces)
			colorAll <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = TRUE)]
			active<-factor(x@values)
			if(length(levels(active))>length(colorAll)){
				cols<-sample(colorAll, length(levels(active)), replace=TRUE)
			}else{
				cols<-sample(colorAll, length(levels(active)), replace=FALSE)
			}
			
			faceColors<-cols[as.numeric(active)]
			
			if(class(actGrid)=="trigrid"){
				
				#in the inner sequence
				#create a source vector as if it was complete
					faceColors2<-rep(NA, length(actGrid@skeleton$uiF))
					names(faceColors2)<-paste("F", 1:length(faceColors2), sep="")
					faceColors2[names(x)]<-faceColors
				
				#order them
					faceColors3<-rep(NA, length(faceColors2))
					faceColors3[actGrid@skeleton$uiF]<-faceColors2
				
				#and get rid of the NAs
				faceColors3<-faceColors3[!is.na(faceColors3)]
			
				
			}
			if(class(actGrid)=="hexagrid"){
				
				tu <- as.numeric(t(actGrid@skeleton$uiF[names(x),]))
				
				empty<-rep(NA, nrow(actGrid@skeleton$f))
				
				fc<-rep(faceColors, each=12)
				
				empty[tu[!is.na(tu)]] <-fc[!is.na(tu)]
				
				noNA<-empty[as.logical(actGrid@skeleton$aSF)]
			
			###	
				# get the subfaces where there is information
			#	f<-as.data.frame(actGrid@skeleton$f[as.logical(actGrid@skeleton$aSF),1:3])
				
				#which outer faces do the subfaces belong?
				aas<- actGrid@skeleton$aSF[as.logical(actGrid@skeleton$aSF)]
				
				#create a vector fro all the total colors (as if the grid was full)
				totCol<-rep(NA,nrow(actGrid@skeleton$uiF))
				names(totCol) <- paste("F", 1:length(totCol), sep="")
				
				# insert the information
				totCol[names(x)]<-faceColors
				
				#reorder the colors to the subfaces
				faceColors3<-totCol[aas]
				
				
			#	
			#	temp<-cbind(f, newCol, stringsAsFactors=FALSE)
			#
			#	temp2<-temp[order(temp[,1]),]
			#	temp2<-unique(temp2)
				
				
				
			}
			faces3d(actGrid,col=rep(faceColors3, each=3),...)
			
			
		}
		
		# when the values are factors!
		if(class(x@values)=="factor"){
			# depending on the number of levels, more color palettes might be useful
			if(length(levels(factor(x@values))) <= 7){
				faces3d(actGrid, col=rep(as.numeric(x@values),each=3), ...)
			
			}
		}	
	#	legend3d()
	}
)
