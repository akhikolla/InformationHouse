
#gridlayer basic class
	#class definition
	gridlayer <- setClass(
		"gridlayer",
		
		slots= c(
			grid= "character",
			tessellation ="numeric",
			gridclass = "character",
			names = "character",
			values= "vector",
			length= "integer"
		)
	
	)

setMethod("show", signature ="gridlayer",
		definition = function(object){
		#	cat(paste(class(object), "of", object@grid ,"with", object@length, class(object@values), "values\n", sep=" "))
		#	cat(object@values, fill=TRUE)
			
			actGrid<-get(object@grid)
			
			
			cat(paste("class        : ", class(object),"\n", sep=""))
			cat(paste("linked grid  : \'", object@grid,"\' (name), ", class(actGrid)," (class), ",
				paste(as.character(object@tessellation), collapse=",")
				, " (tessellation)", "\n", sep=""))
			cat(paste("dimensions   : ", length(object), " (values)", " @ ", "mean edge length: ",round(actGrid@edgeLength[1],2), " km, " ,round(actGrid@edgeLength[2],2), " degrees",  "\n" ,sep=""))
			if(sum(is.na(object@values))==length(object)){
				mx<-NA
				mn<-NA
			}else{
				mx<-max(object@values, na.rm=TRUE)
				mn<-min(object@values, na.rm=TRUE)
			}
			
			cat(paste("values       : ", class(object@values)), "\n")
			cat(paste("max value    : ", mx), "\n")
			cat(paste("min value    : ", mn), "\n")
			cat(paste("missing      : ", sum(is.na(object@values))), "\n")
			
		} 
	)
