#' Add an igraph object to a predefined slot in an icosahedral grid
#'
#' @name newgraph
#		
#' @rdname newgraph
#' @return A new (\code{\link{trigrid}} or \code{\link{hexagrid}}) object with the recalculated graph.
#' @param gridObj (\code{\link{trigrid}}, \code{\link{hexagrid}}) An icosahedral grid.
#' @param ... Arguments passed to the \code{\link{gridgraph}} function.
#' @examples
#' #create a grid
#' g<-trigrid(4, graph=FALSE)
#' g<-newgraph(g)
#' 
#' @exportMethod newgraph
setGeneric(
	name="newgraph",
	package="icosa",
	def=function(gridObj,...){
		standardGeneric("newgraph")
	}

)

#' @rdname newgraph
setMethod(
	"newgraph",
	signature="trigrid",
	definition=function(gridObj,...){
		gridObj@graph<-gridgraph(gridObj,...)
		return(gridObj)
	
	}
)


#' Create or instantiate an \code{\link[igraph:make_graph]{graph}} class graph from the faces of an icosahedral grid
#'
#' The function can be applied to both grids and to \code{\link{facelayer}}-class object of \code{logical} values. The resulting graph will have the characteristics of the original grid (directed/undirected etc.). 
#' @name gridgraph
#' @return The function returns an undirected igraph graph.
#' @param x (\code{\link{trigrid}}, \code{\link{hexagrid}} or \code{\link{facelayer}}) The icosahedral grid or \code{\link{facelayer}}.
#' @param ... Arguments passed to the class specific methods.
#' @rdname gridgraph
#' @exportMethod gridgraph
setGeneric(
		name="gridgraph",
		def=function(x,...){
			standardGeneric("gridgraph")
		}
	
	)
	
#' Create or instantiate an 'igraph' class graph from the faces of a triangular grid
#'
#' @param directed \code{logical} Defaults to \code{FALSE}, creating an undirected graph. If \code{TRUE}, then the graph will be directed.
#' @param distances \code{logical} Defaults to \code{FALSE}. If \code{TRUE}, then the distances between the linked faces will be calculated and will be rendered to the edges as \code{"dist"}.
#' @rdname gridgraph
setMethod(
	f="gridgraph",
	signature="trigrid",
	definition= function(x, directed=FALSE,distances=FALSE){
		# if a graph object already exists in the grid
		if(!suppressWarnings(is.na(x@graph))[1]){
			gridGraph<-x@graph
		}else{
				
			# calculate the outer ordering
				# same format
				boolActFace <- x@skeleton$f[,4]==max(x@skeleton$f[,4])
				replaceSource<-x@skeleton$aF[boolActFace]
			
				# new order
				nOutOldIndex <- x@skeleton$n[x@skeleton$uiF,]
				nOutOldIndex <- nOutOldIndex +1
			
			
			nOut<-x@skeleton$n
			nOut[,1]<-replaceSource[nOutOldIndex[,1]]
			nOut[,2]<-replaceSource[nOutOldIndex[,2]]
			nOut[,3]<-replaceSource[nOutOldIndex[,3]]
			nOut[,4]<-replaceSource[nOutOldIndex[,4]]
			
			# c++ function to create an edgelist 
			edgeList <- .Call(Cpp_icosa_edgeListFromNeighbours_, nOut)
			
			# supress scientific notation!
			options(scipen=999)
			edgeListChar <- matrix(paste("F", edgeList, sep=""), ncol=2)
			options(scipen=0)
			
			# the arguments for the igraph function
		
			# get rid of the double edges
			edgeListChar <- unique(edgeListChar)
			
			
			# order the edges so they are not that messy
			edgeListChar<-edgeListChar[order(edgeListChar[,1]), ]
			
			
			# make a graph from that 
				graphArgs <- c(list(d=edgeListChar), list(directed=FALSE))
				
				gridGraph <- do.call(igraph::graph_from_data_frame, graphArgs)
		}	
		
		# depending on whether the directed argument was specified or not
		if(!is.null(directed)){
			if(directed){
				gridGraph<-igraph::as.directed(gridGraph)
				
			}
		
		}
		
		
		if(distances){
			edgeListChar<-igraph::get.edgelist(gridGraph)
			p0 <- x@faceCenters[edgeListChar[,1],]
			p1 <- x@faceCenters[edgeListChar[,2],]
			weights<- .Call(Cpp_icosa_ArcDistMany_, p0, p1, x@center, x@r)
			igraph::E(gridGraph)$dist <- weights
		}
		
		
		# subset, if necessary
			graph <- igraph::induced_subgraph(gridGraph, v=rownames(x@faces))
		
		return(graph)
})



#' Create or instantiate an 'igraph' class graph from the faces of a penta-hexagonal grid
#' @rdname gridgraph
setMethod(
	f="gridgraph",
	signature="hexagrid",
	definition= function(x, directed=FALSE,distances=FALSE){
	
	# get the edges of the original trigrid
		edgeListChar<-gsub("P", "F",x@skeleton$edgeTri)
		rownames(edgeListChar)<-NULL
	
		# depending on whether the directed argument was specified or not
		if(!is.null(directed)){
			if(directed==TRUE){
				# get rid of the double edges
				edgeList2<-cbind(edgeListChar[,2],edgeListChar[,1])
				edgeListChar<-rbind(edgeListChar,edgeList2)
			}
			
		
		}
		
		# order the edges so they are not that messy
		edgeListChar<-edgeListChar[order(edgeListChar[,1]), ]
		
		# make a graph from that 
			graphArgs <- c(list(d=edgeListChar), list(directed=directed))
			
			gridGraph <- do.call(igraph::graph_from_data_frame, graphArgs)
		
		
			if(distances){
				edgeListChar<-igraph::get.edgelist(gridGraph)
				p0 <- x@faceCenters[edgeListChar[,1],]
				p1 <- x@faceCenters[edgeListChar[,2],]
				weights<- .Call(Cpp_icosa_ArcDistMany_, p0, p1, x@center, x@r)
				igraph::E(gridGraph)$dist <- weights
			}
		
		
		# subset, if necessary
			graph <- igraph::induced_subgraph(gridGraph, v=rownames(x@faces))
		
		return(graph)

})



#' The neighbouring faces of faces in an icosahedral grid
#' 
#' This function will return neighbouring faces of the input faces. 
#' @name vicinity
#' 
#' @param gridObj (\code{\link{trigrid}} or \code{\link{hexagrid}}) Icosahedral grid object. 
#' 
#' @param faces (\code{character}) A vector specifying names of faces. 
#'
#' @param order (\code{numeric}) Passed to the \code{\link[igraph]{ego}} function, an integer value specifying the size of the neighborhood around a face.
#'	
#' @param output (\code{character}) The type of the output. The default \code{"vector"} 
#' 	will give back the names of the faces that adjacent to the faces specified, 
#' 	including themselves. \code{"list"} will return a list.
#'
#' @param self (\code{logical}) Flag indicating whether the input faces should be in the output. For the \code{"list"} output option, the input face names will be
#' omitted only from those character vectors that contain face names that are related to the face in question.
#'
#' @param namedorder (\code{logical}) Should the orders of the neighbouring cells be reported (\code{TRUE}) or just the names of the cells (default, \code{FALSE}).
#' @param ... Arguments passed to the \code{\link[igraph]{ego}} function.
#' @examples
#' g <- trigrid(3)
#' ne <- vicinity(g, c("F4", "F10"))
#' ne
#' 
#' @return A \code{character} vector or a \code{list} of \code{character} vectors.
#' 	
#' @exportMethod vicinity
#' @rdname vicinity
setGeneric(
	name="vicinity",
	def=function(gridObj,faces,...){
		standardGeneric("vicinity")
	}
)

#' @rdname vicinity
setMethod(
	"vicinity",
	signature=c("trigrid","character"),
	definition=function(gridObj,faces, order=1, output="vector",self=TRUE,namedorder=FALSE,  ...){
		#if no @graph found
		if(suppressWarnings(is.na(gridObj@graph)[1])){
			stop("Slot @graph is empty. Use newgraph() to add an igraph respresentation. ")
		}
		
		# get rid of NA entries
		faces<-faces[!is.na(faces)]
		
		# check whether the facein put is actually in the grid
		if(sum(!faces%in%rownames(gridObj@faces))!=0)
			stop("Invalid face name.")
		
		# get the actual vertices
		nList<-igraph::ego(gridObj@graph, order=order, faces,...)
		
		# get only the names of the faces
		newList<-lapply(nList, function(x){x$name})
		
		# should the orders of the neighbours be reported? if yes..
		if(namedorder){
			newStructure <- lapply(newList, function(x){
				temp <- rep(order, length(x))
				names(temp) <- x
				temp
			})

			# rerun the process for all the lower orders iterativels
			for(i in order:1){
				nList<-igraph::ego(gridObj@graph, order=i, faces,...)
		
				# get only the names of the faces
				nL<-lapply(nList, function(x){x$name})

				newStructure <- mapply(FUN=function(x,y){
					y[x] <- i
					return(y)
				},nL, newStructure)

			}

			newList <- lapply(newStructure, function(x){x[1] <- 0; return(x)})

		}


		if(!self & !namedorder){
			newList<-lapply(newList, function(x){x[-1]})
		}

		
		if(output=="list"){
			return(newList)
		}
		
		if(output=="vector"){
			temp <- unlist(newList)
			if(self){
				return(sort(unique(c(temp, faces))))
			}else{
				return(sort(temp[!temp%in%faces]))
			
			}
		}
	}
)
