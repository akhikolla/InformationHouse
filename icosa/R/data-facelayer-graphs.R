
#' Create or instantiate an 'igraph' class graph based on a facelayer
#'
#' @rdname gridgraph
setMethod(
	f="gridgraph",
	signature="facelayer",
	definition=function(x){
		# the grid object of the facelayer
		actGrid <- get(x@grid)
		
		checkLinkedGrid(actGrid, x)
		
		# only the occupied cells should be part of the grid
		# or: where the logical values indicate absence
		if(is.logical(x@values)){
			x@values[x@values==FALSE] <- NA
		}
		
		# these will be the nodes of the graph
		nodes<- x@names[!is.na(x@values)]
		
		# subset the graph to the nodes which are indicated in the facelayer
		actGraph <- igraph::induced_subgraph(actGrid@graph, nodes)
		
		# the values of the facelayer should be also present
		if(!is.logical(x@values)){
			igraph::V(actGraph)$value <- x@values[!is.na(x@values)]
		}
		
		return(actGraph)

	}
)

