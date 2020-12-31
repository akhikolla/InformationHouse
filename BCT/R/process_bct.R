library(stringr)
library(igraph)


# =========================================================================================
#' @useDynLib BCT, .registration=TRUE
#' @importFrom grDevices dev.new graphics.off 
#' @importFrom graphics par text title

# =========================================================================================

#' @title Plot tree with given contexts
#' @description Plots a tree depicting a model with the given set of contexts.
#' @param s vector containing the contexts of the leaves of the desired tree. 
#' @return plot of the desired tree model.
#' @export
#' @seealso \code{\link{BCT}}, \code{\link{draw_models}}
#' 
#' @examples
#' # Construct an example vector:
#' r <- c("a", "ab", "aab", "b", "ba")
#' 
#' show_tree(r)
#' 
#' # If the input contains digits:
#' q <- c(11,1,0)
#' 
#' show_tree(q)
show_tree <- function(s){
  
  edge_size <- 0.1
  vertex_color <- "white"
  vertex_size <- 35
  vertex_frame_color <- "black"
  s <- s[order(unlist(s),decreasing=FALSE)]
  
  # create a graph with one empty-named node and no edges
  # used for show_tree and draw_models functions
  init_g=igraph::graph.empty(n=1)
  igraph::V(init_g)$name=""
  g = init_g
  for(word in c(s)){
    # turns "10100" into c("1","10","101","1010", 10100")
    subwords = stringr::str_sub(word, 1, 1:nchar(word))
    # make a graph long enough to hold all those sub-words + start node
    subg = igraph::graph.lattice(length(subwords)+1,directed=TRUE)
    # set vertex nodes to start node plus sub-words
    igraph::V(subg)$name=c("",subwords)
    # merge *by name* into the existing graph
    g = igraph::graph.union(g, subg)
  }
  
  
  if(Sys.getenv("RSTUDIO")!="1") # if the user uses RStudio or R
    dev.new()
  if(length(s)==1)
    igraph::plot.igraph(init_g,layout=igraph::layout.reingold.tilford, edge.arrow.size = edge_size, vertex.color=vertex_color, vertex.frame.color="white")
  else{
    igraph::plot.igraph(g,layout=igraph::layout.reingold.tilford, edge.arrow.size = edge_size, vertex.color=vertex_color, vertex.frame.color=vertex_frame_color)
  }
}

# =========================================================================================

#' @title Plot the results of the BCT and kBCT functions
#' @description This function plots the models produced by the BCT and kBCT functions.
#' @param lst output of the BCT/kBCT function.
#' @return plots of the BCT/kBCT output models.
#' @export
#' @seealso \code{\link{show_tree}}, \code{\link{BCT}}, \code{\link{kBCT}}
#' 
#' @examples
#' 
#' # Use the pewee dataset as an example:
#' q <- BCT(pewee, 5) # maximum depth of 5
#' 
#' draw_models(q)
#' r <- kBCT(pewee, 5, 3) 
#' 
#' # maximum depth of 5, and k = 3 (top 3 a posteriori most likely models)
#' draw_models(r)
draw_models <- function(lst){
  graphics.off()
  
  if("posterior_odds"%in% colnames(lst[[2]])){ #kbct results -> bct results does not contain the posterior_odds column
    for (i in 1:length(lst[[1]])){
      show_tree(lst[[1]][[i]])
      t <- sprintf("MAP tree %i", i)
      st <- sprintf("Max depth: %i \nNumber of leaves: %i", lst[[2]]['max_depth'][[1]][i], lst[[2]]['number_leaves'][[1]][i])
      st1 <- sprintf("Empty Tree!")
      if(length(lst[[1]][[i]])==1)
        title(main = t, sub = st1) # if not an empty tree
      else
        title(main = t, sub = st) # if an empty tree
      v_par <- par('usr')
      xx <- v_par[1]
      yy <- v_par[4]
      info <- sprintf("log(Prior): %f\nlog(Posterior): %f\nPosterior odds: %f\nBIC: %f\nAIC: %f\nlog(ML): %f", lst[[2]]['log_prior'][[1]][[i]],
                      lst[[2]]['log_posterior'][[1]][[i]], lst[[2]]['posterior_odds'][[1]][[i]], lst[[2]]['BIC'][[1]][[i]]
                      ,lst[[2]]['AIC'][[1]][[i]], lst[[2]]['max_log_lik'][[1]][[i]])   # add the information in the plots
      
      text(x = xx,y = yy, info, adj = 0)
    }
  }
  
  else{ # bct results -> bct results does not contain the posterior_odds column
    show_tree(lst[[1]])
    t <- sprintf("BCT MAP tree")
    st <- sprintf("Max depth: %i \nNumber of leaves: %i", lst[[2]]['max_depth'][[1]], lst[[2]]['number_leaves'][[1]])
    st1 <- sprintf("Empty Tree!")
    if(length(lst[[1]])==1)
      title(main = t, sub = st1)
    else
      title(main = t, sub = st)
    v_par <- par('usr')
    xx <- v_par[1]
    yy <- v_par[4]
    info <- sprintf("log(Prior): %f\nlog(Posterior): %f\nBIC: %f\nAIC: %f\nmax_log_lik: %f", lst[[2]]['log_prior'][[1]],
                    lst[[2]]['log_posterior'][[1]], lst[[2]]['BIC'][[1]]
                    ,lst[[2]]['AIC'][[1]], lst[[2]]['max_log_lik'][[1]])
    
    text(x = xx,y = yy, info, adj = 0)
  }
}


# =========================================================================================

#' @title Sequence generator
#' @description Generates a simulated sequence of data according to a given model and associated parameters. 
#' An initial context of length equal to the maximum depth of the model is first generated uniformly and independently, and it is deleted after the desired  
#' number of samples has been generated. 
#' 
#' @param ct_theta a list containing the contexts that specify a model, and also a parameter vector for each context.
#' @param N length of the sequence to be generated. 
#' @return a simulated sequence as a "character" object
#' @export
#' @seealso \code{\link{BCT}}, \code{\link{kBCT}}, \code{\link{MAP_parameters}}
#' 
#' @examples
#' # Create a list containing contexts and associated parameters.
#' d1 <- list("0" = c(0.2, 0.8), "10" = c(0.9, 0.1), "11" = c(0,1))
#' 
#' # The contexts need to correspond to the leaves of a proper tree. 
#' # The key of each vector is the context. 
#' # For example:
#' # For context "0": P(x_{i+1} = 0 | x_{i} = 0) = 0.2
#' # and P(x_{i+1} = 1 | x_{i} = 0) = 0.8
#' 
#' # If a dataset containing only letters is desired:
#' d2 <- list("ab" = c(0.3, 0.7), "b" = c(0.8, 0.2), "aa" = c(0.5,0.5))
#' 
#' # Generate data from d2
#' gd <- generate_data(d2, 10000) 
#' 
#' # Use the BCT function to find the MAP model
#' BCT(gd, 10) # maximum depth of 10
#' 
#' # or the kBCT function can be used:
#' kBCT(gd, 10, 5) # maximum depth of 10 and top 5 models
generate_data <- function(ct_theta, N){
  contexts <- names(ct_theta)
  
  depth <- 0
  alphabet <- vector()

  # Constructs the alphabet and the maximal depth of the contexts
  for (context in contexts){
    alphabet <- union(alphabet, unlist(strsplit(context, "")))
    if(nchar(context)>depth)
      depth <- nchar(context)
  }
  
  alphabet<- sort(alphabet)
  alphabet_list <- list()
  
  # Constructs a dictionary. Keys: elements of the alphabet. Values: position in the sorted alphabet
  for(i in 1:length(alphabet))
    alphabet_list[alphabet[i]] = i
  
  if(depth == 0){
    x <- sample(alphabet, size = N, replace = TRUE, prob = ct_theta[["-"]])
    
    x <- paste(x, collapse='')
    return(x)
  }
  
  else{ 
    
    x <- integer(N+depth) # First depth elements are used as a random context
    for(i in 1:depth)
      x[i] <- sample(alphabet, size = 1, replace = TRUE)
    
    for (i in (depth+1):(N+depth)) {
      ct <- vector()
      for(j in 1:depth){ 
        ct <- c(ct, x[i-j])
        ct <-paste(ct,collapse="") # for each i, construct the context prior the the character on the ith position
        if(! is.null(ct_theta[[ct]])){
          x[i] <- sample(alphabet, size = 1, replace = TRUE, prob = ct_theta[[ct]])
          break
        }
      }
    }
    output <-  paste(x[(depth+1):(N+depth)], collapse = '')
    return(output)
  }
}

# =========================================================================================


