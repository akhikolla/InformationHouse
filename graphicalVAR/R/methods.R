summary.graphicalVAR <- function(object,...) print(object,...)

# Print method
print.graphicalVAR <- function(x, ...){
  name <- deparse(substitute(x))[[1]]
  if (nchar(name) > 10) name <- "object"
  

  cat("=== graphicalVAR results ===")
  cat("\nNumber of nodes:",nrow(x[['kappa']]),
  "\nNumber of tuning parameters tested:",nrow(x[['path']]),
  "\nEBIC hyperparameter:",x[['gamma']],
  "\nOptimal EBIC:",x[['EBIC']],
      
      
      "\n\nNumber of non-zero Partial Contemporaneous Correlations (PCC):",sum(x[['PCC']][upper.tri(x[['PCC']],diag=FALSE)]!=0) ,
      "\nPCC Sparsity:",mean(x[['PCC']][upper.tri(x[['PCC']],diag=FALSE)]==0) ,
      "\nNumber of PCC tuning parameters tested:",length(unique(x$path$kappa)),
      paste0("\nPCC network stored in ",name,"$PCC"),
      
      "\n\nNumber of non-zero Partial Directed Correlations (PDC):",sum(x[['PDC']][upper.tri(x[['PDC']],diag=FALSE)]!=0) ,
      "\nPDC Sparsity:",mean(x[['PDC']][upper.tri(x[['PDC']],diag=FALSE)]==0) ,
      "\nNumber of PDC tuning parameters tested:",length(unique(x$path$beta)),
      paste0("\nPDC network stored in ",name,"$PDC"),
  
      paste0("\n\nUse plot(",name,") to plot the estimated networks.")
  )
}

# Plot method
plot.graphicalVAR <- plot.gVARmodel <- function(x, include = c("PCC","PDC"), repulsion = 1, horizontal = TRUE, titles = TRUE, sameLayout = TRUE, unweightedLayout = FALSE,...){
  qtitle <-  function (x) 
  {
    text(par("usr")[1] + (par("usr")[2] - par("usr")[1])/40, 
         par("usr")[4] - (par("usr")[4] - par("usr")[3])/40, x, 
         adj = c(0, 1))
  }
  
  if (length(include)>1){
    if (horizontal){
      layout(t(seq_along(include))) 
    } else {
      layout(seq_along(include))
    }
  }
  
  # Choose directed or undirected:
  if (unweightedLayout){
    wPCC <- 1*(x$PCC!=0)
    wPDC <- 1*(x$PDC!=0)
  } else {
    wPCC <- x$PCC
    wPDC <- x$PDC
  }
  
  if (sameLayout & all(c("PCC","PDC") %in% include)){
    Layout <- qgraph::averageLayout(as.matrix(wPCC), as.matrix(wPDC), repulsion=repulsion)
  }
  
  Res <- list()
  
  for (i in seq_along(include)){
    if ("PCC" == include[i]){
      if (sameLayout & all(c("PCC","PDC") %in% include)){

        Res[[i]] <- qgraph::qgraph(x$PCC, layout = Layout, ..., repulsion=repulsion)
      } else {
        L <- qgraph::qgraph(wPCC,DoNotPlot=TRUE,...,repulsion=repulsion)$layout
        Res[[i]] <- qgraph::qgraph(x$PCC, layout = L,..., repulsion=repulsion)
      }

      if (titles){
        qtitle("Partial Contemporaneous Correlations")
      }
    }
    
    if ("PDC" == include[i]){
      if (sameLayout & all(c("PCC","PDC") %in% include)){
        Res[[i]] <- qgraph::qgraph(x$PDC, layout = Layout, ..., repulsion=repulsion, directed=TRUE)
      } else {
        L <- qgraph::qgraph(wPDC,DoNotPlot=TRUE,...,repulsion=repulsion, directed=TRUE)$layout
        Res[[i]] <- qgraph::qgraph(x$PDC,layout=L, ..., repulsion=repulsion, directed=TRUE)
      }
      
      if (titles){
        qtitle("Partial Directed Correlations")
      }
    }
  }
  
  invisible(Res)
}