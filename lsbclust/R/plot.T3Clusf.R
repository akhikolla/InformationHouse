#' Plot Method for Class 'T3Clusf'
#' 
#' Two-dimensional plot method for object of class 'T3Clusf' as output by \code{\link{T3Clusf}}.
#' 
#' @param x An object of class \code{T3Clusf}.
#' @param which An integer vector indicating which item segments to plot.
#' @param arrange Logical indicating whether to arrange the plots on a single page or not
#' @param \dots Additional arguments to \code{\link{theme}}
#' @keywords hplot
#' @method plot T3Clusf
#' @importFrom gridExtra grid.arrange
#' @export
plot.T3Clusf <- function(x, which = seq_len(nclust), arrange = FALSE, ...) {
  
  ## Break if means is NULL
  if (is.null(x$means))
    stop("The cluster means has not been calculated (alpha > 1).")

  ## Avoid globalVariable issues
  
  ## Determine dims and number of clusters
  nr <- nrow(x$means[[1]])
  nc <- ncol(x$means[[1]])
  nclust <- length(x$means)
  
  ## Set up data frame for means
  for (i in which) {
    dfMean <- data.frame(Value = x$means[[i]], Row = , Column = )
    
  }
  dfMeans <- do.call(rbind, x$means)
  dfMeans <- as.data.frame(reshape2::melt(dfMeans, as.is = TRUE))
  colnames(dfMeans) <- c("Row", "Column", "Value")
  
  ## Change/add factor and ensure correct order
  dfMeans$Row <- factor(dfMeans$Row, levels = rev(rownames(x$means[[1]])))
  dfMeans$Column <- factor(dfMeans$Column, levels = colnames(x$means[[1]]))
  dfMeans$Cluster <- rep(rep(seq_len(nclust), each = J), K)
  
  ## Get limits for colour range (symmetric)
  lims <- c(-1, 1) * max(abs(dfMeans$Value))
  
  ## Do plots
  for (i in which) {
    plots[[i]] <- ggplot(data = dfMeans[dfMeans$Cluster == i, ], aes(y = Row, x = Column, fill = Value)) + 
      geom_tile(colour = "white") + 
      scale_fill_gradient2(limits = lims) +
      # ggtitle(paste0("Mean: Interactions ", i, " (", round(100 * nvec[i] / N, 1),  "%)")) +
      theme_bw() + theme(...)
  }
  
  ## Set plot names
  names(plots) <- paste0("plot", seq_len(nclust))
  
  ## Return plots
  if (arrange) do.call(gridExtra::grid.arrange, plots[which])
  else return(plots[which])
}