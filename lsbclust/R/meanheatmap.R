#' Plot Heatmap of A Matrix
#' 
#' Construct a heatmap of a matrix using \pkg{ggplot2}.
#' 
#' @param x Matrix or list of matrices to be plotted
#' @export
#' @importFrom ggplot2 ggplot geom_tile facet_wrap scale_fill_gradient2 theme_bw
#' @examples 
#' set.seed(1)
#' dat <- rlsbclust(ndata = 1, nobs = 100, size = c(6, 6), nclust = c(5, 4, 6, 5))
#' meanheatmap(Map(tcrossprod, dat[[1]]$interactions$C, dat[[1]]$interactions$D))
meanheatmap <- function(x) {
  
  if (is.matrix(x)) {
    ## Convert x to data.frame
    df <- data.frame(Value = c(x), Row = c(row(x)), Column = c(col(x)))
    if (!is.null(rownames(x))) df$Row <- rownames(x)[df$Row]
    if (!is.null(colnames(x))) df$Column <- colnames(x)[df$Column]
    df$Row <- factor(df$Row, levels = rev(seq_len(nrow(x))))
    df$Column <- factor(df$Column)
    
    ## Create plot
    p <- ggplot(data = df, mapping = aes(x = Column, y = Row, fill = Value)) +
      geom_tile() + scale_fill_gradient2() + theme_bw()
    
  } 
  if (is.list(x)) {
    
    ## Convert x to data.frame
    df <- lapply(x, function(x) data.frame(Value = c(x), Row = c(row(x)), Column = c(col(x))))
    df <- Reduce(rbind, df)
    if (!is.null(rownames(x))) df$Row <- rownames(x)[df$Row]
    if (!is.null(colnames(x))) df$Column <- colnames(x)[df$Column]
    df$Row <- factor(df$Row, levels = rev(seq_len(nrow(x[[1]]))))
    df$Column <- factor(df$Column)
    df$Cluster <- rep(seq_along(x), times = sapply(x, function(x) prod(dim(x))))
    
    ## Create plot
    p <- ggplot(data = df, mapping = aes(x = Column, y = Row, fill = Value, group = Cluster)) +
      geom_tile() + facet_wrap(~ Cluster) + scale_fill_gradient2() + theme_bw()
    
  }
  
  ## Return
  return(p)
}