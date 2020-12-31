#' Biplots of 
#' 
#' Construct simple two-dimensional biplots given matrices representing the rows
#' and columns of a two-dimensional matrix using \pkg{ggplot2}.
#' 
#' @param rows A list of matrices representing the rows
#' @param cols A list of matrices representing the columns
#' @importFrom ggplot2 ggplot geom_point aes facet_wrap geom_segment geom_text
#' @importFrom grid arrow unit
#' @export
#' @examples 
#' set.seed(1)
#' dat <- rlsbclust(ndata = 1, nobs = 100, size = c(10, 8), nclust = c(5, 4, 6, 5))
#' meanbiplot(dat[[1]]$interactions$C, dat[[1]]$interactions$D)
meanbiplot <- function(rows, cols) {
  
  ## Construct data frames for the means
  dfRow <- Reduce(rbind, rows)
  colnames(dfRow) <- c("x", "y")
  dfRow <- cbind(dfRow, Cluster = rep(seq_along(rows), each = nrow(rows[[1]])))
  dfRow <- as.data.frame(dfRow)
  dfRow$names <- paste0("R", rep(seq_len(nrow(rows[[1]])), times = length(rows)))
  dfCol <- Reduce(rbind, cols)
  colnames(dfCol) <- c("x", "y")
  dfCol <- cbind(dfCol, Cluster = rep(seq_along(cols), each = nrow(cols[[1]])))
  dfCol <- as.data.frame(dfCol)
  dfCol$xorigin <- rep(0, nrow(dfCol))
  dfCol$yorigin <- rep(0, nrow(dfCol))
  dfCol$names <- paste0("C", rep(seq_len(nrow(cols[[1]])), times = length(cols)))
  
  
  ## Create simple biplots
  p <- ggplot(data = dfRow, mapping = aes(x = x, y = y, group = Cluster, label = names)) + 
    geom_point() +  
    geom_text(nudge_y = 0.15, size = 2.5, colour = "grey40") +
    geom_segment(data = dfCol, mapping = aes(x = xorigin, xend = x, y = yorigin, yend = y), 
                 arrow = grid::arrow(angle = 20, length = grid::unit(0.0175, "npc"))) + 
    geom_text(data = dfCol, mapping = aes(label = names), nudge_y = 0.15, size = 2.5, colour = "grey40") +
    facet_wrap(~ Cluster) + theme_bw()
  
  
  ## Return
  return(p)
}