weighted_node_positions_output_row_names <- function(x, NZ, NP){
  if(is.null(rownames(x)) & is.null(colnames(x))){
    c(paste('r', 1:NZ, sep = ''), paste('c', 1:NP, sep = ''))
  } else if(is.null(rownames(x)) & (!is.null(colnames(x)))) {
    c(paste('r', 1:NZ, sep = ''), colnames(x))
  } else if((!is.null(rownames(x))) & is.null(colnames(x))) {
    c(rownames(x), paste('c', 1:NP, sep = ''))
  } else if((!is.null(rownames(x))) & (!is.null(colnames(x)))) {
    c(rownames(x), colnames(x))
  }
}
