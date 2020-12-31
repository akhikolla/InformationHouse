#' Get Vector with Counts for Positional Attribute.
#' 
#' The return value is an integer vector. The length of the vector is the number of 
#' unique tokens in the corpus / the number of unique ids. The order of the counts
#' corresponds to the number of ids.
#' 
#' @param corpus a CWB corpus
#' @param p_attribute a positional attribute
#' @param registry registry directory
#' @return an integer vector
#' @rdname get_count_vector
#' @export get_count_vector
#' @examples 
#' registry <- use_tmp_registry()
#' y <- get_count_vector(
#'   corpus = "REUTERS", p_attribute = "word",
#'   registry = registry
#'   )
#' df <- data.frame(token_id = 0:(length(y) - 1), count = y)
#' df[["token"]] <- cl_id2str(
#'   "REUTERS", p_attribute = "word",
#'   id = df[["token_id"]], registry = registry
#'   )
#' df <- df[,c("token", "token_id", "count")] # reorder columns
#' df <- df[order(df[["count"]], decreasing = TRUE),]
#' head(df)
get_count_vector <- function(corpus, p_attribute, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_p_attribute(p_attribute = p_attribute, corpus = corpus, registry = registry)
  .get_count_vector(corpus = corpus, p_attribute = p_attribute, registry = registry)
}



#' Perform Count for Vector of IDs.
#' 
#' The return value is a two-column integer matrix. Column one represents the
#' unique ids of the input vector, column two the respective number of
#' occurrences / counts.
#' 
#' @param ids a vector of ids (integer values)
#' @rdname ids_to_count_matrix
#' @examples 
#' ids <- c(1L, 5L, 5L, 7L, 7L, 7L, 7L)
#' ids_to_count_matrix(ids)
#' table(ids) # alternative to get a similar result
ids_to_count_matrix <- function(ids){
  stopifnot(is.integer(ids))
  .ids_to_count_matrix(ids = ids)
}


