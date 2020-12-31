#' Get CBOW Matrix.
#' 
#' Get matrix with moving windows. Negative integer values indicate absence of a
#' token at the respective position.
#' 
#' @param corpus a CWB corpus
#' @param p_attribute a positional attribute
#' @param registry the registry directory
#' @param matrix a matrix
#' @param window window size
#' @rdname get_cbow_matrix
#' @export get_cbow_matrix
#' @examples 
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' 
#' m <- get_region_matrix(
#'   corpus = "REUTERS", s_attribute = "places",
#'   strucs = 0L:5L, registry = registry
#'   )
#' windowsize <- 3L
#' m2 <- get_cbow_matrix(
#'   corpus = "REUTERS", p_attribute = "word",
#'   registry = registry, matrix = m, window = windowsize
#'   )
#' colnames(m2) <- c(-windowsize:-1, "node", 1:windowsize)
get_cbow_matrix <- function(corpus, p_attribute, registry = Sys.getenv("CORPUS_REGISTRY"), matrix, window){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_p_attribute(p_attribute = p_attribute, corpus = corpus, registry = registry)
  check_region_matrix(region_matrix = matrix)
  stopifnot(window >= 1L)
  .get_cbow_matrix(corpus = corpus, p_attribute = p_attribute, registry = registry, matrix = matrix, window = window)
}





