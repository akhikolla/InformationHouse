#' Get Matrix with Regions for Strucs.
#' 
#' The return value is an integer matrix with the left and right corpus positions
#' of the strucs in columns one and two, respectively.
#' 
#' @param corpus a CWB corpus
#' @param s_attribute a structural attribute
#' @param strucs strucs
#' @param registry the registry directory
#' @rdname get_region_matrix
#' @export get_region_matrix
#' @return A matrix with integer values indicating left and right corpus positions
#' (columns 1 and 2, respectively).
#' @examples 
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' y <- get_region_matrix(
#'   corpus = "REUTERS", s_attribute = "id",
#'   strucs = 0L:5L, registry = registry
#'   )
get_region_matrix <- function(corpus, s_attribute, strucs, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_strucs(corpus = corpus, s_attribute = s_attribute, strucs = strucs, registry = registry)
  .get_region_matrix(corpus = corpus, s_attribute = s_attribute, strucs = strucs, registry = registry)
}


#' Get IDs and Counts for Region Matrices.
#' 
#' @param corpus a CWB corpus
#' @param p_attribute a positional attribute
#' @param registry registry directory
#' @param matrix a regions matrix
#' @rdname region_matrix_ops
#' @name region_matrix_ops
#' @export region_matrix_to_ids
#' @examples
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' 
#' # Scenario 1: Get full text for a subcorpus defined by regions
#' m <- get_region_matrix(
#'   corpus = "REUTERS", s_attribute = "places",
#'   strucs = 4L:5L, registry = registry
#'   )
#' ids <- region_matrix_to_ids(
#'   corpus = "REUTERS", p_attribute = "word",
#'   registry = registry, matrix = m
#'   )
#' tokenstream <- cl_id2str(
#'   corpus = "REUTERS", p_attribute = "word",
#'   registry = registry, id = ids
#'   )
#' txt <- paste(tokenstream, collapse = " ")
#' txt
#' 
#' # Scenario 2: Get data.frame with counts for region matrix
#' y <- region_matrix_to_count_matrix(
#'   corpus = "REUTERS", p_attribute = "word",
#'   registry = registry, matrix = m
#'   )
#' df <- as.data.frame(y)
#' colnames(df) <- c("token_id", "count")
#' df[["token"]] <- cl_id2str(
#'   "REUTERS", p_attribute = "word",
#'   registry = registry, id = df[["token_id"]]
#'   )
#' df[order(df[["count"]], decreasing = TRUE),]
#' head(df)
region_matrix_to_ids <- function(corpus, p_attribute, registry = Sys.getenv("CORPUS_REGISTRY"), matrix){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_p_attribute(p_attribute = p_attribute, corpus = corpus, registry = registry)
  check_region_matrix(region_matrix = matrix)
  .region_matrix_to_ids(corpus = corpus, p_attribute = p_attribute, registry = registry, matrix = matrix)
}


#' @rdname region_matrix_ops
#' @export region_matrix_to_count_matrix
region_matrix_to_count_matrix <- function(corpus, p_attribute, registry = Sys.getenv("CORPUS_REGISTRY"), matrix){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_p_attribute(p_attribute = p_attribute, corpus = corpus, registry = registry)
  stopifnot(is.matrix(matrix))
  .region_matrix_to_count_matrix(corpus = corpus, p_attribute = p_attribute, registry = registry, matrix = matrix)
}

