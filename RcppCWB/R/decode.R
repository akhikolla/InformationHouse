#' Decode Structural Attribute.
#' 
#' Get \code{data.frame} with left and right corpus positions (cpos) for structural
#' attributes and values.
#' 
#' Two approaches are implemented: A pure R solution will decode the files directly in 
#' the directory specified by \code{data_dir}. An implementation using Rcpp will use the
#' registry file for \code{corpus} to find the data directory.
#' 
#' @param corpus a CWB corpus
#' @param s_attribute a structural attribute
#' @param data_dir data directory where binary files for corpus are stored
#' @param encoding encoding of the values ("latin-1" or "utf-8")
#' @param registry registry directory
#' @param method character vector, whether to use "R" or "Rcpp" implementation
#' @return A \code{data.frame} with three columns. Column \code{cpos_left} are the start
#' corpus positions of a structural annotation, \code{cpos_right} the end corpus positions.
#' Column \code{value} is the value of the annotation.
#' @export s_attribute_decode
#' @rdname s_attribute_decode
#' @return a character vector
#' @examples 
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' Sys.setenv(CORPUS_REGISTRY = registry)
#' 
#' # pure R implementation (Rcpp implementation fails on Windows in vanilla mode)
#' b <- s_attribute_decode(
#'   data_dir = system.file(package = "RcppCWB", "extdata", "cwb", "indexed_corpora", "reuters"),
#'   s_attribute = "places", method = "R"
#'   )
s_attribute_decode <- function(corpus, data_dir, s_attribute, encoding = NULL, registry = Sys.getenv("CORPUS_REGISTRY"), method = c("R", "Rcpp")){
  if (method == "R"){
    
    if (missing(data_dir)) stop("data_dir needs to be specified to use R method")
    avs_file <- file.path(data_dir, paste(s_attribute, "avs", sep = ".")) # attribute values
    avx_file <- file.path(data_dir, paste(s_attribute, "avx", sep = ".")) # attribute value index
    rng_file <- file.path(data_dir, paste(s_attribute, "rng", sep = ".")) # ranges
    
    avs_file_size <- file.info(avs_file)[["size"]]
    avx_file_size <- file.info(avx_file)[["size"]]
    rng_file_size <- file.info(rng_file)[["size"]]
    
    avs <- readBin(con = avs_file, what = character(), n = avs_file_size)
    rng <- readBin(rng_file, what = integer(), size = 4L, n = rng_file_size / 4, endian = "big")
    avx <- readBin(avx_file, what = integer(), size = 4L, n = avx_file_size / 4, endian = "big")
    
    if (!is.null(encoding)) Encoding(avs) <- encoding
    region_matrix <- matrix(rng, ncol = 2, byrow = TRUE)
    avx_matrix <- matrix(avx, ncol = 2, byrow = TRUE)
    
    colnames(region_matrix) <- c("cpos_left", "cpos_right")
    
    values <- avs[match(avx_matrix[,2], unique(avx_matrix[,2]))]
    
    df <- data.frame(
      cpos_left = region_matrix[,1],
      cpos_right = region_matrix[,2],
      value = values,
      stringsAsFactors = FALSE
    )
  } else if (method == "Rcpp"){
    
    if (missing(corpus)) stop("corpus needs to be specified to use 'Rcpp' method")
    check_registry(registry = registry)
    check_corpus(corpus = corpus, registry = registry)
    check_s_attribute(corpus = corpus, registry = registry, s_attribute = s_attribute)
    values <- .decode_s_attribute(corpus = corpus, s_attribute = s_attribute, registry = registry)
    warning("region matrix can not yet be generated with Rcpp method")
    region_matrix <- NULL
    
    df <- data.frame(
      cpos_left = NA,
      cpos_right = NA,
      value = values,
      stringsAsFactors = FALSE
    )
    
  } else {
    stop("method needs to be either 'R' or 'Rcpp'")
  }
  df
}
