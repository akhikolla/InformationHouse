#' Initialize Corpus Query Processor (CQP).
#' 
#' CQP needs to know where to look for CWB indexed corpora. To initialize CQP,
#' call \code{cqp_initialize}. To reset the registry, use the function
#' \code{cqp_reset_registry}. To get the registry used by CQP, use
#' \code{cqp_get_registry}. To get the initialization status, use
#' \code{cqp_is_initialized}
#' 
#' @param registry the registry directory
#' @export cqp_initialize
#' @rdname cqp_initialize
#' @author Andreas Blaette, Bernard Desgraupes, Sylvain Loiseau
#' @examples
#' cqp_is_initialized() # check initialization status
#' if (!cqp_is_initialized()) cqp_initialize()
#' cqp_is_initialized() # check initialization status (TRUE now?)
#' cqp_get_registry() # get registry dir used by CQP
#' 
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' if (cqp_get_registry() != registry) cqp_reset_registry(registry = registry)
#' cqp_list_corpora() # get list of corpora
cqp_initialize <- function(registry = Sys.getenv("CORPUS_REGISTRY")){
  registry_new <- registry
  # registry # necessary to capture Sys.getenv() assignment
  if (cqp_is_initialized()){
    warning("CQP has already been initialized. Re-initialization is not possible. ",
            "Only resetting registry.")
  } else {
    # workaround to ensure that global registry variable in dynamic
    # library will have 255 characters. Without starting with a very long (fake)
    # initial registry, there is a bug when resetting the registry dir to a dir
    # that is longer than the initial dir
    dummy_superdir <- tempdir()
    dir.create(dummy_superdir, showWarnings = FALSE)
    if (.Platform$OS.type == "windows"){
      dummy_superdir <- normalizePath(dummy_superdir, winslash = "/")
    }
    # the times argument is 247 for Windows compatibility 
    dummy_regdir <- file.path(
      dummy_superdir, 
      paste0(
        rep("x", times = 246 - nchar(dummy_superdir)),
        collapse = ""
      )
    )
    dir.create(dummy_regdir, showWarnings = FALSE)
    Sys.setenv(CORPUS_REGISTRY = dummy_regdir)
    .init_cqp()
  }
  check_registry(registry_new)
  Sys.setenv(CORPUS_REGISTRY = registry_new)
  cqp_reset_registry()
  return( cqp_is_initialized() )
}


#' @export cqp_is_initialized
#' @rdname cqp_initialize
cqp_is_initialized <- function(){
  if (.cqp_get_status() == 0) return(FALSE) else return(TRUE)
}

#' @export cqp_get_registry
#' @rdname cqp_initialize
cqp_get_registry <- function() .cqp_get_registry()

#' @export cqp_reset_registry
#' @rdname cqp_initialize
cqp_reset_registry <- function(registry = Sys.getenv("CORPUS_REGISTRY")){
  registry_dir <- registry
  if (!cqp_is_initialized()){
    warning("cannot reset registry, cqp has not yet been initialized!")
    return( FALSE )
  } else {
    check_registry(registry_dir)
    Sys.setenv(CORPUS_REGISTRY = registry_dir)
    if (nchar(registry_dir) > 255){
      stop("cannot assign new registry: maximum nchar(registry) is 255")
    } else {
      .cqp_set_registry(registry_dir = registry_dir)
      return( TRUE )
    }
  }
}


#' List Available CWB Corpora.
#' 
#' List the corpora described by the registry files in the registry directory
#' that is currently set.
#' 
#' @export cqp_list_corpora
#' @examples
#' if (!cqp_is_initialized()){
#'   registry <- system.file(package = "RcppCWB", "extdata", "cwb", "registry")
#'   cqp_initialize(registry)
#' }
#' cqp_list_corpora()
#' @author Andreas Blaette, Bernard Desgraupes, Sylvain Loiseau
cqp_list_corpora <- function() .cqp_list_corpora()


#' Execute CQP Query and Retrieve Results.
#' 
#' Using CQP queries requires a two-step procedure: At first, you execute a
#' query using \code{cqp_query}. Then, \code{cqp_dump_subcorpus} will return a
#' matrix with the regions of the matches for the query.
#' 
#' The \code{cqp_query} function executes a CQP query. The
#' \code{cqp_subcorpus_size} function returns the number of matches for the CQP
#' query. The \code{cqp_dump_subcorpus} function will return a two-column matrix
#' with the left and right corpus positions of the matches for the CQP query.
#' 
#' @param corpus a CWB corpus
#' @param query a CQP query
#' @param subcorpus subcorpus name
#' @export cqp_query
#' @rdname cqp_query
#' @references 
#' Evert, S. 2005. The CQP Query Language Tutorial. Available online at
#' \url{http://cwb.sourceforge.net/files/CWB_Encoding_Tutorial.pdf}
#' @examples 
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' 
#' if (!cqp_is_initialized()){
#'   cqp_initialize(registry = registry)
#' } else {
#'   if (cqp_get_registry() != registry) cqp_reset_registry(registry)
#' }
#' cqp_query(corpus = "REUTERS", query = '"oil";')
#' cqp_subcorpus_size("REUTERS")
#' cqp_dump_subcorpus("REUTERS")
#' 
#' cqp_query(corpus = "REUTERS", query = '"crude" "oil";')
#' cqp_subcorpus_size("REUTERS", subcorpus = "QUERY")
#' cqp_dump_subcorpus("REUTERS")
#' @author Andreas Blaette, Bernard Desgraupes, Sylvain Loiseau
cqp_query <- function(corpus, query, subcorpus = "QUERY"){
  stopifnot(corpus %in% cqp_list_corpora())
  query <- check_cqp_query(query)
  .cqp_query(corpus = corpus, subcorpus = subcorpus, query = query)
}

#' @export cqp_dump_subcorpus
#' @rdname cqp_query
cqp_dump_subcorpus <- function(corpus, subcorpus = "QUERY"){
  stopifnot(corpus %in% cqp_list_corpora())
  .cqp_dump_subcorpus(paste(corpus, subcorpus, sep = ":"))
}

#' @export cqp_subcorpus_size
#' @rdname cqp_query
cqp_subcorpus_size <- function(corpus, subcorpus = "QUERY"){
  stopifnot(corpus %in% cqp_list_corpora())
  .cqp_subcorpus_size(scorpus = paste(corpus, subcorpus, sep = ":"))
}

#' @export cqp_list_subcorpora
#' @rdname cqp_query
cqp_list_subcorpora <- function(corpus){
  stopifnot(corpus %in% cqp_list_corpora())
  .cqp_list_subcorpora(inCorpus = corpus)
}

