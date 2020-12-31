#' Use Temporary Registry
#' 
#' Use and get temporary registry directory to describe and access the corpora
#' in a package.
#' @param pkg Full path to a package.
#' @export use_tmp_registry
#' @rdname tmp_registry
use_tmp_registry <- function(pkg = system.file(package = "RcppCWB")){
  
  tmp_registry_dir <- file.path(tempdir(), "registry_tmp")
  if (!file.exists(tmp_registry_dir)) dir.create(tmp_registry_dir)
  
  pkg_registry_dir <- file.path(pkg, "extdata", "cwb", "registry")
  pkg_indexed_corpora_dir <- file.path(pkg, "extdata", "cwb", "indexed_corpora")
  
  for (corpus in list.files(pkg_registry_dir)){
    
    registry <- readLines(file.path(pkg_registry_dir, corpus))
    
    home_line_no <- grep("^HOME", registry)
    registry[home_line_no] <- sprintf("HOME \"%s\"", file.path(pkg_indexed_corpora_dir, corpus))
    
    info_line_no <- grep("^INFO", registry)
    registry_info_file <- gsub("^INFO\\s+\"*(.*?)\"*\\s*$", "\\1", registry[info_line_no])
    info_file_new <- file.path(pkg_indexed_corpora_dir, corpus, basename(registry_info_file), fsep = "/")
    registry[info_line_no] <- sprintf("INFO \"%s\"", info_file_new)
    
    writeLines(text = registry, con = file.path(tmp_registry_dir, corpus), sep = "\n")
  }
  
  Sys.setenv("CORPUS_REGISTRY" = tmp_registry_dir)
  if (cqp_is_initialized()) cqp_reset_registry(tmp_registry_dir)
  tmp_registry_dir
}

#' @rdname tmp_registry
#' @export get_tmp_registry
get_tmp_registry <- function() file.path(tempdir(), "registry_tmp")



#' Get Registry Directory Within Package
#' @param pkgname Name of package (character vector)
#' @export get_pkg_registry
get_pkg_registry <- function(pkgname = "RcppCWB") system.file(package = pkgname, "extdata", "cwb", "registry")