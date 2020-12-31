# delayed registration of autoplot, xyplot
vctrs_s3_register = function(generic, class, method = NULL) {
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)
  pieces <- strsplit(generic, "::")[[1]]
  stopifnot(length(pieces) == 2)
  package <- pieces[[1]]
  generic <- pieces[[2]]
  caller <- parent.frame()
  get_method_env <- function() {
    top <- topenv(caller)
    if (isNamespace(top)) {
      asNamespace(environmentName(top))
    }
    else {
      caller
    }
  }
  get_method <- function(method, env) {
    if (is.null(method)) {
      get(paste0(generic, ".", class), envir = get_method_env())
    }
    else {
      method
    }
  }
  method_fn <- get_method(method)
  stopifnot(is.function(method_fn))
  setHook(packageEvent(package, "onLoad"), function(...) {
    ns <- asNamespace(package)
    method_fn <- get_method(method)
    registerS3method(generic, class, method_fn, envir = ns)
  })
  if (!isNamespaceLoaded(package)) {
    return(invisible())
  }
  envir <- asNamespace(package)
  if (exists(generic, envir)) {
    registerS3method(generic, class, method_fn, envir = envir)
  }
  invisible()
}

.onLoad <- function(...) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    vctrs_s3_register("ggplot2::autoplot", "evgram")
  }
  if (requireNamespace("lattice", quietly = TRUE)) {
    vctrs_s3_register("lattice::xyplot", "evgram")
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("gear", libpath)
}

# @rawNamespace if(getRversion() >= "3.6.0") {  # delayed registration
#    S3method(ggplot2::autoplot, evgram)
# } else {
#    S3method(autoplot, evgram)
#    export(autoplot.evgram)
# }
# @rawNamespace if(getRversion() >= "3.6.0") {  # delayed registration
#    S3method(lattice::xyplot, evgram)
# } else {
#    S3method(xyplot, evgram)
#    export(xyplot.evgram)
# }
