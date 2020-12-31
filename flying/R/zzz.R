.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to package flying")
}

.onUnload <- function (libpath) {
  library.dynam.unload("flying", libpath)
}
