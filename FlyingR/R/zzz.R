.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to package FlyingR")
}

.onUnload <- function (libpath) {
  library.dynam.unload("FlyingR", libpath)
}
