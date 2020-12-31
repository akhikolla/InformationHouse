.onAttach = function(...) {
  # Begin Exclude Linting
  greet = paste("# This research was partially supported under NSF grant 1463642")
  packageStartupMessage(greet)
  # End Exclude Linting
}

.onUnload <- function (libpath) {
  library.dynam.unload("smerc", libpath)
}
