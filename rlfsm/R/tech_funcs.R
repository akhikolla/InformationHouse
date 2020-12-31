
# Cleans up after unloading
.onUnload <- function (libpath) {
  library.dynam.unload("rlfsm", libpath)
}
