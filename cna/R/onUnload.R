
# http://markovjumps.blogspot.ch/2011/11/unloading-dynamic-libraries-when.html
.onUnload <- function(libpath) { library.dynam.unload("cna", libpath) }