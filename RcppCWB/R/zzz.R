.onAttach <- function (libname, pkgname) {
  # Previous versions of RcppCWB checked whether paths are set 
  # correctly in the registry files within the package and reverted to
  # a temporary directory only if paths are not set correctly.
  # Using temporary directories is now the standard to avoid issues
  # with staged installs, see:
  # https://developer.r-project.org/Blog/public/2019/02/14/staged-install/index.html
  # 
  # Note: Using use_tmp_registry in a .onAttach() call rather than
  # .onLoad() seems the appropriate solution, because activating the
  # corpora within the package is necessary only when the package is
  # loaded. When it is used as a backend (e.g. by polmineR), the 
  # corpora within the RcppCWB package remain unused, and using the
  # temporary directory is not necessary.
  use_tmp_registry(pkg = file.path(libname, pkgname))
  if (!cqp_is_initialized()) cqp_initialize()
}


.onUnload <- function(libpath) {
  # xml2 pkg served as model
  gc() # trigger finalisers
  library.dynam.unload("RcppCWB", libpath)
}