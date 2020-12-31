.onAttach <- function(libname, pkgname) {
    if (interactive()) {
        packageStartupMessage('landscapeR ', as.character(utils::packageVersion("landscapeR")),' - For help type ?landscapeR or go to: cran.r-project.org/package=landscapeR \n',
                              'License: GPL v3. Centre for Ecology and Hydrology and NERC (UK)\n', domain=NA, appendLF=TRUE)
    }
}

.onUnload <- function (libpath) {
  library.dynam.unload("landscapeR", libpath)
}
