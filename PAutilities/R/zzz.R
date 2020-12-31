if(getRversion() >= "2.15.1") utils::globalVariables(".")

.onAttach <- function(libname, pkgname) {

  packageStartupMessage(
    "package 'PAutilities' was built under R version 3.5.0"
  )

}
