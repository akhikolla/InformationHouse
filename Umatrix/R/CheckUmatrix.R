CheckUmatrix <- function(Umatrix = NULL, shiny = F){
  # Umatrix <- CheckUmatrix(Umxatrix)
  # Checks if the given Umatrix contains errors. Stops program if
  #    an error is found.
  # INPUT
  # Umatrix
  # shiny    is this method called from a shiny application?
  # OUTPUT
  # Umatrix
  # Author FL

  stopShinyOrR <- function(message){
    if(shiny) stopApp(message)
    else stop(message)
  }

  if(is.null(Umatrix)) return(NULL)
  if(!is.matrix(Umatrix)){
    stopShinyOrR("Umatrix is not given as a matrix")
  }
  if(min(Umatrix) < 0){
    stopShinyOrR("Umatrix contains values smaller than 0")
  }
  if(max(Umatrix) > 1){
    stopShinyOrR("Umatrix contains values greater than 1")
  }

  Umatrix
}
