CheckImx <- function(Imx = NULL, shiny = F){
  # Imx <- CheckImx(Imx)
  # Checks if the given Island makes sense. Stops program if an error is found.
  # INPUT
  # Imx
  # shiny    is this method called from a shiny application?
  # OUTPUT
  # Imx
  # Author FL

  stopShinyOrR <- function(message){
    if(shiny) stopApp(message)
    else stop(message)
  }

  if(is.null(Imx)) return(NULL)
  if(!is.matrix(Imx)) stopShinyOrR("Island is not given as a matrix")
  if(sum((Imx!=1)&(Imx!=0)) > 0) stopShinyOrR("Island contains values other than 0 and 1")

  Imx
}
