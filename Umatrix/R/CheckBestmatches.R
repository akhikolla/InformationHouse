CheckBestMatches <- function(BestMatches= NULL, Cls = NULL, shiny = F){
# BestMatches <- CheckBestMatches(BestMatches, Cls)
# Checks if the given BestMatches and their cls data make sense. Stops program if
#    an error is found.
# INPUT
# BestMatches
# Cls
# shiny    is this method called from a shiny application?
# OUTPUT
# BestMatches
# Author FL

  if(!is.null(Cls)){
    if(!is.vector(Cls)){
      if(dim(Cls)[2]>1){
        if(shiny){
          stopApp('Wrong number of dimensions for Cls')
        }else{
          stop('Wrong number of dimensions for Cls')
        }
      }
      else
        Cls = as.vector(Cls)
    }
  }

  stopShinyOrR <- function(message){
    if(shiny) stopApp(message)
    else stop(message)
  }

  if(!is.null(BestMatches)){
    if(!(ncol(BestMatches) %in% c(2,3)))
      stopShinyOrR("BestMatches need to have 2 or 3 columns")
    if(ncol(BestMatches)==2){
      BestMatches = cbind(1:nrow(BestMatches),BestMatches)


      warning('Given BMUs do not have IDs. A column 1:nrow(BestMatches) will be added.')
    }
    if(!is.null(Cls)){
      if(!is.vector(Cls)){
          Cls=as.vector(Cls)
          warning('Cls has to be a vector, Umatrix() called as.vector(). If an follwing error occurs, please change Cls to vector format manually')
      }
      if(!is.numeric(Cls)){
          Cls=as.numeric(Cls)
          warning('Cls has to be numerical, Umatrix() called as.numeric(). If an follwing error occurs, please change Cls to numerical format manually')
      }
      if(nrow(BestMatches) != length(Cls))
        stopShinyOrR("BestMatches need to have the same amount of elements like Cls")
    }
    if(!is.numeric(BestMatches))
      stopShinyOrR('BestMatches has to be numerical, try to fix the problem manually')
    if(min(BestMatches[,c(2,3)]) < 1){
      warning('There are BestMatches with a row or column position below 1. Their position was artificially set to 1!')
      BestMatches[BestMatches[,2] < 1,2] = 1
      BestMatches[BestMatches[,3] < 1,3] = 1
    }
  }

  BestMatches
}
