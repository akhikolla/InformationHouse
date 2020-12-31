upscaleUmatrix <- function(Umatrix, Factor = 2, BestMatches, Imx){
  # scaled = upscaleUmatrix(umx, 2)
  ## upscaleUmatrix skaliert eine Umatrix mittels linearer Interpolation auf ein gr??eres Gitter.
  #
  ## INPUT
  # Umatrix       Eingabe Umatrix die hochskaliert werden soll.
  #
  ## OPTIONAL
  # Factor        Faktor um den die Umatrix pro Achse skaliert werden soll. Standart: 2. Achtung: Die Anzahl der Gitterpunkte w?chst um Faktor Factor^2.
  # Imx           Inselausschnitt der Umatrix, der auf die neue gr??e der Umatrix hochskaliert werden muss.
  #
  ## OUTPUT
  # Liste mit 2 Elementen:
  #   Umatrix     Skalierte Umatrix
  #   Imx         Wenn Imx angegeben wurde: entsprechend skalierte Imx matrix. Ansonsten: NULL
  #
  ## AUTHOR
  # Felix Pape 10/2016
  
  rv = list()
  
  if(Factor %% 1 != 0 || Factor < 2)
    stop("Factor must be a natural number greater 2")
  
  up.vector <- function (x, f = Factor, y = NA) 
  {
    n <- length(x)
    as.vector(rbind(x, matrix(rep(y, (f - 1) * n), nrow = f - 1)))
  }
  
  interpNA <- function (napos, data, f = Factor){
    x = napos[1]
    y = napos[2]
    n = (x-1) %% f
    
    if(x-n+f > dim(data)[1])
      i = 1
    else
      i = x - n + f
    
    l = data[i,y] - data[x-n,y]
    return(l/f * n + data[x-n,y])
  }
  
  minNA <- function (napos, data, f = Factor){
    x = napos[1]
    y = napos[2]
    n = (x-1) %% f
    
    if(x-n+f > dim(data)[1])
      i = 1
    else
      i = x - n + f
    
    return(min(data[i,y],data[x-n,y]))
  }
  
  fillNA <- function(datawithnarows, Factor = Factor, method = interpNA){
    naposs = which(is.na(datawithnarows),arr.ind = T)
    datawithnarows[which(is.na(datawithnarows))] = apply(naposs,1, method, data = datawithnarows, f = Factor)
    return(datawithnarows)
  }
  
  upscale <- function(data, method = interpNA){
    # Add rows
    upsampledrows = apply(data,2,up.vector, f = Factor)
    # Fill rows
    upsampledrows = fillNA(upsampledrows, Factor, method = method)
    # Add columns
    upsampledcols = apply(upsampledrows,1,up.vector, f = Factor)
    # Fill columns
    upsampledcols = fillNA(upsampledcols, Factor, method = method)
    # Transpose, so input columns are columns again (apply with margin 1 transposed the output)
    return(t(upsampledcols))
  }
  
  umx = upscale(Umatrix)
  rv$Umatrix = umx
  
  if(!missing(BestMatches)){
    rv$BestMatches = BestMatches * Factor - (Factor - 1)
  }
  
  if(!missing(Imx)){
    Imx = upscale(Imx,method = minNA)
    rv$Imx = Imx
  }
  
  return(rv)
}