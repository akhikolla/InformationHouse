TopviewTopographicMap <- function(GeneralizedUmatrix,BestMatchingUnits,Cls,ClsColors=NULL,Imx=NULL,Names=NULL, BmSize=6,...) {
  #author: Tim Schreier, Luis Winckelmann, MCT
  udim <- dim(GeneralizedUmatrix)
  #Tiled needed for Imx
  if(!is.null(Imx)){
    Tiled=TRUE
  }
  
  # Error Catching ----
  if (missing(BestMatchingUnits)) {
    BestMatchingUnits = matrix(1, 2, 2)
    warning('BestMatchingUnits are missing.Creating a dummy..')
  }
  if (!is.matrix(BestMatchingUnits))
    stop('Bestmatches have to be a matrix')
  else
    b = dim(BestMatchingUnits)
  
  if (b[2] > 3 | b[2] < 2)
    stop(paste0('Wrong number of Columns of Bestmatches: ', b[2]))
  if (b[2] == 3) {
    BestMatchingUnits = BestMatchingUnits[, 2:3]
  }

  if (missing(Cls))
    Cls = rep(1, b[1])
  
  d = dim(GeneralizedUmatrix)
  if (is.null(d)) {
    stop('GeneralizedUmatrix Dimension is null. Please check Input')
  }
  
 # requireNamespace('matrixStats')
  mini=apply(BestMatchingUnits, 2, min,na.rm=TRUE)
  maxi=apply(BestMatchingUnits, 2, max,na.rm=TRUE)
  #mini = matrixStats::colMins(BestMatchingUnits, na.rm = TRUE)
 # maxi = matrixStats::colMaxs(BestMatchingUnits, na.rm = TRUE)
  if (sum(mini) < 2) {
    stop('Some Bestmatches are below 1 in X or Y/Columns or Lines')
  }
  if (d[1] < maxi[1]) {
    stop(
      paste0(
        'Range of Bestmatches',
        maxi[1],
        ' is higher than Range of GeneralizedUmatrix',
        d[1]
      )
    )
  }
  if (d[2] < maxi[2]) {
    stop(
      paste0(
        'Range of Bestmatches',
        maxi[2],
        ' is higher than Range of GeneralizedUmatrix',
        d[2]
      )
    )
  }
  
  if (!is.vector(Cls)) {
    warning('Cls is not a vector. Calling as.vector()')
    Cls = as.vector(Cls)
  }
  if (!is.numeric(Cls)) {
    warning('Cls is not a numeric Calling as.numeric()')
    Cls = as.numeric(Cls)
  }
  if (sum(!is.finite(Cls)) > 0) {
    warning('Not all values in Cls are finite. Generating nonfiniteclass with value 999')
    Cls[!is.finite(Cls)] = 999
  }

  if (length(Cls) != b[1]) {
    Cls = rep(1, b[1])
    warning(
      paste0(
        'Cls has the length ',
        length(Cls),
        ' which does not equal the number of the BestMatchingUnits: ',
        b[1],
        '. Plotting without Cls.'
      )
    )
  }
  
  #Handle Color ----
  colormap = GeneralizedUmatrix::UmatrixColormap
  
  if (is.null(ClsColors)) {
    ClsColors = GeneralizedUmatrix::DefaultColorSequence
    ClsColors = ClsColors[-5] #green is not visible in plotly
  } else{
    if (length(unique(Cls)) > length(ClsColors)) {
      stop('Length of vector of Clscolor does not match the number of unique Clusters in Cls.')
    }
  }
  
  ## Additional Arguments ----
  dots = list(...)
  #in case of pmatrix
  if (is.null(dots[["Tiled"]]))
    Tiled = FALSE
  else
    Tiled=dots[["Tiled"]]
  if(!is.null(Imx)){
    Tiled=TRUE
  }
  
  #axis with labels
  if (is.null(dots[["ShinyBinding"]]))
    ShinyBinding = FALSE
  else
    ShinyBinding=dots[["ShinyBinding"]]
  
  if (is.null(dots[["ShinyDimension"]]))
    ShinyDimension = 1
  else
    ShinyDimension=dots[["ShinyDimension"]]
  
  if (!is.null(dots[["Session"]]))
    session = dots[["Session"]]
  
  if (is.null(dots[["main"]]))
    main = NULL
  else
    main=dots[["main"]]
  
  #Helper Function ----
  addclass <- function(class,
                       plotbmus,
                       plot,
                       bmu_cols,
                       MarkerSize) {
    
    inds <- which(Cls == class)
    x = as.numeric(plotbmus[inds, 2])
    y = as.numeric(plotbmus[inds, 1])

    plot <- plotly::add_markers(
      plot, x=x, y=y,


    
      marker = list(
        size = MarkerSize,
        color = bmu_cols[class],
        line = list(color = "rgba(80, 80, 80, .8)",
                    width = 1)
      ),
      name = paste("Cluster", class)
    )

    return(plot)
  }
  
  PlotlyUmatrix = function(plotdim,
                           plotumx,
                           colormap,
                           Nrlevels2,
                           plotbmus,
                           class,
                           ClsColors,
                           MarkerSize,
                           ShinyBinding,
                           ShinyDimension,
                           Imx,Cls) {
    

    # configure filter, so that every bestmatch stays in
 
    
    # put Imx on Umatrix and bestmatches if given
    if(!is.null(Imx)){
      if(!is.null(plotbmus)){
        BestMatchesFilter = rep(T,nrow(plotbmus)) # every Bestmatch stays
        
      }
      plotumx[which(Imx == 1)] = 0
      bigImx = Imx
      #print(str(BestMatchesFilter))
      #print(str(plotbmus))
      for(i in 1:nrow(Imx)){
        for(j in 1:ncol(Imx)){
          if(Imx[i,j] == 1){
            #plotumx[i,j] = NA
            if(!is.null(plotbmus)){
              
              BestMatchesFilter[(plotbmus[,1] == i) & (plotbmus[,2] == j)] = F
            }
              
          }
        }
      }
      #BestMatchesFilter = rep(T,nrow(plotbmus)) # every Bestmatch stays
      #print(str(BestMatchesFilter))
      if(!is.null(plotbmus)) plotbmus = plotbmus[BestMatchesFilter,]

      #äprint(str(class))
      #print(class)
      #todo, ohne globale variable loesen
      if(!is.null(Cls)) Cls <<- Cls[BestMatchesFilter]
      
      oceanLine = apply(plotumx, 1, function(x) all(x==0))
      startLine = min(which(!oceanLine),na.rm=T)
      endLine = length(oceanLine) - min(which(rev(!oceanLine)),na.rm=T) + 1
      
      oceanCol = apply(plotumx, 2, function(x) all(x==0))
      startCol = min(which(!oceanCol),na.rm=T)
      endCol = length(oceanCol) - min(which(rev(!oceanCol)),na.rm=T) + 1
   
      if(!is.null(plotbmus)){
        plotbmus <- plotbmus - cbind(rep(startLine-1,nrow(plotbmus)),rep(startCol-1,nrow(plotbmus)))
      }
      plotumx <- plotumx[startLine:endLine,startCol:endCol]
      Imx <- Imx[startLine:endLine,startCol:endCol]
      bigImx <- bigImx[startLine:endLine,startCol:endCol]
      
    }
   


    ax <- list(
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      #showticklabels = FALSE,
      showgrid = FALSE
    )
    ay <- list(
      autorange = "reversed",
      title = "",
      zeroline = FALSE,
      showline = FALSE,
      #showticklabels = FALSE,
      showgrid = FALSE
    )
    
    if (isTRUE(ShinyBinding)) {
      width = (0.95 * as.numeric(ShinyDimension[1]))
      height = udim[1] / udim[2] * (width - 80)
      #print(width)
      plt <- plotly::plot_ly(width = width, height = height * 0.75)
    } else{
      plt <- plotly::plot_ly()
    }
    plt <- plotly::add_contour(
      plt,
      x = 1:plotdim[1],
      y = 1:plotdim[2],
      z = plotumx,
      showscale = FALSE,
      line = list(color = 'black',
                  width = 0.5),
      contours = list(
        start = 0,
        end = 1,
        size = 1 / 15
      ),
      # colors = colorRamp(colormap[c(rep(3, 6),
      #                               seq(
      #                                 from = 4,
      #                                 to = length(colormap) - 30,
      #                                 length.out = ceiling(Nrlevels2 + 1) - 7
      #                               ),
      #                               length(colormap))]),
      colors = colorRamp(colormap[c(  1,2,
                                      seq(
                                      from = 3,
                                      to = length(colormap) - 30,
                                      length.out = ceiling(Nrlevels2 + 1)- 4
                                    ),length(colormap),
                                    length(colormap))]),


      name = "UMatrix"
      # , showscale = FALSE
    )

    for (class in unique(Cls)) {
      
      plt <- addclass(class, plotbmus, plt, ClsColors, MarkerSize)
    }#end add class
 
    plt <- plotly::layout(
      #title <- "Drop plot title here",
      #bgcolor = "rgb(244, 244, 248)",
      plt,
      xaxis = ax,
      yaxis = ay,
      dragmode = 'lasso',
      legend = list(orientation = 'h')
      #, showlegend = FALSE
    )

    #if (isTRUE(ShinyBinding)) {
    #  requireNamespace('shiny')
    #  shiny::updateSelectInput(session,
    #                    "ClsSelect",
    #                    label = "Select Class",
    #                    choices = unique(Cls))
    #}

    return(plt)
  }
  
  if (missing(Cls))
    Cls = rep(1, nrow(BestMatchingUnits))
  
  #Normalizing GeneralizedUmatrix ----
  
 
  quants = quantile(as.vector(GeneralizedUmatrix), c(0.01, 0.5, 0.99))
  minU = quants[1]
  maxU = quants[3]
  GeneralizedUmatrix = (GeneralizedUmatrix - minU)/(maxU -
                                                      minU)
   quants2 = quantile(as.vector(GeneralizedUmatrix), c(0.01, 
                                                      0.5, 0.99))
  minU2 = quants2[1]
  maxU2 = quants2[3]
  HeightScale = round(maxU2 / (2 * max(minU2, 0.05)), 0)
  stretchFactor = sqrt(nrow(GeneralizedUmatrix) ^ 2 + ncol(GeneralizedUmatrix) ^
                         2) / sqrt(50 ^ 2 + 80 ^ 2)
  Nrlevels2 = 2 * HeightScale * stretchFactor
  indMax = which(GeneralizedUmatrix > 1, arr.ind = T)
  indMin = which(GeneralizedUmatrix < 0, arr.ind = T)
  if (length(indMax) > 0)
    GeneralizedUmatrix[indMax] = 1
  if (length(indMin) > 0)
    GeneralizedUmatrix[indMin] = 0

  # GeneralizedUmatrix <- GeneralizedUmatrix * HeightScale * stretchFactor
  if (isTRUE(Tiled)) {
    
    tU <- tileGUM(GeneralizedUmatrix,BestMatchingUnits,Cls)
    GeneralizedUmatrix <- tU$GeneralizedUmatrix
    BestMatchingUnits <- tU$BestMatchingUnits[,2:3] #no key
    Cls <- tU$Cls
    qdim <- udim * 2
    #dmx  <- cbind(GeneralizedUmatrix, GeneralizedUmatrix)
    #qmx  <- rbind(dmx, dmx)
    #dbm  <-
    #  rbind(BestMatchingUnits,
    #        cbind(BestMatchingUnits[, 1], BestMatchingUnits[, 2] + udim[2]))
    #qbm  <- rbind(dbm, cbind(dbm[, 1] + udim[1], dbm[, 2]))
    plotumx <- tU$GeneralizedUmatrix
    plotbmus <- tU$BestMatchingUnits[,2:3] #no key
    plotCls <-  tU$Cls
  } else{
  
    if(is.null(dots[["ExtendBorders"]])){
      #nothing
    }else{
      ExtendBorders=dots$ExtendBorders
      V=ExtendToroidalUmatrix(GeneralizedUmatrix,BestMatchingUnits,ExtendBorders)
      GeneralizedUmatrix=V$Umatrix
      BestMatchingUnits=V$Bestmatches
    }
    plotumx <- GeneralizedUmatrix
    plotbmus <- BestMatchingUnits
    plotCls <- Cls
    qdim <- udim
  }
  
  plotdim <- qdim

  plt=PlotlyUmatrix(
    plotdim,
    plotumx,
    colormap,
    Nrlevels2,
    plotbmus,
    class,
    ClsColors,
    BmSize,
    ShinyBinding,
    ShinyDimension,
    Imx,
    Cls
  )
  if(is.null(main))
    plt=plotly::layout(plt,title = "Topographic Map of Generalized U-Matrix")
  else
    plt=plotly::layout(plt,title = main)
  
  #if (isTRUE(ShinyBinding)) {
  #  PlotR <- plotly::renderPlotly({
  #    plt
  #  })
  #  return(list(Rendered=PlotR,single=plt))
  #} else{
    return(plt)
  #}
 
}#end TopviewTopographicMap
