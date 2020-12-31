plotTopographicMap <- function(GeneralizedUmatrix, BestMatchingUnits, Cls=NULL, ClsColors=NULL,Imx=NULL,Names=NULL, BmSize=0.5,RenderingContourLines=TRUE,...){
# plotTopographicMap(GeneralizedUmatrix, BestMatchingUnits, Cls, Tiled)
# Draws a plot of given GeneralizedUmatrix
# INPUT
# OPTIONAL
# GeneralizedUmatrix(1:Lines,1:Columns)	      GeneralizedUmatrix to be plotted
# BestMatchingUnits(1:n,1:2)		      Positions of BestMatchingUnits to be plotted onto the GeneralizedUmatrix
# Cls(1:n)			            class identifier for the bestmatch at the given point
# ClsColors(1:m)			Vector of colors that will be used to colorize the different classes
# Tiled				                  should the GeneralizedUmatrix be drawn 4times?
# HeightScale                 how high should the mountains be?
# TextureRendering            use the 2D GeneralizedUmatrix Function to create a Texture and use this
              
# author: MT
  #Bei einer Insel muss die Umatrix vervierfacht werden damit funktion funktioniert
 
########################################################################################## 
  #Catch further arguments ----
##########################################################################################
  if(missing(GeneralizedUmatrix)) stop('GeneralizedUmatrix is missing.')
  if(!exists(x = "GeneralizedUmatrix",where = parent.frame())) stop('GeneralizedUmatrix is missing.')
  if(is.null(GeneralizedUmatrix)) stop('GeneralizedUmatrix is missing.')
  
  #
  dots=list(...)
  
  #in case of pmatrix
  if(is.null(dots[["Colormap"]]))
    Colormap=GeneralizedUmatrix::UmatrixColormap
  else
    Colormap=dots$Colormap
  
  #axis with labels
  if(is.null(dots[["ShowAxis"]]))
    ShowAxis=FALSE
  else
    ShowAxis=dots$ShowAxis
  
  if(is.null(dots[["Tiled"]]))
    Tiled=FALSE
  else
    Tiled=dots$Tiled
  
  if(is.null(dots[["Silent"]]))
    Silent=TRUE
  else
    Silent=dots$Silent
  
  #number of contour lines, limit to plot faster
  if(is.null(dots[["NoLevels"]]))
    NoLevels=30
  else
    NoLevels=dots$NoLevels
  
  if(is.null(dots[["main"]]))
    main=NULL
  else
    main=dots$main
  
  if(!is.null(dots[["title"]]))
    main=dots$title
  
  
  if(is.null(dots[["sub"]]))
    sub=NULL
  else
    sub=dots$sub
  
  if(!ShowAxis){
    if(is.null(dots[["xlab"]]))
      xlab=NULL
    else
      xlab=dots$xlab
    
    if(is.null(dots[["ylab"]]))
      ylab=NULL
    else
      ylab=dots$ylab
    
    if(is.null(dots[["zlab"]]))
      zlab=NULL
    else
      zlab=dots$zlab
  }else{
    xlab=ylab=zlab=NULL
  }
 #########################################################################################  
  #check for island ----
##########################################################################################
if(!is.null(Imx))
  Tiled=TRUE

  #TextureRendering=F only with Umatrix package

  # OUTPUT
	if(!requireNamespace("rgl", quietly = T)) warning("Namespace of package rgl could not be loaded. You can try 'TopviewTopographicMap' function as an alternative for plotting.")

  # N#ormalization der GeneralizedUmatrix werte ----
  # Milligan, Copper 1988 A Study of Standadization of Variables in Cluster Analysis,
  # robust Normalization Z_5 :"Z_5 is bounded by 0.0 and 1.0 with at least one observed value at each of these end points"

  quants=quantile(as.vector(GeneralizedUmatrix),c(0.01,0.5,0.99))
  minU=quants[1]
  maxU=quants[3]

  GeneralizedUmatrix=(GeneralizedUmatrix-minU)/(maxU-minU)

  quants2=quantile(as.vector(GeneralizedUmatrix),c(0.01,0.5,0.99))
  minU2=quants2[1]
  maxU2=quants2[3]

  ### Hoehe aus GeneralizedUmatrix schaetzen ----
  #Verhaeltnis zwischen minhoehe/maxHoehe=1/HeightScale
      HeightScale=round(maxU2/(2*max(minU2,0.05)),0) #Factor 2 damit die GeneralizedUmatrix Hoehe nicht zu riesig wird

      # proportionen sind fuer standard GeneralizedUmatrix der groesse 80x50 gedacht. fuer andere groessen linear skalieren
      stretchFactor = sqrt(nrow(GeneralizedUmatrix)^2 + ncol(GeneralizedUmatrix)^2) / sqrt(50^2 + 80^2)

  #MT: Die Level muessen vor der Begrenzung der Werte auf 0 und 1 gesetz werden,
  #aber nachdem der Wertebereich umgeschoben wurde, damit die Dichte in den Grenzbereichen abgeschetzt werden kann

  # proportionen sind fuer standard GeneralizedUmatrix der groesse 80x50 gedacht. fuer andere groessen linear skalieren
  stretchFactor = sqrt(nrow(GeneralizedUmatrix)^2 + ncol(GeneralizedUmatrix)^2) / sqrt(50^2 + 80^2)

  indMax=which(GeneralizedUmatrix>1,arr.ind=T)
  indMin=which(GeneralizedUmatrix<0,arr.ind=T)
  if(length(indMax)>0)
    GeneralizedUmatrix[indMax]=1
  if(length(indMin)>0)
    GeneralizedUmatrix[indMin]=0
##########################################################################################
  #diverse Checks ----
##########################################################################################
  if(missing(BestMatchingUnits)){
    BestMatchingUnits=matrix(1,2,2)
    warning('BestMatchingUnits are missing.Creating a dummy..')
  }
  if (!is.matrix(BestMatchingUnits))
    stop('Bestmatches have to be a matrix')
  else
    b = dim(BestMatchingUnits)
  
  if (b[2] > 3 | b[2] < 2)
    stop(paste0('Wrong number of Columns of Bestmatches: ', b[2]))
  if (b[2] == 2) {
    Points = BestMatchingUnits
    #With Key
    BestMatchingUnits = cbind(1:b[1], BestMatchingUnits)
  } else{
    Points = BestMatchingUnits[, 2:3]
  }
  d = dim(GeneralizedUmatrix)
  if (is.null(d)) {
    stop('GeneralizedUmatrix Dimension is null. Please check Input')
  }
  mini=apply(Points, 2, min,na.rm=TRUE)
  maxi=apply(Points, 2, max,na.rm=TRUE)
  #requireNamespace('matrixStats')
  #mini = matrixStats::colMins(Points, na.rm = TRUE)
  #maxi = matrixStats::colMaxs(Points, na.rm = TRUE)
  if (sum(mini) < 2) {
    stop('Some Bestmatches are below 1 in X or Y/Columns or Lines')
  }
  if (d[1] < maxi[1]) {
    stop(paste0(
      'Range of Bestmatches',
      maxi[1],
      ' is higher than Range of GeneralizedUmatrix',
      d[1]
    ))
  }
  if (d[2] < maxi[2]) {
    stop(paste0(
      'Range of Bestmatches',
      maxi[2],
      ' is higher than Range of GeneralizedUmatrix',
      d[2]
    ))
  }
  ##Cls checks
if(is.null(Cls)) Cls=rep(1,b[1])
if(!is.vector(Cls)){
  warning('Cls is not a vector. Calling as.vector()')
  Cls=as.vector(Cls)
}
  if(!is.numeric(Cls)){
    warning('Cls is not a numeric Calling as.numeric()')
    Cls=as.numeric(Cls)
  }   
  if(sum(!is.finite(Cls))>0){
    warning('Not all values in Cls are finite. Generating nonfiniteclass with value 999')
    Cls[!is.finite(Cls)]=999
  }
  if(length(Cls)!=b[1]){
    Cls=rep(1,b[1])
    warning(paste0('Cls has the length ',length(Cls),' which does not equal the number of the BestMatchingUnits: ',b[1],'. Plotting without Cls.'))
  }
    
if(is.null(ClsColors)){
 ClsColors=GeneralizedUmatrix::DefaultColorSequence
 }else{
	if(length(unique(Cls))!=length(ClsColors)){
		stop('Length of vector of Clscolor does not match the number of unique Clusters in Cls.')
	}
} 

  
  #Settings for Legend ----
  if(!is.null(Names)){
    
    if(!is.character(Names)){
      warning('Names are not a character vector, calling as.character')
      Names=as.character(Names)
    }
    
    # 
     if(!is.null(Names)){
       cls_length=length(unique(Cls))
       if(!length(Names)==cls_length){
         warning('Names are not of length of number of clusters, Renaming Names to Cluster 1.. Cluster n')
         Names=rep("C",cls_length)
         no_vec=seq_len(cls_length)
         Names=paste0(Names,no_vec)
      }
     }
    
    if(is.null(dots[["NamesCex"]]))
      NamesCex=1.5
    else
      NamesCex=dots$NamesCex
    
    if(is.null(dots[["NamesPosition"]]))
      NamesPosition="topright"
    else
      NamesPosition=dots$NamesPosition
    
    if(is.null(dots[["NamesTitle"]]))
      NamesTitle="Clustering"
    else
      NamesTitle=dots$NamesTitle
    
    if(is.null(dots[["NamesPch"]]))
      NamesPch=16
    else
      NamesPch=dots$NamesPch
    
    if(is.null(dots[["NamesColors"]]))
      NamesColors=DataVisualizations::DefaultColorSequence[1:length(Names)]
    else{
      NamesColors=dots$NamesColors
      if(length(NamesColors)!=length(Names)){
        warning('Length of "Names" does not equal length of "NamesColors". Trying to use "ClsColors".')
        if(length(ClsColors)!=length(Names)){
          warning('Length of "Names" does not equal length of "ClsColors". Using default Colors')
          NamesColors=DataVisualizations::DefaultColorSequence[1:length(Names)]
        }else{
          NamesColors=ClsColors
        }
      }
    }
  }
  
  if(isFALSE(Silent)){
    cat("plotToporaphicMaps: Error checks complete.\n")
  }

##########################################################################################
## 4fach Kachelung (tiling) durchfuehren ----
##########################################################################################
  
  if(Tiled){
    tU <- tileGUM(GeneralizedUmatrix,BestMatchingUnits,Cls)
    GeneralizedUmatrix <- tU$GeneralizedUmatrix
    BestMatchingUnits <- tU$BestMatchingUnits #with key
    Cls <- tU$Cls
  }else{
    if(is.null(dots[["ExtendBorders"]])){
      #nothing
    }else{
      ExtendBorders=dots$ExtendBorders
      V=ExtendToroidalUmatrix(GeneralizedUmatrix,BestMatchingUnits[,2:3],ExtendBorders)
      GeneralizedUmatrix=V$Umatrix
      BestMatchingUnits=cbind(BestMatchingUnits[,1],V$Bestmatches)
    }
  }
##########################################################################################
## groessere imx bauen ----
#########
  bigImx = Imx

  ##########
  ## an GeneralizedUmatrix rumrechnen
  ##########
  #Aus showUmatrix3d, package Umatrix
  zcol = cut(GeneralizedUmatrix,128) # divide the Height into 128 different levels for coloring

  lines = seq(1, nrow(GeneralizedUmatrix), len = nrow(GeneralizedUmatrix))
  columns = seq(1, ncol(GeneralizedUmatrix), len = ncol(GeneralizedUmatrix))
  
  Nrlevels2 = 2*HeightScale*stretchFactor #MT Farbintervalle gleich 2*hoehenintervalle, siehe oben
  levelBreaks <- seq(0,1.000001,length.out=(Nrlevels2+1))

##########################################################################################
## remove Heights by island ----
##########################################################################################
  
  if(!is.null(Imx)){#Aus showUmatrix3d, package Umatrix
    GeneralizedUmatrix[which(Imx == 1)] = 0
    #GeneralizedUmatrix[which(bigImx == 1)] = 0

    if(!is.null(BestMatchingUnits)){
      BestMatchingUnitsFilter = rep(T,nrow(BestMatchingUnits)) # every Bestmatch stays
      for(i in 1:nrow(Imx)){
        for(j in 1:ncol(Imx)){
          if(Imx[i,j] == 1){
            if(!is.null(BestMatchingUnits))
              BestMatchingUnitsFilter[(BestMatchingUnits[,2] == i) & (BestMatchingUnits[,3] == j)] = F
          }
        }
      }
      BestMatchingUnits = BestMatchingUnits[BestMatchingUnitsFilter,]
      if(!is.null(Cls)) Cls = Cls[BestMatchingUnitsFilter]
    }
      #### remove ocean around GeneralizedUmatrix
      oceanLine = apply(GeneralizedUmatrix, 1, function(x) all(x==0))
      startLine = min(which(!oceanLine),na.rm=T)
      endLine = length(oceanLine) - min(which(rev(!oceanLine)),na.rm=T) + 1

      oceanCol = apply(GeneralizedUmatrix, 2, function(x) all(x==0))
      startCol = min(which(!oceanCol),na.rm=T)
      endCol = length(oceanCol) - min(which(rev(!oceanCol)),na.rm=T) + 1

      if(!is.null(BestMatchingUnits)){
        BestMatchingUnits <- BestMatchingUnits - cbind(rep(0,nrow(BestMatchingUnits)),startLine-1,startCol-1)
      }
      GeneralizedUmatrix <- GeneralizedUmatrix[startLine:endLine,startCol:endCol]
      Imx <- Imx[startLine:endLine,startCol:endCol]
      bigImx <- bigImx[startLine:endLine,startCol:endCol]
  }




  # change GeneralizedUmatrix into heights
  splittedGeneralizedUmatrix = GeneralizedUmatrix
  for(i in 1:Nrlevels2){
    splittedGeneralizedUmatrix[ (GeneralizedUmatrix >= levelBreaks[i]) & (GeneralizedUmatrix <= levelBreaks[i+1]) ] = levelBreaks[i]
  }
  splittedGeneralizedUmatrix=(floor(splittedGeneralizedUmatrix * length(Colormap)))+1
  color = Colormap[splittedGeneralizedUmatrix]
  #if(!is.null(Imx)) color[which(Imx == 1)] = NA
  if(!is.null(Imx)) color[which(bigImx == 1)] = NA

  z<-GeneralizedUmatrix*HeightScale*stretchFactor #MT: Die Hoehe entspricht den GeneralizedUmatrixwerten, siehe vorne Proportion

  lines = seq(1, nrow(z), len = nrow(z))
  columns = seq(1, ncol(z), len = ncol(z))


  
  if(isFALSE(Silent)){
    cat("plotToporaphicMaps: Heights computed.\n")
  }
##########################################################################################
## GeneralizedUmatrix darstellen ----
##########################################################################################
 #Aus showUmatrix3d, package Umatrix
   rgl::open3d()
   rgl::rgl.viewpoint(0, -30,zoom = 0.8)
  
  if(ShowAxis){
    rgl::material3d(col = "black")

    #if(!TextureRendering)
      rgl::persp3d(x=lines, y=columns, z=z, color=color, aspect=FALSE, lit=F,box=F, texmagfilter="nearest", texminfilter="nearest", texenvmap=TRUE)
   # else
    #  rgl::persp3d(x=lines, y=columns, z=z, col="white", texture="tmpGeneralizedUmatrix.png", textype="rgb",lit=F)
      
      if(isFALSE(Silent)){
        cat("plotToporaphicMaps: persp3d computed.\n")
      }
  }
  else{
    #if(!TextureRendering)
      rgl::surface3d(x=lines, y=columns, z=z, color=color, aspect=FALSE, lit=F)
   # else
     # rgl::surface3d(x=lines, y=columns, z=z, col="white", texture="tmpGeneralizedUmatrix.png", textype="rgb",lit=F)
    if(isFALSE(Silent)){
      cat("plotToporaphicMaps:surface 3d computed. \n")
    }
  }
   #rgl::bgplot3d( suppressWarnings ( fields::image.plot( legend.only=TRUE, legend.args=list(text='legend'), zlim=round(range(z,na.rm = T),1),col=color) )  )
   #rgl::legend3d("topleft", legend="a",col=color, pch=20)
  if(!ShowAxis){
    rgl::title3d(main = main,sub = sub,xlab = xlab,ylab = ylab,zlab = zlab)
  }else{
    rgl::title3d(main = main,sub = sub)
  }
##########################################################################################
## BestMatchingUnits darstellen ----
##########################################################################################
#Aus showUmatrix3d, package Umatrix   
  if(!is.null(BestMatchingUnits)){

    ColorClass = c()
    #for(i in 1:length(unique(Cls))) ColorClass[Cls == sort(unique(Cls))[i]] = i
    for(i in 1:length(unique(Cls))) ColorClass[Cls == unique(Cls)[i]] = i
    
    BestMatchingUnitsHeights <- sapply(1:nrow(BestMatchingUnits), function(x) z[BestMatchingUnits[x,2],BestMatchingUnits[x,3]])

    if(is.list(BestMatchingUnitsHeights))
      BestMatchingUnitsHeights = unlist(BestMatchingUnitsHeights)

    rgl::spheres3d(x=BestMatchingUnits[,2],y = BestMatchingUnits[,3], z = BestMatchingUnitsHeights, col = ClsColors[ColorClass], radius = BmSize)
  
    
    #ToDo: shape selection for multiple clusters (>20 clusters)
    #shapelist3d(shapes, x = 0, y = NULL, z = NULL, size = 1, matrix = NULL, override = TRUE, 
    #            ..., plot = TRUE)
    
    
    if(isFALSE(Silent)){
      cat("plotToporaphicMaps: Bestmatching Unit drawn.\n")
    }
  }

  # konturlinien zeichnen
 # if(!TextureRendering){
   if(RenderingContourLines){
    lines = contourLines(lines,columns,z, nlevels = NoLevels)
      #XYZ=list()
      X=c()
      Y=c()
      Z=c()
     for (i in seq_along(lines)) {
       x <- lines[[i]]$x
       y <- lines[[i]]$y
       z <- rep(lines[[i]]$level, length(x))
       if(isFALSE(Silent)){
         #cat(paste("plotToporaphicMaps: contour line",i," drawn.\n"))
       }
       #rgl::lines3d(x = x, y = y,z =z)
       #XYZ[[i]]=xyz.coords(x=x,y=y,z=z)
       X=c(X,NaN,x)
       Y=c(Y,NaN,y)
       Z=c(Z,NaN,z)
       
     }
      rgl::lines3d(x = X, y = Y,z =Z,alpha=0.75)
      #return(XYZ)
      #lapply(XYZ, function(x) rgl::lines3d(x))#lit=FALSE,point_antialias=FALSE,line_antialias = FALSE,fog=FALSE,depth_test="never",depth_mask=FALSE,texenvmap=FALSE,texmipmap=FALSE))
   }
   if(isFALSE(Silent)){
    cat(paste("plotToporaphicMaps: contour lines computed.\n"))
   }
   
 if(!is.null(Names)){
   if(is.character(Names)){
     if(length(Names)==length(unique(Cls))){
       rgl::legend3d(NamesPosition, legend = paste(Names), 
                     pch = NamesPch, col = NamesColors,
                     cex=NamesCex, 
                     inset=c(0.02),bty="n",title=NamesTitle)
   
     }else{
       warning('Names are not of length of number of clusters and hence ignored.')
     }
   }else{
     warning('Names are not a character vector and hence ignored.')
   }
 }
   

  widgets <- rgl::rglwidget()

  return(invisible(widgets))
}
