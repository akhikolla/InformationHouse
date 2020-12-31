showMatrix3D <- function(Matrix = NULL, BestMatches=NULL, Cls=NULL,Imx=NULL, Toroid=TRUE, HeightScale=NULL,
                          BmSize=0.5, RemoveOcean=T, ColorStyle = "Umatrix",
                          ShowAxis=F, SmoothSlope=F, ClsColors = NULL, FileName = NULL){
# showUmatrix3D(Umatrix, BestMatches, Cls, Toroid)
# Draws a plot of given Umatrix
# INPUT
# OPTIONAL
# Matrix(1:Lines,1:Columns)	      Matrix to be plotted
# BestMatches(1:n,1:2)		      Positions of bestMatches to be plotted onto the Umatrix
# Cls(1:n)			            class identifier for the bestmatch at the given point
# Toroid				                  should the Umatrix be drawn 4times?
# HeightScale                 how high should the mountains be?
# SmoothSlope                 try to increase the Island size, to get smooth slopes around the island

# OUTPUT
# ---
# author: FL
# 1.Editor: MT: Definition der Hoehen im 3Dplot diverse bugfixes
#  requireRpackage("rgl")
  #require(evd)

  if(is.null(ClsColors)) ClsColors = DefaultColorSequence()
	if(!requireNamespace("rgl", quietly = T)) stop("Package Rgl could not be loaded.")

  #########
  #########
  if(is.null(Matrix))
    stop("Matrix needs to be given")
 ## MT: Normalization der Matrix werte
  # Milligan, Copper 1988 A Study of Standadization of Variables in Cluster Analysis,
  # robust Normalization Z_5 :"Z_5 is bounded by 0.0 and 1.0 with at least one observed value at each of these end points"

  quants=quantile(as.vector(Matrix),c(0.01,0.5,0.99))
  minU=quants[1]
  maxU=quants[3]

  ### Hoehe aus Matrix schaetzen
 #Verhaeltnis zwischen minhoehe/maxHoehe=1/HeightScale
  # if(is.null(HeightScale))
  #   HeightScale=round(maxU/max(minU,0.03),0)

  Matrix=(Matrix-minU)/(maxU-minU)

  quants2=quantile(as.vector(Matrix),c(0.01,0.5,0.99))
  minU2=quants2[1]
  maxU2=quants2[3]

  ### Hoehe aus Matrix schaetzen
  #Verhaeltnis zwischen minhoehe/maxHoehe=1/HeightScale
  if(is.null(HeightScale)){
      HeightScale=round(maxU2/(2*max(minU2,0.05)),0) #Factor 2 damit die Matrix Hoehe nicht zu riesig wird
      # proportionen sind fuer standard umatrix der groesse 80x50 gedacht. fuer andere groessen linear skalieren
      stretchFactor = sqrt(nrow(Matrix)^2 + ncol(Matrix)^2) / sqrt(50^2 + 80^2)
  }
  #MT: Die Level muessen vor der Begrenzung der Werte auf 0 und 1 gesetz werden,
  #aber nachdem der Wertebereich umgeschoben wurde, damit die Dichte in den Grenzbereichen abgeschetzt werden kann

  # proportionen sind fuer standard umatrix der groesse 80x50 gedacht. fuer andere groessen linear skalieren
  stretchFactor = sqrt(nrow(Matrix)^2 + ncol(Matrix)^2) / sqrt(50^2 + 80^2)

  indMax=which(Matrix>1,arr.ind=T)
  indMin=which(Matrix<0,arr.ind=T)
  if(length(indMax)>0)
    Matrix[indMax]=1
  if(length(indMin)>0)
    Matrix[indMin]=0


  #######
  ## bestmatches vorbereiten
  #######
  if(!is.null(BestMatches)) if(is.null(Cls)) Cls = rep(1, nrow(BestMatches))
  BestMatches = CheckBestMatches(BestMatches, Cls, shiny=F)
  CheckUmatrix(Matrix, shiny=F)
  CheckImx(Imx, shiny=F)

  ########
  ## tiling durchf?hren
  #########
  if(Toroid){
    tU <- ToroidUmatrix(Matrix,BestMatches,Cls)
    Matrix <- tU$Umatrix
    BestMatches <- tU$BestMatches
    Cls <- tU$Cls
  }


  #########
  ## groessere imx bauen
  #########
  bigImx = Imx

  if(SmoothSlope){
    for(i in 1:10){
      if(!is.null(Imx)){
        tmpImx = bigImx
        for(x in 2:(nrow(Matrix)-1)){
          for(y in 2:(ncol(Matrix)-1)){
            if((Matrix[x-1,y] >= 0.3)&(tmpImx[x-1,y]==0)) bigImx[x,y] = 0
            if((Matrix[x+1,y] >= 0.3)&(tmpImx[x+1,y]==0)) bigImx[x,y] = 0
            if((Matrix[x,y-1] >= 0.3)&(tmpImx[x,y-1]==0)) bigImx[x,y] = 0
            if((Matrix[x,y+1] >= 0.3)&(tmpImx[x,y+1]==0)) bigImx[x,y] = 0
          }
        }
      }
    }
  }


  ##########
  ## an umatrix rumrechnen
  ##########
  zcol = cut(Matrix,128) # divide the Height into 128 different levels for coloring

  lines = seq(1, nrow(Matrix), len = nrow(Matrix))
  columns = seq(1, ncol(Matrix), len = ncol(Matrix))

  if(ColorStyle == "Umatrix") Colormap = c("#3C6DF0","#3C6DF0","#3C6DF0","#006602","#006A02","#006D01","#007101","#007501","#007901","#007C00","#008000","#068103","#118408","#0B8305","#17860A","#1D870D","#228810","#288A12","#2E8B15","#348D18","#398E1A","#3F8F1D","#45911F","#4A9222","#509325","#569527","#5C962A","#61982C","#67992F","#6D9A32","#729C34","#789D37","#7E9F39","#84A03C","#89A13F","#8FA341","#95A444","#9AA547","#A0A749","#A6A84C","#ACAA4E","#B1AB51","#B7AC54","#BDAE56","#C3AF59","#C8B15B","#CEB25E","#CBAF5C","#C8AC59","#C5A957","#C3A654","#C0A352","#BDA050","#BA9D4D","#B7994B","#B49648","#B29346","#AF9044","#AC8D41","#A98A3F","#A6873C","#A3843A","#A08138","#9E7E35","#9B7B33","#987830","#95752E","#92722B","#8F6E29","#8C6B27","#8A6824","#876522","#84621F","#815F1D","#7E5C1B","#7B5918","#795616","#765313","#714E0F","#6C480B","#674307","#6F4D15","#785822","#806230","#896D3E","#91774C","#998159","#A28C67","#AA9675","#B3A183","#BBAB90","#C3B59E","#CCC0AC","#D4CABA","#DDD5C7","#E5DFD5","#E7E1D8","#E9E4DB","#EBE6DE","#ECE8E1","#EEEAE4","#F0EDE7","#F2EFEA","#F4F1ED","#F6F4F0","#F8F6F3","#F9F8F6","#FBFAF9","#FDFDFC","#FFFFFF","#FFFFFF","#FEFEFE","#FEFEFE","#FEFEFE","#FDFDFD","#FDFDFD","#FDFDFD","#FCFCFC","#FCFCFC","#FCFCFC","#FBFBFB","#FBFBFB","#FBFBFB","#FAFAFA","#FAFAFA","#FAFAFA","#F9F9F9","#F9F9F9","#FFFFFF","#FFFFFF")
  else if(ColorStyle == "Pmatrix") Colormap = c("#FFFFFF","#FFFFF7","#FFFFEF","#FFFFE7","#FFFFDF","#FFFFD7","#FFFFCF","#FFFFC7","#FFFFBF","#FFFFB7","#FFFFAF","#FFFFA7","#FFFF9F","#FFFF97","#FFFF8F","#FFFF87","#FFFF80","#FFFF78","#FFFF70","#FFFF68","#FFFF60","#FFFF58","#FFFF50","#FFFF48","#FFFF40","#FFFF38","#FFFF30","#FFFF28","#FFFF20","#FFFF18","#FFFF10","#FFFF08","#FFFF00","#FFFA00","#FFF400","#FFEF00","#FFEA00","#FFE400","#FFDF00","#FFDA00","#FFD400","#FFCF00","#FFCA00","#FFC500","#FFBF00","#FFBA00","#FFB500","#FFAF00","#FFAA00","#FFA500","#FF9F00","#FF9A00","#FF9500","#FF8F00","#FF8A00","#FF8500","#FF8000","#FF7A00","#FF7500","#FF7000","#FF6A00","#FF6500","#FF6000","#FF5A00","#FF5500","#FF5000","#FF4A00","#FF4500","#FF4000","#FF3A00","#FF3500","#FF3000","#FF2B00","#FF2500","#FF2000","#FF1B00","#FF1500","#FF1000","#FF0B00","#FF0500","#FF0000","#FA0000","#F40000","#EF0000","#EA0000","#E40000","#DF0000","#DA0000","#D40000","#CF0000","#CA0000","#C50000","#BF0000","#BA0000","#B50000","#AF0000","#AA0000","#A50000","#9F0000","#9A0000","#950000","#8F0000","#8A0000","#850000","#800000","#7A0000","#750000","#700000","#6A0000","#650000","#600000","#5A0000","#550000","#500000","#4A0000","#450000","#400000","#3A0000","#350000","#300000","#2B0000","#250000","#200000","#1B0000","#150000","#100000","#0B0000","#050000")
  else stop("ColorStyle not found.")


  Nrlevels2 = 2*HeightScale*stretchFactor #MT Farbintervalle gleich 2*hoehenintervalle, siehe oben
  levelBreaks <- seq(0,1.000001,length.out=(Nrlevels2+1))



  ##########
  ## remove Heights by island
  ###########
  if(!is.null(Imx)){
    #Umatrix[which(Imx == 1)] = 0
    Matrix[which(bigImx == 1)] = 0

    if(!is.null(BestMatches)){
      BestMatchesFilter = rep(T,nrow(BestMatches)) # every Bestmatch stays
      for(i in 1:nrow(Imx)){
        for(j in 1:ncol(Imx)){
          if(Imx[i,j] == 1){
            if(!is.null(BestMatches))
              BestMatchesFilter[(BestMatches[,2] == i) & (BestMatches[,3] == j)] = F
          }
        }
      }
      BestMatches = BestMatches[BestMatchesFilter,]
      if(!is.null(Cls)) Cls = Cls[BestMatchesFilter]
    }
      #### remove ocean around Umatrix
    if(RemoveOcean){
      oceanLine = apply(Matrix, 1, function(x) all(x==0))
      startLine = min(which(!oceanLine),na.rm=T)
      endLine = length(oceanLine) - min(which(rev(!oceanLine)),na.rm=T) + 1

      oceanCol = apply(Matrix, 2, function(x) all(x==0))
      startCol = min(which(!oceanCol),na.rm=T)
      endCol = length(oceanCol) - min(which(rev(!oceanCol)),na.rm=T) + 1

      if(!is.null(BestMatches)){
        BestMatches <- BestMatches - cbind(rep(0,nrow(BestMatches)),startLine-1,startCol-1)
      }
      Matrix <- Matrix[startLine:endLine,startCol:endCol]
      Imx <- Imx[startLine:endLine,startCol:endCol]
      bigImx <- bigImx[startLine:endLine,startCol:endCol]
    }
  }

  # change Matrix into 35 fixed values
  splittedMatrix = Matrix
  for(i in 1:Nrlevels2){
    splittedMatrix[ (Matrix >= levelBreaks[i]) & (Matrix <= levelBreaks[i+1]) ] = levelBreaks[i]
  }
  splittedMatrix=(floor(splittedMatrix * length(Colormap)))+1
  color = Colormap[splittedMatrix]
  #if(!is.null(Imx)) color[which(Imx == 1)] = NA
  if(!is.null(Imx)) color[which(bigImx == 1)] = NA

  z<-Matrix*HeightScale*stretchFactor #MT: Die Hoehe entspricht den Umatrixwerten, siehe vorne Proportion

  lines = seq(1, nrow(z), len = nrow(z))
  columns = seq(1, ncol(z), len = ncol(z))



  #########
  ## umatrix darstellen
  #########
  rgl::open3d()
  if(ShowAxis){
    rgl::material3d(col = "black")
    rgl::persp3d(x=lines, y=columns, z=z, color=color, aspect=FALSE, lit=F,box=F, texmagfilter="nearest", texminfilter="nearest", texenvmap=TRUE)
  }
  else{
    rgl::surface3d(x=lines, y=columns, z=z, color=color, aspect=FALSE, lit=F)
  }
  
  if(!is.null(FileName)){
    rgl::writeSTL(FileName)
  }
  
  ###########
  ## bestmatches darstellen
  ###########
  if(!is.null(BestMatches)){
    if(length(ClsColors) < length(DefaultColorSequence()))
      ClsColors = c(ClsColors, DefaultColorSequence()[(length(ClsColors)+1):(length(DefaultColorSequence()))])

    ColorClass = c()
    for(i in 1:length(unique(Cls))) ColorClass[Cls == sort(unique(Cls))[i]] = i

    BestMatchesHeights <- sapply(1:nrow(BestMatches), function(x) z[BestMatches[x,2],BestMatches[x,3]])

    if(is.list(BestMatchesHeights))
      BestMatchesHeights = unlist(BestMatchesHeights)

    bmuIds = rgl::spheres3d(x=BestMatches[,2],BestMatches[,3], BestMatchesHeights, col = ClsColors[ColorClass], radius = BmSize)
  }


  #############
  ## konturlinien zeichnen
  ##############
  lines = contourLines(lines,columns,z, nlevels = 35)
   for (i in seq_along(lines)) {
     x <- lines[[i]]$x
     y <- lines[[i]]$y
     z <- rep(lines[[i]]$level, length(x))
     rgl::lines3d(x, y, z)
   }
}
