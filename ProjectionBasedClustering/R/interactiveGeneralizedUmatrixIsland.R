interactiveGeneralizedUmatrixIsland <- function(Umatrix, Bestmatches=NULL, Cls=NULL){
  # Imx = interactiveGeneralizedUmatrixIsland(Umatrix, Bestmatches, Cls)
  #
  # INPUT
  # Umatrix
  # Bestmatches
  # Cls
  # OUTPUT
  # island          the generated Imx
  # Author: FL, MT

  #requireRpackage("shiny")
  #requireRpackage("shinyjs")
  #requireRpackage("png")
  #requireRpackage("tcltk")
  #FilePath<-getwd() #MT: sonst funktionierts nicht
  DataName = NULL
  UmatrixSize = 1
  requireNamespace("png")
  QuitButtonText = "Quit"
  #######################################################################
  ## Hilfs Funktion
  ######################################################################
  UmxBestmatchesFromIsland <- function(Umatrix, Bestmatches, Imx, Cls= NULL, Toroid=T , RemoveOcean = T){
    # V <- UmxBestmatchesFromIsland(Umx, Bestmatches, Imx, Cls)
    # Cuts out the Imx (Island) from the Umatrix, Bestmatches and Cls
    # INPUT
    # Umatrix
    # Bestmatches
    # Imx
    # Optional
    # Cls
    # Toroid         should the umatrix be Toroid? (was it calculated on a toroid)
    # RemoveOcean   should the umatrix shortened to the actual island?
    # OUTPUT
    # list of
    # Umatrix
    # Bestmatches
    # Cls
    
    # auf vierfache groesse kacheln
    if(Toroid){
      tU <- ToroidUmatrix(Umatrix, Bestmatches, Cls)
      Umatrix <- tU$Umatrix
      Bestmatches <- tU$Bestmatches
      Cls <- tU$Cls
    }
    
    #Bestmatches = CheckBestmatches(Bestmatches, Cls)
    
    # configure filter, so that every bestmatch stays in
    if(!is.null(Bestmatches)){
      BestmatchesFilter = rep(T,nrow(Bestmatches)) # every Bestmatch stays
    }
    
    # put Imx on Umatrix and Bestmatches if given
    if(!is.null(Imx)){
      for(i in 1:nrow(Imx)){
        for(j in 1:ncol(Imx)){
          if(Imx[i,j] == 1){
            Umatrix[i,j] = NA
            if(!is.null(Bestmatches))
              BestmatchesFilter[(Bestmatches[,2] == i) & (Bestmatches[,3] == j)] = F
          }
        }
      }
      
      if(!is.null(Bestmatches)) Bestmatches = Bestmatches[BestmatchesFilter,]
      if((!is.null(Cls)) & (!is.null(Bestmatches))) Cls = Cls[BestmatchesFilter]
    }
    
    #### remove ocean around Umatrix
    if(RemoveOcean){
      oceanLine = !apply(Umatrix, 1, function(x) any(x != -1))
      startLine = min(which(!oceanLine),na.rm=T)
      endLine = length(oceanLine) - min(which(rev(!oceanLine)),na.rm=T) + 1
      
      oceanCol = !apply(Umatrix, 2, function(x) any(x != -1))
      startCol = min(which(!oceanCol),na.rm=T)
      endCol = length(oceanCol) - min(which(rev(!oceanCol)),na.rm=T) + 1
      
      if(!is.null(Bestmatches)){
        Bestmatches <- Bestmatches - cbind(rep(0,nrow(Bestmatches)),startLine-1,startCol-1)
      }
      Umatrix <- Umatrix[startLine:endLine,startCol:endCol]
    }
    
    Umatrix[which(is.na(Umatrix))] = 0
    
    return(list(Umatrix = Umatrix, Bestmatches = Bestmatches, Cls=Cls))
  }
  
  bestUmatrixTranslation <- function(Umatrix, Bestmatches=NULL){
    # rotation = bestUmatrixTranslation(Umatrix, Bestmatches)
    # Get the ideal amount of lines and columns to rotate the Umatrix by, so that the borders
    # are on positions with the greatest Heights and fewest Bestmatches
    # INPUT
    # Umatrix(1:Lines, 1:Columns)	  Umatrix
    # OPTIONAL
    # Bestmatches(1:n,1:2)		  Bestmatches
    # OUTPUT
    # list with
    # lines				  number of lines to rotate up
    # cols				  number of columns to rotate left
    # author: Florian Lerch, Michael Thrun
    
    ncols = ncol(Umatrix)
    nrows = nrow(Umatrix)
    
    # calculate Height for every line and column
    lineHeights = rowSums(Umatrix)
    colHeights = colSums(Umatrix) 
    
    if(!is.null(Bestmatches)){
      # add bm keys if not given
      if(ncol(Bestmatches) == 2) Bestmatches <- cbind(1:nrow(Bestmatches),Bestmatches)
      
      # find intruding Bestmatches on every line and column
      lineIntruders = rep(1,nrows)
      colIntruders = rep(1,ncols)
      for(i in 1:nrows) lineIntruders[i] = length(which(Bestmatches[,2] == i))
      for(i in 1:ncols) colIntruders[i] = length(which(Bestmatches[,3] == i))
      
      lineIntruders2 = rep(1,nrows)
      colIntruders2 = rep(1,ncols)
      for(i in 2:(nrows-1)) lineIntruders2[i] = lineIntruders[(i-1)]+lineIntruders[i]+lineIntruders[(i+1)]
      
      for(i in 2:(ncols-1)) colIntruders2[i] = colIntruders[(i-1)]+colIntruders[i]+colIntruders[(i+1)]
      
      lineIntruders2[1]=lineIntruders[nrows]+lineIntruders[1]+lineIntruders[2]
      lineIntruders2[nrows]=lineIntruders[nrows]+lineIntruders[1]+lineIntruders[nrows-1]
      
      colIntruders2[1]=colIntruders[ncols]+colIntruders[1]+colIntruders[2]
      colIntruders2[ncols]=colIntruders[ncols]+colIntruders[1]+colIntruders[ncols-1]
      
      colIntruders2[colIntruders2==0]=1
      lineIntruders2[lineIntruders2==0]=1
      
      # norm Heights with intruding Bestmatches
      lineHeights = lineHeights / lineIntruders2
      colHeights = colHeights / colIntruders2
    }
    
    # calculate the upper left cutting point
    cutLine = which.max(lineHeights)
    cutCol = which.max(colHeights)
    
    list(lines = cutLine, cols = cutCol)
  }
  createMaskFromPolygon <- function(Height, Width, CutoutPol){
    # Sets every value of the Mask to 0, that is not within CutoutPol or else to 1
    # INPUT
    # Height                            Height of umatrix
    # Width                             Width of umatrix
    # CutoutPol                         dataframe with x and y
    # OUTPUT
    # UmatrixMask[1:lines, 1:cols]      Umatrix Mask
    # author: FL
    
    #library(fields)
    UmatrixMask = matrix(1,nrow=Height, ncol=Width)
    requireNamespace('fields')
    # cut out polygon out of umatrix
    for(i in 1:ncol(UmatrixMask)){
      umatrixRow = cbind(rep(i,nrow(UmatrixMask)), 1:nrow(UmatrixMask))
      umatrixRowFilter = fields::in.poly(umatrixRow, CutoutPol)
      UmatrixMask[(umatrixRowFilter), i] = 0  # set everything to 1 that is within the polygon
    }
    
    UmatrixMask
  }
  
  ToroidUmatrix <- function(Umatrix, Bestmatches = NULL, Cls = NULL){
    # Umatrix <- ToroidUmatrix(Umatrix)
    # make 4 umatrices out of one. Optionally do the same for the Bestmatches and their classes
    # INPUT
    # Umatrix(1:Lines,1:Columns)	      Umatrix to be Toroid
    # OPTIONAL
    # Bestmatches(1:n,1:2)		      Positions of Bestmatches to be plotted onto the Umatrix
    # Cls(1:n)			      class identifier for the bestmatch at the given point
    # OUTPUT
    # list with
    # Umatrix(1:(2*Lines),1:(2*Columns))
    # Bestmatches(1:(4*n),1:2)
    # Cls(1:(4*n)))
    # author: copied out of plotMatrixTopView; FL
    rows=nrow(Umatrix)
    cols=ncol(Umatrix)
    
    Umatrix <- Umatrix[c(1:rows,1:rows),c(1:cols,1:cols)]
    
    bm_keys = NULL
    if(!is.null(Bestmatches)){
      # extract the keys
      if(ncol(Bestmatches) == 3){
        bm_keys = Bestmatches[,1]
        Bestmatches = Bestmatches[,c(2,3)]
      }
      else{
        bm_keys = 1:nrow(Bestmatches)
      }
      
      bmRow <- nrow(Bestmatches)
      Bestmatches <- Bestmatches[rep(1:bmRow,4),]
      Bestmatches[(bmRow+1):(2*bmRow),1] <- Bestmatches[(bmRow+1):(2*bmRow),1]+rows # unterer rechter Quadrant
      Bestmatches[(2*bmRow+1):(3*bmRow),2] <- Bestmatches[(2*bmRow+1):(3*bmRow),2]+cols # oberer linker Quadrant
      Bestmatches[(3*bmRow+1):(4*bmRow),1] <- Bestmatches[(3*bmRow+1):(4*bmRow),1]+rows # oberer rechter Quadrant
      Bestmatches[(3*bmRow+1):(4*bmRow),2] <- Bestmatches[(3*bmRow+1):(4*bmRow),2]+cols # oberer rechter Quadrant
    }
    # If Cls not missing, adjust
    if(!is.null(Cls)){
      Cls <- rep(Cls,4)
    }
    # reattach the keys
    if(!is.null(bm_keys))
      Bestmatches = cbind(rep(bm_keys,4), Bestmatches)
    
    list(Umatrix=Umatrix,Bestmatches=Bestmatches,Cls=Cls)
  } #End ToroidUmatrix
  #########################################################################################################
  
  
  
  #ShowUmatrix_Hlp
  showUmatrix_hlp <- function(Matrix = NULL, Bestmatches=NULL, Cls = NULL, ClsColors=NULL,
                              ColorStyle = "Umatrix", Toroid=TRUE,BmSize=2, DrawLegend=F,
                              FixedRatio=T, CutoutPol=NULL, Nrlevels = NULL, TransparentContours = T,
                              Imx = NULL, Clean=F, RemoveOcean=F, TransparentOcean = FALSE,
                              Title = NULL, BestmatchesLabels = NULL,
                              BestmatchesLabelStyle = c(), BestmatchesShape=19, MarkDuplicatedBestmatches = F, YellowCircle = F){
    
    
    if(is.null(ClsColors)) ClsColors=GeneralizedUmatrix::DefaultColorSequence
    if(ColorStyle == "Umatrix") Colormap = c("#3C6DF0","#3C6DF0","#3C6DF0","#006602","#006A02","#006D01","#007101","#007501","#007901","#007C00","#008000","#068103","#118408","#0B8305","#17860A","#1D870D","#228810","#288A12","#2E8B15","#348D18","#398E1A","#3F8F1D","#45911F","#4A9222","#509325","#569527","#5C962A","#61982C","#67992F","#6D9A32","#729C34","#789D37","#7E9F39","#84A03C","#89A13F","#8FA341","#95A444","#9AA547","#A0A749","#A6A84C","#ACAA4E","#B1AB51","#B7AC54","#BDAE56","#C3AF59","#C8B15B","#CEB25E","#CBAF5C","#C8AC59","#C5A957","#C3A654","#C0A352","#BDA050","#BA9D4D","#B7994B","#B49648","#B29346","#AF9044","#AC8D41","#A98A3F","#A6873C","#A3843A","#A08138","#9E7E35","#9B7B33","#987830","#95752E","#92722B","#8F6E29","#8C6B27","#8A6824","#876522","#84621F","#815F1D","#7E5C1B","#7B5918","#795616","#765313","#714E0F","#6C480B","#674307","#6F4D15","#785822","#806230","#896D3E","#91774C","#998159","#A28C67","#AA9675","#B3A183","#BBAB90","#C3B59E","#CCC0AC","#D4CABA","#DDD5C7","#E5DFD5","#E7E1D8","#E9E4DB","#EBE6DE","#ECE8E1","#EEEAE4","#F0EDE7","#F2EFEA","#F4F1ED","#F6F4F0","#F8F6F3","#F9F8F6","#FBFAF9","#FDFDFC","#FFFFFF","#FFFFFF","#FEFEFE","#FEFEFE","#FEFEFE","#FDFDFD","#FDFDFD","#FDFDFD","#FCFCFC","#FCFCFC","#FCFCFC","#FBFBFB","#FBFBFB","#FBFBFB","#FAFAFA","#FAFAFA","#FAFAFA","#F9F9F9","#F9F9F9","#FFFFFF","#FFFFFF")
    else if(ColorStyle == "Pmatrix") Colormap = c("#FFFFFF","#FFFFF7","#FFFFEF","#FFFFE7","#FFFFDF","#FFFFD7","#FFFFCF","#FFFFC7","#FFFFBF","#FFFFB7","#FFFFAF","#FFFFA7","#FFFF9F","#FFFF97","#FFFF8F","#FFFF87","#FFFF80","#FFFF78","#FFFF70","#FFFF68","#FFFF60","#FFFF58","#FFFF50","#FFFF48","#FFFF40","#FFFF38","#FFFF30","#FFFF28","#FFFF20","#FFFF18","#FFFF10","#FFFF08","#FFFF00","#FFFA00","#FFF400","#FFEF00","#FFEA00","#FFE400","#FFDF00","#FFDA00","#FFD400","#FFCF00","#FFCA00","#FFC500","#FFBF00","#FFBA00","#FFB500","#FFAF00","#FFAA00","#FFA500","#FF9F00","#FF9A00","#FF9500","#FF8F00","#FF8A00","#FF8500","#FF8000","#FF7A00","#FF7500","#FF7000","#FF6A00","#FF6500","#FF6000","#FF5A00","#FF5500","#FF5000","#FF4A00","#FF4500","#FF4000","#FF3A00","#FF3500","#FF3000","#FF2B00","#FF2500","#FF2000","#FF1B00","#FF1500","#FF1000","#FF0B00","#FF0500","#FF0000","#FA0000","#F40000","#EF0000","#EA0000","#E40000","#DF0000","#DA0000","#D40000","#CF0000","#CA0000","#C50000","#BF0000","#BA0000","#B50000","#AF0000","#AA0000","#A50000","#9F0000","#9A0000","#950000","#8F0000","#8A0000","#850000","#800000","#7A0000","#750000","#700000","#6A0000","#650000","#600000","#5A0000","#550000","#500000","#4A0000","#450000","#400000","#3A0000","#350000","#300000","#2B0000","#250000","#200000","#1B0000","#150000","#100000","#0B0000","#050000")
    else stop("ColorStyle not found.")
    
    
    # v <- showUmatrix_hlp()
    # Draws a plot of given Umatrix
    # INPUT
    # OPTIONAL
    # Umatrix(1:Lines,1:Columns)	      Umatrix to be plotted
    # Bestmatches(1:n,1:2)		      Positions of Bestmatches to be plotted onto the Umatrix
    # Cls(1:n)			        Class identifier for the bestmatch at the given point
    # ClsColors(1:m)			Vector of colors that will be used to colorize the different classes
    # Colormap(1:o)				Vector of colors that will be used to colorize the Heights of the Umatrix
    # BmSize                Integer between 0.1 and 5, magnification factor of the points,
    #                       if Bestmatches given
    # DrawLegend
    # FixedRatio            should the ratio of Width and Height be fixed or matched to the window Width and Height
    # CutoutPol(1:n,1:2)     only draws the area within given polygon
    # Nrlevels               nr of breaks that will be done
    # TransparentContours    use half transparent contours. Looks better but is slow
    # Imx                    a Imx (Imx) that will be used to cut out the umatrix
    # Fast                  Faster version that will be drawn using baseplot
    # RemoveOcean             restrict the umatrix to the island, and reduce the ocean to a minimum
    # TransparentOcean        draw the water surrounding the island transparent instead of blue
    # Title                   a title to be printed above the umatrix
    # BestmatchesLabels(1:n)         a vector of strings matching to the Bestmatches, that will be shown above the Bestmatches
    # BestmatchesLabelStyle           a list containing possible modifications to the style of the labels
    # BestmatchesShape                number of the R-Shape to be used for drawing BMU
    # MarkDuplicatedBestmatches
    # YellowCircle            draw a circle around Bestmatches to get a better contrast
    
    # OUTPUT
    # ---()
    # author: Florian Lerch
    # 1.Editor: MT 09/2015 BmSize added, NumberContours added, region=T set
    # 2.Editor: MT 04/2016: Umatrix-Normierung, Definition der Hoehen im topview
    
    #library(reshape2) # necessary for "melt" which converts the Umatrix to a data.frame
    #library(ggplot2)
    #library(fields)
    requireNamespace("reshape2")
    # in einer frueheren Variante wurden Bestmatches mit der groesse 2 wesentlich kleiner gezeichnet. diese proportionen
    # werden hiermit aufrechterhalten.
    BmSize = BmSize / 2
    
    ################
    #### bmuLabel Standard Style ----
    ################
    if(is.null(BestmatchesLabelStyle$color)) BestmatchesLabelStyle$color = "white"
    if(is.null(BestmatchesLabelStyle$fontface)) BestmatchesLabelStyle$fontface = "bold"
    if(is.null(BestmatchesLabelStyle$fill)) BestmatchesLabelStyle$fill = "red"
    if(is.null(BestmatchesLabelStyle$alpha)) BestmatchesLabelStyle$alpha = 0.5
    if(is.null(BestmatchesLabelStyle$angle)) BestmatchesLabelStyle$angle = 0
    if(is.null(BestmatchesLabelStyle$size)) BestmatchesLabelStyle$size = 4
    if(is.null(BestmatchesLabelStyle$nudge_x)) BestmatchesLabelStyle$nudge_x = 0
    if(is.null(BestmatchesLabelStyle$nudge_y)) BestmatchesLabelStyle$nudge_y = 2
    
    if(!Toroid){
      Imx = NULL
      RemoveOcean = FALSE
      TransparentOcean = FALSE
    }
    
    
    if(is.null(Matrix))
      stop("Matrix needs to be given")
    
    if(!is.matrix(Matrix)){
      stop("Matrix has to be of type matrix")
    }
    if(length(dim(Matrix))!=2)
      stop("Matrix has to be a of type matrix, not an array")
    
    if(!is.null(Cls)){
      if(length(ClsColors) < length(GeneralizedUmatrix::DefaultColorSequence))
        ClsColors = c(ClsColors, GeneralizedUmatrix::DefaultColorSequence[(length(ClsColors)+1):(length(GeneralizedUmatrix::DefaultColorSequence))])
      if(max(Cls) > length(ClsColors))
        stop(paste("The amount of given Colors (ClsColors) is not enough for the Cls. The highest Cls Value is",
                   max(Cls),"while the number of ClsColors only reaches to",
                   length(ClsColors)))
    }
    
    ## MT: Normalization der Umatrix werte
    # Milligan, Copper 1988 A Study of Standadization of Variables in Cluster Analysis,
    # robust Normalization Z_5 :"Z_5 is bounded by 0.0 and 1.0 with at least one observed value at each of these end points"
    quants=quantile(as.vector(Matrix),c(0.01,0.5,0.99),na.rm = T)
    # minU=min(Umatrix,na.rm=T)
    # maxU=max(Umatrix,na.rm=T)
    minU=quants[1]
    maxU=quants[3]
    #Verhaeltnis zwischen minhoehe/maxHoehe=1/HeightScale
    
    Matrix=(Matrix-minU)/(maxU-minU)
    
    quants2=quantile(as.vector(Matrix),c(0.01,0.5,0.99),na.rm = T)
    minU2=quants2[1]
    maxU2=quants2[3]
    
    ### Hoehe aus Umatrix schaetzen
    #Verhaeltnis zwischen minhoehe/maxHoehe=1/HeightScale
    if(is.null(Nrlevels)){
      Nrlevels=round(maxU2/max(minU2,0.05),0)
    }
    #MT: Die Level muessen vor der Begrenzung der Werte auf 0 und 1 gesetz werden,
    #aber nachdem der Wertebereich umgeschoben wurde, damit die Dichte in den Grenzbereichen abgeschetzt werden kann
    indMax=which(Matrix>1,arr.ind=T)
    indMin=which(Matrix<0,arr.ind=T)
    
    if(length(indMax)>0)
      Matrix[indMax]=1
    if(length(indMin)>0)
      Matrix[indMin]=0
    
    
    
    # force rownames. this keeps log of the "original position" of the bestmatch within the list.
    # that is necessary because of the way BestmatchesLabels is currently defined.
    if(!is.null(Bestmatches)){
      if(is.null(rownames(Bestmatches))) rownames(Bestmatches) <- 1:nrow(Bestmatches)
    }
    
    # verknuepfe BestmatchesLabels mit BestmatchKey
    if((!is.null(Bestmatches)) & (!is.null(BestmatchesLabels))){
      BestmatchesLabels <- cbind(Bestmatches[,1], BestmatchesLabels)
    }
    
    Nrlevels2=Nrlevels #nicht konturen sondern farbintervalle!
    #lohnt sich nicht eine andere zahl als die die konturenanzahl zu setzen
    # get Toroid values for the Umatrix
    if(Toroid){
      tU <- ToroidUmatrix(Matrix, Bestmatches, Cls)
      Matrix <- tU$Umatrix
      Bestmatches <- tU$Bestmatches
      #ToDo: Wenn ich hier erweitere muss ich auch insel entsprechend verschieben und imx ausgabe entsprechend zurueckverschieben
      #V=ExtendToroidalUmatrix(Matrix,Bestmatches[,2:3],10)
      #Matrix <- V$Umatrix
      #Bestmatches[,2:3] <- V$Bestmatches
      
      Cls <- tU$Cls
    }
    
    # configure filter, so that every bestmatch stays in
    if(!is.null(Bestmatches)){
      BestmatchesFilter = rep(T,nrow(Bestmatches)) # every Bestmatch stays
    }
    
    # put Imx on Umatrix and Bestmatches if given
    if(!is.null(Imx)){
      for(i in 1:nrow(Imx)){
        for(j in 1:ncol(Imx)){
          if(Imx[i,j] == 1){
            Matrix[i,j] = NA
            if(!is.null(Bestmatches))
              BestmatchesFilter[(Bestmatches[,2] == i) & (Bestmatches[,3] == j)] = F
          }
        }
      }
      
      if(!is.null(Bestmatches)) Bestmatches = Bestmatches[BestmatchesFilter,]
      if((!is.null(Cls)) & (!is.null(Bestmatches))) Cls = Cls[BestmatchesFilter]
    }
    
    #### remove ocean around Umatrix
    if(RemoveOcean){
      oceanLine = !apply(Matrix, 1, function(x) any(x != -1))
      startLine = min(which(!oceanLine),na.rm=T)
      endLine = length(oceanLine) - min(which(rev(!oceanLine)),na.rm=T) + 1
      
      oceanCol = !apply(Matrix, 2, function(x) any(x != -1))
      startCol = min(which(!oceanCol),na.rm=T)
      endCol = length(oceanCol) - min(which(rev(!oceanCol)),na.rm=T) + 1
      
      if(!is.null(Bestmatches)){
        Bestmatches <- Bestmatches - cbind(rep(0,nrow(Bestmatches)),startLine-1,startCol-1)
      }
      Matrix <- Matrix[startLine:endLine,startCol:endCol]
    }
    
    if(!TransparentOcean) Matrix[which(is.na(Matrix))] = 0
    
    nrows = nrow(Matrix)
    ncols = ncol(Matrix)
    
    colorx = Colormap
    
    # define level intervals
    # 35 levels
    # upper limit a bit above 1 so that every point inside the Umatrix is below it
    levelBreaks <- seq(0,1.000001,length.out=(Nrlevels2+1))
    
    # change Umatrix into 35 fixed values
    unfixedMatrix <- Matrix
    for(i in 1:Nrlevels2){
      Matrix[ (Matrix >= levelBreaks[i]) & (Matrix <= levelBreaks[i+1]) ] = levelBreaks[i]
    }
    
    
    
    duplicatedBestmatches = matrix(,ncol=2, nrow=0)
    # find duplicated Bestmatches
    if(!is.null(Bestmatches)){
      dup = duplicated(Bestmatches[,1], fromLast = T) | duplicated(Bestmatches[,1])
      duplicatedBestmatches = Bestmatches[dup,,drop=F]
    }
    
    # transparency of the contours
    if(TransparentContours) alpha = 0.5
    else alpha = 1
    
    ###############
    # draw the Umatrix with ggplot
    ###############
    #### convert Umatrix into a dataframe
    unnamedMatrix <- Matrix
    colnames(unnamedMatrix) <- NULL
    rownames(unnamedMatrix) <- NULL
    
    dfMatrix <- reshape2::melt(unnamedMatrix)
    colnames(dfMatrix) <- c("y","x","z")
    
    dfMatrix2 = dfMatrix
    dfMatrix2[is.na(dfMatrix$z),3]=0
    
    MatrixPlot <- ggplot(dfMatrix, aes_string('x','y')) +   # data layer
      #stat_contour(geom="polygon", aes(fill=..level..,z=z, colour=..level..),binWidth=0.05)+
      geom_tile(aes_string(fill = 'z'))+
      
      scale_fill_gradientn(colours=colorx,space='Lab', na.value="transparent")+
      stat_contour(aes_string(z='z'), data=dfMatrix2,bins=Nrlevels,size = 0.25,color='black',alpha=alpha,na.rm = F)+ #80% der Laufzeit wird fuer konturen benoetigt
      ylab("Lines (y)")+xlab("Columns (x)")
    
    #MatrixPlot = MatrixPlot + scale_y_reverse()
    
    
    if(!DrawLegend) MatrixPlot <- MatrixPlot + theme(legend.position="none")
    
    if(FixedRatio)
      MatrixPlot <- MatrixPlot + coord_fixed(1,xlim=c(0.5,ncols+0.5),ylim=c(0.5,nrows+0.5), expand = 0) + scale_y_reverse()
    else
      MatrixPlot <- MatrixPlot + coord_cartesian(xlim=c(0.5,ncols+0.5),ylim=c(nrows+0.5,0.5), expand = 0) + scale_y_reverse()
    
    if(Clean){
      if(is.null(Title)){
        MatrixPlot <- MatrixPlot + theme(
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(0.01,"lines"),
          #axis.text.margin = unit(0,"null"),
          panel.grid = element_blank(),
          axis.title.x=element_blank(),
          plot.margin = unit( c(0.001,0.001,0.001,0.001) , "lines" ),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          
          panel.background = element_blank(),
          #plot.background = element_blank(),
          legend.background = element_blank(),
          panel.border = element_blank(),
          plot.background=element_rect(fill="white", colour=NA)
          
        ) + labs('x'=NULL, 'y'=NULL)
        
        #if(packageVersion("ggplot2") >= 2.2.1) UmatrixPlot <- UmatrixPlot + theme(panel.spacing = unit(0,"null"))
        MatrixPlot <- MatrixPlot +theme(legend.spacing = unit(0,"null"))
      }
      else{
        MatrixPlot <- MatrixPlot + theme(
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.ticks.length = unit(0.01,"lines"),
          #axis.text.margin = unit(0,"null"),
          #panel.grid = element_blank(),
          #legend.spacing = unit(0,"null"),
          axis.title.x=element_blank(),
          #plot.margin = unit( c(0.001,0.001,0.001,0.001) , "lines" ),
          axis.title.y=element_blank(),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          
          panel.background = element_blank(),
          #plot.background = element_blank(),
          legend.background = element_blank(),
          panel.border = element_blank(),
          plot.background=element_rect(fill="white", colour=NA)
          
        ) + labs('x'=NULL, 'y'=NULL)
      }
    }
    
    # add the Bestmatches
    if(!is.null(Bestmatches)){
      if(is.null(Cls)) Cls <- rep(1,nrow(Bestmatches))
      
      Cls = factor(Cls)
      #classes <- sort(na.last=T,c(unique(Cls), setdiff(1:100,unique(Cls))))
      
      uniqueClasses <- sort(na.last=T,unique(Cls))
      
      ClsColors = c("darkgreen", ClsColors)
      names(ClsColors) = 0:(length(ClsColors)-1)
      
      MatrixPlot = MatrixPlot+scale_color_manual(values=ClsColors, name="Clusters")
      #UmatrixPlot = UmatrixPlot+scale_color_manual(c("0"="blue", "1"="green", "2"="red", "3"="yellow", "4"="brown", "5"="black",
      #                                               "6"="magenta", "7" = "white"), name="Clusters")
      
      if(!is.null(BestmatchesShape))
        if(length(BestmatchesShape)>1)
          MatrixPlot = MatrixPlot+scale_shape_manual(values=BestmatchesShape, name="Clusters")
      
      # extract Bestmatches without a class
      d = data.frame(y = Bestmatches[,2], x = Bestmatches[,3], class = Cls)
      #d1 = data.frame(y = Bestmatches[Cls!=0,2], x = Bestmatches[Cls!=0,3], class = Cls[Cls!=0])
      #d2 = data.frame(y = Bestmatches[Cls==0,2], x = Bestmatches[Cls==0,3], class = rep(length(ClsColors),sum(Cls==0)))
      
      if(YellowCircle)
        MatrixPlot <- MatrixPlot + geom_point(data=d, col = "yellow", size = rel(BmSize)*1.4,
                                              shape=BestmatchesShape)
      
      # elemente mit und ohne klasse
      if(length(BestmatchesShape)==1){
        MatrixPlot <- MatrixPlot + geom_point(aes_string(y = 'y', x='x', col='factor(class)'), data=d , size = rel(BmSize),
                                              shape=BestmatchesShape)
        # UmatrixPlot <- UmatrixPlot + geom_point(aes_string(y = 'y', x='x', col='factor(class)'), data=d2 , size = rel(BmSize),
        #                                       shape=BestmatchesShape)
      }
      else{
        MatrixPlot <- MatrixPlot + geom_point(aes_string(y = 'y', x='x', col='factor(class)', shape='factor(class)'), data=d ,
                                              stroke=BmSize)#lwd=3)#size = rel(BmSize), lwd=3)
        # UmatrixPlot <- UmatrixPlot + geom_point(aes_string(y = 'y', x='x', col='factor(class)', shape='factor(class)'), data=d2,
        #                                         stroke=BmSize)#lwd=3)#size = rel(BmSize), lwd=3)
      }
      
      
      
      # add labels
      if(!is.null(BestmatchesLabels)){
        bestmatchLabels = data.frame(y = Bestmatches[,2], x = Bestmatches[,3], 'label' = BestmatchesLabels[as.numeric(rownames(Bestmatches)),2])
        MatrixPlot <- MatrixPlot + geom_label(data=bestmatchLabels,
                                              aes_string(label="label"),
                                              color=BestmatchesLabelStyle$color,
                                              fontface=BestmatchesLabelStyle$fontface,
                                              fill=BestmatchesLabelStyle$fill ,
                                              alpha=BestmatchesLabelStyle$alpha ,
                                              angle=BestmatchesLabelStyle$angle ,
                                              size=BestmatchesLabelStyle$size ,
                                              nudge_x =BestmatchesLabelStyle$nudge_x,
                                              nudge_y =BestmatchesLabelStyle$nudge_y )
        #UmatrixPlot <- UmatrixPlot + geom_text(data=bestmatchLabels, aes(label=label))
      }
      
      # mark duplicates
      if(nrow(duplicatedBestmatches)>0){
        if(MarkDuplicatedBestmatches){
          duplicatedBestmatchesMarker = data.frame(y = duplicatedBestmatches[,2]+1, x = duplicatedBestmatches[,3])
          MatrixPlot <- MatrixPlot + geom_point(data = duplicatedBestmatchesMarker,
                                                shape=24,
                                                color="red",
                                                fill="red",
                                                size=4)
        }
      }
      
    }
    
    return(MatrixPlot + ggtitle(Title) + scale_size_area())  # return the plot

  }# End ShowUmatrix_hlp
  #########################################################################################################
  
  
  ############
  ### Fehler abfangen ----
  ############
  if(missing(Umatrix))
    stop("Missing Umatrix")
  if(!is.null(Cls)){
    if(is.list(Cls)) stop("Cls is not a vector")
    if(nrow(Bestmatches)!=length(Cls))
      stop('The length of the given classification does not equal the row length of the Bestmatches')
  }
  # umatrix normieren
  Umatrix = Umatrix/max(Umatrix)

  if(!is.null(Bestmatches))
    if(is.null(Cls)) Cls=rep(1,nrow(Bestmatches))

  Imx = matrix(0,nrow=nrow(Umatrix)*2,ncol=ncol(Umatrix)*2)
  islandLeft = 1
  islandRight = ncol(Umatrix)*2
  islandTop = 1
  islandBottom = nrow(Umatrix)*2

  Width = islandRight - islandLeft + 1
  Height = islandBottom - islandTop + 1

  # Vorschlag fuer eine Insel berechnen
  idealIsland = NULL
  if(!is.null(Bestmatches))
    idealIsland = bestUmatrixTranslation(Umatrix, Bestmatches)
  ##########
  # Shiny Fenster ----
  ##########
  UmatrixUi = fluidPage(
    useShinyjs(),
    sidebarLayout(position="right",
      mainPanel(
        #div(plotOutput("UmatrixPlot", width=as.character(Width), height=as.character(Height), click = "clickOnUmatrix"), style="margin-left:-70px")
        div(plotOutput("UmatrixPlot", click = "clickOnUmatrix"))
      ),
      div(style="max-width:1150px", # die sidebar ist dann auf 33% davon beschränkt
      sidebarPanel(
        # busy balken
        tags$head(tags$style(type="text/css", "
                             #loadmessage {
                             position: fixed;
                             top: 0px;
                             left: 0px;
                             width: 100%;
                             padding: 5px 0px 5px 0px;
                             text-align: center;
                             font-weight: bold;
                             font-size: 100%;
                             color: #000000;
                             background-color: #CCFF66;
                             z-index: 105;}")),
        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                         tags$div("Loading...",id="loadmessage")),

        fluidRow(
        actionButton("createIsland", "Create Island"),
        actionButton("resetIsland", "Reset Island")),

        fluidRow(
        ## DBT-ONLY
        #actionButton("loadBM", "loadBM"),
        #actionButton("loadCLS", "loadCLS")
        ),

        fluidRow(
          column(12, selectInput("BmSize", "Bestmatchsize:",
                                 c(
                                   "1x" = 1,
                                   "2x" = 2,
                                   "3x" = 3,
                                   "4x" = 4,
                                   "6x" = 6,
                                   "8x" = 8),
                                 selected = 2))
        ),
        fluidRow(
          column(12, selectInput("UmxSize", "Umatrixsize:",
                                 c(
                                   "0.5x" = 0.5,
                                   "1x" = 1,
                                   "1.5x" = 1.5,
                                   "2x" = 2,
                                   "3x" = 3,
                                   "4x" = 4),
                                 selected = UmatrixSize))
        ),

        checkboxInput("ShowWarnings", "Show Warnings", value = TRUE),

        htmlOutput("warnings"),

        br(),

        fluidRow(
          ## DBT-ONLY
          #actionButton("Save", "Save"),
          actionButton("quit", QuitButtonText)
          )
      )
    )))


  #Bestmatches = CheckBestmatches(Bestmatches, Cls, shiny=T)
  #CheckUmatrix(Umatrix, shiny=T)

  outputApp=runApp(list(ui = UmatrixUi, server = function(input, output, session){

    updateSelectInput(session, "UmxSize", selected = UmatrixSize)
    observe({UmatrixSize <<- as.numeric(input$UmxSize)})

    val = reactiveValues()
    val$UmatrixImg = NULL
    val$currentLines = c()
    val$Bestmatches = Bestmatches
    val$Cls = Cls
    val$Imx = matrix(0,nrow=nrow(Umatrix)*2,ncol=ncol(Umatrix)*2)

    val$Width <- Width
    val$Height <- Height

    val$Stretchfactor = 800/Width

    observe({
      if(input$ShowWarnings)
        output$warnings <- renderText(val$warnings)
      else
        output$warnings <- renderText(" ")
    })


    val$warnings = "<font color=\"#00FF00\"><b>No Warnings</b></font>"


    ##########
    # rerender Umatrix if island has changed ----
    ##########
    observe({
      CurrentDirectory <- getwd()
      setwd(tempdir())
      png("tmpUmatrix.png", width= val$Stretchfactor*val$Width,
          height = val$Stretchfactor*val$Height  )


      islandExists = !(sum(val$Imx)==0)


      umx <- showUmatrix_hlp(Umatrix, BmSize = as.numeric(input$BmSize)*2, Bestmatches = val$Bestmatches, Cls = val$Cls, Clean = T,
                        Imx=val$Imx, RemoveOcean = T, TransparentContours = F,
                        Toroid = T, MarkDuplicatedBestmatches=(islandExists&input$ShowWarnings))

      if(!is.null(val$Bestmatches))
        idealIsland = bestUmatrixTranslation(Umatrix, val$Bestmatches)

      if((!is.null(idealIsland))&&(!islandExists)){
        test <- data.frame(x=c(idealIsland$cols,
                               idealIsland$cols,
                               idealIsland$cols+ncol(Umatrix),
                               idealIsland$cols+ncol(Umatrix),
                               idealIsland$cols),
                           y=c(idealIsland$lines,
                               idealIsland$lines+nrow(Umatrix),
                               idealIsland$lines+nrow(Umatrix),
                               idealIsland$lines,
                               idealIsland$lines))

        # zeichne standard insel drueber
        umx = umx + geom_polygon(data=test,
                              aes(x=x, y=y),
                              size=0.5, alpha=0.01,
                              color="yellow", fill=NA)
      }

      print(umx)

      dev.off()
      val$UmatrixImg <- png::readPNG("tmpUmatrix.png")
      setwd(CurrentDirectory)
    })

    ##########
    # plot Umatrix ----
    ##########
    observe({
      output$UmatrixPlot <- renderPlot({
        if(!is.null(val$UmatrixImg)){
          par(mar=c(0,0,0,0))
          par(mai=c(0,0,0,0))
          par(xaxs = 'i',yaxs='i')
          rasterImage(val$UmatrixImg,0,0,1,1)

          lines(val$currentLines[,1], val$currentLines[,2],col="red")
        }
      },
      width = function(){val$Width * as.numeric(input$UmxSize) * val$Stretchfactor},
      height = function(){val$Height * as.numeric(input$UmxSize) * val$Stretchfactor})
    })

    ##########
    # react on click ----
    ##########
    observeEvent(input$clickOnUmatrix,{
      x <- isolate(input$clickOnUmatrix$x)
      y <- isolate(input$clickOnUmatrix$y)

      val$currentLines = rbind(val$currentLines,c(x,y))
    })

    #########
    # react on button: create island ----
    #########
    observeEvent(input$createIsland,{
      # check if there are already enough points to draw a complete polygon
      if(is.null(val$currentLines)) return()
      if(nrow(val$currentLines) < 3) return()

      Polygon = cbind( val$currentLines[,1]*ncol(Umatrix)*2+0.5, (1-val$currentLines[,2])*nrow(Umatrix)*2+0.5 )
      val$Imx <- createMaskFromPolygon(nrow(Umatrix)*2, ncol(Umatrix)*2, CutoutPol = Polygon)
      Imx <<- val$Imx
      val$currentLines = c()
      islandBM = UmxBestmatchesFromIsland(Umatrix, val$Bestmatches, Imx = val$Imx, Toroid = T)$Bestmatches

      if(!is.null(val$Bestmatches)){
        if(length(unique(islandBM[,1])) < nrow(val$Bestmatches)){
          val$warnings = paste0("<font color=\"#FF0000\"><b>",nrow(val$Bestmatches)-length(unique(islandBM[,1]))," Bestmatches are missing!</b></font>")
        }
        else if(length(islandBM[,1]) > nrow(val$Bestmatches)){
          val$warnings = "<font color=\"#FF0000\"><b>Some Bestmatches are duplicated! (marked with red triangles)</b></font>"
        }
        else{
          val$warnings = "<font color=\"#00FF00\"><b>No Warnings</b></font>"
        }
      }
    })

    #######
    # react on button: reset island ----
    #######
    observeEvent(input$resetIsland,{
      val$Imx = matrix(0,nrow=nrow(Umatrix)*2, ncol=ncol(Umatrix)*2)
      val$currentLines = c()
      val$warnings = "<font color=\"#00FF00\"><b>No Warnings</b></font>"
    })

    session$onSessionEnded(function() {
      print("program closed")
      stopApp(Imx)
    })

    observeEvent(input$quit, {
      stopApp(Imx)
    })

  }))

  return(outputApp)
}
