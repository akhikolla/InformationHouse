PlotProjectedPoints <-  function(Points,Cls,BMUorProjected=F,PlotLegend=FALSE,xlab='X',ylab='Y',main="Projected Points",PointSize=2.5){
  # PlotProjectedPoints(ProjectedPoints)
  # PlotProjectedPoints(ProjectedPoints,Cls)
  # plot ProjectedPoints colored by Cls
  #
  # INPUT
  # Points[1:n,1:3] or [1:n,1:2]       the coordinates of the projected Points																				data, BestMatches are also allowed
  #																		 if [1:n,1:3] assuming first column=key
  # Optional
  # Cls[1:n]              integer vector of class identifiers
  # BMUorProjected        Default(FALSE): ProjectedPoints[1:n,Key,x,y]
  #                       TRUE: BestMatches[1:n,Key,Lines,Columns]=BestMatches[1:n,Key,y,x]
  #
  # PlotLegend            lengend of Colors and Classens
  # xlab                  Label of x-Axis
  # ylab                  Label of y-Axis
  # ...                   All arguments from plot außer xlim und ylim
  #
  # author: MT 07/2015
  # uses DefaultColorSymbSequence and DefaultColorSequence
  
  if (!is.matrix(Points))
    stop('ProjectedPoints has to be a matrix')
  
  AnzData = nrow(Points)
  if (missing(Cls)) {
    colnames(Points) = c('X', 'Y')
    requireNamespace('ggplot2')
    ggobject = ggplot(as.data.frame(Points), aes(X, Y)) + geom_point() +
      ggtitle(main) + coord_fixed(ratio = 1)
    
    print(ggobject)
    return(ggobject)
    Cls = rep(1, AnzData)
    DefaultColorSeq = 'black'
  } else{
    if (is.null(Cls)) {
      colnames(Points) = c('X', 'Y')
      ggobject = ggplot(as.data.frame(Points), aes(X, Y)) + geom_point() +
        ggtitle(main) + coord_fixed(ratio = 1)
      
      print(ggobject)
      return(ggobject)
    } else{
      DefaultColorSeq <- ProjectionBasedClustering::DefaultColorSequence
    }
    
  }
  if (sum(Cls) == length(Cls)) {
    colnames(Points) = c('X', 'Y')
    ggobject = ggplot(as.data.frame(Points), aes(X, Y)) + geom_point() +
      ggtitle(main) + coord_fixed(ratio = 1)
    
    print(ggobject)
    return(ggobject)
  }
  c = ncol(Points)
  if (c > 3 | c < 2)
    stop(paste0('Wrong number of Columns of ProjectedPoints: ', c))
  if (c == 3) {
    #With Key
    if (BMUorProjected) {
      Y = Points[, 2]
      X = Points[, 3]
    } else{
      X = Points[, 2]
      Y = Points[, 3]
    }
  } else{
    #Without Key
    if (BMUorProjected) {
      Y = Points[, 1]
      X = Points[, 2]
    } else{
      X = Points[, 1]
      Y = Points[, 2]
    }
  }
  if (length(Cls) < AnzData) {
    warning(paste0(
      "ClassPlot(): too few Classes",
      length(Cls),
      " dummy added",
      AnzData
    ))
    NewCls = rep(1, AnzData) * max(Cls + 1)
    NewCls[1:AnzData] = Cls
    Cls = NewCls
  }
  if (length(Cls) > AnzData) {
    warning("ClassPlot(): Cls too long, shortened")
    Cls = Cls[1:AnzData]
  }
  
  Points = cbind(X, Y)
  colnames(Points) = c('X', 'Y')
  df = as.data.frame(Points)
  #names=DefaultColorSequence()
  
  # noClsL=ClassCount(Cls)
  # noCls=noClsL$numberOfClasses
  #colornames=DefaultColors4Cls(Cls)
  hex_hlp = function(cname) {
    colMat <- grDevices::col2rgb(cname)
    grDevices::rgb(red = colMat[1,] / 255,
                   green = colMat[2,] / 255,
                   blue = colMat[3,] / 255)
  }
  vec = hex_hlp(ProjectionBasedClustering::DefaultColorSequence)
  ind = unique(Cls)
  # print('test')
  # if(length(ind)==1){
  #   df=cbind(df,colornames)
  #   ggobject=ggplot(df, aes(X, Y)) + geom_point(size=PointSize,color='green')
  #     ggtitle(main)
  #     print('test')
  # }else{
  colornames = Cls
  for (i in ind) {
    ind2 = which(Cls == i)
    colornames[ind2] = vec[i]
  }
  # Green    = "#00FF00FF"
  # Red      = "#FF0000FF"
  # Cyan     = "#00FFFFFF"
  # Magenta  = "#FF00FFFF"
  # Yellow   = "#FFFF00FF"
  # Blue     = "#0000FFFF"
  # Black    = "#000000"
  # names=c(Green,Red,Cyan,Magenta,Yellow,Blue,Black,names)
  # colornamesTmp=names[1:noCls]
  # colornames=as.character(Cls)
  # for(i in 1:noCls){
  #   ind=which(Cls==noClsL$uniqueClasses[i])
  #   colornames[ind]=names[i]
  # }
  #print(colornames)
  df = cbind(df, colornames)

  ggobject = ggplot(df, aes(X, Y, colour = colornames)) + geom_point(size =
                                                                       PointSize, colour = colornames) +
    ggtitle(main) + #coord_fixed(ratio = 1)+
    theme(
      panel.background = element_blank(),
      legend.key = element_blank(),
      axis.line = element_line(colour = 'black'),
      axis.title.y = element_text(size = rel(2), angle = 00),
      axis.title.x = element_text(size = rel(2), angle = 00),
      axis.text.x = element_text(size = rel(2)),
      axis.text.y = element_text(size = rel(2)),
      plot.title =  element_text(size = rel(2)),
      legend.text =  element_text(size = rel(2)),
      legend.title =  element_text(size = rel(2)),
      legend.position 	= 'none'
    )+theme(plot.title = element_text(hjust = 0.5),axis.title.x  = element_text(hjust = 0.5),axis.title.y  = element_text(vjust = 0.5))
  
  # }
  #scatter plot to be square shaped
  ggobject=ggobject+coord_equal()+expand_limits(x=ggplot_build(ggobject)$panel$ranges[[1]]$y.range,
                  y=ggplot_build(ggobject)$panel$ranges[[1]]$x.range)
  print(ggobject)
  return(ggobject)
}