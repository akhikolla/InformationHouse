drawPoliticalMap <- function(cls, dirichletTiles, Lines, Columns, Toroid, add=F, alphaV=0.2){ #, lines = NULL, columns = NULL, projectedData = NULL, add=F, alphaV=0.2){
  # draws the generated precalculated dirichletTiles.
  # INPUT
  # cls             classification. Will be used to choose colors for the voronoi cells.
  # dirichletTiles  calculated Voronoi cells
  # add             if true, political map will be drawn onto the already existing plot
  # alphaV          transparancy value.
  # OUTPUT
  # -
  
  #    requireRpackage("deldir")
  
  #if(is.null(dirichletTiles)){
  #  if(is.null(lines) | is.null(columns) | is.null(projectedData))
  #    stop("if dirichletTiles is not given, lines, columns and projectedData are necessary for drawing the political map")
  #}
  
  ccc <- DefaultColorSequence()[cls]
  if(add==T) {
    ccc <- alpha(ccc, alphaV)
    plot(dirichletTiles,fillcol=ccc,close=T, showpoints=T, border=F, add=T)
  }
  else{
    if(Toroid)
      plot(1, type="n", ylim=c(Lines*2,1), xlim=c(1,Columns*2))
    else
      plot(1, type="n", ylim=c(Lines,1), xlim=c(1,Columns))
    
    plot(dirichletTiles,fillcol=ccc,close=T, showpoints=T, border=F, add=T)
  }
}
