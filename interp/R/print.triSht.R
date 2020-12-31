print.triSht<-function(x,...)
{
  if(!inherits(x,"triSht"))
    stop("x must be of class \"triSht\"")
  cat("Delauney triangulation, node and triangle indices:\n")
  cat("triangle: nodes (a,b,c), neighbour triangles [i,j,k] \n")
  for (i in 1:x$nt)
    {
      cat(i,": (",x$trlist[i,"i1"],",",x$trlist[i,"i2"],",",x$trlist[i,"i3"],"), [",x$trlist[i,"j1"],",",x$trlist[i,"j2"],",",x$trlist[i,"j3"],"]\n",sep="")
  }
  cat("boundary nodes: ", x$chull, "\n", sep=" ")
}


