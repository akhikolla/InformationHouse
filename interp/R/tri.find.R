tri.find<-function(tri.obj,x,y)
{
  if(!inherits(tri.obj,"triSht"))
      stop("tri.obj must be of class \"triSht\"")

  ans <- triFind(tri.obj$nt, tri.obj$x, tri.obj$y,
                 tri.obj$trlist[,"i1"], tri.obj$trlist[,"i2"], tri.obj$trlist[,"i3"],
                 x,y)

  ans
}
