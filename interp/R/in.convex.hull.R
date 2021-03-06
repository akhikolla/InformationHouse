
in.convex.hull<-function(tri.obj,x,y,eps=1E-16,strict=TRUE)
{
  if(!inherits(tri.obj,"triSht"))
    stop("tri.obj must be of class \"triSht\"")
  if(length(x)!=length(y))
    stop("x and y must be of same length")
  n<-length(x)
  if(n==0)
    stop("length of x (resp. y) is 0")
  inhull <- inHull(tri.obj,x,y,eps)
  if(!strict)
      inhull <- inhull | onHull(tri.obj,x,y,eps)
  inhull
}


