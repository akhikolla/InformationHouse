convex.hull<-function(tri.obj, plot.it=FALSE, add=FALSE,...)
{
  if(!inherits(tri.obj,"triSht"))
    stop("tri.obj must be of class \"triSht\"")

  ret<-list(x=tri.obj$x[tri.obj$chull],
            y=tri.obj$y[tri.obj$chull],
            i=tri.obj$chull)
  if(plot.it)
    {
      if (!add)
        {
          plot.new()
          plot.window(range(ret$x), range(ret$y), "")
        }
      lines(cbind(ret$x,ret$x[1]),cbind(ret$y,ret$y[1]), ...)
      invisible(ret)
    }
  else
    ret
}
