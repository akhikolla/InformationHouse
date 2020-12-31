summary.triSht<-function(object, ...)
{
  if(!inherits(object,"triSht"))
    stop("object must be of class \"triSht\"")

  ans<-list(n=object$n,
            na=object$narcs,
            nb=object$nchull,
            nt=object$nt,
            call=object$call)
  class(ans)<-"summary.triSht"
  ans
}
