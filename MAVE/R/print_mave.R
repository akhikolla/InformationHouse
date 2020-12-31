#'@method print mave
#'@export
print.mave<-function(x,...){
  cat('\nCall:\n',paste(deparse(x$call),sep="\n",collapse="\n"),"\n\n",sep="")
  arg='central space'
  if(substr(x$method,1,2)!='CS') arg='central mean space'
  cat(arg, ' of dimensions ', 1:x$max.dim, ' are computed\n\n')

  #len=length(x$which.dim)
  #if(sum('cv'==names(x))==0){
  #for(i in 1:min(len,3)){
  #  cat(mave.dir(x,dim=x$which.dim[i]))
  #}
  #cat('(only the first 3 ',arg,' are displayed,
  #      to display the space of dimension k, call mave.dir(dr,dim) or object$dir[[k]])\n')
  #}
}


