#'@method print mave.dim
#'@export

print.mave.dim<-function(x,...){
  cat('\nCall:\n',paste(deparse(x$call),sep="\n",collapse="\n"),"\n\n",sep="")
  arg='central space'
  if(substr(x$method,1,2)!='CS') arg='central mean space'
  p=length(x$dir)
  pp=1
  cv = c(x$cv0,x$cv)
  idx = which(!is.infinite(cv)) - 1
  cat('The cross-validation is run on dimensions of',idx,'\n')
  while(pp<=length(idx)){
    np=min(pp+10,length(idx))
    cat('Dimension\t')
    for(i in pp:np) cat(idx[i],'\t')
    cat('\n')
    cat('CV-value\t')
    for(i in pp:np) cat(round(cv[idx[i]+1],2),'\t')
    cat('\n')
    pp=np+1
  }
  cat('\n')
  d=which.min(cv)-1
  cat('The selected dimension of ',arg,' is ',d)
  cat('\n\n')
  #cat(paste('The matrix of the best',arg,'selected by cross-validation is of',d,'dimensions,','\nwhich is given below\n',sep=' '))
  #print(x$dir[[d]])
}
