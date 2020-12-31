M3d2list<-function(BB,xnames){
	p=dim(BB)[1]
	B=list()
	if(is.null(xnames)){
	  xnames = 1:p
	}
	for(i in 1:dim(BB)[2]){
		B[[i]] = as.matrix(BB[,1:i,i])
		colnames(B[[i]])<-paste('dir',1:i,sep='')
    rownames(B[[i]])<-xnames
	}
	return(B)
}

