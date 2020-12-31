networkflow_odd<-function(a,b,C,threads){
  a<-rbind(a,1)
  b<-rbind(b,1)
  #This is an atrocity. However, the fastest I could come up with, since rbind is horribly slow. 
  C<-unname(as.matrix(rbindlist(list(data.frame(C),data.frame(matrix(max(C)+1,1,dim(C)[2]))))))
  C<-cbind(C,matrix(max(C)*10^16,dim(C)[1],1))
  C[dim(C)[1],dim(C)[2]]<-0
  res<-networkflow(a,b,C,threads)
  res$frame[,3][res$frame[,1]==length(a)]<-0
  res$plan<-res$plan[1:(length(a)-1),1:length(b)-1]
  res$potential<-res$potential[c(-length(a),-(length(a)+length(b)))]
  return(res)
}