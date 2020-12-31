triangles<-function(tri.obj){
  if(!inherits(tri.obj,"triSht"))
    stop("tri.obj must be of class \"triSht\"")

  ret<-tri.obj$trlist
  colnames(ret)<-c("node1","node2","node3","tr1","tr2","tr3","arc1","arc2","arc3")
  ret
}

arcs<-function(tri.obj){
  if(!inherits(tri.obj,"triSht"))
    stop("tri.obj must be of class \"triSht\"")

  ret<-cbind(tri.obj$arcs[,"from"],tri.obj$arcs[,"to"])
  colnames(ret)<-c("from","to")
  ret
}

area<-function(tri.obj){
    if(!inherits(tri.obj,"triSht"))
        stop("tri.obj must be of class \"triSht\"")

    ret<-tri.obj$cclist[,"area"]

    ret
}
