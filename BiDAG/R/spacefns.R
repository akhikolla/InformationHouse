newspaceskel<-function(n,startspace,currspace,softlimit,hardlimit,posterior,blacklist,
                       MCMCtrace=NULL,mergetype="skeleton") {
  switch(mergetype,
         "dag" = { 
           mdag<-dag.threshold(MCMCtrace,pbarrier=posterior,pdag=FALSE)
           newadj<-1*(!blacklist&(startspace|mdag))
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-(1*(!blacklist&mdag))[,toomanyneib]}
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-currspace[,toomanyneib]}
         },
         "cpdag" = { 
           mcp<-dag.threshold(n,MCMCtrace,pbarrier=posterior,pdag=TRUE)
           newadj<-1*(!blacklist&(startspace|mcp))
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             mdag<-dag.threshold(MCMCtrace,pbarrier=posterior,pdag=FALSE)
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|mdag)))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-mdag[,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
         },
         "skeleton" = { 
           mskel<-1*(dag.threshold(MCMCtrace,pbarrier=posterior,pdag=FALSE)|t(dag.threshold(MCMCtrace,pbarrier=posterior,pdag=FALSE)))
           newadj<-1*(!blacklist&(startspace|mskel))
           toomanyneib<-which(apply(newadj,2,sum)>4)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|dag.threshold(MCMCtrace,pbarrier=posterior,pdag=TRUE))))[,toomanyneib]
           }
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             mdag<-dag.threshold(MCMCtrace,pbarrier=posterior,pdag=FALSE)
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|mdag)))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-mdag[,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
         }
  )
  return(newadj)
}

newspacemap<-function(n,startspace,currspace,softlimit,hardlimit,blacklist, 
                      maxN,MCMCtrace=NULL,mergetype="skeleton",accum) {
  switch(mergetype,
         "dag" = { 
           maxdag<-MCMCtrace[[maxN]]
           newadj<-1*(!blacklist&(startspace|maxdag))
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-(1*(!blacklist&maxdag))[,toomanyneib]}
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-currspace[,toomanyneib]}
           if(accum) {
             newnewadj<-1*(newadj|currspace)
             toomanyneib<-which(apply(newnewadj,2,sum)>hardlimit)
             if(length(toomanyneib)>0){newnewadj[,toomanyneib]<-newadj[,toomanyneib]}
             newadj<-newnewadj
           }
         },
         "cpdag" = { 
           maxcp<-dagadj2cpadj(MCMCtrace[[maxN]])
           newadj<-1*(!blacklist&(startspace|maxcp))
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|MCMCtrace[[maxN]])))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-MCMCtrace[[maxN]][,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
           if(accum) {
             newnewadj<-1*(newadj|currspace)
             toomanyneib<-which(apply(newnewadj,2,sum)>hardlimit)
             if(length(toomanyneib)>0){newnewadj[,toomanyneib]<-newadj[,toomanyneib]}
             newadj<-newnewadj
           }
         },
         "skeleton" = { 
           maxskel<-1*(MCMCtrace[[maxN]]|t(MCMCtrace[[maxN]]))
           newadj<-1*(!blacklist&(startspace|maxskel))
           toomanyneib<-which(apply(newadj,2,sum)>7)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|dagadj2cpadj(MCMCtrace[[maxN]]))))[,toomanyneib]
           }
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|MCMCtrace[[maxN]])))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-MCMCtrace[[maxN]][,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
           if(accum) {
             newnewadj<-1*(newadj|currspace)
             toomanyneib<-which(apply(newnewadj,2,sum)>hardlimit)
             if(length(toomanyneib)>0){newnewadj[,toomanyneib]<-newadj[,toomanyneib]}
             newadj<-newnewadj
           }
         }
  )
  return(newadj)
}

definestartspace<-function(alpha,param,cpdag=FALSE,algo="pc") {
  if(is.null(alpha)) {alpha<-0.1}
  
  local_type <- param$type
  if(local_type=="usr") {
    if(param$pctesttype%in%c("bde","bge","bdecat")) {
      local_type<-param$pctesttype
    } 
  }
  
  if(param$DBN){
    if(param$stationary) {
    othersliceskel <- definestartspace(alpha,param$otherslices,cpdag=FALSE,algo="pc")
    firstsliceskel <- definestartspace(alpha,param$firstslice,cpdag=FALSE,algo="pc")
    startspace <- othersliceskel
    startspace[param$intstr$rows,param$intstr$cols] <- 1*(startspace[param$intstr$rows,param$intstr$cols] | firstsliceskel[param$intstr$rows,param$intstr$cols])
    } else {
      skels<-list()
      skels[[1]]<-definestartspace(alpha,param$paramsets[[1]],cpdag=FALSE,algo="pc")
      startspace<-skels[[1]]
      for(i in 2:length(param$paramsets)) {
        skels[[i]]<-definestartspace(alpha,param$paramsets[[i]],cpdag=FALSE,algo="pc")
        startspace<-1*(skels[[i]]|startspace)
      }
}
  } else { # otherwise use old versions
    
    if(local_type=="bde") {
      if(cpdag){
        pc.skel<-pc(suffStat = list(d1=param$d1,d0=param$d0,data=param$data),
                    indepTest = weightedbinCItest, alpha = alpha, labels = colnames(param$data),
                    verbose = FALSE)
        
      } else {
        pc.skel<-pcalg::skeleton(suffStat = list(d1=param$d1,d0=param$d0,data=param$data),
                                 indepTest = weightedbinCItest, alpha = alpha, labels = colnames(param$data),
                                 verbose = FALSE)
      }
    } else if(local_type=="bdecat") {
      if(cpdag){
        pc.skel<-pc(suffStat = param,
                    indepTest = weightedcatCItest, alpha = alpha, labels = colnames(param$data),
                    verbose = FALSE)
        
      } else {
        pc.skel<-pcalg::skeleton(suffStat = param,
                                 indepTest = weightedcatCItest, alpha = alpha, labels = colnames(param$data),
                                 verbose = FALSE)
      }
    } else if(local_type=="bge") {
      if(is.null(param$weightvector)) {
        cormat<-cor(param$data)
        N<-nrow(param$data)
      } else { N<-sum(param$weightvector)
      cormat<-cov.wt(param$data,wt=param$weightvector,cor=TRUE)$cor}
      if(cpdag){
        pc.skel<-pcalg::pc(suffStat = list(C = cormat, n = N),
                           indepTest = gaussCItest,
                           alpha=alpha,labels=colnames(param$data),skel.method="stable",verbose = FALSE)
      } else {
        pc.skel<-pcalg::skeleton(suffStat = list(C = cormat, n = N),
                                 indepTest = gaussCItest,
                                 alpha=alpha,labels=colnames(param$data),method="stable",verbose = FALSE)
      }
    } else if (local_type=="usr") {

        if(cpdag){
        pc.skel<-pc(suffStat = param,
                    indepTest = usrCItest, alpha = alpha, labels = colnames(param$data),
                    verbose = FALSE)
        
      } else {
        pc.skel<-pcalg::skeleton(suffStat = param,
                                 indepTest = usrCItest, alpha = alpha, labels = colnames(param$data),
                                 verbose = FALSE)
      }
    }
    
    g<-pc.skel@graph
    startspace<-1*(graph2m(g))
  }
  return(startspace)
}


