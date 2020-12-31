orderMCMCmain<-function(param,iterations,stepsave,MAP=TRUE, posterior=0.5,
                        startorder=c(1:n),moveprobs,plus1=FALSE,chainout=TRUE,
                        scoreout=FALSE,startspace=NULL,blacklist=NULL,gamma=1,verbose=FALSE,alpha=NULL,
                        hardlimit=ifelse(plus1,15,22),
                        cpdag=FALSE,addspace=NULL,scoretable=NULL){
  result<-list()
  maxobj<-list()
  MCMCtraces<-list()
  
  n<-param$n
  nsmall<-param$nsmall
  matsize<-ifelse(param$DBN,n+nsmall,n)
  
  #defining startorder and updatenodes 
  if(!param$DBN) { 
    
    if(param$bgn!=0) {
      updatenodes<-c(1:n)[-param$bgnodes]
    } else { 
      updatenodes<-c(1:n)
    }
    
  } else { #for DBNs startorder is defined in main.R
    updatenodes<-c(1:nsmall)
  }
  maxorder<-startorder
  
  #creating blacklist objects
  if (is.null(blacklist)) {
    blacklist<-matrix(0,nrow=matsize,ncol=matsize)
  }
  diag(blacklist)<-1
  if(!is.null(param$bgnodes)) {
    for(i in param$bgnodes) {
      blacklist[,i]<-1
    }
  }
  blacklistparents<-list()
  for  (i in 1:matsize) {
    blacklistparents[[i]]<-which(blacklist[,i]==1)
  }
  
  #defining startspace and startskeleton
  if (is.null(startspace)) {
    startspace<-definestartspace(alpha,param,cpdag=cpdag,algo="pc")
  } 
  startskeleton<-1*(startspace&!blacklist)
  
  if(!is.null(addspace)) { startskel<-1*((addspace|startskeleton)&!blacklist)
  } else {startskel<-startskeleton }
  
  print(paste("maximum parent set size is", max(apply(startskel,2,sum))))
  if(max(apply(startskel,2,sum))>hardlimit) {
    stop("the size of maximal parent set is higher that the hardlimit; redifine the search space or increase the hardlimit!")
  }
  
  tablestart<-Sys.time()
  #computing score tables
  ptab<-listpossibleparents.PC.aliases(startskel,isgraphNEL=FALSE,n,updatenodes)
  
  if (verbose) {
    print("skeleton ready")
   flush.console()
  }
  
  if (plus1==FALSE) {
  
    parenttable<-ptab$parenttable # basic parenttable without plus1 lists
    aliases<-ptab$aliases #aliases for each node since all nodes in parent tables are named as 1,2,3,4. not real parent names
    numberofparentsvec<-ptab$numberofparentsvec
    numparents<-ptab$numparents
    rowmaps<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
    scoretable<-scorepossibleparents.alias(parenttable,aliases,n,param,updatenodes,rowmaps,
                                         numparents,numberofparentsvec)
    posetparenttable<-poset(parenttable,numberofparentsvec,rowmaps,n,updatenodes)
  
  if(MAP==TRUE){
    maxmatrices<-posetscoremax(posetparenttable,scoretable,numberofparentsvec,rowmaps,
                               n,plus1lists=NULL,updatenodes)
    tableend<-Sys.time()
    MCMCresult<-orderMCMCbasemax(n,nsmall,startorder,iterations,stepsave,moveprobs,parenttable,
                                 scoretable,aliases,numparents,rowmaps,maxmatrices,
                                 numberofparentsvec,gamma=gamma,bgnodes=param$bgnodes,matsize=matsize)

    mcmcend<-Sys.time()
  } else {
    bannedscore<-poset.scores(posetparenttable,scoretable,numberofparentsvec,rowmaps,
                              n,plus1lists=NULL,ptab$numparents,updatenodes=updatenodes)
    tableend<-Sys.time()
    MCMCresult<-orderMCMCbase(n,nsmall,startorder,iterations,stepsave,moveprobs,parenttable,
                              scoretable,aliases,numparents,rowmaps,
                              bannedscore,numberofparentsvec,gamma=gamma,
                              bgnodes=param$bgnodes,matsize=matsize)
    mcmcend<-Sys.time()
    
  }
  
} else {
  parenttable<-ptab$parenttable # basic parenttable without plus1 lists
  aliases<-ptab$aliases #aliases for each node since all nodes in parent tables are done as 1,2,3,4... not real parent names
  numberofparentsvec<-ptab$numberofparentsvec
  numparents<-ptab$numparents
  plus1lists<-PLUS1(matsize,aliases,updatenodes,blacklistparents)
  rowmaps<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)
  if(is.null(scoretable)) {
    scoretable<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,param,updatenodes,
                                           rowmaps,numparents,numberofparentsvec) 
    }
    posetparenttable<-poset(parenttable,numberofparentsvec,rowmaps,n,updatenodes)
  
  if(MAP==TRUE){
    maxmatrices<-posetscoremax(posetparenttable,scoretable,numberofparentsvec,
                               rowmaps,n,plus1lists=plus1lists,updatenodes)
  } else {
    bannedscore<-poset.scores(posetparenttable,scoretable,ptab$numberofparentsvec,rowmaps,
                              n,plus1lists=plus1lists,ptab$numparents,updatenodes)
  }

  if(verbose) {
    print(paste("score tables calculated, MCMC plus1 starts"))
    flush.console()
  }

  tableend<-Sys.time()
    
  #running MCMC scheme   
  if(MAP) {
    MCMCresult<-orderMCMCplus1max(n,nsmall,startorder,iterations,stepsave,moveprobs,parenttable,
                                      scoretable,aliases,numparents,rowmaps,plus1lists,maxmatrices,numberofparentsvec,
                                      gamma=gamma,bgnodes=param$bgnodes,matsize=matsize)
    
    } else {
    MCMCresult<-orderMCMCplus1(n,nsmall,startorder,iterations,stepsave,moveprobs,parenttable,
                                   scoretable,aliases,numparents,rowmaps,plus1lists,
                                   bannedscore,numberofparentsvec,gamma=gamma,
                                   bgnodes=param$bgnodes,matsize=matsize)
    }
    
   mcmcend<-Sys.time()
  
}

  #defining result object
  if(chainout) {
    if(param$DBN) {
      MCMCchain<-lapply(MCMCresult$incidence,function(x)DBNtransform(x,param=param))
      MCMCtraces$incidence<-MCMCchain
      MCMCtraces$orders<-lapply(MCMCresult$orders,order2var,varnames=param$firstslice$labels)
    } else {
      MCMCtraces$incidence<-lapply(MCMCresult$incidence,function(x)assignLabels(x,param$labels))
      MCMCtraces$orders<-lapply(MCMCresult$orders,order2var,varnames=param$labels)
    }
    MCMCtraces$DAGscores<-MCMCresult$DAGscores
    MCMCtraces$orderscores<-MCMCresult$orderscores
  } 
  
  
  maxobj<-storemaxMCMC(MCMCresult,param)
  maxN<-which.max(unlist(MCMCresult$DAGscores))
  maxobj$reach<-maxN
  
  if (scoreout){
    if(chainout){output<-4}
    else{output<-3}
  } else {
    if(chainout) {output<-2}
    else {output<-1}
  }
  
  result$max<-maxobj
  attr(result$max,"class")<-"MCMCmax"
  result$info<-list()
  tabletime<-tableend-tablestart
  if(units(tabletime)=="mins") {
    tabletime<-as.numeric(tabletime*60)
  }
  mcmctime<-mcmcend-tableend
  if(units(mcmctime)=="mins") {
    mcmctime<-as.numeric(mcmctime*60)
  }
  result$info$times<-c(tabletime,mcmctime)
  
  switch(as.character(output),
         "1"={ # return only maximum DAG and order
           return(result)
         },
         "2"={ # return all MCMC all saved MCMC steps: incidence, DAGscore, orderscore and order and max result
           result$chain<-MCMCtraces
           attr(result$chain,"class")<-"MCMCtrace"
           attr(result,"class")<-"MCMCres"
         },
         "3"={ # return max DAG, order, last search space incidence and all scoretables
           result$space<-list()
           result$space$adjacency<-startskel
           result$space$scoretable<-scoretable
           attr(result$space,"class")<-"MCMCspace"
         },
         "4"={ # return all MCMC all saved MCMC steps,max result,last search space and scoretables
           result$space<-list()
           result$space$adjacency<-startskel
           result$space$scoretable<-scoretable
           attr(result$space,"class")<-"MCMCspace"
           result$chain<-MCMCtraces
           attr(result$chain,"class")<-"MCMCtrace"
         }
  )
  return(result)
  
}