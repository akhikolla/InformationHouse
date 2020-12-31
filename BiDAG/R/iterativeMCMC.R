iterativeMCMCplus1<-function(param,iterations,stepsave,plus1it=NULL,MAP=TRUE, posterior=0.5,
                             startorder=NULL,moveprobs,softlimit=9,hardlimit=14,chainout=FALSE,
                             scoreout=FALSE,startspace=NULL,blacklist=NULL,gamma=1,verbose=FALSE,alpha=NULL,
                             cpdag=FALSE,mergecp="skeleton",addspace=NULL,scoretable=NULL,
                             accum) {
  n<-param$n
  nsmall<-param$nsmall
  matsize<-ifelse(param$DBN,n+nsmall,n)
  
  
  objsizes<-list()
  maxlist<-list()
  maxobj<-list()
  updatenodeslist<-list()
  MCMCtraces<-list()
  
  if(!param$DBN) {
    if(param$bgn!=0) {
      updatenodes<-c(1:n)[-param$bgnodes]
    } else {
      updatenodes<-c(1:n)
    }
  } else {
    updatenodes<-c(1:nsmall)
  }
  
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
  
  
  starttimeall<-Sys.time()
  maxorder<-startorder
  

  ptab<-listpossibleparents.PC.aliases(startskel,isgraphNEL=FALSE,n,updatenodes)
  
  
  if (verbose) {
    print("skeleton ready")
    flush.console()
  }
    ##################################
    #calculating initial score tables#
    ##################################
  
    parenttable<-ptab$parenttable # basic parenttable without plus1 lists
    aliases<-ptab$aliases #aliases for each node since all nodes in parent tables are done as 1,2,3,4... not real parent names
    numberofparentsvec<-ptab$numberofparentsvec
    numparents<-ptab$numparents
    plus1lists<-PLUS1(matsize,aliases,updatenodes,blacklistparents)
    rowmaps<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)

    if(is.null(scoretable)) {
      scoretable<-scorepossibleparents.PLUS1(parenttable=parenttable,plus1lists=plus1lists,
                                             n=n,param=param,updatenodes=updatenodes,
                                             rowmaps,numparents,numberofparentsvec) }
      posetparenttable<-poset(parenttable,numberofparentsvec,rowmaps,n,updatenodes)
    
    if(MAP==TRUE){
      maxmatrices<-posetscoremax(posetparenttable,scoretable,numberofparentsvec,
                                 rowmaps,n,plus1lists=plus1lists,updatenodes)
      endtimetab<-Sys.time()
    } else {
      bannedscore<-poset.scores(posetparenttable,scoretable,ptab$numberofparentsvec,rowmaps,
                                n,plus1lists=plus1lists,ptab$numparents,updatenodes)
      endtimetab<-Sys.time()
    }
    oldadj<-startskeleton
    
  ###########################
  #MCMC chain is constructed#  
  ###########################
    
  i<-1  
  if(is.null(plus1it)) plus1it<-100
  while (length(updatenodes)>0 & i<=plus1it){
        starttimeit<-Sys.time()
        if(i>1){
          newptab<-listpossibleparents.PC.aliases(newadj,isgraphNEL=FALSE,n,updatenodes)
          parenttable[updatenodes]<-newptab$parenttable[updatenodes] # basic parenttable without plus1 lists
          aliases[updatenodes]<-newptab$aliases[updatenodes] #aliases for each node since all nodes in parent tables are done as 1,2,3,4... not real parent names
          
          numberofparentsvec[updatenodes]<-newptab$numberofparentsvec[updatenodes]
          numparents[updatenodes]<-newptab$numparents[updatenodes]
          newplus1lists<-PLUS1(matsize,aliases,updatenodes,blacklistparents)
          plus1lists$mask[updatenodes]<- newplus1lists$mask[updatenodes]
          plus1lists$parents[updatenodes]<- newplus1lists$parents[updatenodes]
          plus1lists$aliases[updatenodes]<- newplus1lists$aliases[updatenodes]
          rowmaps[updatenodes]<-parentsmapping(parenttable,numberofparentsvec,n,updatenodes)[updatenodes]
          scoretable[updatenodes]<-scorepossibleparents.PLUS1(parenttable,plus1lists,n,param,updatenodes,rowmaps,numparents,numberofparentsvec)[updatenodes]
          posetparenttable[updatenodes]<-poset(parenttable,numberofparentsvec,rowmaps,n,updatenodes)[updatenodes]
          
          if (MAP) {
            newmaxmatrices<-posetscoremax(posetparenttable,scoretable,numberofparentsvec,
                                          rowmaps,n,plus1lists=plus1lists,updatenodes)
            maxmatrices$maxmatrix[updatenodes]<- newmaxmatrices$maxmatrix[updatenodes]
            maxmatrices$maxrow[updatenodes]<- newmaxmatrices$maxrow[updatenodes]
            #objsizes$sumscore<-object.size(maxmatrices)
            
          } else {
            newbannedscore<-poset.scores(posetparenttable,scoretable,numberofparentsvec,rowmaps,
                                         n,plus1lists=plus1lists,numparents,updatenodes)
            bannedscore[updatenodes]<-newbannedscore[updatenodes]
            #objsizes$sumscore<-object.size(bannedscore)
          }
          
          if(verbose) {
            print(paste("MCMC plus1 iteration",i))
            flush.console()
          }
        } else {
          if(verbose) {
            print(paste("score tables completed, MCMC plus1 starts"))
            flush.console()
          }
        }
        if(MAP) {
          starttimemcmc<-Sys.time()
          MCMCresult<-orderMCMCplus1max(n,nsmall,startorder,iterations,stepsave,moveprobs,parenttable,
                                        scoretable,aliases,numparents,rowmaps,plus1lists,maxmatrices,numberofparentsvec,
                                        gamma=gamma,bgnodes=param$bgnodes,matsize=matsize)
          endtimeit<-Sys.time()
          if(verbose) {
            print(endtimeit-starttimeit)
            flush.console()
          }
        } else {
          starttimemcmc<-Sys.time()
          MCMCresult<-orderMCMCplus1(n,nsmall,startorder,iterations,stepsave,moveprobs,parenttable,
                                     scoretable,aliases,numparents,rowmaps,plus1lists,
                                     bannedscore,numberofparentsvec,gamma=gamma,bgnodes=param$bgnodes,
                                     matsize=matsize)
          endtimeit<-Sys.time()
          if(verbose){
            print(endtimeit-starttimeit)
            flush.console()
          }
        }
        
        if(chainout) {
          if(param$DBN) {
            MCMCtraces$incidence[[i]]<-lapply(MCMCresult$incidence,function(x)DBNtransform(x,param=param))
            MCMCtraces$orders[[i]]<-lapply(MCMCresult$orders,order2var,varnames=param$firstslice$labels)
          } else {
            MCMCtraces$incidence[[i]]<-lapply(MCMCresult$incidence,function(x)assignLabels(x,param$labels))
            MCMCtraces$orders[[i]]<-lapply(MCMCresult$orders,order2var,varnames=param$labels)
          }
          MCMCtraces$DAGscores[[i]]<-MCMCresult$DAGscores
          MCMCtraces$orderscores[[i]]<-MCMCresult$orderscores
        } 
        
        maxobj<-storemaxMCMC(MCMCresult,param)
        maxlist[[i]]<-maxobj
        
        maxN<-which.max(unlist(MCMCresult$DAGscores))

        if(i>1){
          if (maxobj$score>maxscore){
            maxDAG<-maxobj$DAG
            maxorder<-maxobj$order
            maxscore<-maxobj$score
            maxit<-i
          } 
        } else {
          maxDAG<-maxobj$DAG
          maxscore<-maxobj$score
          maxorder<-maxobj$order
          maxit<-1
        }
        if (MAP) {
          newadj<-newspacemap(n,startskeleton,oldadj,softlimit,hardlimit,blacklist,
                               maxN=maxN,MCMCtrace=MCMCresult[[1]],mergetype=mergecp,
                               accum=accum)
        } else {
          newadj<-newspaceskel(n,startskeleton,oldadj,softlimit,hardlimit,posterior,
                                blacklist,MCMCtrace=MCMCresult[[1]],mergetype=mergecp)
        }
        updatenodes<-which(apply(newadj==oldadj,2,all)==FALSE)
        updatenodeslist[[i]]<-updatenodes
        oldadj<-newadj
        startorder<-c(MCMCresult$orders[[maxN]],param$bgnodes)
        i<-i+1
      }
      
    
    addedge<-sum(newadj)-sum(startskeleton)
    endtimeall<-Sys.time()
    result<-list()
    if (scoreout){
      if(chainout){output<-4}
      else{output<-3}
    } else {
      if(chainout) {output<-2}
      else {output<-1}
    }
    result$maxtrace<-maxlist
    result$max<-maxlist[[length(maxlist)]]
    if(param$DBN) {
      result$startspace<-DBNtransform(startskeleton,param)
      result$endspace<-DBNtransform(oldadj,param)
      } else {
          result$startspace<-startskeleton
          result$endspace<-oldadj
      }
    switch(as.character(output),
           "1"={ # return only maximum DAG and order
                 # do not need to do anything else
           },
           "2"={ # return all MCMC all saved MCMC steps: incidence, DAGscore, orderscore and order and max result
             if (i==1) {
               result$chain<-MCMCresult
               attr(result$chain,"class")<-"MCMCtrace"
               return(result)
             } else {
               result$chain<-MCMCtraces
               attr(result$chain,"class")<-"MCMCtracemult"}
           },
           "3"={ # return max DAG, order, last search space incidence and all scoretables
             result$scoretable<-scoretable
           },
           "4"={ # return all MCMC all saved MCMC steps,max result,last search space and scoretables
             result$scoretable<-scoretable
             if(i==1) {
               result$chain<-MCMCresult
               attr(result$chain,"class")<-"MCMCtrace"
             } else {
               result$chain<-MCMCtraces
               attr(result$chain,"class")<-"MCMCtracemult"
             }
           }
    )
  
    return(result)
}
