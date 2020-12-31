#implements order MCMC algorithm with a defined searchspace, MAP version

orderMCMCbasemax<-function(n,nsmall,startorder,iterations,stepsave,moveprobs,parenttable,scoretable,aliases,numparents,
                           rowmaps,maxmatrices,numberofparentsvec,gamma=1,bgnodes,matsize) {
  
  #n - number of nodes (background included)
  #nsmall - number of nodes excluding background
  #matsize - number of rows/columns in adjacency matrix 
  
  if(!is.null(bgnodes)) {
    mainnodes<-c(1:matsize)[-bgnodes]
  } else {mainnodes<-c(1:matsize)}
  
  currentpermy<-startorder #starting order represented as a permutation
  currentorderscores<-orderscoreBasemax(matsize,startorder[1:nsmall],c(1:nsmall),parenttable,aliases,numparents,
                                        rowmaps,scoretable,maxmatrices,currentpermy) #starting score
  currenttotallogscore<-sum(currentorderscores$totscores[mainnodes]) #log total score of all DAGs in the starting order
  currentDAG<-samplescoreplus1.max(matsize=matsize,samplenodes=mainnodes,scores=currentorderscores,plus1lists=NULL,
                                   maxmatrices=maxmatrices,scoretable=scoretable,
                                   parenttable,numberofparentsvec,aliases) #score of a single DAG sampled from the starting order
  L1 <- list() # stores the adjacency matrix of a DAG sampled from the orders
  L2 <- list() # stores its log BGe score
  L3 <- list() # stores the log BGe score of the entire order
  L4 <- list() # stores the orders as permutations
  
  zlimit<- floor(iterations/stepsave) + 1 # number of outer iterations
  length(L1) <- zlimit
  length(L2) <- zlimit
  length(L3) <- zlimit
  length(L4) <- zlimit
  
  L1[[1]]<-currentDAG$incidence #starting DAG adjacency matrix
  L2[[1]]<-currentDAG$logscore #starting DAG score
  L3[[1]]<-currenttotallogscore #starting order score
  L4[[1]]<-currentpermy[1:nsmall] #starting order
  
  moveprobsstart<-moveprobs
  
  for (z in 2:zlimit){ #the MCMC chain loop with 'iteration' steps is in two parts
    for (count in 1:stepsave){ #since we only save the results to the lists each 'stepsave'
      
      
      chosenmove<-sample.int(4,1,prob=moveprobs)
      if(chosenmove<4){  # if it is 3 then we stay still
        
        proposedpermy<-currentpermy #sample a new order by swapping two elements
        switch(as.character(chosenmove),
               "1"={ # swap any two elements at random
                 sampledelements<-sample.int(nsmall,2,replace=FALSE) #chosen at random
               },
               "2"={ # swap any adjacent elements
                 k<-sample.int(nsmall-1,1) #chose the smallest at random
                 sampledelements<-c(k,k+1)
               },
               "3"={ # find best position for 1 element
                 sampledpos<-sample.int(nsmall,1)
               },
               "4"={ # stay still
               },
               {# if neither is chosen, we have a problem
                 print('The move sampling has failed!')
               }) #end switch
        
        if(chosenmove<3){
          scorepositions<-c(min(sampledelements):max(sampledelements))
          proposedpermy[sampledelements]<-currentpermy[rev(sampledelements)] #proposed new order
          rescorenodes<-proposedpermy[scorepositions] #we only need to rescore these nodes between the swapped elements to speed up the calculation
          proposedorderrescored<-orderscoreBasemax(matsize,rescorenodes,scorepositions,parenttable,aliases,numparents,rowmaps,scoretable,maxmatrices,proposedpermy)
          proposedtotallogscore<-currenttotallogscore-sum(currentorderscores$totscores[rescorenodes])+sum(proposedorderrescored$totscores[rescorenodes]) #and the new log total score by updating only the necessary nodes
          scoreratio<-exp((proposedtotallogscore-currenttotallogscore)*gamma) #acceptance probability
          
          if(runif(1)<scoreratio){ #Move accepted then set the current order and scores to the proposal
            currentpermy<-proposedpermy
            currentorderscores$therow[rescorenodes]<-proposedorderrescored$therow[rescorenodes]
            currentorderscores$totscores[rescorenodes]<-proposedorderrescored$totscores[rescorenodes]
            currenttotallogscore<-proposedtotallogscore
          } 
        } else {
           neworder<-positionscorebasemax(matsize,nsmall,currentorderscores,sampledpos,currentpermy,aliases,
                                          rowmaps,numparents,scoretable,maxmatrices)
           currentorderscores<-neworder$score
           currentpermy<-neworder$order
           currenttotallogscore<-sum(currentorderscores$totscores)
        }
      }
    }
    currentDAG<-samplescoreplus1.max(n,mainnodes,currentorderscores,plus1lists=NULL,maxmatrices,scoretable,
                                     parenttable,numberofparentsvec,aliases)
    L1[[z]]<-currentDAG$incidence #store adjacency matrix of a sampled DAG each 'stepsave'
    L2[[z]]<-currentDAG$logscore #and log score of a sampled DAG
    L3[[z]]<-currenttotallogscore #and the current order score
    L4[[z]]<-currentpermy[1:nsmall] #and store current order
  }
  result<-list()
  result$incidence<-L1
  result$DAGscores<-L2
  result$orderscores<-L3
  result$orders<-L4
  return(result)
}

