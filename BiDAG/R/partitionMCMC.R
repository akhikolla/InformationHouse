#implements partition MCMC scheme for structure learning problem, searches in plus1 neighbourhood of a search space defined
#by aliases and scoretables 

partitionMCMCplus1<-function(n,nsmall,startpermy,startparty,iterations,stepsave,parenttable,scoretable,scoretab,
                             aliases,scoresneeded,scoresallowed,plus1lists, rowmapsneeded,rowmapsallowed,
                             needednodetable,numberofparentsvec,numberofpartitionparentsvec,needednodebannedrow,
                             neededparentsvec,moveprobs,bgnodes,matsize) {
  
  #n - number of nodes (background included)
  #nsmall - number of nodes excluding background
  #matsize - number of rows/columns in adjacency matrix 
  
  if(!is.null(bgnodes)) {
    mainnodes<-c(1:matsize)[-bgnodes]
  } else {mainnodes<-c(1:matsize)}
  
  
  currentpermy<-startpermy[1:nsmall] #starting permutation
  currentparty<-startparty #starting partition
  currentposy<-parttolist(nsmall,currentparty) #create a list of which nodes are in which partition element

  currentpartitionscores<-partitionscoreplus1(n,nsmall,currentpermy,c(1:nsmall),bgnodes,parenttable,aliases,scoretab,plus1lists,rowmapsneeded,
                                              rowmapsallowed,scoresneeded,scoresallowed,
                                              currentpermy,currentparty,currentposy) #starting score of all DAGs compatible with the starting permutation and partition

  currenttotallogscore<-sum(currentpartitionscores$totscores[mainnodes]) #log total score of all DAGs in the starting partition and permutation

  currentDAG<-samplescore.partition.plus1(matsize,mainnodes,currentpartitionscores,scoretable,
                                          scoresallowed,scoresneeded,scoretab,
                                          parenttable,needednodetable,
                                          numberofparentsvec,needednodebannedrow,
                                          numberofpartitionparentsvec,plus1lists) #log score of a single sampled DAG from the partition and permutation

  L1 <- list() #stores the adjacency matrix of a DAG sampled from the partition and permutation
  L2 <- list() #stores the log BGe score of a DAG sampled from the partition and permutation
  L3 <- list() #stores the log BGe score of the entire partition following the permutation
  L4 <- list() #stores the permutations
  L5 <- list() #stores the partitions
  zlimit<- floor(iterations/stepsave) + 1 # number of outer iterations
  length(L1) <- zlimit
  length(L2) <- zlimit
  length(L3) <- zlimit
  length(L4) <- zlimit
  length(L5) <- zlimit

  L1[[1]]<-currentDAG$incidence #starting DAG adjacency matrix
  L2[[1]]<-currentDAG$logscore #starting DAG score
  L3[[1]]<-currenttotallogscore #starting partition score
  L4[[1]]<-currentpermy #starting permutation
  L5[[1]]<-currentparty #starting partition

  # Set some flags for when we need to recalculate neighbourhoods
  permdiffelemflag<-1
  permneighbourflag<-1

  partstepflag<-1
  partjoinholeflag<-1

  for (z in 2:zlimit){ #the MCMC chain loop with 'iteration' steps is in two parts
    count<-1
    while (count <=stepsave){ #since we only save the results to the lists each 'stepsave'
        chosenmove<-sample.int(5,1,prob=moveprobs) # sample what type of move
      if(chosenmove<3){  # if it is <3 then we swap two elements

        if(length(currentparty)>1){ # if the partition only has one element then we cannot move
          switch(as.character(chosenmove),
                 "1"={ # swap any two elements from diffent partition elements
                   if(permdiffelemflag>0){ # do we need to recalculate the neighbourhood?
                     permdiffelemposs<-parttopermdiffelemposs(nsmall,currentparty)
                     permdiffelemflag<-0
                   }
                   temp<-swapdiffelementnodes(nsmall,currentparty,currentposy,currentpermy,permdiffelemposs)
                   proposedpermy<-temp[[1]]
                   rescorenodes<-temp[[2]]
                   scorepositions<-temp[[3]]
                 },
                 "2"={ # swap any elements in adjacent partition elements
                   if(permneighbourflag>0){ # do we need to recalculate the neighbourhood?
                     permneighbourposs<-parttopermneighbourposs(nsmall,currentparty)
                     permneighbourflag<-0
                   }
                   temp<-swapadjacentnodes(nsmall,currentparty,currentposy,currentpermy,permneighbourposs)
                   proposedpermy<-temp[[1]]
                   rescorenodes<-temp[[2]]
                   scorepositions<-temp[[3]]
                 },
                 {# if neither is chosen, we have a problem
                   print('The move sampling has failed!')
                 })

          proposedpartitionrescored<-partitionscoreplus1(n,nsmall,rescorenodes,scorepositions,bgnodes,parenttable,aliases,scoretab,plus1lists,rowmapsneeded,rowmapsallowed,scoresneeded,scoresallowed,
                                                         proposedpermy,currentparty,currentposy) #their scores
            proposedtotallogscore<-currenttotallogscore-sum(currentpartitionscores$totscores[rescorenodes])+sum(proposedpartitionrescored$totscores[rescorenodes]) #and the new log total score by updating only the necessary nodes
            scoreratio<-exp(proposedtotallogscore-currenttotallogscore) #acceptance probability
            count<-count+1
            if(runif(1)<scoreratio){ #Move accepted then set the current permutation and scores to the proposal
              currentpermy<-proposedpermy
              currentpartitionscores$neededrow[rescorenodes]<-proposedpartitionrescored$neededrow[rescorenodes]
              currentpartitionscores$allowedrow[rescorenodes]<-proposedpartitionrescored$allowedrow[rescorenodes]
              currentpartitionscores$plus1allowedlists[rescorenodes]<-proposedpartitionrescored$plus1allowedlists[rescorenodes]
              currentpartitionscores$plus1neededlists[rescorenodes]<-proposedpartitionrescored$plus1neededlists[rescorenodes]
              currentpartitionscores$totscores[rescorenodes]<-proposedpartitionrescored$totscores[rescorenodes]
              currenttotallogscore<-proposedtotallogscore
            }
          }
      } else if(chosenmove<5) { # we move in the space of partitions

        switch(as.character(chosenmove),
               "3"={ # we split a partition element or join one
                 if(partstepflag>0){ # do we need to recalculate the neighbourhood?
                   currentpartstepposs<-partysteps(nsmall,currentparty)
                   currentpartstepnbhood<-sum(currentpartstepposs)
                   partstepflag<-0
                 }
                 temp<-partitionsplitorjoin(nsmall,currentparty,currentposy,currentpermy,currentpartstepposs)
                 proposedparty<-temp[[1]]
                 proposedposy<-temp[[2]]
                 proposedpermy<-temp[[3]]
                 rescorenodes<-temp[[4]]
                 proposedpartstepposs<-temp[[5]]
                 scorepositions<-temp[[6]]
                 proposedpartstepnbhood<-sum(proposedpartstepposs)
               },
               "4"={ # we move a single node into another partition element or into a new one
                 if(partjoinholeflag>0){ # do we need to recalculate the neighbourhood?
                   currentpartjoinposs<-partyjoin(nsmall,currentparty,currentposy)
                   currentpartjoinnbhood<-sum(currentpartjoinposs)
                   currentpartholeposs<-partyhole(nsmall,currentparty,currentposy)
                   currentpartholenbhood<-sum(currentpartholeposs)
                   partjoinholeflag<-0
                 }

                 joinorhole<-sample.int(2,1,prob=c(currentpartjoinnbhood,currentpartholenbhood)) # choose the type of move

                 switch(as.character(joinorhole),
                        "1"={ # we join the node to another partition element
                          temp<-joinnode(nsmall,currentparty,currentposy,currentpermy,currentpartjoinposs)
                        },
                        "2"={ # we place the node in a new partition element
                          temp<-holenode(nsmall,currentparty,currentposy,currentpermy,currentpartholeposs)
                        },
                        {# if nothing is chosen, we have a problem
                          print('The move sampling has failed!')
                        })

                 proposedparty<-temp[[1]]
                 proposedposy<-temp[[2]]
                 proposedpermy<-temp[[3]]
                 rescorenodes<-temp[[4]]
                 scorepositions<-temp[[5]]

                 # these neighbourhoods should be updated for efficiency
                 proposedpartjoinposs<-partyjoin(nsmall,proposedparty,proposedposy)
                 proposedpartjoinnbhood<-sum(proposedpartjoinposs)
                 proposedpartholeposs<-partyhole(nsmall,proposedparty,proposedposy)
                 proposedpartholenbhood<-sum(proposedpartholeposs)
               },
               {# if nothing is chosen, we have a problem
                 print('The move sampling has failed!')
               })

        proposedpartitionrescored<-partitionscoreplus1(n,nsmall,rescorenodes,scorepositions,bgnodes,parenttable,aliases,scoretab,plus1lists,rowmapsneeded,rowmapsallowed,scoresneeded,scoresallowed,
                                                       proposedpermy,proposedparty,proposedposy) #only rescore the necessary nodes
        proposedtotallogscore<-currenttotallogscore-sum(currentpartitionscores$totscores[rescorenodes])+sum(proposedpartitionrescored$totscores[rescorenodes]) #and calculate the new log total score by updating only the necessary nodes
        count<-count+1
        scoreratio<-exp(proposedtotallogscore-currenttotallogscore)

        switch(as.character(chosenmove),
               "3"={ # we split a partition element or joined one
                 scoreratio<-scoreratio*(currentpartstepnbhood/proposedpartstepnbhood) # neighbourhood correction
               },
               "4"={# we moved a single node
                 scoreratio<-scoreratio*((currentpartjoinnbhood+currentpartholenbhood)/(proposedpartjoinnbhood+proposedpartholenbhood)) # neighbourhood correction
               },
               {# if nothing is chosen, we have a problem
                 print('The move sampling has failed!')
               })

        if(runif(1)<scoreratio){ #Move accepted then set the current partition and scores to the proposal
          currentpermy<-proposedpermy
          currentparty<-proposedparty
          currentposy<-proposedposy
          currentpartitionscores$neededrow[rescorenodes]<-proposedpartitionrescored$neededrow[rescorenodes]
          currentpartitionscores$allowedrow[rescorenodes]<-proposedpartitionrescored$allowedrow[rescorenodes]
          currentpartitionscores$plus1allowedlists[rescorenodes]<-proposedpartitionrescored$plus1allowedlists[rescorenodes]
          currentpartitionscores$plus1neededlists[rescorenodes]<-proposedpartitionrescored$plus1neededlists[rescorenodes]
          currentpartitionscores$totscores[rescorenodes]<-proposedpartitionrescored$totscores[rescorenodes]

          currenttotallogscore<-proposedtotallogscore

          permdiffelemflag<-1 # need to recalculate the permutation possibilities
          permneighbourflag<-1 # in principle these could be updated instead

          switch(as.character(chosenmove),
                 "3"={ # we split a partition element or joined one
                   partjoinholeflag<-1
                   currentpartstepposs<-proposedpartstepposs
                   currentpartstepnbhood<-proposedpartstepnbhood
                 },
                 "4"={# we made a different partition move?
                   partstepflag<-1
                   currentpartjoinposs<-proposedpartjoinposs
                   currentpartjoinnbhood<-proposedpartjoinnbhood
                   currentpartholeposs<-proposedpartholeposs
                   currentpartholenbhood<-proposedpartholenbhood
                 },
                 {# if nothing is chosen, we have a problem
                   print('The move sampling has failed!')
                 })
        }
        }
    }
    currentDAG<-samplescore.partition.plus1(matsize,mainnodes,currentpartitionscores,
                                            scoretable,scoresallowed,scoresneeded,scoretab,parenttable,needednodetable,numberofparentsvec,
                                            needednodebannedrow,numberofpartitionparentsvec,plus1lists)
    L1[[z]]<-currentDAG$incidence #store adjacency matrix of a sampled DAG
    L2[[z]]<-currentDAG$logscore #and its log score
    L3[[z]]<-currenttotallogscore #store the current total partition score
    L4[[z]]<-currentpermy[1:nsmall] #store current permutation each 'stepsave'
    L5[[z]]<-currentparty #store current partition each 'stepsave'
  }

  result<-list()
  result$incidence<-L1
  result$DAGscores<-L2
  result$partitionscores<-L3
  result$order<-L4
  result$partition<-L5
  return(result)
}
