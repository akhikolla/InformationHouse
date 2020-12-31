partitionscoreplus1<-function(n,nsmall,scorenodes,scorepositions,bgnodes,parenttable,aliases,scoretab,plus1lists,rowmapsneeded,
                              rowmapsallowed,scoresneeded,scoresallowed,permy,party,posy){
  
  partitionscores<-rep(0,n)
  scoresvec<-vector("numeric")
  allowedrow<-vector("numeric")
  scoretabrow<-vector("numeric")
  neededrow<-vector("numeric")
  m<-length(party)
  plus1neededlists<-list() #lists where score is taken from allowedscore, so only plus1 lists where top element is required
  plus1allowedlists<-list()  #lists where score is taken from neededscore, i.e. NOT_LISTS where top element is required
  k<-1
  for (i in scorenodes){
    tablesize<-dim(parenttable[[i]]) # just to remove some arguments
    position<-scorepositions[k]
    partyelement<-posy[position]
    if(partyelement==m){# no parents are allowed
      if(is.null(bgnodes)) {
        partitionscores[i]<-scoretab[[i]][1,1] #changed
        scoretabrow[i]<-1
        neededrow[i]<-0
        allowedrow[i]<-0 
      } else {
        neededrow[i]<--1
        allowednodes<-bgnodes
        bannednodes<-permy[1:nsmall] 
        bannedpool<-which(aliases[[i]]%in%bannednodes) 
        
        if(length(bannedpool)==tablesize[2]){ #no parents from the main table allowed
          allowedrow[i]<-nrow(parenttable[[i]]) #choose last row, so only plus1 parent is allowed and required
        }
        else { 
          allowedrow[i]<-rowmapsallowed[[i]]$backwards[sum(2^bannedpool)/2+1] 
        }
        plus1neededlists[[i]]<-c(1,which(plus1lists$parents[[i]]%in%allowednodes)+1)
        scoresvec<-scoresallowed[[i]][allowedrow[i],plus1neededlists[[i]]] 
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
    }
    else if(tablesize[1]==0) { #no parents in PC-table, 1 row in scoretables, then only plus1 lists with requirednodes are good
      nextpartyelement<-partyelement+1 #new
      lastbannedindex<-sum(party[1:partyelement]) #new
      requirednodes<-permy[(lastbannedindex+1):(lastbannedindex+party[nextpartyelement])] #new
#      allowedneeded<-intersect(requirednodes,aliases[[i]])
      allowedrow[i]<-1
      neededrow[i]<-0
      plus1neededlists[[i]]<-which(plus1lists$parents[[i]]%in%requirednodes)+1
      scoresvec<-scoresallowed[[i]][neededrow[i],plus1neededlists[[i]]]
      maxallowed<-max(scoresvec)
      partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
    } else  {
      lastbannedindex<-sum(party[1:partyelement]) #new
      bannednodes<-permy[1:lastbannedindex] #new
      nextpartyelement<-partyelement+1 #new
      lastreqindex<-lastbannedindex+party[nextpartyelement] #new
      requirednodes<-permy[(lastbannedindex+1):lastreqindex] #new
      if(lastreqindex==nsmall) { #new
        allowednotrequired<-bgnodes #new
      } else {
        allowednotrequired<-c(permy[(lastreqindex+1):nsmall], bgnodes) #new
      }
      bannedpool<-which(aliases[[i]]%in%bannednodes)
      allowedneeded<-which(aliases[[i]]%in%requirednodes)
      
      if(length(bannedpool)==tablesize[2]){ #no parents from the main table allowed
        neededrow[i]<-0 #in the main table no rows are allowed, so no score from neededscores is taken
        allowedrow[i]<-nrow(parenttable[[i]]) #choose last row, so only plus1 parent is allowed and required
        plus1neededlists[[i]]<-which(plus1lists$parents[[i]]%in%requirednodes)+1 #plus1 lists for nodes which are in required set
        scoresvec<-scoresallowed[[i]][allowedrow[i],plus1neededlists[[i]]]
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
      else if(length(bannedpool)<tablesize[2]&&length(allowedneeded)==0){ #no nodes from main tables are in required nodes set,
        #but some are in the allowed set
        neededrow[i]<-0 #in the main table no rows are needed, so no score from neededscores is taken
        allowedrow[i]<-rowmapsallowed[[i]]$backwards[sum(2^bannedpool)/2+1] 
        plus1neededlists[[i]]<-which(plus1lists$parents[[i]]%in%requirednodes)+1
        scoresvec<-scoresallowed[[i]][allowedrow[i],plus1neededlists[[i]]] 
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
      else if (length(aliases[[i]])==1){ #only 1 parent in possible parent set in the main table and it is required
        neededrow[i]<-1
        allowedrow[i]<-1
        plus1neededlists[[i]]<-c(which(plus1lists$parents[[i]]%in%requirednodes)+1)
        #top element does not enter the sum
        #first list is added because we do not want to add "no parents" case to the score
        plus1allowedlists[[i]]<-c(1,which(plus1lists$parents[[i]]%in%allowednotrequired)+1) 
        scoresvec<-scoresneeded[[i]][neededrow[i],plus1allowedlists[[i]]]
        scoresvec<-c(scoresvec,scoresallowed[[i]][allowedrow[i],plus1neededlists[[i]]])
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
      else {
        allowedrow[i]<-rowmapsallowed[[i]]$backwards[sum(2^bannedpool)/2+1]
        if (length(bannedpool)==0){
          neededrow[i]<-rowmapsneeded[[i]]$backwards[2*sum(3^allowedneeded)/3+1]
        } else {
          neededrow[i]<-rowmapsneeded[[i]]$backwards[sum(3^bannedpool)/3+2*sum(3^allowedneeded)/3+1]}
        plus1neededlists[[i]]<-which(plus1lists$parents[[i]]%in%requirednodes)+1
        plus1allowedlists[[i]]<-c(1,which(plus1lists$parents[[i]]%in%allowednotrequired)+1)
        scoresvec<-scoresneeded[[i]][neededrow[i],plus1allowedlists[[i]]]
        scoresvec<-c(scoresvec,scoresallowed[[i]][allowedrow[i],plus1neededlists[[i]]])
        maxallowed<-max(scoresvec)
        partitionscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
    }
    k<-k+1
  }
  
  scores<-list()
  scores$neededrow<-neededrow
  scores$allowedrow<-allowedrow
  scores$plus1neededlists<-plus1neededlists
  scores$plus1allowedlists<-plus1allowedlists
  scores$totscores<-partitionscores
  return(scores)
  
}



