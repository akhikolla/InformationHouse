#scores a single order base version (base neighbourhood)
orderscoreBase<-function(n,scorenodes,scorepositions,parenttable,aliases,numparents,rowmaps,
                         scoretable,scoresmatrices,permy) {
  
  orderscores<-vector("double",n)
  therows<-vector("integer",n)
  k<-1
  for (i in scorenodes){
    position<-scorepositions[k]
    if (position==n) {  #no parents allowed, i.e. only first row, only first list
      orderscores[i]<-scoretable[[i]][1,1]
      therows[i]<-c(2^numparents[i])
    }
    else {
      bannednodes<-permy[1:position]
      allowednodes<-permy[(position+1):n]
      bannedpool<-which(aliases[[i]]%in%bannednodes)
      if (numparents[i]==0||length(bannedpool)==0) { #all parents allowed or no parents in the parent table
        therows[i]<-c(1)
      }
      else {
        therows[i]<-rowmaps[[i]]$backwards[sum(2^bannedpool)/2+1]
      }
      orderscores[i]<-scoresmatrices[[i]][therows[i],1]
    }
    k<-k+1
  }
  scores<-list()
  scores$therow<-therows
  scores$totscores<-orderscores
  return(scores)
  
}

#scores a single order base version (plus1 neighbourhood)

orderscorePlus1<-function(n,scorenodes,scorepositions,parenttable,aliases,numparents,
                          rowmaps,plus1lists,scoretable,scoresmatrices,permy) {
  
  orderscores<-vector("double",n)
  allowedscorelists<-vector("list",n)
  therows<-vector("integer",n)
  k<-1
  for (i in scorenodes){
    position<-scorepositions[k]
    if (position==n) {  #no parents allowed, i.e. only first row, only first list
      orderscores[i]<-scoretable[[i]][[1]][1,1]
      allowedscorelists[[i]]<-c(1)
      therows[i]<-c(2^numparents[i])
    }
    else {
      bannednodes<-permy[1:position]
      allowednodes<-permy[(position+1):n]
      bannedpool<-which(aliases[[i]]%in%bannednodes)
      if (numparents[i]==0||length(bannedpool)==0) { #all parents allowed
        therows[i]<-c(1)
      }
      else {
        therows[i]<-rowmaps[[i]]$backwards[sum(2^bannedpool)/2+1]
      }
      allowedscorelists[[i]]<-c(1,which(plus1lists$parents[[i]]%in%allowednodes)+1)
      scoresvec<-scoresmatrices[[i]][therows[i],allowedscorelists[[i]]]
      maxallowed<-max(scoresvec)
      orderscores[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
    }
    k<-k+1
  }
  scores<-list()
  scores$therow<-therows
  scores$allowedlists<-allowedscorelists
  scores$totscores<-orderscores
  return(scores)
  
  
}

#returns an order and its BGe/BDe logscore such that is sampled from all orders (according to their scores), obtained via
#putting a given node in every position from 1 to n with all other nodes fixed  (plus1 neighbourhood)

positionscorePlus1<-function(n,nsmall,currentscore,positionx,permy,aliases,rowmaps,plus1lists,numparents,
                             scoretable,scoresmatrices) {
  vectorx<-vector(length=nsmall) #scores of node in each place in the order
  vectorall<-vector(length=nsmall) #scores for each node when x takes its position
  base<-currentscore$totscores
  totalall<-vector(length=nsmall) #result vector with log scores of all orders
  totalall[positionx]<-sum(base) #totallogscore of initial permutation
  x<-permy[positionx] #node for which we search max/sample position
  vectorx[positionx]<-base[x] #its score in current permy
  allowedlistx<-list()
  allowedlisty<-list()
  therowx<-vector()
  therowy<-vector()
  
  if (positionx>1) {
    rightpart<-permy[-positionx]
    for (i in (positionx-1):1) {
      
      nodey<-permy[i] #node which changes position with nodex
      
      #row with the new score of nodex
      if (numparents[x]==0||i==1) {
        therowx[i]<-c(1)
      } else {
        bannedpool<-which(aliases[[x]]%in%permy[1:(i-1)])
        if (length(bannedpool)==0) {
          therowx[i]<-c(1)
        } else {
          therowx[i]<-rowmaps[[x]]$backwards[sum(2^bannedpool)/2+1]
        }
      }
      
      allowedlistx[[i]]<-c(1,which(plus1lists$parents[[x]]%in%c(permy[(i+1):n],nodey))+1)
      scoresvec<-scoresmatrices[[x]][therowx[i],allowedlistx[[i]]]
      maxallowed<-max(scoresvec)
      vectorx[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      
      newpos<-i+1 #new position of nodey
      if(newpos==n) { #no parents allowed, i.e. only first row, only first list
        vectorall[nodey]<-scoretable[[nodey]][[1]][1,1]
        allowedlisty[[nodey]]<-c(1)
        therowy[nodey]<-c(2^numparents[nodey])
      } else {
        bannedpool<-which(aliases[[nodey]]%in%c(permy[1:(i-1)],x))
        if (numparents[nodey]==0||length(bannedpool)==0) {
          therowy[nodey]<-c(1)
        }
        else {
          therowy[nodey]<-rowmaps[[nodey]]$backwards[sum(2^bannedpool)/2+1]
        }
        if (newpos==n) {allowedlisty[[nodey]]<-c(1)} else {
          allowedlisty[[nodey]]<-c(1,which(plus1lists$parents[[nodey]]%in%rightpart[(newpos):(n-1)])+1)
        }
        
        scoresvec<-scoresmatrices[[nodey]][therowy[nodey],allowedlisty[[nodey]]]
        maxallowed<-max(scoresvec)
        vectorall[nodey]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
      
      totalall[i]<-totalall[i+1]-vectorx[i+1]+vectorx[i]-base[nodey]+vectorall[nodey]
    }
  }
  
  if (positionx < nsmall) {
    
    for (i in (positionx+1):nsmall) {
      nodey<-permy[i]
      
      if (numparents[x]==0) { #there is only 1 row in the score table
        therowx[i]<-c(1)
      } else if (i==n) { #no parents allowed, i.e. only first row, only first list
        vectorall[x]<-scoretable[[x]][[1]][1,1]
        therowx[i]<-c(2^numparents[x])
      } else {
        bannedpool<-which(aliases[[x]]%in%permy[1:i])
        if (length(bannedpool)==0) {
          therowx[i]<-c(1)
        } else {
          therowx[i]<-rowmaps[[x]]$backwards[sum(2^bannedpool)/2+1]
        }
      }
      if (i==n) { allowedlistx[[i]]<-c(1)} else {
        allowedlistx[[i]]<-c(1,which(plus1lists$parents[[x]]%in%permy[(i+1):n])+1)
      }
      
      scoresvec<-scoresmatrices[[x]][therowx[i],allowedlistx[[i]]]
      maxallowed<-max(scoresvec)
      vectorx[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      
      
      newpos<-i-1
      if (newpos==1) {
        therowy[nodey]<-c(1)
      } else {
        bannedpool<-which(aliases[[nodey]]%in%((permy[1:newpos])[-positionx]))
        if (numparents[nodey]==0||length(bannedpool)==0) {
          therowy[nodey]<-c(1)
        } else {
          therowy[nodey]<-rowmaps[[nodey]]$backwards[sum(2^bannedpool)/2+1]
        }
      }
      
      allowedlisty[[nodey]]<-c(1,which(plus1lists$parents[[nodey]]%in%c(permy[i:n],x))+1)
      scoresvec<-scoresmatrices[[nodey]][therowy[nodey],allowedlisty[[nodey]]]
      maxallowed<-max(scoresvec)
      vectorall[nodey]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      
      totalall[i]<-totalall[i-1]-vectorx[i-1]+vectorx[i]-base[nodey]+vectorall[nodey]
    }
  }
  
  totalmax<-max(totalall)
  allscore<-totalmax+log(sum(exp(totalall-totalmax)))
  maxi<-sample.int(nsmall,1,prob=exp(totalall-allscore))
  res<-list()
  
  if (maxi==positionx) {
    res$score<-currentscore
    res$order<-permy
    res$tot<-totalall[positionx]
    return(res)
  } else if (maxi>positionx) {
    newscore<-currentscore
    newscore$therow[x]<-therowx[maxi]
    newscore$allowedlists[[x]]<-allowedlistx[[maxi]]
    newscore$totscores[x]<-vectorx[maxi]
    updatenodes<-permy[(positionx+1):maxi]
    newscore$therow[updatenodes]<-therowy[updatenodes]
    newscore$allowedlists[updatenodes]<-allowedlisty[updatenodes]
    newscore$totscores[updatenodes]<-vectorall[updatenodes]
    res$score<-newscore
    res$order<-movenode(permy,positionx,maxi,n)
    res$tot<-totalall[maxi]
    return(res)
  } else {
    newscore<-currentscore
    newscore$therow[x]<-therowx[maxi]
    newscore$allowedlists[[x]]<-allowedlistx[[maxi]]
    newscore$totscores[x]<-vectorx[maxi]
    updatenodes<-permy[maxi:(positionx-1)]
    newscore$therow[updatenodes]<-therowy[updatenodes]
    newscore$allowedlists[updatenodes]<-allowedlisty[updatenodes]
    newscore$totscores[updatenodes]<-vectorall[updatenodes]
    res$score<-newscore
    res$order<-movenode(permy,positionx,maxi,n)
    res$tot<-totalall[maxi]
    return(res)
  }
  
  
}

#returns an order and its BGe/BDe logscore such that is sampled from all orders (according to their scores), obtained via
#putting a given node in every position from 1 to n with all other nodes fixed  (plus1 neighbourhood)

positionscorebase<-function(n,nsmall,currentscore,positionx,permy,aliases,rowmaps,numparents,
                            scoretable,scoresmatrices) {
  vectorx<-vector(length=nsmall) #scores of node in each place in the order
  vectorall<-vector(length=nsmall) #scores for each node when x takes its position
  base<-currentscore$totscores
  totalall<-vector(length=nsmall) #result vector with log scores of all orders
  totalall[positionx]<-sum(base) #totallogscore of initial permutation
  x<-permy[positionx] #node for which we search max/sample position
  vectorx[positionx]<-base[x] #its score in current permy
  therowx<-vector()
  therowy<-vector()
  
  if (positionx>1) {
    rightpart<-permy[-positionx]
    for (i in (positionx-1):1) {
      nodey<-permy[i]
      if (numparents[x]==0||i==1) {
        therowx[i]<-c(1)
      } else {
        bannedpool<-which(aliases[[x]]%in%permy[1:(i-1)])
        if (length(bannedpool)==0) {
          therowx[i]<-c(1)
        } else {
          therowx[i]<-rowmaps[[x]]$backwards[sum(2^bannedpool)/2+1]
        }
      }
      
      scoresvec<-scoresmatrices[[x]][therowx[i],1]
      maxallowed<-max(scoresvec)
      vectorx[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      
      newpos<-i+1 #new position of nodey
      if(newpos==n) { #no parents allowed, i.e. only first row, only first list
        vectorall[nodey]<-scoretable[[nodey]][1,1]
        therowy[nodey]<-c(2^numparents[nodey])
      } else {
        bannedpool<-which(aliases[[nodey]]%in%c(permy[1:(i-1)],x))
        if (numparents[nodey]==0||length(bannedpool)==0) {
          therowy[nodey]<-c(1)
        }
        else {
          therowy[nodey]<-rowmaps[[nodey]]$backwards[sum(2^bannedpool)/2+1]
        }
        
        scoresvec<-scoresmatrices[[nodey]][therowy[nodey],1]
        maxallowed<-max(scoresvec)
        vectorall[nodey]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      }
      
      totalall[i]<-totalall[i+1]-vectorx[i+1]+vectorx[i]-base[nodey]+vectorall[nodey]
    }
  }
  
  if (positionx < nsmall) {
    
    for (i in (positionx+1):nsmall) {
      nodey<-permy[i]
      
      if (numparents[x]==0) { #there is only 1 row in the score table
        therowx[i]<-c(1)
      } else if (i==n) { #no parents allowed, i.e. only first row, only first list
        vectorall[x]<-scoretable[[x]][1,1]
        therowx[i]<-c(2^numparents[x])
      } else {
        bannedpool<-which(aliases[[x]]%in%permy[1:i])
        if (length(bannedpool)==0) {
          therowx[i]<-c(1)
        } else {
          therowx[i]<-rowmaps[[x]]$backwards[sum(2^bannedpool)/2+1]
        }
      }
      
      scoresvec<-scoresmatrices[[x]][therowx[i],1]
      maxallowed<-max(scoresvec)
      vectorx[i]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      
      
      newpos<-i-1
      if (newpos==1) {
        therowy[nodey]<-c(1)
      } else {
        bannedpool<-which(aliases[[nodey]]%in%((permy[1:newpos])[-positionx]))
        if (numparents[nodey]==0||length(bannedpool)==0) {
          therowy[nodey]<-c(1)
        } else {
          therowy[nodey]<-rowmaps[[nodey]]$backwards[sum(2^bannedpool)/2+1]
        }
      }
      
      scoresvec<-scoresmatrices[[nodey]][therowy[nodey],1]
      maxallowed<-max(scoresvec)
      vectorall[nodey]<-maxallowed+log(sum(exp(scoresvec-maxallowed)))
      
      totalall[i]<-totalall[i-1]-vectorx[i-1]+vectorx[i]-base[nodey]+vectorall[nodey]
    }
  }
  
  totalmax<-max(totalall)
  allscore<-totalmax+log(sum(exp(totalall-totalmax)))
  maxi<-sample.int(nsmall,1,prob=exp(totalall-allscore))
  res<-list()
  
  if (maxi==positionx) {
    res$score<-currentscore
    res$order<-permy
    res$tot<-totalall[positionx]
    return(res)
  } else if (maxi>positionx) {
    newscore<-currentscore
    newscore$therow[x]<-therowx[maxi]
    newscore$totscores[x]<-vectorx[maxi]
    updatenodes<-permy[(positionx+1):maxi]
    newscore$therow[updatenodes]<-therowy[updatenodes]
    newscore$totscores[updatenodes]<-vectorall[updatenodes]
    res$score<-newscore
    res$order<-movenode(permy,positionx,maxi,n)
    res$tot<-totalall[maxi]
    return(res)
  } else {
    newscore<-currentscore
    newscore$therow[x]<-therowx[maxi]
    newscore$totscores[x]<-vectorx[maxi]
    updatenodes<-permy[maxi:(positionx-1)]
    newscore$therow[updatenodes]<-therowy[updatenodes]
    newscore$totscores[updatenodes]<-vectorall[updatenodes]
    res$score<-newscore
    res$order<-movenode(permy,positionx,maxi,n)
    res$tot<-totalall[maxi]
    return(res)
  }
  
}
