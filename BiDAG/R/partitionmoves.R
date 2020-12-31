# swap only nodes in different partition elements # move 1

swapdiffelementnodes<-function(n,currentparty,currentposy,currentpermy,permpossibs){
  
  scorelength<-length(permpossibs) 
  selectedelement<-sample.int(scorelength,1,prob=permpossibs) # sample an element
  
  leftnodes<-which(currentposy==selectedelement) #nodes in the left partition element
  maxleftnodes<-max(leftnodes)
  sampledleftnode<-propersample(leftnodes)
  sampledrightnode<-propersample(c((maxleftnodes+1):n)) #chose the right node from the remaining ones
  sampledelements<-c(sampledleftnode,sampledrightnode) #the sampled pair
  rightnodes<-which(currentposy==currentposy[sampledrightnode]) #the remaining nodes on the right
  minrightnodes<-min(rightnodes)
  
  centralnodes<-c()
  if((maxleftnodes+1)<minrightnodes){ #the nodes inbetween
    centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
  }
  scorepositions<-c(leftnodes,centralnodes,sampledrightnode) #all these need to be rescored
  if(selectedelement>1){ #as well as the partition element further left, if it exists
    scorepositions<-c(which(currentposy==(selectedelement-1)),scorepositions)
  }
  
  proposedpermy<-currentpermy #create the new permutation
  proposedpermy[sampledelements]<-currentpermy[rev(sampledelements)] #by swapping the elements
  
  rescorenodes<-proposedpermy[scorepositions] #map to the correct labels
  
  return(list(proposedpermy,rescorenodes,scorepositions))
  
}

# swap any two adjacent nodes # move 2

swapadjacentnodes<-function(n,currentparty,currentposy,currentpermy,permpossibs){
  
  scorelength<-length(permpossibs) 
  selectedelement<-sample.int(scorelength,1,prob=permpossibs) # sample an element
  
  leftnodes<-which(currentposy==selectedelement) #nodes in the left partition element
  rightnodes<-which(currentposy==(selectedelement+1)) #and in the right
  sampledelements<-c(propersample(leftnodes),propersample(rightnodes)) #sample a pair
  
  scorepositions<-c(leftnodes,sampledelements[2]) #need to rescore the left partition element and just the selected node in the right partition
  
  if(selectedelement>1){ #as well as the partition element further left, if it exists
    scorepositions<-c(which(currentposy==(selectedelement-1)),scorepositions)
  }
  
  proposedpermy<-currentpermy #create the new permutation
  proposedpermy[sampledelements]<-currentpermy[rev(sampledelements)] #by swapping the elements
  
  rescorenodes<-proposedpermy[scorepositions] #map to the correct labels
  
  return(list(proposedpermy,rescorenodes,scorepositions))
  
}

# Split or join adjacent partition elements # move 3

partitionsplitorjoin<-function(n,currentparty,currentposy,currentpermy,partpossibs){
  
  scorelength<-length(partpossibs) 
  sampledelement<-sample.int(scorelength,1,prob=partpossibs) # sample an element
  
  # vectors to be updated later
  proposedposy<-currentposy
  proposedpartposs<-partpossibs
  
  if(currentposy[sampledelement]==currentposy[sampledelement+1]){ #if the next node is in the same partition element we split
    involvednodes<-which(currentposy==currentposy[sampledelement]) #the nodes in the partition element to be split
    newleftelementsize<-sampledelement-min(involvednodes)+1 #the size of the new partition element on the left
    newrightelementsize<-length(involvednodes)-newleftelementsize #and the right
    proposedparty<-append(currentparty,newrightelementsize,currentposy[sampledelement]) #create the proposed partition
    proposedparty[currentposy[sampledelement]]<-newleftelementsize
    # update the list to be slightly more efficient
    proposedposy[(sampledelement+1):n]<-currentposy[(sampledelement+1):n]+1
    
    # sample which nodes go to the new left and right partition elements # note that involvednodes must have at least two elements so `sample' should work fine
    newleftnodes<-sample(involvednodes,newleftelementsize) 
    newrightnodes<-involvednodes[-which(involvednodes%in%newleftnodes)]
    proposedpermy<-currentpermy #and update the proposed permutation according to how the nodes are distributed 
    proposedpermy[(sampledelement-newleftelementsize+1):sampledelement]<-currentpermy[newleftnodes]
    proposedpermy[(sampledelement+1):(sampledelement+newrightelementsize)]<-currentpermy[newrightnodes]
    
    # update the neighbourhood of the new left element
    proposedpartposs[(sampledelement-newleftelementsize+1):sampledelement]<-choose(newleftelementsize,1:newleftelementsize)
    
    scorepositions<-c((sampledelement-newleftelementsize+1):sampledelement) #we only need to rescore the nodes in the new left partition element
    if(currentposy[sampledelement]>1){ # and the partition element to the left, if it exists
      scorepositions<-c(which(currentposy==(currentposy[sampledelement]-1)),scorepositions)
    }
    
    if(currentposy[sampledelement]<length(currentparty)){ #update the neighbourhood of the right element
      proposedpartposs[(sampledelement+1):(sampledelement+newrightelementsize)]<-choose(newrightelementsize,1:newrightelementsize)
    } else { # if the partition element was last then we cut the binomial coeffiecients short
      if(newrightelementsize>1){
        proposedpartposs[(sampledelement+1):(sampledelement+newrightelementsize-1)]<-choose(newrightelementsize,1:(newrightelementsize-1))
      }
    }
    rescorenodes<-proposedpermy[scorepositions] #finally map to the correct labels
  }
  
  if(currentposy[sampledelement]<currentposy[sampledelement+1]){ #otherwise join partition elements
    newelementsize<-currentparty[currentposy[sampledelement]]+currentparty[currentposy[sampledelement+1]]
    proposedparty<-currentparty[-currentposy[sampledelement]] #create the proposed partition
    proposedparty[currentposy[sampledelement]]<-newelementsize
    # update the list to be slightly more efficient 
    proposedposy[(sampledelement+1):n]<-currentposy[(sampledelement+1):n]-1
    proposedpermy<-currentpermy #now changing the permutation has no effect
    
    scorepositions<-which(currentposy==currentposy[sampledelement]) #we only need to rescore the nodes from the old element on the left and the element further left
    newelementstart<-min(scorepositions) #store the start of the new element
    if(proposedposy[sampledelement]>1){ #include the partition element to the left, if it exists
      scorepositions<-c(which(proposedposy==(proposedposy[sampledelement]-1)),scorepositions)
    }
    if(proposedposy[sampledelement]<length(proposedparty)){ #update the neighbourhood 
      proposedpartposs[(newelementstart):(newelementstart+newelementsize-1)]<-choose(newelementsize,1:newelementsize)
    } else { # if the new element is last we cut the binomial coeffiecients short
      proposedpartposs[(newelementstart):(newelementstart+newelementsize-2)]<-choose(newelementsize,1:(newelementsize-1)) # the new element must have a size of at least two
    }
    rescorenodes<-proposedpermy[scorepositions] #again map to the correct labels
  }
  
  return(list(proposedparty,proposedposy,proposedpermy,rescorenodes,proposedpartposs,
              scorepositions))
  
}

# move a node to a different partition element (following further restrictions) # move 4 part 1

joinnode<-function(n,currentparty,currentposy,currentpermy,joinpossibs){
  m<-length(currentparty)
  scorelength<-length(joinpossibs) 
  nodetomove<-sample.int(scorelength,1,prob=joinpossibs) # sample an element
  
  nodeselement<-currentposy[nodetomove]
  elementtomoveto<-sample.int(joinpossibs[nodetomove],1) # sample where to go to
  
  if(elementtomoveto==nodeselement){ #if element selected is the same, replace by m which is not otherwise sampled
    elementtomoveto<-m
  } else if (elementtomoveto==(nodeselement+1)){
    if(currentparty[nodeselement]==1){ #if partition element has size one
      if(currentparty[nodeselement+1]==1){ #if next partition element is also size one
        elementtomoveto<-m-1 #replace not allowed sample by 'm-1' which also not be otherwise sampled
      }
    }
  }
  
  proposedparty<-currentparty # update the partition
  proposedparty[elementtomoveto]<-currentparty[elementtomoveto]+1
  if(currentparty[nodeselement]>1){
    proposedparty[nodeselement]<-currentparty[nodeselement]-1
  } else {
    proposedparty<-proposedparty[-nodeselement]
  }
  proposedposy<-parttolist(n,proposedparty) #should be updated for efficiency
  
  if(elementtomoveto<nodeselement){
    leftnodes<-which(currentposy==elementtomoveto)
    newleftnodes<-which(proposedposy==elementtomoveto)
    maxleftnodes<-max(leftnodes)
    rightnodes<-which(currentposy==nodeselement) # this includes the node which is moved
    minrightnodes<-min(rightnodes)
    proposedpermy<-currentpermy[-nodetomove]
    proposedpermy<-append(proposedpermy,currentpermy[nodetomove],maxleftnodes)
    centralnodes<-c()
    newcentralnodes<-c()
    
    if((maxleftnodes+1)<minrightnodes){
      centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      newcentralnodes<-centralnodes+1
    }
    newscorepositions<-c(newleftnodes,newcentralnodes) #all these need to be rescored
    if(elementtomoveto>1){ #as well as the partition element further left, if it exists
      newscorepositions<-c(which(currentposy==(elementtomoveto-1)),newscorepositions)
    }
  } else {
    leftnodes<-which(currentposy==nodeselement) # this includes the node which is moved
    if(length(leftnodes)>1) {
      newleftnodes<-which(proposedposy==nodeselement)
    } else {
      newleftnodes<-c()
    }
    maxleftnodes<-max(leftnodes)
    rightnodes<-which(currentposy==elementtomoveto)
    minrightnodes<-min(rightnodes)
    proposedpermy<-currentpermy[-nodetomove]
    newnodepos<-min(rightnodes)-1
    proposedpermy<-append(proposedpermy,currentpermy[nodetomove],newnodepos-1)
    centralnodes<-c()
    newcentralnodes<-c()
    
    if((maxleftnodes+1)<minrightnodes){
      centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      newcentralnodes<-centralnodes-1
    }
    newscorepositions<-c(newleftnodes,newcentralnodes,newnodepos)
    if(nodeselement>1){ #as well as the partition element further left, if it exists
      newscorepositions<-c(which(currentposy==(nodeselement-1)),newscorepositions)
    }
  }
  
  newrescorenodes<-proposedpermy[newscorepositions]

  return(list(proposedparty,proposedposy,proposedpermy,newrescorenodes,newscorepositions))
  
}


# move a node to a new partition element (following further restrictions) # move 4 part 2

holenode<-function(n,currentparty,currentposy,currentpermy,holepossibs){
  m<-length(currentparty)
  
  scorelength<-length(holepossibs) 
  nodetomove<-sample.int(scorelength,1,prob=holepossibs) # sample an element
  
  nodeselement<-currentposy[nodetomove]
  holetomoveto<-sample.int(holepossibs[nodetomove],1) # sample where to go to
  
  if(currentparty[nodeselement]==1){
    if(holetomoveto>=nodeselement){
      holetomoveto<-holetomoveto+1 #counts hole when partition element removed
      if(nodeselement<m){
        if(currentparty[nodeselement+1]==1){ #if next partition element is also size one
          holetomoveto<-holetomoveto+1 #shift one hole further along
        }
      }
    }
    proposedparty<-currentparty[-nodeselement]
    proposedparty<-append(proposedparty,1,holetomoveto-1)
    proposedposy<-parttolist(n,proposedparty) #should be updated for efficiency
    
    if(holetomoveto<nodeselement){
      leftnodes<-which(proposedposy==holetomoveto) # this should be the node which is moved
      maxleftnodes<-max(leftnodes)
      rightnodes<-which(proposedposy==nodeselement) # these were to the left of the moved node
      minrightnodes<-min(rightnodes)
      proposedpermy<-currentpermy[-nodetomove]
      proposedpermy<-append(proposedpermy,currentpermy[nodetomove],leftnodes-1)
      
      centralnodes<-c()
      if((maxleftnodes+1)<minrightnodes){
        centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      }
      scorepositions<-c(leftnodes,centralnodes,rightnodes) #all these need to be rescored
      if(holetomoveto>1){ #as well as the partition element further left, if it exists
        scorepositions<-c(which(currentposy==(holetomoveto-1)),scorepositions)
      }
    } else {
      leftnodes<-which(proposedposy==nodeselement) # these were to the right of the moved node
      maxleftnodes<-max(leftnodes)
      rightnodes<-which(proposedposy==holetomoveto) # this should be the node that is moved
      minrightnodes<-min(rightnodes)
      proposedpermy<-currentpermy[-nodetomove]
      proposedpermy<-append(proposedpermy,currentpermy[nodetomove],rightnodes-1)
      
      centralnodes<-c()
      if((maxleftnodes+1)<minrightnodes){
        centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      }
      scorepositions<-c(leftnodes,centralnodes,rightnodes) #all these need to be rescored
      if(nodeselement>1){ #as well as the partition element further left, if it exists
        scorepositions<-c(which(currentposy==(nodeselement-1)),scorepositions)
      }
    }
  }  else {
    if(currentparty[nodeselement]==2){
      if(holetomoveto==nodeselement){
        holetomoveto<-m+1
      }
    }
    proposedparty<-currentparty
    proposedparty[nodeselement]<-currentparty[nodeselement]-1
    proposedparty<-append(proposedparty,1,holetomoveto-1)
    proposedposy<-parttolist(n,proposedparty) #should be updated for efficiency
    
    if(holetomoveto<=nodeselement){
      leftnodes<-which(proposedposy==holetomoveto) # this is the node which is moved
      maxleftnodes<-max(leftnodes)
      rightnodes<-which(proposedposy==(nodeselement+1)) # these are the nodes left over
      minrightnodes<-min(rightnodes)
      proposedpermy<-currentpermy[-nodetomove]
      proposedpermy<-append(proposedpermy,currentpermy[nodetomove],leftnodes-1)
      
      centralnodes<-c()
      if((maxleftnodes+1)<minrightnodes){
        centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      }
      scorepositions<-c(leftnodes,centralnodes) #all these need to be rescored
      if(holetomoveto>1){ #as well as the partition element further left, if it exists
        scorepositions<-c(which(currentposy==(holetomoveto-1)),scorepositions)
      }
    } else {
      leftnodes<-which(proposedposy==nodeselement) # these are the nodes left over
      maxleftnodes<-max(leftnodes)
      rightnodes<-which(proposedposy==holetomoveto) # this is the node which is moved
      minrightnodes<-min(rightnodes)
      proposedpermy<-currentpermy[-nodetomove]
      proposedpermy<-append(proposedpermy,currentpermy[nodetomove],rightnodes-1)
      
      centralnodes<-c()
      if((maxleftnodes+1)<minrightnodes){
        centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      }
      scorepositions<-c(leftnodes,centralnodes,rightnodes) #all these need to be rescored
      if(nodeselement>1){ #as well as the partition element further left, if it exists
        scorepositions<-c(which(currentposy==(nodeselement-1)),scorepositions)
      }
    }
    
  }
  
  rescorenodes<-proposedpermy[scorepositions] #map to the correct labels
  
  return(list(proposedparty,proposedposy,proposedpermy,rescorenodes,scorepositions))
  
}








# swap only nodes in different partition elements # move 1
swapdiffelementnodes.old<-function(n,currentparty,currentposy,currentpermy,permpossibs){
  
  scorelength<-length(permpossibs) 
  selectedelement<-sample.int(scorelength,1,prob=permpossibs) # sample an element
  
  leftnodes<-which(currentposy==selectedelement) #nodes in the left partition element
  maxleftnodes<-max(leftnodes)
  sampledleftnode<-propersample(leftnodes)
  sampledrightnode<-propersample(c((maxleftnodes+1):n)) #chose the right node from the remaining ones
  sampledelements<-c(sampledleftnode,sampledrightnode) #the sampled pair
  rightnodes<-which(currentposy==currentposy[sampledrightnode]) #the remaining nodes on the right
  minrightnodes<-min(rightnodes)
  
  centralnodes<-c()
  if((maxleftnodes+1)<minrightnodes){ #the nodes inbetween
    centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
  }
  rescorenodes<-c(leftnodes,centralnodes,sampledrightnode) #all these need to be rescored
  if(selectedelement>1){ #as well as the partition element further left, if it exists
    rescorenodes<-c(which(currentposy==(selectedelement-1)),rescorenodes)
  }
  
  proposedpermy<-currentpermy #create the new permutation
  proposedpermy[sampledelements]<-currentpermy[rev(sampledelements)] #by swapping the elements
  
  rescorenodes<-proposedpermy[rescorenodes] #map to the correct labels
  scorepositions<-which(proposedpermy%in%rescorenodes)
  rescorenodes<-proposedpermy[scorepositions]
  
  
  return(list(proposedpermy,rescorenodes,scorepositions))
  
}

# swap any two adjacent nodes # move 2

swapadjacentnodes.old<-function(n,currentparty,currentposy,currentpermy,permpossibs){
  
  scorelength<-length(permpossibs) 
  selectedelement<-sample.int(scorelength,1,prob=permpossibs) # sample an element
  
  leftnodes<-which(currentposy==selectedelement) #nodes in the left partition element
  rightnodes<-which(currentposy==(selectedelement+1)) #and in the right
  sampledelements<-c(propersample(leftnodes),propersample(rightnodes)) #sample a pair
  
  rescorenodes<-c(leftnodes,sampledelements[2]) #need to rescore the left partition element and just the selected node in the right partition
  
  if(selectedelement>1){ #as well as the partition element further left, if it exists
    rescorenodes<-c(which(currentposy==(selectedelement-1)),rescorenodes)
  }
  
  proposedpermy<-currentpermy #create the new permutation
  proposedpermy[sampledelements]<-currentpermy[rev(sampledelements)] #by swapping the elements
  
  rescorenodes<-proposedpermy[rescorenodes] #map to the correct labels
  scorepositions<-which(proposedpermy%in%rescorenodes)
  rescorenodes<-proposedpermy[scorepositions]
  
  return(list(proposedpermy,rescorenodes,scorepositions))
  
}

# Split or join adjacent partition elements # move 3

partitionsplitorjoin.old<-function(n,currentparty,currentposy,currentpermy,partpossibs){
  
  scorelength<-length(partpossibs) 
  sampledelement<-sample.int(scorelength,1,prob=partpossibs) # sample an element
  
  # vectors to be updated later
  proposedposy<-currentposy
  proposedpartposs<-partpossibs
  
  if(currentposy[sampledelement]==currentposy[sampledelement+1]){ #if the next node is in the same partition element we split
    involvednodes<-which(currentposy==currentposy[sampledelement]) #the nodes in the partition element to be split
    newleftelementsize<-sampledelement-min(involvednodes)+1 #the size of the new partition element on the left
    newrightelementsize<-length(involvednodes)-newleftelementsize #and the right
    proposedparty<-append(currentparty,newrightelementsize,currentposy[sampledelement]) #create the proposed partition
    proposedparty[currentposy[sampledelement]]<-newleftelementsize
    # update the list to be slightly more efficient
    proposedposy[(sampledelement+1):n]<-currentposy[(sampledelement+1):n]+1
    
    # sample which nodes go to the new left and right partition elements # note that involvednodes must have at least two elements so `sample' should work fine
    newleftnodes<-sample(involvednodes,newleftelementsize) 
    newrightnodes<-involvednodes[-which(involvednodes%in%newleftnodes)]
    proposedpermy<-currentpermy #and update the proposed permutation according to how the nodes are distributed 
    proposedpermy[(sampledelement-newleftelementsize+1):sampledelement]<-currentpermy[newleftnodes]
    proposedpermy[(sampledelement+1):(sampledelement+newrightelementsize)]<-currentpermy[newrightnodes]
    
    # update the neighbourhood of the new left element
    proposedpartposs[(sampledelement-newleftelementsize+1):sampledelement]<-choose(newleftelementsize,1:newleftelementsize)
    
    rescorenodes<-c((sampledelement-newleftelementsize+1):sampledelement) #we only need to rescore the nodes in the new left partition element
    if(currentposy[sampledelement]>1){ # and the partition element to the left, if it exists
      rescorenodes<-c(which(currentposy==(currentposy[sampledelement]-1)),rescorenodes)
    }
    
    if(currentposy[sampledelement]<length(currentparty)){ #update the neighbourhood of the right element
      proposedpartposs[(sampledelement+1):(sampledelement+newrightelementsize)]<-choose(newrightelementsize,1:newrightelementsize)
    } else { # if the partition element was last then we cut the binomial coeffiecients short
      if(newrightelementsize>1){
        proposedpartposs[(sampledelement+1):(sampledelement+newrightelementsize-1)]<-choose(newrightelementsize,1:(newrightelementsize-1))
      }
    }
    rescorenodes<-proposedpermy[rescorenodes] #finally map to the correct labels
  }
  
  if(currentposy[sampledelement]<currentposy[sampledelement+1]){ #otherwise join partition elements
    newelementsize<-currentparty[currentposy[sampledelement]]+currentparty[currentposy[sampledelement+1]]
    proposedparty<-currentparty[-currentposy[sampledelement]] #create the proposed partition
    proposedparty[currentposy[sampledelement]]<-newelementsize
    # update the list to be slightly more efficient 
    proposedposy[(sampledelement+1):n]<-currentposy[(sampledelement+1):n]-1
    proposedpermy<-currentpermy #now changing the permutation has no effect
    
    rescorenodes<-which(currentposy==currentposy[sampledelement]) #we only need to rescore the nodes from the old element on the left and the element further left
    newelementstart<-min(rescorenodes) #store the start of the new element
    if(proposedposy[sampledelement]>1){ #include the partition element to the left, if it exists
      rescorenodes<-c(which(proposedposy==(proposedposy[sampledelement]-1)),rescorenodes)
    }
    if(proposedposy[sampledelement]<length(proposedparty)){ #update the neighbourhood 
      proposedpartposs[(newelementstart):(newelementstart+newelementsize-1)]<-choose(newelementsize,1:newelementsize)
    } else { # if the new element is last we cut the binomial coeffiecients short
      proposedpartposs[(newelementstart):(newelementstart+newelementsize-2)]<-choose(newelementsize,1:(newelementsize-1)) # the new element must have a size of at least two
    }
    rescorenodes<-proposedpermy[rescorenodes] #again map to the correct labels
  }
  
  scorepositions<-which(proposedpermy%in%rescorenodes)
  rescorenodes<-proposedpermy[scorepositions]
  return(list(proposedparty,proposedposy,proposedpermy,rescorenodes,proposedpartposs,
              scorepositions))
  
}

# move a node to a different partition element (following further restrictions) # move 4 part 1

joinnode.old<-function(n,currentparty,currentposy,currentpermy,joinpossibs){
  m<-length(currentparty)
  scorelength<-length(joinpossibs) 
  nodetomove<-sample.int(scorelength,1,prob=joinpossibs) # sample an element
  
  nodeselement<-currentposy[nodetomove]
  elementtomoveto<-sample.int(joinpossibs[nodetomove],1) # sample where to go to
  
  if(elementtomoveto==nodeselement){ #if element selected is the same, replace by m which is not otherwise sampled
    elementtomoveto<-m
  } else if (elementtomoveto==(nodeselement+1)){
    if(currentparty[nodeselement]==1){ #if partition element has size one
      if(currentparty[nodeselement+1]==1){ #if next partition element is also size one
        elementtomoveto<-m-1 #replace not allowed sample by 'm-1' which also not be otherwise sampled
      }
    }
  }
  
  proposedparty<-currentparty # update the partition
  proposedparty[elementtomoveto]<-currentparty[elementtomoveto]+1
  if(currentparty[nodeselement]>1){
    proposedparty[nodeselement]<-currentparty[nodeselement]-1
  } else {
    proposedparty<-proposedparty[-nodeselement]
  }
  proposedposy<-parttolist(n,proposedparty) #should be updated for efficiency
  
  if(elementtomoveto<nodeselement){
    leftnodes<-which(currentposy==elementtomoveto)
    maxleftnodes<-max(leftnodes)
    rightnodes<-which(currentposy==nodeselement) # this includes the node which is moved
    minrightnodes<-min(rightnodes)
    proposedpermy<-currentpermy[-nodetomove]
    proposedpermy<-append(proposedpermy,currentpermy[nodetomove],maxleftnodes)
    
    centralnodes<-c()
    if((maxleftnodes+1)<minrightnodes){
      centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
    }
    rescorenodes<-c(leftnodes,centralnodes,nodetomove) #all these need to be rescored
    if(elementtomoveto>1){ #as well as the partition element further left, if it exists
      rescorenodes<-c(which(currentposy==(elementtomoveto-1)),rescorenodes)
    }
  } else {
    leftnodes<-which(currentposy==nodeselement) # this includes the node which is moved
    maxleftnodes<-max(leftnodes)
    rightnodes<-which(currentposy==elementtomoveto)
    minrightnodes<-min(rightnodes)
    proposedpermy<-currentpermy[-nodetomove]
    proposedpermy<-append(proposedpermy,currentpermy[nodetomove],min(rightnodes)-2)
    
    centralnodes<-c()
    if((maxleftnodes+1)<minrightnodes){
      centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
    }
    rescorenodes<-c(leftnodes,centralnodes) #all these need to be rescored - the moved node is included in leftnodes
    if(nodeselement>1){ #as well as the partition element further left, if it exists
      rescorenodes<-c(which(currentposy==(nodeselement-1)),rescorenodes)
    }
  }
  
  rescorenodes<-currentpermy[rescorenodes] #map to the correct labels
  scorepositions<-which(proposedpermy%in%rescorenodes)
  rescorenodes<-proposedpermy[scorepositions]
  
  return(list(proposedparty,proposedposy,proposedpermy,rescorenodes,scorepositions))
  
}

# move a node to a new partition element (following further restrictions) # move 4 part 2

holenode.old<-function(n,currentparty,currentposy,currentpermy,holepossibs){
  m<-length(currentparty)
  
  scorelength<-length(holepossibs) 
  nodetomove<-sample.int(scorelength,1,prob=holepossibs) # sample an element
  
  nodeselement<-currentposy[nodetomove]
  holetomoveto<-sample.int(holepossibs[nodetomove],1) # sample where to go to
  
  if(currentparty[nodeselement]==1){
    if(holetomoveto>=nodeselement){
      holetomoveto<-holetomoveto+1 #counts hole when partition element removed
      if(nodeselement<m){
        if(currentparty[nodeselement+1]==1){ #if next partition element is also size one
          holetomoveto<-holetomoveto+1 #shift one hole further along
        }
      }
    }
    proposedparty<-currentparty[-nodeselement]
    proposedparty<-append(proposedparty,1,holetomoveto-1)
    proposedposy<-parttolist(n,proposedparty) #should be updated for efficiency
    
    if(holetomoveto<nodeselement){
      leftnodes<-which(proposedposy==holetomoveto) # this should be the node which is moved
      maxleftnodes<-max(leftnodes)
      rightnodes<-which(proposedposy==nodeselement) # these were to the left of the moved node
      minrightnodes<-min(rightnodes)
      proposedpermy<-currentpermy[-nodetomove]
      proposedpermy<-append(proposedpermy,currentpermy[nodetomove],leftnodes-1)
      
      centralnodes<-c()
      if((maxleftnodes+1)<minrightnodes){
        centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      }
      rescorenodes<-c(leftnodes,centralnodes,rightnodes) #all these need to be rescored
      if(holetomoveto>1){ #as well as the partition element further left, if it exists
        rescorenodes<-c(which(currentposy==(holetomoveto-1)),rescorenodes)
      }
    } else {
      leftnodes<-which(proposedposy==nodeselement) # these were to the right of the moved node
      maxleftnodes<-max(leftnodes)
      rightnodes<-which(proposedposy==holetomoveto) # this should be the node that is moved
      minrightnodes<-min(rightnodes)
      proposedpermy<-currentpermy[-nodetomove]
      proposedpermy<-append(proposedpermy,currentpermy[nodetomove],rightnodes-1)
      
      centralnodes<-c()
      if((maxleftnodes+1)<minrightnodes){
        centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      }
      rescorenodes<-c(leftnodes,centralnodes,rightnodes) #all these need to be rescored
      if(nodeselement>1){ #as well as the partition element further left, if it exists
        rescorenodes<-c(which(currentposy==(nodeselement-1)),rescorenodes)
      }
    }
  }  else {
    if(currentparty[nodeselement]==2){
      if(holetomoveto==nodeselement){
        holetomoveto<-m+1
      }
    }
    proposedparty<-currentparty
    proposedparty[nodeselement]<-currentparty[nodeselement]-1
    proposedparty<-append(proposedparty,1,holetomoveto-1)
    proposedposy<-parttolist(n,proposedparty) #should be updated for efficiency
    
    if(holetomoveto<=nodeselement){
      leftnodes<-which(proposedposy==holetomoveto) # this is the node which is moved
      maxleftnodes<-max(leftnodes)
      rightnodes<-which(proposedposy==(nodeselement+1)) # these are the nodes left over
      minrightnodes<-min(rightnodes)
      proposedpermy<-currentpermy[-nodetomove]
      proposedpermy<-append(proposedpermy,currentpermy[nodetomove],leftnodes-1)
      
      centralnodes<-c()
      if((maxleftnodes+1)<minrightnodes){
        centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      }
      rescorenodes<-c(leftnodes,centralnodes) #all these need to be rescored
      if(holetomoveto>1){ #as well as the partition element further left, if it exists
        rescorenodes<-c(which(currentposy==(holetomoveto-1)),rescorenodes)
      }
    } else {
      leftnodes<-which(proposedposy==nodeselement) # these are the nodes left over
      maxleftnodes<-max(leftnodes)
      rightnodes<-which(proposedposy==holetomoveto) # this is the node which is moved
      minrightnodes<-min(rightnodes)
      proposedpermy<-currentpermy[-nodetomove]
      proposedpermy<-append(proposedpermy,currentpermy[nodetomove],rightnodes-1)
      
      centralnodes<-c()
      if((maxleftnodes+1)<minrightnodes){
        centralnodes<-c((maxleftnodes+1):(minrightnodes-1))
      }
      rescorenodes<-c(leftnodes,centralnodes,rightnodes) #all these need to be rescored
      if(nodeselement>1){ #as well as the partition element further left, if it exists
        rescorenodes<-c(which(currentposy==(nodeselement-1)),rescorenodes)
      }
    }
    
  }
  
  rescorenodes<-proposedpermy[rescorenodes] #map to the correct labels
  scorepositions<-which(proposedpermy%in%rescorenodes)
  rescorenodes<-proposedpermy[scorepositions]
  return(list(proposedparty,proposedposy,proposedpermy,rescorenodes,
              scorepositions))
  
}



