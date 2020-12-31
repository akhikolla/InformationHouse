#' Score IBD length.
#' 
#' \code{ibd.length} returns the total length of IBD segemnt between two haplotypes.
#' 
#' This function works with output from \code{\link{sim.recomb}}. 
#' 
#' @param inher.hap1,inher.hap2 numeric matrix.
#' @param startpos,endpos non-negative number.
#' @return A non-negative number representing the length of IBD segment in Haldane centiMorgan.
#' @examples 
#' # a simple pedigree with sibling marriage
#' pedigree = as.character(rep(1, 5))
#' member = as.character(c(11, 12, 21, 22, 31))
#' sex = as.numeric(c(1, 2, 1, 2, 1))
#' father = as.character(c(NA, NA, 11, 11, 21))
#' mother = as.character(c(NA, NA, 12, 12, 22))
#' pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
#' inheritance = sim.recomb(pedinfo, 100)
#' 
#' # IBD length between the two haplotypes of inbred individual 31
#' ibd.length(inheritance[[9]], inheritance[[10]])
#' @export
ibd.length = function(inher.hap1, inher.hap2, startpos = NULL, endpos = NULL){
  n1 = nrow(inher.hap1)
  n2 = nrow(inher.hap2)

  if(is.null(startpos)){
    startpos = 0
  }

  if(is.null(endpos)){
    endpos = inher.hap1[n1,2]
  }
  
  if(startpos < 0 || endpos > inher.hap1[n1,2] || startpos > endpos){
    stop("Check starting and ending genetic positions.")
  }
  
  index1 = 1
  index2 = 1

  relatedness = 0
  ibd = 0

  while(index1 <= n1 && index2 <= n2){
    # start from specified starting position
    if(inher.hap1[index1,2] < startpos){
      index1 = index1 + 1
      next
    }else if(inher.hap2[index2,2] < startpos){
      index2 = index2 + 1
      next
    }
    
    # get current ibd status
    if(inher.hap1[index1,1] == inher.hap2[index2,1]){
      ibd = 1
    }else{
      ibd = 0
    }

    if(inher.hap1[index1,2] >= endpos && inher.hap2[index2,2] >= endpos){
      relatedness = relatedness + ibd * (endpos - startpos)
      break
    }else{
      if(inher.hap1[index1,2] < inher.hap2[index2,2]){
        relatedness = relatedness + ibd * (inher.hap1[index1,2] - startpos)
        startpos = inher.hap1[index1,2]
        index1 = index1 + 1
      }else{
        relatedness = relatedness + ibd * (inher.hap2[index2,2] - startpos)
        startpos = inher.hap2[index2,2]
        index2 = index2 + 1
      }
    }
  }
  return(as.numeric(relatedness))
}


#' Score IBD proportion.
#' 
#' \code{ibd.proportion} returns the proportion of IBD sharing between two haplotypes of the same individual or two individuals.
#' 
#' When only one individual index is supplied, \code{ibd.proportion} returns the realized inbreeding coefficient of the individual. When two individual indices are supplied, \code{ibd.proportion} returns the realized relatedness of the two individuals.
#' 
#' @param inheritance list of matrices.
#' @param ind1index,ind2index positive integer.
#' @param startpos,endpos non-negative number.
#' @return A value between 0 and 1 representing the proportion of IBD segment.
#' @examples 
#' # a simple pedigree with sibling marriage
#' pedigree = as.character(rep(1, 5))
#' member = as.character(c(11, 12, 21, 22, 31))
#' sex = as.numeric(c(1, 2, 1, 2, 1))
#' father = as.character(c(NA, NA, 11, 11, 21))
#' mother = as.character(c(NA, NA, 12, 12, 22))
#' pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
#' inheritance = sim.recomb(pedinfo, 100)
#' 
#' # realized inbreeding of inbred child
#' get.pedindex(pedinfo, "31")
#' ibd.proportion(inheritance, 5)
#' 
#' # realized relatedness between individual 21 and 22 (parents of inbred child)
#' get.pedindex(pedinfo, c("21", "22"))
#' ibd.proportion(inheritance, 3, 4)
#' @export
ibd.proportion = function(inheritance, ind1index, ind2index = NULL, startpos = NULL, endpos = NULL){
  if(is.null(startpos)){
    startpos = 0
  }
  
  if(is.null(endpos)){
    endpos = tail(c(inheritance[[1]]), 1)
  }
  seglength = endpos - startpos
  
  if(is.null(ind2index)){
    proportion = ibd.length(inheritance[[2*ind1index-1]], inheritance[[2*ind1index]], startpos, endpos)/seglength
  }else{
    proportion = (ibd.length(inheritance[[2*ind1index-1]], inheritance[[2*ind2index-1]], startpos, endpos) + ibd.length(inheritance[[2*ind1index-1]], inheritance[[2*ind2index]], startpos, endpos) + ibd.length(inheritance[[2*ind1index]], inheritance[[2*ind2index-1]], startpos, endpos) + ibd.length(inheritance[[2*ind1index]], inheritance[[2*ind2index]], startpos, endpos)) / 2 / seglength
  }
  return(proportion)
}


#' Score IBD state.
#' 
#' \code{fgl2ibd} determines pairwise IBD state given the four founder genome labels of two individuals at a marker.
#' 
#' IBD states take value from 1 to 15, which represent the indices of the underlying IBD states from 1111 to 1234 in lexicographical order. E.g., output 1 means IBD state 1111, output 2 means IBD state 1112 etc. Recoding in, e.g., Jacquard order, can be obtained using \code{\link{recode.ibd}}.
#' 
#' @param fgl1p,fgl1m,fgl2p,fgl2m positive integer, represents founder genome label.
#' @return A value between 1 and 15 representing index of IBD state in lexicographical order.
#' @export
#' @examples 
#' fgl2ibd(1, 1, 1, 1)
#' fgl2ibd(1, 2, 1, 2)
#' fgl2ibd(3, 4, 5, 6)
#' fgl2ibd(4, 5, 4, 4)
fgl2ibd = function(fgl1p, fgl1m, fgl2p, fgl2m){
  if(fgl1p == fgl1m){
    if(fgl2p == fgl2m){
      if(fgl1p ==fgl2p){
        ibdstate = 1 # 1111
      }else{
        ibdstate = 4 # 1122
      }
    }else{
      if(fgl1p == fgl2p){
        ibdstate = 2 # 1112
      }else if(fgl1p == fgl2m){
        ibdstate = 3 # 1121
      }else{
        ibdstate = 5 # 1123
      }
    }
  }else{
    if(fgl2p == fgl2m){
      if(fgl1p == fgl2p){
        ibdstate = 6 # 1211
      }else if(fgl1m == fgl2p){
        ibdstate = 10 # 1222
      }else{
        ibdstate = 14 # 1233
      }
    }else{
      if(fgl1p == fgl2p){
        if(fgl1m == fgl2m){
          ibdstate = 7 # 1212
        }else{
          ibdstate = 8 # 1213
        }
      }else if(fgl1p == fgl2m){
        if(fgl1m == fgl2p){
          ibdstate = 9 # 1221
        }else{
          ibdstate = 12 # 1231
        }
      }else{
        if(fgl1m == fgl2p){
          ibdstate = 11 # 1223
        }else if(fgl1m == fgl2m){
          ibdstate = 13 # 1232
        }else{
          ibdstate = 15 # 1234
        }
      }
    }
  }
  return(ibdstate)
}

#' Score pairwise relatedness.
#' 
#' \code{fgl2relatedness} determines pairwise relatedness given the four founder genome labels of two individuals at a marker.
#' 
#' @param fgl1p,fgl1m,fgl2p,fgl2m positive integer, represents founder genome label.
#' @return A value in [0, 0.5, 1, 2] representing local relatedness coefficient.
#' @export
#' @examples 
#' fgl2relatedness(1, 1, 1, 1)
#' fgl2relatedness(1, 2, 1, 2)
#' fgl2relatedness(1, 2, 1, 3)
#' fgl2relatedness(3, 4, 5, 6)
#' fgl2relatedness(4, 5, 4, 4)
fgl2relatedness = function(fgl1p, fgl1m, fgl2p, fgl2m){
  relatedness = ((fgl1p == fgl2p) + (fgl1p == fgl2m) + (fgl1m == fgl2p) + (fgl1m == fgl2m))/2
  return(relatedness)
}


#' Score IBD sharing by segment.
#' 
#' \code{ibd.segment} determines the starting and endping genetic positions of segments with different amount of pairwise IBD sharing.
#' 
#' When only index of one individual is supplied, IBD sharing status for each segment is coded as 0 (not IBD) or 1 (IBD) between the two haplotypes of the individual.
#' 
#' When indices of two individuals are supplied, IBD sharing status for each segment is either in relatedness (default) or lexicographical order of IBD state, where recoding can be done using \code{\link{recode.ibd}}.
#' 
#' @param inheritance list of numeric matrices.
#' @param ind1index,ind2index positive integer, represents index of individual in pedigree.
#' @param relatedness logical, determines coding of IBD information.
#' @return A dataframe of three variables. \code{ibd} represents IBD sharing status of a segment, and \code{startpos}/\code{endpos} represents starting/ending genetic position of the segment.
#' @export
#' @examples 
#' # a simple pedigree with sibling marriage
#' pedigree = as.character(rep(1, 5))
#' member = as.character(c(11, 12, 21, 22, 31))
#' sex = as.numeric(c(1, 2, 1, 2, 1))
#' father = as.character(c(NA, NA, 11, 11, 21))
#' mother = as.character(c(NA, NA, 12, 12, 22))
#' pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
#' inheritance = sim.recomb(pedinfo, 100)
#' 
#' # IBD segments between the two haplotypes of the inbred individual
#' ibd.segment(inheritance, 5)
#' 
#' # IBD segments between the two full sibs
#' ibd.segment(inheritance, 3, 4) # relatedness
#' ibd.segment(inheritance, 3, 4, relatedness = FALSE) # lexicographical order of IBD state
ibd.segment = function(inheritance, ind1index, ind2index = NULL, relatedness = TRUE){
  ibd = NULL
  startpos = 0
  endpos = NULL
  
  if(is.null(ind2index)){ # between two haplotypes
    N = nrow(inheritance[[2*ind1index-1]]) + nrow(inheritance[[2*ind1index]])
    ibd = ifelse(inheritance[[2*ind1index-1]][1, 1] == inheritance[[2*ind1index]][1, 1], 1, 0)
    recomb.index = c(1, 1)
    
    data = list(inheritance[[2*ind1index-1]], inheritance[[2*ind1index]])
    
    while(sum(recomb.index) < N){
      min.index = which.min(c(data[[1]][recomb.index[1], 2], data[[2]][recomb.index[2], 2]))
      min.recomb = data[[min.index]][recomb.index[min.index], 2]
      if(data[[1]][recomb.index[1],2] == min.recomb){
        recomb.index[1] =  recomb.index[1] + 1
      }
      if(data[[2]][recomb.index[2],2] == min.recomb){
        recomb.index[2] =  recomb.index[2] + 1
      }
      ibd.temp = ifelse((data[[1]][recomb.index[1], 1] == data[[2]][recomb.index[2], 1]), 1, 0)
      if(ibd.temp != tail(ibd, 1)){
        ibd = c(ibd, ibd.temp)
        endpos = c(endpos, data[[min.index]][recomb.index[min.index]-1, 2])
        startpos = c(startpos, data[[min.index]][recomb.index[min.index]-1, 2])
      }
    }
    endpos = c(endpos, inheritance[[2*ind1index-1]][recomb.index[1], 2])
    return(data.frame(ibd, startpos, endpos))
  }else{ # between two individuals
    N = nrow(inheritance[[2*ind1index-1]]) + nrow(inheritance[[2*ind1index]]) + nrow(inheritance[[2*ind2index-1]]) + nrow(inheritance[[2*ind2index]])
    recomb.index = rep(1, 4)
    
    data = list(inheritance[[2*ind1index-1]], inheritance[[2*ind1index]], inheritance[[2*ind2index-1]], inheritance[[2*ind2index]])
    
    if(relatedness){ # output local relatedness
      relatedness = fgl2relatedness(data[[1]][1, 1], data[[2]][1, 1], data[[3]][1, 1], data[[4]][1, 1])
      while(sum(recomb.index) < N){
        min.index = which.min(c(data[[1]][recomb.index[1], 2], data[[2]][recomb.index[2], 2], data[[3]][recomb.index[3], 2], data[[4]][recomb.index[4], 2]))
        min.recomb = data[[min.index]][recomb.index[min.index], 2]
        recomb.index[min.index] = recomb.index[min.index] + 1
        while(min(c(data[[1]][recomb.index[1], 2], data[[2]][recomb.index[2], 2], data[[3]][recomb.index[3], 2], data[[4]][recomb.index[4], 2])) == min.recomb){
          min.index = which.min(c(data[[1]][recomb.index[1], 2], data[[2]][recomb.index[2], 2], data[[3]][recomb.index[3], 2], data[[4]][recomb.index[4], 2]))
          recomb.index[min.index] = recomb.index[min.index] + 1
        }
        relatedness.temp = fgl2relatedness(inheritance[[2*ind1index-1]][recomb.index[1], 1], inheritance[[2*ind1index]][recomb.index[2], 1], inheritance[[2*ind2index-1]][recomb.index[3], 1], inheritance[[2*ind2index]][recomb.index[4], 1])
        if(relatedness.temp != tail(relatedness, 1)){
          relatedness = c(relatedness, relatedness.temp)
          endpos = c(endpos, data[[min.index]][recomb.index[min.index]-1, 2])
          startpos = c(startpos, data[[min.index]][recomb.index[min.index]-1, 2])
        }
      }
      endpos = c(endpos, inheritance[[2*ind1index-1]][recomb.index[1], 2])
      return(data.frame(relatedness, startpos, endpos))
    }else{ # output ibd state
      ibd = fgl2ibd(data[[1]][1, 1], data[[2]][1, 1], data[[3]][1, 1], data[[4]][1, 1])
      while(sum(recomb.index) < N){
        min.index = which.min(c(data[[1]][recomb.index[1], 2], data[[2]][recomb.index[2], 2], data[[3]][recomb.index[3], 2], data[[4]][recomb.index[4], 2]))
        min.recomb = data[[min.index]][recomb.index[min.index], 2]
        recomb.index[min.index] = recomb.index[min.index] + 1
        while(min(c(data[[1]][recomb.index[1], 2], data[[2]][recomb.index[2], 2], data[[3]][recomb.index[3], 2], data[[4]][recomb.index[4], 2])) == min.recomb){
          min.index = which.min(c(data[[1]][recomb.index[1], 2], data[[2]][recomb.index[2], 2], data[[3]][recomb.index[3], 2], data[[4]][recomb.index[4], 2]))
          recomb.index[min.index] = recomb.index[min.index] + 1
        }
        ibd.temp = fgl2ibd(inheritance[[2*ind1index-1]][recomb.index[1], 1], inheritance[[2*ind1index]][recomb.index[2], 1], inheritance[[2*ind2index-1]][recomb.index[3], 1], inheritance[[2*ind2index]][recomb.index[4], 1])
        if(ibd.temp != tail(ibd, 1)){
          ibd = c(ibd, ibd.temp)
          endpos = c(endpos, data[[min.index]][recomb.index[min.index]-1, 2])
          startpos = c(startpos, data[[min.index]][recomb.index[min.index]-1, 2])
        }
      }
      endpos = c(endpos, inheritance[[2*ind1index-1]][recomb.index[1], 2])
      return(data.frame(ibd, startpos, endpos))
    }
  }
}


#' Score IBD sharing at a list of marker positions.
#' 
#' \code{ibd.marker} determines pairwise IBD sharing at marker positions.
#' 
#' When only index of one individual is supplied, IBD sharing status at each marker is coded as 0 (not IBD) or 1 (IBD) between the two haplotypes of the individual.
#' 
#' When indices of two individuals are supplied, IBD sharing status at each marker is either in relatedness (default) or lexicographical order of IBD state, where recoding can be done using \code{\link{recode.ibd}}.
#' 
#' @param inheritance list of numeric matrices.
#' @param marker numeric vector.
#' @param ind1index,ind2index positive integer, represents index of individual in pedigree.
#' @param relatedness logical, determines coding of IBD information.
#' @return A numeric vector of IBD sharing status at the list of marker positions.
#' @export
#' @examples 
#' # a simple pedigree with sibling marriage
#' pedigree = as.character(rep(1, 5))
#' member = as.character(c(11, 12, 21, 22, 31))
#' sex = as.numeric(c(1, 2, 1, 2, 1))
#' father = as.character(c(NA, NA, 11, 11, 21))
#' mother = as.character(c(NA, NA, 12, 12, 22))
#' pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
#' inheritance = sim.recomb(pedinfo, 100)
#' nsnp = 10
#' marker = sort(runif(nsnp, 0, 100))
#' 
#' # IBD at markers between the two haplotypes of the inbred individual
#' ibd.marker(inheritance, marker, 5)
#' 
#' # IBD at markers between the two full sibs, with different IBD coding
#' ibd.marker(inheritance, marker, 3, 4) # relatedness
#' ibd.marker(inheritance, marker, 3, 4, relatedness = FALSE) # lexicographical order of IBD state
ibd.marker = function(inheritance, marker, ind1index, ind2index = NULL, relatedness = TRUE){
  L = tail(c(inheritance[[1]]),1)
  if(min(marker) < 0 || max(marker) > L){
    stop("Marker positon out of bounds.")
  }
  
  nsnp = length(marker)
  output = rep(0, nsnp)
  
  if(relatedness){
    segment = ibd.segment(inheritance, ind1index, ind2index)
  }else{
    segment = ibd.segment(inheritance, ind1index, ind2index, relatedness = FALSE)
  }
  
  current.index = 1
  current.ibd = segment[1, 1]
  current.recomb = segment[1, 3]
  start.mindex = 1
  end.mindex = 1
  L = segment[nrow(segment), 3]
  
  if(current.recomb == L){
    output = rep(current.ibd, nsnp)
  }
  
  while(current.recomb < L && end.mindex <= nsnp){
    if(marker[end.mindex] > current.recomb){
      if(start.mindex != end.mindex){
        output[start.mindex:(end.mindex-1)] = rep(current.ibd, (end.mindex-start.mindex))
        start.mindex = end.mindex
      }
      current.index = current.index + 1
      current.ibd = segment[current.index, 1]
      current.recomb = segment[current.index, 3]
    }else{
      end.mindex = end.mindex + 1
    }
  }
  output[start.mindex:nsnp] = rep(current.ibd, (nsnp-start.mindex+1))
  return(output)
}
