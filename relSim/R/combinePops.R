combinePops = function(pop1, pop2, frac = 1){
  
  ## check if the populations have the same locus sets
  
  if(pop1$nLoci == pop2$nLoci){
    if(all(pop1$Freqs$loci %in% pop2$Freqs$loci)){ ## note this doesn't deal with loci in different order.
      
      ## firstly turn each population into a matrix
      
      p1 = matrix(pop1$profiles, nrow = pop1$nProfiles, ncol = 2 * pop1$nLoci)
      p2 = matrix(pop2$profiles, nrow = pop2$nProfiles, ncol = 2 * pop2$nLoci)
      
      ## determine the locus order and permute the columns accordingly
      m = match(pop2$Freqs$loci, pop2$Freqs$loci)
      p2 = p2[,m]
      
      if(frac != 1){
        if(frac <= 0 || frac > 1){
          stop("The fraction sampled from population 2 must be between 0 and 1")
        }
        
        ## From Thomas Lumley
        sainte_lague = function(votes, nseats) {
          nparties = length(votes)
          denominators = 2 * (1:nseats) - 1
          quotients = outer(votes, denominators, "/")
          last = sort(quotients, decreasing = TRUE)[nseats]
          clear = rowSums(quotients > last)
          borderline = rowSums(quotients == last)
          borderline[sample(which(borderline > 0), sum(borderline) - (nseats -
                                                                        sum(clear)))] = 0
          total = clear + borderline
          error = votes - sum(votes) * total / nseats
          rval = rep(1:nparties, clear + borderline)
          
          return(list(sizes = table(rval), rval = rval, error = error))
        }
        
        
        id = matrix(1:nrow(p2), nrow = pop2$ns, byrow = TRUE)
      }
      
      p = rbind(p1, p2)
      
      pop = list(profiles = as.vector(t(p)), nProfiles = nrow(p),  nSubpops = pop1$ns + pop2$ns,
                 nLoci = pop1$nLoci, theta = 0, Freqs = p1$Freqs)
      class(pop) = "population"
      
      return(pop)
      
      
    }else{
      stop("Each population must have the same locus set")
    }
  }else{
    stop("Each population must have an equal number of loci")
  }
  
  
  
}