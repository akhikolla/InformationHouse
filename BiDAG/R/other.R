# $Id: combinations.R 1083 2007-03-23 22:53:00Z warnes $
#

##
## From email by Brian D Ripley <ripley@stats.ox.ac.uk> to r-help
## dated Tue, 14 Dec 1999 11:14:04 +0000 (GMT) in response to
## Alex Ahgarin <datamanagement@email.com>.  Original version was
## named "subsets" and was Written by Bill Venables.  
##

combinations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed=FALSE) {
  if(mode(n) != "numeric" || length(n) != 1 
     || n < 1 || (n %% 1) != 0) stop("bad value of n") 
  if(mode(r) != "numeric" || length(r) != 1 
     || r < 1 || (r %% 1) != 0) stop("bad value of r") 
  if(!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    v <- unique(sort(v))
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
    { 
      if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
          if(n == 1) matrix(v, 1, r) else
            rbind( cbind(v[1], Recall(n, r-1, v)),
                   Recall(n-1, r, v[-1]))
    }
  else
    sub <- function(n, r, v)
    { 
      if(r == 0) v0 else
        if(r == 1) matrix(v, n, 1) else
          if(r == n) matrix(v, 1, n) else
            rbind(cbind(v[1], Recall(n-1, r-1, v[-1])),
                  Recall(n-1, r, v[-1]))
    }
  sub(n, r, v[1:n])
}

##
## Original version by Bill Venables and cited by by Matthew
## Wiener (mcw@ln.nimh.nih.gov) in an email to R-help dated
## Tue, 14 Dec 1999 09:11:32 -0500 (EST) in response to
## Alex Ahgarin <datamanagement@email.com>
##
##


permutations <- function(n, r, v = 1:n, set = TRUE, repeats.allowed=FALSE)
{
  if(mode(n) != "numeric" || length(n) != 1 
     || n < 1 || (n %% 1) != 0) stop("bad value of n") 
  if(mode(r) != "numeric" || length(r) != 1 
     || r < 1 || (r %% 1) != 0) stop("bad value of r") 
  if(!is.atomic(v) || length(v) < n) 
    stop("v is either non-atomic or too short")
  if( (r > n) & repeats.allowed==FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if(set) {
    v <- unique(sort(v))
    if (length(v) < n) stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  ## Inner workhorse
  if(repeats.allowed)
    sub <- function(n, r, v)
    {
      if(r==1) matrix(v,n,1) else
        if(n==1) matrix(v,1,r) else
        {
          inner  <-  Recall(n, r-1, v)
          cbind( rep( v, rep(nrow(inner),n)  ),
                 matrix( t(inner), ncol=ncol(inner), nrow=nrow(inner) * n ,
                         byrow=TRUE )
          )
        }
    }
  else
    sub <- function(n, r, v)
    {
      if(r==1) matrix(v,n,1) else
        if(n==1) matrix(v,1,r) else
        {
          X  <-  NULL
          for(i in 1:n)
            X  <-  rbind( X, cbind( v[i], Recall(n-1, r - 1, v[-i])))
          X
        }
    }
  
  sub(n, r, v[1:n])
}


# This function samples an element from a vector properly
propersample <- function(x){if(length(x)==1) x else sample(x,1)} 


compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

is.subset<-function(v1,v2){
  issubset<-TRUE
  l<-length(v1)
  for  (i in 1:l){
    if (is.na(v1[i])) {break }
    else  {issubset<-v1[i] %in% v2}
  }
  return(issubset)
}

orderbgn<-function(permy,bgnodes) {
  if(!is.null(bgnodes)) {
    movenodes<-which(permy%in%bgnodes)
    newpermy<-permy[-movenodes]
    return(c(newpermy,bgnodes))
  } else {
    return(permy)
  }
}

storemaxMCMC<-function(MCMCres,param) {
  maxobj<-list()
  maxN<-which.max(unlist(MCMCres[[2]]))
  if(param$DBN) {
    maxobj$DAG<-DBNtransform(MCMCres[[1]][[maxN]],param)
    maxobj$DAGorig<-MCMCres[[1]][[maxN]]
    maxobj$order<-order2var(MCMCres[[4]][[maxN]],param$firstslice$labels)
  } else {
    maxobj$DAG<-MCMCres[[1]][[maxN]]
    colnames(maxobj$DAG)<-param$labels
    rownames(maxobj$DAG)<-param$labels
    maxobj$order<-order2var(MCMCres[[4]][[maxN]],param$labels)
  }
  maxobj$score<-MCMCres[[2]][[maxN]]
  attr(maxobj,"class")<-"MCMCmax"
  return(maxobj)
}
assignLabels<-function(adjacency,nodelabels){
  colnames(adjacency)<-nodelabels
  rownames(adjacency)<-nodelabels
  return(adjacency)
}

defcolrange<-function(value) {
  if(value>0.8) {
    return(5)
  } else if (value>0.6) {
    return(4)
  } else if (value>0.4) {
    return(3)
  } else if (value>0.2) {
    return(2)
  } else {
    return(1)
  }
}


#checking startorder, if NULL generating random order of right length
checkstartorder<-function(order,varnames,mainnodes,bgnodes,
                          DBN=FALSE, split=FALSE) {
  
  matsize<-length(varnames) #maximum length of the startorder
  nsmall<-length(mainnodes)#minimum length of the startorder
  errortext<-"ok"
  errorflag<-0
  bgn<-length(bgnodes)
  n<-nsmall+bgn
  lo<-length(order)
  
  error1<-"startorder should contain either variable names or their respective indices in the data object!"
  
  error2<-"the variables (or indices) in the startorder should be similar to variable names (or their indices in the data object)!"
  
  error3<-"DBN samestruct=FALSE should have following format: vector of length equal to
  at 2*nsmall or 2*nsmall+bgn, where nsmall=number of dynamic variables and bgn=number of static variables.
  The first half of the order should contain indices (bgn+1):(nsmall+bgn) representing the order of 
  variables in the initial time slice, the second half represents the order of variables
  in any other time slice should contain indices (nsmall+bgn+1):(2*nsmall+bgn)"
  

  if(!DBN) { #usual BN
    if(is.null(order)) {
      order<-c(sample(mainnodes,nsmall,replace=FALSE),bgnodes)
    } else {
      
      if(all(is.character(order))) {#convert to indices
        order<-match(order,varnames)
      }
      if(any(is.na(order))) {
        errortext<-error2
        errorflag<-2
      } else if(!all(is.numeric(order)))  {
        errortext<-error1
        errorflag<-1
      } else if (!setequal(order,mainnodes) & !setequal(order,c(1:n))) {
        errortext<-error2
        errorflag<-2
      } else {#startorder is defined correctly
        if(lo==n) {
          order<-orderbgn(order,bgnodes)#need to put bgnodes to the end
        }
        else {
          order<-c(order,bgnodes)#or attach them if they are missing
        }
      }
    }
  } else {#DBN
    if(split) {#DBN: initial and transition structure (internal edges) are different
      
      if(is.null(order)) {#need to define initial and transition order
        
        order<-list()
        
        #we change indices of variables turning into internal representation for DBNs
        order$init<-c(sample(mainnodes-bgn,nsmall,replace=FALSE),bgnodes+nsmall)
        order$trans<-c(sample(mainnodes-bgn,nsmall,replace=FALSE),1:n+nsmall)
        
      } else {
        
        if(all(is.character(order))) {#convert to indices
          order<-match(order,varnames)
        }
        
        mainall<-c(mainnodes,mainnodes+nsmall)
        
        if(any(is.na(order))) {
            errortext<-error2
            errorflag<-2
          } else if(!all(is.numeric(order)))  {
            errortext<-error1
            errorflag<-1
          } else if (!setequal(order,c(1:matsize)) & !setequal(order,mainall)) {
            errortext<-error2
            errorflag<-2
           } else { 
            
            order.init<-order[1:nsmall]-bgn #get order for initial structure
            order.trans<-order[1:nsmall+nsmall]-nsmall-bgn #get order for transition structure
            if(!setequal(order.init,c(1:nsmall)) | !setequal(order.trans,c(1:nsmall))) {
              errortext<-error3
              errorflag<-3
             } else {#startorder is defined correctly
                order<-list()
                order$init<-c(order.init,bgnodes+nsmall)
                order$trans<-c(order.trans,1:n+nsmall)
            }
        }
      }
    }  else {#DBN: initial and transition structure (internal edges) are the same
      if(is.null(order)) {
        
        #we change indices of variables turning into internal representation
        order<-c(sample(mainnodes,nsmall,replace=FALSE)-bgn,bgnodes+nsmall)
        
      } else {
        
        if(all(is.character(order)))  {
          order<-match(order,varnames)
        }
        
        if (!setequal(order,mainnodes) & !setequal(order,c(1:n))) {
          errortext<-error2
          errorflag<-2
        } else {#startorder is defined correctly
            if(lo==n) {
              order<-orderbgn(order,bgnodes)
              print(order)
              #we change indices of variables turning into internal representation
              if(bgn>0) {
                order<-c(order[1:nsmall]-bgn,order[1:bgn+nsmall]+nsmall)
              }
            } else {
              order<-c(order-bgn,bgnodes+nsmall)
            }
        }
      }
    }
    
  }
  
  res<-list()
  res$errorflag<-errorflag
  res$order<-order
  res$errortext<-errortext
  
  return(res)
}


#transform numeric order into varnames
order2var<-function(order,varnames) {
  return(varnames[order])
}






