predict.fast = function(object, data= NULL, quantiles= c(0.1,0.5,0.9),obs=1,...) {
  
  origObs = object$origObs
  nobs = length(origObs)
  origNodes = object$origNodes
  ntree = object$num.trees
  thres = 5*.Machine$double.eps
  filterednodes = rep(0, nobs*ntree)
  z = matrix(nrow=nobs, ncol=ntree)
  newnodes = matrix(nrow = nobs, ncol = ntree)
  newindex = matrix(0, nrow = nobs, ncol = ntree)
  z = apply(origNodes, 2, function(x) order(x, stats::rnorm(length(x)))) #ordering the nodes with randomization
  newnodes = sapply(seq(ncol(z)), function(x) origNodes[z[, x], x])
  
  if(is.null(data)){
    weightvec = rep(0, nobs*nobs)
    quant = matrix(nrow=nobs,ncol=length(quantiles))
    result = Findweightsinbagfast(as.double(as.vector(origNodes)),
                                  as.double(as.vector(newnodes)),
                                  as.double(filterednodes),
                                  as.integer(as.vector(z)),
                                  as.integer(as.vector(newindex)),
                                  as.integer(as.vector(unlist(t(as.data.frame(object$inbag))))),
                                  as.double(weightvec),
                                  as.integer(nobs),
                                  as.integer(ntree),
                                  as.double(thres),
                                  as.integer(obs))
    
    } else {
    nnew = nrow(data)
    weightvec = rep(0, nobs*nnew)
    quant = matrix(nrow = nrow(data), ncol = length(quantiles))
    nodes = getnodes(object, data)
    result = Findweightsfast(as.double(as.vector(newnodes)),
                             as.double(as.vector(nodes)),
                             as.double(filterednodes),
                             as.integer(as.vector(z)),
                             as.integer(as.vector(newindex)),
                             as.double(weightvec),
                             as.integer(nobs),
                             as.integer(nnew),
                             as.integer(ntree),
                             as.double(thres),
                             as.integer(obs))
  }
  
  weights = matrix(result, nrow = nobs)
  
  ord = order(origObs)
  origObs = origObs[ord]
  weights = weights[ord, , drop = FALSE]
  cumweights = apply(weights, 2, cumsum)
  cumweights = sweep(cumweights, 2, as.numeric(cumweights[nobs,]), FUN = "/")
  
  for (qc in 1:length(quantiles)){
    larg = cumweights<quantiles[qc]
    wc = apply(larg, 2, sum)+1
    ind1 = which(wc<1.1) 
    indn1 = which(wc>1.1)
    quant[ind1,qc] = rep(origObs[1], length(ind1))
    quantmax = origObs[wc[indn1]]
    quantmin = origObs[wc[indn1]-1]
    weightmax = cumweights[cbind(wc[indn1], indn1)]
    weightmin = cumweights[cbind(wc[indn1]-1, indn1)]
    factor = numeric(length(indn1))
    indz = weightmax-weightmin<10^(-10)
    factor[indz] = 0.5
    factor[!indz] = (quantiles[qc]-weightmin[!indz])/(weightmax[!indz]-weightmin[!indz])
    quant[indn1,qc] = quantmin + factor* (quantmax-quantmin)
  }
  colnames(quant) = paste("quantile=", quantiles)
  return(quant)
}