#'Estimating posterior probabilities of single edges
#'
#'This function estimates the posterior probabilities of edges by averaging over a sample of DAGs
#'obtained via an MCMC scheme.
#'
#'@param MCMCchain list of square matrices with elements in \code{\{0,1\}} and representing adjacency matrices of a sample of DAGs obtained via an MCMC scheme
#'(objects of classes 'MCMCtrace' or 'MCMCres' are also valid data types)
#'@param pdag logical, if TRUE (FALSE by default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin (optional) number between \code{0} and \code{1}, indicates the percentage of the samples which will be discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@param endstep (optional) number between \code{0} and \code{1}; 1 by default 
#'@return a square matrix with dimensions equal to the number of variables; each entry \code{[i,j]} is an estimate of the posterior probability of the edge from node \code{i} to node \code{j}
#'@examples
#'Bostonscore<-scoreparameters(14, "bge", Boston)
#'\dontrun{
#'samplefit<-orderMCMC(Bostonscore, iterations=25000,chainout=TRUE)
#'edgesposterior<-edges.posterior(samplefit, pdag=TRUE, burnin=0.2)
#'}
#'@export
edges.posterior<-function(MCMCchain,pdag=FALSE,burnin=0.2,endstep=1) {
  if(!is.matrix(MCMCchain[[1]])) {
    if(class(MCMCchain)=="MCMCres") {
        MCMCchain<-MCMCchain$chain$incidence
    } else if (class(MCMCchain)=="MCMCtrace") {
       MCMCchain<-MCMCchain$incidence
    } else {
      stop("bad input format! no list of matrices found")
    }
  }
  if(endstep==1) {
  endstep<-length(MCMCchain)
  } else {
    endstep<-ceiling(length(MCMCchain)*endstep)
  }
  startstep<-max(as.integer(burnin*endstep),1)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    return(Reduce('+', cpdags)/(endstep-startstep+1))
  } else {
    return(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1))
  }
}


#'Estimating a graph corresponding to a posterior probability threshold
#'
#'This function constructs a directed graph (not necessarily acyclic) including all edges with a posterior probability above a certain threshold.  The posterior probability is evaluated as the Monte Carlo estimate from a sample of DAGs obtained via an MCMC scheme.
#'
#'@param MCMCchain list of adjacency matrices with dimensions equal to n and elements in \code{\{0,1\}}, representing a sample of DAGs from an MCMC scheme 
#'(objects of classes 'MCMCtrace' or 'MCMCres' are also valid data types)
#'@param pbarrier threshold such that only edges with a higher posterior probability will be retained in the directed graph summarising the sample of DAGs
#'@param pdag logical, if TRUE (FALSE by default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin (optional)  number between \code{0} and \code{1}, indicates the percentage of the samples which will be  the discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@return a square matrix with dimensions equal to the number of variables representing the adjacency matrix of the directed graph summarising the sample of DAGs
#'@examples
#'Bostonscore<-scoreparameters(14, "bge", Boston)
#'\dontrun{
#'orderfit<-orderMCMC(Bostonscore, MAP=FALSE, iterations=25000, chainout=TRUE)
#'hdag<-dag.threshold(orderfit, pbarrier=0.9)
#'}
#'@export
dag.threshold<-function(MCMCchain, pbarrier, pdag=FALSE, burnin=0.2) {
  if(!is.matrix(MCMCchain[[1]])) {
    if(class(MCMCchain)=="MCMCres") {
      MCMCchain<-MCMCchain$chain$incidence
    } else if (class(MCMCchain)=="MCMCtrace") {
      MCMCchain<-MCMCchain$incidence
    } else {
      stop("bad input format! no list of matrices found")
    }
  }
  varlabels<-colnames(MCMCchain[[1]])
  n<-nrow(MCMCchain[[1]])
  incidence<-matrix(rep(0, n*n), nrow=n, ncol=n)
  endstep<-length(MCMCchain)
  startstep<-max(as.integer(burnin*endstep),1)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    incidence[which(Reduce('+', cpdags)/(endstep-startstep+1)>pbarrier)]<-1
  } else {
  incidence[which(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1)>pbarrier)]<-1
  }
  colnames(incidence)<-varlabels
  rownames(incidence)<-varlabels
  return(incidence)
}


#'Performance assessment of iterative MCMC scheme against a known Bayesian network
#'
#'This function calculates the number of true and false positives, the true positive rate, the structural Hamming distance and score for each iteration in the search procedure implemented in the function \code{\link{iterativeMCMC}}.
#'
#'@param MCMCmult an object which of class \code{MCMCmult} (output of the function \code{\link{iterativeMCMC}})
#'@param truedag ground truth DAG which generated the data used in the search procedure; represented by an object of class \code{\link[graph]{graphNEL}}
#@param sample logical (FALSE by default), indicates if \code{MCMCmult} contains sample or maximum score DAGs
#'@param cpdag logical, if TRUE (FALSE by default) all DAGs in the \code{MCMCmult} are first converted to their respective equivalence class (CPDAG) before the averaging if parameter \code{sample} set to TRUE
#'@param pbarrier threshold such that only edges with a higher posterior probability will be retained in the directed graph summarising the sample of DAGs at each iteration from \code{MCMCmult} if parameter \code{sample} set to TRUE
#'@param trans logical, for DBNs indicates if model comparions are performed for transition structure; when \code{trans} equals FALSE the comparison is performed for initial structures of estimated models and the ground truth DBN; for usual BNs the parameter is disregarded
#'@return A matrix with the number of rows equal to the number of elements in \code{MCMCmult}, and 5 columns reporting for 
#'the maximally scoring DAG uncovered at each iteration (or for a summary over the sample of DAGs if \code{sample} parameter set to TRUE) 
#'the number of true positive edges (`TP'), the number of false positive edges (`FP'), 
#'the true positive rate (`TPR'), the structural Hamming distance (`SHD') and the score of the DAG (`score'). 
#'Note that the maximum estimated DAG as well as the true DAG are first converted to 
#'the corresponding equivalence class (CPDAG) when calculating the SHD.
#' @examples
#' gsim.score<-scoreparameters(100, "bge", gsim)
#' \dontrun{
#' MAPestimate<-iterativeMCMC(gsim.score)
#' iterations.check(MAPestimate, gsimmat)
#' }
#'@export
iterations.check<-function(MCMCmult, truedag, cpdag=TRUE, pbarrier=0.5,trans=TRUE) {
  TP<-vector()
  FP<-vector()
  TPR<-vector()
  SHD<-vector()
  SC<-vector()

    if(MCMCmult$info$DBN) { #we need to extract either transition or initial structure
     
       if(trans) {
        if(!is.matrix(truedag)) truedag<-graph2m(truedag)
        truedag<-m2graph(DBNcut(truedag,n.dynamic=MCMCmult$info$nsmall,n.static=MCMCmult$info$bgn))
        trueskeleton<-graph2skeleton(truedag)
      } else {
        if(!is.matrix(truedag)) truedag<-graph2m(truedag)
        truedag<-m2graph(DBNinit(truedag,n.dynamic=MCMCmult$info$nsmall,n.static=MCMCmult$info$bgn))
        trueskeleton<-graph2skeleton(truedag)
      }
      
      if(MCMCmult$info$split) {
        if(trans) {
          MCMCmult$maxtrace<-MCMCmult$trans$maxtrace
          MCMCmult$trans$maxtrace<-NULL
        } else {
          MCMCmult$maxtrace<-MCMCmult$init$maxtrace
          MCMCmult$init$maxtrace<-NULL
        }
      } else {
          newtrace<-lapply(MCMCmult$maxtrace,function(x)x$DAG)
          if(trans) {
            newtrace<-lapply(newtrace,DBNcut,n.dynamic=MCMCmult$info$nsmall,n.static=MCMCmult$info$bgn)
            for(i in 1:length(newtrace)) {
              MCMCmult$maxtrace[[i]]$DAG<-newtrace[[i]]
            }
          } else {
            newtrace<-lapply(newtrace,DBNinit,n.dynamic=MCMCmult$info$nsmall,n.static=MCMCmult$info$bgn)
            for(i in 1:length(newtrace)) {
              MCMCmult$maxtrace[[i]]$DAG<-newtrace[[i]]
            }
            }
      }
    } else {
      if(is.matrix(truedag)) truedag<-m2graph(truedag)
      trueskeleton<-graph2skeleton(truedag)
    }
  
  numedges<-sum(trueskeleton)
  chainl<-length(MCMCmult$maxtrace)
  
    for (j in  1:chainl) {
      SC[j]<-MCMCmult$maxtrace[[j]]$score
      maxadj<-MCMCmult$maxtrace[[j]]$DAG
      estskelmcmc<-graph2skeleton(maxadj) #output is an upper triangular matrix, so works for DBNs too
      diffmcmc<-estskelmcmc-trueskeleton
      mc.dag<-m2graph(maxadj)
      cpmcmc<-pcalg::dag2cpdag(mc.dag)
      truecp<-pcalg::dag2cpdag(truedag)
      TP[j]<-numedges-sum(diffmcmc<0)
      FP[j]<-sum(diffmcmc>0)
      SHD[j]<-pcalg::shd(truecp, cpmcmc)
    }
    TPR<-TP/numedges
    result<-cbind(TP, FP, TPR, SHD, SC)
    colnames(result)<-c("TP", "FP", "TPR", "SHD", "score")
    return(result) 
}


#'Performance assessment of sampling algorithms against a known Bayesian network
#'
#'This function calculates the number of true and false positives and the structural Hamming distance between a ground truth DAG and a directed graph summarising a sample of DAGs obtained from an MCMC scheme, as the posterior probability threshold is varied
#'
#'@param MCMCchain an object of class MCMCres, representing the output of structure sampling function \code{\link{partitionMCMC}} or \code{\link{orderMCMC}} (the latter when parameter \code{chainout}=TRUE)
#'@param truedag ground truth DAG which generated the data used in the search procedure; represented by an object of class  \code{\link[graph]{graphNEL}}
#'@param pbarrier (optional) a vector of numeric values between 0 and 1, defining posterior probabilities according to which the edges of assessed structures are drawn, please note very low barriers can lead to very dense structures; by default 
#'\eqn{pbarrier=c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)}
#'@param pdag logical, if TRUE (default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin (optional)  number between \code{0} and \code{1}, indicates the percentage of the samples which will be  the discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@param trans logical, for DBNs indicates if model comparions are performed for transition structure; when \code{trans} equals FALSE the comparison is performed for initial structures of estimated models and the ground truth DBN; for usual BNs the parameter is disregarded
#'@return A matrix with the number of rows equal to the number of posterior thresholds tested, and 4 columns reporting for each  thresholded directed graphs the number of true positive edges (`TP'), the number of false positive edges (`FP'), the structural Hamming distance (`SHD') and the posterior threshold
#' @examples
#' gsim.score<-scoreparameters(100, "bge", gsim)
#' \dontrun{
#' mapest<-iterativeMCMC(gsim.score)
#' ordersample<-orderMCMC(gsim.score, MAP=FALSE, startspace=mapest$endspace)
#' sample.check(ordersample, gsimmat)
#' }
#'@export
sample.check<-function(MCMCchain, truedag, pbarrier=c(0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2),
                                         pdag=TRUE, burnin=0.2, trans=TRUE) {
  if(is.matrix(truedag)) truedag<-m2graph(truedag)
  MCMCmatlist<-MCMCchain$chain$incidence
  n<-nrow(MCMCmatlist[[1]])
  truecp<-pcalg::dag2cpdag(truedag)
  if(MCMCchain$info$DBN) {
    if(trans==TRUE) {
        print("comparison is performed for transition structures")
        trueadj<-DBNcut(graph2m(truedag),MCMCchain$info$nsmall,MCMCchain$info$bgn)
        truedag<-m2graph(trueadj)
        truecpadj<-graph2m(truecp)
        truecpadj<-DBNcut(truecpadj,MCMCchain$info$nsmall,MCMCchain$info$bgn)
        truecp<-m2graph(truecpadj) 
      } else {
          print("comparison is performed for initial structures")
          trueadj<-DBNinit(graph2m(truedag),MCMCchain$info$nsmall,MCMCchain$info$bgn)
          truedag<-m2graph(trueadj)
          truecpadj<-DBNinit(graph2m(truecp),MCMCchain$info$nsmall,MCMCchain$info$bgn)
          truecp<-m2graph(truecpadj)
          n<-MCMCchain$info$nsmall+MCMCchain$info$bgn
      }
  }
  trueskeleton<-graph2skeleton(truedag)
  results<-matrix(ncol=6,nrow=length(pbarrier))
  results[,6]<-pbarrier
  numedges<-sum(trueskeleton)
  endstep<-length(MCMCmatlist)
  startstep<-max(as.integer(burnin*endstep),1)
  
  if(pdag) {
  dags<-lapply(MCMCmatlist[startstep:endstep],dagadj2cpadj) #first convert every DAG in the sample to equivalence class
  } else {dags<-MCMCmatlist[startstep:endstep]}
  
  if(MCMCchain$info$DBN){
    if(trans) {
    dags<-lapply(dags,DBNcut,n.dynamic=MCMCchain$info$nsmall,n.static=MCMCchain$info$bgn)
    trueskeleton<-DBNcut(trueskeleton,MCMCchain$info$nsmall,MCMCchain$info$bgn)
    } else {
      dags<-lapply(dags,DBNinit,n.dynamic=MCMCchain$info$nsmall,n.static=MCMCchain$info$bgn)
      trueskeleton<-DBNinit(trueskeleton,MCMCchain$info$nsmall,MCMCchain$info$bgn)
    }
  }
  postprobmat<-Reduce('+', dags)/(endstep-startstep+1)
 
  for (p in 1:length(pbarrier)) {
    sampledag<-matrix(0, nrow=n,ncol=n)
    sampledag[which(postprobmat>pbarrier[p])]<-1 #average over graphs
    sampleest<-graph2skeleton(sampledag)
    diffmcmc<-sampleest-trueskeleton
    mc.dag<-m2graph(sampledag)

    results[p,1]<-numedges-sum(diffmcmc<0)
    results[p,2]<-sum(diffmcmc>0)
    results[p,4]<-results[p,1]/numedges
    results[p,5]<-results[p,2]/numedges
    
    
    if(pdag) {
        results[p,3]<-pcalg::shd(truecp, mc.dag)
    } else {
      results[p,3]<-pcalg::shd(truedag, mc.dag)
    }
  }
  colnames(results)<-c("TP", "FP", "SHD","TPR","FPRn","post.thr.")
  return(results)
}
