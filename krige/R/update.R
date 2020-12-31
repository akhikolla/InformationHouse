#################### UPDATE MCMC RUN FOR KRIGE OBJECT ##########################
### LAST UPDATE: 11/08/2020; Le Bao

#' Update Model
#' 
#' Update the Markov chain associated with the \code{metropolis.krige} model
#' 
#' @param object An \code{krige} object from the \code{metropolis.krige} function.
#' @param n.iter Number of iterations for the update run.
#' @param n.burnin The number of burnin iterations. Defaults to 0.
#' @param combine Logical value indicate whether to combine the update run with 
#'   original output.
#' @param \dots Additional arguments passed to \code{update} methods. Not supported for 
#'   \code{krige} objects.
#'   
#' @return An object of class \code{krige} that includes the output MCMC matrix
#'   of sampled values from the posterior distribution as well as the record of 
#'   function arguments, model frame, acceptance rates, and standing parameters.  
#'   
#' @details Since MCMC calculations typically need to run relatively long. This 
#'   function can continue the MCMC calculations by \code{metropolis.krige()}. 
#'   The parameters of the original model and the estimated results from the previous
#'   run are passed through the \code{krige} object. 
#'   
#'   As geospatial estimation can be computationally taxing, the users may want to 
#'     preserve more iterations of the posterior samples. The \code{combine} argument 
#'     can be used to indicate whether combine the updated run with previous results. 
#'     This includes both the posterior samples and acceptance rates.
#' 
#' @examples
#' \dontrun{
#' # Summarize Data
#' summary(ContrivedData)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' M <- 100
#' #M<-10000
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'    data = ContrivedData, n.iter = M, range.tol = 0.05)
#' summary(contrived.run)
#' 
#' # Update run
#' contrived.run2 <- update(contrived.run, n.iter = M, combine = TRUE)
#' 
#' summary(contrived.run2)
#' #plot(contrived.run2)
#' }
#' 
#' @importFrom stats rbeta rgamma runif vcov var formula update
#' @importFrom Rcpp evalCpp
#' 
#' @export
update.krige <- function (object, n.iter, n.burnin=0, combine=FALSE, ...){
  if (!inherits(object, "krige")) stop("The input object is not a 'krige' object") 
  if(n.iter <= 0) stop("'n.iter' must be greater than 0.")
  if(n.iter%%1 != 0) stop("'n.iter' must be an integer.")
  if (n.burnin < 0) {
    n.burnin <- 0; warning("The burn-in period is negative. 'n.burnin = 0' is used.")
  } else if (n.burnin >= n.iter) {stop("The number of iterations is less than the burn-in period.")}
  
  #INITIALIZE INVISIBLE VARIABLES
  X<-y<-b.cand.var<-nugget.tune<-psill.tune<-max.distance<-beta.A<-beta.B<-dist.mat<-
    powered.exp<-easting<-northing<-spatial.share<-range.share<-range.tol<-beta.var<-
    err.var<-min.distance<-b.tune<-NULL
  #CALL
  object$call$n.iter <- n.iter
  if ("progress.bar" %in% names(object$call)) {
    progress.bar <- object$call$progress.bar
  } else {progress.bar="message"}
  
  #NUMBER OF ITERATIONS
  new.iter <- n.iter
  n.iter <- ifelse(combine==TRUE, new.iter+1, new.iter)
  
  #ATTACH KRIGE OBJECT
  ## Iterations
  old.iter <- object$n.iter
  ## Priors
  for(i in seq(1, length(object$priors))){
    mpriors <- paste0(names(object$priors)[i])
    assign(mpriors, object$priors[[i]])
  }
  
  ## Model matrix
  for(i in seq(1, length(object$model.data.list))){
    mdata <- paste0(names(object$model.data.list)[i])
    assign(mdata, object$model.data.list[[i]])
  }
  
  ## Standing parameters
  for(i in seq(1, length(object$standing.parameter))){
    mpara <- paste0(names(object$standing.parameter)[i])
    assign(mpara, object$standing.parameter[[i]])
  }
  if (!is.null(dist.mat)) {distance.matrix <- TRUE} else {
    dist.mat <- k_distmat(cbind(easting,northing)); distance.matrix <- FALSE}
    
  #CREATE OUTPUT MATRIX
  old.mcmc.mat <- object$mcmc.mat
  mcmc.mat<-matrix(NA,nrow=n.iter,ncol=3+ncol(X))
  dimnames(mcmc.mat)[[2]]<-dimnames(old.mcmc.mat)[[2]]
  mcmc.mat[1,] <- object$end
  
  #INITIALIZE ACCEPTANCE RATES
  if (combine == FALSE){
    accepted.beta<-0; accepted.nugget<-0; accepted.decay<-0; accepted.psill<-0
  } #else {
  #  for(i in seq(1, length(object$acceptance.rate))){
  #    mar <- paste0(names(object$acceptance.rate)[i])
  #    assign(mar, object$acceptance.rate[[i]])
  #  }
  #}
  
  #START OF MARKOV CHAIN
  for (i in 1:(nrow(mcmc.mat)-1))  {
    ### Progress Bar
    if (!progress.bar==FALSE) progBar(iter=i, total=n.iter, progress.bar=progress.bar, 
                                      nugget=mcmc.mat[i,1], decay=mcmc.mat[i,2], 
                                      partial.sill=mcmc.mat[i,3], interac = interactive())
    
    old.beta<-mcmc.mat[i,4:ncol(mcmc.mat)]
    old.var<-var(y-X%*%old.beta)
    beta<-simp.mvrnorm(n=1,mu=old.beta,Sigma=b.cand.var)
    local.share<-mcmc.mat[i,3]/(mcmc.mat[i,1]+mcmc.mat[i,3])
    
    local.tau2.shape<-1+nugget.tune/(1-local.share)
    local.sigma2.shape<-1+psill.tune/local.share
    tau2<-1/rgamma(1,shape=local.tau2.shape,rate=nugget.tune*old.var)
    sigma2<-1/rgamma(1,shape=local.sigma2.shape,rate=psill.tune*old.var)
    
    phi<-1/(max.distance*rbeta(1,shape1=beta.A,shape2=beta.B))
    if(is.null(local.Sigma)) local.Sigma<-ifelse(dist.mat>0, mcmc.mat[i,3]*exp(-abs(mcmc.mat[i,2]*dist.mat)^powered.exp), mcmc.mat[i,1]+mcmc.mat[i,3])
    current <- krige.posterior(mcmc.mat[i,1],mcmc.mat[i,2],mcmc.mat[i,3],old.beta,
                               y,X,easting,northing,semivar.exp=powered.exp,p.spatial.share=spatial.share,
                               p.range.share=range.share,p.range.tol=range.tol,p.beta.var=beta.var,tot.var=err.var,local.Sigma=local.Sigma,max.distance=max.distance)
    candidate.beta <- krige.posterior(mcmc.mat[i,1],mcmc.mat[i,2],mcmc.mat[i,3],beta,
                                      y,X,easting,northing,semivar.exp=powered.exp,p.spatial.share=spatial.share,
                                      p.range.share=range.share,p.range.tol=range.tol,p.beta.var=beta.var,tot.var=err.var,local.Sigma=local.Sigma,max.distance=max.distance)
    candidate.nugget <- krige.posterior(tau2,mcmc.mat[i,2],mcmc.mat[i,3],old.beta,
                                        y,X,easting,northing,semivar.exp=powered.exp,p.spatial.share=spatial.share,
                                        p.range.share=range.share,p.range.tol=range.tol,p.beta.var=beta.var,tot.var=err.var,local.Sigma=NULL,max.distance=max.distance)
    candidate.decay <- krige.posterior(mcmc.mat[i,1],phi,mcmc.mat[i,3],old.beta,
                                       y,X,easting,northing,semivar.exp=powered.exp,p.spatial.share=spatial.share,
                                       p.range.share=range.share,p.range.tol=range.tol,p.beta.var=beta.var,tot.var=err.var,local.Sigma=NULL,max.distance=max.distance)
    candidate.psill <- krige.posterior(mcmc.mat[i,1],mcmc.mat[i,2],sigma2,old.beta,
                                       y,X,easting,northing,semivar.exp=powered.exp,p.spatial.share=spatial.share,
                                       p.range.share=range.share,p.range.tol=range.tol,p.beta.var=beta.var,tot.var=err.var,local.Sigma=NULL,max.distance=max.distance)
    a.beta<-exp(candidate.beta-current)
    a.nugget<-exp(candidate.nugget-current)
    a.decay<-exp(candidate.decay-current)
    a.psill<-exp(candidate.psill-current)
    #if (is.na(a.beta) || is.na(a.nugget) || is.na(a.decay) || is.na(a.psill)) {print(paste("Undefined acceptance ratio at iteration ",i)); a.beta<-a.nugget<-a.decay<-a.psill<-0}
    if (a.beta > runif(1)) {
      accepted.beta <- accepted.beta + 1
      mcmc.mat[(i+1),4:ncol(mcmc.mat)] <- beta
    } else mcmc.mat[(i+1),4:ncol(mcmc.mat)] <- mcmc.mat[i,4:ncol(mcmc.mat)]
    if (a.nugget > runif(1)) {
      accepted.nugget<-accepted.nugget + 1
      mcmc.mat[(i+1),1] <- tau2
      local.Sigma<-NULL
    } else mcmc.mat[(i+1),1] <- mcmc.mat[i,1]
    if (a.decay > runif(1)) {
      accepted.decay<-accepted.decay+1
      mcmc.mat[(i+1),2] <- phi
      local.Sigma<-NULL
    } else mcmc.mat[(i+1),2] <- mcmc.mat[i,2]
    if (a.psill > runif(1)) {
      accepted.psill<-accepted.psill + 1
      mcmc.mat[(i+1),3] <- sigma2
      local.Sigma<-NULL
    } else mcmc.mat[(i+1),3] <- mcmc.mat[i,3]
    
    ### Progress Bar
    if (!progress.bar==FALSE & i==n.iter-1) 
      progBar(iter=i+1, total=n.iter, progress.bar=progress.bar, 
              nugget=mcmc.mat[i,1], decay=mcmc.mat[i,2], 
              partial.sill=mcmc.mat[i,3], interac = interactive())
  }
  
  if (combine==TRUE) mcmc.mat <- rbind(old.mcmc.mat, mcmc.mat[-1,])
  beta.rate<-accepted.beta/(nrow(mcmc.mat)-1)
  tau2.rate<-accepted.nugget/(nrow(mcmc.mat)-1)
  phi.rate<-accepted.decay/(nrow(mcmc.mat)-1)
  sigma2.rate<-accepted.psill/(nrow(mcmc.mat)-1)
  ar.rate <- list(beta.rate = beta.rate, tau2.rate = tau2.rate, 
                  phi.rate = phi.rate, sigma2.rate = sigma2.rate)
  mcmc.mat2 <- ifelse(n.burnin > 0, burnin.matrix(mcmc.mat, n.burnin = n.burnin), mcmc.mat)
  if (distance.matrix == FALSE) dist.mat <- NULL
  standing <- list(dist.mat = dist.mat, max.distance = max.distance, 
                   min.distance = min.distance, b.cand.var = b.cand.var, 
                   err.var = err.var, beta.A = beta.A, beta.B = beta.B, 
                   local.Sigma = local.Sigma, accepted.beta=accepted.beta, 
                   accepted.nugget=accepted.nugget, accepted.decay=accepted.decay, 
                   accepted.psill=accepted.psill, powered.exp=powered.exp)
  krige.out <- list(call = object$call, 
                    formula = object$formula,
                    coords = object$coords,
                    n.iter = nrow(mcmc.mat),
                    n.burnin = n.burnin,
                    init.ols = object$init.ols,
                    priors = list(spatial.share = spatial.share,range.share = range.share, 
                                  beta.var = beta.var,range.tol = range.tol, 
                                  b.tune = b.tune,nugget.tune = nugget.tune, 
                                  psill.tune = psill.tune),
                    data = object$data,
                    model.data.list = list(y = y, X = X, easting = easting, northing = northing),
                    standing.parameter = standing,
                    acceptance.rate = ar.rate, 
                    start = mcmc.mat[1,],
                    end = mcmc.mat[n.iter,],
                    mcmc.mat = mcmc.mat)
  class(krige.out) <- "krige"
  krige.out
  }
