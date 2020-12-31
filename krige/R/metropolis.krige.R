#SAMPLER: PERFORMING METROPOLIS-HASTINGS SAMPLING FOR A LINEAR MODEL SPECIFIED #
#OVER POINT-REFERENCED GEOSPATIAL DATA
### LAST UPDATE: 11/08/2020; Le Bao

#' Sampling Technique Using Metropolis-Hastings
#'
#' This function performs Metropolis-Hastings sampling for a linear model specified 
#' over point-referenced geospatial data. It returns MCMC iterations, with which 
#' results of the geospatial linear model can be summarized.
#'
#' @param formula An object of class \code{\link[stats:formula]{formula}} (or one 
#'   that can be coerced to that classes). A symbolic description of the model to 
#'   be fitted. Alternatively, the model can be specified in \code{y} (a vector of
#'   the outcome variable) and \code{X} (a matrix of explanatory variables).
#' @param coords A matrix of coordinates for all observations or a vector of variable 
#'   names indicating the coordinates variables in the data. Alternatively, the 
#'   coordinates can also be specified seperately using \code{east} and \code{north}.
#' @param data An data frame containing the variables in the model. 
#' @param powered.exp This exponent, which must be greater than 0 and less than 
#'   or equal to 2, specifies a powered exponential correlation structure for the 
#'   data. One widely used specification is setting this to 1, which yields an 
#'   exponential correlation structure. Another common specification is setting 
#'   this to 2 (the default), which yields a Gaussian correlation structure.
#' @param n.iter Number of MCMC iterations (defaults to 100).
#' @param n.burnin Number of iterations that will be discarded for burnin (warmup).
#'   The number of burnin should not be larger than \code{n.iter} and the default is 0. 
#' @param y Alternative specification for the outcome variable that is used in the 
#'   kriging model. If formula is used, this argument will be suppressed.
#' @param X Alternative specification for the matrix of explanatory variables used 
#'   in the kriging model. Different forms of the variables such as transformations 
#'   and interactions also need to be specified accordingly beforehand.
#' @param east Alternative specification for the vector of eastings for all observations.
#' @param north Alternative specification for the vector of northing for all observations.
#' @param na.action A function which indicates what should happen when the data 
#'   contain NAs. The default is "na.fail." Another possible value is "na.omit."
#' @param spatial.share Prior for proportion of unexplained variance that is spatial 
#'   in nature. Must be greater than 0 and less than 1. Defaults to an even split, 
#'   valued at 0.5.
#' @param range.share Prior for the effective range term, as a proportion of the 
#'   maximum distance in the data. Users should choose the proportion of distance 
#'   at which they think the spatial correlation will become negligible. Must be 
#'   greater than 0. Values greater than 1 are permitted, but users should recognize 
#'   that this implies that meaningful spatial correlation would persist outside 
#'   of the convex hull of data. Defaults to half the maximum distance, valued at 0.5.
#' @param beta.var Prior for the variance on zero-meaned normal priors on the 
#'   regression coefficients. Must be greater than 0. Defaults to 10.
#' @param range.tol Tolerance term for setting the effective range. At the distance 
#'   where the spatial correlation drops below this term, it is judged that the 
#'   effective range has been met. The default value is the commonly-used 0.05. 
#'   Must be greater than 0 and less than 1.
#' @param b.tune Tuning parameter for candidate generation of regression coefficients 
#'   that must be greater than 0. A value of 1 means that draws will be based on 
#'   the variance-covariance matrix of coefficients from OLS. Larger steps are taken 
#'   for values greater than 1, and smaller steps are taken for values from 0 to 1. 
#'   Defaults to 1.0.
#' @param nugget.tune Tuning parameter for candidate generation of the nugget term 
#'   (\code{tau2}) that must be greater than 0. A value of 1 means that draws will 
#'   be based on the typical variance of an inverse gamma distribution. \emph{Smaller} 
#'   steps are taken for values \emph{greater} than 1, and \emph{larger} steps are 
#'   taken for \emph{decimal} values from 0 to 1. Defaults to 10.0.
#' @param psill.tune Tuning parameter for candidate generation of the partial sill 
#'   term (\code{sigma2}) that must be greater than 0. A value of 1 means that draws 
#'   will be based on the typical variance of an inverse gamma distribution. 
#'   \emph{Smaller} steps are taken for values \emph{greater} than 1, and \emph{larger} 
#'   steps are taken for \emph{decimal} values from 0 to 1. Defaults to 1.0.
#' @param distance.matrix Logical value indicates whether to save the distance matrix 
#'   in the output object. Saving distance matrix can save time for furthur use such as
#'   in \code{update()} function but may results in larger file size. Defaults to \code{FALSE}.
#' @param progress.bar Types of progress bar. The default is "message" and will 
#'   report variance terms. Other possible values are "TRUE" (simple percentage) 
#'   and "FALSE" (suppress the progress bar).
#' @param accept.rate.warning Logical values indicating whether to show the warnings 
#'   when the acceptance rates are too high or too low. Defaults to \code{TRUE}.
#' 
#' @return An object of class \code{krige} that includes the output MCMC matrix
#'   of sampled values from the posterior distribution as well as the record of 
#'   function arguments, model frame, acceptance rates, and standing parameters.  
#'   
#' @details Analysts should use this function if they want to estimate a linear 
#'   regression model in which each observation can be located at points in geographic 
#'   space. That is, each observation is observed for a set of coordinates in eastings 
#'   & northings or longitude & latitude. 
#'   
#'   Researchers must specify their model in the following manner: \code{formula}
#'   should be a symbolic description of the model to be fitted; it is similar to 
#'   \code{R} model syntax as used in \code{lm()}. In addition, a matrix of 
#'   coordinates must be specified for the geospatial model in \code{coords}. \code{coords}
#'   should be a matrix with two columns that specify west-east and north-south 
#'   coordinates, respectively (ideally eastings and northings but possibly longitude 
#'   and latitude). It can also be a vector of strings indicating the variables names 
#'   of the coordinates in the \code{data}. \code{data} should be a data frame 
#'   containing the variables in the model including both the formula and coordinates 
#'   (if only the names are provided). Alternatively, users can also specify the 
#'   variables using \code{y}, \code{X}, \code{east}, and \code{north} for outcome, 
#'   explanatory, west-east coordinates, and north-south coordinates variables, 
#'   respectively. This alternative specification is compatible with the one used 
#'   in the early version of this package.
#'   
#'   \code{n.iter} is the number of iterations to sample from the posterior distribution 
#'   using the Metropolis-Hastings algorithm. This defaults to 100 iterations, but 
#'   many more iterations would normally be preferred. \code{n.burnin} is set to 0
#'   by default to preserve all the iterations since the kriging model usually takes 
#'   a relatively long time to run. Users can set a number for burnin or use \code{burnin}
#'   function afterwards to discard the burnin period. The output of the function 
#'   prints the proportion of candidate values for the coefficients and for the 
#'   variance terms accepted by the Metropolis-Hastings algorithm. Particularly 
#'   low or high acceptance rates respectively may indicate slow mixing (requiring 
#'   more iterations) or a transient state (leading to nonconvergence), so additional 
#'   messages will print for extreme acceptance rates. Users may want to adjust the 
#'   tuning parameters \code{b.tune}, \code{nugget.tune}, or \code{psill.tune}, 
#'   or perhaps the tolerance parameter \code{range.tol} if the acceptance rate 
#'   is too high or too low.
#'   
#'   The function returns a "krige" list object including the output MCMC matrix
#'   of sampled values from the posterior distribution as well as the record of 
#'   function arguments, model frame, acceptance rates, and standing parameters. 
#'   Users can use the generic \code{summary} function to summarize the results or
#'   extract the elements of the object for further use.  
#'   
#' @references 
#'   Jeff Gill. 2020. Measuring Constituency Ideology Using Bayesian Universal Kriging. 
#'     \emph{State Politics & Policy Quarterly}. \code{doi:10.1177/1532440020930197}
#'   
#' @examples
#' \dontrun{
#' # Summarize example data
#' summary(ContrivedData)
#' 
#' # Initial OLS model
#' contrived.ols<-lm(y~x.1+x.2,data=ContrivedData)
#' # summary(contrived.ols)
#' 
#' # Set seed
#' set.seed(1241060320)
#' 
#' #For simple illustration, we set to few iterations.
#' #In this case, a 10,000-iteration run converges to the true parameters.
#' #If you have considerable time and hardware, delete the # on the next line.
#' #10,000 iterations took 39 min. with 8 GB RAM & a 1.5 GHz Quad-Core processor.
#' M <- 100
#' #M<-10000
#' 
#' contrived.run <- metropolis.krige(y ~ x.1 + x.2, coords = c("s.1","s.2"), 
#'    data = ContrivedData, n.iter = M, n.burnin=20, range.tol = 0.05)
#' # Alternatively, use burnin() after estimation  
#' #contrived.run <- burnin(contrived.run, n.burnin=20)
#' 
#' # Summarize the results and examine results against true coefficients	
#' summary(contrived.run)
#' (TRUTH<-c(0.5,2.5,0.5,0,1,2))
#' 
#' # Extract the MCMC matrix of the posterior distribution
#' contrived.run.mat <- mcmc.samples(contrived.run)
#' head(contrived.run.mat)
#' 
#' # Diagnostics
#' geweke(contrived.run, early.prop=0.5)
#' heidel.welch(contrived.run)
#' 
#' # Semivariogram
#' ### Semivariance
#' semivariance(contrived.run)
#' ### Plot semivariogram
#' semivariogram(contrived.run)
#' ### Alternatively, use generic plot() on a krige object
#' plot(contrived.run)
#' }
#'
#' @importFrom stats rbeta rgamma runif vcov var formula resid model.matrix model.response model.frame lm
#' @importFrom Rcpp evalCpp
#' 
#' @export

metropolis.krige <- function(formula,coords,data,n.iter=100,powered.exp=2,n.burnin=0,
                             y,X,east,north,na.action="na.fail",spatial.share=0.5,
                             range.share=0.5, beta.var=10,range.tol=0.05,b.tune=1.0,
                             nugget.tune=10.0, psill.tune=1.0, distance.matrix=FALSE,
                             progress.bar="message", accept.rate.warning=TRUE){
  
  # ERROR CHECKS
  if(n.iter <= 0) stop("'n.iter' must be greater than 0.")
  if(n.iter%%1 != 0) stop("'n.iter' must be an integer.")
  if(powered.exp<=0 | powered.exp>2) stop("powered.exp must be greater than 0 and less than or equal to 2.")
  if(spatial.share<=0 | spatial.share>=1) stop("spatial.share must be between 0 and 1.")
  if(range.share<=0) stop("range.share must be greater than 0.")
  if(range.tol<=0 | range.tol>=1) stop("p.range.tol must be between 0 and 1.")
  if(beta.var<=0) stop("beta.var must be greater than 0.")
  if(b.tune<=0) stop("b.tune must be greater than 0.")
  if(nugget.tune<=0) stop("nugget.tune must be greater than 0.")
  if(psill.tune<=0) stop("psill.tune must be greater than 0.")
  if (n.burnin < 0) {
    n.burnin <- 0; warning("The burn-in period is negative. 'n.burnin = 0' is used.")
  } else if (n.burnin >= n.iter) {stop("The number of iterations is less than the burn-in period.")}
  
  # IMPUT DATA
  ## Model frame
  cl <- match.call()
  if (missing(formula) || is.null(formula)) formula <- y ~ X - 1
  if (missing(coords) || is.null(coords)) coords <- cbind(east, north)
  if (missing(data) || is.null(data)) data <- model.frame(formula)
  if (is.character(coords) & length(coords) == 2){
    coords.names <- coords
    east <- data[coords[1]]; north <- data[coords[2]]
    coords <- cbind(east, north)
    colnames(coords) <- coords.names
  } 
  if (is.null(colnames(coords))) colnames(coords) <- c("east", "north")
  data <- as.data.frame(cbind(model.frame(formula,data), coords))
  allvar <- update(formula, paste("~ . +",paste(colnames(coords), collapse=" + ")))
  data <- model.frame(formula = allvar, data = data, na.action = na.action)
  coords <- unlist(colnames(coords))
  
  # DATA
  mf <- model.frame(formula = formula, na.action = na.action, 
                    data = data)
  y <- model.response(mf, type = "numeric")
  X <- model.matrix(formula, data = data)
  if(is.null(dimnames(X)[[2]])) {
    mynames<-rep(".",dim(X)[2])
    for(k in 1:dim(X)[2]) mynames[k]<-paste("x.",k,sep="")
    dimnames(X)[[2]]<-mynames
  }
  easting <- as.vector(unlist(data[coords[1]]))
  northing <- as.vector(unlist(data[coords[2]]))

  #DEFINE STANDING PARAMETERS
  dist.mat<-k_distmat(cbind(easting,northing))
  max.distance<-max(dist.mat)
  min.distance<-min(dist.mat[dist.mat>0])
  init.ols <- lm(formula, data=data) # STARTING VALUES
  b.cand.var<-b.tune*vcov(init.ols) 
  err.var<-var(resid(init.ols))
  beta.A<-round(100*range.share/((-log(range.tol))^(1/powered.exp)))
  beta.B<-round(100*(1-(range.share/((-log(range.tol))^(1/powered.exp)))))
  attr(init.ols$terms, ".Environment") <- NULL
  
  #CREATE OUTPUT MATRIX
  mcmc.mat<-matrix(NA,nrow=n.iter,ncol=3+ncol(X))
  dimnames(mcmc.mat)[[2]]<-c("tau2","phi","sigma2",dimnames(X)[[2]])
  start.decay<-((-log(range.tol))^(1/powered.exp))/(max(dist.mat)*range.share)
  mcmc.mat[1,]<-c((1-spatial.share)*err.var,start.decay,spatial.share*err.var,init.ols$coef) 
  local.Sigma<-NULL
  
  #INITIALIZE ACCEPTANCE RATES
  accepted.beta<-0
  accepted.nugget<-0
  accepted.decay<-0
  accepted.psill<-0
  
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
    if(is.null(local.Sigma)) local.Sigma<-ifelse(dist.mat>0, 
                                                 mcmc.mat[i,3]*exp(-abs(mcmc.mat[i,2]*dist.mat)^powered.exp),
                                                 mcmc.mat[i,1]+mcmc.mat[i,3])
    current <- krige.posterior(mcmc.mat[i,1],mcmc.mat[i,2],mcmc.mat[i,3],old.beta,
                               y,X,easting,northing,semivar.exp=powered.exp,
                               p.spatial.share=spatial.share,p.range.share=range.share,
                               p.range.tol=range.tol,p.beta.var=beta.var,tot.var=err.var,
                               local.Sigma=local.Sigma,max.distance=max.distance)
    candidate.beta <- krige.posterior(mcmc.mat[i,1],mcmc.mat[i,2],mcmc.mat[i,3],
                                      beta,y,X,easting,northing,semivar.exp=powered.exp,
                                      p.spatial.share=spatial.share,p.range.share=range.share,
                                      p.range.tol=range.tol,p.beta.var=beta.var,tot.var=err.var,
                                      local.Sigma=local.Sigma,max.distance=max.distance)
    candidate.nugget <- krige.posterior(tau2,mcmc.mat[i,2],mcmc.mat[i,3],old.beta,
                                        y,X,easting,northing,semivar.exp=powered.exp,
                                        p.spatial.share=spatial.share,p.range.share=range.share,
                                        p.range.tol=range.tol,p.beta.var=beta.var,
                                        tot.var=err.var,local.Sigma=NULL,
                                        max.distance=max.distance)
    candidate.decay <- krige.posterior(mcmc.mat[i,1],phi,mcmc.mat[i,3],old.beta,
                                       y,X,easting,northing,semivar.exp=powered.exp,
                                       p.spatial.share=spatial.share,p.range.share=range.share,
                                       p.range.tol=range.tol,p.beta.var=beta.var,
                                       tot.var=err.var,local.Sigma=NULL,
                                       max.distance=max.distance)
    candidate.psill <- krige.posterior(mcmc.mat[i,1],mcmc.mat[i,2],sigma2,old.beta,
                                       y,X,easting,northing,semivar.exp=powered.exp,
                                       p.spatial.share=spatial.share,p.range.share=range.share,
                                       p.range.tol=range.tol,p.beta.var=beta.var,
                                       tot.var=err.var,local.Sigma=NULL,
                                       max.distance=max.distance)
    a.beta<-exp(candidate.beta-current)
    a.nugget<-exp(candidate.nugget-current)
    a.decay<-exp(candidate.decay-current)
    a.psill<-exp(candidate.psill-current)

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
    
    if (!progress.bar==FALSE & i==n.iter-1) 
      progBar(iter=i+1, total=n.iter, progress.bar=progress.bar, 
              nugget=mcmc.mat[i,1], decay=mcmc.mat[i,2], 
              partial.sill=mcmc.mat[i,3], interac = interactive())
  }
  
  beta.rate<-accepted.beta/(nrow(mcmc.mat)-1)
  tau2.rate<-accepted.nugget/(nrow(mcmc.mat)-1)
  phi.rate<-accepted.decay/(nrow(mcmc.mat)-1)
  sigma2.rate<-accepted.psill/(nrow(mcmc.mat)-1)
  ar.rate <- list(beta.rate = beta.rate, tau2.rate = tau2.rate, 
                  phi.rate = phi.rate, sigma2.rate = sigma2.rate)
  if (n.burnin > 0) {mcmc.mat2 <- burnin.matrix(mcmc.mat, n.burnin = n.burnin)
   } else {mcmc.mat2 <- mcmc.mat}
  if (distance.matrix == FALSE) dist.mat <- NULL
  standing <- list(dist.mat = dist.mat, max.distance = max.distance, 
                   min.distance = min.distance, b.cand.var = b.cand.var, 
                   err.var = err.var, beta.A = beta.A, beta.B = beta.B, 
                   local.Sigma = local.Sigma, accepted.beta=accepted.beta, 
                   accepted.nugget=accepted.nugget, accepted.decay=accepted.decay, 
                   accepted.psill=accepted.psill, powered.exp=powered.exp)
  krige.out <- list(call = cl, 
                    formula = formula,
                    coords = coords,
                    n.iter = n.iter,
                    n.burnin = n.burnin,
                    init.ols = init.ols,
                    priors = list(spatial.share = spatial.share, 
                                  range.share = range.share, beta.var = beta.var, 
                                  range.tol = range.tol, b.tune = b.tune, 
                                  nugget.tune = nugget.tune, psill.tune = psill.tune),
                    data = data,
                    model.data.list = list(y = y, X = X, easting = easting, 
                                           northing = northing),
                    standing.parameter = standing,
                    acceptance.rate = ar.rate, 
                    start = mcmc.mat[1,],
                    end = mcmc.mat[nrow(mcmc.mat),],
                    mcmc.mat = mcmc.mat2)
  class(krige.out) <- "krige"
  krige.out
}