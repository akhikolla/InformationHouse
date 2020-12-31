#' @include generics.R empiricalCopula.R
NULL

.Cort = setClass(Class = "Cort", contains = "empiricalCopula",
                        slots = c(p_value_for_dim_red = "numeric",
                                  number_max_dim = "numeric",
                                  min_node_size = "numeric",
                                  verbose_lvl="numeric",
                                  vols = "numeric",
                                  f = "numeric",
                                  p = "numeric",
                                  a = "matrix",
                                  b = "matrix"),
                        validity = function(object) {
                          errors <- c()
                          if(!(length(dim(object@data))==2)){
                            errors <- c(errors, "data should be a matrix")
                          }
                          if(!(all(object@data >= 0) & all(1 > object@data))){
                            errors <- c(errors,"data should be a matrix inside 0,1")
                          }
                          if(object@number_max_dim > ncol(object@data)){
                            errors <- c(errors, "the maximum number of splitting dimensions should be smaller than the number of dimensons!")
                          }
                          if(object@number_max_dim < 2){
                            errors <- c(errros, "splits cannot be done in dimmensions smaller than 2.")
                          }
                          if (length(errors) == 0)
                            TRUE else errors
                        })

#' Cort class
#'
#' This class implements the CORT algorithm to a fit a multivariate copula using piece constant density. Given a dataset `x`, the function will produce an estimator for the copula of this dataset
#' that is tree-shaped, by recursive partitioning of the unit hypercube. the `min_node_size` parameter controls the stopping conditions for the splitting procedure. Once the space is splitted,
#' we ran a quadratic solver, which options can be tweaked via the `osqp_options` parameter, to ensure that the weights respect the copula conditions.
#'
#'
#' Once the model is fitted, it can be used through the classical (r/d/p/v)Copula functions to compute, respectively, random number generations, the density, the cdf and the volume function of the copula.
#'
#'
#' See O. Laverny, E. Masiello, V. Maume-Deschamps and D. Rullière (2020) for the details of this density estimation procedure, and `vignettes(package='cort')` for examples of usecases.
#'
#'
#' @param x The data, must be provided as a matrix with each row as an observation.
#' @param p_value_for_dim_red a p_value for the localized dimension reduction test
#' @param min_node_size The minimum number of observation available in a leaf to initialize a split.
#' @param pseudo_data set to True if you are already providing data on the copula space.
#' @param number_max_dim The maximum number of dimension a split occurs in. Defaults to be all of the dimensions.
#' @param slsqp_options options for nloptr::slsqp to find breakpoints : you can change defaults.
#' @param verbose_lvl numeric. set the verbosity. 0 for no output and bigger you set it the most output you get.
#' @param N The number of bootstrap samples for p_values computations.
#' @param osqp_options options for the weights optimization. You can pass a call to osqp::osqpSettings, or NULL for defaults.
#' @param force_grid Set to TRUE to force breakpoints to be on the n-checkerboard grid.
#'
#' @name Cort-Class
#' @title Cort copulas
#' @rdname Cort-Class
#'
#' @return An instance of the `Cort` S4 class. The object represent the fitted copula and can be used through several methods to query classical (r/d/p/v)Copula methods, constraint influence, etc.
#' Beside returning some inputted parameters, notable slots are :
#'
#' \itemize{
#'  \item{`data` }{Your original data}
#'  \item{`dim` }{The dimension of problem, number of columns of your dataset}
#'  \item{`f` }{The empirical frequency in the leaves}
#'  \item{`p` }{The fitted probabilities of each leaf}
#'  \item{`a` }{Minimum points of leaves}
#'  \item{`b` }{Maximum points of leaves}
#'  \item{`vols` }{Volume of the leaves}
#' }
#'
#' More details about these slots can be found in the reference.
#'
#' @export
#'
#' @references
#' \insertRef{laverny2020}{cort}
#'
#' @examples
#' (Cort(LifeCycleSavings[,1:3]))
Cort = function(x,
                p_value_for_dim_red=0.75,
                min_node_size=1,
                pseudo_data=FALSE,
                number_max_dim=NULL,
                verbose_lvl=1,
                slsqp_options = NULL,
                osqp_options = NULL,
                N = 999,
                force_grid=FALSE) {

  # Coerce the data :
  data= as.matrix(x)
  row.names(data) = NULL
  colnames(data) = NULL
  if(!pseudo_data){
    data = apply(data,2,rank,ties.method="random")/(nrow(data)+1)
  }

  # initialisation
  d              = ncol(data)
  n_obs          = nrow(data)
  vols = f = p   = 1
  a              = matrix(0,ncol=d,nrow=1)
  b              = matrix(1,ncol=d,nrow=1)
  number_max_dim = min(number_max_dim,d)
  dd             = list(1:n_obs) # index of the data points in each leave
  ss             = list(1:d) # dimensions of splits in each leave

  # Deal with solver parameters :
  DEFAULT_SLQP_OPTIONS = list(stopval = -Inf, xtol_rel = 1e-4, maxeval = 100000, ftol_rel = 1e-6, ftol_abs = 1e-6)
  slsqp_options = purrr::map2(DEFAULT_SLQP_OPTIONS,names(DEFAULT_SLQP_OPTIONS),function(i,n){
    ifelse(is.null(slsqp_options[[n]]),i,slsqp_options[[n]])
  })
  if(is.null(osqp_options)){
    if (verbose_lvl>1) {
      osqp_options = osqp::osqpSettings(max_iter = 100000L, eps_abs = 0.000001, eps_rel = 0.000001, eps_prim_inf = 0.000001, eps_dual_inf = 0.000001, verbose = TRUE)
    } else {
      osqp_options = osqp::osqpSettings(max_iter = 100000L, eps_abs = 0.000001, eps_rel = 0.000001, eps_prim_inf = 0.000001, eps_dual_inf = 0.000001, verbose = FALSE)
    }
  }

  if(verbose_lvl>0) {cat("Splitting...\n")}

  # Loop until there are no more splittable leaves :
  while(TRUE){

    are_splittables = purrr::map_lgl(1:nrow(a),~(length(dd[[.x]])>min_node_size) && (length(ss[[.x]])>1))
    if(!any(are_splittables)){
      break # This is the only way out of the while loop
    }

    if(verbose_lvl>0){
      cat("\n    ",sum(are_splittables),"leaves to split...")
      if(verbose_lvl>1){cat("\n")}
    }

    leaves_to_remove = numeric()
    for(i_leaf in which(are_splittables)){

      spltd = ss[[i_leaf]]
      di    = dd[[i_leaf]]

      if(verbose_lvl>1){
        cat(paste0("        Leaf with ",length(di)," points.\n"))
        if(verbose_lvl>2){
          verb_df            = data.frame(min = a[i_leaf,], max = b[i_leaf,])
          verb_df$bp         = rep(NaN,d)
          verb_df$p_value    = rep(NaN,d)
          verb_df$action     = rep("",d)
          verb_df$reason     = rep("",d)
          row.names(verb_df) = paste0("             ",1:d)
        }
      }

      # Randomize splitting dimensions :

      if(length(spltd) > number_max_dim){
        random_dims    = resample(x =spltd, size=number_max_dim, replace=FALSE)
        non_taken_dims = spltd[!(spltd %in% random_dims)]
        spltd   = random_dims
      } else{
        non_taken_dims = numeric()
      }

      if(verbose_lvl>2){
        verb_df$action[non_taken_dims] = "Dissmissed"
        verb_df$reason[non_taken_dims] = "Randomly"
      }

      d_split = length(spltd)
      if(d_split>1){

        # Compute prerequisites for the optimisation of the breakpoint
        n = length(di)
        leaf_a = a[i_leaf,spltd]
        leaf_b = b[i_leaf,spltd]
        z = (t(data[di,spltd,drop=FALSE]) - leaf_a)/(leaf_b-leaf_a) # d*n
        bin_repr = sapply(1:(2^d_split),number2binary,d_split)

        #Launche optimisation routine for the breakpoint :
        optimizer = nloptr::slsqp(
          x0 = rowMeans(z),
          fn = lossFunc, # coded in Rcpp
          lower = rep(0,d_split),
          upper = rep(1,d_split),
          nl.info=verbose_lvl>3,
          control=slsqp_options,
          bin_repr = bin_repr,
          z = z)

        # Get the breakpoints and the final splitting :
        z_bp = optimizer$par
        bp = leaf_a + z_bp*(leaf_b-leaf_a)

        if(force_grid){
          bp = round(bp*2*(n_obs+1))/(n_obs+1)/2
          z_bp = (bp-leaf_a)/(b-leaf_a)
        }

        # Are we going to keep the breakpoint in every splitting dimensions ?
        z_min                     = z_bp*bin_repr
        z_max                     = z_bp^(1-bin_repr)
        p_values                  = cortMonteCarlo(z,z_min,z_max,as.integer(N))
        p_values[is.na(p_values)] = 0 # really usefull ?
        p_val_too_big             = p_values > p_value_for_dim_red
        threshold                 = 1/ifelse(force_grid,(n_obs+1)^2,min((length(di)+1)^2,1000))
        close_to_bound            = (z_bp< threshold) + (z_bp > 1-threshold) > 0
        to_be_removed             = p_val_too_big+close_to_bound>0

        if(verbose_lvl>2){
          verb_df$bp[spltd]                     = bp
          verb_df$p_value[spltd]                = p_values
          verb_df$action[spltd[p_val_too_big]]  = "Removed"
          verb_df$reason[spltd[p_val_too_big]]  = "Independence test"
          verb_df$action[spltd[close_to_bound]] = "Removed"
          verb_df$reason[spltd[close_to_bound]] = "Close to boundary"
        }

        spltd = spltd[!to_be_removed]
        bp    = bp[!to_be_removed]

        if(length(spltd) < 2){
          if(verbose_lvl>2){
            verb_df$action[spltd] = "Dissmissed"
            verb_df$reason[spltd] = "Split dim < 2"
          }
          ss[[i_leaf]] = c(spltd,non_taken_dims)
        } else {
          # NOW WE FINALY SPLIT
          leaves_to_remove = c(leaves_to_remove,i_leaf)
          if(verbose_lvl>2){verb_df$action[spltd] = "Splitted"}

          # remove the breakpoint from the data points if it's one of them :
          are_the_breakpoint  = (colSums(t(data[di,spltd,drop=FALSE]) == bp) == length(spltd))
          if(any(are_the_breakpoint)){
            if(verbose_lvl>4){cat("Be carrefull, we are splitting on a point.\n")}
            di = di[!are_the_breakpoint]
          }

          # construct new information for new leaves :
          d_split       = length(spltd)
          D             = 2^d_split
          bin_repr      = sapply(1:D,number2binary,d_split)
          new_a         = a[rep(i_leaf,D),]
          new_b         = b[rep(i_leaf,D),]
          new_a[,spltd] = t(bp^bin_repr * a[i_leaf,spltd]^(1-bin_repr))
          new_b[,spltd] = t(bp^(1-bin_repr) * b[i_leaf,spltd]^bin_repr)

          # store new edges and new split dims :
          a = rbind(a,new_a)
          b = rbind(b,new_b)
          ss = c(ss,rep(list(c(spltd,non_taken_dims)),D))
          # store new points assignements :
          for (i in 1:D){
            dd = c(dd,list(di[colSums((t(data[di, spltd, drop=FALSE]) >= new_a[i,spltd]) *
                                                  (t(data[di, spltd, drop=FALSE]) <  new_b[i,spltd]))==d_split]))
          }
        }
      }

      if(verbose_lvl>2) {
        cat(toString.data.frame(verb_df,digits=8))
        cat("\n\n")
      }
    }
    # remove information from the splitted leaves :
    if(length(leaves_to_remove) >= 1){
      a = a[-leaves_to_remove,]
      b = b[-leaves_to_remove,]
      ss = ss[-leaves_to_remove]
      dd = dd[-leaves_to_remove]
    }
  }

  if(verbose_lvl>0) {cat("\nEnforcing constraints...\n")}
  n_leaves = nrow(a)
  vols     = apply(b-a,1,prod)
  f        = purrr::map_dbl(dd,~length(.x))/nrow(data)
  eval_pts = purrr::map(1:d,~unique((b+a)[,.x]/2))
  dims     = unlist(purrr::map2(eval_pts,1:d,function(x,y){rep(y,length(x))}))
  F_vec    = unlist(eval_pts)
  lambdas  = pmin(pmax((F_vec - t(a[,dims,drop=FALSE]))/(t(b[,dims,drop=FALSE] - a[,dims,drop=FALSE])),0),1)
  model    = osqp::osqp(P=diag(1/vols),
                        q=-f/vols,
                        A=rbind(lambdas,rep(1,n_leaves),diag(n_leaves)),
                        l=c(F_vec,1,rep(0,n_leaves)),
                        u=c(F_vec,1,rep(Inf,n_leaves)),
                        pars=osqp_options)
  model$WarmStart(x=f)
  rez = model$Solve()$x
  p = pmax(rez,0)/sum(pmax(rez,0)) # correction for small negative weights.

  # remove unnecessary names :
  row.names(a) <- NULL
  row.names(b) <- NULL
  names(vols) <- NULL

  if(verbose_lvl>0){cat("Done !\n")}
  return(.Cort(data = data, p_value_for_dim_red = p_value_for_dim_red,
               number_max_dim = number_max_dim, min_node_size = min_node_size,
               verbose_lvl=verbose_lvl, dim = d, vols = vols, f = f, p = p, a = a, b = b))
}

setMethod(f = "show", signature = c(object = "Cort"), definition = function(object){
  cat(paste0("Cort copula model: ",nrow(object@data),"x",ncol(object@data),"-dataset and ",
             nrow(object@a)," leaves."))
})


#' @describeIn rCopula-methods Method for the class Cort
setMethod(f = "rCopula", signature = c(n = "numeric", copula = "Cort"), definition = function(n, copula) {
  # we sample the boxes and then sample from them.
  sampled_indexes = resample(1:nrow(copula@a),size=n,prob = copula@p,replace=TRUE)
  matrix(runif(n*copula@dim,min = as.vector(copula@a[sampled_indexes,]),max=as.vector(copula@b[sampled_indexes,])),nrow=n,ncol=copula@dim)
})

#' @describeIn pCopula-methods Method for the class Cort
setMethod(f = "pCopula", signature = c(u = "matrix",  copula = "Cort"), definition = function(u, copula) {
  # The implementation is in Rcpp.
  pCort(copula@a,
        copula@b,
        copula@p,
        u)
})

#' @describeIn dCopula-methods Method for the class Cort
setMethod(f = "dCopula", signature = c(u = "matrix",  copula="Cort"),   definition = function(u, copula) {
  # The implementation is in Rcpp
  dCort(copula@a,
        copula@b,
        copula@p/copula@vols,
        u)
})

#' @describeIn biv_rho-methods Method for the class Cort
setMethod(f = "biv_rho", signature = c(copula="Cort"),   definition = function(copula) {
  # The implementation is in Rcpp.
  bivRho(copula@a,
         copula@b,
         copula@p)
})

#' @describeIn biv_tau-methods Method for the class Cort
setMethod(f = "biv_tau", signature = c(copula="Cort"),   definition = function(copula) {
  # The implmeentation is in Rcpp
  bivTau(copula@a,
         copula@b,
         copula@p)
})

#' @describeIn loss-methods Method for the class Cort
setMethod(f = "loss", signature = c(object="Cort"),   definition = function(object) {
  return(sum((object@p^2/2 - object@p*object@f)/object@vols))
})

#' @describeIn constraint_infl-methods Method for the class Cort
setMethod(f = "constraint_infl", signature = c(object="Cort"),   definition = function(object) {
  return(sum(((object@p - object@f)^2/object@vols)))
})

#' @describeIn quad_norm-methods Method for the class Cort
setMethod(f = "quad_norm", signature = c(object="Cort"),   definition = function(object) {
  return(sum(object@p^2/object@vols))
})

#' @describeIn quad_prod_with_data-methods Method for the class Cort
setMethod(f = "quad_prod_with_data", signature = c(object="Cort"),   definition = function(object) {
  return(sum(object@f * object@p/object@vols))
})

#' @describeIn quad_prod-methods Method for the class Cort
setMethod(f = "quad_prod", signature = c(object="Cort",other_tree = "Cort"),   definition = function(object,other_tree) {

  # The implementation is in Rcpp
  return(quadProd(
    object@a,
    object@b,
    object@p/object@vols,
    other_tree@a,
    other_tree@b,
    other_tree@p/other_tree@vols
  ))
})

#' @param M the number of simulations
#'
#' @describeIn kendall_func-methods Method for the class Cort
setMethod(f = "kendall_func", signature = c(object="Cort"),   definition = function(object,t,M=1000) {

  # Simulate a bunch of random variables :
  m = length(t)
  rng = pCopula(rCopula(M,object),object)
  dim(rng) = c(M,1)
  dim(t) = c(1,m)
  colMeans(rng[,rep(1,m),drop=FALSE]<=t[rep(1,M),,drop=FALSE])

})

#' @describeIn project_on_dims-methods Method for the class Cort
setMethod(f = "project_on_dims", signature = c(object="Cort"),   definition = function(object,dims) {

  # The implementation is in Rcpp
  cpp_result = projectOnTwoDims(a=t(object@a),
                                b=t(object@b),
                                p=object@p,
                                f=object@f,
                                dims=dims,
                                data = object@data[,dims],
                                kern = object@p/object@vols)

  object@dim = 2
  object@data = object@data[,dims]
  object@f = cpp_result$f
  object@p = cpp_result$p
  object@a = cpp_result$a
  object@b = cpp_result$b
  object@vols = cpp_result$vols
  return(object)
})

#' @export
plot.Cort <- function(x,...){

  d = ncol(x@data)
  dd = do.call(rbind,unlist(purrr::map(1:(d-1),function(i){
    purrr::map((i+1):(d),function(j){
      proj = project_on_dims(x,c(i,j))
      rbind(
        data.frame(xmin = proj@a[,1],
                   ymin = proj@a[,2],
                   xmax = proj@b[,1],
                   ymax = proj@b[,2],
                   weight = proj@p,
                   volume = proj@vols,
                   dim_x = rep(i,nrow(proj@a)),
                   dim_y = rep(j,nrow(proj@a))),
        data.frame(ymin = proj@a[,1],
                   xmin = proj@a[,2],
                   ymax = proj@b[,1],
                   xmax = proj@b[,2],
                   weight = proj@p,
                   volume = proj@vols,
                   dim_y = rep(i,nrow(proj@a)),
                   dim_x = rep(j,nrow(proj@a)))
      )
    })}),recursive=FALSE))
    dd$col = log(1+dd$weight/dd$volume)
    dd$col = dd$col/max(dd$col)

    opar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))
    graphics::par(mfrow=c(d,d),
      mai = c(0,0,0,0),
      oma = c(3,3,5,3))

  for (i in 1:d){
    for (j in 1:d){
      graphics::plot(c(0,1),c(0,1),type="n",xlab="",ylab="",xaxt='n',yaxt='n')
      if(i != j){
        xx = dd[(dd$dim_x == i)&(dd$dim_y==j),]
        graphics::rect(xx$ymin,xx$xmin,xx$ymax,xx$xmax,col=grDevices::gray(1-xx$col),border=NA,density=NA)
        graphics::points(x@data[,j],x@data[,i],cex=0.5,col="red")
      }
    }
  }
}





























