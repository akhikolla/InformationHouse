#### PPass function - brings it all together

#' Person Assessment function
#' 
#' Estimate Person Paramters and calculate Person Fit in one step to gain resonse pattern assessment. Submit a data.frame which contains item responses, or an fitted model (Rasch Model and Partial Credit Model are supported) of the \code{eRm} package.
#' 
#' @param \ldots Submit arguments to the underlying functions: \code{PP_4pl}, \code{PP_gpcm} and \code{PPall} (see documentation files) for person parameter estimation.
#' 
#' @rdname PPass
#' 
#' @example R/.example_ppass.R
#' 
#' @export
#'

PPass <- function(...) UseMethod("PPass")

# ---------------------------------------------------------------------



#' @param respdf A data.frame which contains the items, and perhaps other informations. Each row is a person related resonse patter. Each column denotes a variable.
#' 
#' @param items A numeric (integer) vector which indicates the positions of the items in the data.frame (\code{respdf}). If \code{items} = 'all', \bold{all columns} are treated as items.
#' 
#' @param mod Choose your data generating model. This argument switches between the three person parameter estimating functions \code{PP_4pl}, \code{PP_gpcm} and \code{PPall}.
#' 
#' @param fitindices A character vector which denotes the fit indices to compute.
#' 
#' @details PPass fuses Person Parameter estimation and Person Fit computation into a single function.
#' 
#' @return The original data.frame and
#' 
#' \itemize{
#' 
#' \item The Person Parameter estimates incl. Standard Errors (2 columns)
#' 
#' \item Person Fit Indices you chose (1 or more)
#' 
#' }
#' 
#' @author Manuel Reif, Jan Steinfeld
#' 
#' @method PPass default
#' 
#' @rdname PPass
#' 
#' @export
#' 
#' @seealso \link{PP_4pl}, \link{PP_gpcm}, \link{PPall}, \link{Pfit}
#' 

PPass.default <- function(respdf, items="all", mod=c("1PL","2PL","3PL","4PL","PCM","GPCM","MIXED"), fitindices= c("lz","lzstar","infit","outfit"), ...)
{
  
  # catch additional arguments
  all_pts <- list(respdf,...)
  fitindices <- match.arg(fitindices, several.ok = TRUE)  
  
  ########### ESTIMATE PERSON PARAMETERS ###############################  
  
  ## checks concering the input  
  stopifnot(is.data.frame(respdf)) # muss ein df sein als input
  stopifnot((is.numeric(items) & all(items >= 1)) | all(items == "all")) # indices oder alles
  
  ## create matrix, keep the rest
  if(all(items == "all")) # all variables are items
  {
    respm <- as.matrix(respdf)
  } else {
    
    rest  <- respdf[ , -items, drop=FALSE]
    respm <- as.matrix(respdf[ , items, drop=FALSE])
  }
  
  # collect the arguments and change data argument
  all_pts$respdf <- respm
  names(all_pts)[grep("respdf",names(all_pts))] <- "respm"
  args_4pl  <- setdiff(names(formals(PP_4pl)), all_pts)
  args_4pl <- all_pts[names(all_pts)%in%args_4pl]
  args_gpcm <- setdiff(names(formals(PP_gpcm)), all_pts)
  args_gpcm <- all_pts[names(all_pts)%in%args_gpcm]
  args_all <- setdiff(names(formals(PPall)), all_pts)
  args_all <- all_pts[names(all_pts)%in%args_all]
  args_pfit <- setdiff(names(formals(Pfit)), all_pts)
  args_pfit <- all_pts[args_pfit]
  args_pfit$fitindices <- fitindices
  # check if first element is character
  if(is.character(respm[1,1])) stop("At least one response is of type character!\n")
  
  if(mod %in% c("1PL","2PL","3PL","4PL"))
  {
    pp_est <- do.call(PP_4pl,args_4pl)
  } else if(mod %in% c("PCM", "GPCM"))
  {
    pp_est <- PP_gpcm(respm, ...)
    
  } else if("MIXED"){ # mixed
    pp_est <- PPall(respm, ...)
  }
  
  ########### CALCULATE PERSON FIT ###############################
  # first we have to store the estimated person parameter
  args_pfit[[2]] <- pp_est
  names(args_pfit)[2] <- "pp"
  # this is needed to make shure, that no 'NA' arguments are in the list
  args_pfit <- args_pfit[names(args_pfit)%in%names(formals(Pfit))] 
  fit_calc <- do.call(Pfit,args_pfit)
  # cbind all pers fits together
  fit_calc <- do.call(cbind,fit_calc)
  ########### PUT IT ALL TOGETHER ############################### 
  #out <- list("personparameter"=pp_est,"personfit"=fit_calc)
  out <- data.frame(pp_est[["resPP"]][["resPP"]], fit_calc)
  return(out)
  
}





# ================= for eRm input =================================


#' @param RMobj A fitted Rasch Model (\code{RM()}) object which stems from the \code{eRm} package.
#' 
#' @rdname PPass
#' 
#' @export
#' 
#' @method PPass Rm
#' 
#' 
PPass.Rm <- function(RMobj, fitindices= c("lz","lzstar","infit","outfit"), ...)
{
  
  # catch additional arguments
  all_pts <- list(RMobj,...)
  fitindices <- match.arg(fitindices, several.ok = TRUE)  
  
  # ---------------------------------------------------
  # collect the arguments and change data argument
  all_pts$RMobj <- RMobj$X
  names(all_pts)[grep("RMobj",names(all_pts))] <- "respm"
  # --------------
  args_4pl  <- setdiff(names(formals(PP_4pl)), all_pts)
  args_4pl <- all_pts[names(all_pts)%in%args_4pl]
  # --------------
  args_gpcm <- setdiff(names(formals(PP_gpcm)), all_pts)
  args_gpcm <- all_pts[names(all_pts)%in%args_gpcm]
  # --------------
  args_pfit <- setdiff(names(formals(Pfit)), all_pts)
  args_pfit <- all_pts[args_pfit]
  args_pfit$fitindices <- fitindices
  # ---------------------------------------------------
  
  # geht leider nicht anders weil sowohl PCM als auch RM die Klassen Rm als auch eRm haben.
  if(RMobj$model == "RM")
  {
    args_4pl <- append( args_4pl,list("thres"=RMobj$betapar * (-1)) )
   
    # pp_est <- PP_4pl(respm=RMobj$X, thres=RMobj$betapar * (-1), ...)
    pp_est <- do.call(PP_4pl,args_4pl)
    
    
  } else if(RMobj$model == "PCM")
  {
    
    # get threshold parameters:
    
    # create threshold matrix:
    #tps <- eRm::thresholds(RMobj)$threshpar
    
    #len1        <- apply(RMobj$X, 2, function(x) length(unique(x))-1)
    #itemss      <- unlist(lapply(1:length(len1), function(x) rep(x, each=len1[x])))
    #wohin_zeile <- unlist(lapply(len1, function(x) 1:x))
    
    #names(wohin_zeile) <- NULL
    
    #thres <- matrix(NA,ncol=ncol(RMobj$X), nrow=max(len1))
    
    #for(i in 1:length(tps))
    #{
    #  thres[wohin_zeile[i], itemss[i]] <- tps[i]
    #}
    tps <- eRm::thresholds(RMobj)$threshtable[[1]]
    tps <- t(tps)
    tps[1,] <- 0
    thres <- tps
    # thres <- rbind(0,thres)
    slopes <- rep(1,ncol(thres))
    
    ########### ESTIMATE PERSON PARAMETERS ###############################  
    args_gpcm <- append( args_gpcm,list("thres"=thres, slopes=slopes) )
    # pp_est <- PP_gpcm(respm=RMobj$X, thres=thres, slopes=slopes, ...)
    pp_est <- do.call(PP_gpcm,args_gpcm)
    
  } else {
    
    stop("I don't know this model!")
    
  }
  
  
  ########### CALCULATE PERSON FIT ###############################  
  # first we have to store the estimated person parameter
  args_pfit[[2]] <- pp_est
  names(args_pfit)[2] <- "pp"
  
  # this is needed to make shure, that no 'NA' arguments are in the list
  args_pfit <- args_pfit[names(args_pfit)%in%names(formals(Pfit))] 
  fit_calc <- do.call(Pfit,args_pfit)
  # fit_calc <- Pfit(respm=RMobj$X,pp=pp_est,fitindices=fitindices,...)
  # cbind all pers fits together
  fit_calc <- do.call(cbind,fit_calc)
  ########### PUT IT ALL TOGETHER ############################### 
  #out <- list("personparameter"=pp_est,"personfit"=fit_calc)
  
  # if(class(pp_est)[1] == "gpcm")
  #   {
  #   out <- data.frame(pp_est[["resPP"]][["resPP"]], fit_calc)
  #   } else {
  #          out <- data.frame(pp_est[["resPP"]][["resPP"]], fit_calc)  
  #          }
  out <- data.frame(pp_est[["resPP"]][["resPP"]], fit_calc) 
  return(out)
}



