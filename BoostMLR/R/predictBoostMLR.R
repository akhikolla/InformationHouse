predictBoostMLR    <-  function(Object,
                                 x,
                                 tm,
                                 id,
                                 y,
                                 M,
                                 importance = FALSE,
                                 eps = 1e-5,
                                 ...)
{
  
  user.option <- list(...)
  dt_Add <- is.hidden.predict.dt_Add(user.option)
  
  importance_Coef <- is.hidden.importance_Coef(user.option)
  
  if(missing(tm) && missing(x)){
    stop("tm and x both missing")
  }
  
  if(!missing(tm) && missing(id)){
    stop("id is missing")
  }
  
  CrossSectional <- FALSE
  if(missing(tm) && !missing(x) ){
    if(!missing(id)){
      if(!(length(sort(unique(id))) == nrow(x)) ){
        stop("tm is missing")
      }
    }else
    {
      id <- 1:nrow(x)
    }
    tm <- rep(0, length(x))
    CrossSectional <- TRUE
  }
  
  if(missing(x) && !missing(tm)){
    x_miss <- TRUE
    All_RawX <- TRUE
    x <- cbind(rep(1,length(tm)))
    if(length(sort(unique(id))) == nrow(x) ){
      CrossSectional <- TRUE
    }
  }else
  {
    x_miss <- FALSE
  }
  
  Time_Varying <- Object$Time_Varying
  
  if(!missing(tm) && Time_Varying == FALSE && CrossSectional == FALSE){
    if(x_miss){
      x <- tm
    }else
    {
      x <- x #cbind(x,tm)
    }
  }
  
  if (any(is.na(id))) {
    stop("missing values encountered in id: remove observations with missing values")
  }
  
  if (!missing(y)) {
    if ( any(is.na(y)) ) {
      #stop("missing values encountered in y: remove observations with missing values")
    }
    testFlag <- TRUE
    y_Names <- colnames(y)
  }
  else{
    testFlag <- FALSE
    L <- Object$Grow_Object$Dimensions$L
    y <- matrix(0, nrow =  nrow(x),ncol = L)
    y_Names <- paste("y",1:L,sep="")
  }
  
  if(!is.matrix(y)){
    y <- data.matrix(y)
  }
  
  Time_Unmatch <- Object$Grow_Object$Time_Unmatch
  N <- nrow(x)
  
  if(!is.null(dt_Add)){
    if(!is.list(dt_Add)){
      stop("dt_Add must be a list")
    }
    K_Add <- length(dt_Add)
    nullObj <- lapply(1:K_Add,function(kk){
      nc_K_Add <- ncol(dt_Add[[kk]])
      if(nc_K_Add != 3){
        stop("Each element of dt_Add must be a dataset with 3 columns arrange in order of id, time, x")
      }
      NULL
    })
    
    Ord_id_tm <- Order_Time(ID = id,Time = tm)
    id  <- id[Ord_id_tm]
    tm  <- tm[Ord_id_tm]
    x   <- x[Ord_id_tm,,drop = FALSE]
    y   <- y[Ord_id_tm,,drop = FALSE]
    x_Add_New <- matrix(NA,nrow = N,ncol = K_Add)
    x_Names_Add <- rep(NA,K_Add)
    Time_Add_New <- matrix(NA,nrow = N,ncol = K_Add)
    Time_Names_Add <- rep(NA,K_Add)
    for(kk in 1:K_Add){
      Ord_id_tm_Add <- Order_Time(ID = dt_Add[[kk]][,1],Time =  dt_Add[[kk]][,2])
      dt_Add[[kk]] <- dt_Add[[kk]][Ord_id_tm_Add,,drop = FALSE]
      id_Add <- dt_Add[[kk]][,1]
      x_Names_Add[kk] <- names(dt_Add[[kk]][,3,drop = FALSE])
      Time_Names_Add[kk] <- names(dt_Add[[kk]][,2,drop = FALSE])
      if(any(is.na(id_Add))){
        stop("Missing values observed for id in dt_Add")
      }
      unq_id_Add <- unique(id_Add)
      n_Add <- length(unq_id_Add)
      nullObj <- unlist(lapply(1:n_Add,function(i){
        Which_id <- which(unq_id_Add[i] == id)
        ni <- length(Which_id)
        if(ni > 0){
          Which_id_Add <- which(id_Add == unq_id_Add[i])
          ni_Add <- length(Which_id_Add)
          tm_Add <- dt_Add[[kk]][Which_id_Add,2]
          x_Add <- dt_Add[[kk]][Which_id_Add,3]
          for(j in 1:ni){
            for(jj in 1:ni_Add){
              if((!is.na(tm_Add[jj]) && !is.na(tm[Which_id[j]]))){
                if(tm_Add[jj] <=  tm[Which_id[j]]){
                  x_Add_New[Which_id[j], kk] <<- x_Add[jj]
                  Time_Add_New[Which_id[j], kk] <<- tm_Add[jj]
                }
              }
            }
          }
        }
        NULL  
      }))
    }
    colnames(x_Add_New) <- x_Names_Add
    x <- cbind(x,x_Add_New)
    colnames(Time_Add_New) <- Time_Names_Add
  } else
  {
    Time_Add_New <- matrix(0,nrow = N,ncol = 1)
    colnames(Time_Add_New) <- "Time_Add"
  }
  
  if(!is.matrix(x)){
    x <- data.matrix(x)
  }
  
  x_Names <- colnames(x)
  
  K <- ncol(x)
  
  if(is.null(x_Names)){
    x_Names <- paste("x",1:K,sep="")
  }
  
  if(!identical(x_Names , Object$x_Names) ){
    stop("Covariate from grow and predict function are not matching")
  }
  
  if(missing(M)){
    M <- Object$Grow_Object$Regulate$M  
  }
  
  L <- ncol(y)
  
  if(is.null(y_Names)){
    y_Names <- paste("y",1:L,sep="")
  }
  
  H  <- Object$Grow_Object$Dimensions$H
  Dk <- Object$Grow_Object$Dimensions$Dk
  
  x_Mean      <- Object$Grow_Object$Data$x_Mean
  x_Std_Error <- Object$Grow_Object$Data$x_Std_Error
  y_Mean      <- Object$Grow_Object$Data$y_Mean
  y_Std_Error <- Object$Grow_Object$Data$y_Std_Error
  
  unq_tm   <- Object$Grow_Object$Index$unq_tm
  unq_x    <- Object$Grow_Object$Index$unq_x
  
  Bt   <- Object$Grow_Object$BS$Bt
  Bx   <- Object$Grow_Object$BS$Bx
  
  nu  <- Object$Grow_Object$Regulate$nu
  
  Beta             <- Object$Grow_Object$Beta_Estimate$Beta
  Beta_Hat_List    <- Object$Grow_Object$Beta_Estimate$Beta_Hat_List
  
  UseRaw <- Object$UseRaw
  vimpFlag <- (importance == TRUE && testFlag == TRUE)
  vimpFlag_Coef <- (importance_Coef == TRUE && testFlag == TRUE)
  
  obj_C <- predict_BoostMLR_C(x,
                              tm,
                              id,
                              y,
                              x_Mean,
                              x_Std_Error,
                              y_Mean,
                              y_Std_Error,
                              K,
                              L,
                              H,
                              Dk,
                              unq_tm,
                              unq_x,
                              Bt,
                              Bx,
                              UseRaw,
                              Time_Add_New,
                              Time_Unmatch,
                              Beta,
                              Beta_Hat_List,
                              testFlag,
                              M,
                              nu,
                              Time_Varying,
                              vimpFlag,
                              vimpFlag_Coef,
                              eps)
  
  Error_Rate <- obj_C$Error_Rate
  colnames(Error_Rate) <- y_Names
  vimp <- obj_C$vimp
  vimp_Coef <- obj_C$vimp_Coef

  if(vimpFlag){
    names(vimp) <- y_Names
    for(l in 1:L){
      rownames(vimp[[l]]) <- x_Names
      if(H == 1){
        vimp[[l]] <- vimp[[l]][,1,drop = FALSE]
      }
      if(H == 1){
        colnames(vimp[[l]]) <- "Main_Eff"
      } else {
        colnames(vimp[[l]]) <- c("Main_Eff",paste("Int_Eff.",1:H,sep=""))
      }
    }
  }
  
  if(vimpFlag_Coef){
    names(vimp_Coef) <- y_Names
    for(l in 1:L){
      rownames(vimp_Coef[[l]]) <- x_Names
      if(H == 1){
        vimp_Coef[[l]] <- vimp_Coef[[l]][,1,drop = FALSE]
      }
      if(H == 1){
        colnames(vimp_Coef[[l]]) <- "Main_Eff"
      } else {
        colnames(vimp_Coef[[l]]) <- c("Main_Eff",paste("Int_Eff.",1:H,sep=""))
      }
    }
  }
  
  mu <- obj_C$Org_mu
  colnames(mu) <- y_Names
  
  if(testFlag)
  {
    mu_Mopt <- obj_C$Org_mu_Mopt
    colnames(mu_Mopt) <- y_Names
  } else
  {
    mu_Mopt <- NA
  }
  
  Pred_Object <- obj_C$Pred_Object
  Pred_Object$Dimensions = obj_C$Dimensions
  Pred_Object$Index = obj_C$Index
  Pred_Object$BS = obj_C$BS
  Pred_Object$UseRaw = UseRaw
  Pred_Object$Time_Varying = Time_Varying
  Pred_Object$Beta_Hat_List = Beta_Hat_List
  
  obj <- list(Data = obj_C$Data,
              x_Names = x_Names,
              y_Names = y_Names,
              mu = mu,
              mu_Mopt = mu_Mopt,
              Error_Rate = Error_Rate,
              Mopt = obj_C$Mopt,
              nu = nu,
              rmse = obj_C$rmse,
              vimp = vimp,
              vimp_Coef = vimp_Coef,
              Pred_Object = Pred_Object)
  class(obj) <- c("BoostMLR", "predict")
  invisible(obj)
}