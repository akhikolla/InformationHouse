BoostMLR     <- function(x,
                         tm,
                         id,
                         y,
                         Time_Varying = TRUE,
                         BS_Time = TRUE,
                         nknots_t = 10,
                         d_t = 3,
                         All_RawX = TRUE,
                         RawX_Names,
                         nknots_x = 7,
                         d_x = 3,
                         M = 200,
                         nu = 0.05,
                         Mod_Grad = TRUE,
                         Shrink = FALSE,
                         VarFlag = TRUE,
                         lower_perc = 0.25,
                         upper_perc = 0.75,
                         NLambda = 100,
                         Verbose = TRUE,
                         ...)
{
  if(missing(y)){
    stop("y is missing")
  }

  if(is.data.frame(y)){
    y <- data.matrix(y)
  }
  
  if(!is.vector(y) && !is.matrix(y) ){
    stop("y must be a vector or matrix")
  } else 
    {
    if(is.vector(y) && is.atomic(y)) {
      y <- matrix(y,ncol = 1)  
    }
  }
  
  if(missing(tm) && missing(x)){
    stop("tm and x both missing")
  }

  if(!missing(tm) && missing(id)){
    stop("id is missing")
  }

  CrossSectional <- FALSE
  if(missing(tm) && !missing(x) ){
    if(!missing(id)){
      if(!(length(sort(unique(id))) == nrow(y)) ){
        stop("tm is missing")
      }
    } else
    {
      id <- 1:nrow(y)
    }
    tm <- rep(0, nrow(y))
    nknots_t <-  1
    d_t <- 0
    CrossSectional <- TRUE
    VarFlag <- FALSE
  }
  
  if(missing(x) && !missing(tm)){
    x_miss <- TRUE
    All_RawX <- TRUE
    x <- cbind(rep(1,nrow(y)))
    if(length(sort(unique(id))) == nrow(y) ){
      CrossSectional <- TRUE
      VarFlag <- FALSE
    }
  }else
  {
    x_miss <- FALSE
  }

  if(nknots_t <= 0 && BS_Time == TRUE){
    stop("nknots_t must be > 0")
  }
  
  if(All_RawX == FALSE && any(nknots_x <= 0) ){
    stop("nknots_x must be > 0")
  }
  
  if(missing(tm) && Time_Varying == TRUE){
    Time_Varying <- FALSE
  }
  
  if(!BS_Time){
    nknots_t <-  1
    d_t <- 0
  }
  
  #if(!missing(tm) && Time_Varying == FALSE && CrossSectional == FALSE){
  #  if(x_miss){
  #    x <- tm
  #  }
  #  nknots_t <-  1
  #  d_t <- 0
  #}
  
  if ( any(is.na(id))  ) {  
    stop("missing values encountered in id: remove observations with missing values")
  }
  
  Time_Unmatch <- rep(FALSE,ncol(x))
  N <- nrow(y)
  
 user.option <- list(...)
 Lambda_Scale <- 1
 rho <- is.hidden.rho(user.option)
 phi <- is.hidden.phi(user.option)
 Lambda_Ridge <- 0
 Ridge_Penalty <- FALSE
  
 dt_Add <- is.hidden.dt_Add(user.option) 
  
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
    Time_Unmatch <- c(Time_Unmatch,rep(TRUE, ncol(x_Add_New)))
    colnames(Time_Add_New) <- Time_Names_Add
  } else
  {
    Time_Add_New <- matrix(0,nrow = N,ncol = 1)
    colnames(Time_Add_New) <- "Time_Add"
  }
  
  if(!is.matrix(x)){
    x <- data.matrix(x)
  }
  
  K <- ncol(x)
  L <- ncol(y)

  x_Names <- colnames(x)
  y_Names <- colnames(y)

  if(is.null(x_Names)){
    x_Names <- paste("x",1:K,sep="")
  }
  
  if(is.null(y_Names)){
    y_Names <- paste("y",1:L,sep="")
  }

  if( ( is.null(rho) && !is.null(phi) ) || ( !is.null(rho) && is.null(phi) )  ){
    stop("rho or phi is null")
  }
  
  if( !is.null(rho) && !is.null(phi) ){
    VarFlag <- FALSE
    if(! ( (length(rho) == 1 || length(rho) == L) && (length(phi) == 1 || length(phi) == L) ) ){
      stop("Length of rho and phi must be 1 or L")
    }
    
    if( any(rho < -1) || any(rho > 1) || any(phi <= 0) ){
      stop("rho should be between -1 and 1 and phi should be > 0")
    }
    
    if(length(rho) == 1 && length(phi) == 1){
      rho <- rep(rho,L)
      phi <- rep(phi,L)
    }
    Rho <- matrix(rep(rho,M),nrow = M,byrow = TRUE)
    Phi <- matrix(rep(phi,M),nrow = M,byrow = TRUE)
  }
  else {
    rho <- rep(0,L)
    phi <- rep(1,L)
    Rho <- matrix(rep(rho,M),nrow = M,byrow = TRUE)
    Phi <- matrix(rep(phi,M),nrow = M,byrow = TRUE)
  }
  
  unq_x <- lapply(1:K,function(k){
    sort_unique_C_NA(x[,k,drop = TRUE])
  })
  
  n_unq_x <- unlist(lapply(1:K,function(k){
    length(unq_x[[k]])
  }))
  
  if( !(length(nknots_x) == 1 || length(nknots_x) == K) ){
    stop( paste("Length of nknots_x should be either 1 or",K,sep = " ") )
  }
  
  if( !(length(d_x) == 1 || length(d_x) == K) ){
    stop( paste("Length of d_x should be either 1 or",K,sep = " ") )
  }
  
  if(length(nknots_x) == 1){
    nknots_x <- rep(nknots_x,K)
  }
  
  if(length(d_x) == 1){
    d_x <- rep(d_x,K)
  }
  
  Unique_Limit <- unlist(lapply(1:K,function(k){
    nknots_x[k] + d_x[k]
  }))
  
  Categorical_x <- x_Names[which(unlist(lapply(1:K,function(k){
    (n_unq_x[k] < Unique_Limit[k])
  })))]
  
  if(length(Categorical_x) > 0){ 
    if(missing(RawX_Names)){
      RawX_Names <- Categorical_x
    }else
    {
      RawX_Names <- union(RawX_Names,Categorical_x)
    }
  }
  
  if(All_RawX){
    UseRaw <- rep(TRUE,K)
  }else
  {
    if(missing(RawX_Names) ){
      UseRaw <- rep(FALSE,K)
    }else
    {
      UseRaw <- (!is.na(match(x_Names,RawX_Names)))
      if(all(UseRaw == FALSE)){
        stop("RawX.names do not match with the variables names from x")
      }
    }
  }

  ProcessedData <- DataProcessing_C(x,
                                    y,
                                    id, 
                                    tm,
                                    x_miss)
  
  Org_x       <- ProcessedData$Data$Org_x
  Org_y       <- ProcessedData$Data$Org_y
  id          <- ProcessedData$Data$id
  tm          <- ProcessedData$Data$tm
  x           <- ProcessedData$Data$x
  y           <- ProcessedData$Data$y
  
  if(!is.matrix(Org_x)){
    Org_x <- data.matrix(Org_x)
  }
  if(!is.matrix(Org_y)){
    Org_y <- data.matrix(Org_y)
  }
  if(!is.matrix(x)){
    x <- data.matrix(x)
  }
  if(!is.matrix(y)){
    y <- data.matrix(y)
  }
  
  x_Mean      <- ProcessedData$Data$x_Mean
  x_Std_Error <- ProcessedData$Data$x_Std_Error
  y_Mean      <- ProcessedData$Data$y_Mean
  y_Std_Error <- ProcessedData$Data$y_Std_Error
  
  unq_id   <- ProcessedData$Index$unq_id
  id_index <- ProcessedData$Index$id_index
  
  n  <- ProcessedData$Dimensions$n
  K  <- ProcessedData$Dimensions$K
  L  <- ProcessedData$Dimensions$L
  ni <- ProcessedData$Dimensions$ni
  N  <- ProcessedData$Dimensions$N

  if(all(ni == 1)){
    VarFlag <- FALSE
  }
  
  H <- nknots_t + d_t
  unq_tm <- sort_unique_C_NA(tm)
  n_unq_tm <- length(unq_tm)
  if( (n_unq_tm < H) && (BS_Time == TRUE) ){
    H <- n_unq_tm
  }

  if(CrossSectional == TRUE || Time_Varying == FALSE){
    Bt <- cbind(rep(1,n_unq_tm));
  }else
  {
    if(BS_Time){
      Bt <- bs(x = unq_tm,df = H,degree = d_t,intercept = TRUE)
    }else
    {
      Bt <-  cbind(unq_tm)
    }
  }

  unq_x <- lapply(1:K,function(k){
    if(UseRaw[k]){
      NA
    }else
    {
      sort_unique_C_NA(x[,k,drop = TRUE])
    }
  })

  n_unq_x <- unlist(lapply(1:K,function(k){
    if(UseRaw[k]){
      NA
    }else
    {
      length(unq_x[[k]])
    }
  }))
  
  Dk <- unlist(lapply(1:K,function(k){
    if(UseRaw[k]){
     out <- 1
    }else
    {
      Dk_Temp <- nknots_x[k] + d_x[k]
      if(n_unq_x[k] < Dk_Temp){
        out <- n_unq_x[k]
      }else
      {
        out <- Dk_Temp  
      }
    }
   out  
  }))
  
  count <- 0
  Bx <- lapply(1:K,function(k){
    if(UseRaw[k]){
      if(Time_Unmatch[k]){
        count <<- count + 1
        x[,k,drop = FALSE]*Time_Add_New[ , count, drop = FALSE]
      }
      else {
        x[,k,drop = FALSE]
      }
    }else
    {
        bs(x = unq_x[[k]],df = Dk[k],degree = d_x[k],intercept = TRUE)
    }
  })

  Bx_Scale <- lapply(1:K,function(k){
    if(UseRaw[k]){
      1
    }else
    {
      rep(1,Dk[k])
    }
  })
  
  Lambda_Ridge_Vec <- unlist(lapply(1:K,function(k){
    Lambda_Ridge
  }))
  
  if(M < 10){
    Verbose <- FALSE
  }

  obj_C  <- BoostMLR_C(Org_x,
                       Org_y,
                       id,
                       tm,
                       x,
                       y,
                       x_Mean,
                       x_Std_Error,
                       y_Mean,
                       y_Std_Error,
                       n,
                       K,
                       L,
                       H,
                       Dk,
                       ni,
                       N,
                       unq_id,
                       unq_tm,
                       unq_x,
                       id_index,
                       Bt,
                       Bx,
                       Bx_Scale,
                       Time_Add_New,
                       Time_Unmatch,
                       nu,
                       M,
                       Mod_Grad,
                       UseRaw,
                       Lambda_Ridge_Vec,
                       Ridge_Penalty,
                       Shrink,
                       lower_perc,
                       upper_perc,
                       Lambda_Scale,
                       NLambda,
                       VarFlag,
                       rho,
                       phi,
                       Verbose)
  
  
  Tm_Beta <- lapply(1:obj_C$Dimensions$L,function(l){
    Out <- matrix(unlist(lapply(1:obj_C$Dimensions$K,function(k){
      if(!UseRaw[k]){
        rep(NA, obj_C$Dimensions$N)
      }else
      {
        Reduce("+",lapply(1:obj_C$Dimensions$H,function(h){
          unlist(lapply(1:obj_C$Dimensions$n,function(i){
            obj_C$Beta_Estimate$Tm_Beta_C[[k]][[1]][[h]][[l]][[i]]
          }))
        }))
      }
    })),ncol = obj_C$Dimensions$K,byrow = FALSE)
    colnames(Out) <- x_Names
    Out
  })

  if(Time_Varying == FALSE){
    Tm_Beta <- lapply(1:obj_C$Dimensions$L,function(l){
      Tm_Beta[[l]][1,,drop = TRUE]
    })
  }
  
  names(Tm_Beta) <- y_Names
  
  Beta_Estimate <- obj_C$Beta_Estimate
  Beta_Estimate$Tm_Beta <- Tm_Beta
  
  Rho <- Phi <- matrix(NA,nrow = M,ncol = L)
  colnames(Phi) <- y_Names
  colnames(Rho) <- y_Names
  Error_Rate <- obj_C$Error_Rate
  colnames(Error_Rate) <- y_Names
  
  if(FALSE){
  if(VarFlag){
    NullObj <- lapply(1:L,function(l){
      lapply(1:M,function(m){
        
        Residual_Data <- data.frame(y = (obj_C$Data$Org_y[,l] - obj_C$mu_List[[m]][,l]) ,tm = obj_C$Data$tm, id = obj_C$Data$id, obj_C$Data$Org_x)
        gls.obj <- tryCatch({gls(y ~ ., data = Residual_Data,
                                 correlation = corCompSymm(form = ~ 1 | id))},
                            error = function(ex){NULL})
        if (is.null(gls.obj)) {
          gls.obj <- tryCatch({gls(y ~ 1, data = Residual_Data,
                                   correlation = corCompSymm(form = ~ 1 | id))},
                              error = function(ex){NULL})
        }
        
        if (!is.null(gls.obj)) {
          phi_Temp <- gls.obj$sigma^2
          Phi[m,l] <<- ifelse(phi_Temp == 0,1,phi_Temp)
          rho_Temp <- as.numeric(coef(gls.obj$modelStruct$corStruc, unconstrained = FALSE))
          Rho[m,l] <<- max(min(0.999, rho_Temp, na.rm = TRUE), -0.999)
          Result <- c(phi_Temp,rho_Temp)
        }
        NULL
      })
      NULL
    })
  }
  }  
  
  x       <- obj_C$Data$Org_x
  y       <- obj_C$Data$Org_y
  id      <- obj_C$Data$id
  tm      <- obj_C$Data$tm
  M       <- obj_C$Regulate$M
  nu      <- obj_C$Regulate$nu
  mu      <- obj_C$mu_List[[M]] 
  
  if(VarFlag){
    phi <- obj_C$Phi[M,]
    rho <- obj_C$Rho[M,]
  }
  
  Grow_Object <- list(Data = obj_C$Data,
                      Dimensions = obj_C$Dimensions,
                      Index = obj_C$Index,
                      BS = obj_C$BS,
                      Regulate = obj_C$Regulate,
                      Beta_Estimate = Beta_Estimate,
                      mu = obj_C$mu,
                      mu_List = obj_C$mu_List,
                      mu_zero = obj_C$mu_zero,
                      Vec_zero = obj_C$Vec_zero,
                      Mod_Grad = Mod_Grad,
                      phi = phi,
                      rho = rho,
                      Time_Unmatch = Time_Unmatch,
                      Time_Add_New = if(is.null(dt_Add)) NULL else Time_Add_New)

  
  Variable_Select = (obj_C$Variable_Select + 1)
  Response_Select = (obj_C$Response_Select + 1)
  
  Variable_Select[Variable_Select == 0] <- NA
  Response_Select[Response_Select == 0] <- NA
  
  Phi <- obj_C$Phi
  Rho <- obj_C$Rho
  colnames(Phi) <- colnames(Rho) <- y_Names
  
  obj <- list(x  = x,
              id = id,
              tm = tm,
              y  = y,
              UseRaw = UseRaw,
              x_Names = x_Names,
              y_Names = y_Names,
              M = M,
              nu = nu,
              Tm_Beta = Tm_Beta,
              mu = mu,
              Error_Rate = Error_Rate,
              Variable_Select = Variable_Select,
              Response_Select = Response_Select,
              VarFlag = VarFlag,
              Time_Varying = Time_Varying,
              Phi = Phi,
              Rho = Rho,
              Lambda_List = obj_C$Lambda_List,
              Grow_Object = Grow_Object)
  
  class(obj) <- c("BoostMLR")
  invisible(obj)
}
