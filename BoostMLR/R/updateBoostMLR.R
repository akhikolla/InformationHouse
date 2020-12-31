updateBoostMLR <- function(Object,M_Add,Verbose = TRUE,...){
  
  n  <- Object$Grow_Object$Dimensions$n
  K  <- Object$Grow_Object$Dimensions$K
  L  <- Object$Grow_Object$Dimensions$L
  H  <- Object$Grow_Object$Dimensions$H
  Dk <- Object$Grow_Object$Dimensions$Dk
  ni <- Object$Grow_Object$Dimensions$ni
  N  <- Object$Grow_Object$Dimensions$N
  
  Org_x       <- Object$Grow_Object$Data$Org_x
  Org_y       <- Object$Grow_Object$Data$Org_y
  id          <- Object$Grow_Object$Data$id
  tm          <- Object$Grow_Object$Data$tm
  x           <- Object$Grow_Object$Data$x
  y           <- Object$Grow_Object$Data$y
  x_Mean      <- Object$Grow_Object$Data$x_Mean
  x_Std_Error <- Object$Grow_Object$Data$x_Std_Error
  y_Mean      <- Object$Grow_Object$Data$y_Mean
  y_Std_Error <- Object$Grow_Object$Data$y_Std_Error
  
  unq_id   <- Object$Grow_Object$Index$unq_id
  unq_tm   <- Object$Grow_Object$Index$unq_tm
  unq_x    <- Object$Grow_Object$Index$unq_x
  id_index <- Object$Grow_Object$Index$id_index
  tm_index <- Object$Grow_Object$Index$tm_index
  x_index  <- Object$Grow_Object$Index$x_index
  
  Bt   <- Object$Grow_Object$BS$Bt
  Bx   <- Object$Grow_Object$BS$Bx
  Bt_H <- Object$Grow_Object$BS$Bt_H
  Bx_K <- Object$Grow_Object$BS$Bx_K
  Bxt  <- Object$Grow_Object$BS$Bxt
  Bx_Scale  <- Object$Grow_Object$BS$Bx_Scale
  
  nu               <- Object$Grow_Object$Regulate$nu
  M                <- Object$Grow_Object$Regulate$M
  Lambda_Ridge_Vec <- Object$Grow_Object$Regulate$Lambda_Ridge_Vec
  Shrink           <- Object$Grow_Object$Regulate$Shrink
  Ridge_Penalty    <- Object$Grow_Object$Regulate$Ridge_Penalty
  Lambda_Scale     <- Object$Grow_Object$Regulate$Lambda_Scale
  NLambda          <- Object$Grow_Object$Regulate$NLambda
  lower_perc       <- Object$Grow_Object$Regulate$lower_perc
  upper_perc       <- Object$Grow_Object$Regulate$upper_perc
  Mod_Grad         <- Object$Grow_Object$Mod_Grad
  
  Error_Rate        <- Object$Error_Rate
  Variable_Select   <- (Object$Variable_Select - 1)
  Response_Select   <- (Object$Response_Select - 1)
  mu_List           <- Object$Grow_Object$mu_List
  mu                <- Object$Grow_Object$mu
  mu_zero           <- Object$Grow_Object$mu_zero
  Vec_zero          <- Object$Grow_Object$Vec_zero
  UseRaw            <- Object$UseRaw
  VarFlag           <- Object$VarFlag
  Time_Varying      <- Object$Time_Varying 
  x_Names           <- Object$x_Names
  y_Names           <- Object$y_Names
  Lambda_List       <- Object$Lambda_List
  Phi               <- Object$Phi
  Rho               <- Object$Rho
  phi               <- Object$Grow_Object$phi
  rho               <- Object$Grow_Object$rho
  
  Beta                 <- Object$Grow_Object$Beta_Estimate$Beta  
  Beta_Hat_List        <- Object$Grow_Object$Beta_Estimate$Beta_Hat_List
  Beta_Hat_List_Iter   <- Object$Grow_Object$Beta_Estimate$Beta_Hat_List_Iter
  Sum_Beta_Hat_List    <- Object$Grow_Object$Beta_Estimate$Sum_Beta_Hat_List
  lower_Beta_Hat_Noise <- Object$Grow_Object$Beta_Estimate$lower_Beta_Hat_Noise
  upper_Beta_Hat_Noise <- Object$Grow_Object$Beta_Estimate$upper_Beta_Hat_Noise
  List_Trace_Bxt_gm    <- Object$Grow_Object$Beta_Estimate$List_Trace_Bxt_gm
  
  M_New <- M + M_Add
  
  if(M_Add < 10){
    Verbose <- FALSE
  }
  
  obj_C <- update_BoostMLR_C(Org_x,
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
                             tm_index,
                             x_index,
                             Bt,
                             Bx,
                             Bt_H,
                             Bx_K,
                             Bxt,
                             Bx_Scale,
                             nu,
                             M,
                             M_New,
                             UseRaw,
                             Shrink,
                             Ridge_Penalty,
                             Lambda_Ridge_Vec,
                             Lambda_Scale,
                             NLambda,
                             lower_perc,
                             upper_perc,
                             Lambda_List,
                             mu,
                             mu_List,
                             mu_zero,
                             Vec_zero,
                             Error_Rate,
                             Variable_Select,
                             Response_Select,
                             Beta_Hat_List,
                             Sum_Beta_Hat_List,
                             Beta,
                             Beta_Hat_List_Iter,
                             lower_Beta_Hat_Noise,
                             upper_Beta_Hat_Noise,
                             List_Trace_Bxt_gm,
                             Mod_Grad,
                             VarFlag,
                             phi,
                             rho,
                             Phi,
                             Rho,
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
  
  Rho <- Phi <- matrix(NA,nrow = M_New,ncol = L)
  colnames(Phi) <- y_Names
  colnames(Rho) <- y_Names
  Error_Rate <- obj_C$Error_Rate
  colnames(Error_Rate) <- y_Names
  
  if(FALSE){
  if(VarFlag){
    Rho[1:M, ] <- Object$Rho
    Phi[1:M, ] <- Object$Phi
    
    NullObj <- lapply(1:L,function(l){
      lapply( (M+1) : M_New,function(m){
        
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
                      Time_Unmatch = Object$Grow_Object$Time_Unmatch,
                      Time_Add_New = Object$Grow_Object$Time_Add_New)
  
  
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
  
  class(obj) <- c("BoostMLR", "update")
  invisible(obj)
}
