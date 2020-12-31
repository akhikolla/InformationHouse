vimp.BoostMLR <- function(Object,xvar.names = NULL,joint = FALSE){
  
   if(sum(inherits(Object, c("BoostMLR", "predict"), TRUE) == c(1, 2)) != 2) {
    stop("This function only works for objects of class `(BoostMLR, predict)'")
  }

x <- Object$Data$Org_x
P <- ncol(x)
x_Names <- Object$x_Names
y_Names <- Object$y_Names

if(is.null(xvar.names)){
  vimp_set <- (1:P) - 1
  x_Names <- x_Names
  if(joint){
    x_Names <- "joint_vimp"
  }
}else
{
  n.x.names <- length(xvar.names)
  vimp_set <- (match(xvar.names,x_Names)) - 1
  if(any(is.na( vimp_set ))){
    stop("xvar.names do not match with variable names from original data")
  }
  
  if(joint){
    x_Names <- "joint_vimp"
  }else
  {
    x_Names <- xvar.names
  }
}

if(joint) {
  p <- 1
}
else {
  p <- length(vimp_set)  
}

Org_x <- Object$Data$Org_x
Org_y <- Object$Data$Org_y
id <- Object$Data$id
tm <- Object$Data$tm
x_Mean <- Object$Data$x_Mean
x_Std_Error <- Object$Data$x_Std_Error
y_Mean <- Object$Data$y_Mean
y_Std_Error <- Object$Data$y_Std_Error

n <- Object$Pred_Object$Dimensions$n
ni <- Object$Pred_Object$Dimensions$ni
N <- Object$Pred_Object$Dimensions$N
L <- Object$Pred_Object$Dimensions$L
K <- Object$Pred_Object$Dimensions$K
Dk <- Object$Pred_Object$Dimensions$Dk
H <- Object$Pred_Object$Dimensions$H
n_unq_tm <- Object$Pred_Object$Dimensions$n_unq_tm
  
unq_id <- Object$Pred_Object$Index$unq_id
unq_tm <- Object$Pred_Object$Index$unq_tm
unq_x <- Object$Pred_Object$Index$unq_x
id_index <- Object$Pred_Object$Index$id_index
tm_index <- Object$Pred_Object$Index$tm_index
x_index <- Object$Pred_Object$Index$x_index
unq_x_New <- Object$Pred_Object$Index$unq_x_New
Index_Bt <- Object$Pred_Object$Index$Index_Bt

Bt <- Object$Pred_Object$BS$Bt
Bxt <- Object$Pred_Object$BS$Bxt
Bx_K <- Object$Pred_Object$BS$Bx_K
Bt_H <- Object$Pred_Object$BS$Bt_H
Bx <- Object$Pred_Object$BS$Bx

UseRaw <- Object$Pred_Object$UseRaw
Beta_Hat_List <- Object$Pred_Object$Beta_Hat_List
Mopt <- Object$Mopt
rmse <- Object$rmse
nu <- Object$nu

Vec_zero <- Object$Pred_Object$Vec_zero
mu_zero_vec <- Object$Pred_Object$mu_zero_vec

Time_Varying <- Object$Pred_Object$Time_Varying

obj_C <- vimp_BoostMLR_C(Org_x,
                         Org_y,
                         tm,
                         id,
                         x_Mean,
                         x_Std_Error,
                         y_Mean,
                         y_Std_Error,
                         n,
                         ni,
                         N,
                         L,
                         K,
                         p,
                         H,
                         Dk,
                         n_unq_tm,
                         UseRaw,
                         id_index,
                         tm_index,
                         unq_x_New,
                         Index_Bt,
                         vimp_set,
                         joint,
                         Bt,
                         Bt_H,
                         Bx,
                         Bxt,
                         Bx_K,
                         Beta_Hat_List,
                         Mopt,
                         nu,
                         rmse,
                         Time_Varying,
                         Vec_zero,
                         mu_zero_vec)

vimp <- obj_C$vimp
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

obj <- vimp
invisible(obj)

}
