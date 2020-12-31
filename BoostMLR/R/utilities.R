#-----------------------------------------------------------------------------------
# Find the index of the data based on increasing value of time for a given subject
Order_Time <- function(ID,Time){
  if(is.numeric(ID)){
    unq_id <- sort(unique(ID))
  }else
  {
    unq_id <- unique(ID)
  }
  n <- length(unq_id)
  ID_index <- unlist(lapply(1:n,function(i){
    Temp_index <- which(ID == unq_id[i])
    Temp_time <- Time[Temp_index]
    Ord_time <- order(Temp_time)
    Temp_index <- Temp_index[Ord_time]
    Temp_index
  }))
  return(ID_index)
}



#--------------------------------------------------
# Plot beta coefficient on the original scale

plotBeta <- function(Object,
                     N_Row = 1,
                     N_Column = 1,
                     plot.it = TRUE,
                     path_saveplot = NULL,
                     Verbose = TRUE){
  Raw_Variables <- which(Object$UseRaw == TRUE)
  K <- length(Raw_Variables)
  if(K == 0){
    stop("Variables not provided for plotting")
  }
  x_Names_plot <- Object$x_Names[Raw_Variables]
  L <- Object$Grow_Object$Dimensions$L
  tm <- Object$tm
  Ord_tm <- order(tm)
  if(plot.it){
    if(is.null(path_saveplot)){
      path_saveplot <- tempdir()
    }
    pdf(file = paste(path_saveplot,"/","Beta_Estimate.pdf",sep=""),width = 14,height = 14)
    for(l in 1:L){
      count <- 0
      oldpar <- par("mfrow", "mar")
      on.exit(par(oldpar)) 
      par(mfrow = c(N_Row,N_Column))
      for(k in Raw_Variables){
        count <- count + 1
        plot(tm[Ord_tm], Object$Tm_Beta[[ l ]][,k,drop = TRUE][Ord_tm],type = "l",xlab = "Time",ylab = "Beta Estimate",main = x_Names_plot[count],lwd = 3)   
      }
    }
    dev.off()
    if(Verbose){
      cat("Plot will be saved at:",path_saveplot,sep = "")
    }
  }
}

#--------------------------------------------------
# Plot mu and y on the original scale
plotMu <- function(Object,smooth_plot = TRUE,overall_mean = FALSE,y_range = NULL,plot.it = TRUE,path_saveplot = NULL,Verbose = TRUE){
  mu         <- Object$mu
  n          <- Object$Grow_Object$Dimensions$n
  id_index   <- Object$Grow_Object$Index$id_index
  y          <- Object$y
  tm         <- Object$tm
  L          <- Object$Grow_Object$Dimensions$L
  y_Names    <- Object$y_Names
  
  if(plot.it){
    if(is.null(path_saveplot)){
      path_saveplot <- tempdir()
    }
    pdf(file = paste(path_saveplot,"/","PredictedMu.pdf",sep=""),width = 14,height = 14)
    for(l in 1:L){
      if( is.null(y_range) ){
        y_range    <- range( c(mu[,l,drop = TRUE], y[,l,drop = TRUE])  )
      }
      plot(tm, y[,l,drop = TRUE],type = "n",ylim = y_range,xlab = "Time",ylab = "Predicted Y",main = y_Names[l])
      if(!overall_mean){
        for(i in 1:n){
          tm_subject <- tm[ id_index[[i]] ]
          Ord_tm     <- order(tm_subject)
          mu_subject <- mu[ id_index[[i]], l, drop = TRUE ]
          y_subject <- y[ id_index[[i]], l, drop = TRUE ]
          lines(tm_subject[Ord_tm], y_subject[Ord_tm],type = "p",lwd = 1,col = "gray")
          if(smooth_plot){
            fit.lo <- tryCatch({loess(mu_subject ~ tm_subject,span = 0.75)},
                               error = function(ex){NULL})
            if(!is.null(fit.lo)){
              predict.lo <- tryCatch({predict(fit.lo,newdata = tm_subject)},
                                     error = function(ex){NULL})
            }
            if( !is.null(predict.lo) ){
              lines(tm_subject[Ord_tm],predict.lo[Ord_tm],lwd = 1,col = 2,type = "l")
            }else
            {
              lines(tm_subject[Ord_tm], mu_subject[Ord_tm],type = "l",lwd = 1,col = 2)
            }
          }else
          {
            lines(tm_subject[Ord_tm], mu_subject[Ord_tm],type = "l",lwd = 1,col = 2)  
          }
        }
      }else
      {
        tm_unique <- sort(unique(tm))[unique(as.integer(quantile(1: length(unique(tm)) ,probs = seq(0,0.9,length.out = 15))))]
        tm.q <- cut(tm,breaks = tm_unique,labels = tm_unique[-1],include.lowest = TRUE)
        y_Mean <- tapply(y[,l,drop = TRUE],tm.q,mean,na.rm = TRUE)
        mu_Mean <- tapply(mu[,l,drop = TRUE],tm.q,mean,na.rm = TRUE)
        tm_Mean <- as.numeric(levels(sort(unique(tm.q))))
        lo.mu <- lowess(tm_Mean,mu_Mean)$y
        lines(tm_Mean,y_Mean,type = "p",pch = 19)
        lines(tm_Mean,lo.mu,type = "l",lwd = 3,col = 1)
      }
    }
    dev.off()
    if(Verbose){
      cat("Plot will be saved at:",path_saveplot,sep = "")
    }
  }
}


is.hidden.dt_Add <-  function (user.option) {
  if (is.null(user.option$dt_Add)) {
    NULL
  }
  else {
    user.option$dt_Add
  }
}


is.hidden.predict.dt_Add <-  function (user.option) {
  if (is.null(user.option$dt_Add)) {
    NULL
  }
  else {
    user.option$dt_Add
  }
}

is.hidden.rho <-  function (user.option) {
  if (is.null(user.option$rho)) {
    NULL
  }
  else {
    user.option$rho
  }
}


is.hidden.phi <-  function (user.option) {
  if (is.null(user.option$phi)) {
    NULL
  }
  else {
    user.option$phi
  }
}


is.hidden.prob_min <-  function (user.option) {
  if (is.null(user.option$prob_min)) {
    0.1
  }
  else {
    user.option$prob_min
  }
}


is.hidden.prob_max <-  function (user.option) {
  if (is.null(user.option$prob_max)) {
    0.9
  }
  else {
    user.option$prob_max
  }
}


is.hidden.importance_Coef <-  function (user.option) {
  if (is.null(user.option$importance_Coef)) {
    FALSE
  }
  else {
    TRUE
  }
}
