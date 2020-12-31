plot.Variable_Select <- function(Variable_Select,K,plot.it = FALSE,N_R = 2,N_C = 2){
  H <- ncol(Variable_Select)
  M <- nrow(Variable_Select)
  List_Proportion_Variable_Select <- lapply(1:H,function(h){
    prop.var.full <- matrix(unlist(lapply(1:K,function(i){
      unlist(lapply(1:M,function(j){
        sum(Variable_Select[1:j,h] == i,na.rm = TRUE)/j
      }))
    })),nrow = M,byrow = FALSE)
    prop.var.full <- cbind(prop.var.full[,1:4],rowMeans(prop.var.full[,5:K]) )
  })
  if(plot.it){
    oldpar <- par("mfrow", "mar")
    on.exit(par(oldpar))  
    par(mfrow=c(N_R,N_C))
    for(h in 1:H){
      plot(range(1:M),range(List_Proportion_Variable_Select[[h]]),type = "n",xlab = "m",ylab = "Proportion of Selected Variables")
      K <- ncol(List_Proportion_Variable_Select[[h]])
      for(i in 1:K){
        if(i <= 4){
          lines(1:M,List_Proportion_Variable_Select[[h]][,i],type = "l",lwd = 3,col = i)  
        }else
        {
          lines(1:M,List_Proportion_Variable_Select[[h]][,i],type = "l",lwd = 3,col = "gray")
        }
      }      
    }
  }
  List_Proportion_Variable_Select
}