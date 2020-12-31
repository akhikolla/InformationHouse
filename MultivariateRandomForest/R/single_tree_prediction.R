
single_tree_prediction <- function(Single_Model,X_test,Variable_number){


  Y_pred=matrix(  0*(1:nrow(X_test)*Variable_number)  ,nrow=nrow(X_test),  ncol=Variable_number)

  for (k in 1:nrow(X_test)){
    xt=X_test[k, ]
    i=1
    Result_temp=predicting(Single_Model,i,xt,Variable_number)
    Y_pred[k,]=unlist(Result_temp)

  }
  #Y_pred1=unlist(Y_pred, recursive = TRUE)
  #Y_pred1=matrix(Y_pred1,nrow=nrow(X_test))
  return(Y_pred)
}
