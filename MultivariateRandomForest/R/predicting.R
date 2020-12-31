
predicting <- function(Single_Model,i,X_test,Variable_number){

  Result=NULL

  if(length(Single_Model[[i]])==5){
    feature_no=Single_Model[[i]][[3]]
    feature_value=X_test[feature_no]
    if(feature_value<Single_Model[[i]][[4]]){  #feature value less than threshold value
      #i=i*2+1
      Result=predicting(Single_Model,Single_Model[[i]][[5]][1],X_test,Variable_number)
    }else{                                    #feature value greater than threshold value
      #i=i*2+2
      Result=predicting(Single_Model,Single_Model[[i]][[5]][2],X_test,Variable_number)
    }
  }else{
    Result=matrix(  0*Variable_number,  ncol=Variable_number)
    if (Variable_number>1){
      for (jj in 1:Variable_number){
        Result[,jj]=mean(Single_Model[[i]][[1]][,jj])
      }
    }else {
      for (jj in 1:Variable_number){
        Result[,jj]=mean(unlist(Single_Model[[i]][[1]]))
      }
    }

  }
  return(Result)
}
