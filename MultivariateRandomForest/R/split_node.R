
split_node <- function(X,Y,m_feature,Index,i,model,min_leaf,Inv_Cov_Y,Command){
  ii=NULL
  Index_left=NULL
  Index_right=NULL
  if(length(Index)>min_leaf){ #create problem with 2
    ff2 = ncol(X) # number of features
    ff =sort(sample(ff2, m_feature)) #randomly taken 10 features, for each splits vary
    Result = splitt2(X,Y,m_feature,Index,Inv_Cov_Y,Command, ff)
    Index_left=Result[[1]]
    Index_right=Result[[2]]
    if(i==1){
      Result[[5]]=c(2,3)
    }else{
      j=1
      while (length(model[[j]])!=0){
        j=j+1
      }
      Result[[5]]=c(j,j+1)
    }

    model[[i]]=Result
    k=i
    i=1 #maybe unnecessary
    while (length(model[[i]])!=0){
      i=i+1
    }
    model[[Result[[5]][1]]]=Result[[1]]
    model[[Result[[5]][2]]]=Result[[2]]

    model=split_node(X,Y,m_feature,Index_left,model[[k]][[5]][1],model,min_leaf,Inv_Cov_Y,Command)

    model=split_node(X,Y,m_feature,Index_right,model[[k]][[5]][2],model,min_leaf,Inv_Cov_Y,Command)


  }else{
    ii[[1]]=matrix(Y[Index,],ncol=ncol(Y))
    model[[i]]=ii
  }


  return(model)
}
