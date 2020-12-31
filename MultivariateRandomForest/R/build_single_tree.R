
build_single_tree <- function(X, Y, m_feature, min_leaf,Inv_Cov_Y,Command){
  NN=round(nrow(X)/min_leaf)*100
  model=rep( list(NULL), NN )
  i=1
  Index=1:nrow(X)

  model=split_node(X,Y,m_feature,Index,i,model,min_leaf,Inv_Cov_Y,Command)
  return(model)
}
