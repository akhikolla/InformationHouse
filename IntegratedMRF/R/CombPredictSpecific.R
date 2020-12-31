CombPredictSpecific <- function(finalX,finalY_train,Cell,finalY_train_cell,finalY_test_cell,n_tree,m_feature,min_leaf,Coeff){
  Serial=NULL
  for (p in length(Cell):1){
    nk=combn(1:length(Cell),p)
    sk=length(Serial)
    for (q in 1:dim(nk)[2]){
      Serial[[sk+q]]=nk[ ,q]
    }
  }
  ##
  Common_cell_dataset=NULL
  for (q in 1:length(Cell)){
    Common_cell_dataset[[q]]=intersect(finalY_train_cell,Cell[[q]])
  }
  Variable_number=ncol(finalY_train)
  if (ncol(finalY_train)==1){
    Command=1
  }else if (ncol(finalY_train)>1){
    Command=2
  }

  final_dataset=NULL
  finalY_dataset=NULL
  for (q in 1:length(Cell)){
    Cell_ind=match(Common_cell_dataset[[q]],Cell[[q]])
    final_dataset[[q]]=finalX[[q]][Cell_ind, ]
    final_dataset[[q]]=matrix(as.numeric(final_dataset[[q]]),nrow = dim(final_dataset[[q]])[1], ncol = dim(final_dataset[[q]])[2])

    Cell_ind_Y=match(Common_cell_dataset[[q]],finalY_train_cell)
    finalY_dataset[[q]]=matrix(finalY_train[Cell_ind_Y,],ncol=Variable_number)
    if (class(min_leaf)=="character" || min_leaf%%1!=0 || min_leaf<1 || min_leaf>nrow(finalY_dataset[[q]])) stop('Minimum leaf number can not be fractional or negative integer or string or greater than number of samples')
  }
  if (class(n_tree)=="character" || n_tree%%1!=0 || n_tree<1) stop('Number of trees in the forest can not be fractional or negative integer or string')
  if (class(m_feature)=="character" || m_feature%%1!=0 || m_feature<1) stop('Number of randomly selected features considered for a split can not be fractional or negative integer or string')

  Y_hat_test=NULL#matrix(rep(0,length(Cell)*length(finalY_test_cell)),nrow=length(finalY_test_cell))
  final_test=NULL
  Common_cell_test=NULL
  for (q in 1:length(Cell)){
    Common_cell_test=intersect(finalY_test_cell,Cell[[q]])
    Cell_ind=match(Common_cell_test,Cell[[q]])
    final_test[[q]]=finalX[[q]][Cell_ind, ]
    final_test[[q]]=matrix(as.numeric(final_test[[q]]),nrow = dim(final_test[[q]])[1], ncol = dim(final_test[[q]])[2])

    Y_hat_test[[q]] = build_forest_predict(final_dataset[[q]],finalY_dataset[[q]], n_tree, m_feature, min_leaf, final_test[[q]])
  }

  Coeff2=rep( list(NULL), length(Serial) )
  for (q in 1:length(Serial)){
    Temp2=Coeff[Serial[q][[1]]]
    Coeff2[[q]]=Temp2/sum(Temp2)
  }

  Common_cell=NULL
  Final_test=NULL
  Drug_sensitivity_cell_test2=finalY_test_cell
  for (q in 1:length(Serial)){
    D=NULL
    D=Drug_sensitivity_cell_test2
    for (check in 1:length(Serial[[q]])){
      D=intersect(D,Cell[Serial[[q]][check]][[1]])
    }
    Common_cell[[q]]=D
    Drug_sensitivity_cell_test2=setdiff(Drug_sensitivity_cell_test2,D)

    F_test=matrix(rep(0,length(Common_cell[[q]])*Variable_number),ncol=Variable_number)
    match_ind=NULL
    W=Serial[[q]]
    for (RR in 1:Variable_number){
      for (R in 1:length(W)){
        match_ind[[R]]=match(Common_cell[[q]],Common_cell_test[[W[R]]])
        F_test[,RR]=F_test[,RR]+Coeff2[[q]][R]*matrix(Y_hat_test[[W[R]]][match_ind[[R]],RR],ncol=1)
      }
    }

    Final_test[[q]]=F_test
  }
  Final_result=NULL
  for (q in 1:length(Serial)){
    if (length(Common_cell[[q]])>0){
      Final_result=rbind(Final_result,cbind(matrix(Common_cell[[q]],ncol=1),matrix(Final_test[[q]],ncol=Variable_number)))
    }
  }
  match_ind=match(finalY_test_cell,Final_result[,1])
  Final_prediction=Final_result[match_ind,]
  return(Final_prediction)
}
