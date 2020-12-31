Combination <- function(finalX,finalY_train,Cell,finalY_train_cell,n_tree,m_feature,min_leaf,Confidence_Level){

  Serial=NULL
  for (p in length(Cell):1){
    nk=combn(1:length(Cell),p)
    sk=length(Serial)
    for (q in 1:dim(nk)[2]){
      Serial[[sk+q]]=nk[ ,q]
    }
  }
  ##
  Common_cell_train=finalY_train_cell
  for (q in 1:length(Cell)){
    Common_cell_train=intersect(Common_cell_train,Cell[[q]])
  }

  final=NULL
  MM_feature=rep(0,length(Cell))
  for (q in 1:length(Cell)){
    Cell_ind=match(Common_cell_train,Cell[[q]])
    final[[q]]=finalX[[q]][Cell_ind, ]
    MM_feature[q]=ncol(final[[q]])
  }
  finalY=NULL
  ia6=match(Common_cell_train,finalY_train_cell)
  finalY=matrix(finalY_train[ia6,],ncol=ncol(finalY_train))
  Variable_number=ncol(finalY)
  if (Variable_number>1){
    Command=2
  }else if(Variable_number==1){
    Command=1
  }
  if (class(n_tree)=="character" || n_tree%%1!=0 || n_tree<1) stop('Number of trees in the forest can not be fractional or negative integer or string')
  if (class(m_feature)=="character" || m_feature%%1!=0 || m_feature<1 || m_feature>min(MM_feature)) stop('Number of randomly selected features considered for a split can not be fractional or negative integer or string or greater than minimum number of features')
  if (class(min_leaf)=="character" || min_leaf%%1!=0 || min_leaf<1 || min_leaf>nrow(finalY)) stop('Minimum leaf number can not be fractional or negative integer or string or greater than number of samples')
  if (class(Confidence_Level)=="character" || Confidence_Level>100 || Confidence_Level<1) stop('Confidence Interval can not be negative integer or string or greater than 100')

  ################################## BSP ###############################
  ptm1=proc.time()
  if (nrow(finalY)<50){
    N=floor(0.75*nrow(finalY))
  }else if (nrow(finalY)>=50 && nrow(finalY)<101){
    N=floor(nrow(finalY)/2)
  }else if (nrow(finalY)>=101){
    N=floor(nrow(finalY)/3)
  }

  Y_hat_BSP=NULL
  for (q in 1:length(Cell)){
    Y_hat_BSP[[q]]=matrix(rep(0,length(finalY)),ncol=ncol(finalY))
  }
  bootsam_FF=NULL
  Index=NULL
  Index=1:nrow(finalY)
  theta <- function(x){x}
  results <- bootstrap::bootstrap(Index,N,theta) #no indics, gives number
  bootsam=results$thetastar

  Store=rep( list(NULL), length(Cell) )
  for (q in 1:length(Cell)){
    Store[[q]]=rep( list(NULL), nrow(finalY) )
  }
  Store_Jack=rep( list(NULL), Variable_number )
  for (RR in 1:Variable_number){
    Store_Jack[[RR]]=rep( list(NULL), length(Serial))
    for (q in 1:length(Serial)){
      Store_Jack[[RR]][[q]]=rep( list(NULL), nrow(finalY) )
    }
  }
  BSP_error_alll_mae=rep(list(NULL), Variable_number)
  BSP_error_alll_mse=rep(list(NULL), Variable_number)
  BSP_error_alll_corr=rep(list(NULL), Variable_number)
  for (S in 1:Variable_number){
    BSP_error_alll_mae[[S]]=matrix(rep(0,length(Serial)*N),ncol=N)
    BSP_error_alll_mse[[S]]=matrix(rep(0,length(Serial)*N),ncol=N)
    BSP_error_alll_corr[[S]]=matrix(rep(0,length(Serial)*N),ncol=N)
  }
  for (FF in 1:N){
    bootsam_FF=bootsam[,FF]
    Index_FF=unique(bootsam_FF)
    Index_pred=setdiff(Index, Index_FF)

    finalY_bsp=matrix(finalY[bootsam_FF,],ncol=Variable_number)
    finalY_bsp_pred=matrix(finalY[Index_pred,],ncol=Variable_number)

    final_genome=NULL
    finalX_bsp=NULL
    finalX_bsp_pred=NULL
    finalY_pred=NULL

    for (q in 1:length(Cell)){
      finalX_bsp[[q]]=final[[q]][bootsam_FF,]
      finalX_bsp_pred[[q]]=final[[q]][Index_pred,]
      finalY_pred[[q]]=build_forest_predict(finalX_bsp[[q]], finalY_bsp, n_tree, m_feature, min_leaf, finalX_bsp_pred[[q]])
      for (R in 1:length(Index_pred)){
        Store[[q]][[Index_pred[R]]]=rbind(Store[[q]][[Index_pred[R]]],finalY_pred[[q]][R,])
      }
    }

    for (RR in 1:Variable_number){
      for (S in 1:length(Serial)){
        final_genome_BSP1=NULL
        W=Serial[[S]]
        for (q in 1:length(Serial[[S]])){
          final_genome_BSP1=cbind(final_genome_BSP1,matrix(finalY_pred[[W[q]]][,RR],ncol=1))
        }
        BSP_temp_coeff1=matrix(limSolve::lsei(A=final_genome_BSP1, B=matrix(finalY_bsp_pred[,RR],ncol=1), E=rep(1,dim(final_genome_BSP1)[2]), F=1)$X, ncol=1)
        BSP_error_alll_mae[[RR]][S,FF]=mean(abs(final_genome_BSP1%*%BSP_temp_coeff1-finalY_bsp_pred[,RR]))
        BSP_error_alll_mse[[RR]][S,FF]=mean((final_genome_BSP1%*%BSP_temp_coeff1-finalY_bsp_pred[,RR])^2)
        BSP_error_alll_corr[[RR]][S,FF]=stats::cor(final_genome_BSP1%*%BSP_temp_coeff1,finalY_bsp_pred[,RR])
        Err_bsp_temp=abs(final_genome_BSP1%*%BSP_temp_coeff1-finalY_bsp_pred[,RR])
        for (R in 1:length(Index_pred)){
          Store_Jack[[RR]][[S]][[Index_pred[R]]]=rbind(Store_Jack[[RR]][[S]][[Index_pred[R]]],Err_bsp_temp[R,])
        }
      }
    }
  }
  Error_BSP_Jack=rep( list(NULL), Variable_number )
  for (RR in 1:Variable_number){
    Error_BSP_Jack[[RR]]=matrix(rep(0,length(Serial)*nrow(finalY)),ncol=nrow(finalY))
    for (S in 1:length(Serial)){
      for (FFF in 1:nrow(finalY)){
        Error_BSP_Jack[[RR]][S,FFF]=mean(Store_Jack[[RR]][[S]][[FFF]])
      }
    }
  }
  BSP_error_all_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  BSP_error_all_mse=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  BSP_error_all_corr=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  BSP_errors=NULL
  for (RR in 1:Variable_number){
    for (S in 1:length(Serial)){
      BSP_error_all_mae[S,RR]=mean(BSP_error_alll_mae[[RR]][S,])
      BSP_error_all_mse[S,RR]=mean(BSP_error_alll_mse[[RR]][S,])
      BSP_error_all_corr[S,RR]=mean(BSP_error_alll_corr[[RR]][S,])
    }
    BSP_errors=cbind( BSP_errors,rbind(BSP_error_all_mae[1,RR],BSP_error_all_mse[1,RR],BSP_error_all_corr[1,RR]))
  }


  final_genome_BSP=NULL
  for (RR in 1:Variable_number){
    final_genome_BSP[[RR]]=matrix(rep(0,length(Cell)*nrow(finalY)),nrow=nrow(finalY))
  }
  for (RR in 1:Variable_number){
    for (q in 1:length(Cell)){
      for (FF in 1:nrow(finalY)){
        final_genome_BSP[[RR]][FF,q]=mean(Store[[q]][[FF]][,RR])
      }
    }
  }
  final_BSP=matrix(rep(0,length(finalY)),ncol=Variable_number)
  for (RR in 1:Variable_number){
    Store2=rep( list(NULL), nrow(finalY) )
    for (FF in 1:N){
      bootsam_FF=bootsam[,FF]
      Index_FF=unique(bootsam_FF)
      Index_pred=setdiff(Index, Index_FF)

      finalY_bsp1=matrix(finalY[bootsam_FF,RR],ncol=1)
      finalY_bsp1_pred=matrix(finalY[Index_pred,RR],ncol=1)

      BSP_temp_coeff=matrix(limSolve::lsei(A=final_genome_BSP[[RR]][bootsam_FF,], B=finalY_bsp1, E=rep(1,dim(final_genome_BSP[[RR]])[2]), F=1)$X, ncol=1)
      final_BSP_index=final_genome_BSP[[RR]][Index_pred,]%*%BSP_temp_coeff

      for (R in 1:length(Index_pred)){
        Store2[[Index_pred[R]]]=c(Store2[[Index_pred[R]]],final_BSP_index[R,])
      }
    }
    for (FF in 1:nrow(finalY)){
      final_BSP[FF,RR]=mean(Store2[[FF]])
    }
  }
  ## Find Combination Weight
  BSP_coeff=rep(list(NULL), length(Serial))
  for (S in 1:length(Serial)){
    BSP_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }

  for (RR in 1:Variable_number){
    for (q in 1:length(Cell)){
      Y_hat_BSP[[q]]=final_genome_BSP[[RR]][,q]
    }
    for (S in 1:length(Serial)){
      final_genome_BSP1=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_BSP1=cbind(final_genome_BSP1,matrix(Y_hat_BSP[[W[q]]],ncol=1))
      }
      BSP_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_BSP1, B=finalY[,RR], E=rep(1,dim(final_genome_BSP1)[2]), F=1)$X, ncol=1)
    }
  }
  ptm2=proc.time()-ptm1
  message("Elapsed Time for Bootstrap Error Estimation is ", ptm2[[3]])
  ####################### Resubstitution Error######################
  ptm1 <- proc.time()
  finalY_pred_resub=NULL
  final_genome_resub=rep(list(NULL), Variable_number)
  for (q in 1:length(Cell)){
    finalY_pred_resub[[q]]=build_forest_predict(final[[q]], finalY, n_tree, m_feature, min_leaf, final[[q]])
    for (RR in 1:Variable_number){
      final_genome_resub[[RR]]=cbind(final_genome_resub[[RR]],finalY_pred_resub[[q]][,RR])
    }
  }
  Resub_errors=NULL
  Resub_error_mae=NULL
  Resub_error_mse=NULL
  Resub_corr=NULL
  final_resub=matrix(rep(0,length(finalY)),ncol=Variable_number)
  for (RR in 1:Variable_number){
    error=error_calculation(final_genome_resub[[RR]],matrix(finalY[,RR],ncol=1))
    Resub_error_mae[RR]=error[[2]]
    Resub_error_mse[RR]=error[[3]]
    Resub_corr[RR]=error[[4]]
    final_resub[,RR]=error[[1]]
  }
  Resub_errors=rbind(Resub_error_mae,Resub_error_mse,Resub_corr)

  Resub_coeff=rep(list(NULL), length(Serial))
  Y_hat_resub=NULL
  for (S in 1:length(Serial)){
    Resub_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }
  Error_resub_Jack=rep( list(NULL), Variable_number )
  Resub_error_all_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  Resub_error_all_mse=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  Resub_error_all_corr=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  for (RR in 1:Variable_number){
    Error_resub_Jack[[RR]]=matrix(rep(0,length(Serial)*nrow(finalY)),ncol=nrow(finalY))
    for (q in 1:length(Cell)){
      Y_hat_resub[[q]]=final_genome_resub[[RR]][,q]
    }
    for (S in 1:length(Serial)){
      final_genome_Resub1=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_Resub1=cbind(final_genome_Resub1,matrix(Y_hat_resub[[W[q]]],ncol=1))
      }
      Resub_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_Resub1, B=finalY[,RR], E=rep(1,dim(final_genome_Resub1)[2]), F=1)$X, ncol=1)
      Resub_error_all_mae[S,RR]=mean(abs(final_genome_Resub1%*%Resub_coeff[[S]][,RR]-finalY[,RR]))
      Resub_error_all_mse[S,RR]=mean((final_genome_Resub1%*%Resub_coeff[[S]][,RR]-finalY[,RR])^2)
      Resub_error_all_corr[S,RR]=cor(final_genome_Resub1%*%Resub_coeff[[S]][,RR],finalY[,RR])
      Error_resub_Jack[[RR]][S,]=abs(final_genome_Resub1%*%Resub_coeff[[S]][,RR]-finalY[,RR])
    }
  }
  ptm2=proc.time()-ptm1
  message("Elapsed Time for Resubstitution Error Estimation is ", ptm2[[3]])
  ####################### 0.632+BSP Error######################
  ptm1 <- proc.time()
  BSP632_coeff=rep(list(NULL), length(Serial))
  Y_hat_BSP632=NULL
  for (S in 1:length(Serial)){
    BSP632_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }
  BSP632_error_all_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  BSP632_error_all_mse=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  BSP632_error_all_corr=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  R_hat=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  W_hat=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  Gamma=rep(0,Variable_number)
  for (RR in 1:Variable_number){
    Gamma_all=NULL
    for (YY in 1:nrow(finalY)){
      Gamma_all=c(Gamma_all,abs(finalY[YY,RR]-finalY[,RR]))
    }
    Gamma[RR]=mean(Gamma_all)
    for (S in 1:length(Serial)){
      R_hat[S,RR]=(BSP_error_all_mae[S,RR]-Resub_error_all_mae[S,RR])/(Gamma[RR]-Resub_error_all_mae[S,RR])
      W_hat[S,RR]=0.632/(1-0.368*R_hat[S,RR])
      BSP632_error_all_mae[S,RR]=W_hat[S,RR]*BSP_error_all_mae[S,RR]+(1-W_hat[S,RR])*Resub_error_all_mae[S,RR]
      BSP632_error_all_mse[S,RR]=W_hat[S,RR]*BSP_error_all_mse[S,RR]+(1-W_hat[S,RR])*Resub_error_all_mse[S,RR]
      BSP632_error_all_corr[S,RR]=W_hat[S,RR]*BSP_error_all_corr[S,RR]+(1-W_hat[S,RR])*Resub_error_all_corr[S,RR]
      final_genome_BSP6321=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_BSP6321=cbind(final_genome_BSP6321,matrix(W_hat[S,RR]*final_genome_BSP[[RR]][,W[q]]+(1-W_hat[S,RR])*final_genome_resub[[RR]][,W[q]],ncol=1))
      }
      BSP632_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_BSP6321, B=finalY[,RR], E=rep(1,dim(final_genome_BSP6321)[2]), F=1)$X, ncol=1)
    }
  }
  BSP632_errors=rbind(BSP632_error_all_mae[1,],BSP632_error_all_mse[1,],BSP632_error_all_corr[1,])
  ptm2=proc.time()-ptm1
  message("Elapsed Time for 0.632+ Bootstrap Error Estimation is ", ptm2[[3]])
  ################################ 0.632BSP error Jackkniffe Confidence Interval ####################
  ptm1 <- proc.time()
  Low_confidence=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  High_confidence=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  mean_D=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  for (RR in 1:Variable_number){
    for (S in 1:length(Serial)){
      err=0.632*Error_BSP_Jack[[RR]][S,]+0.368*Error_resub_Jack[[RR]][S,]
      #err=abs((0.632*final_BSP[,RR]+0.368*final_resub[,RR])-finalY[,RR])
      n=length(err)
      s=sqrt((1/(n*(n-1)))*sum((err-mean(err))^2))
      alpha=1-(Confidence_Level/100)
      z_alpha=stats::qnorm((1+2*alpha/2-1)/2) # z_alpha=sqrt(2)*erfinv(2*alpha/2-1)
      z_1alpha=stats::qnorm((1+2*(1-alpha/2)-1)/2) # z_1alpha=sqrt(2)*erfinv(2*(1-alpha/2)-1)
      mean_D[S,RR]=mean(err)

      Low_confidence[S,RR]=max(0,(mean(err)-(s*abs(z_alpha))))# Z_a/2 error function or, quantile function
      High_confidence[S,RR]=mean(err)+(s*(z_1alpha))
    }
  }
  ptm2=proc.time()-ptm1
  message("Elapsed Time for Jackknife Confidence Interval Calculation is ", ptm2[[3]])
  ################################ Leave One out error #######################################
  ptm1 <- proc.time()
  Y_hat_LOO=NULL
  for (q in 1:length(Cell)){
    Y_hat_LOO[[q]]=matrix(rep(0,length(finalY)),ncol=Variable_number)
  }

  for (q in 1:length(Cell)){
    for (FF in 1:nrow(finalY)){
      Index2=1:nrow(finalY)
      Index2=setdiff(Index2,FF)
      Y1=matrix(finalY[Index2,],ncol=Variable_number)
      ##
      X1=final[[q]][Index2,]
      Xt=matrix(final[[q]][FF,],nrow=1)
      Y_hat_LOO[[q]][FF,] = build_forest_predict(X1,Y1, n_tree, m_feature, min_leaf, Xt)
    }
  }
  LOO_errors=NULL
  LOO_error_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  LOO_error_mse=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  LOO_corr=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)

  LOO_coeff=rep(list(NULL), length(Serial))
  for (S in 1:length(Serial)){
    LOO_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }
  for (RR in 1:Variable_number){
    for (S in 1:length(Serial)){
      final_genome_LOO=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_LOO=cbind(final_genome_LOO,matrix(Y_hat_LOO[[W[q]]][,RR],ncol=1))
      }
      LOO_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_LOO, B=finalY[,RR], E=rep(1,dim(final_genome_LOO)[2]), F=1)$X, ncol=1)
      LOO_error_mae[S,RR]=mean(abs(final_genome_LOO%*%LOO_coeff[[S]][,RR]-finalY[,RR]))
      LOO_error_mse[S,RR]=mean((final_genome_LOO%*%LOO_coeff[[S]][,RR]-finalY[,RR])^2)
      LOO_corr[S,RR]=stats::cor(final_genome_LOO%*%LOO_coeff[[S]][,RR],finalY[,RR])
    }
  }
  LOO_errors=rbind(LOO_error_mae[1,],LOO_error_mse[1,],LOO_corr[1,])
  ptm2=proc.time()-ptm1
  message("Elapsed Time for Leave-one-out Error Estimation is ", ptm2[[3]])

  ################################ N_fold cross validation error #######################################
  ptm1 <- proc.time()
  if (nrow(finalY)<50){
    N_Fold=10
  }else if (nrow(finalY)>=50 && nrow(finalY)<101){
    N_Fold=10
  }else if (nrow(finalY)>=101){
    N_Fold=5
  }
  Y_hat_Nfold=rep(list(NULL),length(Cell))
  Result_CV=CrossValidation(final[[1]],finalY,N_Fold)
  FoldedIndex=Result_CV[[5]]
  Folded=NULL
  for (Fold in 1:N_Fold){
    Folded=c(Folded,FoldedIndex[[Fold]])
  }
  tmp=sort(Folded,index.return=TRUE)
  s=tmp$x; idx=tmp$ix
  for (q in 1:length(Cell)){
    for (FF in 1:N_Fold){
      Index2=1:nrow(finalY)
      Index2=setdiff(Index2,FoldedIndex[[FF]])
      Y1=matrix(finalY[Index2,],ncol=Variable_number)
      ##
      X1=final[[q]][Index2,]
      Xt=final[[q]][FoldedIndex[[FF]],]
      Y_hat_Nfold[[q]] = rbind(Y_hat_Nfold[[q]],build_forest_predict(X1,Y1, n_tree, m_feature, min_leaf, Xt))
    }
    Y_hat_Nfold[[q]]=matrix(Y_hat_Nfold[[q]][idx,],ncol=Variable_number)
  }
  Nfold_errors=NULL
  Nfold_error_mae=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  Nfold_error_mse=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)
  Nfold_corr=matrix(rep(0,length(Serial)*Variable_number),ncol=Variable_number)

  Nfold_coeff=rep(list(NULL), length(Serial))
  for (S in 1:length(Serial)){
    Nfold_coeff[[S]]=matrix(rep(0,length(Serial[[S]])*Variable_number),ncol=Variable_number)
  }
  for (RR in 1:Variable_number){
    for (S in 1:length(Serial)){
      final_genome_Nfold=NULL
      W=Serial[[S]]
      for (q in 1:length(Serial[[S]])){
        final_genome_Nfold=cbind(final_genome_Nfold,matrix(Y_hat_Nfold[[W[q]]][,RR],ncol=1))
      }
      Nfold_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_Nfold, B=finalY[,RR], E=rep(1,dim(final_genome_Nfold)[2]), F=1)$X, ncol=1)
      Nfold_error_mae[S,RR]=mean(abs(final_genome_Nfold%*%Nfold_coeff[[S]][,RR]-finalY[,RR]))
      Nfold_error_mse[S,RR]=mean((final_genome_Nfold%*%Nfold_coeff[[S]][,RR]-finalY[,RR])^2)
      Nfold_corr[S,RR]=cor(final_genome_Nfold%*%Nfold_coeff[[S]][,RR],finalY[,RR])
    }
  }
  Nfold_errors=rbind(Nfold_error_mae[1,],Nfold_error_mse[1,],Nfold_corr[1,])
  ptm2=proc.time()-ptm1
  message("Elapsed Time for N-fold cross validation Error Estimation is ", ptm2[[3]])

  Result=NULL
  Result[[1]]=BSP_coeff
  Result[[2]]=Nfold_coeff
  Result[[3]]=BSP632_coeff
  Result[[4]]=LOO_coeff
  Result[[5]]=rbind(BSP_errors,Nfold_errors,BSP632_errors,LOO_errors)
  Result[[6]]=rep( list(NULL), 2 )
  Result[[6]][[1]]=Low_confidence
  Result[[6]][[2]]=High_confidence
  Result[[7]]=BSP_error_all_mae
  Result[[8]]=Nfold_error_mae
  Result[[9]]=BSP632_error_all_mae
  Result[[10]]=LOO_error_mae
  #### Figures ######
  Error_Matrix1=rbind(Result[[7]],Result[[8]],Result[[9]],Result[[10]])
  Error_Matrix2=matrix(c(rep("BSP_error",length(Serial)),rep("N_fold_error",length(Serial)),rep("0.632+ BSP_error",length(Serial)),rep("LOO_error",length(Serial))),ncol=1)

  SS=matrix(rep(0,length(Serial)*length(Cell)),ncol=length(Cell))
  for (i in 1:length(Serial)){
    for (j in 1:length(Cell)){
      if (length(which(j==Serial[[i]]))==1){
        SS[i,j]=1
      }
    }
  }
  SSC=transform(paste0(SS[,1],SS[,2]))
  if (length(Cell)>2) {
    for (k in 3:length(Cell)){
      SSC=transform(paste0(SSC[,1],SS[,k]))
    }
  }
  Error_type=NULL
  Error_Matrix3=cbind(rep(matrix(SSC[,1],ncol=1),4))
  Error_Matrix=cbind(Error_Matrix2,Error_Matrix3,Error_Matrix1)
  colnames(Error_Matrix) <- c("Error_type","Combination",rep("Response",ncol(finalY_train)))
  for (i in 1:ncol(finalY_train)){
    Error_Matrix[,i+2]=matrix(format(round(as.numeric(Error_Matrix[,i+2]), 4), nsmall = 4),ncol=1)
  }
  #  library(ggplot2)
  plot_list = list()
  for (i in 1:ncol(finalY_train)){
    qq=ggplot2::ggplot(data.frame(Error_Matrix), ggplot2::aes(x=Combination,y=Error_Matrix[,(i+2)], fill=Error_type))+ ggplot2::geom_bar( stat="identity",position="dodge")+ ggplot2::labs(title = "Different Error Estimation for Integrated Models")+ ggplot2::xlab("Combination Index of Integrated Models")+ggplot2::ylab("Mean Absolute Error")+ ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust = 0))
    plot_list[[i]] = qq
  }
  for (i in 1:ncol(finalY_train)){
    tiff(paste("Error_Estimation_Response_",i,".tiff",sep=""), width = 7, height = 10, units = 'in', res = 700)
    print(plot_list[[i]])
    dev.off()
  }

  return(Result)
}
