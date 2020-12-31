#' Integrated Prediction of Testing samples from integrated RF or MRF model
#'
#' Generates Random Forest or Multivariate Random Forest model for each subtype of dataset and predicts testing samples using the generated models.
#' Subsequently, the prediction for different subtypes of dataset are combined using the Combination weights generated from
#' Integrated Model which is based on Bootstrap error estimate
#'
#' @param finalX List of Matrices where each matrix represent a specific data subtype (such as genomic characterizations for
#' drug sensitivity prediction). Each subtype can have different types of features. For example, if there are three subtypes containing
#'  100, 200 and 250 features respectively,  finalX will be a list containing 3 matrices of sizes M x 100, M x 200 and M x 250
#'  where M is the number of Samples.
#' @param finalY_train A M x T matrix of output features for training samples, where M is number of samples and
#' T is the number of output features. The dataset is assumed to contain no missing values. If there are missing values, an imputation method
#' should be applied before using the function. A function 'Imputation' is included within the package.
#' @param Cell It contains a list of samples (the samples can be represented either numerically by indices or by names) for each data subtype.
#' For the example of 3 data subtypes, it will be a list containing 3 arrays where each array contains the sample information for each data subtype.
#' @param finalY_train_cell Cell lines of output features for training samples
#' @param finalY_test_cell Cell lines of output features for testing samples
#' @param n_tree number of trees in the forest, which must be positive integer
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node, which must be positive integer
#' @param min_leaf minimum number of samples in the leaf node, which must be positive integer and less than or equal to M (number of training samples)
#'
#' @return Final Prediction of testing samples based on provided testing sample names.
#' @details
#' Input matrix and output response of training samples have been used to build Random Forest or Multivariate Random Forest model for each subtype of
#' a dataset. These models are used to calculate prediction of
#' testing samples for each subtype separately. Subsequently Combination Weights are used to integrate the predictions from data subtypes.
#'
#' Combination Weight Generation: For M x N dataset, N number of bootstrap sampling sets are considered. For each bootstrap sampling set and each subtype, a Random Forest (RF)
#' or, Multivariate Random Forest (MRF) model is generated, which is used for calculating the prediction performance for out-of-bag samples.
#' The prediction performance for each dataset subtypes is based on the averaging over different bootstrap training sets.
#' The combination weights (regression coefficients) for each combination of subtypes are generated using least Square Regression from the
#' individual subtype predictions and used to integrate the predictions from data subtypes.
#'
#' The specific set of combination weights to be used for testing samples will depend on the number of data subtypes available
#' for the testing samples. Note that not all subtype information maybe available for all samples.
#' As an example with three data subtypes, a testing sample with all subtype data available will use
#' the combination weights corresponding to Serial [1 2 3] where as if subtype 3 is not available, the function will
#' use the combination weights corresponding to Serial [1 2].
#' @examples
#' library(IntegratedMRF)
#' data(Dream_Dataset)
#' Tree=1
#' Feature=1
#' Leaf=10
#' finalX=Dream_Dataset[[1]]
#' Cell=Dream_Dataset[[2]]
#' Y_train_Dream=Dream_Dataset[[3]]
#' Y_train_cell=Dream_Dataset[[4]]
#' Y_test=Dream_Dataset[[5]]
#' Y_test_cell=Dream_Dataset[[6]]
#' Drug=c(1,2,3)
#' Y_train_Drug=matrix(Y_train_Dream[,Drug],ncol=length(Drug))
#' IntegratedPrediction(finalX,Y_train_Drug,Cell,Y_train_cell,Y_test_cell,Tree,Feature,Leaf)
#'
#' @importFrom utils combn
#' @importFrom stats lsfit
#' @export

IntegratedPrediction <- function(finalX,finalY_train,Cell,finalY_train_cell,finalY_test_cell,n_tree,m_feature,min_leaf){
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
  for (q in 1:length(Cell)){
    Cell_ind=match(Common_cell_train,Cell[[q]])
    final[[q]]=finalX[[q]][Cell_ind, ]
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
  if (class(m_feature)=="character" || m_feature%%1!=0 || m_feature<1) stop('Number of randomly selected features considered for a split can not be fractional or negative integer or string')
  if (class(min_leaf)=="character" || min_leaf%%1!=0 || min_leaf<1 || min_leaf>nrow(finalY)) stop('Minimum leaf number can not be fractional or negative integer or string or greater than number of samples')

  ################################## BSP ###############################
  if (nrow(finalY)<50){
    N=floor(0.75*nrow(finalY))
  }else if (nrow(finalY)>=50 && nrow(finalY)<101){
    N=floor(nrow(finalY)/2)
  }else if (nrow(finalY)>=101){
    N=floor(nrow(finalY)/5)
  }

  Y_hat_BSP=NULL
  for (q in 1:length(Cell)){
    Y_hat_BSP[[q]]=matrix(rep(0,length(finalY)),ncol=ncol(finalY))
  }
  bootsam_FF=NULL
  Index=NULL
  Index=1:nrow(finalY)
  #   library(bootstrap)
  theta <- function(x){x}
  results <- bootstrap::bootstrap(Index,N,theta) #no indics, gives number
  bootsam=results$thetastar

  Store=rep( list(NULL), length(Cell) )
  for (q in 1:length(Cell)){
    Store[[q]]=rep( list(NULL), nrow(finalY) )
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
      #final_genome=cbind(final_genome,finalY_pred[[q]])
      for (R in 1:length(Index_pred)){
        Store[[q]][[Index_pred[R]]]=rbind(Store[[q]][[Index_pred[R]]],finalY_pred[[q]][R,])
      }
    }
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
#       BB=stats::lsfit(final_genome_BSP1, finalY[,RR], wt = NULL, intercept = FALSE, tolerance = 1e-07)
#       BSP_coeff[[S]][,RR]=unname(BB$coefficients)
      BSP_coeff[[S]][,RR]=matrix(limSolve::lsei(A=final_genome_BSP1, B=finalY[,RR], E=rep(1,dim(final_genome_BSP1)[2]), F=1)$X, ncol=1)
    }
  }
  ############################
  Common_cell_dataset=NULL
  for (q in 1:length(Cell)){
    Common_cell_dataset[[q]]=intersect(finalY_train_cell,Cell[[q]])
  }

  final_dataset=NULL
  finalY_dataset=NULL
  for (q in 1:length(Cell)){
    Cell_ind=match(Common_cell_dataset[[q]],Cell[[q]])
    final_dataset[[q]]=finalX[[q]][Cell_ind, ]
    final_dataset[[q]]=matrix(as.numeric(final_dataset[[q]]),nrow = dim(final_dataset[[q]])[1], ncol = dim(final_dataset[[q]])[2])

    Cell_ind_Y=match(Common_cell_dataset[[q]],finalY_train_cell)
    finalY_dataset[[q]]=matrix(finalY_train[Cell_ind_Y,],ncol=Variable_number)
  }

  Y_hat_test=NULL
  final_test=NULL
  Common_cell_test=NULL
  for (q in 1:length(Cell)){
    Common_cell_test[[q]]=intersect(finalY_test_cell,Cell[[q]])
    Cell_ind=match(Common_cell_test[[q]],Cell[[q]])
    final_test[[q]]=finalX[[q]][Cell_ind, ]
    final_test[[q]]=matrix(as.numeric(final_test[[q]]),nrow = dim(final_test[[q]])[1], ncol = dim(final_test[[q]])[2])

    Y_hat_test[[q]] = build_forest_predict(final_dataset[[q]],finalY_dataset[[q]], n_tree, m_feature, min_leaf, final_test[[q]])
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
        F_test[,RR]=F_test[,RR]+BSP_coeff[[q]][R,RR]*matrix(Y_hat_test[[W[R]]][match_ind[[R]],RR],ncol=1)
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
