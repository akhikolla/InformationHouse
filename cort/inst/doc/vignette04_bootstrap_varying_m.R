## ----setup, include = FALSE---------------------------------------------------
library(cort)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7
)

## ----"getting data and parameters"--------------------------------------------
set.seed(1)
df <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
d = ncol(df)
n = nrow(df)
nb_replicates = 20 # number of replication of the resampling.
nb_fold = 5 # "k" value for the k-fold resampling method.
nb_cop = 50 # the number of m-randomized checkerboard copulas we will test.
pairs(df,lower.panel = NULL)

## ----"fitting_functions"------------------------------------------------------

make_k_fold_samples <- function(data,k = nb_fold,n_repeat = 1){
  
  sapply(1:n_repeat,function(n){
      # shuffle data : 
      data <- data[sample(nrow(data)),]
    
      # create k folds : 
      folds <- cut(seq(1,nrow(data)),breaks=k,labels=FALSE)
      
      sapply(1:k,function(i){
        test_index <- which(folds==i,arr.ind=TRUE)
        return(list(train = data[-test_index,], 
                    test = data[test_index,]))
      },simplify=FALSE)
  })
}

build_random_m <- function(how_much=1,dim = d,nrow_data){
  t(sapply(1:how_much,function(i){
    m_pos = (2:nrow_data)[nrow_data%%(2:nrow_data)==0]
    sample(m_pos,d,replace=TRUE)
  }))
}

build_all_checkerboard <- function(sampled_data,m){
  lapply(sampled_data,function(d){
    apply(m,1,function(m_){
      cbCopula(x = d$train,m = m_,pseudo=TRUE)
    })
  })
}

samples <- make_k_fold_samples(df,k=nb_fold,n_repeat=nb_replicates)
rand_m <- build_random_m(nb_cop,d,nrow(samples[[1]]$train))
cops <- build_all_checkerboard(samples,rand_m)


## -----------------------------------------------------------------------------
pEmpCop <- function(points,data=df){
  sapply(1:nrow(points),function(i){
    sum(colSums(t(data) <= points[i,]) == d)
  }) / nrow(data)
}


## -----------------------------------------------------------------------------
error <- function(cop,i,j){
  test <- samples[[i]]$test
  return(sum((pCopula(test,cop) - pEmpCop(test))^2))
  
}

errors <- sapply(1:(nb_replicates*nb_fold),function(i){
  sapply(1:nb_cop,function(j){
   error(cops[[i]][[j]],i,j)
  })
})

rmse_by_model <- sqrt(rowMeans(errors))
plot(rmse_by_model)


## -----------------------------------------------------------------------------
rand_m

## -----------------------------------------------------------------------------
convex_combination <- ConvexCombCopula(unlist(cops,recursive=FALSE),alpha = rep(1/rmse_by_model,nb_replicates*nb_fold))
simu = rCopula(1000,convex_combination)
pairs(simu,lower.panel = NULL,cex=0.5)

