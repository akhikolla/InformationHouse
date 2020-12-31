#' This function calculates the prediction value of k* nearest neighbors algorithm.
#' @param Label vectors of the known labels of the samples.
#' @param Distance vectors of the distance between the target sample we want to predict and the other samples.
#' @param L_C parameter of k* nearest neighbors algorithm.
#' @return the prediction value(pred) and the weight of the samples(alpha).
#' @note This algorithm is based on Anava and Levy(2017).
#' @export
#' @examples
#'   library(ksNN)
#'   set.seed(1)
#'
#'   #make the nonlinear regression problem
#'   X<-runif(100)
#'   Y<-X^6-3*X^3+5*X^2+2
#'
#'   suffle<-order(rnorm(length(X)))
#'   X<-X[suffle]
#'   Y<-Y[suffle]
#'
#'   test_X<-X[1]
#'   test_Y<-Y[1]
#'
#'   train_X<-X[-1]
#'   train_Y<-Y[-1]
#'
#'   Label<-train_Y
#'   Distance<-sqrt((test_X-train_X)^2)
#'
#'   pred_ksNN<-ksNN(Label,Distance,L_C=1)
#'
#'   #the predicted value with k*NN
#'   pred_ksNN$pred
#'
#'   #the 'true' value
#'   test_Y
ksNN<-function(Label,Distance,L_C=1){
  
  #sort the distance in ascending order
  ord<-order(Distance)
  Label<-Label[ord]
  Distance<-Distance[ord]

  Distances <- L_C*Distance
  n <- length(Distances)
  k <- 0

  #calculate Lambda
  lambda<-Distances[1] + 1

  for (k in 1:n) {
    if(lambda[k] > Distances[k]){
      beta <- sum(Distances[1:k])
      beta2 <- sum(Distances[1:k]^2)
      tmp <- (1/k)*(beta + sqrt(k + (beta^2 - k*beta2 )))
      lambda<-c(lambda,tmp)
    }
    else{
      break
    }
  }

  #calculate Alpha
  alpha <-numeric(0)
  sum_lambda <- 0
  i <- 0
  tail_lambda <- lambda[length(lambda)]

  while(i<n){
    if(Distances[i+1] < tail_lambda){
      tmp = tail_lambda - Distances[i+1]
      sum_lambda = sum_lambda + tmp
      alpha<-c(alpha,tmp)
      i = i + 1
    }else{
      break
    }
  }

  alpha <- alpha/sum_lambda
  pred <- sum(alpha*Label[1:length(alpha)])
  return(list(pred=pred,alpha=alpha))
}
