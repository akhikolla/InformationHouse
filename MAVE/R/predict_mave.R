#' Make predictions based on the dimension reduction space
#'
#' This method make predictions based the reduced dimension of data using \code{\link{mars}} function.
#'
#' @aliases predict.mave
#' @aliases predict.mave.dim
#' @rdname predict.mave
#'
#' @param object the object of class 'mave'
#' @param newx Matrix of the new data to be predicted
#' @param dim the dimension of central space or central mean space. The matrix of the original data will be
#' multiplied by the matrix of dimension reduction directions of given dimension. Then the prediction will be
#' made based on the data of given dimensions. The value of dim should be given when the class of the
#' argument dr is mave. When the class of the argument dr is mave.dim and dim is not given, the
#' function will return the basis matrix of CS or CMS of dimension selected by \code{\link{mave.dim}}
#' @param ... further arguments passed to \code{\link{mars}} function such as degree.
#' @return the prediced response of the new data
#' @seealso \code{\link{mave}} for computing the dimension reduction space and \code{\link{mave.dim}} for
#' estimating the dimension of the dimension reduction space
#' @examples
#'
#' X = matrix(rnorm(10000),1000,10)
#' beta1 = as.matrix(c(1,1,1,1,0,0,0,0,0,0))
#' beta2 = as.matrix(c(0,0,0,1,1,1,1,1,0,0))
#' err = as.matrix(rnorm(1000))
#' Y = X%*%beta1+X%*%beta2+err
#'
#' train = sample(1:1000)[1:500]
#' x.train = X[train,]
#' y.train = as.matrix(Y[train])
#' x.test = X[-train,]
#' y.test = as.matrix(Y[-train])
#'
#' dr = mave(y.train~x.train, method = 'meanopg')
#'
#' yp = predict(dr,x.test,dim=3,degree=2)
#' #mean error
#' mean((yp-y.test)^2)
#'
#' dr.dim = mave.dim(dr)
#'
#' yp = predict(dr.dim,x.test,degree=2)
#' #mean error
#' mean((yp-y.test)^2)
#'
#'@method predict mave
#'@method predict mave.dim
#'@export

predict.mave<-function(object, newx, dim, ...){

  n0 = nrow(object$x)
  n1 = nrow(newx)

  x.train.mave <- mave.data(object, object$x, dim)
  x.test.mave <- mave.data(object, newx, dim)

  y.pred = matrix(0,n1,ncol(object$y))

  for(k in 1:ncol(object$y)){
    thresh0 <- 0.5/sqrt(n0)
    thresh <- 0.5/sqrt(n1)
    yk.pred <- matrix(Inf, n1, 1)
    y.ub <- max(object$y[,k]) + 0.5*(max(object$y[,k])-min(object$y[,k]))
    y.lb <- min(object$y[,k]) - 0.5*(max(object$y[,k])-min(object$y[,k]))

    idx <- which((yk.pred>y.ub)|(yk.pred<y.lb))
    counter <- 0
    while(length(idx)>0){
      if(any(substr(rep("thresh", 5), 1, 2:6) %in% names(list(...)))){
        model <- mda::mars(x.train.mave,object$y[,k],...)
      }
      else{
        model <- mda::mars(x.train.mave,object$y[,k],thresh=thresh,...)
      }
      yk.pred[idx] <- predict(model,x.test.mave[idx,,drop=F])
      idx <- which( (yk.pred>y.ub) | (yk.pred<y.lb) )
      thresh = thresh + thresh0

      counter <- counter+1
      if(counter>100){
        break
      }
    }
    y.pred[,k] = yk.pred
  }
  return(y.pred)
}

#' @rdname predict.mave
#' @export
predict.mave.dim<-function(object,newx,dim='dim.min',...){
  if(dim=='dim.min'){
    dim <- which.min(object$cv)
    y.pred <- predict.mave(object,newx,dim,...)
  }else{
    y.pred <- predict.mave(object,newx,dim,...)
  }
  return(y.pred)

}
