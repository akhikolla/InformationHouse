`summary.selvarmix` <- 
  function(object, ...)
  {
    x <- object
    if (class(x) != "selvarmix") {
      stop(paste(sQuote("x"), sep = ""), " must be of class ", 
           paste(dQuote("selvarmix"), sep = ""), sep = "")
    }
    if(length(x)==2)
      for(i in 1:2)
      {
        cat("Criterion:", x[[i]]$criterion, "\n")
        cat("Criterion value:", x[[i]]$criterionValue,"\n")
        cat("Number of clusters:", x[[i]]$nbcluster,"\n")
        cat("Gaussian mixture model:", x[[i]]$model,"\n")
        cat("Regression covariance model:", x[[i]]$rmodel,"\n")
        cat("Independent covariance model:", x[[i]]$imodel,"\n")
        cat("The SRUW model:\n")
        cat(" S:", x[[i]]$S,"\n")
        cat(" R:", x[[i]]$R,"\n")
        cat(" U:", x[[i]]$U,"\n")
        cat(" W:", x[[i]]$W,"\n")
      }
    else{
      if(is.null(x$error))
      {
        cat("Criterion:", x$criterion, "\n")
        cat("Criterion value:", x$criterionValue,"\n")
        cat("Number of clusters:", x$nbcluster,"\n")
        cat("Gaussian mixture model:", x$model,"\n")
        cat("Regression covariance model:", x$rmodel,"\n")
        cat("Independent covariance model:", x$imodel,"\n")
        cat("The SRUW model:\n")
        cat(" S:", x$S,"\n")
        cat(" R:", x$R,"\n")
        cat(" U:", x$U,"\n")
        cat(" W:", x$W,"\n")
      }
      else
      {
        cat("Criterion:", x$criterion, "\n")
        cat("Criterion value:", x$criterionValue,"\n")
        cat("Number of clusters:", x$nbcluster,"\n")
        cat("Gaussian mixture model:", x$model,"\n")
        cat("Prediction error:", round(x$error, 2),"\n")
        cat("Regression covariance model:", x$rmodel,"\n")
        cat("Independent covariance model:", x$imodel,"\n")
        cat("The SRUW model:\n")
        cat(" S:", x$S,"\n")
        cat(" R:", x$R,"\n")
        cat(" U:", x$U,"\n")
        cat(" W:", x$W,"\n")
        
        
      }  
    }
  }