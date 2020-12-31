selvarmix <- function(...) UseMethod("selvarmix")
print.selvarmix <- 
  function(object)
  {
    x <- object
    if (class(x) != "selvarmix") {
      stop(paste(sQuote("x"), sep = ""), " must be of class ", 
           paste(dQuote("selvarmix"), sep = ""), sep = "")
    }
    
    if(length(x)==2)
      for(i in 1:2)
      {
        print(x[[i]]$parameters)
        cat("Regression parameters:\n")
        print(x[[i]]$regparameters)
      }  
    else
    {
      print(x$parameters)
      cat("Regression parameters:\n")
      print(x$regparameters)
    }
  }

