print.npmresult = function(x, bylocus = FALSE, digits = 3, ...){
    if(bylocus){
      Loci = rownames(x$raw)
      nLoci = length(Loci)
      nContrib = length(x$summary) + 1
      
      tbl = matrix(0, nrow = nLoci, ncol = nContrib - 1)
      colnames(tbl) = paste(paste("Pr(", nContrib, "->", sep = ""), 1:(nContrib - 1), ")", sep = "")
      
      for(i in 1:(nContrib - 1)){
        tbl[,i] = signif(x$byNC[[i]]$byLoc, digits)
      }
      
      tbl = rbind(tbl, signif(x$summary, digits))
      
      rownames(tbl) = c(Loci, "Product")    
      print(tbl, ...)
      cat("\n")
    }else{  
      nContrib = length(x$summary) + 1
      names(x$summary) = paste(paste("Pr(", nContrib, "->", sep = ""), 1:(nContrib - 1), ")", sep = "")
      print(signif(x$summary, digits, ...))
    }
}
