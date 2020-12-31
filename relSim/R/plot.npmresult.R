#' @importFrom graphics barplot boxplot box
plot.npmresult = function(x, bylocus = FALSE, ...){
    if(bylocus){
      Loci = rownames(x$raw)
      nLoci = length(Loci)
      nContrib = length(x$summary) + 1
      
      tbl = matrix(0, nrow = nLoci, ncol = nContrib - 1)
      colnames(tbl) = paste(paste("Pr(", nContrib, "->", sep = ""), 1:(nContrib - 1), ")", sep = "")
      
      for(i in 1:(nContrib - 1)){
        tbl[,i] = x$byNC[[i]]$byLoc
      }
      
      tbl = rbind(tbl, x$summary)
      
      rownames(tbl) = c(Loci, "Product")    
      barplot(tbl, beside = T, ...)
      box()
    }else{  
      nContrib = length(x$summary) + 1
      names(x$summary) = paste(paste("Pr(", nContrib, "->", sep = ""), 1:(nContrib - 1), ")", sep = "")
      barplot(x$summary, ...)
      box()
    }
}
