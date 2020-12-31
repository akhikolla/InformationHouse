#' @importFrom utils getFromNamespace
#' @importFrom xtable xtable
xtable.npmresult = function(x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
                            display = NULL, bylocus = FALSE, ...){
    if(bylocus){
      Loci = rownames(x$raw)
      nLoci = length(Loci)
      nContrib = length(x$summary) + 1
      
      tbl = matrix(0, nrow = nLoci, ncol = nContrib - 1)
      colnames(tbl) = paste(paste("$\\Pr\\left(", nContrib, "\\rightarrow ", sep = ""), 1:(nContrib - 1), "\\right$)", sep = "")
      
      for(i in 1:(nContrib - 1)){
        tbl[,i] = signif(x$byNC[[i]]$byLoc, digits)
      }
      
      tbl = rbind(tbl, signif(x$summary, digits))
      
      rownames(tbl) = c(Loci, "Product")  
      
      xtdf = getFromNamespace("xtable.data.frame", "xtable")
      
      return(xtdf(data.frame(tbl, check.names = FALSE), 
                               caption = caption, label = label, align = align, digits = digits, 
                               display = display, ...))
    }else{  
      nContrib = length(x$summary) + 1
      names(x$summary) = paste(paste("$\\Pr\\left(", nContrib, "\\rightarrow ", sep = ""), 1:(nContrib - 1), "\\right$)", sep = "")
      return(xtable(data.frame(x$summary, check.names = FALSE), 
                               caption = caption, label = label, align = align, digits = digits, 
                               display = display,   ...))
    }
}
