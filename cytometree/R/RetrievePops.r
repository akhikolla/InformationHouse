#' Retrieve cell populations found using Annotation.
#' 
#'@param AnnotationObj An object of class Annotation.
#'
#'@param phenotypes A list containing at least one element of class 
#'matrix describing a sought phenotype. Each matrix should have two
#'columns where the name of a used marker is associated to a value
#'chosen between 0, 1 and 2. 0 translates to -, 1 to + and 2 to ++.
#' 
#'@return A \code{list} of two elements.
#'\itemize{
#'\item{\code{phenotypesinfo}}{ A \code{list} containing informations 
#'about sought populations.}
#'\item{\code{Mergedleaves}}{ The partitioning of the set of n cells
#'with potentially merged leaves.}
#'}
#'
#'@importFrom methods is
#'
#'@author Chariff Alkhassim, Boris Hejblum
#'
#'@export 
RetrievePops <-function(AnnotationObj, phenotypes)
{
  if(!methods::is(AnnotationObj, "Annotation")){
    stop("AnnotationObj must be of class Annotation.")
  }
  if(methods::is(phenotypes, "list")){
    if(!methods::is(phenotypes[[1]], "matrix")){
      stop("Elements of phenotypes should be of class matrix.") 
    }
  }
  else{
    stop("phenotypes should be of class list.")
  }
  labels <- AnnotationObj$labels
  labelmerge <- labels
  combinations <- AnnotationObj$combinations
  colnames_combinations <- colnames(combinations)
  maxlab <- max(labels)
  L <- length(phenotypes)
  outlist <- list()
  Prop <- combinations[, c("prop")]
  for(l in 1:L){
    outlist[[l]] <- list()
    temp <- phenotypes[[l]]
    inputest <- temp[, 1] %in% colnames_combinations
    if(any(!(inputest))){
      logicalinds <- as.logical(1-inputest)
      if(sum(logicalinds) > 1){
        wstr <- paste("markers",paste(c(temp[logicalinds, 1]), 
                                      collapse=", "),
                      "are not in CytomeTree.", sep = " ")
      }else{
        wstr <- paste("marker",paste(c(temp[logicalinds,1]), 
                                      collapse=", "),
                      "is not in CytomeTree.", sep = " ")
      }
      stop(wstr)
    }
    tempcombinations <- combinations[, temp[,1], drop=FALSE]
    selected_pop <- apply(tempcombinations, 1, 
                        FUN = function(x,y){ 
                          x == y
                        }, y = as.numeric(temp[, 2]))
    if(is.null(dim(selected_pop))){
      scores <- 1*selected_pop
    }else{
      scores <- rowSums(t(selected_pop))
    }
    tempres <- which(scores == nrow(temp))
    if(!length(tempres)){
      outlist[[l]][["phenotype"]] <- apply(temp[, 1:2, drop=FALSE], 1, paste, collapse="-") 
      outlist[[l]][["proportion"]] <- NA
      outlist[[l]][["Mergedlabels"]] <- NA
      outlist[[l]][["Newlabel"]] <- NA
    }
    else {
      leaves <- combinations[,c("leaves")][tempres]
      outlist[[l]][["phenotype"]] <- apply(temp[,1:2, drop=FALSE],1,paste,collapse="-") 
      outlist[[l]][["proportion"]] <- sum(Prop[tempres])
      if(length(tempres) > 1) {
        outlist[[l]][["Mergedlabels"]] <- leaves
        outlist[[l]][["Newlabel"]] <- maxlab + 1
        labelmerge[labels%in%leaves] <- maxlab + 1
        maxlab <- maxlab + 1
      } else {
        outlist[[l]][["label"]] <- leaves
      }
    }
  }
  return(list("phenotypesinfo" = outlist, "Mergedleaves" = labelmerge))
}
