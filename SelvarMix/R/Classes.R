selvarmix <- function(bestModel)  {
  if(length(bestModel)==2)
  {
    output <- vector(mode="selvarmix",length = 2)
    for(el in 1:2)
    {
      object <- list(S=bestModel[[cl]]$S, 
                     R=bestModel[[cl]]$R, 
                     U=bestModel[[cl]]$U,
                     W=bestModel[[cl]]$W,  
                     criterionValue=bestModel[[cl]]$criterionValue, 
                     criterion=bestModel[[cl]]$criterion, 
                     model=bestModel[[cl]]$model, 
                     nbCluster=bestModel[[cl]]$nbCluster, 
                     partition=bestModel[[cl]]$partition, 
                     proba=bestModel[[cl]]$proba)
      class(object) <- "selvarmix"
      output[el] <-object 
    }
  }else
  {
    output <- list(S=bestModel[[1]]$S, 
                   R=bestModel[[1]]$R, 
                   U=bestModel[[1]]$U,
                   W=bestModel[[1]]$W,  
                   criterionValue=bestModel[[1]]$criterionValue, 
                   criterion=bestModel[[1]]$criterion, 
                   model=bestModel[[1]]$model, 
                   nbCluster=bestModel[[1]]$nbCluster, 
                   partition=bestModel[[1]]$partition, 
                   proba=bestModel[[1]]$proba)
    class(output) <- "selvarmix"
  }
  return (output)
}