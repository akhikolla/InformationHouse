variable_importance_measure <- function(Model_VIM,NumVariable){
  NumRepeatation=matrix(rep(0,NumVariable),nrow=1)
  for (kk in 1:length(Model_VIM)){
    if (length(Model_VIM[[kk]])==5){
      NumRepeatation[1, Model_VIM[[kk]][3][[1]]]=NumRepeatation[1, Model_VIM[[kk]][3][[1]]]+1
    }
  }
  return(NumRepeatation)
}