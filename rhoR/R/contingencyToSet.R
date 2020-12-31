contingencyToSet = function(TP, FP, FN, TN){
  setLength = TP + FP + FN + TN;

  gold1s = TP + FN;
  gold0s = setLength - gold1s;

  gold = c(rep(1, gold1s), rep(0, gold0s));
  silver = c(rep(1, TP),rep(0, gold1s - TP), rep(1, FP), rep(0, gold0s - FP));

  return(cbind(gold, silver));
};
