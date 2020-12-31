###
# Precision Recall Set
#
# This creates a set with a given recall and precision
# @param precision This is the precision of the given set
# @param recall This is the recall of the given set
# @param length This is the length of the set
# @param baserate This is the baserate of the set
# @return returns a set with the given parameters
###
prset = function(precision, recall, length, baserate) {
  #generates a set given the different 4 quadrants
  gold1s = max(round(baserate * length , 0), 1); # TP + FN
  gold0s = length - gold1s; # FP + FN
  TP = max(round(gold1s * recall , 0) , 1); # TP
  FP = min(gold0s, max(round(TP * (1 - precision) / precision , 0), 1)); # FP
  gold = c(rep(1, gold1s), rep(0, gold0s));
  silver = c(rep(1, TP),rep(0, gold1s - TP), rep(1, FP), rep(0, gold0s - FP));

  return(cbind(gold, silver));
}

pr.ct = function(precision, recall, length, baserate) {
  #generates a set given the different 4 quadrants
  gold1s = max(round(baserate * length , 0), 1); # TP + FN
  gold0s = length - gold1s; # FP + FN
  TP = max(round(gold1s * recall , 0) , 1); # TP
  FP = min(gold0s, max(round(TP * (1 - precision) / precision , 0), 1)); # FP
  
  ct = matrix(0, ncol = 2, nrow = 2)
  ct[1,1] = TP
  ct[1,2] = gold1s - TP
  ct[2,1] = FP
  ct[2,2] = gold0s - FP
  return(ct);
}
