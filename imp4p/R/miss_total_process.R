#########################
#
# Estimation du processus des valeurs manquantes
#
#########################
#F_na, F_tot, F_mcar et F_mnar sont suppos?s ?tre ?valu?s aux points d'abscisses abs
miss.total.process=function(abs,pi_na,F_na,F_tot){
  prob.miss=NULL
  for (x in 1:(length(abs)-1)){
      prob.miss=c(prob.miss,pi_na*(F_na[x+1]-F_na[x])/(F_tot[x+1]-F_tot[x]))
  }
  return(list(abs=abs[1:(length(abs)-1)],p=prob.miss))
}

