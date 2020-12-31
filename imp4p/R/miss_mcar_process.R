miss.mcar.process=function(abs,pi_mcar,F_tot,F_na){
  prob.miss.mcar=NULL
  for (x in 1:(length(abs)-1)){
       prob.miss.mcar=c(prob.miss.mcar,pi_mcar*(F_tot[x+1]-F_tot[x])/(F_na[x+1]-F_na[x]));
  }
  return(list(abs=abs[1:(length(abs)-1)],p=prob.miss.mcar))
}

