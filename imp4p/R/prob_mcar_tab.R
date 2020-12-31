
#Compute a matrix of probabilities to be MCAR

prob.mcar.tab=function(tab.u,res){
  nri=nrow(tab.u);
  nci=ncol(tab.u);
  prob=matrix(0,nri,nci);
  for (j in 1:nci){
    prob[,j]=prob.mcar(b.u=tab.u[,j], absc=res$abs.mod, pi.na=res$pi.na[j],pi.mcar=res$pi.mcar[j], F.tot=res$F.tot[,j], F.obs=res$F.obs[,j]);
  }
  return(prob)
}
