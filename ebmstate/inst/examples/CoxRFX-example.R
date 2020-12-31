# Fit an empirical Bayes Cox model using
# simulated, illness-death data from 250
# patients ('mstate_data_sample').

#load simulated data
data("mstate_data_sample")

# Set class of ‘mstate_data_sample’
class(mstate_data_sample)<-c("data.frame","msdata")

# add transition matrix as attribute
tmat<-mstate::transMat(x=list(c(2,3),c(4),c(),c()),
      names=c("health","illness","death",
     "death_after_illness"))
attr(mstate_data_sample,"trans")<-tmat 

# expand covariates by transition:
covariates.expanded<-mstate::expand.covs(
      mstate_data_sample,
      covs=names(mstate_data_sample)
      [!names(mstate_data_sample)%in%c("id","from",
      "to","trans","Tstart","Tstop","time","status",
      "strata")],append=FALSE)


# argument ‘Z’ of coxrfx
Z<-data.frame(covariates.expanded,
   trans=mstate_data_sample$trans,
   strata=mstate_data_sample$trans)

# argument ‘surv’ for a non-homogeneous 
# Markov model
surv<-survival::Surv(mstate_data_sample$Tstart,
           mstate_data_sample$Tstop,
           mstate_data_sample$status)

# argument ‘groups’ of coxrfx
groups<-paste0(rep("group", ncol(Z)-2),c("_1","_2","_3"))

#fit random effects model
coxrfx_object<-CoxRFX(Z,surv,groups,tmat)

#show point estimates
summary(coxrfx_object)


