# Compute cumulative hazard rates
# under a (pre-estimated) empirical Bayes Cox
# model.

#load simulated data (illness-death model,
#500 patients) and estimated empirical
# Bayes Cox model
data("mstate_data_sample")
data("coxrfx_object_sample")

# Make objects 'surv' and 'Z'
# with the data used in the estimation

#outcome data
surv<-coxrfx_object_sample$surv

#covariate data
Z<-coxrfx_object_sample$Z

# Build a data frame 'patient_data'
# with the covariate values for which 
# cumulative hazards are to be computed
# (patient 1 covariate values in this case).
# 'patient_data' must have one row for each
# transition in the model 
# and the same columns as 'Z'. The assignment
# of transitions to strata (made in the 'strata'
# column) must follow the original model in
# 'coxrfx_object_sample'.

patient_data<-mstate_data_sample[mstate_data_sample$id==1,
   ,drop=FALSE][rep(1,3),]
patient_data$strata<-patient_data$trans<-1:3
patient_data<-mstate::expand.covs(patient_data,
   covs=names(patient_data)[!names(patient_data)%in%
   c("id","from","to","trans","Tstart","Tstop","time",
   "to","trans","Tstart","Tstop","time","status",
   "strata")],append=TRUE)

# compute cumulative hazards
msfit_object<-msfit_generic(coxrfx_object_sample,
                            patient_data,
                            coxrfx_object_sample$tmat)

# show estimates
print(msfit_object)
