\donttest{
# Compute transition probabilities
# from an object with (pre-estimated) 
# cumulative hazard rates.

#load object with estimated
#cumulative hazard rates 
data("msfit_object_sample")

#compute transition probabilities
probtrans_object<-probtrans_ebmstate("health",
   msfit_object_sample,"Markov")
}
