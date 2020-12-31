##################################
##  PREPARES ARGUMENTS FOR USE  ##
##################################

# CONTAINS
# get.varType() - grabs first part of string, sep specifiable  

# Anjana Commented: Moved the calculation of exSpArgs for the hold attributes to createExpSpace
# getArguments<-function(attPerturb=attPerturb, attHold=attHold, exSpArgs=exSpArgs){
#   
#   # make attSel 
#   attSel=c(attPerturb,attHold)
# 
#   if(!is.character(exSpArgs)) {
#     
#     # CHECKS FOR BOUNDS
#     boundVars<-sapply(names(exSpArgs$bounds),get.varType,USE.NAMES=FALSE,sep="_")
#     attVars<-sapply(attSel,get.varType,USE.NAMES = FALSE)
#     boundNames<-names(exSpArgs$bounds)
#     
#     # Checks that bounds are specified for each perturbed attribute
#     if(!is.null(exSpArgs$samp)){
#       if(length(exSpArgs$samp)!=length(attPerturb)){
#         stop("Samp needs to be specified for each attPerturb")
#       }
#     }
#     
#     # Code for creating a default set of historical bounds
#     # The bounds are set to '0' for attributes of Temperature (absolute values) & '1' for other variables
#     n=length(attVars)                       
#     boundsdefault=vector(length = n)
#     for (i in 1:n) {
#       if(attVars[i]=="Temp"){
#         boundsdefault[i]=0
#       } else {
#         boundsdefault[i]=1
#       }
#     }
#     tmp=as.list(boundsdefault)
#     names(tmp)=attSel
#     
#     # this fills in samp/bounds
#     # A default list of exSpArgs - specific to the particular case setup
#     exSpArgsdefault=list(type="regGrid",
#                          samp=rep(1,n),
#                          bounds=tmp)
#     
#     # Combine user specified and defaults
#     exSpArgs=modifyList(exSpArgsdefault,exSpArgs)
#     exSpArgs$samp=c(exSpArgs$samp,rep(1,length(attHold)))
#   
#   }
# 
#   return(exSpArgs)
# }
