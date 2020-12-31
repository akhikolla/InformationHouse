
# stats
setMethod("Math", signature="gridlayer", 
definition=function(x){
	methods::callGeneric(x@values)
}
)

setMethod("Summary", signature="gridlayer", 
definition=function(x){
	methods::callGeneric(x@values)
}
)


setMethod("Ops", signature=c("gridlayer", "numeric"),
definition=function(e1,e2){
	methods::callGeneric(e1@values,e2)
}
)


setMethod("Ops", signature=c("gridlayer", "gridlayer"),
definition=function(e1,e2){
	methods::callGeneric(e1@values,e2@values)
}
)

# basic statistics
setMethod("mean", signature="gridlayer", 
definition=function(x,...){
	methods::callGeneric(x@values,...)
}
)

# basic statistics
setMethod("sd", signature="gridlayer", 
definition=function(x,na.rm=FALSE){
	methods::callGeneric(x@values,na.rm=na.rm)
}
)

# basic statistics
setMethod("summary", signature="gridlayer", 
definition=function(object,...){
	methods::callGeneric(object@values,...)
}
)

# basic statistics
setMethod("quantile", signature="gridlayer", 
definition=function(x,...){
	methods::callGeneric(x@values,...)
}
)

# basic statistics
setMethod("median", signature="gridlayer", 
definition=function(x,na.rm=FALSE){
	methods::callGeneric(x@values,na.rm=na.rm)
}
)
