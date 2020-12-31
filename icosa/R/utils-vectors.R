collapse<-function(vect){
	oldClass<-class(vect)
	factorized<-factor(vect)
	newFactor<-factor(.Call(Cpp_icosa_Collapse_, factorized))
	levels(newFactor)<-levels(factorized)
	newVect<-as.character(newFactor)
	class(newVect) <- oldClass
	return(newVect)
}
