extrapolate <- function (extrapolation,betamatrix,lambda,nbeta){
  if (extrapolation=="quadratic"){
    extrapomodel<-apply(betamatrix,MARGIN=2, FUN=function(x){
      lambda2<-lambda^2
      model<-lm(x~lambda+lambda2)
      newdata<-data.frame(lambda=-1,lambda2=1)
      betaresults<-predict(model,newdata)
    })
    coefs <- extrapomodel[1:nbeta]
  } else if (extrapolation=="linear"){
    extrapomodel<-apply(betamatrix,MARGIN=2, FUN=function(x){
      model<-lm(x~lambda)
      newdata<-data.frame(lambda=-1)
      betaresults<-predict(model,newdata)
    })
    coefs <- extrapomodel[1:nbeta]
  } else if (extrapolation=="both"){
    extrapomodel1<-apply(betamatrix,MARGIN=2, FUN=function(x){
      lambda2<-lambda^2
      model<-lm(x~lambda+lambda2,na.action=na.omit)
      newdata<-data.frame(lambda=-1,lambda2=1)
      betaresults<-predict(model,newdata)
    })
    
    extrapomodel2<-apply(betamatrix,MARGIN=2, FUN=function(x){
      model<-lm(x~lambda)
      newdata<-data.frame(lambda=-1)
      betaresults<-predict(model,newdata)
    })
    
    coefs <- c(extrapomodel1[1:nbeta],extrapomodel2[1:nbeta])
  }
  
  return(coefs)
}
  


fastapproach<-function(family,cpp){
  if (missing(family)) return(list(scorefun=score.modifieduser,fast=cpp))
  if (cpp==TRUE){
    if (family$family=="gaussian"){
      if (family$link=="identity") {return(list(scorefun=scoregaussiancpp,fast=TRUE))}
    }
    if (family$family=="binomial"){
      if (family$link=="logit") {return(list(scorefun=scorelogitcpp,fast=TRUE))}
      if (family$link=="probit") {return(list(scorefun=scoreprobitcpp,fast=TRUE))}
      if (family$link=="cloglog") {return(list(scorefun=scorecloglogcpp,fast=TRUE))}
    }
    if (family$family=="poisson"){
      if (family$link=="sqrt") {return(list(scorefun=scorepoisqrtcpp,fast=TRUE))}
      return(list(scorefun=scorepoilogcpp,fast=TRUE))
    }
  }
  return(list(scorefun=score.modifiedglm,fast=FALSE))
}

Getalpha<-function(ValidationData,mismodel=NULL,mis.var,mis.true,err.var){
  if (length(mismodel)[1]!=2) stop("The number of responses in misclassification model is incorrect. Use | to seperate responses.")
  if (length(mismodel)[2]>2) stop("The number of responses model should not be greater than 2.")
  repsname <- all.vars(mismodel)[1:2]
  if (!all(repsname %in% c("qi","pi"))) stop("The specification of reponse in misclassification should be pi and qi.")
  if (length(mismodel)[2]==1){
    for (i in 1:2){
      if (repsname[i]=="pi") {pimodel<- formula(mismodel,lhs=i,rhs=1)}
       else {qimodel<- formula(mismodel,lhs=i,rhs=1)}
    }
  } else {
    for (i in 1:2){
      if (repsname[i]=="pi") {pimodel<- formula(mismodel,lhs=i,rhs=i)}
      else {qimodel<- formula(mismodel,lhs=i,rhs=i)}
    }
  }
  
  ValidationZ1<-ValidationData[ValidationData[,mis.true]==1,]
  ValidationZ1.data<-data.frame(pi=1-(ValidationZ1[,mis.var]==ValidationZ1[,mis.true])*1,ValidationZ1)
  ValidationZ0<-ValidationData[ValidationData[,mis.true]==0,]
  ValidationZ0.data<-data.frame(qi=1-(ValidationZ0[,mis.var]==ValidationZ0[,mis.true])*1,ValidationZ0)
  
  alphahat1<-glm(pimodel,family=binomial,data=ValidationZ1.data)
  alphahat2<-glm(qimodel,family=binomial,data=ValidationZ0.data)
  return(list(alphahat1=alphahat1,alphahat2=alphahat2))
}

imputeX<-function(maindata,err.true,err.var,lambda,Sigma_ehat,nsize,repeated,repind){
  if (!repeated){
  e_ib<-mvrnorm(n=nsize,mu=rep(0,length=dim(Sigma_ehat)[1]),Sigma=Sigma_ehat)
    return(as.matrix(maindata[,err.true])+e_ib*sqrt(lambda))} else {
      Wbi<-lapply(1:length(repind),FUN=function(i){
        if (is.null(names(repind))) repvar<-repind[[i]] else repvar<-repind[[err.var[i]]]
        return(repeatgenerate(as.matrix(maindata[,repvar]),lambda))
      })
      return(matrix(unlist(Wbi),ncol=length(repind),byrow=T))
    }
}

score.glm<-function(beta,Y,DataM,weight,offset,linkinv,var,mueta){
  results<-lapply(1:dim(DataM)[2],FUN=function(i){
    S<-apply(cbind(Y,DataM),MARGIN=1,function(x){
      eta<- matrix(beta,nrow=1) %*% matrix(as.numeric(x[2:length(x)]),ncol=1)
      mu<-linkinv(eta)
      return(weight[i]*mueta(eta)/(var(mu))*(x[1]-mu)* x[i+1])})
    return(S)}
  )
  return(matrix(unlist(results),ncol=dim(DataM)[2]))
}

score.modifiedglm<-function(beta,Y,DataM,DataM0,DataM1,phat,qhat,weight,offset,linkinv,var,mueta){
  ScoreZ0<-score.glm(beta,Y,DataM0,weight,offset,linkinv,var,mueta)
  ScoreZ1<-score.glm(beta,Y,DataM1,weight,offset,linkinv,var,mueta)
  Gere<-lapply(1:(dim(DataM)[2]),FUN = function(i){
    value<-((1-DataM[,dim(DataM)[2]])*(ScoreZ0[,i]*(1-phat)-ScoreZ1[,i]*qhat)-DataM[,dim(DataM)[2]]*(ScoreZ0[,i]*phat-ScoreZ1[,i]*(1-qhat)))/(1-phat-qhat)
    return(sum(value))
  })
  return(unlist(Gere))
}


score.modifieduser<-function(beta,Y,DataM,DataM0,DataM1,phat,qhat,weight,offset,sfun){
  ScoreZ0<-sfun(beta,Y,DataM0,weight,offset)
  ScoreZ1<-sfun(beta,Y,DataM1,weight,offset)
  Gere<-lapply(1:(dim(DataM)[2]),FUN = function(i){
    value<-((1-DataM[,dim(DataM)[2]])*(ScoreZ0[,i]*(1-phat)-ScoreZ1[,i]*qhat)-DataM[,dim(DataM)[2]]*(ScoreZ0[,i]*phat-ScoreZ1[,i]*(1-qhat)))/(1-phat-qhat)
    return(sum(value))
  })
  return(unlist(Gere))
}