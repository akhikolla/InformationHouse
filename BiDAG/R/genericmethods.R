#this script contains definition of classes in BiDAG
#and all summary and print functions for these classes

setClass("MCMCmax",
         slots = c(DAG = "matrix",
                   score = "numeric",
                   reach = "integer",
                   order= "vector"))
setClass("MCMCtrace",
         slots = c(incidence = "list",
                   DAGscores = "list",
                   partitionscores = "list",
                   orderscores = "list",
                   order= "list",
                   partiton="list"))
setClass("MCMCspace",
         slots = c(adjacency = "matrix",
                   scoretables = "list"))

setClass("MCMCres",
         slots = c(max = "MCMCmax",
                   trace = "MCMCtrace",
                   space = "MCMCspace",
                   info="list"))

setClass("MCMCmult",
         slots = c(max = "MCMCmax",
                   maxtrace = "list",
                   trace = "MCMCtrace",
                   space = "MCMCspace",
                   info="list"))

#summary methods
#setGeneric("summary")

summary.default <- function(object, ...) { 
  base::summary(object, ...)
}

#'@method summary MCMCres
#'@export
summary.MCMCres <- function(object, ...) {
  cat(paste("algorithm:",object$info$algo,"\n"))
  cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
  cat(paste("initial search space:",object$info$spacealgo,"\n"))
  cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
  cat(paste("number of sampling steps:",object$info$samplesteps,"\n"))
  cat(paste("score of the maximum score DAG:",round(object$max$score,2),"\n"))
  cat(paste("maximum reached at:",object$max$reach, "step out of", object$info$samplesteps,"\n"))
  
  if(is.null(object$chain)) {
    cat("MCMC trace was not saved","\n")
  } else {
    cat(paste("MCMC trace contains:",length(object$chain$incidence),"saved steps","\n"))
  }
  if(!is.null(object$space)) {
    cat("additional objects returned: the search space and the score tables","\n")
  }
}
setMethod("summary", "MCMCres", summary.MCMCres)

#'@method summary MCMCmult
#'@export
summary.MCMCmult <- function(object, ...) {
  cat("MCMC settings\n")
  cat(paste("algorithm:",object$info$algo,"\n"))
  cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
  cat(paste("initial search space:",object$info$spacealgo,"\n"))
  cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
  cat(paste("number of sample steps:",object$info$samplesteps,"\n"))
  cat("\n")
  
  cat("   Search space optimization\n")
  cat(paste("number of space expansion steps:",length(object$maxtrace),"\n"))
  cat(paste("number of edges in the intial search space:",sum(object$startspace),"\n"))
  cat(paste("number of added edges:",sum(object$endspace)-sum(object$startspace),"\n"))
  cat("\n")
  
  cat("   MAP estimate\n")
  cat(paste("score of the maximum score DAG:",round(object$max$score,2),"\n"))
  cat("\n")
  
  cat("   Additional objects\n")
  if(is.null(object$chain)) {
  } else {
    cat(paste("MCMC trace, contains ",length(object$maxtrace)," x ",length(object$chain$incidence[[1]]),"saved steps","\n"))
  }
  if(!is.null(object$scoretable)) {
    cat("score tables","\n")
  }
}
setMethod("summary", "MCMCmult", summary.MCMCmult)

# setMethod("summary", "MCMCres",
#           function(object, ...) {
#             cat(paste("algorithm:",object$info$algo,"\n"))
#             cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
#             cat(paste("initial search space:",object$info$spacealgo,"\n"))
#             cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
#             cat(paste("number of sampling steps:",object$info$samplesteps,"\n"))
#             cat(paste("score of the maximum score DAG:",round(object$max$score,2),"\n"))
#             cat(paste("maximum reached at:",object$max$reach, "step out of", object$info$samplesteps,"\n"))
#             
#             if(is.null(object$chain)) {
#               cat("MCMC trace was not saved","\n")
#             } else {
#               cat(paste("MCMC trace contains:",length(object$chain$incidence),"saved steps","\n"))
#             }
#             if(!is.null(object$space)) {
#               cat("additional objects returned: the search space and the score tables","\n")
#             }
#           })


# setMethod("summary", "MCMCmult",
#           function(object, ...) {
#             cat("MCMC settings\n")
#             cat(paste("algorithm:",object$info$algo,"\n"))
#             cat(paste("sample/MAP: ",object$info$sampletype,"\n"))
#             cat(paste("initial search space:",object$info$spacealgo,"\n"))
#             cat(paste("number of MCMC iterations:",object$info$iterations,"\n"))
#             cat(paste("number of sample steps:",object$info$samplesteps,"\n"))
#             cat("\n")
#             
#             cat("   Search space optimization\n")
#             cat(paste("number of space expansion steps:",length(object$maxtrace),"\n"))
#             cat(paste("number of edges in the intial search space:",sum(object$startspace),"\n"))
#             cat(paste("number of added edges:",sum(object$endspace)-sum(object$startspace),"\n"))
#             cat("\n")
#             
#             cat("   MAP estimate\n")
#             cat(paste("score of the maximum score DAG:",round(object$max$score,2),"\n"))
#             cat("\n")
#             
#             cat("   Additional objects\n")
#             if(is.null(object$chain)) {
#             } else {
#               cat(paste("MCMC trace, contains ",length(object$maxtrace)," x ",length(object$chain$incidence[[1]]),"saved steps","\n"))
#             }
#             if(!is.null(object$scoretable)) {
#               cat("score tables","\n")
#             }
#           })


#PRINT methods



#'@export
#
print.MCMCmax <-function(x,...){
  
  cat("Number of edges in the maximum score DAG: ", sum(x$DAG), "\n")
  cat("Maximum DAG score: ", x$score, "\n")
  cat("Order: ", x$order, "\n")
  
}

#'@export
#a generic method (does not need description) print for scoreparameters class
print.scoreparameters <-function(x,...){
  if (x$type=="bge") {
    cat("Score type: BGe", "\n")
  } else if (x$type=="bde") {
    cat("Score type: BDe", "\n")
    cat("Prior pseudo counts:", x$chi,"\n")
    cat("Edge penalization factor:", x$pf,"\n")
  } else if (x$type=="bdecat") {
    cat("Score type: BDe for categorical data", "\n")
    cat("Prior pseudo counts:", x$chi,"\n")
    cat("Edge penalization factor:", x$pf,"\n")
  } else if (x$type=="dbn") {
    cat("Score type: DBN", "\n")
    cat("DBN score type: ", "\n")
    } else {
    cat("Score type: user defined", "\n")
  }
  if(is.null(x$weightvector)) {
    cat("Data is not weighted\n")
  } else {
    "Data is weighted according to the weight vector"
  }
  cat("Score constant vector:", x$scoreconstvec,"\n")
}

#'@export
#a standard print method for function compareDAGs
print.compDAGs<-function(x,...) {
  cat("TP: ", x[2], "\n")
  cat("TPR: ", x[5], "\n")
  cat("SHD: ", x[1], "\n")
  cat("FP: ", x[3], "\n")
  cat("FN: ", x[4], "\n")
}

#'@export
#a standard print method for function MCMCtrace
print.MCMCtrace <-function(x,...){
  cat("Number of saved steps: ", length(x$incidence), "\n")
}


#plot functions

# setMethod("plot", "MCMCmult",
#           function(x, ...) {
#             x<-x$chain
#             if(is.null(x)) {
#               stop("no saved MCMC steps found! set chainout=TRUE")
#             }
#             
#             
#             col5<-c("#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84")
#             ncol<-length(col5)
#             nchains<-length(x$DAGscores)
#             scorevecmin<-unlist(x$DAGscores[[1]])
#             scorevecprevmax<-unlist(x$DAGscores[[nchains-1]])
#             minprev<-min(scorevecprevmax)
#             scorevecmax<-unlist(x$DAGscores[[nchains]])
#             vecl<-length(scorevecmin)
#             scoremax<-max(scorevecmax)
#             scoremin<-min(scorevecmin)
#             scoremaxmin<-min(scorevecmax)
#             colvec<-assigncolor(nchains,ncol)
#             par(mfrow=c(1,2))
#             par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
#             plot(scorevecmin,type="l",xlab="iteration",ylab="logscore",
#                  ylim=c(scoremin,scoremax),main="DAG scores: all expansion iterations",cex.main=1.2,
#                  col=col5[colvec[1]])
#             for(i in 2:nchains) {
#               lines(unlist(x$DAGscores[[i]]),col=col5[colvec[i]])
#             }
#             par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
#             plot(c(scorevecprevmax,scorevecmax),type="l",xlab="iteration",ylab="logscore",
#                  ylim=c(minprev,scoremax),main="last iterations",cex.main=1.2,col=col5[ncol])
#           })

setMethod("plot", "MCMCmult",
          function(x, ...) {
            x<-x$chain
            if(is.null(x)) {
              stop("no saved MCMC steps found! set chainout=TRUE")
            }
            
            
            nchains<-length(x$DAGscores)
            scorevecmin<-unlist(x$DAGscores[[1]])
            scorevecprevmax<-unlist(x$DAGscores[[nchains-1]])
            minprev<-min(scorevecprevmax)
            scorevecmax<-unlist(x$DAGscores[[nchains]])
            vecl<-length(scorevecmin)
            scoremax<-max(scorevecmax)
            scoremin<-min(scorevecmin)
            scoremaxmin<-min(scorevecmax)
            # Add extra space to right of plot area; change clipping to figure
            par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
            
            #par(mfrow=c(1,2))
            #par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
            plot(scorevecmin,type="l",xlab="iteration",ylab="DAG logscores",
                 ylim=c(scoremin,scoremax),main="DAG scores: all expansion iterations",cex.main=1.2,
                 col=nchains)
            for(i in 2:nchains) {
              lines(unlist(x$DAGscores[[i]]),col=nchains-i+1)
            }
            legend("topright", inset=c(-0.2,0), legend=paste(1:nchains, " "), col=c(nchains:1), lty=rep(1,nchains), title="expansion \niteration:",bty="n")
            
          })

assigncolor<-function(nit,ncol) {
  colind<-1:(ncol-1)
  ncolsmall<-ncol-1
  nitsmall<-nit-2
  colvec<-vector()
  
  if (nit<=ncol+1) {
    if(nit<3){
      return(rep(ncol,nit))
    } else {
      return(c(tail(colind,nitsmall),ncol,ncol))
    }
  } else {
    rmndr <- nitsmall %%  ncolsmall
    div<-nitsmall %/% ncolsmall
    if(rmndr==0) firstadd<-ncol else firstadd<-tail(colind,rmndr)[1]
    for(i in 1:ncolsmall) {
      if(i>=firstadd) {
        colvec<-c(colvec,rep(colind[i],div+1))
      } else {
        colvec<-c(colvec,rep(colind[i],div))
      }
    }
  }
  return(c(colvec,ncol,ncol))
}

setMethod("plot", "MCMCres",
          function(x, ...,burnin=0.2) {
            x<-x$chain
            if(is.null(x)) {
              stop("no saved MCMC steps found! set chainout=TRUE")
            }
            scorevec<-unlist(x$DAGscores)
            vecl<-length(scorevec)
            burnin<-ceiling(vecl*burnin)
            score20<-min(scorevec[burnin:vecl])
            scoremax<-max(scorevec)
            scoremin<-min(scorevec)
            
            par(mfrow=c(1,2))
            par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
            plot(scorevec,type="l",col="#0c2c84",xlab="iteration",ylab="logscore",
                 ylim=c(scoremin,scoremax+(scoremax-scoremin)*0.02),main="DAG logscores",cex.main=1)
            par(mar = c(2, 2, 2, 2)) # Set the margin on all sides to 2
            plot(x=c(burnin:vecl),y=scorevec[burnin:vecl],type="l",col="#0c2c84",xlab="iteration",ylab="logscore",
                 ylim=c(score20,scoremax),main="DAG logscores excluding burn-in",cex.main=1)
            
          })



