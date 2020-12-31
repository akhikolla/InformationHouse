\name{glmscore}
\alias{glmscore}
\title{
Score Value in Generalized Linear Model
}
\description{
This function can be used to calculate the value of score function. This function can serve as a tool to assist the users in debugging their self-defined score function.
}
\usage{
glmscore(beta, Y, data, weight, offset, family)
}

\arguments{
  \item{beta}{
 a vector specifying the value of parameters.
}
  \item{Y}{
 a vector specifying the response.
}
  \item{data}{
 a matrix or data frame of covariates. If this is for score function debugging purpose, see detail for the rearrangement of the dataset.
}
  \item{weight}{
 a vector of weight.
}
  \item{offset}{
 a vector of offset.
}
  \item{family}{
 a \code{\link[stats]{glm}} ``family'' object.
}

}

\details{
In a general setting, Y can be treated as a response and data as covariates. 

This function helps to debug the user specified score function. The data are needed to be specially arranged such that the score function output can be successfully passed into \code{augSIMEX} function in a correct order. The data should be a matrix or data frame of a combination with error-prone covariates (in the order of \code{err.var} as in \code{\link[augSIMEX]{augSIMEX}}), other covariates and binary variable that is prone to misclassification. The intercept, if considered, should be included the category of ``other covariates''.
}
\value{
A matrix of n times m, where n is the number of observations and m is the number of parameters. Each entry in the matrix represents the calculated score value for subject \eqn{i} on parameter \eqn{j}. The vector of row sum will be the score value vector.
}
\references{
McCullagh P. Generalized linear models[J]. European Journal of Operational Research, 1984, 16(3): 285-292.
}
\author{
Qihuang Zhang and Grace Y. Yi.
}

\seealso{
\code{\link[stats]{glm}}
}
\examples{
### The user specified function to be checked. (logit link in binomial family)
scorefunction=function(beta,Y,data,weight,offset){
  results<-lapply(1:dim(data)[2],
                  FUN=function(i){
                    S<-lapply(1:dim(data)[1],function(x){
                      eta<- matrix(beta,nrow=1) %*% matrix(as.numeric(data[x,]),ncol=1)-offset[x]
                      return(weight[x]*Y[x]*data[x,i]-weight[x]*exp(eta)/(1+exp(eta))*data[x,i])})
                    return(S)}
  )
  return(matrix(unlist(results),ncol=dim(data)[2]))
}

data(ToyUni)

### Data need to rearranged. See detail.
nsize<-length(ToyUni$Main[,"Y"])
data.in.score<-data.frame(intercept=1,X=ToyUni$Main[,"Xstar"],
                          W=ToyUni$Main[,"W"],Z=ToyUni$Main[,"Zstar"])

## compare. The results should be identical.
glmscore(rep(0,4),ToyUni$Main[,"Y"],data.in.score,rep(1,nsize),
         rep(0,nsize),family=binomial(link=logit))
scorefunction(rep(0,4),ToyUni$Main[,"Y"],data.in.score,rep(1,nsize),rep(0,nsize))
}

\keyword{function}

