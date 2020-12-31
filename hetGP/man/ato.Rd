\name{ato}
\alias{ato}
\alias{X}
\alias{Z}
\alias{Xa}
\alias{Xtrain}
\alias{Xtrain.out}
\alias{Ztrain}
\alias{Ztrain.out}
\alias{train}
\alias{mult}
\alias{reps}
\alias{kill}
\alias{Zm}
\alias{Zv}
\alias{out}
\alias{Xtest}
\alias{Ztest}
\alias{nc}
\alias{ato.a}
\alias{Xa}
\alias{Za}
\alias{out.a}
\docType{data}
\title{ Assemble To Order (ATO) Data and Fits }
\description{
  A batch design-evaluated ATO data set, random partition into training and
  testing, and fitted \pkg{hetGP} model; similarly a sequentially designed
  adaptive horizon data set, and associated fitted \pkg{hetGP} model 
}
\usage{data(ato)}
\format{
  Calling \code{data(ato)} causes the following objects to be loaded into the namespace.

  \describe{
    \item{\code{X}}{ 2000x8 \code{matrix} of inputs coded from 1,...,20 to the unit 8-cube; original inputs can be recreated as \code{X*19 + 1} }
    \item{\code{Z}}{ 2000x10 \code{matrix} of normalized outputs obtained in ten replicates at each of the 2000 inputs \code{X}.  Original outputs can be obtained as \code{Z*sqrt(Zv) + Zm} }
    \item{\code{Zm}}{ scalar mean used to normalize \code{Z} }
    \item{\code{Zv}}{ scalar variance used to normalize \code{Z} }
    \item{\code{train}}{ vector of 1000 rows of \code{X} and \code{Z} selected for training }
    \item{\code{Xtrain}}{ 1000x8 \code{matrix} obtained as a random partition of \code{X} }
    \item{\code{Ztrain}}{ length 1000 list of vectors containing the selected (replicated) observations at each row of \code{Xtrain} }
    \item{\code{mult}}{ the length of each entry of \code{Ztrain}; same as \code{unlist(lapply(Ztrain, length))} }
    \item{\code{kill}}{ a \code{logical} vector indicating which rows of \code{Xtrain} for which all replicates of \code{Z} are selected for \code{Ztrain};  same as \code{mult == 10}}
    \item{\code{Xtrain.out}}{ 897x8 \code{matrix} comprised of the subset of \code{X} where not all replicates are selected for training; i.e., those for which \code{kill == FALSE} }
    \item{\code{Ztrain.out}}{ \code{list} of length 897 containing the replicates of \code{Z} not selected for \code{Ztrain} }
    \item{\code{nc}}{ \code{noiseControl} argument for \code{mleHetGP} call }
    \item{\code{out}}{ \code{mleHetGP} model based on \code{Xtrain} and \code{Ztrain} using \code{noiseControl=nc} }
    \item{\code{Xtest}}{ 1000x8 \code{matrix} containing the other partition of \code{X} of locations not selected for training }
    \item{\code{Ztest}}{ 1000x10 \code{matrix} of responses from the partition of \code{Z} not selected for training }
    \item{\code{ato.a}}{ 2000x9 \code{matrix} of sequentially designed inputs (8) and outputs (1) obtained under an adaptive horizon scheme }
    \item{\code{Xa}}{ 2000x8 matrix of coded inputs from \code{ato.a} as \code{(ato.a[,1:8]-1)/19} }
    \item{\code{Za}}{ length 2000 vector of outputs from \code{ato.a} as \code{(ato.a[,9] - Zm)/sqrt(Zv)} }
    \item{\code{out.a}}{ \code{mleHetGP} model based on \code{Xa} and \code{Za} using \code{noiseControl=nc} }
  }
}
\details{
  The assemble to order (ATO) simulator  (Hong, Nelson, 2006) is a queuing
  simulation targeting inventory management scenarios.  The setup is as follows.
  A company manufactures \eqn{m}{m} products.  Products are built from base parts
  called items, some of which are \dQuote{key} in that the product cannot be built
  without them.  If a random request comes in for a product that is missing a
  key item, a replenishment order is executed, and is filled after a random
  period.  Holding items in inventory is expensive, so there is a balance
  between inventory costs and revenue. Hong & Nelson built a
  \code{Matlab} simulator for this setup, which was subsequently
  reimplemented by Xie, et al., (2012).  

  Binois, et al (2018a) describe an out-of-sample experiment based on this
  latter implementation in its default (Hong & Nelson) setting, specifying
  item cost structure, product makeup (their items) and revenue, distribution
  of demand and replenishment time, under target stock vector inputs \eqn{b \in
  \{1,\dots,20\}^8}{b=(1,...,20)^8} for eight items.  They worked with 2000
  random uniform input locations (\code{X}), and ten replicate responses at
  each location (\code{Z}). The partition of 1000 training data points
  (\code{Xtrain} and
  \code{Ztrain}) and 1000 testing (\code{Xtest} and \code{Ztest}) sets
  provided here is an example of one that was used for the Monte Carlo
  experiment in that paper.  The elements \code{Xtrain.out} and
  \code{Ztrain.out} comprise of replicates from the training inputs which were
  not used in training, so may be used for out-of-sample testing.  For more
  details on how the partitions were build, see the code in the examples
  section below.

  Binois, et al (2018b) describe an adaptive lookahead horizon scheme for
  building a sequential design (\code{Xa}, \code{Za}) of size 2000 whose
  predictive performance, via proper scores, is almost as good as the
  approximately 5000 training data sites in each of the Monte Carlo
  repetitions described above.  The example code below demonstrates this via
  out-of-sample predictions on \code{Xtest} (measured against \code{Ztest})
  when \code{Xtrain} and \code{Ztrain} are used compared to those from
  \code{Xa} and \code{Za}.
}

\note{ The \code{mleHetGP} output objects were build with
\code{return.matrices=FALSE} for more compact storage.  Before these objects
can be used for calculations, e.g., prediction or design, these covariance
matrices need to be rebuilt with \code{\link{rebuild}}.  The generic
\code{predict} method will call \code{\link{rebuild}} automatically, 
however, some of the other methods will not, and it is often more
efficient to call \code{\link{rebuild}} once at the outset, rather
than for every subsequent \code{predict} call }

\author{ 
Mickael Binois, \email{mbinois@mcs.anl.gov}, and
Robert B. Gramacy, \email{rbg@vt.edu}
}

\references{
  Hong L., Nelson B. (2006), Discrete optimization via simulation using COMPASS. Operations Research, 54(1), 115-129.

  Xie J., Frazier P., Chick S. (2012). Assemble to Order Simulator. \url{http://simopt.org/wiki/index.php?title=Assemble_to_Order&oldid=447}.

  M. Binois, J. Huang, R. Gramacy, M. Ludkovski (2018a), Replication or exploration? Sequential design for stochastic simulation experiments,
  arXiv preprint arXiv:1710.03206.

  M. Binois, Robert B. Gramacy, M. Ludkovski (2018b), Practical heteroskedastic Gaussian process modeling for large simulation experiments,
  arXiv preprint arXiv:1611.05902.
}

\examples{
data(ato)

\dontrun{
##
## the code below was used to create the random partition 
##

## recover the data in its original form
X <- X*19+1
Z <- Z*sqrt(Zv) + Zm

## code the inputs and outputs; i.e., undo the transformation 
## above
X <- (X-1)/19
Zm <- mean(Z)
Zv <- var(as.vector(Z))
Z <- (Z - Zm)/sqrt(Zv)

## random training and testing partition
train <- sample(1:nrow(X), 1000)
Xtrain <- X[train,]
Xtest <- X[-train,]
Ztest <- as.list(as.data.frame(t(Z[-train,])))
Ztrain <- Ztrain.out <- list()
mult <- rep(NA, nrow(Xtrain))
kill <- rep(FALSE, nrow(Xtrain))
for(i in 1:length(train)) {
  reps <- sample(1:ncol(Z), 1)
  w <- sample(1:ncol(Z), reps)
  Ztrain[[i]] <- Z[train[i],w]
  if(reps < 10) Ztrain.out[[i]] <- Z[train[i],-w]
  else kill[i] <- TRUE
  mult[i] <- reps
}

## calculate training locations and outputs for replicates not
## included in Ztrain
Xtrain.out <- Xtrain[!kill,]
Ztrain.out <- Ztrain[which(!kill)]

## fit hetGP model
out <- mleHetGP(X=list(X0=Xtrain, Z0=sapply(Ztrain, mean), mult=mult),
  Z=unlist(Ztrain), lower=rep(0.01, ncol(X)), upper=rep(30, ncol(X)),
  covtype="Matern5_2", noiseControl=nc, known=list(beta0=0), 
  maxit=100000, settings=list(return.matrices=FALSE))

##
## the adaptive lookahead design is read in and fit as 
## follows
##
Xa <- (ato.a[,1:8]-1)/19
Za <- ato.a[,9]
Za <- (Za - Zm)/sqrt(Zv)

## uses nc defined above
out.a <- mleHetGP(Xa, Za, lower=rep(0.01, ncol(X)), 
  upper=rep(30, ncol(X)), covtype="Matern5_2", known=list(beta0=0), 
  noiseControl=nc, maxit=100000, settings=list(return.matrices=FALSE))
}

##
## the following code duplicates a predictive comparison in
## the package vignette
##

## first using the model fit to the train partition (out)
out <- rebuild(out)

## predicting out-of-sample at the test sights
phet <- predict(out, Xtest)
phets2 <- phet$sd2 + phet$nugs
mhet <- as.numeric(t(matrix(rep(phet$mean, 10), ncol=10)))
s2het <- as.numeric(t(matrix(rep(phets2, 10), ncol=10)))
sehet <- (unlist(t(Ztest)) - mhet)^2
sc <- - sehet/s2het - log(s2het)
mean(sc)

## predicting at the held-out training replicates
phet.out <- predict(out, Xtrain.out)
phets2.out <- phet.out$sd2 + phet.out$nugs
s2het.out <- mhet.out <- Ztrain.out
for(i in 1:length(mhet.out)) {
  mhet.out[[i]] <- rep(phet.out$mean[i], length(mhet.out[[i]]))
  s2het.out[[i]] <- rep(phets2.out[i], length(s2het.out[[i]]))
}
mhet.out <- unlist(t(mhet.out))
s2het.out <- unlist(t(s2het.out))
sehet.out <- (unlist(t(Ztrain.out)) - mhet.out)^2
sc.out <- - sehet.out/s2het.out - log(s2het.out)
mean(sc.out)

\dontrun{
## then using the model trained from the "adaptive" 
## sequential design, with comparison from the "batch" 
## one above, using the scores function
out.a <- rebuild(out.a)
sc.a <- scores(out.a, Xtest = Xtest, Ztest = Ztest)
c(batch=mean(sc), adaptive=sc.a)

## an example of one iteration of sequential design

  Wijs <- Wij(out.a$X0, theta=out.a$theta, type=out.a$covtype)
  h <- horizon(out.a, Wijs=Wijs)
  control = list(tol_dist=1e-4, tol_diff=1e-4, multi.start=30, maxit=100)
  opt <- IMSPE_optim(out.a, h, Wijs=Wijs, control=control)
  opt$par
}
}

\seealso{ \code{\link{bfs}}, \code{\link{sirEval}}, \code{link{rebuild}}, 
  \code{\link{horizon}}, \code{\link{IMSPE_optim}}, \code{\link{mleHetGP}}, 
  \code{vignette("hetGP")} }

\keyword{datasets}
\keyword{datagen}
