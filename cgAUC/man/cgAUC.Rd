\name{cgAUC}
\alias{cgAUC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate AUC when gold standard is continuous with large variables.}
\description{
The cgAUC can calculate the AUC-type measure of Obuchowski(2006) when gold standard is continuous, and find the optimal linear combination of variables with respect to this measure.}
\usage{
cgAUC(x, z, h, delta = 1, auto = FALSE, tau = 1, scale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The potential variables. It is a matrix with column of values of a variables. It should be standardized in this application.}
  \item{z}{
The gold standard variable. It should be standardized.}
  \item{h}{
The parameter controls the window width of smoothing function.}
  \item{delta}{
The parameter be used in TGDM. The default value is one.}
  \item{auto}{
Find the optimal delta in TGDN using cross-validation. If the auto is TRUE. The default is FALSE.}
  \item{tau}{
The parameter used in TGDM. The default value is one.}
  \item{scale}{
Scaling data when scale = 1, no scaling data when scale = 0. The default value is 1.}
}
\details{
In this package, we use the TGDM to find the optimal linear combination of variables in order to maximize the AUC-type measure. Before using this function, all of variables, including gold standard variable, should be standardized first. Below are parameters used in the algorithm:}
\value{
\item{Rev}{When Rev = 0 means l * 1; otherwise, l * -1.}
\item{l}{The estimate of coefficients for the optimal linear combination of variables.}
\item{theta.sh.h.p}{The estimate of the theta of Chang(2012) for the optimal linear combination of variables.}
\item{theta.sh.h.p.var}{The estimate of variance for the theta of Chang(2012).}
\item{cntin.ri}{The estimate of the theta of Chang(2012) for each single vaiable.}
\item{theta.h.p}{The estimate of the theta of Obuchowski(2006) for the optimal linear combination of variables.}
\item{theta.h.p.var}{The estimate of variance for the theta of Obuchowski(2006).}
\item{dscrt.ri}{The estimate of the theta of Obuchowski(2006) for each single vaiable.}
\item{delta}{The value of delta.}

}
\references{
Chang, YCI. Maximizing an ROC type measure via linear combination of markers when the gold reference is continuous. Statistics in Medicine 2012.\cr
Obuchowski NA. An ROC-type measure of diagnostic accuracy when the gold standard is continuous-scale. Statistics in Medicine 2006; 25:481--493.\cr
Obuchowski N. Estimating and comparing diagnostic tests accuracy when the gold standard is not binary. Statistics in Medicine 2005; 20:3261--3278.\cr
Friedman JH, Popescu BE. Gradient directed regularization for linear regression and classification. Technical Report, Department of Statistics, Stanford University, 2004.\cr
}
\author{
Yu-chia Chang}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

# n = 100; p = 5;
# r.x = matrix(rnorm(n * p), , p) # raw data
# r.z = r.x[ ,1] + rnorm(n) # gold standard
# x = scale(r.x) # standardized of raw data
# z = scale(r.z) # standardized of gold standard
# h = n^(-1 / 2)
# t1 = cgAUC(r.x, r.z, h, delta = 1, auto = FALSE, tau = 1, scale = 1) # the delta be constant
# t1
# t2 = cgAUC(r.x, r.z, h, delta = 1, auto = TRUE, tau = 1, scale = 1) # the delta be variable
# t2

## The function is currently defined as
function (x, z, h, delta = 1, auto = FALSE, tau = 1)
{
	x = scale(x)
	z = scale(z)
	conv = FALSE
	n = dim(x)[1]
	p = dim(x)[2]
	cntin.ri = dscrt.ri = rep(0, p)
	id = diag(p)
	for (i in 1:p) {
		dscrt.ri[i] = dscrt(x, z, id[i, ])$theta.h.p
		cntin.ri[i] = cntin(x, z, id[i, ], h)$theta.sh.h.p
	}
	beta.i = ifelse(cntin.ri > 0.5, 1, -1)
	dscrt.ri = ifelse(dscrt.ri > 0.5, dscrt.ri, (1 - dscrt.ri))
	cntin.ri = ifelse(cntin.ri > 0.5, cntin.ri, (1 - cntin.ri))
	y = x * matrix(beta.i, n, p, byrow = TRUE)
	max.x = which(cntin.ri == max(cntin.ri))
	theta.sh.h.p = 0
	l = id[max.x, ]
	while (conv == FALSE) {
		d.l = d.theta.sh.h.p(y, z, l, h)
		max.d.l = max(d.l)
		ind.d.l = ifelse(d.l >= (tau * max.d.l), 1, 0) * d.l
		if (auto == TRUE) {
			delta = optimal.delta(y, z, l, h, ind.d.l)
		}
		l = l + delta * ind.d.l
		l = l/max(l)
		theta.temp = cntin(y, z, l, h)$theta.sh.h.p
		ifelse(abs(theta.temp - theta.sh.h.p) < 1e-04, conv <- TRUE, conv <- FALSE)
		theta.sh.h.p = theta.temp
	}
	optimal.dscrt = dscrt(y, z, l)
	theta.sh.h.p.var = cntin(y, z, l, h)$var
	l = l * beta.i
	return(list(l = l, theta.sh.h.p = theta.sh.h.p, theta.sh.h.p.var = theta.sh.h.p.var,
				cntin.ri = cntin.ri, theta.h.p = optimal.dscrt$theta.h.p,
				theta.h.p.var = optimal.dscrt$var, dscrt.ri = dscrt.ri,
				delta = delta))
}
## The function is currently defined as
function (x, z, h, delta = 1, auto = FALSE, tau = 1) 
{
    x = scale(x)
    z = scale(z)
    conv = FALSE
    n = dim(x)[1]
    p = dim(x)[2]
    cntin.ri = dscrt.ri = rep(0, p)
    id = diag(p)
    for (i in 1:p) {
        dscrt.ri[i] = dscrt(x, z, id[i, ])$theta.h.p
        cntin.ri[i] = cntin(x, z, id[i, ], h)$theta.sh.h.p
    }
    beta.i = ifelse(cntin.ri > 0.5, 1, -1)
    dscrt.ri = ifelse(dscrt.ri > 0.5, dscrt.ri, (1 - dscrt.ri))
    cntin.ri = ifelse(cntin.ri > 0.5, cntin.ri, (1 - cntin.ri))
    y = x * matrix(beta.i, n, p, byrow = TRUE)
    max.x = which(cntin.ri == max(cntin.ri))
    theta.sh.h.p = 0
    l = id[max.x, ]
    while (conv == FALSE) {
        d.l = d.theta.sh.h.p(y, z, l, h)
        max.d.l = max(d.l)
        ind.d.l = ifelse(d.l >= (tau * max.d.l), 1, 0) * d.l
        if (auto == TRUE) {
            delta = optimal.delta(y, z, l, h, ind.d.l)
        }
        l = l + delta * ind.d.l
        l = l/max(l)
        theta.temp = cntin(y, z, l, h)$theta.sh.h.p
        ifelse(abs(theta.temp - theta.sh.h.p) < 1e-04, conv <- TRUE, 
            conv <- FALSE)
        theta.sh.h.p = theta.temp
    }
    optimal.dscrt = dscrt(y, z, l)
    theta.sh.h.p.var = cntin(y, z, l, h)$var
    l = l * beta.i
    return(list(l = l, theta.sh.h.p = theta.sh.h.p, theta.sh.h.p.var = theta.sh.h.p.var, 
        cntin.ri = cntin.ri, theta.h.p = optimal.dscrt$theta.h.p, 
        theta.h.p.var = optimal.dscrt$var, dscrt.ri = dscrt.ri, 
        delta = delta))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cgAUC }
