
context("Increment functions checks")

#### Parameter setup
m<-25; M<-60; N<-2^12-M
alpha<-1.8; H<-0.8; sigma<-1.8
k<-2; p<-0.3; p_prime<-0.1
t1<-1; t2<-2
###################

List<-path(N,m,M,alpha,H,sigma,freq='L',disable_X=FALSE,levy_increments=NULL,seed=NULL)
X<-List$lfsm


###################
test_that("increment(s) should return a numeric", {
    expect_is(increment(r=2, i=k+5, k, X), "numeric")
    expect_is(increments(2, k, X), "numeric")
})

test_that("increment(s) should return an error here", {
    expect_error(increment(r=2, i=4, k=length(X), X))
    expect_error(increments(2, k, NULL))
})


###################
test_that("Computation of increments is consistent", {
    expect_equal(X[2]-X[1], increment(r=1, i=1, k=1, X))
    expect_equal(increments(r=3, k=1, X)[1], increment(r=3, i=3, k=1, X))
    expect_equal(increments(r=3, k=1, X), increment(r=3, i=seq(1*3, (length(X)-1)), k=1, X))
})


X_1=c(1,4,3,6,8,5,3,5,8,5,1,8,6)
r=1; k=1
n <- length(X_1) - 1
DeltaX = increment(seq(r*k, n), path = X_1, k = k, r = r)
test_that("Computation of increments is consistent", {
    expect_equal(DeltaX, increments(k=k,r=r,X_1))
})


r=2; k=1
DeltaX = increment(seq(r*k, n), path = X_1, k = k, r = r)
sum(DeltaX == increments(k=k,r=r,X_1)) == length(DeltaX)
test_that("Computation of increments is consistent", {
    expect_equal(DeltaX, increments(k=k,r=r,X_1))
})

r=2; k=2
DeltaX = increment(seq(r*k, n), path = X_1, k = k, r = r)
test_that("Computation of increments is consistent", {
    expect_equal(DeltaX, increments(k=k,r=r,X_1))
})


###################



###################



###################


