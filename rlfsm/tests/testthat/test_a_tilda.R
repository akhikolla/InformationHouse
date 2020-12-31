
context("Work of a_tilde auxiliary functions")

#### Parameter setup
m<-25; M<-60; N<-2^12-M
alpha<-1.8; H<-0.8; sigma<-1.8
k<-2; p<-0.3; p_prime<-0.1
t1<-1; t2<-2
###################


test_that("a_tilda should produce a number", {
    expect_is(a_tilda(N,m,M,alpha,H), "numeric")
    expect_is(a_tilda_R(N,m,M,alpha,H), "numeric")
})

rr_C <- a_tilda(N,m,M,alpha,H)
rr_R <- a_tilda_R(N,m,M,alpha,H)

test_that("Diff versions of a_tilda should produce the same outcome", {
    expect_equivalent(rr_C, rr_R)
})

N=1; M=1; m=2
rr_C <- a_tilda(N,m,M,alpha,H)
rr_R <- a_tilda_R(N,m,M,alpha,H)
test_that("Diff versions of a_tilda should produce the same outcome", {
    expect_equivalent(rr_C, rr_R)
})

N=10; M=2; m=1
rr_C <- a_tilda(N,m,M,alpha,H)
rr_R <- a_tilda_R(N,m,M,alpha,H)
test_that("Diff versions of a_tilda should produce the same outcome", {
    expect_equivalent(rr_C, rr_R)
})
