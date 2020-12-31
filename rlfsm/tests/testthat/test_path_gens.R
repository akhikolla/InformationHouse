
context("Work of sample path generators and auxiliary functions")

#### Parameter setup
m<-25; M<-60; N<-2^12-M
alpha<-1.8; H<-0.8; sigma<-1.8
k<-2; p<-0.3; p_prime<-0.1
t1<-1; t2<-2
###################

List<-path(N,m,M,alpha,H,sigma,freq='L',disable_X=FALSE,levy_increments=NULL,seed=NULL)
X<-List$lfsm

###################
#### Check correctness of scaling with sigma ####
List01<-path(N=2^10-600,m=256,M=600,alpha=1.8,H=0.8,
             sigma=0.1,freq='L',disable_X=FALSE,seed=3)
lfsm01<-List01$lfsm

List10<-path(N=2^10-600,m=256,M=600,alpha=1.8,H=0.8,
             sigma=10,freq='L',disable_X=FALSE,seed=3)
lfsm10<-List10$lfsm

test_that("Check correctness of scaling with sigma", {
    expect_equal(round(lfsm01*100,8), round(lfsm10,8))
})

###################


library(stabledist)
Xpath<-path(N,m,M,alpha,H,sigma,freq='L',
            disable_X=FALSE,levy_increments=NULL,seed=2)$lfsm
set.seed(2)
Fast_path<-rlfsm:::path_fast(N,m,M,alpha,H,sigma,freq='L')

test_that("fast version of the path equals the regular one", {
    expect_equal(Xpath, Fast_path)
})

###################


test_that("all path generator should return either a list or a vector", {
    expect_is(List, 'list')
    expect_is(Fast_path, 'numeric')
})

###################


## Tests from MO18
# We test 2 different passes that we go through generating oLm and lfsm
seed<-3
### Checks levy increments
# Creating levy motion
levyIncrems<-path(N=N, m=m, M=M, alpha, H, sigma, freq='L',
                  disable_X=T, levy_increments=NULL, seed=seed)

# Creating lfsm based on the levy motion
lfsm_full<-path(m=m, M=M, alpha=alpha, H=H, sigma=sigma, freq='L',
                disable_X=F, levy_increments=levyIncrems$levy_increments,
                seed=seed)

sum(levyIncrems$levy_increments==lfsm_full$levy_increments)==length(lfsm_full$levy_increments)

test_that("Supplied Levy increments should be the same as in output", {
    expect_equal(levyIncrems$levy_increments, lfsm_full$levy_increments)
})


### Checks lfsm
lfsm_full2<-path(N=N, m=m, M=M, alpha, H, sigma, freq='L',
                 disable_X=F, levy_increments=NULL, seed=seed)
sum(lfsm_full$lfsm==lfsm_full2$lfsm)==length(lfsm_full$lfsm)


test_that("Seeded pregenerated increments supplied to path produce the same path as seeded generator", {
    expect_equal(lfsm_full$lfsm, lfsm_full2$lfsm)
})
###################


#### New parameter setup
m<-25; M<-60; N<-10
alpha<-1.8; H<-0.8; sigma<-1.8
k<-2; p<-0.3; p_prime<-0.1
t1<-1; t2<-2
###################


List<-path(N,m,M,alpha,H,sigma,freq='L',disable_X=FALSE,levy_increments=NULL,seed=NULL)
Fast_path<-rlfsm:::path_fast(N,m,M,alpha,H,sigma,freq='L')

test_that("Paths work with different M and N", {
    expect_is(List, 'list')
    expect_is(Fast_path, 'numeric')
})
###################

#### New parameter setup
m<-25; M<-60; N<-10
alpha<-0.8; H<-0.4; sigma<-0.2
k<-2; p<-0.3; p_prime<-0.1
t1<-1; t2<-2
###################


List<-path(N,m,M,alpha,H,sigma,freq='L',disable_X=FALSE,levy_increments=NULL,seed=NULL)
Fast_path<-rlfsm:::path_fast(N,m,M,alpha,H,sigma,freq='L')

test_that("Paths work with negative H-1/alpha", {
    expect_is(List, 'list')
    expect_is(Fast_path, 'numeric')
})
###################


List<-path(N,m,M,alpha,H,sigma,freq='H',disable_X=FALSE,levy_increments=NULL,seed=NULL)
Fast_path<-rlfsm:::path_fast(N,m,M,alpha,H,sigma,freq='H')

test_that("Paths work with freq='H'", {
    expect_is(List, 'list')
    expect_is(Fast_path, 'numeric')
})
