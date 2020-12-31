
context("Statistic-related functions checks")

#### Parameter setup
m<-25; M<-60; N<-2^12-M
alpha<-1.8; H<-0.8; sigma<-1.8
k<-2; p<-0.3; p_prime<-0.1
t1<-1; t2<-2
###################

List<-path(N,m,M,alpha,H,sigma,freq='L',disable_X=FALSE,levy_increments=NULL,seed=NULL)
X<-List$lfsm



###################

test_that("estimators should return a numeric", {
    expect_is(m_pk(k,p,alpha,H,sigma), "numeric")
    expect_is(R_hl(p,k,X), "numeric")
    expect_is(H_hat(p,k,X), "numeric")
    expect_is(phi(t1,k,X,freq="L"), "numeric")
    expect_is(alpha_hat(t1,t2,k,X,freq="L"), "numeric")
    expect_is(alpha_hat(t1,t2,k,X,H=H,freq="H"), "numeric")
    expect_is(sigma_hat(t1,k,X,alpha,H,freq="H"), "numeric")
})


###################


Y<-paths(N_var=10,parallel=F,N=N,m=m,M=M,
         alpha=alpha,H=H,sigma=sigma,freq='L',
         disable_X=FALSE,levy_increments=NULL)
Rtr_st_low <- Retrieve_stats(paths=Y,true_val=sigma,Est=sigma_hat,t1=t1,k=2,alpha=alpha,H=H,freq="L")

test_that("Retrieve_stats in low frequency should return a list", {
    expect_is(Rtr_st_low, "list")
})

###################

Y<-paths(N_var=10,parallel=F,N=N,m=m,M=M,
         alpha=alpha,H=H,sigma=sigma,freq='L',
         disable_X=FALSE,levy_increments=NULL)
Rtr_st_H <- Retrieve_stats(paths=Y,true_val=sigma,Est=sigma_hat,t1=t1,k=2,alpha=alpha,H=H,freq="L")

test_that("Retrieve_stats in high frequency should return a list", {
    expect_is(Rtr_st_H, "list")
})


H_hat_bv    <- Retrieve_stats(paths=Y,true_val=H,Est=H_hat,p=0.3,k=2)
alph_hat_bv <- Retrieve_stats(paths=Y,true_val=alpha,Est=alpha_hat,t1=t1,t2=t2,k=2,H=H,freq='H')
si_hat_bv   <- Retrieve_stats(paths=Y,true_val=sigma,Est=sigma_hat,t1=t1,k=2,alpha=alpha,H=H,freq='H')

test_that("Retrieve_stats in high frequency should return a list", {
    expect_is(H_hat_bv, "list")
    expect_is(alph_hat_bv, "list")
    expect_is(si_hat_bv, "list")
})



###################

N_alpha <- Norm_alpha(h_kr,alpha=1.8,k=2,r=1,H=0.8,l=0)

test_that("Norm_alpha should return a value", {
    expect_is(N_alpha, "list")
    expect_is((unlist(N_alpha[1])), "numeric")
    expect_is((unlist(N_alpha[2])), "numeric")
})
###################

###################

###################
