
context("Graphical functions checks")

registerDoParallel(2)
#### Parameter setup
m<-45; M<-60
p<-.4;
t1<-1; t2<-2; k<-2

NmonteC<-1e2
S<-1e3*(1:5)
alpha<-.8; H<-0.8; sigma<-0.3

############################################################################

H_hat_f <- function(p,k,path) {hh<-H_hat(p,k,path); list(H=hh)}
theor_3_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
                             m=m,M=M,alpha=alpha,H=H,
                             sigma=sigma,H_hat_f,
                             p=p,k=k)
###################

l_plot_vb <- Plot_vb(theor_3_1_H_clt)

suppressWarnings({
    l_plot_dens<- Plot_dens(par_vec=c('H'), MC_data=theor_3_1_H_clt, Nnorm=1e7)
})

# Test if the objects returned are of the correct sizes
test_that("Graphics works with H-estimatior only", {
    expect_gt(length(l_plot_vb),expected=1)
    expect_equal(length(l_plot_dens), 1)
})

############################################################################

NmonteC<-5e1
alpha<-1.8; H<-0.8; sigma<-0.3


# How to plot empirical density

theor_3_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
                             m=m,M=M,alpha=alpha,H=H,
                             sigma=sigma,ContinEstim,
                             t1=t1,t2=t2,p=p,k=k)

###################

l_plot_vb <- Plot_vb(theor_3_1_H_clt)

suppressWarnings({
    l_plot_dens<- Plot_dens(par_vec=c('alpha', 'H'), MC_data=theor_3_1_H_clt, Nnorm=1e7)
})

# Test if the objects returned are of the correct sizes
test_that("Graphics works with ContinEstim estimatior", {
    expect_gt(length(l_plot_vb),expected=1)
    expect_equal(length(l_plot_dens), 2)
})
