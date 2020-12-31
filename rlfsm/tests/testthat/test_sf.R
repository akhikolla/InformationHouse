
context("sf function checks")

#### Parameter setup
m<-25; M<-60; N<-2^12-M
alpha<-1.8; H<-0.8; sigma<-1.8
k<-2; p<-0.3; p_prime<-0.1
t1<-1; t2<-2
###################

List<-path(N,m,M,alpha,H,sigma,freq='L',disable_X=FALSE,levy_increments=NULL,seed=NULL)
X<-List$lfsm

X_1=c(1,4,3,6,8,5,3,5,8,5,1,8,6)

###################
test_that("sf should return a numeric", {
    expect_is(sf(path=X,f=abs,k=2,r=1,H=0.7,freq='H'), "numeric")
    expect_is(sf(path=X,f=cos,k=2,r=1,H=0.7,freq='H'), "numeric")
})


###################
sf_old<-function(path,f,k,r,H,freq,...){

    n<-length(path)-1 # -1 because the scalling factor is 1/(n-1)
    vector<-seq(r*k,n,by=1) # is passed as i to increment
    v1<-sapply(vector,FUN=increment,r=r,k=k,path=path)

    if(freq=='L') v2<-sapply(v1,FUN=f,...) else{
        if(freq=='H') v2<-sapply((n^H)*v1,FUN=f,...) else{
            stop('Parameter freq could take either "H" or "L" values')
        }
    }
    sum(v2)/(length(v2))
}


test_that("Different sfs should return the same", {
    expect_equal(sf(path=X_1,f=cos,k=2,r=1,H=NA,freq='L'),
                 sf_old(path=X_1,f=cos,k=2,r=1,H=NA,freq='L'))
    expect_equal(sf(path=X_1,f=cos,k=3,r=2,H=NA,freq='L'),
                 sf_old(path=X_1,f=cos,k=3,r=2,H=NA,freq='L'))
    expect_equal(sf(path=X_1,f=cos,k=3,r=2,H=0.8,freq='H'),
                 sf_old(path=X_1,f=cos,k=3,r=2,H=0.8,freq='H'))
})
