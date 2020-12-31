library(matrixcalc)

ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
    dnx <- dimnames(X)
    if (is.null(dnx)) 
        dnx <- vector("list", 2)
    s <- svd(X)
    nz <- s$d > tol * s$d[1]
    structure(if (any(nz)) 
        s$v[, nz] %*% (t(s$u[, nz])/s$d[nz])
        else X, dimnames = dnx[2:1])
}

#One iteration to run Newton Raphson to get optimal weights
Inv<-function(M,I) {
    check_sin<-is.singular.matrix(M)
    if(check_sin==FALSE) inv<-try(solve(M),TRUE) else inv<-try(solve(M+I),TRUE)
    if(length(inv)==1) {
        M<-.00001*M
        check_sin<-is.singular.matrix(M)
        if(check_sin==FALSE) inv<-solve(M) else inv<-solve(M+I)
        inv<-.00001*inv
    }
    inv
}

D_weight<-function(W,T,X,d,q) {
    M3<-upinfor(W,T[1,],X,3)
    M4<-upinfor(W,T[2,],X,4)
    M5<-upinfor(W,T[3,],X,5)
    I3<-10^-10*diag(3)
    I4<-10^-10*diag(4)
    I5<-10^-10*diag(5)
    inv3<-Inv(M3,I3)
    inv4<-Inv(M4,I4)
    inv5<-Inv(M5,I5)
    f1<-D_weight_1(q,W,T[1,],T[2,],T[3,],X,inv3,inv4,inv5)
    f2<-D_weight_2(q,W,T[1,],T[2,],T[3,],X,inv3,inv4,inv5)
    newweight<-W-d*(f1%*%ginv(f2))
    newweight
}

DD_weight<-function(W,T,X,d,I4,I5,order) {
    M4<-upinfor(W,T,X,order-1)
    M5<-upinfor(W,T,X,order)
    inv1<-Inv(M4,I4)
    inv<-Inv(M5,I5)
    f1<-DD_weight_1(W,T,X,inv,inv1,order)
    f2<-DD_weight_2(W,T,X,inv,inv1,order)
    newweight<-W-d*(f1%*%ginv(f2))
    newweight
}

c_weight<-function(W,T,X,d,p,order,UB,I) {
    M<-upinfor(W,T,X,order)
    if(exp(UB)<999) {
        inv<-ginv(M)
    } else {
        inv<-Inv(M,I)
    }
    f1<-c_weight_1(W,T,X,inv,p,order)
    f2<-c_weight_2(W,T,X,inv,p,order)
    newweight<-W-d*(f1%*%ginv(f2))
    newweight
}

#Newton Raphson method to get optimal weights for given design points "X"
#S_weight<-function(X,T,e1,q,p,order) {
S_weight<-function(X,T,e1,f,...){
    diff<-10
    W<-rep(1/length(X),length(X)-1)
    while(diff>e1) {
         d<-1
         #NW<-D_weight(W,T,X,d,q)
         #NW<-D_weight(W,T,X,d,...)
         #NW<-c_weight(W,T,X,d,p,order)
         #NW<-c_weight(W,T,X,d,...)
         NW<-f(W,T,X,d,...)
         minW<-min(min(NW),1-sum(NW))
         while(minW<0 & d>.0001) {
               d<-d/2
               #NW<-D_weight(W,T,X,d,q)
               #NW<-D_weight(W,T,X,d,...)
               #NW<-c_weight(W,T,X,d,p,order)
               #NW<-c_weight(W,T,X,d,...)
               NW<-f(W,T,X,d,...)
               minW<-min(min(NW),1-sum(NW))
         }
         NW<-c(NW,1-sum(NW))
         n<-length(NW)
         minW<-min(NW)
         diff<-max(abs(W-NW[1:n-1]))
         if (abs(minW)<.0000001||minW<0) {
               for(i in 1:n)
                   if (NW[i]==minW) NW[i]<-0         
         }
         D<-rbind(X,NW)
         for (i in 1:n) 
               if (D[2,i]==0) D[,i]<-NA
         X<-D[1,]
         W<-D[2,]
         X<-na.omit(X)
         W<-na.omit(W)
         W<-W[1:(length(X)-1)]
    }
    W<-c(W,1-sum(W))
    D<-rbind(X,W)
    D
}




