#include <Rcpp.h>

using namespace Rcpp;
//' Summation of diagonal elements in a matrix
//' @param M A numeric matrix
//' @export
// [[Rcpp::export]]

double SDM(NumericMatrix M)
{
    double sum=0;
    for(int i=0;i<M.nrow();i++)
        sum=sum+M(i,i);

    return sum;
}
//' Transpose of a matrix
//' @param M A numeric matrix
//' @export
// [[Rcpp::export]]

NumericMatrix Trans(NumericMatrix M)
{
    int nrow=M.ncol(),ncol=M.nrow();
    NumericMatrix N(nrow,ncol);

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++)
        N(i,j)=M(j,i);

    return N;
}
//' Matrix subtraction
//' @param M1 A numeric matrix
//' @param M2 A numeric matrix
//' @export
// [[Rcpp::export]]

NumericMatrix Minus(NumericMatrix M1, NumericMatrix M2)
{
    int nrow=M1.nrow(),ncol=M1.ncol();
    NumericMatrix M3(nrow,ncol);

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++)
            M3(i,j)=M1(i,j)-M2(i,j);

    return M3;
}
//' Matrix addition
//' @param M1 A numeric matrix
//' @param M2 A numeric matrix
//' @export
// [[Rcpp::export]]

NumericMatrix Plus(NumericMatrix M1, NumericMatrix M2)
{
    int nrow=M1.nrow(),ncol=M1.ncol();
    NumericMatrix M3(nrow,ncol);

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++)
            M3(i,j)=M1(i,j)+M2(i,j);

    return M3;
}
//' Matrix multiplication
//' @param M1 A numeric matrix
//' @param M2 A numeric matrix
//' @export
// [[Rcpp::export]]

NumericMatrix Multiple(NumericMatrix M1, NumericMatrix M2)
{
    int nrow=M1.nrow(),ncol=M2.ncol(),P=M1.ncol();
    NumericMatrix M3(nrow,ncol);
    double sum;

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++) {
            sum=0;
            for(int k=0;k<P;k++)
                sum=sum+M1(i,k)*M2(k,j);
            M3(i,j)=sum;
        }

    return M3;
}
//' Multiply a constant to a matrix
//' @param s Numeric
//' @param M A numeric matrix
//' @export
// [[Rcpp::export]]

NumericMatrix sMultiple(double s, NumericMatrix M)
{
    int nrow=M.nrow(),ncol=M.ncol();
    NumericMatrix N(nrow,ncol);

    for(int i=0;i<nrow;i++)
        for(int j=0;j<ncol;j++)
            N(i,j)=s*M(i,j);

    return N;
}
//' Obtain a information matrix at a single design point
//' @param T A numeric vector. Model parameter values
//' @param x numeric. A single design point(dose level)
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

NumericMatrix infor(NumericVector T, double x, int order)
{
    NumericMatrix M(order,1);

    if(order==3) {
        M(0,0)=1/(1+exp(T[1]*(T[2]-x)));
        M(1,0)=-T[0]*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),2)*(T[2]-x);
        M(2,0)=-T[1]*T[0]*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),2);
    }
    if(order==4) {
        M(0,0)=1/pow(1+exp(T[1]*(T[2]-x)),T[4]);
        M(1,0)=-(T[0]-T[3])*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),1+T[4])*T[4]*(T[2]-x);
        M(2,0)=-T[1]*(T[0]-T[3])*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),1+T[4])*T[4];
        M(3,0)=1-1/pow(1+exp(T[1]*(T[2]-x)),T[4]);
    }
    if(order==5) {
        M(0,0)=1/pow(1+exp(T[1]*(T[2]-x)),T[4]);
        M(1,0)=-(T[0]-T[3])*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),1+T[4])*T[4]*(T[2]-x);
        M(2,0)=-T[1]*(T[0]-T[3])*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),1+T[4])*T[4];
        M(3,0)=1-1/pow(1+exp(T[1]*(T[2]-x)),T[4]);
        M(4,0)=-(T[0]-T[3])/pow(1+exp(T[1]*(T[2]-x)),T[4])*log(1+exp(T[1]*(T[2]-x)));
    }

    return Multiple(M,Trans(M));
}
//' Obtain normalized Fisher information matrix for the 3,4,5PL models
//' @param W A numeric vector. The first K-1 weights for a given design
//' @param T A numeric vector. Model parameter values
//' @param X A numeric vector. K design points for a given design
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

NumericMatrix upinfor(NumericVector W, NumericVector T, NumericVector X, int order)
{
    NumericMatrix last_infor=infor(T,X[X.size()-1],order);
    NumericMatrix infor_new=sMultiple(1-sum(W),last_infor);

    for(int i=0;i<X.size()-1;i++)
        infor_new=Plus(infor_new,sMultiple(W[i],infor(T,X[i],order)));

    return infor_new;
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x A numeric
//' @param xl A numeric
//' @param inv A numeric matrix.
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

double smalld1(NumericVector T, double x, double xl, NumericMatrix inv, int order)
{
    NumericMatrix M=Multiple(inv,Minus(infor(T,x,order),infor(T,xl,order)));
    double s=SDM(M);

    return s;
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x1 A numeric
//' @param x2 A numeric
//' @param xl A numeric
//' @param inv A numeric matrix
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

double smalldd1(NumericVector T, double x1, double x2, double xl, NumericMatrix inv, int order)
{
    NumericMatrix M=Multiple(Multiple(Multiple(inv,Minus(infor(T,x2,order),infor(T,xl,order))),inv),Minus(infor(T,x1,order),infor(T,xl,order)));
    double s=SDM(M);

    return -s;
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x A numeric
//' @param xl A numeric
//' @param inv A numeric matrix
//' @param inv1 A numeric matrix
//' @param order numeric. The number of model paraemters
// [[Rcpp::export]]

double d11(NumericVector T, double x, double xl, NumericMatrix inv, NumericMatrix inv1, int order)
{
    NumericMatrix M1=Multiple(inv,Minus(infor(T,x,order),infor(T,xl,order)));
    NumericMatrix M2=Multiple(inv1,Minus(infor(T,x,order-1),infor(T,xl,order-1)));

    return SDM(M1)-SDM(M2);
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x1 A numeric
//' @param x2 A numeric
//' @param xl A numeric
//' @param inv A numeric matrix
//' @param inv1 A numeric matrix
//' @param order numeric. The number of model paraemters
// [[Rcpp::export]]

double dd11(NumericVector T, double x1, double x2, double xl, NumericMatrix inv, NumericMatrix inv1, int order)
{
    NumericMatrix M1=Multiple(Multiple(Multiple(inv,Minus(infor(T,x2,order),infor(T,xl,order))),inv),Minus(infor(T,x1,order),infor(T,xl,order)));
    NumericMatrix M2=Multiple(Multiple(Multiple(inv1,Minus(infor(T,x2,order-1),infor(T,xl,order-1))),inv1),Minus(infor(T,x1,order-1),infor(T,xl,order-1)));

    return SDM(M2)-SDM(M1);
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x A numeric
//' @param order numeric. The number of model paraemters
// [[Rcpp::export]]

NumericMatrix f(NumericVector T, double x, int order)
{
    NumericMatrix M(order,1);

    if(order==3) {
        M(0,0)=1/(1+exp(T[1]*(T[2]-x)));
        M(1,0)=-T[0]*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),2)*(T[2]-x);
        M(2,0)=-T[1]*T[0]*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),2);
    }
    if(order==4) {
        M(0,0)=1/pow(1+exp(T[1]*(T[2]-x)),T[4]);
        M(1,0)=-(T[0]-T[3])*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),1+T[4])*T[4]*(T[2]-x);
        M(2,0)=-T[1]*(T[0]-T[3])*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),1+T[4])*T[4];
        M(3,0)=1-1/pow(1+exp(T[1]*(T[2]-x)),T[4]);
    }
    if(order==5) {
        M(0,0)=1/pow(1+exp(T[1]*(T[2]-x)),T[4]);
        M(1,0)=-(T[0]-T[3])*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),1+T[4])*T[4]*(T[2]-x);
        M(2,0)=-T[1]*(T[0]-T[3])*exp(T[1]*(T[2]-x))/pow(1+exp(T[1]*(T[2]-x)),1+T[4])*T[4];
        M(3,0)=1-1/pow(1+exp(T[1]*(T[2]-x)),T[4]);
        M(4,0)=-(T[0]-T[3])/pow(1+exp(T[1]*(T[2]-x)),T[4])*log(1+exp(T[1]*(T[2]-x)));
    }
    return M;
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x A numeric
//' @param inv A numeric matrix
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

double smallds1(NumericVector T, double x, NumericMatrix inv, int order)
{
    return Multiple(Multiple(Trans(f(T,x,order)),inv),f(T,x,order))(0,0)/order;
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x A numeric
//' @param inv A numeric matrix
//' @param inv1 A numeric matrix
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

double ds11(NumericVector T, double x, NumericMatrix inv, NumericMatrix inv1, int order)
{
    double s=Multiple(Multiple(Trans(f(T,x,order)),inv),f(T,x,order))(0,0);
    double t=Multiple(Multiple(Trans(f(T,x,order-1)),inv1),f(T,x,order-1))(0,0);

    return s-t;
}
//' EDp
//' @param T A numeric vector. Model parameter values
//' @param p A numeric
//' @export
// [[Rcpp::export]]

double Dp(NumericVector T, double p)
{
    double r;
    r=T[2]-1/T[1]*log(pow(1/p,1/T[4])-1);

    return r;
}
//' Partial derivative of the EDp w.r.t the model parameters
//' @param T A numeric vector. Model parameter values
//' @param p A numeric
//' @export
// [[Rcpp::export]]

NumericMatrix g(NumericVector T, double p)
{
    NumericMatrix M(T.size(),1);

    M[0]=0;
    M[1]=pow(T[1],-2)*log(pow(1/p,1/T[4])-1);
    M[2]=1;
    M[3]=0;
    M[4]=1/(T[1]*pow(T[4],2)*(pow(1/p,1/T[4])-1))*pow(1/p,1/T[4])*log(1/p);

    return M;
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x A numeric
//' @param xl A numeric
//' @param inv A numeric matrix
//' @param p A numeric
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

double D1(NumericVector T, double x, double xl, NumericMatrix inv, double p, int order)
{
    double s=Multiple(Multiple(Multiple(Multiple(Trans(g(T,p)),inv),Minus(infor(T,x,order),infor(T,xl,order))),inv),g(T,p))(0,0);

    return -s;
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x1 A numeric
//' @param x2 A numeric
//' @param xl A numeric
//' @param inv A numeric matrix
//' @param p A numeric
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

double DD1(NumericVector T, double x1, double x2, double xl, NumericMatrix inv, double p, int order)
{
    double s=Multiple(Multiple(Multiple(Multiple(Multiple(Multiple(Trans(g(T,p)),inv),Minus(infor(T,x2,order),infor(T,xl,order))),inv),Minus(infor(T,x1,order),infor(T,xl,order))),inv),g(T,p))(0,0);
    double t=Multiple(Multiple(Multiple(Multiple(Multiple(Multiple(Trans(g(T,p)),inv),Minus(infor(T,x1,order),infor(T,xl,order))),inv),Minus(infor(T,x2,order),infor(T,xl,order))),inv),g(T,p))(0,0);

    return s+t;
}
//' Matrix calculation
//' @param T A numeric vector. Model parameter values
//' @param x A numeric
//' @param inv A numeric matrix
//' @param p A numeric
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

double DS1(NumericVector T,double x,NumericMatrix inv,double p, int order)
{
    double s=Multiple(Multiple(Trans(f(T,x,order)),inv),g(T,p))(0,0);
    double t=Multiple(Multiple(Trans(g(T,p)),inv),f(T,x,order))(0,0);
    double k=pow(Multiple(Multiple(Trans(g(T,p)),inv),g(T,p))(0,0),-1);

    return s*t*k;
}
//' Matrix calculation
//' @param q A numeric vector
//' @param W A numeric vector
//' @param T1 A numeric vector
//' @param T2 A numeric vector
//' @param T3 A numeric vector
//' @param X A numeric vector
//' @param inv1 A numeric matrix
//' @param inv2 A numeric matrix
//' @param inv3 A numeric matrix
//' @export
// [[Rcpp::export]]

NumericVector D_weight_1(NumericVector q, NumericVector W, NumericVector T1, NumericVector T2, NumericVector T3, NumericVector X, NumericMatrix inv1, NumericMatrix inv2, NumericMatrix inv3)
{
    int p=W.size(),k=X.size();
    NumericVector f1(p);

    for(int i=0;i<p;i++)
        f1[i]=q[0]*smalld1(T1,X[i],X[k-1],inv1,3)/3+q[1]*smalld1(T2,X[i],X[k-1],inv2,4)/4+q[2]*smalld1(T3,X[i],X[k-1],inv3,5)/5;

    return f1;
}
//' Matrix calculation
//' @param q A numeric vector
//' @param W A numeric vector
//' @param T1 A numeric vector
//' @param T2 A numeric vector
//' @param T3 A numeric vector
//' @param X A numeric vector
//' @param inv1 A numeric matrix
//' @param inv2 A numeric matrix
//' @param inv3 A numeric matrix
//' @export
// [[Rcpp::export]]

NumericVector D_weight_2(NumericVector q, NumericVector W, NumericVector T1, NumericVector T2, NumericVector T3, NumericVector X, NumericMatrix inv1, NumericMatrix inv2, NumericMatrix inv3)
{
    int p=W.size(),k=X.size();
    NumericMatrix f2(p,p);

    for(int i=0;i<p;i++)
        for(int j=0;j<p;j++)
            f2(i,j)=q[0]*smalldd1(T1,X[i],X[j],X[k-1],inv1,3)/3+q[1]*smalldd1(T2,X[i],X[j],X[k-1],inv2,4)/4+q[2]*smalldd1(T3,X[i],X[j],X[k-1],inv3,5)/5;

    return f2;
}
//' Matrix calculation
//' @param W A numeric vector
//' @param T A numeric vector
//' @param X A numeric vector
//' @param inv A numeric matrix
//' @param p A numeric
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

NumericVector c_weight_1(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, double p, int order)
{
    int p2=W.size(),k=X.size();
    NumericVector f1(p2);

    for(int i=0;i<p2;i++)
        f1[i]=D1(T,X[i],X[k-1],inv,p,order);

    return f1;
}
//' Matrix calculation
//' @param W A numeric vector
//' @param T A numeric vector
//' @param X A numeric vector
//' @param inv A numeric matrix
//' @param p A numeric
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

NumericMatrix c_weight_2(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, double p, int order)
{
    int p2=W.size(),k=X.size();
    NumericMatrix f2(p2,p2);

    for(int i=0;i<p2;i++)
        for(int j=0;j<p2;j++)
            f2(i,j)=DD1(T,X[i],X[j],X[k-1],inv,p,order);

    return f2;
}
//' Matrix calculation
//' @param W A numeric vector
//' @param T A numeric vector
//' @param X A numeric vector
//' @param inv A numeric matrix
//' @param inv1 A numeric matrix
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

NumericVector DD_weight_1(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, NumericMatrix inv1, int order)
{
    int p2=W.size(),k=X.size();
    NumericVector f1(p2);

    for(int i=0;i<p2;i++)
        f1[i]=d11(T,X[i],X[k-1],inv,inv1,order);

    return f1;
}
//' Matrix calculation
//' @param W A numeric vector
//' @param T A numeric vector
//' @param X A numeric vector
//' @param inv A numeric matrix
//' @param inv1 A numeric matrix
//' @param order numeric. The number of model paraemters
//' @export
// [[Rcpp::export]]

NumericMatrix DD_weight_2(NumericVector W, NumericVector T, NumericVector X, NumericMatrix inv, NumericMatrix inv1, int order)
{
    int p2=W.size(),k=X.size();
    NumericMatrix f2(p2,p2);

    for(int i=0;i<p2;i++)
        for(int j=0;j<p2;j++)
            f2(i,j)=dd11(T,X[i],X[j],X[k-1],inv,inv1,order);

    return f2;
}
