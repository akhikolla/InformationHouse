// $Id: chol_down.cpp 225 2020-06-27 14:31:29Z lao $
#include <math.h>
#include <Rcpp.h>
#include "outlier.h"

using namespace Rcpp;


/*
\prg\r\bin\R.exe --vanilla
library(Rcpp)
compileAttributes("Benchmarking", verbose=TRUE)
*/
// g++ -c -Ofast -Ic:\prg\R\library\Rcpp\include -Ic:\prg\R\include



// Find Cholesky faktorisering af symmetrisk positiv definit matrix A;
// resultat er en nedre trekantsmatrix.
// Inputmatricen A ændres ikke.
// Særlig beregnet til brug i outlier.ap derfor ingen kontrol af input.
// Der er ingen test af noget så funktionen kan give underlige
// resultater hvis 'A' ikke er symetrisk og positiv definit.
// [[Rcpp::export]]
NumericMatrix chol_LO(const NumericMatrix A)
{
    NumericMatrix L(A.nrow(), A.ncol());
    L(0,0) = sqrt(A(0,0));
    for ( int i=0; i < A.nrow(); ++i)  {
        double s = 0;
        for (int k=0; k < i; ++k) {
            s += L(i,k)*L(i,k);
        }
        L(i,i) = sqrt(A(i,i) - s);
        for (int j=i+1; j < A.nrow(); ++j) {
            s = 0;
            for (int k=0; k < i; ++k) {
                s += L(j,k)*L(i,k);
            }
            L(j,i) = (A(j,i) - s)/L(i,i);
        } // for j
    }  // for i
    return(L);
}  // chol_LO



// Opdater Cholesky faktorisering ved en rang 1 ændring enten
// ved plus eller ved minus.
// Input L og v ændres IKKE undervejs og er derfor IKKE ændret ved retur.
// [[Rcpp::export]]
NumericMatrix chol_downdate(const NumericMatrix L, const NumericVector v)
{
    int m = v.size() - 1;
    // NumericMatrix A(m+1, m+1);
    // NumericVector u(m+1);
    NumericMatrix A = L + 0.0;
    NumericVector u = v + 0.0;
    for (int i=0; i < m; ++i)  {
        double r = sqrt(A(i,i)*A(i,i) - u[i]*u[i]);
        double s = u[i] / A(i,i);
        double t = r / A(i,i);
        A(i,i) = r;
        for (int j=i+1; j < m+1; ++j)  {
            A(j, i) = (A(j,i) - s * u[j]) / t;
            u[j] = t * u[j] - s * A(j,i);
        }
    }
    A(m,m) = sqrt(A(m,m)*A(m,m) - u[m]*u[m]);
    return A;
}  // chol_down...

// [[Rcpp::export]]
NumericMatrix chol_downdate2(const NumericMatrix L, const NumericVector v)
{
    int m = v.size();
    NumericMatrix A(m, m);
    NumericVector u(m);
    // A = L + 0.0;
    // u = v + 0.0;
    for (int i=0; i < m; ++i)  {
        double r = sqrt(L(i,i)*L(i,i) - v[i]*v[i]);
        double s = v[i] / L(i,i);
        double t = r / L(i,i);
        A(i,i) = r;
        for (int j=i+1; j < m+1; ++j)  {
            A(j, i) = (L(j,i) - s * v[j]) / t;
            u[j] = t * v[j] - s * A(j,i);
        }
    }
    A(m,m) = sqrt(L(m,m)*L(m,m) - u[m]*u[m]);
    return A;
}  // chol_down...




// Beregn determinant ud fra Cholesky faktorisering 'L' der
// dowddates med rang 1.
// Bemærk L bliver IKKE ændret i funktionen!!!
// [[Rcpp::export]]
double det_chol_downdate(const NumericMatrix L, const NumericVector v)
{
    NumericMatrix LL = L + 0;
    LL = chol_downdate(LL, v);
    // Beregn determinanten som produktet af diagonalelementerne i anden
    double det = 1;
    for (int i=0; i<v.size(); ++i)  {
        det *= LL(i,i);
    }
    det = det * det;
    return(det);
}  // chol_down...




// [[Rcpp::export]]
NumericMatrix chol_update(NumericMatrix L, NumericVector v)
{
    int m = v.size() - 1;
    double r, s, t;
    for (int i=0; i < m; ++i)  {
        r = sqrt(L(i,i)*L(i,i) + v(i)*v(i));
        s = v(i) / L(i,i);
        t = r / L(i,i);
        L(i,i) = r;
        for (int j=i+1; j<m+1; ++j)  {
            L(j,i) = (L(j,i) + s * v[j]) / t;
            v(j) = t * v(j) - s * L(j,i);
        }
    }
    L(m,m) = sqrt(L(m,m)*L(m,m) + v(m)*v(m));
    return L;

    // Beregn determinanten som produktet af diagonalelementerne i anden
    double det = 1;
    for (int i=0; i<m+1; ++i)  {
        det *= L(i,i);
    }
    det = det * det;
    return(det);
}  // chol_updateCC()



// Invers af en nedre trekantmatrix; input overskrives ikke
// [[Rcpp::export]]
NumericMatrix inverse_LO(const NumericMatrix L)
{
    int n = L.nrow();
    NumericMatrix A(L.nrow(), L.ncol());
    // std::fill(A.begin(), A.end(), 0);
    for (int i=0; i<n; ++i) {
        for (int j=0; j<i; ++j) {
            double s = 0;
            for (int k=0; k<i; ++k) {
               s += L(i,k) * A(k,j);
            }
            A(i,j) = - s / L(i,i);
        }
        A(i,i) = 1/L(i,i);
    }
    return A;
} // solve_LO



// Matrix produkt: x'x
// [[Rcpp::export]]
NumericMatrix matProdT_LO(const NumericMatrix X)
{
     NumericMatrix XX(X.ncol(), X.ncol());
     for (int i=0; i < X.ncol(); ++i) {
        for (int j=0; j < X.ncol(); ++j) {
            double s = 0;
            for (int h=0; h < X.nrow(); ++h) {
                s += X(h,i) * X(h,j);
            }
            XX(i,j) = s;
        }
    }
    return XX;
}



// Invers af symmetrisk positiv definit matrix
// [[Rcpp::export]]
NumericMatrix inverse_spd(const NumericMatrix A, bool lower_triangel=false)
{
    NumericMatrix L;
    if (lower_triangel)
        L = A;
    else
        L = chol_LO(A);  // Cholesky dekomposition
    NumericMatrix Linv = inverse_LO(L);
    NumericMatrix Ainv = matProdT_LO(Linv);
    return Ainv;
} // solve_LO



// Løser ligningssystemet Ly=d hver L er nedre trekantmatrix
// [[Rcpp::export]]
NumericVector solve_LO(const NumericMatrix L, const NumericVector d)
{
    int n = d.size();
    NumericVector y(L.ncol());
    y[0] = d[0]/L(0,0);
    for (int i=1; i<n; ++i) {
        double s = 0;
        for (int j=0; j<i; ++j) {
            s += L(i,j) * y[j];
        }
        y[i] = (d[i] - s)/L(i,i);
    }
    return(y);
}  // solve_LO



// Downdate determinant via Matrix Determinant Lemma hvor
// A er invers matrix
// [[Rcpp::export]]
double det_downdate(const NumericMatrix A, const NumericVector v, const double det)
{
    int n = v.size();
    double s = 0;
    for(int i=0; i<n; ++i)  {
        for (int j=0; j<n; ++j)  {
            s += v[j] * A(j,i) * v[i];
        }
    }
    double res = (1 - s) * det;
    return res;
}
