// $Id: outlier.cpp 231 2020-07-26 20:56:26Z lao $
#include <math.h>
#include <Rcpp.h>
#include <vector>
#include "outlier.h"


// #define RCPP_DO_BOUNDS_CHECKS
using namespace Rcpp;
using namespace std;


/*
\prg\r\bin\R.exe --vanilla
library(Rcpp)
compileAttributes("Benchmarking", verbose=TRUE)
*/
// g++ -c -Ofast -Ic:\prg\R\library\Rcpp\include -Ic:\prg\R\include


// Determinant af trekantmatrix; produktet af diagonalelementerne
double determinant(NumericMatrix A)
{
    double det = 1.0;
    for(int i=0; i < A.nrow(); ++i) det *= A(i,i);
    det = det * det;
    return det;
}


#if 0
#define PR {for (int i=0; i < (int)del.size(); ++i) Rprintf(" %i", del[i]); Rprintf("\n");}
#else
#define PR {;}
#endif


// Outlier.ap hoveddel
// [[Rcpp::export]]
void outlierCpp(const int K, const int R, NumericMatrix xy,
        NumericMatrix ratio, NumericMatrix imat, NumericVector rmin)
{
    NumericMatrix XY = matProdT_LO(xy);
    NumericMatrix L = chol_LO(XY);
    double S = determinant(L);

    bool lav_ny_chol {false};
    int last = ratio.nrow();

    NumericMatrix Ainv(XY.ncol(), XY.ncol());
    double detA;
    double rr;
    NumericVector RX(last);
    double rrlast;
    int e, h;
    int nmmp1;
    vector<int> del(1);  // Mængde af indeks der skal udelades, antal firms
    vector<int> j(1);

    for (int r=1; r <= R; ++r) {
        if (r>1) {
            del.resize(r);
            L = chol_downdate(L, xy(r-1-1,_));
            Ainv = inverse_spd(L, true);
            detA = determinant(L);
        } else {
            Ainv = inverse_spd(L, true);
            detA = S;
        }
        RX.fill(DBL_MAX);
        rrlast = DBL_MAX;

        e = 0;
        h = r;
        del.resize(r);
        for (int i=0; i < r; ++i) del[i] = i+1;

        rr = det_downdate(Ainv, xy(del[r-1]-1,_), detA);
        PR;
        rr = rr / S;
        RX[0] = rr;
        rrlast = RX[last-1];
        nmmp1 = K - r + 1;

        while(del[0] < nmmp1) {
            if (e < K - h) {
                h = 1;
                e = del[r-1];
                j.resize(1);
                j[0] = 1;
            } else {
                e = del.at(r-h-1);
                ++h;
                j.resize(h);
                for (int i = 0; i < h; ++i)  j[i] = i+1;
                lav_ny_chol = true;
            }
            for ( int i=1; i <= h; ++i) del[r - h + j[i-1] - 1] = e + j[i-1];
            if (lav_ny_chol)  {
                NumericMatrix A(xy.ncol(), xy.ncol());
                double aa;
                for (int j=0; j < (int)XY.ncol(); ++j) {
                    for (int i=j; i < (int)XY.nrow(); ++i) {
                        aa = 0.;
                        for (int k=0; k < r-1; ++k) {
                            aa += xy.at(del[k]-1, i) * xy.at(del[k]-1, j);
                        }
                        A.at(i,j) = XY(i,j) - aa;
                        A(j,i) = A(i,j);
                    }
                }
                NumericMatrix Ln = chol_LO(A);
                Ainv = inverse_spd(Ln, true);
                detA = determinant(Ln);
                lav_ny_chol = false;
            }  // if (lav_ny_chol)
            rr = det_downdate(Ainv, xy(del[r-1]-1,_), detA);
            PR;
            rr = rr / S;
            if (rr < rrlast) {
                if (rr < min(RX)) {
                    for (int i=0; i < r; ++i) {
                        imat(r-1, i) = (double) del[i];
                    }
                }
                RX[last-1] = rr;
                RX = RX.sort();
                rrlast = RX[last-1];
            }
            rmin[r-1] = min(RX);
            for (int i=0; i < ratio.nrow(); ++i) {
                double Rratio = log(RX[i]/min(RX));
                ratio(i,r-1) = Rratio;
            }
        }  // while
    }  // for (r)
}  // outlierCpp
