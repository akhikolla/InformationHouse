// $Id: outlier.h 225 2020-06-27 14:31:29Z lao $
#ifndef NORCPP
    #include <Rcpp.h>

    using namespace Rcpp;
#else
    #define Rcpp
    #include <vector>
#endif

using namespace std;


#ifdef NORCPP
namespace Rcpp {
class NumericMatrix {
    public:
        NumericMatrix() = default;
        NumericMatrix(int r, int c) {nrows=r; ncols=c; elem.resize(r*c);};
        NumericMatrix(NumericMatrix&&) = default;
        NumericMatrix& operator=(NumericMatrix&&) = default;
        ~NumericMatrix() = default;
        int nrow() { return nrows;};
        int ncol() { return ncols;};
        double at(int& i, int& j) {
            if (i<0||nrows<=i || j<0 || ncols<=j) throw out_of_range("NumericMatrix::at()");
            return elem[i + j*nrows];
        }
    private:
        int nrows;
        int ncols;
        vector<double> elem;
};

class NumericVector : std::vector<double> {
    public:
        NumericVector() = default;
        NumericVector(int s) {elem.resize(s);};
        ~NumericVector() = default;
        double at(int i) {
            if (i<0 || (int)elem.size()<=i) throw out_of_range("NumericVector::at()");
            return elem[i];
        }
    private:
        vector<double> elem;
};
};
#endif



Rcpp::NumericMatrix chol_LO(const Rcpp::NumericMatrix A);
Rcpp::NumericMatrix matProdT_LO(const Rcpp::NumericMatrix X);

Rcpp::NumericMatrix chol_downdate(const Rcpp::NumericMatrix L, const Rcpp::NumericVector v);
Rcpp::NumericMatrix inverse_spd(const Rcpp::NumericMatrix A, bool lower_triangel);
double det_downdate(const Rcpp::NumericMatrix A, const Rcpp::NumericVector v, const double det);

