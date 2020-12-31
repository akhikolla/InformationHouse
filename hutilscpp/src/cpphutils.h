#ifndef cpphutils_H
#define cpphutils_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DoubleVector do_range_dbl(NumericVector x, double halt_if_min = 1, double halt_if_max = -1);

// [[Rcpp::export]]
DoubleVector do_range_int(IntegerVector x, int halt_if_min = 1, int halt_if_max = -1);

bool single_ox_x1_x2(int x, int oix, int x1, int x2);
bool single_ox_x1_x2(double x, int oix, double x1, double x2);

bool do_in_int(int x, IntegerVector table);

bool do_is_safe2int(double x);
int type_safe2int(double x);
IntegerVector dblTable2int(DoubleVector table);

// switch(operator,
//        "!=" = 1L,
//        "==" = 2L,
//        ">=" = 3L,
//        "<=" = 4L,
//        ">"  = 5L,
//        "<"  = 6L,
//        "%in%" = 7L,
//        "%between%" = 8L,
//        "%(between)%" = 9L,
//        "%]between[%" = 10L,
//        0L)

// const int op = !(eq || gt || lt) ? 0 : (eq ? (gt ? 2 : (lt ? 3 : 1)) : (gt ? 4 : 5));
// != == >= <=  >  <
//  0  1  2  3  4  5
const int OP_NE = 1;
const int OP_EQ = 2;
const int OP_GE = 3;
const int OP_LE = 4;
const int OP_GT = 5;
const int OP_LT = 6;
const int OP_IN = 7;
const int OP_BW = 8;
const int OP_BO = 9;
const int OP_BC = 10;


// texParse
constexpr char SPACE = 32;
constexpr char OPEN_BRACE = 123;
constexpr char STOP_BRACE = 125;
constexpr char OPEN_SQBRK = 91;
constexpr char STOP_SQBRK = 93;



#endif
