#include <Rcpp.h>
#include <cmath>

#ifndef INCLUDED_INTERFACE_H
#define INCLUDED_INTERFACE_H

using namespace Rcpp;
using namespace std;

inline int random_integer(int exc_max)
{
	return (int)(unif_rand() * (double)exc_max) % exc_max;
}

Rcpp::IntegerVector sample_int(int n, int r);

RcppExport SEXP R_es(SEXP S, SEXP r);

RcppExport SEXP R_gset(
	SEXP N,
	SEXP S,
	SEXP r,
	SEXP min_its,
	SEXP max_its,
	SEXP signif,
	SEXP log_dismiss
);

#endif
