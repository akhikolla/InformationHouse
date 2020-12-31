#include "Interface.h"

using namespace Rcpp;
using namespace std;

IntegerVector sample_int(int n, int r) {
	IntegerVector result(r);
	LogicalVector still_in(n);
	for (int i = 0; i < n; i++)
		still_in[i] = true;

	for (int i = 0; i < r; i++) {
		do {
			result[i] = random_integer(n);
		}
		while (!still_in[result[i]]);
		still_in[result[i]] = false;
	}

	return result;
}

// [[Rcpp::export]]
double es_raw(IntegerVector S, NumericVector r) {
	double N_r = 0.0;
	int N = r.length();
	int n = S.length();
	for (int i = 0; i < n; i++) {
		N_r += r[S[i]];
	}
	double max_hit = 0.0;
	double p_hit = 0.0;

	for (int i = 0; i < n; i++) {
		p_hit += r[S[i]]/N_r;
		double p_miss = ((double)(S[i] - i))/((double)N-(double)n);
		if ((p_hit - p_miss) > max_hit) max_hit = p_hit - p_miss;
	}

	return max_hit;
}

// [[Rcpp::export]]
double gset_raw(
	int N,
	IntegerVector S,
	NumericVector r,
	int min_its,
	int max_its,
	double signif,
	double log_dismiss
) {
	int n = S.length();
	double max_hit = es_raw(S, r);
	int as_sim = 0;
	int samples = 0;

	do {
		IntegerVector samp = sample_int(N, n);
		sort(samp.begin(), samp.end());
		samples++;
		as_sim += (int)(es_raw(samp, r) >= max_hit);
	}
	while (samples < min_its || ((R::pnorm((double)as_sim, (double)samples*signif, sqrt((double)samples*signif*(1.0-signif)),false,true) > log_dismiss) && (samples < max_its)));

	return (double)(as_sim + 1) / (double)(samples + 1);
}
