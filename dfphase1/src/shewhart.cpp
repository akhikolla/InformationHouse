#include <Rcpp.h>
using namespace Rcpp;

namespace {
    inline double ggmedian(int n, double *x, double *w) {
	int half = n/2;
	std::copy(x, x+n, w);
	std::nth_element(w, w+half, w+n);
	if (n == 2*half) {
	    return (*std::max_element(w,w+half)+w[half])/2;
	} else {
	    return w[half];
	}
    } 

    struct Comparator
    {
	Comparator(const double *x) : data(x) {}
	bool operator()(int left, int right) const { return data[left] < data[right]; }
	const double *data;
    };

    inline void ggrank(int n, double *x, double *r, int *a) {
	int i, j, s;
	double d;
	for (i=0; i<n; i++) a[i] = i;
	std::sort(a,a+n,Comparator(x));
	for (i=0; i<n;) {
	    for (s=i, j=i+1; (j<n)&&(x[a[j]]<=x[a[i]]); j++) s += j;
	    for (d=1+((double)s)/(j-i); i<j; i++) r[a[i]] = d;
	}   
    }

    int isample(const int n) {return floor(n*unif_rand());}

    inline void ggcolmeans(NumericMatrix x, NumericVector xbar) {
	int n=x.nrow(), m=x.ncol(), i, j;
	double a, *xi=x.begin();
	for (i=0; i<m; i++, xi+=n) {
	    for (j=0, a=0.0; j<n; j++) a += xi[j];
	    xbar[i] = a/n;
	}
    }

    inline void ggcolxbars(NumericMatrix x, NumericVector xbar, NumericVector s) {
	int n=x.nrow(), m=x.ncol(), i, j;
	double a, b, c, *xi=x.begin(),
	    c4 = M_SQRT2*exp(R::lgammafn(0.5*n)-R::lgammafn(0.5*(n-1)));
	for (i=0; i<m; i++, xi+=n) {
	    for (j=0, a=b=0.0; j<n; j++) {
		a += c = xi[j];
		b += c*c;
	    }
	    xbar[i] = a/n;
	    s[i] = sqrt(b-a*a/n)/c4;
	}
    }	

    // length(w)=m if !aggr_with_mean
    inline void horsexbars(NumericMatrix x, bool aggr_with_mean,
			   NumericVector xbar, NumericVector s, NumericVector est,
			   NumericVector w) {
	ggcolxbars(x,xbar,s);
	if (aggr_with_mean) {
	    est[0] = mean(xbar);
	    est[1] = mean(s);
	} else {
	    est[0] = ggmedian(xbar.size(),xbar.begin(),w.begin()); 
	    est[1] = ggmedian(s.size(),s.begin(),w.begin());
	}
    }

    inline void horserank(NumericMatrix x, NumericVector l, NumericVector s,
			  NumericMatrix r, NumericVector v, IntegerVector a) {
	int i, j, n = x.nrow(), m = x.ncol(), N=m*n, MN=m*(n-1);
	double aq, bq, cq, mq, sq, med = ggmedian(N,x.begin(),r.begin()),
	    c4 = exp(R::lgammafn(0.5*MN)-R::lgammafn(0.5*(MN-1)))*sqrt(2.0/(MN-1.0));
	ggrank(N, x.begin(), r.begin(), a.begin());
	ggcolmeans(r, l);
	l = (l-(N+1)/2.0) / sqrt((m-1)*(N+1)/12.0);
	for (i=0; i<N ; i++) x[i] = fabs(x[i]-med);
	ggrank(N, x.begin(), r.begin(), a.begin());
	for (i=0, mq=sq=0.0; i<m; i++) {
	    for (j=0, aq=bq=0; j<n; j++) {
		aq += cq = r(j,i)*r(j,i);
		bq += cq*cq;
	    }
	    aq /= n;
	    s[i] = aq;
	    mq += aq ;
	    sq += (bq-n*aq*aq)/(n-1);
	}
	mq /= m;
	sq = sqrt(sq/m)/c4;
	s = sqrt(static_cast<double>(n))*(s-mq)/sq;
    }

}

// [[Rcpp::export]]
List ggxbars(NumericMatrix x, bool aggr_with_mean, int L) {
    int i, n=x.nrow(), m = x.ncol();
    double sn = sqrt(static_cast<double>(n));
    NumericVector xb(m), s(m), est(2), w(aggr_with_mean?0:m);
    NumericMatrix xx=clone(x), pstat(3,L);
    for (i=0; i<L; i++) {
	std::random_shuffle(xx.begin(),xx.end(),isample);
	horsexbars(xx, aggr_with_mean, xb, s, est, w);
	xb = abs(xb-est[0]);
	pstat(0,i) = sn*(*std::max_element(xb.begin(),xb.end()))/est[1];
	pstat(1,i) = -(*std::min_element(s.begin(),s.end()))/est[1];
	pstat(2,i) = (*std::max_element(s.begin(),s.end()))/est[1];
    }
    horsexbars(x, aggr_with_mean, xb, s, est, w);
    return List::create(_["Xbar"]=xb,_["S"]=s,_["center"]=est[0],
			_["scale"]=est[1],_["perm"]=pstat);
}


// [[Rcpp::export]]
NumericMatrix ggxbarsall(int n, int m, bool aggr_with_mean, int rep) {
    int i;
    double sn = sqrt(static_cast<double>(n));
    NumericVector xb(m), s(m), est(2), w(aggr_with_mean?0:m);
    NumericMatrix x(n,m), stat(3,rep);
    for (i=0; i<rep; i++) {
	std::generate(x.begin(),x.end(),norm_rand);
	horsexbars(x, aggr_with_mean, xb, s, est, w);
	xb = abs(xb-est[0]);
	stat(0,i) = sn*(*std::max_element(xb.begin(),xb.end()))/est[1];
	stat(1,i) = -(*std::min_element(s.begin(),s.end()))/est[1];
	stat(2,i) = (*std::max_element(s.begin(),s.end()))/est[1];
    }
    return stat;
}


// [[Rcpp::export]]
List ggrank(NumericMatrix x) {
    int n=x.nrow(), m=x.ncol(), N=m*n;
    NumericMatrix xx=clone(x), r(n,m);
    IntegerVector a(N);
    NumericVector l(m), s(m), v(m);
    horserank(xx, l, s, r, v, a);
    return List::create(_["lRank"]=l,_["sRank"]=s);
}

// [[Rcpp::export]]
NumericMatrix ggrankall(int n, int m, int rep) {
    NumericVector l(m), s(m), v(m);
    NumericMatrix x(n,m), r(n,m), stat(2,rep);
    IntegerVector a(m*n);
    while (rep--) {
	std::generate(x.begin(),x.end(),unif_rand);
	horserank(x, l, s, r, v, a);
	l = abs(l);
	s = abs(s);
	stat(0,rep) = *std::max_element(l.begin(),l.end());
	stat(1,rep) = *std::max_element(s.begin(),s.end());
    }
    return stat;
}
