/*
 
  This file is part of BWStest.
  
  BWStest is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  BWStest is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with BWStest.  If not, see <http://www.gnu.org/licenses/>.

  Baumgartner-Weiss-Schindler 2-sample test of equal distributions.
 
  see also: 

  + W. Baumgartner, P. Weiss, H. Schindler, 'A Nonparametric Test for the General Two-Sample Problem', 
    Biometrics, Vol. 54, No. 3 (Sep., 1998), pp. 1129-1135. http://doai.io/10.2307/2533862
  + M. Neuhauser, 'Exact Tests Based on the Baumgartner-Weiss-Schindler Statistic--A Survey', 
    Statistical Papers, Vol 46 (2005), pp. 1-30. http://doai.io/10.1007/BF02762032
  + M. Neuhauser, 'One-Sided Two-Sample and Trend Tests Based on a Modified Baumgartner-Weiss-Schindler 
    Statistic', J. Nonparametric Statistics, Vol 13 (2001) pp 729-739. http://doai.io/10.1080/10485250108832874
  + H. Murakami, 'K-Sample Rank Test Based on Modified Baumgartner Statistic and its Power Comparison', 
    J. Jpn. Comp. Statist. Vol 19 (2006), pp. 1-13. http://doai.io/10.1080/00949655.2010.551516


  Created: 2016.04.06
  Copyright: Steven E. Pav, 2016
  Author: Steven E. Pav <steven@corecast.io>
  Comments: Steven E. Pav
*/

#ifndef __DEF_MURAKAMI_STAT__
#define __DEF_MURAKAMI_STAT__

#include <math.h>
#include "util.h"

// 2FIX: deal with NA/NAN here?
#define MAX(a,b) ((a>b)? (a):(b))
#define MIN(a,b) ((a<b)? (a):(b))

#endif /* __DEF_MURAKAMI_STAT__ */

#include <Rcpp.h>
using namespace Rcpp;

// a weird interface, but basically take the sample_set
// and compute a bunch of the Murakami B_fx statistics,
// for f=0,1,2,3,4,5, set by the 'flavor'. We increment
// through a bunch of sample_sets for efficiency.
//
// later make this back to a template to squeeze out the last bits of efficiency?
//template <int flavor>
//NumericVector murakami_pre_B(const size_t N,const size_t nx,IntegerVector parts,const size_t numits) {
NumericVector murakami_pre_B(const size_t N,const size_t nx,IntegerVector parts,const size_t numits,const int flavor) {
	if ((flavor < 0) || (flavor > 5)) { stop("unsupported flavor."); }
	int iii,jjj,nnn;
	double evx,vvx,Np1,nxp1;
	double nonce,npart,dpart,bplus;
	// preallocate
	NumericVector B1(numits);
	const int ny = (int)N - (int)nx;

	Np1  = (double)N + 1.0;
	nxp1 = (double)nx + 1.0;

	switch(flavor) {
		case 0:
		case 2:
			evx = ((double)N) / ((double)nx);
			vvx = ((double)ny) * evx;
			break;
		default:
			evx = Np1 / nxp1;
			vvx = ((double)ny * Np1) / ((double)nx + 2.0);
			break;
	}

	//if ((flavor == 0) || (flavor == 2)) {
		//evx = (double)N / (double)nx;
		//vvx = (double)ny * evx;
	//}
	//if ((flavor == 1) || (flavor == 3) || (flavor == 4) || (flavor == 5)) {
		//evx = Np1 / nxp1;
		//vvx = (double)ny * Np1 / ((double)nx + 2.0);
	//}
	for (jjj=0;jjj<numits;jjj++) {
		// for more on partitions nonsense, see also:
		// http://howardhinnant.github.io/combinations.html
		// note this will not scale to the multi-class problem very well...
		if (jjj > 0) {
			// increment
			iii = nx - 1;
			while ((iii > 0) && (parts(iii) == N - (nx - 1) + iii)) {
				iii--;
			}
			parts(iii)++;
			for (;iii < (nx - 1);iii++) {
				parts(iii+1) = parts(iii)+1;
			}
		}
		// just to be sure we have prealloc'ed with zeros:
		B1(jjj) = 0.0;
		for (iii=1;iii<=nx;iii++) {
			nnn = parts(iii-1);
			npart = ((double)nnn - evx * (double)iii);
			nonce = ((double)iii) / nxp1;
			dpart = nonce * (1.0 - nonce) * vvx;

			switch(flavor) {
				case 0:
				case 1:
					bplus = npart * npart / dpart;
					break;
				case 2:
					bplus = (npart * std::abs(npart)) / dpart;
					break;
				case 3:
					bplus = npart * npart / (dpart * dpart);
					break;
				case 4:
					bplus = std::abs(npart) / (dpart * dpart);
					break;
				case 5:
					bplus = npart * npart / log(dpart);
					break;
			}
			B1(jjj) += bplus / ((double) nx);
		}
	}
	return B1;
}

//// from when it was templated...
//// same interface, not templated, supposedly faster because of it
//NumericVector murakami_B(const size_t N,const size_t nx,IntegerVector parts,const size_t numits,int flavor) {
	//NumericVector B1;
	//switch(flavor) {
		//case 0:
			//B1 = murakami_pre_B<0>(N,nx,parts,numits);
			//break;
		//case 1:
			//B1 = murakami_pre_B<1>(N,nx,parts,numits);
			//break;
		//case 2:
			//B1 = murakami_pre_B<2>(N,nx,parts,numits);
			//break;
		//case 3:
			//B1 = murakami_pre_B<3>(N,nx,parts,numits);
			//break;
		//case 4:
			//B1 = murakami_pre_B<4>(N,nx,parts,numits);
			//break;
		//case 5:
			//B1 = murakami_pre_B<5>(N,nx,parts,numits);
			//break;
	//}
	//return B1;
//}


/*
double murakami_raw(int N,NumericVector x,int flavor=0) {
	const size_t nx = (size_t)x.size();
	IntegerVector parts(nx);
	for (int iii=0;iii<nx;iii++) { parts(iii) = (int)x(iii); }
	NumericVector B1 = murakami_pre_B(N,nx,parts,1,flavor);
	return B1(0);
}
*/

NumericVector murakami_many_B(const int N,const int nx,int flavor) {
	IntegerVector parts(nx);
	int iii;
	for (iii=0;iii<nx;iii++) {
		// our partitions are 1-based.
		parts(iii)=iii+1;
	}
	size_t numits = (size_t)Combinatorics::bincoef[N][nx];
	return murakami_pre_B(N,nx,parts,numits,flavor);
}

//' @title
//' Compute Murakami's test statistic.
//'
//' @description
//'
//' Compute one of the modified Baumgartner-Weiss-Schindler test statistics proposed
//' by Murakami, or Neuhauser.
//'
//' @details
//'
//' Given vectors \eqn{X} and \eqn{Y}, computes \eqn{B_{jX}} and \eqn{B_{jY}} 
//' for some \eqn{j} as described by Murakami and by Neuhauser, returning either their 
//' their average or their average distance.
//' The test statistics approximate the weighted square norm of the
//' difference in CDFs of the two distributions. 
//'
//' The test statistic is based only on the ranks of the input. If the same
//' monotonic transform is applied to both vectors, the result should be unchanged.
//'
//' The various \sQuote{flavor}s of test statistic are:
//' \describe{
//' \item{0}{The statistic of Baumgartner-Weiss-Schindler.}
//' \item{1}{Murakami's \eqn{B_1} statistic, from his 2006 paper.}
//' \item{2}{Neuhauser's difference statistic, denoted by Murakami as \eqn{B_2} in his 
//' 2012 paper.}
//' \item{3}{Murakami's \eqn{B_3} statistic, from his 2012 paper.}
//' \item{4}{Murakami's \eqn{B_4} statistic, from his 2012 paper.}
//' \item{5}{Murakami's \eqn{B_5} statistic, from his 2012 paper, with a log weighting.}
//' }
//'
//' @param x a vector of the first sample.
//' @param y a vector of the second sample.
//' @param flavor which \sQuote{flavor} of test statistic. 
//' @param nx the length of \code{x}, the first sample.
//' @param ny the length of \code{y}, the second sample.
//'
//' @return The BWS test statistic, \eqn{B_j}. For \code{murakami_stat_perms}, a vector of
//' the test statistics for \emph{all} permutations of the input.
//' @note \code{NA} and \code{NaN} are not yet dealt with!
//' @seealso \code{\link{bws_stat}}.
//' @examples
//'
//' set.seed(1234)
//' x <- runif(1000)
//' y <- runif(100)
//' bval <- murakami_stat(x,y,1)
//'
//' \dontrun{
//' nx <- 6
//' ny <- 5
//' # monte carlo
//' set.seed(1234)
//' repli <- replicate(3000,murakami_stat(rnorm(nx),rnorm(ny),0L))
//' # under the null, perform the permutation test:
//' allem <- murakami_stat_perms(nx,ny,0L)
//' plot(ecdf(allem)) 
//' lines(ecdf(repli),col='red') 
//' }
//'
//' @template etc
//' @template ref-bws
//' @template ref-modtests
//' @rdname murakami_stat
//' @export
// [[Rcpp::export]]
double murakami_stat(NumericVector x,NumericVector y,int flavor=0) {
	// 2FIX: add na_rm?
	// put into the form to be consumed by the other function...
	NumericVector sortx = clone(x); std::sort(sortx.begin(), sortx.end());
	NumericVector sorty = clone(y); std::sort(sorty.begin(), sorty.end());
	IntegerVector G = full_rank<NumericVector, double>(sorty, sortx);
	IntegerVector H = full_rank<NumericVector, double>(sortx, sorty);

	const size_t nx = (size_t)x.size();
	const size_t ny = (size_t)y.size();
	const size_t N = nx + ny;

	NumericVector B1 = murakami_pre_B(N,nx,G,1,flavor);
	NumericVector B2 = murakami_pre_B(N,ny,H,1,flavor);
	double B;

	switch(flavor) {
		case 2:
			B = 0.5 * (B2(0) - B1(0));
			break;
		default:
			B = 0.5 * (B2(0) + B1(0));
			break;
	}
	return B;
}
//' @rdname murakami_stat
//' @export
// [[Rcpp::export]]
NumericVector murakami_stat_perms(int nx, int ny,int flavor=0) {
	const int N = nx + ny;
	if (N > Combinatorics::BINCOEF_MAX_ORD) { stop("N too large"); }

	const bool IS_SYMMETRIC = (nx == ny);
	NumericVector B1,B2,B;
	B1 = murakami_many_B(N,nx,flavor);
	int nlen=B1.size();
	B = NumericVector(nlen);
	int iii;
	if (IS_SYMMETRIC) {
		// combine them;
		if (flavor == 2) {
			for (iii=0;iii < nlen;iii++) {
				B(iii) = 0.5 * (B1(nlen-iii-1) - B1(iii));
			}
		} else {
			for (iii=0;iii < nlen;iii++) {
				B(iii) = 0.5 * (B1(nlen-iii-1) + B1(iii));
			}
		}
	} else {
		B2 = murakami_many_B(N,ny,flavor);
		// combine them;
		if (flavor == 2) {
			for (iii=0;iii < nlen;iii++) {
				B(iii) = 0.5 * (B2(nlen-iii-1) - B1(iii));
			}
		} else {
			for (iii=0;iii < nlen;iii++) {
				B(iii) = 0.5 * (B2(nlen-iii-1) + B1(iii));
			}
		}
	}
	return B;
}

// for speed comparison, timings of the new method defined here:
//    a <- murakami_all_B(5, 5, 1)     24.12     25.35     29.67     26.77     32.91     51.08   100
//    b <- murakami_all_B(7, 7, 1)    366.92    370.64    420.82    379.81    456.58   1235.38   100
//   cc <- murakami_all_B(9, 9, 1)   7106.97   7192.79   7588.08   7304.08   7566.63  10396.69   100
// ee <- murakami_all_B(11, 11, 1) 125301.28 126489.78 135681.89 127329.02 149259.48 183445.64   100
// and the old method based on partitions::setparts and whatnot:
// aa <- murakami_BB(5, 5, 1)     646.9     719.4     835.2     818.4     885.1    1636   100
// bb <- murakami_BB(7, 7, 1)   20391.0   20632.6   21556.7   20869.3   21993.5   26766   100
// cc <- murakami_BB(9, 9, 1) 1374527.5 1405699.7 1474427.9 1449461.6 1497241.3 2010725   100
// I never got up to 11 & 11 with the old method.
// the takeaway is that we are much faster now. yay.
// 
//

/*

		 3/2 split:
 
		 123 45
		 124 35
		 125 34
		 134 25
		 135 24
		 145 23
		 234 15
		 235 14
		 245 13
		 345 12


		 2/3 split:

		 12 345
		 13 245
		 14 235
		 15 234
		 23 145
		 24 135
		 25 134
		 34 124
		 35 124
		 45 123

		 require(partitions)
		 M <- partitions::setparts(c(9,6)) 
		 allem <- sapply(seq_len(ncol(M)),function(ccc) {
		 		xv <- which(M[,ccc] == 1)
		 		yv <- which(M[,ccc] == 2)
				murakami_stat(xv,yv,flavor=2)
		 })

		 require(partitions)
		 M <- partitions::setparts(c(3,2)) 
		 allem <- sapply(seq_len(ncol(M)),function(ccc) {
		 		xv <- which(M[,ccc] == 1)
		 		yv <- which(M[,ccc] == 2)
				murakami_stat(xv,yv,flavor=2)
		 })

		 sort(allem)

		 fooz <- c(murakami_stat(c(1,2,3),c(4,5),flavor=2),
		 murakami_stat(c(1,2,4),c(3,5),flavor=2),
		 murakami_stat(c(1,2,5),c(3,4),flavor=2),
		 murakami_stat(c(1,3,4),c(2,5),flavor=2),
		 murakami_stat(c(1,3,5),c(2,4),flavor=2),
		 murakami_stat(c(1,4,5),c(2,3),flavor=2),
		 murakami_stat(c(2,3,4),c(1,5),flavor=2),
		 murakami_stat(c(2,3,5),c(1,4),flavor=2),
		 murakami_stat(c(2,4,5),c(1,3),flavor=2),
		 murakami_stat(c(3,4,5),c(1,2),flavor=2))

	sort(fooz) - sort(allem)
	sort(murakami_stat_perms(3,2,flavor=2)) 



 */


//for vim modeline: (do not edit)
// vim:nowrap:ts=2:sw=2:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
