
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

  Created: 2016.04.24
  Copyright: Steven E. Pav, 2016
  Author: Steven E. Pav <steven@corecast.io>
  Comments: Steven E. Pav
*/

#ifndef __DEF_BWS_UTIL__
#define __DEF_BWS_UTIL__

#include <Rcpp.h>
using namespace Rcpp;

// return number of elements in sortx[lo:hi) <= val
// or lo - 1 if none.
// these are C indexed, so you should have
// lo = 0 and hi = length(sortx) to do the
// whole vector.
template <typename T,typename D,bool comp_with_eq>
int bin_search_lstar(T sortx, D val, const int lo, const int hi) {
	int kidx;
	int ilo=lo;
	int ihi=hi;

	if (ilo < 0) { stop("out of bounds"); }

	while (ilo < ihi) {
		kidx = (ilo + ihi) / 2;

		if (comp_with_eq) {
			if (val <= sortx[kidx]) {
				ihi = kidx;
			} else {
				ilo = kidx + 1;
			}
		} else {
			if (val < sortx[kidx]) {
				ihi = kidx;
			} else {
				ilo = kidx + 1;
			}
		}
	}
	return ilo - 1;
}

// for each val in y, 
// return number of elements in sortx[lo:hi) <= val
// or lo - 1 if none.
// these are C indexed, so you should have
// lo = 0 and hi = length(sortx) to do the
// whole vector.
template <typename T,typename D,bool comp_with_eq>
IntegerVector zip_index_lstar(T sortx, T refy, const int lo, const int hi) {
	int kidx;
	int xidx, yidx, lastv;
	int ynel = refy.length();
	IntegerVector retv(ynel);

	if (ynel == 1) {
		retv[0] = bin_search_lstar<T,D,comp_with_eq>(sortx, refy[0], lo, hi);
	} else {
		if (lo < 0) { stop("out of bounds"); }
		xidx = lo;
		yidx = 0;
		if (comp_with_eq) {

			while ((xidx < hi) && (yidx < ynel)) {
				if (sortx[xidx] <= refy[yidx]) {
					xidx++;
				} else {
					retv[yidx] = xidx - 1;
					yidx++;
				}
			}

		} else {

			while ((xidx < hi) && (yidx < ynel)) {
				if (sortx[xidx] < refy[yidx]) {
					xidx++;
				} else {
					retv[yidx] = xidx - 1;
					yidx++;
				}
			}

		} 
		lastv = xidx - 1;
		while (yidx < ynel) {
			retv[yidx] = lastv;
			yidx++;
		}
	}
	return retv;
}

//2FIX: this is a stupid hacky way to compute this.
//smarter would be to modify zip_index_lstar to 
//do the addition for you?
//ah, but ties, right.
// given sorted x and sorted y, for each element in y
// find the number of elements in union(x,y) less than
// or equal to it. that is
// retv[i] = # { z in union(x,y) | z <= y[i] } for 1 <= i <= length(y)
template <typename T,typename D>
IntegerVector full_rank(T sortx, T sorty) {
	IntegerVector retv;
	retv = 2 + zip_index_lstar<T, D, true>(sortx, sorty, 0, sortx.length()) +
		zip_index_lstar<T, D, true>(sorty, sorty, 0, sorty.length());
	return retv;
}

/*
 
set.seed(1234)
x <- rnorm(1000)
y <- rnorm(1000)
stopifnot(all(unlist(lapply(sort(x),function(anx) { sum(c(x,y) <= anx) })) == fool(x,y)))

//// [[Rcpp::export]]
//IntegerVector fool(NumericVector x,NumericVector y) {
	//NumericVector sortx = clone(x); std::sort(sortx.begin(), sortx.end());
	//NumericVector sorty = clone(y); std::sort(sorty.begin(), sorty.end());
	//IntegerVector G = full_rank<NumericVector, double>(sorty, sortx);
	//return G;
//}

*/


// namespaces FTW
// http://stackoverflow.com/a/20794407/164611

// preallocate the Binomial Coefficients, for efficiency

// this is R code used to generate C code. [ducks]

//MAXORD <- 32
//refv <- matrix(0L,nrow=MAXORD,ncol=MAXORD)
//for (iii in seq(1,nrow(refv))) {
	 //refv[iii,1] = 1L; refv[iii,iii] = 1L; 
	 //if (iii > 2) { for (jjj in seq(2,iii-1)) { refv[iii,jjj] = refv[iii-1,jjj-1] + refv[iii-1,jjj]; } }
 //} 
//cat(sprintf('namespace Combinatorics\n{\nconst int BINCOEF_MAX_ORD = %d;\nconst int bincoef[%d][%d] = {',ncol(refv)-1,ncol(refv),ncol(refv)),
		//paste0('\n{ ',lapply(seq_len(nrow(refv)),function(rn) { paste0(sprintf('%9s',as.character(refv[rn,])),collapse=', ') }),'},'),
		//'};\n}\n\n',file='/tmp/binc.txt')

#define MAX_ORD 31
namespace Combinatorics
{
const int BINCOEF_MAX_ORD = 31;
const int bincoef[32][32] = { 
{         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         2,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         3,         3,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         4,         6,         4,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         5,        10,        10,         5,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         6,        15,        20,        15,         6,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         7,        21,        35,        35,        21,         7,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         8,        28,        56,        70,        56,        28,         8,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,         9,        36,        84,       126,       126,        84,        36,         9,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        10,        45,       120,       210,       252,       210,       120,        45,        10,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        11,        55,       165,       330,       462,       462,       330,       165,        55,        11,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        12,        66,       220,       495,       792,       924,       792,       495,       220,        66,        12,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        13,        78,       286,       715,      1287,      1716,      1716,      1287,       715,       286,        78,        13,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        14,        91,       364,      1001,      2002,      3003,      3432,      3003,      2002,      1001,       364,        91,        14,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        15,       105,       455,      1365,      3003,      5005,      6435,      6435,      5005,      3003,      1365,       455,       105,        15,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        16,       120,       560,      1820,      4368,      8008,     11440,     12870,     11440,      8008,      4368,      1820,       560,       120,        16,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        17,       136,       680,      2380,      6188,     12376,     19448,     24310,     24310,     19448,     12376,      6188,      2380,       680,       136,        17,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        18,       153,       816,      3060,      8568,     18564,     31824,     43758,     48620,     43758,     31824,     18564,      8568,      3060,       816,       153,        18,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        19,       171,       969,      3876,     11628,     27132,     50388,     75582,     92378,     92378,     75582,     50388,     27132,     11628,      3876,       969,       171,        19,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        20,       190,      1140,      4845,     15504,     38760,     77520,    125970,    167960,    184756,    167960,    125970,     77520,     38760,     15504,      4845,      1140,       190,        20,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        21,       210,      1330,      5985,     20349,     54264,    116280,    203490,    293930,    352716,    352716,    293930,    203490,    116280,     54264,     20349,      5985,      1330,       210,        21,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        22,       231,      1540,      7315,     26334,     74613,    170544,    319770,    497420,    646646,    705432,    646646,    497420,    319770,    170544,     74613,     26334,      7315,      1540,       231,        22,         1,         0,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        23,       253,      1771,      8855,     33649,    100947,    245157,    490314,    817190,   1144066,   1352078,   1352078,   1144066,    817190,    490314,    245157,    100947,     33649,      8855,      1771,       253,        23,         1,         0,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        24,       276,      2024,     10626,     42504,    134596,    346104,    735471,   1307504,   1961256,   2496144,   2704156,   2496144,   1961256,   1307504,    735471,    346104,    134596,     42504,     10626,      2024,       276,        24,         1,         0,         0,         0,         0,         0,         0,         0}, 
{         1,        25,       300,      2300,     12650,     53130,    177100,    480700,   1081575,   2042975,   3268760,   4457400,   5200300,   5200300,   4457400,   3268760,   2042975,   1081575,    480700,    177100,     53130,     12650,      2300,       300,        25,         1,         0,         0,         0,         0,         0,         0}, 
{         1,        26,       325,      2600,     14950,     65780,    230230,    657800,   1562275,   3124550,   5311735,   7726160,   9657700,  10400600,   9657700,   7726160,   5311735,   3124550,   1562275,    657800,    230230,     65780,     14950,      2600,       325,        26,         1,         0,         0,         0,         0,         0}, 
{         1,        27,       351,      2925,     17550,     80730,    296010,    888030,   2220075,   4686825,   8436285,  13037895,  17383860,  20058300,  20058300,  17383860,  13037895,   8436285,   4686825,   2220075,    888030,    296010,     80730,     17550,      2925,       351,        27,         1,         0,         0,         0,         0}, 
{         1,        28,       378,      3276,     20475,     98280,    376740,   1184040,   3108105,   6906900,  13123110,  21474180,  30421755,  37442160,  40116600,  37442160,  30421755,  21474180,  13123110,   6906900,   3108105,   1184040,    376740,     98280,     20475,      3276,       378,        28,         1,         0,         0,         0}, 
{         1,        29,       406,      3654,     23751,    118755,    475020,   1560780,   4292145,  10015005,  20030010,  34597290,  51895935,  67863915,  77558760,  77558760,  67863915,  51895935,  34597290,  20030010,  10015005,   4292145,   1560780,    475020,    118755,     23751,      3654,       406,        29,         1,         0,         0}, 
{         1,        30,       435,      4060,     27405,    142506,    593775,   2035800,   5852925,  14307150,  30045015,  54627300,  86493225, 119759850, 145422675, 155117520, 145422675, 119759850,  86493225,  54627300,  30045015,  14307150,   5852925,   2035800,    593775,    142506,     27405,      4060,       435,        30,         1,         0}, 
{         1,        31,       465,      4495,     31465,    169911,    736281,   2629575,   7888725,  20160075,  44352165,  84672315, 141120525, 206253075, 265182525, 300540195, 300540195, 265182525, 206253075, 141120525,  84672315,  44352165,  20160075,   7888725,   2629575,    736281,    169911,     31465,      4495,       465,        31,         1} };
}

#endif /* __DEF_BWS_UTIL__ */

//for vim modeline: (do not edit)
// vim:nowrap:ts=2:sw=2:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
