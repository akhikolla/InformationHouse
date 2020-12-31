
/*
 

  Created: 2017.03.19
  Copyright: Steven E. Pav, 2017
  Author: Steven E. Pav <steven@gilgamath.com>
  Comments: Steven E. Pav
*/

#include <Rcpp.h>
using namespace Rcpp;

/*
 https://github.com/RcppCore/Rcpp/issues/636
*/

void R_init_BWStest(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}

//for vim modeline: (do not edit)
// vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
