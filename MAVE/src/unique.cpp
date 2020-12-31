/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * unique.cpp
 *
 * Code generation for function 'unique'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "unique.h"

/* Function Definitions */
void count_nonfinites(const emxArray_real_T *b, int nb, int *nMInf, int *nFinite,
                      int *nPInf, int *nNaN)
{
  int k;
  k = 0;
  while ((k + 1 <= nb) && rtIsInf(b->data[k]) && (b->data[k] < 0.0)) {
    k++;
  }

  *nMInf = k;
  k = nb;
  while ((k >= 1) && rtIsNaN(b->data[k - 1])) {
    k--;
  }

  *nNaN = nb - k;
  while ((k >= 1) && rtIsInf(b->data[k - 1]) && (b->data[k - 1] > 0.0)) {
    k--;
  }

  *nPInf = (nb - k) - *nNaN;
  *nFinite = k - *nMInf;
}

/* End of code generation (unique.cpp) */
