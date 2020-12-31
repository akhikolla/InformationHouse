/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * qrsolve.cpp
 *
 * Code generation for function 'qrsolve'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "qrsolve.h"

/* Function Definitions */
int rankFromQR(const emxArray_real_T *A)
{
  int r;
  int minmn;
  int maxmn;
  double tol;
  r = 0;
  if (A->size[0] < A->size[1]) {
    minmn = A->size[0];
    maxmn = A->size[1];
  } else {
    minmn = A->size[1];
    maxmn = A->size[0];
  }

  if (minmn > 0) {
    tol = (double)maxmn * std::abs(A->data[0]) * 2.2204460492503131E-16;
    while ((r < minmn) && (std::abs(A->data[r + A->size[0] * r]) >= tol)) {
      r++;
    }
  }

  return r;
}

/* End of code generation (qrsolve.cpp) */
