/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * norm.cpp
 *
 * Code generation for function 'norm'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "norm.h"

/* Function Definitions */
double norm(const emxArray_real_T *x)
{
  double y;
  double scale;
  int k;
  double absxk;
  double t;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = 0.0;
    if (x->size[0] == 1) {
      y = std::abs(x->data[0]);
    } else {
      scale = 2.2250738585072014E-308;
      for (k = 1; k <= x->size[0]; k++) {
        absxk = std::abs(x->data[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = 1.0 + y * t * t;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

/* End of code generation (norm.cpp) */
