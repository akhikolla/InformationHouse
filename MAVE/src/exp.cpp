/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * exp.cpp
 *
 * Code generation for function 'exp'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "exp.h"

/* Function Definitions */
void b_exp(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[0];
  for (k = 0; k + 1 <= nx; k++) {
    x->data[k] = std::exp(x->data[k]);
  }
}

/* End of code generation (exp.cpp) */
