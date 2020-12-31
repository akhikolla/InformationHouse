/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sum.cpp
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "sum.h"
#include "combine_vector_elements.h"

/* Function Definitions */
void b_sum(const emxArray_real_T *x, emxArray_real_T *y)
{
  b_combine_vector_elements(x, y);
}

double c_sum(const emxArray_real_T *x)
{
  return c_combine_vector_elements(x);
}

double d_sum(const emxArray_real_T *x)
{
  double y;
  int k;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= x->size[0]; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

void sum(const emxArray_real_T *x, emxArray_real_T *y)
{
  combine_vector_elements(x, y);
}

/* End of code generation (sum.cpp) */
