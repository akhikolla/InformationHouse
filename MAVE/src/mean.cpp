/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * mean.cpp
 *
 * Code generation for function 'mean'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "mean.h"
#include "combine_vector_elements.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_mean(const emxArray_real_T *x, emxArray_real_T *y)
{
  int b_x;
  int i3;
  int loop_ub;
  combine_vector_elements(x, y);
  b_x = x->size[1];
  i3 = y->size[0];
  emxEnsureCapacity((emxArray__common *)y, i3, sizeof(double));
  loop_ub = y->size[0];
  for (i3 = 0; i3 < loop_ub; i3++) {
    y->data[i3] /= (double)b_x;
  }
}

double c_mean(const emxArray_real_T *x)
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

  y /= (double)x->size[0];
  return y;
}

double d_mean(const emxArray_real_T *x)
{
  double b_x;
  b_x = c_combine_vector_elements(x);
  return b_x / (double)x->size[1];
}

void mean(const emxArray_real_T *x, emxArray_real_T *y)
{
  int b_y;
  int c_y;
  int b_x;
  b_combine_vector_elements(x, y);
  b_y = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)y, b_y, sizeof(double));
  b_y = y->size[0];
  c_y = y->size[1];
  b_x = x->size[0];
  c_y *= b_y;
  for (b_y = 0; b_y < c_y; b_y++) {
    y->data[b_y] /= (double)b_x;
  }
}

/* End of code generation (mean.cpp) */
