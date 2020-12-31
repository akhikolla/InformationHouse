/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * abs.cpp
 *
 * Code generation for function 'abs'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "abs.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_abs(const emxArray_boolean_T *x, emxArray_real_T *y)
{
  int n;
  unsigned int uv3[2];
  int k;
  for (n = 0; n < 2; n++) {
    uv3[n] = (unsigned int)x->size[n];
  }

  n = y->size[0] * y->size[1];
  y->size[0] = (int)uv3[0];
  y->size[1] = (int)uv3[1];
  emxEnsureCapacity((emxArray__common *)y, n, sizeof(double));
  n = x->size[0] * x->size[1];
  for (k = 0; k + 1 <= n; k++) {
    y->data[k] = x->data[k];
  }
}

void c_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  int n;
  unsigned int uv4[2];
  int k;
  for (n = 0; n < 2; n++) {
    uv4[n] = (unsigned int)x->size[n];
  }

  n = y->size[0] * y->size[1];
  y->size[0] = (int)uv4[0];
  y->size[1] = (int)uv4[1];
  emxEnsureCapacity((emxArray__common *)y, n, sizeof(double));
  n = x->size[0] * x->size[1];
  for (k = 0; k + 1 <= n; k++) {
    y->data[k] = std::abs(x->data[k]);
  }
}

void d_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int x_idx_0;
  int k;
  x_idx_0 = (unsigned int)x->size[0];
  k = y->size[0];
  y->size[0] = (int)x_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  for (k = 0; k + 1 <= x->size[0]; k++) {
    y->data[k] = std::abs(x->data[k]);
  }
}

/* End of code generation (abs.cpp) */
