/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * power.cpp
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "power.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  unsigned int a_idx_0;
  int k;
  a_idx_0 = (unsigned int)a->size[0];
  k = y->size[0];
  y->size[0] = (int)a_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  for (k = 0; k + 1 <= a->size[0]; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

void power(const emxArray_real_T *a, emxArray_real_T *y)
{
  int n;
  unsigned int uv1[2];
  int k;
  for (n = 0; n < 2; n++) {
    uv1[n] = (unsigned int)a->size[n];
  }

  n = y->size[0] * y->size[1];
  y->size[0] = (int)uv1[0];
  y->size[1] = (int)uv1[1];
  emxEnsureCapacity((emxArray__common *)y, n, sizeof(double));
  n = a->size[0] * a->size[1];
  for (k = 0; k + 1 <= n; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

/* End of code generation (power.cpp) */
