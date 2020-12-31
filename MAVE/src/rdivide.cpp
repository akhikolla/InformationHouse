/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * rdivide.cpp
 *
 * Code generation for function 'rdivide'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "rdivide.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void b_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
               emxArray_real_T *z)
{
  int i12;
  int loop_ub;
  i12 = z->size[0] * z->size[1];
  z->size[0] = x->size[0];
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i12, sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i12 = 0; i12 < loop_ub; i12++) {
    z->data[i12] = x->data[i12] / y->data[i12];
  }
}

void rdivide(const emxArray_real_T *y, emxArray_real_T *z)
{
  int i6;
  int loop_ub;
  i6 = z->size[0];
  z->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)z, i6, sizeof(double));
  loop_ub = y->size[0];
  for (i6 = 0; i6 < loop_ub; i6++) {
    z->data[i6] = 1.0 / y->data[i6];
  }
}

/* End of code generation (rdivide.cpp) */
