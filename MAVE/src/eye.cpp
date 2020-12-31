/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * eye.cpp
 *
 * Code generation for function 'eye'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "eye.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void eye(double varargin_1, emxArray_real_T *I)
{
  double t;
  int k;
  int loop_ub;
  if (varargin_1 < 0.0) {
    t = 0.0;
  } else {
    t = varargin_1;
  }

  k = I->size[0] * I->size[1];
  I->size[0] = (int)t;
  I->size[1] = (int)t;
  emxEnsureCapacity((emxArray__common *)I, k, sizeof(double));
  loop_ub = (int)t * (int)t;
  for (k = 0; k < loop_ub; k++) {
    I->data[k] = 0.0;
  }

  if ((int)t > 0) {
    for (k = 0; k + 1 <= (int)t; k++) {
      I->data[k + I->size[0] * k] = 1.0;
    }
  }
}

/* End of code generation (eye.cpp) */
