/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * quantile.cpp
 *
 * Code generation for function 'quantile'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "quantile.h"
#include "MAVEfast_emxutil.h"
#include "prctile.h"
#include "linspace.h"

/* Function Definitions */
void quantile(const emxArray_real_T *x, const emxArray_real_T *p,
              emxArray_real_T *y)
{
  emxArray_real_T *p100;
  int i9;
  double delta;
  int loop_ub;
  emxInit_real_T1(&p100, 2);
  if ((p->size[1] == 1) && (p->data[0] == std::floor(p->data[0])) && (p->data[0]
       > 1.0)) {
    delta = 100.0 / (p->data[0] + 1.0);
    linspace(delta, 100.0 - delta, p->data[0], p100);
  } else {
    i9 = p100->size[0] * p100->size[1];
    p100->size[0] = 1;
    p100->size[1] = p->size[1];
    emxEnsureCapacity((emxArray__common *)p100, i9, sizeof(double));
    loop_ub = p->size[0] * p->size[1];
    for (i9 = 0; i9 < loop_ub; i9++) {
      p100->data[i9] = 100.0 * p->data[i9];
    }
  }

  prctile(x, p100, y);
  emxFree_real_T(&p100);
}

/* End of code generation (quantile.cpp) */
