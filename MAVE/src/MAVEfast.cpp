/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * MAVEfast.cpp
 *
 * Code generation for function 'MAVEfast'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "MAVEfast_emxutil.h"
#include "power.h"
#include "eye.h"
#include "mldivide.h"
#include "repmat.h"
#include "exp.h"
#include "sort1.h"
#include "diag.h"
#include "eig.h"
#include "kron.h"
#include "sum.h"
#include "strcmp.h"
#include "mean.h"
#include "std.h"
#include "norm.h"
#include "floor.h"
#include "abs.h"
#include "unique.h"
#include "sortIdx.h"
#include "quantile.h"
#include "rdivide.h"
#include "sqrt.h"
#include "upper.h"
#include "bsxfun.h"

/* Function Declarations */
static void dcorrVS(const emxArray_real_T *X, const emxArray_real_T *Y,
                    emxArray_real_T *sele);
static double rt_powd_snf(double u0, double u1);
static void unifD(const emxArray_real_T *x, double m, emxArray_real_T *I);

/* Function Definitions */
static void dcorrVS(const emxArray_real_T *X, const emxArray_real_T *Y,
                    emxArray_real_T *sele)
{
  emxArray_real_T *x2;
  emxArray_real_T *ztemp;
  emxArray_real_T *b_x2;
  int i0;
  int br;
  emxArray_real_T *b;
  emxArray_real_T *varargin_1;
  int ar;
  int k;
  int i1;
  unsigned int uv0[2];
  int m;
  int n;
  int c;
  int ic;
  int ib;
  int ia;
  emxArray_real_T *r0;
  double A;
  int b_b[1];
  emxArray_real_T c_b;
  emxArray_real_T *d_b;
  emxArray_real_T *b_c;
  double dcovYY;
  emxArray_real_T *a;
  emxArray_real_T *x;
  emxArray_real_T *b_a;
  emxArray_real_T *c_a;
  emxArray_real_T *c_x2;
  emxArray_real_T *b_x;
  emxArray_int32_T *iidx;
  int d_a[1];
  double dcovXX;
  emxInit_real_T(&x2, 1);
  emxInit_real_T1(&ztemp, 2);
  emxInit_real_T1(&b_x2, 2);

  /* DCORR computes the distance correlation between two random variables */
  /*  X and Y. Rows represent the examples, and columns the variables. */
  /*  Date: 10-03-2015 */
  /*  Author: Paolo Inglese (paolo.ingls@gmail.com) */
  /*  Reference: http://en.wikipedia.org/wiki/Distance_correlation */
  /* b = pdist2(Y,Y); %conversion of pdist2 not supported */
  power(Y, ztemp);
  sum(ztemp, x2);
  i0 = b_x2->size[0] * b_x2->size[1];
  b_x2->size[0] = 1;
  b_x2->size[1] = x2->size[0];
  emxEnsureCapacity((emxArray__common *)b_x2, i0, sizeof(double));
  br = x2->size[0];
  for (i0 = 0; i0 < br; i0++) {
    b_x2->data[b_x2->size[0] * i0] = x2->data[i0];
  }

  emxInit_real_T1(&b, 2);
  emxInit_real_T1(&varargin_1, 2);
  bsxfun(x2, b_x2, varargin_1);
  i0 = b->size[0] * b->size[1];
  b->size[0] = Y->size[1];
  b->size[1] = Y->size[0];
  emxEnsureCapacity((emxArray__common *)b, i0, sizeof(double));
  br = Y->size[0];
  emxFree_real_T(&b_x2);
  for (i0 = 0; i0 < br; i0++) {
    ar = Y->size[1];
    for (i1 = 0; i1 < ar; i1++) {
      b->data[i1 + b->size[0] * i0] = Y->data[i0 + Y->size[0] * i1];
    }
  }

  if ((Y->size[1] == 1) || (b->size[0] == 1)) {
    i0 = ztemp->size[0] * ztemp->size[1];
    ztemp->size[0] = Y->size[0];
    ztemp->size[1] = b->size[1];
    emxEnsureCapacity((emxArray__common *)ztemp, i0, sizeof(double));
    br = Y->size[0];
    for (i0 = 0; i0 < br; i0++) {
      ar = b->size[1];
      for (i1 = 0; i1 < ar; i1++) {
        ztemp->data[i0 + ztemp->size[0] * i1] = 0.0;
        c = Y->size[1];
        for (n = 0; n < c; n++) {
          ztemp->data[i0 + ztemp->size[0] * i1] += Y->data[i0 + Y->size[0] * n] *
            b->data[n + b->size[0] * i1];
        }
      }
    }
  } else {
    k = Y->size[1];
    uv0[0] = (unsigned int)Y->size[0];
    uv0[1] = (unsigned int)b->size[1];
    i0 = ztemp->size[0] * ztemp->size[1];
    ztemp->size[0] = (int)uv0[0];
    ztemp->size[1] = (int)uv0[1];
    emxEnsureCapacity((emxArray__common *)ztemp, i0, sizeof(double));
    m = Y->size[0];
    i0 = ztemp->size[0] * ztemp->size[1];
    emxEnsureCapacity((emxArray__common *)ztemp, i0, sizeof(double));
    br = ztemp->size[1];
    for (i0 = 0; i0 < br; i0++) {
      ar = ztemp->size[0];
      for (i1 = 0; i1 < ar; i1++) {
        ztemp->data[i1 + ztemp->size[0] * i0] = 0.0;
      }
    }

    if ((Y->size[0] == 0) || (b->size[1] == 0)) {
    } else {
      c = Y->size[0] * (b->size[1] - 1);
      n = 0;
      while ((m > 0) && (n <= c)) {
        i0 = n + m;
        for (ic = n; ic + 1 <= i0; ic++) {
          ztemp->data[ic] = 0.0;
        }

        n += m;
      }

      br = 0;
      n = 0;
      while ((m > 0) && (n <= c)) {
        ar = -1;
        i0 = br + k;
        for (ib = br; ib + 1 <= i0; ib++) {
          if (b->data[ib] != 0.0) {
            ia = ar;
            i1 = n + m;
            for (ic = n; ic + 1 <= i1; ic++) {
              ia++;
              ztemp->data[ic] += b->data[ib] * Y->data[ia];
            }
          }

          ar += m;
        }

        br += k;
        n += m;
      }
    }
  }

  i0 = varargin_1->size[0] * varargin_1->size[1];
  emxEnsureCapacity((emxArray__common *)varargin_1, i0, sizeof(double));
  n = varargin_1->size[0];
  c = varargin_1->size[1];
  br = n * c;
  for (i0 = 0; i0 < br; i0++) {
    varargin_1->data[i0] -= 2.0 * ztemp->data[i0];
  }

  for (i0 = 0; i0 < 2; i0++) {
    uv0[i0] = (unsigned int)varargin_1->size[i0];
  }

  i0 = ztemp->size[0] * ztemp->size[1];
  ztemp->size[0] = (int)uv0[0];
  ztemp->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)ztemp, i0, sizeof(double));
  i0 = b->size[0] * b->size[1];
  b->size[0] = (int)uv0[0];
  b->size[1] = (int)uv0[1];
  emxEnsureCapacity((emxArray__common *)b, i0, sizeof(double));
  n = ztemp->size[0] * ztemp->size[1];
  for (k = 0; k + 1 <= n; k++) {
    A = varargin_1->data[k];
    if (!(A > 0.0)) {
      A = 0.0;
    }

    b->data[k] = A;
  }

  emxInit_real_T1(&r0, 2);
  b_sqrt(b);

  /* 2018-05-16 */
  mean(b, r0);
  b_mean(b, x2);
  b_bsxfun(r0, x2, ztemp);
  n = b->size[0];
  c = b->size[1];
  b_b[0] = n * c;
  c_b = *b;
  c_b.size = (int *)&b_b;
  c_b.numDimensions = 1;
  A = c_mean(&c_b);
  i0 = b->size[0] * b->size[1];
  emxEnsureCapacity((emxArray__common *)b, i0, sizeof(double));
  br = b->size[1];
  for (i0 = 0; i0 < br; i0++) {
    ar = b->size[0];
    for (i1 = 0; i1 < ar; i1++) {
      b->data[i1 + b->size[0] * i0] = (b->data[i1 + b->size[0] * i0] -
        ztemp->data[i1 + ztemp->size[0] * i0]) + A;
    }
  }

  emxInit_real_T1(&d_b, 2);
  i0 = d_b->size[0] * d_b->size[1];
  d_b->size[0] = b->size[0];
  d_b->size[1] = b->size[1];
  emxEnsureCapacity((emxArray__common *)d_b, i0, sizeof(double));
  br = b->size[0] * b->size[1];
  for (i0 = 0; i0 < br; i0++) {
    d_b->data[i0] = b->data[i0] * b->data[i0];
  }

  emxInit_real_T(&b_c, 1);
  b_sum(d_b, r0);
  dcovYY = c_sum(r0);
  i0 = b_c->size[0];
  b_c->size[0] = X->size[1];
  emxEnsureCapacity((emxArray__common *)b_c, i0, sizeof(double));
  ib = 0;
  emxFree_real_T(&d_b);
  emxInit_real_T1(&a, 2);
  emxInit_real_T(&x, 1);
  emxInit_real_T1(&b_a, 2);
  emxInit_real_T1(&c_a, 2);
  emxInit_real_T1(&c_x2, 2);
  emxInit_real_T1(&b_x, 2);
  while (ib <= X->size[1] - 1) {
    /* a = pdist2(X(:,i), X(:,i)); % matrix euclidean distances */
    br = X->size[0];
    i0 = x->size[0];
    x->size[0] = br;
    emxEnsureCapacity((emxArray__common *)x, i0, sizeof(double));
    for (i0 = 0; i0 < br; i0++) {
      x->data[i0] = X->data[i0 + X->size[0] * ib];
    }

    b_power(x, x2);
    i0 = c_x2->size[0] * c_x2->size[1];
    c_x2->size[0] = 1;
    c_x2->size[1] = x2->size[0];
    emxEnsureCapacity((emxArray__common *)c_x2, i0, sizeof(double));
    br = x2->size[0];
    for (i0 = 0; i0 < br; i0++) {
      c_x2->data[c_x2->size[0] * i0] = x2->data[i0];
    }

    bsxfun(x2, c_x2, ztemp);
    i0 = b_x->size[0] * b_x->size[1];
    b_x->size[0] = x->size[0];
    b_x->size[1] = x->size[0];
    emxEnsureCapacity((emxArray__common *)b_x, i0, sizeof(double));
    br = x->size[0];
    for (i0 = 0; i0 < br; i0++) {
      ar = x->size[0];
      for (i1 = 0; i1 < ar; i1++) {
        b_x->data[i0 + b_x->size[0] * i1] = x->data[i0] * x->data[i1];
      }
    }

    i0 = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = ztemp->size[0];
    varargin_1->size[1] = ztemp->size[1];
    emxEnsureCapacity((emxArray__common *)varargin_1, i0, sizeof(double));
    br = ztemp->size[1];
    for (i0 = 0; i0 < br; i0++) {
      ar = ztemp->size[0];
      for (i1 = 0; i1 < ar; i1++) {
        varargin_1->data[i1 + varargin_1->size[0] * i0] = ztemp->data[i1 +
          ztemp->size[0] * i0] - 2.0 * b_x->data[i1 + b_x->size[0] * i0];
      }
    }

    for (i0 = 0; i0 < 2; i0++) {
      uv0[i0] = (unsigned int)varargin_1->size[i0];
    }

    i0 = ztemp->size[0] * ztemp->size[1];
    ztemp->size[0] = (int)uv0[0];
    ztemp->size[1] = (int)uv0[1];
    emxEnsureCapacity((emxArray__common *)ztemp, i0, sizeof(double));
    i0 = a->size[0] * a->size[1];
    a->size[0] = (int)uv0[0];
    a->size[1] = (int)uv0[1];
    emxEnsureCapacity((emxArray__common *)a, i0, sizeof(double));
    n = ztemp->size[0] * ztemp->size[1];
    for (k = 0; k + 1 <= n; k++) {
      A = varargin_1->data[k];
      if (!(A > 0.0)) {
        A = 0.0;
      }

      a->data[k] = A;
    }

    b_sqrt(a);

    /* 2018-05-16 */
    mean(a, r0);
    b_mean(a, x2);
    b_bsxfun(r0, x2, ztemp);
    c = a->size[0];
    n = a->size[1];
    d_a[0] = c * n;
    c_b = *a;
    c_b.size = (int *)&d_a;
    c_b.numDimensions = 1;
    A = c_mean(&c_b);
    i0 = a->size[0] * a->size[1];
    emxEnsureCapacity((emxArray__common *)a, i0, sizeof(double));
    br = a->size[1];
    for (i0 = 0; i0 < br; i0++) {
      ar = a->size[0];
      for (i1 = 0; i1 < ar; i1++) {
        a->data[i1 + a->size[0] * i0] = (a->data[i1 + a->size[0] * i0] -
          ztemp->data[i1 + ztemp->size[0] * i0]) + A;
      }
    }

    i0 = c_a->size[0] * c_a->size[1];
    c_a->size[0] = a->size[0];
    c_a->size[1] = a->size[1];
    emxEnsureCapacity((emxArray__common *)c_a, i0, sizeof(double));
    br = a->size[0] * a->size[1];
    for (i0 = 0; i0 < br; i0++) {
      c_a->data[i0] = a->data[i0] * b->data[i0];
    }

    b_sum(c_a, r0);
    A = c_sum(r0);
    i0 = b_a->size[0] * b_a->size[1];
    b_a->size[0] = a->size[0];
    b_a->size[1] = a->size[1];
    emxEnsureCapacity((emxArray__common *)b_a, i0, sizeof(double));
    br = a->size[0] * a->size[1];
    for (i0 = 0; i0 < br; i0++) {
      b_a->data[i0] = a->data[i0] * a->data[i0];
    }

    b_sum(b_a, r0);
    dcovXX = c_sum(r0);
    b_c->data[ib] = std::sqrt(A / std::sqrt(dcovXX * dcovYY));
    ib++;
  }

  emxFree_real_T(&b_x);
  emxFree_real_T(&c_x2);
  emxFree_real_T(&c_a);
  emxFree_real_T(&b_a);
  emxFree_real_T(&r0);
  emxFree_real_T(&x);
  emxFree_real_T(&ztemp);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxInit_int32_T(&iidx, 1);
  sort(b_c, iidx);
  i0 = x2->size[0];
  x2->size[0] = iidx->size[0];
  emxEnsureCapacity((emxArray__common *)x2, i0, sizeof(double));
  br = iidx->size[0];
  emxFree_real_T(&b_c);
  for (i0 = 0; i0 < br; i0++) {
    x2->data[i0] = iidx->data[i0];
  }

  emxFree_int32_T(&iidx);
  i0 = sele->size[0] * sele->size[1];
  sele->size[0] = 1;
  sele->size[1] = x2->size[0];
  emxEnsureCapacity((emxArray__common *)sele, i0, sizeof(double));
  br = x2->size[0];
  for (i0 = 0; i0 < br; i0++) {
    sele->data[i0] = x2->data[i0];
  }

  emxFree_real_T(&x2);
}

static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = std::abs(u0);
    d1 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = 1.0;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

static void unifD(const emxArray_real_T *x, double m, emxArray_real_T *I)
{
  int p;
  emxArray_real_T *y;
  double k;
  int b_k;
  int i;
  int loop_ub;
  emxArray_real_T *di;
  int xoffset;
  boolean_T empty_non_axis_sizes;
  int nrows;
  int result;
  emxArray_real_T *Xremain;
  emxArray_int32_T *Ii;
  emxArray_int32_T *b_Ii;
  emxArray_int32_T *iidx;
  emxArray_real_T *b_x;
  emxArray_boolean_T *b;
  emxArray_real_T *b_Xremain;
  emxArray_real_T *b_y;
  emxArray_real_T *c_Xremain;
  emxArray_real_T *d_Xremain;
  emxArray_real_T *b_I;
  emxArray_real_T *e_Xremain;
  double s;
  int nrowx;
  int ncolx;

  /*  m is the number of points where teh gradients are calculated */
  m = std::floor(m);
  p = x->size[1];
  emxInit_real_T1(&y, 2);
  if (std::ceil((double)x->size[0] / m) == 1.0) {
    if (x->size[0] < 1) {
      b_k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)y, b_k, sizeof(double));
    } else {
      b_k = x->size[0];
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)((double)b_k - 1.0) + 1;
      emxEnsureCapacity((emxArray__common *)y, i, sizeof(double));
      loop_ub = (int)((double)b_k - 1.0);
      for (b_k = 0; b_k <= loop_ub; b_k++) {
        y->data[y->size[0] * b_k] = 1.0 + (double)b_k;
      }
    }

    b_k = I->size[0] * I->size[1];
    I->size[0] = 1;
    I->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)I, b_k, sizeof(double));
    loop_ub = y->size[0] * y->size[1];
    for (b_k = 0; b_k < loop_ub; b_k++) {
      I->data[b_k] = y->data[b_k];
    }
  } else {
    k = std::floor((double)x->size[0] / m);
    if (x->size[0] < 1) {
      b_k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)y, b_k, sizeof(double));
    } else {
      b_k = x->size[0];
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)((double)b_k - 1.0) + 1;
      emxEnsureCapacity((emxArray__common *)y, i, sizeof(double));
      loop_ub = (int)((double)b_k - 1.0);
      for (b_k = 0; b_k <= loop_ub; b_k++) {
        y->data[y->size[0] * b_k] = 1.0 + (double)b_k;
      }
    }

    emxInit_real_T(&di, 1);
    b_k = di->size[0];
    di->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)di, b_k, sizeof(double));
    loop_ub = y->size[1];
    for (b_k = 0; b_k < loop_ub; b_k++) {
      di->data[b_k] = y->data[y->size[0] * b_k];
    }

    if (!((x->size[0] == 0) || (x->size[1] == 0))) {
      xoffset = x->size[0];
    } else if (!(di->size[0] == 0)) {
      xoffset = di->size[0];
    } else {
      xoffset = x->size[0];
      if (!(xoffset > 0)) {
        xoffset = 0;
      }
    }

    empty_non_axis_sizes = (xoffset == 0);
    if (empty_non_axis_sizes || (!((x->size[0] == 0) || (x->size[1] == 0)))) {
      nrows = x->size[1];
    } else {
      nrows = 0;
    }

    if (empty_non_axis_sizes || (!(di->size[0] == 0))) {
      result = 1;
    } else {
      result = 0;
    }

    emxInit_real_T1(&Xremain, 2);
    b_k = Xremain->size[0] * Xremain->size[1];
    Xremain->size[0] = xoffset;
    Xremain->size[1] = nrows + result;
    emxEnsureCapacity((emxArray__common *)Xremain, b_k, sizeof(double));
    for (b_k = 0; b_k < nrows; b_k++) {
      for (i = 0; i < xoffset; i++) {
        Xremain->data[i + Xremain->size[0] * b_k] = x->data[i + xoffset * b_k];
      }
    }

    for (b_k = 0; b_k < result; b_k++) {
      for (i = 0; i < xoffset; i++) {
        Xremain->data[i + Xremain->size[0] * (b_k + nrows)] = di->data[i +
          xoffset * b_k];
      }
    }

    b_k = I->size[0] * I->size[1];
    I->size[0] = 0;
    I->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)I, b_k, sizeof(double));
    emxInit_int32_T(&Ii, 1);
    emxInit_int32_T(&b_Ii, 1);
    emxInit_int32_T(&iidx, 1);
    emxInit_real_T1(&b_x, 2);
    emxInit_boolean_T(&b, 2);
    emxInit_real_T1(&b_Xremain, 2);
    emxInit_real_T1(&b_y, 2);
    emxInit_real_T1(&c_Xremain, 2);
    emxInit_real_T1(&d_Xremain, 2);
    emxInit_real_T1(&b_I, 2);
    emxInit_real_T1(&e_Xremain, 2);
    while (Xremain->size[0] > 1) {
      if (1 > p) {
        loop_ub = 0;
        result = 0;
      } else {
        loop_ub = p;
        result = p;
      }

      b_k = d_Xremain->size[0] * d_Xremain->size[1];
      d_Xremain->size[0] = 1;
      d_Xremain->size[1] = loop_ub;
      emxEnsureCapacity((emxArray__common *)d_Xremain, b_k, sizeof(double));
      for (b_k = 0; b_k < loop_ub; b_k++) {
        d_Xremain->data[d_Xremain->size[0] * b_k] = Xremain->data[Xremain->size
          [0] * b_k];
      }

      repmat(d_Xremain, (double)Xremain->size[0], b_x);
      loop_ub = Xremain->size[0];
      b_k = c_Xremain->size[0] * c_Xremain->size[1];
      c_Xremain->size[0] = loop_ub;
      c_Xremain->size[1] = result;
      emxEnsureCapacity((emxArray__common *)c_Xremain, b_k, sizeof(double));
      for (b_k = 0; b_k < result; b_k++) {
        for (i = 0; i < loop_ub; i++) {
          c_Xremain->data[i + c_Xremain->size[0] * b_k] = Xremain->data[i +
            Xremain->size[0] * b_k] - b_x->data[i + b_x->size[0] * b_k];
        }
      }

      power(c_Xremain, b_x);
      sum(b_x, di);
      c_sort(di, iidx);
      b_k = b_Ii->size[0];
      b_Ii->size[0] = iidx->size[0];
      emxEnsureCapacity((emxArray__common *)b_Ii, b_k, sizeof(int));
      loop_ub = iidx->size[0];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        b_Ii->data[b_k] = iidx->data[b_k];
      }

      s = Xremain->size[0];
      if (k < s) {
        s = k;
      }

      if (1.0 > s) {
        loop_ub = 0;
      } else {
        loop_ub = (int)s;
      }

      b_k = Ii->size[0];
      Ii->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)Ii, b_k, sizeof(int));
      for (b_k = 0; b_k < loop_ub; b_k++) {
        Ii->data[b_k] = b_Ii->data[b_k];
      }

      if (1 > p) {
        result = 0;
      } else {
        result = p;
      }

      b_k = b_x->size[0] * b_x->size[1];
      b_x->size[0] = Ii->size[0];
      b_x->size[1] = result;
      emxEnsureCapacity((emxArray__common *)b_x, b_k, sizeof(double));
      for (b_k = 0; b_k < result; b_k++) {
        nrows = Ii->size[0];
        for (i = 0; i < nrows; i++) {
          b_x->data[i + b_x->size[0] * b_k] = Xremain->data[(Ii->data[i] +
            Xremain->size[0] * b_k) - 1];
        }
      }

      b_k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = result;
      emxEnsureCapacity((emxArray__common *)y, b_k, sizeof(double));
      if ((loop_ub == 0) || (result == 0)) {
        b_k = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, b_k, sizeof(double));
        result = y->size[1];
        for (b_k = 0; b_k < result; b_k++) {
          y->data[y->size[0] * b_k] = 0.0;
        }
      } else {
        for (i = 0; i + 1 <= result; i++) {
          xoffset = i * loop_ub;
          s = b_x->data[xoffset];
          for (b_k = 2; b_k <= loop_ub; b_k++) {
            s += b_x->data[(xoffset + b_k) - 1];
          }

          y->data[i] = s;
        }
      }

      if (1 > p) {
        result = 0;
      } else {
        result = p;
      }

      b_k = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = y->size[1];
      emxEnsureCapacity((emxArray__common *)b_y, b_k, sizeof(double));
      nrows = y->size[0] * y->size[1];
      for (b_k = 0; b_k < nrows; b_k++) {
        b_y->data[b_k] = y->data[b_k] / (double)loop_ub;
      }

      repmat(b_y, (double)loop_ub, b_x);
      b_k = b_Xremain->size[0] * b_Xremain->size[1];
      b_Xremain->size[0] = Ii->size[0];
      b_Xremain->size[1] = result;
      emxEnsureCapacity((emxArray__common *)b_Xremain, b_k, sizeof(double));
      for (b_k = 0; b_k < result; b_k++) {
        loop_ub = Ii->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_Xremain->data[i + b_Xremain->size[0] * b_k] = Xremain->data
            [(Ii->data[i] + Xremain->size[0] * b_k) - 1] - b_x->data[i +
            b_x->size[0] * b_k];
        }
      }

      power(b_Xremain, b_x);
      sum(b_x, di);
      c_sort(di, iidx);
      b_k = di->size[0];
      di->size[0] = iidx->size[0];
      emxEnsureCapacity((emxArray__common *)di, b_k, sizeof(double));
      loop_ub = iidx->size[0];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        di->data[b_k] = iidx->data[b_k];
      }

      s = Xremain->data[(b_Ii->data[(int)di->data[0] - 1] + Xremain->size[0] * p)
        - 1];
      if (!((I->size[0] == 0) || (I->size[1] == 0))) {
        nrows = I->size[1];
      } else {
        nrows = 0;
      }

      b_k = b_I->size[0] * b_I->size[1];
      b_I->size[0] = 1;
      b_I->size[1] = nrows + 1;
      emxEnsureCapacity((emxArray__common *)b_I, b_k, sizeof(double));
      for (b_k = 0; b_k < nrows; b_k++) {
        b_I->data[b_I->size[0] * b_k] = I->data[b_k];
      }

      b_I->data[b_I->size[0] * nrows] = s;
      b_k = I->size[0] * I->size[1];
      I->size[0] = b_I->size[0];
      I->size[1] = b_I->size[1];
      emxEnsureCapacity((emxArray__common *)I, b_k, sizeof(double));
      loop_ub = b_I->size[1];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        result = b_I->size[0];
        for (i = 0; i < result; i++) {
          I->data[i + I->size[0] * b_k] = b_I->data[i + b_I->size[0] * b_k];
        }
      }

      b_k = iidx->size[0];
      iidx->size[0] = Ii->size[0];
      emxEnsureCapacity((emxArray__common *)iidx, b_k, sizeof(int));
      loop_ub = Ii->size[0];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        iidx->data[b_k] = Ii->data[b_k];
      }

      nrowx = Xremain->size[0];
      ncolx = Xremain->size[1];
      if (iidx->size[0] == 1) {
        nrows = Xremain->size[0] - 1;
        for (result = 0; result + 1 <= ncolx; result++) {
          for (i = iidx->data[0]; i < nrowx; i++) {
            Xremain->data[(i + Xremain->size[0] * result) - 1] = Xremain->data[i
              + Xremain->size[0] * result];
          }
        }
      } else {
        b_k = b->size[0] * b->size[1];
        b->size[0] = 1;
        b->size[1] = Xremain->size[0];
        emxEnsureCapacity((emxArray__common *)b, b_k, sizeof(boolean_T));
        loop_ub = Xremain->size[0];
        for (b_k = 0; b_k < loop_ub; b_k++) {
          b->data[b_k] = false;
        }

        for (b_k = 1; b_k <= iidx->size[0]; b_k++) {
          b->data[iidx->data[b_k - 1] - 1] = true;
        }

        xoffset = 0;
        for (b_k = 1; b_k <= b->size[1]; b_k++) {
          xoffset += b->data[b_k - 1];
        }

        nrows = Xremain->size[0] - xoffset;
        i = 0;
        for (b_k = 1; b_k <= nrowx; b_k++) {
          if ((b_k > b->size[1]) || (!b->data[b_k - 1])) {
            for (result = 0; result + 1 <= ncolx; result++) {
              Xremain->data[i + Xremain->size[0] * result] = Xremain->data[(b_k
                + Xremain->size[0] * result) - 1];
            }

            i++;
          }
        }
      }

      if (1 > nrows) {
        loop_ub = 0;
      } else {
        loop_ub = nrows;
      }

      xoffset = Xremain->size[1];
      b_k = e_Xremain->size[0] * e_Xremain->size[1];
      e_Xremain->size[0] = loop_ub;
      e_Xremain->size[1] = xoffset;
      emxEnsureCapacity((emxArray__common *)e_Xremain, b_k, sizeof(double));
      for (b_k = 0; b_k < xoffset; b_k++) {
        for (i = 0; i < loop_ub; i++) {
          e_Xremain->data[i + e_Xremain->size[0] * b_k] = Xremain->data[i +
            Xremain->size[0] * b_k];
        }
      }

      b_k = Xremain->size[0] * Xremain->size[1];
      Xremain->size[0] = e_Xremain->size[0];
      Xremain->size[1] = e_Xremain->size[1];
      emxEnsureCapacity((emxArray__common *)Xremain, b_k, sizeof(double));
      loop_ub = e_Xremain->size[1];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        result = e_Xremain->size[0];
        for (i = 0; i < result; i++) {
          Xremain->data[i + Xremain->size[0] * b_k] = e_Xremain->data[i +
            e_Xremain->size[0] * b_k];
        }
      }
    }

    emxFree_real_T(&e_Xremain);
    emxFree_real_T(&b_I);
    emxFree_real_T(&d_Xremain);
    emxFree_real_T(&c_Xremain);
    emxFree_real_T(&b_y);
    emxFree_real_T(&b_Xremain);
    emxFree_boolean_T(&b);
    emxFree_real_T(&b_x);
    emxFree_int32_T(&iidx);
    emxFree_int32_T(&b_Ii);
    emxFree_int32_T(&Ii);
    emxFree_real_T(&di);
    emxFree_real_T(&Xremain);
  }

  emxFree_real_T(&y);
}

void MAVEfast(emxArray_real_T *x, const emxArray_real_T *y, emxArray_char_T
              *method, double max_dim, double screen, emxArray_real_T *BB,
              emxArray_real_T *ky, emxArray_real_T *BBvs, emxArray_real_T *idx,
              emxArray_real_T *C)
{
  int p0;
  int n;
  double m;
  emxArray_real_T *sele;
  int i13;
  int i14;
  double p;
  int loop_ub;
  emxArray_real_T *b_x;
  int b_idx;
  int br;
  emxArray_real_T *meanky;
  emxArray_real_T *a;
  emxArray_real_T *b_y;
  int k;
  int vstride;
  int b_m;
  emxArray_real_T *KY;
  emxArray_real_T *c_y;
  double pp;
  int ic;
  emxArray_real_T *V;
  emxArray_creal_T *Vc;
  emxArray_creal_T *Dc;
  int ar;
  int ib;
  int ia;
  emxArray_real_T *b_Dc;
  emxArray_real_T *D;
  emxArray_real_T *b_D;
  emxArray_boolean_T *c_x;
  int trueCount;
  emxArray_int32_T *iidx;
  emxArray_real_T *b;
  emxArray_real_T *ss;
  emxArray_real_T *d_x;
  emxArray_char_T *b_method;
  emxArray_real_T *ky1;
  emxArray_real_T *ky2;
  emxArray_real_T *DD;
  emxArray_real_T *qx;
  emxArray_real_T *yi;
  emxArray_real_T *yj;
  emxArray_real_T *U;
  emxArray_real_T *dxij;
  emxArray_int32_T *r1;
  emxArray_boolean_T *r2;
  emxArray_real_T *b_a;
  emxArray_int32_T *c_idx;
  cell_wrap_0 reshapes[2];
  emxArray_real_T *varargin_2;
  cell_wrap_0 b_reshapes[2];
  emxArray_real_T *c_D;
  emxArray_real_T *c_Dc;
  emxArray_boolean_T *b_yi;
  emxArray_real_T *b_meanky;
  emxArray_real_T *b_V;
  emxArray_real_T *c_V;
  emxArray_real_T *b_ky;
  emxArray_real_T *b_KY;
  double e_x;
  emxArray_real_T *r3;
  int exitg1;
  int j;
  emxArray_real_T *onexi;
  int exponent;
  boolean_T empty_non_axis_sizes;
  emxArray_real_T *B;
  emxArray_real_T *SEQ;
  boolean_T exitg2;
  int iter_ip;
  emxArray_real_T *b_BB;
  emxArray_real_T *Ifast;
  emxArray_real_T *Vfast;
  static const double c_a[51] = { 1.0, 0.8, 0.64000000000000012,
    0.51200000000000012, 0.40960000000000008, 0.32768000000000008,
    0.2621440000000001, 0.20971520000000007, 0.16777216000000006,
    0.13421772800000006, 0.10737418240000006, 0.08589934592000005,
    0.06871947673600004, 0.054975581388800043, 0.043980465111040035,
    0.03518437208883203, 0.028147497671065624, 0.0225179981368525,
    0.018014398509482003, 0.014411518807585602, 0.011529215046068483,
    0.0092233720368547871, 0.00737869762948383, 0.0059029581035870641,
    0.0047223664828696518, 0.0037778931862957215, 0.0030223145490365774,
    0.0024178516392292619, 0.0019342813113834097, 0.0015474250491067279,
    0.0012379400392853823, 0.00099035203142830582, 0.00079228162514264483,
    0.00063382530011411582, 0.0005070602400912927, 0.00040564819207303422,
    0.00032451855365842736, 0.00025961484292674189, 0.00020769187434139353,
    0.00016615349947311485, 0.00013292279957849188, 0.00010633823966279351,
    8.5070591730234822E-5, 6.8056473384187861E-5, 5.4445178707350287E-5,
    4.3556142965880231E-5, 3.4844914372704188E-5, 2.7875931498163353E-5,
    2.2300745198530684E-5, 1.7840596158824548E-5, 1.4272476927059638E-5 };

  int ii_data[1];
  emxArray_real_T *xfast;
  emxArray_real_T *xij;
  emxArray_real_T *xk;
  double Ip[51];
  int di_data[1];
  emxArray_real_T *abi;
  signed char tmp_data[51];
  emxArray_real_T *dd;
  emxArray_real_T *dc;
  emxArray_real_T *d_D;
  double ndbl;
  emxArray_creal_T *R;
  double apnd;
  emxArray_int32_T *r4;
  cell_wrap_0 c_reshapes[2];
  int tmp_size[2];
  emxArray_real_T *b_dd;
  emxArray_real_T *b_abi;
  emxArray_real_T *d_a;
  emxArray_real_T *d_V;
  double b_tmp_data[51];
  emxArray_real_T *b_Vfast;
  emxArray_real_T *b_xfast;
  emxArray_real_T *e_a;
  emxArray_real_T *e_V;
  emxArray_real_T *b_B;
  emxArray_int32_T *b_Ifast;
  emxArray_int32_T *c_Ifast;
  emxArray_real_T *b_dc;
  emxArray_real_T *c_abi;
  emxArray_int32_T *d_Ifast;
  emxArray_int32_T *e_Ifast;
  emxArray_real_T *BB0;
  double ip;
  int K;
  emxArray_real_T *c_B;
  emxArray_real_T *d_B;
  boolean_T guard1 = false;
  int iter;
  double d_idx;
  double h2;
  double cdiff;
  double absa;
  double absb;
  unsigned int U_idx_0;

  /*  */
  p0 = x->size[1];
  n = x->size[0];
  m = screen;
  emxInit_real_T1(&sele, 2);
  if (x->size[1] > screen) {
    dcorrVS(x, y, sele);
  } else {
    if (x->size[1] < 1) {
      i13 = sele->size[0] * sele->size[1];
      sele->size[0] = 1;
      sele->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)sele, i13, sizeof(double));
    } else {
      i13 = x->size[1];
      i14 = sele->size[0] * sele->size[1];
      sele->size[0] = 1;
      sele->size[1] = (int)((double)i13 - 1.0) + 1;
      emxEnsureCapacity((emxArray__common *)sele, i14, sizeof(double));
      loop_ub = (int)((double)i13 - 1.0);
      for (i13 = 0; i13 <= loop_ub; i13++) {
        sele->data[sele->size[0] * i13] = 1.0 + (double)i13;
      }
    }

    m = x->size[1];
  }

  p = x->size[1];
  if (m < p) {
    p = m;
  }

  if (1.0 > m) {
    loop_ub = 0;
  } else {
    loop_ub = (int)m;
  }

  emxInit_real_T1(&b_x, 2);
  b_idx = x->size[0];
  i13 = b_x->size[0] * b_x->size[1];
  b_x->size[0] = b_idx;
  b_x->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)b_x, i13, sizeof(double));
  for (i13 = 0; i13 < loop_ub; i13++) {
    for (i14 = 0; i14 < b_idx; i14++) {
      b_x->data[i14 + b_x->size[0] * i13] = x->data[i14 + x->size[0] * ((int)
        sele->data[i13] - 1)];
    }
  }

  i13 = x->size[0] * x->size[1];
  x->size[0] = b_x->size[0];
  x->size[1] = b_x->size[1];
  emxEnsureCapacity((emxArray__common *)x, i13, sizeof(double));
  loop_ub = b_x->size[1];
  for (i13 = 0; i13 < loop_ub; i13++) {
    br = b_x->size[0];
    for (i14 = 0; i14 < br; i14++) {
      x->data[i14 + x->size[0] * i13] = b_x->data[i14 + b_x->size[0] * i13];
    }
  }

  emxFree_real_T(&b_x);
  if (1.0 > m) {
    loop_ub = 0;
  } else {
    loop_ub = (int)m;
  }

  i13 = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)idx, i13, sizeof(double));
  for (i13 = 0; i13 < loop_ub; i13++) {
    idx->data[idx->size[0] * i13] = sele->data[i13];
  }

  emxInit_real_T1(&meanky, 2);
  emxInit_real_T1(&a, 2);
  mean(x, meanky);
  repmat(meanky, (double)n, a);
  i13 = x->size[0] * x->size[1];
  emxEnsureCapacity((emxArray__common *)x, i13, sizeof(double));
  loop_ub = x->size[1];
  for (i13 = 0; i13 < loop_ub; i13++) {
    br = x->size[0];
    for (i14 = 0; i14 < br; i14++) {
      x->data[i14 + x->size[0] * i13] -= a->data[i14 + a->size[0] * i13];
    }
  }

  emxInit_real_T1(&b_y, 2);
  i13 = a->size[0] * a->size[1];
  a->size[0] = x->size[1];
  a->size[1] = x->size[0];
  emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
  loop_ub = x->size[0];
  for (i13 = 0; i13 < loop_ub; i13++) {
    br = x->size[1];
    for (i14 = 0; i14 < br; i14++) {
      a->data[i14 + a->size[0] * i13] = x->data[i13 + x->size[0] * i14];
    }
  }

  if ((a->size[1] == 1) || (x->size[0] == 1)) {
    i13 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = a->size[0];
    b_y->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
    loop_ub = a->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = x->size[1];
      for (i14 = 0; i14 < br; i14++) {
        b_y->data[i13 + b_y->size[0] * i14] = 0.0;
        vstride = a->size[1];
        for (b_m = 0; b_m < vstride; b_m++) {
          b_y->data[i13 + b_y->size[0] * i14] += a->data[i13 + a->size[0] * b_m]
            * x->data[b_m + x->size[0] * i14];
        }
      }
    }
  } else {
    k = a->size[1];
    vstride = a->size[0];
    b_idx = x->size[1];
    i13 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = vstride;
    b_y->size[1] = b_idx;
    emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
    b_m = a->size[0];
    i13 = b_y->size[0] * b_y->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
    loop_ub = b_y->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = b_y->size[0];
      for (i14 = 0; i14 < br; i14++) {
        b_y->data[i14 + b_y->size[0] * i13] = 0.0;
      }
    }

    if ((a->size[0] == 0) || (x->size[1] == 0)) {
    } else {
      b_idx = a->size[0] * (x->size[1] - 1);
      vstride = 0;
      while ((b_m > 0) && (vstride <= b_idx)) {
        i13 = vstride + b_m;
        for (ic = vstride; ic + 1 <= i13; ic++) {
          b_y->data[ic] = 0.0;
        }

        vstride += b_m;
      }

      br = 0;
      vstride = 0;
      while ((b_m > 0) && (vstride <= b_idx)) {
        ar = 0;
        i13 = br + k;
        for (ib = br; ib + 1 <= i13; ib++) {
          if (x->data[ib] != 0.0) {
            ia = ar;
            i14 = vstride + b_m;
            for (ic = vstride; ic + 1 <= i14; ic++) {
              ia++;
              b_y->data[ic] += x->data[ib] * a->data[ia - 1];
            }
          }

          ar += b_m;
        }

        br += k;
        vstride += b_m;
      }
    }
  }

  emxInit_real_T1(&KY, 2);
  emxInit_real_T1(&c_y, 2);
  pp = rt_powd_snf((double)n, 5.0);
  eye(p, KY);
  i13 = c_y->size[0] * c_y->size[1];
  c_y->size[0] = b_y->size[0];
  c_y->size[1] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)c_y, i13, sizeof(double));
  loop_ub = b_y->size[0] * b_y->size[1];
  for (i13 = 0; i13 < loop_ub; i13++) {
    c_y->data[i13] = b_y->data[i13] / (double)n + KY->data[i13] / pp;
  }

  emxInit_real_T1(&V, 2);
  emxInit_creal_T(&Vc, 2);
  emxInit_creal_T(&Dc, 2);
  eig(c_y, Vc, Dc);
  i13 = V->size[0] * V->size[1];
  V->size[0] = Vc->size[0];
  V->size[1] = Vc->size[1];
  emxEnsureCapacity((emxArray__common *)V, i13, sizeof(double));
  loop_ub = Vc->size[0] * Vc->size[1];
  emxFree_real_T(&c_y);
  for (i13 = 0; i13 < loop_ub; i13++) {
    V->data[i13] = Vc->data[i13].re;
  }

  emxInit_real_T1(&b_Dc, 2);
  i13 = b_Dc->size[0] * b_Dc->size[1];
  b_Dc->size[0] = Dc->size[0];
  b_Dc->size[1] = Dc->size[1];
  emxEnsureCapacity((emxArray__common *)b_Dc, i13, sizeof(double));
  loop_ub = Dc->size[0] * Dc->size[1];
  for (i13 = 0; i13 < loop_ub; i13++) {
    b_Dc->data[i13] = Dc->data[i13].re;
  }

  emxInit_real_T(&D, 1);
  emxInit_real_T(&b_D, 1);
  diag(b_Dc, D);
  i13 = b_D->size[0];
  b_D->size[0] = D->size[0];
  emxEnsureCapacity((emxArray__common *)b_D, i13, sizeof(double));
  loop_ub = D->size[0];
  emxFree_real_T(&b_Dc);
  for (i13 = 0; i13 < loop_ub; i13++) {
    b_D->data[i13] = D->data[i13];
  }

  emxInit_boolean_T1(&c_x, 1);
  rdivide(b_D, D);
  pp = rt_powd_snf((double)n, 5.0);
  i13 = c_x->size[0];
  c_x->size[0] = D->size[0];
  emxEnsureCapacity((emxArray__common *)c_x, i13, sizeof(boolean_T));
  loop_ub = D->size[0];
  emxFree_real_T(&b_D);
  for (i13 = 0; i13 < loop_ub; i13++) {
    c_x->data[i13] = (D->data[i13] >= pp);
  }

  vstride = c_x->size[0] - 1;
  trueCount = 0;
  for (ar = 0; ar <= vstride; ar++) {
    if (c_x->data[ar]) {
      trueCount++;
    }
  }

  emxInit_int32_T(&iidx, 1);
  i13 = iidx->size[0];
  iidx->size[0] = trueCount;
  emxEnsureCapacity((emxArray__common *)iidx, i13, sizeof(int));
  b_idx = 0;
  for (ar = 0; ar <= vstride; ar++) {
    if (c_x->data[ar]) {
      iidx->data[b_idx] = ar + 1;
      b_idx++;
    }
  }

  loop_ub = iidx->size[0];
  for (i13 = 0; i13 < loop_ub; i13++) {
    D->data[iidx->data[i13] - 1] = 0.0;
  }

  emxInit_real_T1(&b, 2);
  d_sqrt(D);
  b_diag(D, b);
  if ((V->size[1] == 1) || (b->size[0] == 1)) {
    i13 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = V->size[0];
    b_y->size[1] = b->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
    loop_ub = V->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = b->size[1];
      for (i14 = 0; i14 < br; i14++) {
        b_y->data[i13 + b_y->size[0] * i14] = 0.0;
        vstride = V->size[1];
        for (b_m = 0; b_m < vstride; b_m++) {
          b_y->data[i13 + b_y->size[0] * i14] += V->data[i13 + V->size[0] * b_m]
            * b->data[b_m + b->size[0] * i14];
        }
      }
    }
  } else {
    k = V->size[1];
    vstride = V->size[0];
    b_idx = b->size[1];
    i13 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = vstride;
    b_y->size[1] = b_idx;
    emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
    b_m = V->size[0];
    i13 = b_y->size[0] * b_y->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
    loop_ub = b_y->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = b_y->size[0];
      for (i14 = 0; i14 < br; i14++) {
        b_y->data[i14 + b_y->size[0] * i13] = 0.0;
      }
    }

    if ((V->size[0] == 0) || (b->size[1] == 0)) {
    } else {
      b_idx = V->size[0] * (b->size[1] - 1);
      vstride = 0;
      while ((b_m > 0) && (vstride <= b_idx)) {
        i13 = vstride + b_m;
        for (ic = vstride; ic + 1 <= i13; ic++) {
          b_y->data[ic] = 0.0;
        }

        vstride += b_m;
      }

      br = 0;
      vstride = 0;
      while ((b_m > 0) && (vstride <= b_idx)) {
        ar = 0;
        i13 = br + k;
        for (ib = br; ib + 1 <= i13; ib++) {
          if (b->data[ib] != 0.0) {
            ia = ar;
            i14 = vstride + b_m;
            for (ic = vstride; ic + 1 <= i14; ic++) {
              ia++;
              b_y->data[ic] += b->data[ib] * V->data[ia - 1];
            }
          }

          ar += b_m;
        }

        br += k;
        vstride += b_m;
      }
    }
  }

  i13 = b->size[0] * b->size[1];
  b->size[0] = V->size[1];
  b->size[1] = V->size[0];
  emxEnsureCapacity((emxArray__common *)b, i13, sizeof(double));
  loop_ub = V->size[0];
  for (i13 = 0; i13 < loop_ub; i13++) {
    br = V->size[1];
    for (i14 = 0; i14 < br; i14++) {
      b->data[i14 + b->size[0] * i13] = V->data[i13 + V->size[0] * i14];
    }
  }

  emxInit_real_T1(&ss, 2);
  if ((b_y->size[1] == 1) || (b->size[0] == 1)) {
    i13 = ss->size[0] * ss->size[1];
    ss->size[0] = b_y->size[0];
    ss->size[1] = b->size[1];
    emxEnsureCapacity((emxArray__common *)ss, i13, sizeof(double));
    loop_ub = b_y->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = b->size[1];
      for (i14 = 0; i14 < br; i14++) {
        ss->data[i13 + ss->size[0] * i14] = 0.0;
        vstride = b_y->size[1];
        for (b_m = 0; b_m < vstride; b_m++) {
          ss->data[i13 + ss->size[0] * i14] += b_y->data[i13 + b_y->size[0] *
            b_m] * b->data[b_m + b->size[0] * i14];
        }
      }
    }
  } else {
    k = b_y->size[1];
    vstride = b_y->size[0];
    b_idx = b->size[1];
    i13 = ss->size[0] * ss->size[1];
    ss->size[0] = vstride;
    ss->size[1] = b_idx;
    emxEnsureCapacity((emxArray__common *)ss, i13, sizeof(double));
    b_m = b_y->size[0];
    i13 = ss->size[0] * ss->size[1];
    emxEnsureCapacity((emxArray__common *)ss, i13, sizeof(double));
    loop_ub = ss->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = ss->size[0];
      for (i14 = 0; i14 < br; i14++) {
        ss->data[i14 + ss->size[0] * i13] = 0.0;
      }
    }

    if ((b_y->size[0] == 0) || (b->size[1] == 0)) {
    } else {
      b_idx = b_y->size[0] * (b->size[1] - 1);
      vstride = 0;
      while ((b_m > 0) && (vstride <= b_idx)) {
        i13 = vstride + b_m;
        for (ic = vstride; ic + 1 <= i13; ic++) {
          ss->data[ic] = 0.0;
        }

        vstride += b_m;
      }

      br = 0;
      vstride = 0;
      while ((b_m > 0) && (vstride <= b_idx)) {
        ar = 0;
        i13 = br + k;
        for (ib = br; ib + 1 <= i13; ib++) {
          if (b->data[ib] != 0.0) {
            ia = ar;
            i14 = vstride + b_m;
            for (ic = vstride; ic + 1 <= i14; ic++) {
              ia++;
              ss->data[ic] += b->data[ib] * b_y->data[ia - 1];
            }
          }

          ar += b_m;
        }

        br += k;
        vstride += b_m;
      }
    }
  }

  i13 = a->size[0] * a->size[1];
  a->size[0] = x->size[0];
  a->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i13 = 0; i13 < loop_ub; i13++) {
    a->data[i13] = x->data[i13];
  }

  emxInit_real_T1(&d_x, 2);
  if ((x->size[1] == 1) || (ss->size[0] == 1)) {
    i13 = d_x->size[0] * d_x->size[1];
    d_x->size[0] = x->size[0];
    d_x->size[1] = ss->size[1];
    emxEnsureCapacity((emxArray__common *)d_x, i13, sizeof(double));
    loop_ub = x->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = ss->size[1];
      for (i14 = 0; i14 < br; i14++) {
        d_x->data[i13 + d_x->size[0] * i14] = 0.0;
        vstride = x->size[1];
        for (b_m = 0; b_m < vstride; b_m++) {
          d_x->data[i13 + d_x->size[0] * i14] += x->data[i13 + x->size[0] * b_m]
            * ss->data[b_m + ss->size[0] * i14];
        }
      }
    }

    i13 = x->size[0] * x->size[1];
    x->size[0] = d_x->size[0];
    x->size[1] = d_x->size[1];
    emxEnsureCapacity((emxArray__common *)x, i13, sizeof(double));
    loop_ub = d_x->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = d_x->size[0];
      for (i14 = 0; i14 < br; i14++) {
        x->data[i14 + x->size[0] * i13] = d_x->data[i14 + d_x->size[0] * i13];
      }
    }
  } else {
    k = x->size[1];
    vstride = x->size[0];
    b_idx = ss->size[1];
    i13 = x->size[0] * x->size[1];
    x->size[0] = vstride;
    x->size[1] = b_idx;
    emxEnsureCapacity((emxArray__common *)x, i13, sizeof(double));
    b_m = a->size[0];
    i13 = x->size[0] * x->size[1];
    emxEnsureCapacity((emxArray__common *)x, i13, sizeof(double));
    loop_ub = x->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = x->size[0];
      for (i14 = 0; i14 < br; i14++) {
        x->data[i14 + x->size[0] * i13] = 0.0;
      }
    }

    if ((a->size[0] == 0) || (ss->size[1] == 0)) {
    } else {
      b_idx = a->size[0] * (ss->size[1] - 1);
      vstride = 0;
      while ((b_m > 0) && (vstride <= b_idx)) {
        i13 = vstride + b_m;
        for (ic = vstride; ic + 1 <= i13; ic++) {
          x->data[ic] = 0.0;
        }

        vstride += b_m;
      }

      br = 0;
      vstride = 0;
      while ((b_m > 0) && (vstride <= b_idx)) {
        ar = 0;
        i13 = br + k;
        for (ib = br; ib + 1 <= i13; ib++) {
          if (ss->data[ib] != 0.0) {
            ia = ar;
            i14 = vstride + b_m;
            for (ic = vstride; ic + 1 <= i14; ic++) {
              ia++;
              x->data[ic] += ss->data[ib] * a->data[ia - 1];
            }
          }

          ar += b_m;
        }

        br += k;
        vstride += b_m;
      }
    }
  }

  emxFree_real_T(&d_x);
  emxInit_char_T(&b_method, 2);
  i13 = b_method->size[0] * b_method->size[1];
  b_method->size[0] = method->size[0];
  b_method->size[1] = method->size[1];
  emxEnsureCapacity((emxArray__common *)b_method, i13, sizeof(char));
  loop_ub = method->size[0] * method->size[1];
  for (i13 = 0; i13 < loop_ub; i13++) {
    b_method->data[i13] = method->data[i13];
  }

  upper(b_method, method);

  /* which_dim=sort(which_dim,'descend'); */
  /* if(which_dim(1)~=p) */
  /*     which_dim=[p,which_dim]; */
  /* end */
  /* len_which_dim=length(which_dim); */
  emxFree_char_T(&b_method);
  emxInit_real_T1(&ky1, 2);
  emxInit_real_T1(&ky2, 2);
  emxInit_real_T1(&DD, 2);
  emxInit_real_T1(&qx, 2);
  emxInit_real_T1(&yi, 2);
  emxInit_real_T1(&yj, 2);
  emxInit_real_T1(&U, 2);
  emxInit_real_T(&dxij, 1);
  emxInit_int32_T(&r1, 1);
  emxInit_boolean_T(&r2, 2);
  emxInit_real_T1(&b_a, 2);
  emxInit_int32_T1(&c_idx, 2);
  emxInitMatrix_cell_wrap_0(reshapes);
  emxInit_real_T1(&varargin_2, 2);
  emxInitMatrix_cell_wrap_0(b_reshapes);
  emxInit_real_T1(&c_D, 2);
  emxInit_real_T1(&c_Dc, 2);
  emxInit_boolean_T(&b_yi, 2);
  emxInit_real_T1(&b_meanky, 2);
  emxInit_real_T1(&b_V, 2);
  emxInit_real_T1(&c_V, 2);
  emxInit_real_T1(&b_ky, 2);
  emxInit_real_T1(&b_KY, 2);
  if (b_strcmp(method) || c_strcmp(method)) {
    i13 = ky->size[0] * ky->size[1];
    ky->size[0] = y->size[0];
    ky->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)ky, i13, sizeof(double));
    loop_ub = y->size[0] * y->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      ky->data[i13] = y->data[i13];
    }

    i13 = ky1->size[0] * ky1->size[1];
    ky1->size[0] = y->size[0];
    ky1->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)ky1, i13, sizeof(double));
    loop_ub = y->size[0] * y->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      ky1->data[i13] = y->data[i13];
    }

    i13 = ky2->size[0] * ky2->size[1];
    ky2->size[0] = 1;
    ky2->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)ky2, i13, sizeof(double));
    ky2->data[0] = 1.0;
    i13 = DD->size[0] * DD->size[1];
    DD->size[0] = 1;
    DD->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)DD, i13, sizeof(double));
    DD->data[0] = 1.0;
    i13 = C->size[0] * C->size[1];
    C->size[0] = 0;
    C->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)C, i13, sizeof(double));
  } else {
    i13 = KY->size[0] * KY->size[1];
    KY->size[0] = 0;
    KY->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)KY, i13, sizeof(double));
    for (k = 0; k < y->size[1]; k++) {
      loop_ub = y->size[0];
      i13 = D->size[0];
      D->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)D, i13, sizeof(double));
      for (i13 = 0; i13 < loop_ub; i13++) {
        D->data[i13] = y->data[i13 + y->size[0] * k];
      }

      if (99 < n) {
        b_idx = 99;
      } else {
        b_idx = n;
      }

      e_x = std::floor(rt_powd_snf((double)n, 0.6));
      if ((b_idx > e_x) || rtIsNaN(e_x)) {
        pp = b_idx;
      } else {
        pp = e_x;
      }

      if (rtIsNaN(pp)) {
        i13 = meanky->size[0] * meanky->size[1];
        meanky->size[0] = 1;
        meanky->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
        meanky->data[0] = rtNaN;
      } else if (pp < 1.0) {
        i13 = meanky->size[0] * meanky->size[1];
        meanky->size[0] = 1;
        meanky->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
      } else if (rtIsInf(pp) && (1.0 == pp)) {
        i13 = meanky->size[0] * meanky->size[1];
        meanky->size[0] = 1;
        meanky->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
        meanky->data[0] = rtNaN;
      } else {
        i13 = meanky->size[0] * meanky->size[1];
        meanky->size[0] = 1;
        meanky->size[1] = (int)(pp - 1.0) + 1;
        emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
        loop_ub = (int)(pp - 1.0);
        for (i13 = 0; i13 <= loop_ub; i13++) {
          meanky->data[meanky->size[0] * i13] = 1.0 + (double)i13;
        }
      }

      i13 = b_meanky->size[0] * b_meanky->size[1];
      b_meanky->size[0] = 1;
      b_meanky->size[1] = meanky->size[1];
      emxEnsureCapacity((emxArray__common *)b_meanky, i13, sizeof(double));
      loop_ub = meanky->size[0] * meanky->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        b_meanky->data[i13] = meanky->data[i13] / (pp + 1.0);
      }

      quantile(D, b_meanky, qx);
      i13 = b_a->size[0] * b_a->size[1];
      b_a->size[0] = 1;
      b_a->size[1] = qx->size[1];
      emxEnsureCapacity((emxArray__common *)b_a, i13, sizeof(double));
      loop_ub = qx->size[0] * qx->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        b_a->data[i13] = qx->data[i13];
      }

      vstride = qx->size[1];
      sortIdx(qx, c_idx);
      b_idx = qx->size[1];
      i13 = qx->size[0] * qx->size[1];
      qx->size[0] = 1;
      qx->size[1] = b_idx;
      emxEnsureCapacity((emxArray__common *)qx, i13, sizeof(double));
      for (br = 0; br + 1 <= vstride; br++) {
        qx->data[br] = b_a->data[c_idx->data[br] - 1];
      }

      count_nonfinites(qx, vstride, &br, &b_idx, &trueCount, &ar);
      vstride = -1;
      if (br > 0) {
        vstride = 0;
      }

      b_idx += br;
      while (br + 1 <= b_idx) {
        e_x = qx->data[br];
        do {
          exitg1 = 0;
          br++;
          if (br + 1 > b_idx) {
            exitg1 = 1;
          } else {
            pp = std::abs(e_x / 2.0);
            if ((!rtIsInf(pp)) && (!rtIsNaN(pp))) {
              if (pp <= 2.2250738585072014E-308) {
                pp = 4.94065645841247E-324;
              } else {
                frexp(pp, &exponent);
                pp = std::ldexp(1.0, exponent - 53);
              }
            } else {
              pp = rtNaN;
            }

            if ((std::abs(e_x - qx->data[br]) < pp) || (rtIsInf(qx->data[br]) &&
                 rtIsInf(e_x) && ((qx->data[br] > 0.0) == (e_x > 0.0)))) {
              empty_non_axis_sizes = true;
            } else {
              empty_non_axis_sizes = false;
            }

            if (!empty_non_axis_sizes) {
              exitg1 = 1;
            }
          }
        } while (exitg1 == 0);

        vstride++;
        qx->data[vstride] = e_x;
      }

      if (trueCount > 0) {
        vstride++;
        qx->data[vstride] = qx->data[b_idx];
      }

      br = b_idx + trueCount;
      for (j = 1; j <= ar; j++) {
        vstride++;
        qx->data[vstride] = qx->data[(br + j) - 1];
      }

      if (1 > vstride + 1) {
        i13 = 0;
      } else {
        i13 = vstride + 1;
      }

      i14 = qx->size[0] * qx->size[1];
      qx->size[1] = i13;
      emxEnsureCapacity((emxArray__common *)qx, i14, sizeof(double));
      b_repmat(D, (double)i13, a);
      i13 = yi->size[0] * yi->size[1];
      yi->size[0] = a->size[0];
      yi->size[1] = a->size[1];
      emxEnsureCapacity((emxArray__common *)yi, i13, sizeof(double));
      loop_ub = a->size[0] * a->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        yi->data[i13] = a->data[i13];
      }

      repmat(qx, (double)n, a);
      i13 = yj->size[0] * yj->size[1];
      yj->size[0] = a->size[0];
      yj->size[1] = a->size[1];
      emxEnsureCapacity((emxArray__common *)yj, i13, sizeof(double));
      loop_ub = a->size[0] * a->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        yj->data[i13] = a->data[i13];
      }

      i13 = b_yi->size[0] * b_yi->size[1];
      b_yi->size[0] = yi->size[0];
      b_yi->size[1] = yi->size[1];
      emxEnsureCapacity((emxArray__common *)b_yi, i13, sizeof(boolean_T));
      loop_ub = yi->size[0] * yi->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        b_yi->data[i13] = (yi->data[i13] - yj->data[i13] < 0.0);
      }

      b_abs(b_yi, yi);
      i13 = ky->size[0] * ky->size[1];
      ky->size[0] = yi->size[0];
      ky->size[1] = yi->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i13, sizeof(double));
      loop_ub = yi->size[0] * yi->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        ky->data[i13] = yi->data[i13];
      }

      mean(ky, meanky);
      repmat(meanky, (double)n, a);
      i13 = ky->size[0] * ky->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i13, sizeof(double));
      b_idx = ky->size[0];
      vstride = ky->size[1];
      loop_ub = b_idx * vstride;
      for (i13 = 0; i13 < loop_ub; i13++) {
        ky->data[i13] -= a->data[i13];
      }

      b_mean(ky, D);
      b_repmat(D, (double)ky->size[1], a);
      i13 = yi->size[0] * yi->size[1];
      yi->size[0] = a->size[0];
      yi->size[1] = a->size[1];
      emxEnsureCapacity((emxArray__common *)yi, i13, sizeof(double));
      loop_ub = a->size[0] * a->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        yi->data[i13] = a->data[i13];
      }

      i13 = ky1->size[0] * ky1->size[1];
      ky1->size[0] = ky->size[0];
      ky1->size[1] = ky->size[1];
      emxEnsureCapacity((emxArray__common *)ky1, i13, sizeof(double));
      loop_ub = ky->size[0] * ky->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        ky1->data[i13] = ky->data[i13] - yi->data[i13];
      }

      mean(ky1, meanky);
      repmat(meanky, (double)n, a);
      i13 = ky1->size[0] * ky1->size[1];
      emxEnsureCapacity((emxArray__common *)ky1, i13, sizeof(double));
      b_idx = ky1->size[0];
      vstride = ky1->size[1];
      loop_ub = b_idx * vstride;
      for (i13 = 0; i13 < loop_ub; i13++) {
        ky1->data[i13] -= a->data[i13];
      }

      if (!((ky->size[0] == 0) || (ky->size[1] == 0))) {
        trueCount = ky->size[0];
      } else if (!((ky1->size[0] == 0) || (ky1->size[1] == 0))) {
        trueCount = ky1->size[0];
      } else {
        trueCount = ky->size[0];
        if (!(trueCount > 0)) {
          trueCount = 0;
        }

        if (ky1->size[0] > trueCount) {
          trueCount = ky1->size[0];
        }
      }

      empty_non_axis_sizes = (trueCount == 0);
      if (empty_non_axis_sizes || (!((ky->size[0] == 0) || (ky->size[1] == 0))))
      {
        b_idx = ky->size[1];
      } else {
        b_idx = 0;
      }

      if (empty_non_axis_sizes || (!((ky1->size[0] == 0) || (ky1->size[1] == 0))))
      {
        vstride = ky1->size[1];
      } else {
        vstride = 0;
      }

      i13 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
      reshapes[1].f1->size[0] = trueCount;
      reshapes[1].f1->size[1] = vstride;
      emxEnsureCapacity((emxArray__common *)reshapes[1].f1, i13, sizeof(double));
      loop_ub = trueCount * vstride;
      for (i13 = 0; i13 < loop_ub; i13++) {
        reshapes[1].f1->data[i13] = ky1->data[i13];
      }

      i13 = U->size[0] * U->size[1];
      U->size[0] = trueCount;
      U->size[1] = b_idx;
      emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
      loop_ub = trueCount * b_idx;
      for (i13 = 0; i13 < loop_ub; i13++) {
        U->data[i13] = ky->data[i13];
      }

      i13 = b_ky->size[0] * b_ky->size[1];
      b_ky->size[0] = trueCount;
      b_ky->size[1] = b_idx + reshapes[1].f1->size[1];
      emxEnsureCapacity((emxArray__common *)b_ky, i13, sizeof(double));
      for (i13 = 0; i13 < b_idx; i13++) {
        for (i14 = 0; i14 < trueCount; i14++) {
          b_ky->data[i14 + b_ky->size[0] * i13] = ky->data[i14 + trueCount * i13];
        }
      }

      loop_ub = reshapes[1].f1->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = reshapes[1].f1->size[0];
        for (i14 = 0; i14 < br; i14++) {
          b_ky->data[i14 + b_ky->size[0] * (i13 + b_idx)] = reshapes[1].f1->
            data[i14 + reshapes[1].f1->size[0] * i13];
        }
      }

      i13 = ky->size[0] * ky->size[1];
      ky->size[0] = b_ky->size[0];
      ky->size[1] = b_ky->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i13, sizeof(double));
      loop_ub = b_ky->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = b_ky->size[0];
        for (i14 = 0; i14 < br; i14++) {
          ky->data[i14 + ky->size[0] * i13] = b_ky->data[i14 + b_ky->size[0] *
            i13];
        }
      }

      /*         clear yi ky1; */
      i13 = varargin_2->size[0] * varargin_2->size[1];
      varargin_2->size[0] = U->size[0];
      varargin_2->size[1] = U->size[1] + vstride;
      emxEnsureCapacity((emxArray__common *)varargin_2, i13, sizeof(double));
      loop_ub = U->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = U->size[0];
        for (i14 = 0; i14 < br; i14++) {
          varargin_2->data[i14 + varargin_2->size[0] * i13] = U->data[i14 +
            U->size[0] * i13];
        }
      }

      for (i13 = 0; i13 < vstride; i13++) {
        for (i14 = 0; i14 < trueCount; i14++) {
          varargin_2->data[i14 + varargin_2->size[0] * (i13 + U->size[1])] =
            ky1->data[i14 + trueCount * i13];
        }
      }

      if (!((KY->size[0] == 0) || (KY->size[1] == 0))) {
        trueCount = KY->size[0];
      } else if (!((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0))) {
        trueCount = ky->size[0];
      } else {
        trueCount = KY->size[0];
        if (!(trueCount > 0)) {
          trueCount = 0;
        }

        if (ky->size[0] > trueCount) {
          trueCount = ky->size[0];
        }
      }

      empty_non_axis_sizes = (trueCount == 0);
      if (empty_non_axis_sizes || (!((KY->size[0] == 0) || (KY->size[1] == 0))))
      {
        b_idx = KY->size[1];
      } else {
        b_idx = 0;
      }

      if (empty_non_axis_sizes || (!((varargin_2->size[0] == 0) ||
            (varargin_2->size[1] == 0)))) {
        vstride = varargin_2->size[1];
      } else {
        vstride = 0;
      }

      i13 = b_reshapes[1].f1->size[0] * b_reshapes[1].f1->size[1];
      b_reshapes[1].f1->size[0] = trueCount;
      b_reshapes[1].f1->size[1] = vstride;
      emxEnsureCapacity((emxArray__common *)b_reshapes[1].f1, i13, sizeof(double));
      loop_ub = trueCount * vstride;
      for (i13 = 0; i13 < loop_ub; i13++) {
        b_reshapes[1].f1->data[i13] = varargin_2->data[i13];
      }

      i13 = b_KY->size[0] * b_KY->size[1];
      b_KY->size[0] = trueCount;
      b_KY->size[1] = b_idx + b_reshapes[1].f1->size[1];
      emxEnsureCapacity((emxArray__common *)b_KY, i13, sizeof(double));
      for (i13 = 0; i13 < b_idx; i13++) {
        for (i14 = 0; i14 < trueCount; i14++) {
          b_KY->data[i14 + b_KY->size[0] * i13] = KY->data[i14 + trueCount * i13];
        }
      }

      loop_ub = b_reshapes[1].f1->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = b_reshapes[1].f1->size[0];
        for (i14 = 0; i14 < br; i14++) {
          b_KY->data[i14 + b_KY->size[0] * (i13 + b_idx)] = b_reshapes[1]
            .f1->data[i14 + b_reshapes[1].f1->size[0] * i13];
        }
      }

      i13 = KY->size[0] * KY->size[1];
      KY->size[0] = b_KY->size[0];
      KY->size[1] = b_KY->size[1];
      emxEnsureCapacity((emxArray__common *)KY, i13, sizeof(double));
      loop_ub = b_KY->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = b_KY->size[0];
        for (i14 = 0; i14 < br; i14++) {
          KY->data[i14 + KY->size[0] * i13] = b_KY->data[i14 + b_KY->size[0] *
            i13];
        }
      }
    }

    i13 = a->size[0] * a->size[1];
    a->size[0] = KY->size[1];
    a->size[1] = KY->size[0];
    emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
    loop_ub = KY->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = KY->size[1];
      for (i14 = 0; i14 < br; i14++) {
        a->data[i14 + a->size[0] * i13] = KY->data[i13 + KY->size[0] * i14];
      }
    }

    if ((a->size[1] == 1) || (KY->size[0] == 1)) {
      i13 = C->size[0] * C->size[1];
      C->size[0] = a->size[0];
      C->size[1] = KY->size[1];
      emxEnsureCapacity((emxArray__common *)C, i13, sizeof(double));
      loop_ub = a->size[0];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = KY->size[1];
        for (i14 = 0; i14 < br; i14++) {
          C->data[i13 + C->size[0] * i14] = 0.0;
          vstride = a->size[1];
          for (b_m = 0; b_m < vstride; b_m++) {
            C->data[i13 + C->size[0] * i14] += a->data[i13 + a->size[0] * b_m] *
              KY->data[b_m + KY->size[0] * i14];
          }
        }
      }
    } else {
      k = a->size[1];
      vstride = a->size[0];
      b_idx = KY->size[1];
      i13 = C->size[0] * C->size[1];
      C->size[0] = vstride;
      C->size[1] = b_idx;
      emxEnsureCapacity((emxArray__common *)C, i13, sizeof(double));
      b_m = a->size[0];
      i13 = C->size[0] * C->size[1];
      emxEnsureCapacity((emxArray__common *)C, i13, sizeof(double));
      loop_ub = C->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = C->size[0];
        for (i14 = 0; i14 < br; i14++) {
          C->data[i14 + C->size[0] * i13] = 0.0;
        }
      }

      if ((a->size[0] == 0) || (KY->size[1] == 0)) {
      } else {
        b_idx = a->size[0] * (KY->size[1] - 1);
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          i13 = vstride + b_m;
          for (ic = vstride; ic + 1 <= i13; ic++) {
            C->data[ic] = 0.0;
          }

          vstride += b_m;
        }

        br = 0;
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          ar = 0;
          i13 = br + k;
          for (ib = br; ib + 1 <= i13; ib++) {
            if (KY->data[ib] != 0.0) {
              ia = ar;
              i14 = vstride + b_m;
              for (ic = vstride; ic + 1 <= i14; ic++) {
                ia++;
                C->data[ic] += KY->data[ib] * a->data[ia - 1];
              }
            }

            ar += b_m;
          }

          br += k;
          vstride += b_m;
        }
      }
    }

    c_abs(C, a);
    i13 = r2->size[0] * r2->size[1];
    r2->size[0] = a->size[0];
    r2->size[1] = a->size[1];
    emxEnsureCapacity((emxArray__common *)r2, i13, sizeof(boolean_T));
    loop_ub = a->size[0] * a->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      r2->data[i13] = (a->data[i13] < 1.0E-16);
    }

    vstride = r2->size[0] * r2->size[1] - 1;
    trueCount = 0;
    for (ar = 0; ar <= vstride; ar++) {
      if (r2->data[ar]) {
        trueCount++;
      }
    }

    i13 = r1->size[0];
    r1->size[0] = trueCount;
    emxEnsureCapacity((emxArray__common *)r1, i13, sizeof(int));
    b_idx = 0;
    for (ar = 0; ar <= vstride; ar++) {
      if (r2->data[ar]) {
        r1->data[b_idx] = ar + 1;
        b_idx++;
      }
    }

    loop_ub = r1->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      C->data[r1->data[i13] - 1] = 1.0E-16;
    }

    /* special change: after conversion 0 is used as INF */
    eig(C, Vc, Dc);

    /* change for C */
    i13 = V->size[0] * V->size[1];
    V->size[0] = Vc->size[0];
    V->size[1] = Vc->size[1];
    emxEnsureCapacity((emxArray__common *)V, i13, sizeof(double));
    loop_ub = Vc->size[0] * Vc->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      V->data[i13] = Vc->data[i13].re;
    }

    /* change for C */
    /* change for C */
    /* clear C; */
    i13 = c_Dc->size[0] * c_Dc->size[1];
    c_Dc->size[0] = Dc->size[0];
    c_Dc->size[1] = Dc->size[1];
    emxEnsureCapacity((emxArray__common *)c_Dc, i13, sizeof(double));
    loop_ub = Dc->size[0] * Dc->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      c_Dc->data[i13] = Dc->data[i13].re;
    }

    diag(c_Dc, D);
    d_abs(D, dxij);
    sort(dxij, iidx);
    i13 = D->size[0];
    D->size[0] = dxij->size[0];
    emxEnsureCapacity((emxArray__common *)D, i13, sizeof(double));
    loop_ub = dxij->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      D->data[i13] = dxij->data[i13];
    }

    b_idx = 2;
    if (dxij->size[0] != 1) {
      b_idx = 1;
    }

    if (b_idx <= 1) {
      i13 = dxij->size[0];
    } else {
      i13 = 1;
    }

    if ((!(dxij->size[0] == 0)) && (i13 > 1)) {
      vstride = 1;
      k = 1;
      while (k <= b_idx - 1) {
        vstride *= dxij->size[0];
        k = 2;
      }

      for (j = 0; j + 1 <= vstride; j++) {
        for (k = 1; k < i13; k++) {
          D->data[j + k * vstride] += D->data[j + (k - 1) * vstride];
        }
      }
    }

    pp = d_sum(dxij);
    i13 = c_x->size[0];
    c_x->size[0] = D->size[0];
    emxEnsureCapacity((emxArray__common *)c_x, i13, sizeof(boolean_T));
    loop_ub = D->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      c_x->data[i13] = (D->data[i13] / pp >= 0.99);
    }

    k = c_x->size[0];
    if (1 < k) {
      k = 1;
    }

    b_idx = 0;
    vstride = 1;
    exitg2 = false;
    while ((!exitg2) && (vstride <= c_x->size[0])) {
      if (c_x->data[vstride - 1]) {
        b_idx = 1;
        ii_data[0] = vstride;
        exitg2 = true;
      } else {
        vstride++;
      }
    }

    if (k == 1) {
      if (b_idx == 0) {
        k = 0;
      }
    } else {
      k = !(1 > b_idx);
    }

    for (i13 = 0; i13 < k; i13++) {
      di_data[i13] = ii_data[i13];
    }

    b_idx = V->size[0];
    i13 = b_V->size[0] * b_V->size[1];
    b_V->size[0] = b_idx;
    b_V->size[1] = iidx->size[0];
    emxEnsureCapacity((emxArray__common *)b_V, i13, sizeof(double));
    loop_ub = iidx->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      for (i14 = 0; i14 < b_idx; i14++) {
        b_V->data[i14 + b_V->size[0] * i13] = V->data[i14 + V->size[0] *
          (iidx->data[i13] - 1)];
      }
    }

    i13 = V->size[0] * V->size[1];
    V->size[0] = b_V->size[0];
    V->size[1] = b_V->size[1];
    emxEnsureCapacity((emxArray__common *)V, i13, sizeof(double));
    loop_ub = b_V->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = b_V->size[0];
      for (i14 = 0; i14 < br; i14++) {
        V->data[i14 + V->size[0] * i13] = b_V->data[i14 + b_V->size[0] * i13];
      }
    }

    if (1 > di_data[0]) {
      loop_ub = 0;
    } else {
      loop_ub = di_data[0];
    }

    if (1 > di_data[0]) {
      br = 0;
    } else {
      br = di_data[0];
    }

    b_idx = V->size[0];
    i13 = c_V->size[0] * c_V->size[1];
    c_V->size[0] = b_idx;
    c_V->size[1] = br;
    emxEnsureCapacity((emxArray__common *)c_V, i13, sizeof(double));
    for (i13 = 0; i13 < br; i13++) {
      for (i14 = 0; i14 < b_idx; i14++) {
        c_V->data[i14 + c_V->size[0] * i13] = V->data[i14 + V->size[0] * i13];
      }
    }

    i13 = V->size[0] * V->size[1];
    V->size[0] = c_V->size[0];
    V->size[1] = c_V->size[1];
    emxEnsureCapacity((emxArray__common *)V, i13, sizeof(double));
    br = c_V->size[1];
    for (i13 = 0; i13 < br; i13++) {
      vstride = c_V->size[0];
      for (i14 = 0; i14 < vstride; i14++) {
        V->data[i14 + V->size[0] * i13] = c_V->data[i14 + c_V->size[0] * i13];
      }
    }

    if ((KY->size[1] == 1) || (V->size[0] == 1)) {
      i13 = U->size[0] * U->size[1];
      U->size[0] = KY->size[0];
      U->size[1] = V->size[1];
      emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
      br = KY->size[0];
      for (i13 = 0; i13 < br; i13++) {
        vstride = V->size[1];
        for (i14 = 0; i14 < vstride; i14++) {
          U->data[i13 + U->size[0] * i14] = 0.0;
          b_idx = KY->size[1];
          for (b_m = 0; b_m < b_idx; b_m++) {
            U->data[i13 + U->size[0] * i14] += KY->data[i13 + KY->size[0] * b_m]
              * V->data[b_m + V->size[0] * i14];
          }
        }
      }
    } else {
      k = KY->size[1];
      vstride = KY->size[0];
      b_idx = V->size[1];
      i13 = U->size[0] * U->size[1];
      U->size[0] = vstride;
      U->size[1] = b_idx;
      emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
      b_m = KY->size[0];
      i13 = U->size[0] * U->size[1];
      emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
      br = U->size[1];
      for (i13 = 0; i13 < br; i13++) {
        vstride = U->size[0];
        for (i14 = 0; i14 < vstride; i14++) {
          U->data[i14 + U->size[0] * i13] = 0.0;
        }
      }

      if ((KY->size[0] == 0) || (V->size[1] == 0)) {
      } else {
        b_idx = KY->size[0] * (V->size[1] - 1);
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          i13 = vstride + b_m;
          for (ic = vstride; ic + 1 <= i13; ic++) {
            U->data[ic] = 0.0;
          }

          vstride += b_m;
        }

        br = 0;
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          ar = 0;
          i13 = br + k;
          for (ib = br; ib + 1 <= i13; ib++) {
            if (V->data[ib] != 0.0) {
              ia = ar;
              i14 = vstride + b_m;
              for (ic = vstride; ic + 1 <= i14; ic++) {
                ia++;
                U->data[ic] += V->data[ib] * KY->data[ia - 1];
              }
            }

            ar += b_m;
          }

          br += k;
          vstride += b_m;
        }
      }
    }

    i13 = D->size[0];
    D->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)D, i13, sizeof(double));
    for (i13 = 0; i13 < loop_ub; i13++) {
      D->data[i13] = dxij->data[i13];
    }

    d_sqrt(D);

    /* ky1 = bsxfun(@(x,c)x./c, U, s'); */
    i13 = c_D->size[0] * c_D->size[1];
    c_D->size[0] = 1;
    c_D->size[1] = D->size[0];
    emxEnsureCapacity((emxArray__common *)c_D, i13, sizeof(double));
    loop_ub = D->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      c_D->data[c_D->size[0] * i13] = D->data[i13];
    }

    repmat(c_D, (double)U->size[0], a);
    b_rdivide(U, a, ky1);
    b_diag(D, a);
    i13 = b->size[0] * b->size[1];
    b->size[0] = V->size[1];
    b->size[1] = V->size[0];
    emxEnsureCapacity((emxArray__common *)b, i13, sizeof(double));
    loop_ub = V->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = V->size[1];
      for (i14 = 0; i14 < br; i14++) {
        b->data[i14 + b->size[0] * i13] = V->data[i13 + V->size[0] * i14];
      }
    }

    if ((a->size[1] == 1) || (b->size[0] == 1)) {
      i13 = ky2->size[0] * ky2->size[1];
      ky2->size[0] = a->size[0];
      ky2->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)ky2, i13, sizeof(double));
      loop_ub = a->size[0];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = b->size[1];
        for (i14 = 0; i14 < br; i14++) {
          ky2->data[i13 + ky2->size[0] * i14] = 0.0;
          vstride = a->size[1];
          for (b_m = 0; b_m < vstride; b_m++) {
            ky2->data[i13 + ky2->size[0] * i14] += a->data[i13 + a->size[0] *
              b_m] * b->data[b_m + b->size[0] * i14];
          }
        }
      }
    } else {
      k = a->size[1];
      vstride = a->size[0];
      b_idx = b->size[1];
      i13 = ky2->size[0] * ky2->size[1];
      ky2->size[0] = vstride;
      ky2->size[1] = b_idx;
      emxEnsureCapacity((emxArray__common *)ky2, i13, sizeof(double));
      b_m = a->size[0];
      i13 = ky2->size[0] * ky2->size[1];
      emxEnsureCapacity((emxArray__common *)ky2, i13, sizeof(double));
      loop_ub = ky2->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = ky2->size[0];
        for (i14 = 0; i14 < br; i14++) {
          ky2->data[i14 + ky2->size[0] * i13] = 0.0;
        }
      }

      if ((a->size[0] == 0) || (b->size[1] == 0)) {
      } else {
        b_idx = a->size[0] * (b->size[1] - 1);
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          i13 = vstride + b_m;
          for (ic = vstride; ic + 1 <= i13; ic++) {
            ky2->data[ic] = 0.0;
          }

          vstride += b_m;
        }

        br = 0;
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          ar = 0;
          i13 = br + k;
          for (ib = br; ib + 1 <= i13; ib++) {
            if (b->data[ib] != 0.0) {
              ia = ar;
              i14 = vstride + b_m;
              for (ic = vstride; ic + 1 <= i14; ic++) {
                ia++;
                ky2->data[ic] += b->data[ib] * a->data[ia - 1];
              }
            }

            ar += b_m;
          }

          br += k;
          vstride += b_m;
        }
      }
    }

    if ((ky1->size[1] == 1) || (ky2->size[0] == 1)) {
      i13 = ky->size[0] * ky->size[1];
      ky->size[0] = ky1->size[0];
      ky->size[1] = ky2->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i13, sizeof(double));
      loop_ub = ky1->size[0];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = ky2->size[1];
        for (i14 = 0; i14 < br; i14++) {
          ky->data[i13 + ky->size[0] * i14] = 0.0;
          vstride = ky1->size[1];
          for (b_m = 0; b_m < vstride; b_m++) {
            ky->data[i13 + ky->size[0] * i14] += ky1->data[i13 + ky1->size[0] *
              b_m] * ky2->data[b_m + ky2->size[0] * i14];
          }
        }
      }
    } else {
      k = ky1->size[1];
      vstride = ky1->size[0];
      b_idx = ky2->size[1];
      i13 = ky->size[0] * ky->size[1];
      ky->size[0] = vstride;
      ky->size[1] = b_idx;
      emxEnsureCapacity((emxArray__common *)ky, i13, sizeof(double));
      b_m = ky1->size[0];
      i13 = ky->size[0] * ky->size[1];
      emxEnsureCapacity((emxArray__common *)ky, i13, sizeof(double));
      loop_ub = ky->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = ky->size[0];
        for (i14 = 0; i14 < br; i14++) {
          ky->data[i14 + ky->size[0] * i13] = 0.0;
        }
      }

      if ((ky1->size[0] == 0) || (ky2->size[1] == 0)) {
      } else {
        b_idx = ky1->size[0] * (ky2->size[1] - 1);
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          i13 = vstride + b_m;
          for (ic = vstride; ic + 1 <= i13; ic++) {
            ky->data[ic] = 0.0;
          }

          vstride += b_m;
        }

        br = 0;
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          ar = 0;
          i13 = br + k;
          for (ib = br; ib + 1 <= i13; ib++) {
            if (ky2->data[ib] != 0.0) {
              ia = ar;
              i14 = vstride + b_m;
              for (ic = vstride; ic + 1 <= i14; ic++) {
                ia++;
                ky->data[ic] += ky2->data[ib] * ky1->data[ia - 1];
              }
            }

            ar += b_m;
          }

          br += k;
          vstride += b_m;
        }
      }
    }

    /* clear U V; */
    i13 = b->size[0] * b->size[1];
    b->size[0] = ky2->size[1];
    b->size[1] = ky2->size[0];
    emxEnsureCapacity((emxArray__common *)b, i13, sizeof(double));
    loop_ub = ky2->size[0];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = ky2->size[1];
      for (i14 = 0; i14 < br; i14++) {
        b->data[i14 + b->size[0] * i13] = ky2->data[i13 + ky2->size[0] * i14];
      }
    }

    if ((ky2->size[1] == 1) || (b->size[0] == 1)) {
      i13 = DD->size[0] * DD->size[1];
      DD->size[0] = ky2->size[0];
      DD->size[1] = b->size[1];
      emxEnsureCapacity((emxArray__common *)DD, i13, sizeof(double));
      loop_ub = ky2->size[0];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = b->size[1];
        for (i14 = 0; i14 < br; i14++) {
          DD->data[i13 + DD->size[0] * i14] = 0.0;
          vstride = ky2->size[1];
          for (b_m = 0; b_m < vstride; b_m++) {
            DD->data[i13 + DD->size[0] * i14] += ky2->data[i13 + ky2->size[0] *
              b_m] * b->data[b_m + b->size[0] * i14];
          }
        }
      }
    } else {
      k = ky2->size[1];
      vstride = ky2->size[0];
      b_idx = b->size[1];
      i13 = DD->size[0] * DD->size[1];
      DD->size[0] = vstride;
      DD->size[1] = b_idx;
      emxEnsureCapacity((emxArray__common *)DD, i13, sizeof(double));
      b_m = ky2->size[0];
      i13 = DD->size[0] * DD->size[1];
      emxEnsureCapacity((emxArray__common *)DD, i13, sizeof(double));
      loop_ub = DD->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = DD->size[0];
        for (i14 = 0; i14 < br; i14++) {
          DD->data[i14 + DD->size[0] * i13] = 0.0;
        }
      }

      if ((ky2->size[0] == 0) || (b->size[1] == 0)) {
      } else {
        b_idx = ky2->size[0] * (b->size[1] - 1);
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          i13 = vstride + b_m;
          for (ic = vstride; ic + 1 <= i13; ic++) {
            DD->data[ic] = 0.0;
          }

          vstride += b_m;
        }

        br = 0;
        vstride = 0;
        while ((b_m > 0) && (vstride <= b_idx)) {
          ar = 0;
          i13 = br + k;
          for (ib = br; ib + 1 <= i13; ib++) {
            if (b->data[ib] != 0.0) {
              ia = ar;
              i14 = vstride + b_m;
              for (ic = vstride; ic + 1 <= i14; ic++) {
                ia++;
                DD->data[ic] += b->data[ib] * ky2->data[ia - 1];
              }
            }

            ar += b_m;
          }

          br += k;
          vstride += b_m;
        }
      }
    }
  }

  emxFree_real_T(&b_KY);
  emxFree_real_T(&b_ky);
  emxFree_real_T(&c_V);
  emxFree_real_T(&b_V);
  emxFree_real_T(&b_meanky);
  emxFree_boolean_T(&b_yi);
  emxFree_real_T(&c_Dc);
  emxFree_real_T(&c_D);
  emxFreeMatrix_cell_wrap_0(b_reshapes);
  emxFree_real_T(&varargin_2);
  emxFreeMatrix_cell_wrap_0(reshapes);
  emxFree_int32_T(&c_idx);
  emxFree_real_T(&b_a);
  emxFree_boolean_T(&c_x);
  emxFree_boolean_T(&r2);
  emxFree_int32_T(&r1);
  emxFree_real_T(&yj);
  emxFree_real_T(&yi);
  emxFree_real_T(&qx);
  if (p == 0.0) {
    i13 = BBvs->size[0] * BBvs->size[1] * BBvs->size[2];
    BBvs->size[0] = 1;
    BBvs->size[1] = 1;
    BBvs->size[2] = 1;
    emxEnsureCapacity((emxArray__common *)BBvs, i13, sizeof(double));
    BBvs->data[0] = 1.0;
    i13 = BB->size[0] * BB->size[1] * BB->size[2];
    BB->size[0] = 1;
    BB->size[1] = 1;
    BB->size[2] = 1;
    emxEnsureCapacity((emxArray__common *)BB, i13, sizeof(double));
    BB->data[0] = 1.0;
  } else {
    i13 = BB->size[0] * BB->size[1] * BB->size[2];
    BB->size[0] = (int)p;
    BB->size[1] = (int)p;
    BB->size[2] = (int)p;
    emxEnsureCapacity((emxArray__common *)BB, i13, sizeof(double));
    loop_ub = (int)p * (int)p * (int)p;
    for (i13 = 0; i13 < loop_ub; i13++) {
      BB->data[i13] = 0.0;
    }

    emxInit_real_T1(&r3, 2);
    eye(p, r3);
    loop_ub = r3->size[1];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = r3->size[0];
      for (i14 = 0; i14 < br; i14++) {
        BB->data[(i14 + BB->size[0] * i13) + BB->size[0] * BB->size[1] * ((int)p
          - 1)] = r3->data[i14 + r3->size[0] * i13];
      }
    }

    emxFree_real_T(&r3);
    emxInit_real_T1(&onexi, 2);
    i13 = onexi->size[0] * onexi->size[1];
    onexi->size[0] = n;
    onexi->size[1] = (int)(p + 1.0);
    emxEnsureCapacity((emxArray__common *)onexi, i13, sizeof(double));
    loop_ub = n * (int)(p + 1.0);
    for (i13 = 0; i13 < loop_ub; i13++) {
      onexi->data[i13] = 1.0;
    }

    emxInit_real_T1(&B, 2);
    eye(p, B);
    i13 = U->size[0] * U->size[1];
    U->size[0] = (int)p;
    U->size[1] = (int)p;
    emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
    loop_ub = (int)p * (int)p;
    for (i13 = 0; i13 < loop_ub; i13++) {
      U->data[i13] = 0.0;
    }

    if (p > 1.0) {
      if ((p < max_dim) || rtIsNaN(max_dim)) {
        max_dim = p;
      }

      emxInit_real_T1(&SEQ, 2);
      if (p == p) {
        i13 = SEQ->size[0] * SEQ->size[1];
        SEQ->size[0] = 1;
        SEQ->size[1] = (int)std::floor(-(1.0 - p)) + 1;
        emxEnsureCapacity((emxArray__common *)SEQ, i13, sizeof(double));
        loop_ub = (int)std::floor(-(1.0 - p));
        for (i13 = 0; i13 <= loop_ub; i13++) {
          SEQ->data[SEQ->size[0] * i13] = p - (double)i13;
        }
      } else {
        b_idx = (int)std::floor((1.0 - p) / -1.0 + 0.5) - 1;
        vstride = ((int)p - b_idx) - 1;
        if (std::abs(1.0 - (double)vstride) < 4.4408920985006262E-16 * (double)
            (int)p) {
          b_idx++;
          vstride = 1;
        } else if (1 - vstride > 0) {
          vstride = (int)p - b_idx;
        } else {
          b_idx++;
        }

        i13 = SEQ->size[0] * SEQ->size[1];
        SEQ->size[0] = 1;
        SEQ->size[1] = b_idx + 1;
        emxEnsureCapacity((emxArray__common *)SEQ, i13, sizeof(double));
        SEQ->data[0] = p;
        if (b_idx + 1 > 1) {
          SEQ->data[b_idx] = vstride;
          br = b_idx / 2;
          for (k = 1; k < br; k++) {
            SEQ->data[k] = p + -(double)k;
            SEQ->data[b_idx - k] = (double)vstride - (-(double)k);
          }

          if (br << 1 == b_idx) {
            SEQ->data[br] = (p + (double)vstride) / 2.0;
          } else {
            SEQ->data[br] = p + -(double)br;
            SEQ->data[br + 1] = (double)vstride - (-(double)br);
          }
        }
      }

      if (p > 20.0) {
        trueCount = 0;
        for (ar = 0; ar < 51; ar++) {
          pp = c_a[ar] * p;
          if (pp >= max_dim) {
            trueCount++;
          }

          Ip[ar] = pp;
        }

        b_idx = 0;
        for (ar = 0; ar < 51; ar++) {
          if (Ip[ar] >= max_dim) {
            tmp_data[b_idx] = (signed char)(ar + 1);
            b_idx++;
          }
        }

        if (rtIsNaN(max_dim)) {
          i13 = meanky->size[0] * meanky->size[1];
          meanky->size[0] = 1;
          meanky->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
          meanky->data[0] = rtNaN;
        } else if (max_dim < 1.0) {
          i13 = meanky->size[0] * meanky->size[1];
          meanky->size[0] = 1;
          meanky->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
        } else if (rtIsInf(max_dim) && (max_dim == 1.0)) {
          i13 = meanky->size[0] * meanky->size[1];
          meanky->size[0] = 1;
          meanky->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
          meanky->data[0] = rtNaN;
        } else if (std::floor(max_dim) == max_dim) {
          i13 = meanky->size[0] * meanky->size[1];
          meanky->size[0] = 1;
          meanky->size[1] = (int)std::floor(-(1.0 - max_dim)) + 1;
          emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
          loop_ub = (int)std::floor(-(1.0 - max_dim));
          for (i13 = 0; i13 <= loop_ub; i13++) {
            meanky->data[meanky->size[0] * i13] = max_dim - (double)i13;
          }
        } else {
          ndbl = std::floor((1.0 - max_dim) / -1.0 + 0.5);
          apnd = max_dim + -ndbl;
          if (max_dim > 1.0) {
            pp = max_dim;
          } else {
            pp = 1.0;
          }

          if (std::abs(1.0 - apnd) < 4.4408920985006262E-16 * pp) {
            ndbl++;
            apnd = 1.0;
          } else if (1.0 - apnd > 0.0) {
            apnd = max_dim + -(ndbl - 1.0);
          } else {
            ndbl++;
          }

          b_idx = (int)ndbl;
          i13 = meanky->size[0] * meanky->size[1];
          meanky->size[0] = 1;
          meanky->size[1] = b_idx;
          emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof(double));
          if (b_idx > 0) {
            meanky->data[0] = max_dim;
            if (b_idx > 1) {
              meanky->data[b_idx - 1] = apnd;
              br = (b_idx - 1) / 2;
              for (k = 1; k < br; k++) {
                meanky->data[k] = max_dim + -(double)k;
                meanky->data[(b_idx - k) - 1] = apnd - (-(double)k);
              }

              if (br << 1 == b_idx - 1) {
                meanky->data[br] = (max_dim + apnd) / 2.0;
              } else {
                meanky->data[br] = max_dim + -(double)br;
                meanky->data[br + 1] = apnd - (-(double)br);
              }
            }
          }
        }

        tmp_size[0] = 1;
        tmp_size[1] = trueCount;
        for (i13 = 0; i13 < trueCount; i13++) {
          b_tmp_data[i13] = Ip[tmp_data[i13] - 1];
        }

        b_floor(b_tmp_data, tmp_size);
        i13 = SEQ->size[0] * SEQ->size[1];
        SEQ->size[0] = 1;
        SEQ->size[1] = tmp_size[1] + meanky->size[1];
        emxEnsureCapacity((emxArray__common *)SEQ, i13, sizeof(double));
        loop_ub = tmp_size[1];
        for (i13 = 0; i13 < loop_ub; i13++) {
          SEQ->data[SEQ->size[0] * i13] = b_tmp_data[tmp_size[0] * i13];
        }

        loop_ub = meanky->size[1];
        for (i13 = 0; i13 < loop_ub; i13++) {
          SEQ->data[SEQ->size[0] * (i13 + tmp_size[1])] = meanky->data
            [meanky->size[0] * i13];
        }
      }

      iter_ip = 0;
      emxInit_real_T1(&Ifast, 2);
      emxInit_real_T1(&Vfast, 2);
      emxInit_real_T1(&xfast, 2);
      emxInit_real_T1(&xij, 2);
      emxInit_real_T1(&xk, 2);
      emxInit_real_T1(&abi, 2);
      emxInit_real_T1(&dd, 2);
      emxInit_real_T(&dc, 1);
      emxInit_real_T1(&d_D, 2);
      emxInit_creal_T1(&R, 1);
      emxInit_int32_T1(&r4, 2);
      emxInitMatrix_cell_wrap_0(c_reshapes);
      emxInit_real_T1(&b_dd, 2);
      emxInit_real_T1(&b_abi, 2);
      emxInit_real_T1(&d_a, 2);
      emxInit_real_T1(&d_V, 2);
      emxInit_real_T1(&b_Vfast, 2);
      emxInit_real_T1(&b_xfast, 2);
      emxInit_real_T1(&e_a, 2);
      emxInit_real_T(&e_V, 1);
      emxInit_real_T1(&b_B, 2);
      emxInit_int32_T1(&b_Ifast, 2);
      emxInit_int32_T1(&c_Ifast, 2);
      emxInit_real_T(&b_dc, 1);
      emxInit_real_T1(&c_abi, 2);
      emxInit_int32_T1(&d_Ifast, 2);
      emxInit_int32_T1(&e_Ifast, 2);
      exitg2 = false;
      while ((!exitg2) && (iter_ip <= SEQ->size[1] - 1)) {
        ip = SEQ->data[iter_ip];
        b_idx = (int)std::floor((double)n / 2.0);
        pp = SEQ->data[iter_ip];
        m = b_idx;
        if (pp < m) {
          m = pp;
        }

        K = ((SEQ->data[iter_ip] <= 10.0) << 2) + (SEQ->data[iter_ip] < 20.0) *
          (SEQ->data[iter_ip] > 10.0);
        if (d_strcmp(method)) {
          K = 0;
        }

        if (1.0 > SEQ->data[iter_ip]) {
          loop_ub = 0;
        } else {
          loop_ub = (int)SEQ->data[iter_ip];
        }

        b_idx = B->size[0];
        i13 = b_B->size[0] * b_B->size[1];
        b_B->size[0] = b_idx;
        b_B->size[1] = loop_ub;
        emxEnsureCapacity((emxArray__common *)b_B, i13, sizeof(double));
        for (i13 = 0; i13 < loop_ub; i13++) {
          for (i14 = 0; i14 < b_idx; i14++) {
            b_B->data[i14 + b_B->size[0] * i13] = B->data[i14 + B->size[0] * i13];
          }
        }

        i13 = B->size[0] * B->size[1];
        B->size[0] = b_B->size[0];
        B->size[1] = b_B->size[1];
        emxEnsureCapacity((emxArray__common *)B, i13, sizeof(double));
        loop_ub = b_B->size[1];
        for (i13 = 0; i13 < loop_ub; i13++) {
          br = b_B->size[0];
          for (i14 = 0; i14 < br; i14++) {
            B->data[i14 + B->size[0] * i13] = b_B->data[i14 + b_B->size[0] * i13];
          }
        }

        for (iter = 0; iter <= K; iter++) {
          if (d_strcmp(method)) {
            i13 = V->size[0] * V->size[1];
            V->size[0] = y->size[0];
            V->size[1] = y->size[1];
            emxEnsureCapacity((emxArray__common *)V, i13, sizeof(double));
            loop_ub = y->size[0] * y->size[1];
            for (i13 = 0; i13 < loop_ub; i13++) {
              V->data[i13] = y->data[i13];
            }

            m = 1.0;
          } else if ((x->size[1] == 1) || (B->size[0] == 1)) {
            i13 = V->size[0] * V->size[1];
            V->size[0] = x->size[0];
            V->size[1] = B->size[1];
            emxEnsureCapacity((emxArray__common *)V, i13, sizeof(double));
            loop_ub = x->size[0];
            for (i13 = 0; i13 < loop_ub; i13++) {
              br = B->size[1];
              for (i14 = 0; i14 < br; i14++) {
                V->data[i13 + V->size[0] * i14] = 0.0;
                vstride = x->size[1];
                for (b_m = 0; b_m < vstride; b_m++) {
                  V->data[i13 + V->size[0] * i14] += x->data[i13 + x->size[0] *
                    b_m] * B->data[b_m + B->size[0] * i14];
                }
              }
            }
          } else {
            k = x->size[1];
            vstride = x->size[0];
            b_idx = B->size[1];
            i13 = V->size[0] * V->size[1];
            V->size[0] = vstride;
            V->size[1] = b_idx;
            emxEnsureCapacity((emxArray__common *)V, i13, sizeof(double));
            b_m = x->size[0];
            i13 = V->size[0] * V->size[1];
            emxEnsureCapacity((emxArray__common *)V, i13, sizeof(double));
            loop_ub = V->size[1];
            for (i13 = 0; i13 < loop_ub; i13++) {
              br = V->size[0];
              for (i14 = 0; i14 < br; i14++) {
                V->data[i14 + V->size[0] * i13] = 0.0;
              }
            }

            if ((x->size[0] == 0) || (B->size[1] == 0)) {
            } else {
              b_idx = x->size[0] * (B->size[1] - 1);
              vstride = 0;
              while ((b_m > 0) && (vstride <= b_idx)) {
                i13 = vstride + b_m;
                for (ic = vstride; ic + 1 <= i13; ic++) {
                  V->data[ic] = 0.0;
                }

                vstride += b_m;
              }

              br = 0;
              vstride = 0;
              while ((b_m > 0) && (vstride <= b_idx)) {
                ar = 0;
                i13 = br + k;
                for (ib = br; ib + 1 <= i13; ib++) {
                  if (B->data[ib] != 0.0) {
                    ia = ar;
                    i14 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i14; ic++) {
                      ia++;
                      V->data[ic] += B->data[ib] * x->data[ia - 1];
                    }
                  }

                  ar += b_m;
                }

                br += k;
                vstride += b_m;
              }
            }
          }

          if (100 < n) {
            b_idx = 100;
          } else {
            b_idx = n;
          }

          e_x = std::floor(rt_powd_snf((double)n, 0.5));
          if ((b_idx > e_x) || rtIsNaN(e_x)) {
            d_idx = b_idx;
          } else {
            d_idx = e_x;
          }

          unifD(V, d_idx, Ifast);
          i13 = b_Ifast->size[0] * b_Ifast->size[1];
          b_Ifast->size[0] = Ifast->size[0];
          b_Ifast->size[1] = Ifast->size[1];
          emxEnsureCapacity((emxArray__common *)b_Ifast, i13, sizeof(int));
          loop_ub = Ifast->size[1];
          for (i13 = 0; i13 < loop_ub; i13++) {
            br = Ifast->size[0];
            for (i14 = 0; i14 < br; i14++) {
              b_Ifast->data[i14 + b_Ifast->size[0] * i13] = (int)Ifast->data[i14
                + Ifast->size[0] * i13];
            }
          }

          exponent = Ifast->size[0] * Ifast->size[1];
          loop_ub = V->size[1];
          i13 = Vfast->size[0] * Vfast->size[1];
          Vfast->size[0] = exponent;
          Vfast->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)Vfast, i13, sizeof(double));
          for (i13 = 0; i13 < loop_ub; i13++) {
            for (i14 = 0; i14 < exponent; i14++) {
              Vfast->data[i14 + Vfast->size[0] * i13] = V->data[(b_Ifast->
                data[i14] + V->size[0] * i13) - 1];
            }
          }

          i13 = c_Ifast->size[0] * c_Ifast->size[1];
          c_Ifast->size[0] = Ifast->size[0];
          c_Ifast->size[1] = Ifast->size[1];
          emxEnsureCapacity((emxArray__common *)c_Ifast, i13, sizeof(int));
          loop_ub = Ifast->size[1];
          for (i13 = 0; i13 < loop_ub; i13++) {
            br = Ifast->size[0];
            for (i14 = 0; i14 < br; i14++) {
              c_Ifast->data[i14 + c_Ifast->size[0] * i13] = (int)Ifast->data[i14
                + Ifast->size[0] * i13];
            }
          }

          exponent = Ifast->size[0] * Ifast->size[1];
          loop_ub = x->size[1];
          i13 = xfast->size[0] * xfast->size[1];
          xfast->size[0] = exponent;
          xfast->size[1] = loop_ub;
          emxEnsureCapacity((emxArray__common *)xfast, i13, sizeof(double));
          for (i13 = 0; i13 < loop_ub; i13++) {
            for (i14 = 0; i14 < exponent; i14++) {
              xfast->data[i14 + xfast->size[0] * i13] = x->data[(c_Ifast->
                data[i14] + x->size[0] * i13) - 1];
            }
          }

          exponent = Ifast->size[0] * Ifast->size[1];
          b_std(V, meanky);
          pp = 1.2 * d_mean(meanky);
          pp /= rt_powd_snf((double)n, 1.0 / (m + 4.0));
          h2 = 2.0 * pp * pp;
          if (e_strcmp(method) || d_strcmp(method) || c_strcmp(method) || (ip >
               5.0)) {
            i13 = U->size[0] * U->size[1];
            U->size[0] = (int)p;
            U->size[1] = (int)p;
            emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
            loop_ub = (int)p * (int)p;
            for (i13 = 0; i13 < loop_ub; i13++) {
              U->data[i13] = 0.0;
            }

            for (j = 0; j < exponent; j++) {
              i13 = xij->size[0] * xij->size[1];
              xij->size[0] = n;
              xij->size[1] = (int)p;
              emxEnsureCapacity((emxArray__common *)xij, i13, sizeof(double));
              for (ar = 0; ar < (int)p; ar++) {
                loop_ub = x->size[0] - 1;
                i13 = e_Ifast->size[0] * e_Ifast->size[1];
                e_Ifast->size[0] = Ifast->size[0];
                e_Ifast->size[1] = Ifast->size[1];
                emxEnsureCapacity((emxArray__common *)e_Ifast, i13, sizeof(int));
                br = Ifast->size[1];
                for (i13 = 0; i13 < br; i13++) {
                  vstride = Ifast->size[0];
                  for (i14 = 0; i14 < vstride; i14++) {
                    e_Ifast->data[i14 + e_Ifast->size[0] * i13] = (int)
                      Ifast->data[i14 + Ifast->size[0] * i13];
                  }
                }

                e_x = x->data[(e_Ifast->data[j] + x->size[0] * ar) - 1];
                for (i13 = 0; i13 <= loop_ub; i13++) {
                  xij->data[i13 + xij->size[0] * ar] = x->data[i13 + x->size[0] *
                    ar] - e_x;
                }
              }

              i13 = dxij->size[0];
              dxij->size[0] = n;
              emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
              for (i13 = 0; i13 < n; i13++) {
                dxij->data[i13] = 0.0;
              }

              for (ar = 0; ar < (int)m; ar++) {
                loop_ub = V->size[0];
                i13 = d_Ifast->size[0] * d_Ifast->size[1];
                d_Ifast->size[0] = Ifast->size[0];
                d_Ifast->size[1] = Ifast->size[1];
                emxEnsureCapacity((emxArray__common *)d_Ifast, i13, sizeof(int));
                br = Ifast->size[1];
                for (i13 = 0; i13 < br; i13++) {
                  vstride = Ifast->size[0];
                  for (i14 = 0; i14 < vstride; i14++) {
                    d_Ifast->data[i14 + d_Ifast->size[0] * i13] = (int)
                      Ifast->data[i14 + Ifast->size[0] * i13];
                  }
                }

                pp = V->data[(d_Ifast->data[j] + V->size[0] * ar) - 1];
                i13 = e_V->size[0];
                e_V->size[0] = loop_ub;
                emxEnsureCapacity((emxArray__common *)e_V, i13, sizeof(double));
                for (i13 = 0; i13 < loop_ub; i13++) {
                  e_V->data[i13] = V->data[i13 + V->size[0] * ar] - pp;
                }

                b_power(e_V, D);
                i13 = dxij->size[0];
                emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
                loop_ub = dxij->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  dxij->data[i13] += D->data[i13];
                }
              }

              i13 = D->size[0];
              D->size[0] = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)D, i13, sizeof(double));
              loop_ub = dxij->size[0];
              for (i13 = 0; i13 < loop_ub; i13++) {
                D->data[i13] = dxij->data[i13];
              }

              d_sort(D);
              if ((h2 > D->data[(int)(2.0 * m) - 1]) || rtIsNaN(D->data[(int)
                   (2.0 * m) - 1])) {
                pp = h2;
              } else {
                pp = D->data[(int)(2.0 * m) - 1];
              }

              i13 = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
              loop_ub = dxij->size[0];
              for (i13 = 0; i13 < loop_ub; i13++) {
                dxij->data[i13] = -dxij->data[i13] / pp;
              }

              b_exp(dxij);
              loop_ub = xij->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = xij->size[0];
                for (i14 = 0; i14 < br; i14++) {
                  onexi->data[i14 + onexi->size[0] * i13] = xij->data[i14 +
                    xij->size[0] * i13];
                }
              }

              b_repmat(dxij, p + 1.0, a);
              i13 = xk->size[0] * xk->size[1];
              xk->size[0] = onexi->size[1];
              xk->size[1] = onexi->size[0];
              emxEnsureCapacity((emxArray__common *)xk, i13, sizeof(double));
              loop_ub = onexi->size[0];
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = onexi->size[1];
                for (i14 = 0; i14 < br; i14++) {
                  xk->data[i14 + xk->size[0] * i13] = onexi->data[i13 +
                    onexi->size[0] * i14] * a->data[i13 + a->size[0] * i14];
                }
              }

              /* abi = inv(xk*onexi+eye(p+1)/n)*(xk*ky1); */
              if ((xk->size[1] == 1) || (onexi->size[0] == 1)) {
                i13 = a->size[0] * a->size[1];
                a->size[0] = xk->size[0];
                a->size[1] = onexi->size[1];
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                loop_ub = xk->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = onexi->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    a->data[i13 + a->size[0] * i14] = 0.0;
                    vstride = xk->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      a->data[i13 + a->size[0] * i14] += xk->data[i13 + xk->
                        size[0] * b_m] * onexi->data[b_m + onexi->size[0] * i14];
                    }
                  }
                }
              } else {
                k = xk->size[1];
                vstride = xk->size[0];
                b_idx = onexi->size[1];
                i13 = a->size[0] * a->size[1];
                a->size[0] = vstride;
                a->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                b_m = xk->size[0];
                i13 = a->size[0] * a->size[1];
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                loop_ub = a->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = a->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    a->data[i14 + a->size[0] * i13] = 0.0;
                  }
                }

                if ((xk->size[0] == 0) || (onexi->size[1] == 0)) {
                } else {
                  b_idx = xk->size[0] * (onexi->size[1] - 1);
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    i13 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i13; ic++) {
                      a->data[ic] = 0.0;
                    }

                    vstride += b_m;
                  }

                  br = 0;
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (onexi->data[ib] != 0.0) {
                        ia = ar;
                        i14 = vstride + b_m;
                        for (ic = vstride; ic + 1 <= i14; ic++) {
                          ia++;
                          a->data[ic] += onexi->data[ib] * xk->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    vstride += b_m;
                  }
                }
              }

              eye(p + 1.0, KY);
              if ((xk->size[1] == 1) || (ky1->size[0] == 1)) {
                i13 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = xk->size[0];
                b_y->size[1] = ky1->size[1];
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                loop_ub = xk->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = ky1->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    b_y->data[i13 + b_y->size[0] * i14] = 0.0;
                    vstride = xk->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      b_y->data[i13 + b_y->size[0] * i14] += xk->data[i13 +
                        xk->size[0] * b_m] * ky1->data[b_m + ky1->size[0] * i14];
                    }
                  }
                }
              } else {
                k = xk->size[1];
                vstride = xk->size[0];
                b_idx = ky1->size[1];
                i13 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = vstride;
                b_y->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                b_m = xk->size[0];
                i13 = b_y->size[0] * b_y->size[1];
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                loop_ub = b_y->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = b_y->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    b_y->data[i14 + b_y->size[0] * i13] = 0.0;
                  }
                }

                if ((xk->size[0] == 0) || (ky1->size[1] == 0)) {
                } else {
                  b_idx = xk->size[0] * (ky1->size[1] - 1);
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    i13 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i13; ic++) {
                      b_y->data[ic] = 0.0;
                    }

                    vstride += b_m;
                  }

                  br = 0;
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (ky1->data[ib] != 0.0) {
                        ia = ar;
                        i14 = vstride + b_m;
                        for (ic = vstride; ic + 1 <= i14; ic++) {
                          ia++;
                          b_y->data[ic] += ky1->data[ib] * xk->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    vstride += b_m;
                  }
                }
              }

              i13 = e_a->size[0] * e_a->size[1];
              e_a->size[0] = a->size[0];
              e_a->size[1] = a->size[1];
              emxEnsureCapacity((emxArray__common *)e_a, i13, sizeof(double));
              loop_ub = a->size[0] * a->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                e_a->data[i13] = a->data[i13] + KY->data[i13] / (double)n;
              }

              mldivide(e_a, b_y, abi);
              if ((DD->size[0] == 0) || (DD->size[1] == 0)) {
                b_idx = 0;
              } else {
                trueCount = DD->size[0];
                b_idx = DD->size[1];
                if (trueCount > b_idx) {
                  b_idx = trueCount;
                }
              }

              if (b_idx == 1) {
                /* 2018-05-16 */
                eye((double)abi->size[1], a);
                pp = DD->data[0];
                i13 = DD->size[0] * DD->size[1];
                DD->size[0] = a->size[0];
                DD->size[1] = a->size[1];
                emxEnsureCapacity((emxArray__common *)DD, i13, sizeof(double));
                loop_ub = a->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = a->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    DD->data[i14 + DD->size[0] * i13] = a->data[i14 + a->size[0]
                      * i13] * pp;
                  }
                }
              }

              i13 = abi->size[1];
              if ((i13 == 1) || (DD->size[0] == 1)) {
                loop_ub = abi->size[1];
                i13 = c_abi->size[0] * c_abi->size[1];
                c_abi->size[0] = (int)p;
                c_abi->size[1] = loop_ub;
                emxEnsureCapacity((emxArray__common *)c_abi, i13, sizeof(double));
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = (int)p;
                  for (i14 = 0; i14 < br; i14++) {
                    c_abi->data[i14 + c_abi->size[0] * i13] = abi->data[i14 +
                      abi->size[0] * i13];
                  }
                }

                i13 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = c_abi->size[0];
                b_y->size[1] = DD->size[1];
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                loop_ub = c_abi->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = DD->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    b_y->data[i13 + b_y->size[0] * i14] = 0.0;
                    vstride = c_abi->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      b_y->data[i13 + b_y->size[0] * i14] += c_abi->data[i13 +
                        c_abi->size[0] * b_m] * DD->data[b_m + DD->size[0] * i14];
                    }
                  }
                }
              } else {
                i13 = abi->size[1];
                b_idx = DD->size[1];
                i14 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = (int)p;
                b_y->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)b_y, i14, sizeof(double));
                i14 = b_y->size[0] * b_y->size[1];
                emxEnsureCapacity((emxArray__common *)b_y, i14, sizeof(double));
                loop_ub = b_y->size[1];
                for (i14 = 0; i14 < loop_ub; i14++) {
                  br = b_y->size[0];
                  for (b_m = 0; b_m < br; b_m++) {
                    b_y->data[b_m + b_y->size[0] * i14] = 0.0;
                  }
                }

                if (DD->size[1] != 0) {
                  b_idx = (int)p * (DD->size[1] - 1);
                  for (vstride = 0; vstride <= b_idx; vstride += (int)p) {
                    i14 = vstride + (int)p;
                    for (ic = vstride; ic + 1 <= i14; ic++) {
                      b_y->data[ic] = 0.0;
                    }
                  }

                  br = 0;
                  for (vstride = 0; vstride <= b_idx; vstride += (int)p) {
                    ar = 0;
                    i14 = br + i13;
                    for (ib = br; ib + 1 <= i14; ib++) {
                      if (DD->data[ib] != 0.0) {
                        ia = ar;
                        b_m = vstride + (int)p;
                        for (ic = vstride; ic + 1 <= b_m; ic++) {
                          ia++;
                          b_y->data[ic] += DD->data[ib] * abi->data[(ia - 1) %
                            (int)p + abi->size[0] * ((ia - 1) / (int)p)];
                        }
                      }

                      ar += (int)p;
                    }

                    br += i13;
                  }
                }
              }

              loop_ub = abi->size[1];
              i13 = b->size[0] * b->size[1];
              b->size[0] = loop_ub;
              b->size[1] = (int)p;
              emxEnsureCapacity((emxArray__common *)b, i13, sizeof(double));
              br = (int)p;
              for (i13 = 0; i13 < br; i13++) {
                for (i14 = 0; i14 < loop_ub; i14++) {
                  b->data[i14 + b->size[0] * i13] = abi->data[i13 + abi->size[0]
                    * i14];
                }
              }

              if ((b_y->size[1] == 1) || (b->size[0] == 1)) {
                i13 = a->size[0] * a->size[1];
                a->size[0] = b_y->size[0];
                a->size[1] = b->size[1];
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                loop_ub = b_y->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = b->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    a->data[i13 + a->size[0] * i14] = 0.0;
                    vstride = b_y->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      a->data[i13 + a->size[0] * i14] += b_y->data[i13 +
                        b_y->size[0] * b_m] * b->data[b_m + b->size[0] * i14];
                    }
                  }
                }
              } else {
                k = b_y->size[1];
                vstride = b_y->size[0];
                b_idx = b->size[1];
                i13 = a->size[0] * a->size[1];
                a->size[0] = vstride;
                a->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                b_m = b_y->size[0];
                i13 = a->size[0] * a->size[1];
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                loop_ub = a->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = a->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    a->data[i14 + a->size[0] * i13] = 0.0;
                  }
                }

                b_idx = b_y->size[0] * (b->size[1] - 1);
                for (vstride = 0; vstride <= b_idx; vstride += b_m) {
                  i13 = vstride + b_m;
                  for (ic = vstride; ic + 1 <= i13; ic++) {
                    a->data[ic] = 0.0;
                  }
                }

                br = 0;
                for (vstride = 0; vstride <= b_idx; vstride += b_m) {
                  ar = 0;
                  i13 = br + k;
                  for (ib = br; ib + 1 <= i13; ib++) {
                    if (b->data[ib] != 0.0) {
                      ia = ar;
                      i14 = vstride + b_m;
                      for (ic = vstride; ic + 1 <= i14; ic++) {
                        ia++;
                        a->data[ic] += b->data[ib] * b_y->data[ia - 1];
                      }
                    }

                    ar += b_m;
                  }

                  br += k;
                }
              }

              i13 = U->size[0] * U->size[1];
              emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
              b_idx = U->size[0];
              vstride = U->size[1];
              loop_ub = b_idx * vstride;
              for (i13 = 0; i13 < loop_ub; i13++) {
                U->data[i13] += a->data[i13];
              }
            }

            eig(U, Vc, Dc);
            i13 = B->size[0] * B->size[1];
            B->size[0] = Vc->size[0];
            B->size[1] = Vc->size[1];
            emxEnsureCapacity((emxArray__common *)B, i13, sizeof(double));
            loop_ub = Vc->size[0] * Vc->size[1];
            for (i13 = 0; i13 < loop_ub; i13++) {
              B->data[i13] = Vc->data[i13].re;
            }

            /* Change for C */
            c_diag(Dc, R);
            i13 = dxij->size[0];
            dxij->size[0] = R->size[0];
            emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
            loop_ub = R->size[0];
            for (i13 = 0; i13 < loop_ub; i13++) {
              dxij->data[i13] = R->data[i13].re;
            }

            sort(dxij, iidx);
            i13 = D->size[0];
            D->size[0] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)D, i13, sizeof(double));
            loop_ub = iidx->size[0];
            for (i13 = 0; i13 < loop_ub; i13++) {
              D->data[i13] = iidx->data[i13];
            }

            loop_ub = B->size[0];
            i13 = U->size[0] * U->size[1];
            U->size[0] = loop_ub;
            U->size[1] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
            br = iidx->size[0];
            for (i13 = 0; i13 < br; i13++) {
              for (i14 = 0; i14 < loop_ub; i14++) {
                U->data[i14 + U->size[0] * i13] = B->data[i14 + B->size[0] *
                  (iidx->data[i13] - 1)];
              }
            }

            if (1.0 > ip) {
              loop_ub = 0;
            } else {
              loop_ub = (int)ip;
            }

            br = B->size[0];
            i13 = B->size[0] * B->size[1];
            B->size[0] = br;
            B->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)B, i13, sizeof(double));
            for (i13 = 0; i13 < loop_ub; i13++) {
              for (i14 = 0; i14 < br; i14++) {
                B->data[i14 + B->size[0] * i13] = U->data[i14 + U->size[0] * i13];
              }
            }
          } else {
            i13 = dd->size[0] * dd->size[1];
            dd->size[0] = (int)(m * p);
            dd->size[1] = (int)(m * p);
            emxEnsureCapacity((emxArray__common *)dd, i13, sizeof(double));
            loop_ub = (int)(m * p) * (int)(m * p);
            for (i13 = 0; i13 < loop_ub; i13++) {
              dd->data[i13] = 0.0;
            }

            i13 = dc->size[0];
            dc->size[0] = (int)(m * p);
            emxEnsureCapacity((emxArray__common *)dc, i13, sizeof(double));
            loop_ub = (int)(m * p);
            for (i13 = 0; i13 < loop_ub; i13++) {
              dc->data[i13] = 0.0;
            }

            i13 = d_D->size[0] * d_D->size[1];
            d_D->size[0] = (int)m;
            d_D->size[1] = (int)m;
            emxEnsureCapacity((emxArray__common *)d_D, i13, sizeof(double));
            loop_ub = (int)m * (int)m;
            for (i13 = 0; i13 < loop_ub; i13++) {
              d_D->data[i13] = 0.0;
            }

            /* Change for C */
            for (j = 0; j < exponent; j++) {
              loop_ub = x->size[1];
              i13 = b_xfast->size[0] * b_xfast->size[1];
              b_xfast->size[0] = 1;
              b_xfast->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_xfast, i13, sizeof(double));
              for (i13 = 0; i13 < loop_ub; i13++) {
                b_xfast->data[b_xfast->size[0] * i13] = xfast->data[j +
                  xfast->size[0] * i13];
              }

              repmat(b_xfast, (double)n, xij);
              i13 = xij->size[0] * xij->size[1];
              xij->size[0] = x->size[0];
              xij->size[1] = x->size[1];
              emxEnsureCapacity((emxArray__common *)xij, i13, sizeof(double));
              loop_ub = x->size[0] * x->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                xij->data[i13] = x->data[i13] - xij->data[i13];
              }

              loop_ub = V->size[1];
              i13 = b_Vfast->size[0] * b_Vfast->size[1];
              b_Vfast->size[0] = 1;
              b_Vfast->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_Vfast, i13, sizeof(double));
              for (i13 = 0; i13 < loop_ub; i13++) {
                b_Vfast->data[b_Vfast->size[0] * i13] = Vfast->data[j +
                  Vfast->size[0] * i13];
              }

              repmat(b_Vfast, (double)n, a);
              i13 = d_V->size[0] * d_V->size[1];
              d_V->size[0] = V->size[0];
              d_V->size[1] = V->size[1];
              emxEnsureCapacity((emxArray__common *)d_V, i13, sizeof(double));
              loop_ub = V->size[0] * V->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                d_V->data[i13] = V->data[i13] - a->data[i13];
              }

              power(d_V, a);
              sum(a, dxij);
              i13 = D->size[0];
              D->size[0] = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)D, i13, sizeof(double));
              loop_ub = dxij->size[0];
              for (i13 = 0; i13 < loop_ub; i13++) {
                D->data[i13] = dxij->data[i13];
              }

              d_sort(D);
              if ((h2 > D->data[(int)(2.0 * m) - 1]) || rtIsNaN(D->data[(int)
                   (2.0 * m) - 1])) {
                pp = h2;
              } else {
                pp = D->data[(int)(2.0 * m) - 1];
              }

              i13 = dxij->size[0];
              emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
              loop_ub = dxij->size[0];
              for (i13 = 0; i13 < loop_ub; i13++) {
                dxij->data[i13] = -dxij->data[i13] / pp;
              }

              b_exp(dxij);
              b_repmat(dxij, p + 1.0, U);
              if ((xij->size[1] == 1) || (B->size[0] == 1)) {
                i13 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = xij->size[0];
                b_y->size[1] = B->size[1];
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                loop_ub = xij->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = B->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    b_y->data[i13 + b_y->size[0] * i14] = 0.0;
                    vstride = xij->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      b_y->data[i13 + b_y->size[0] * i14] += xij->data[i13 +
                        xij->size[0] * b_m] * B->data[b_m + B->size[0] * i14];
                    }
                  }
                }
              } else {
                k = xij->size[1];
                vstride = xij->size[0];
                b_idx = B->size[1];
                i13 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = vstride;
                b_y->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                b_m = xij->size[0];
                i13 = b_y->size[0] * b_y->size[1];
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                loop_ub = b_y->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = b_y->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    b_y->data[i14 + b_y->size[0] * i13] = 0.0;
                  }
                }

                if ((xij->size[0] == 0) || (B->size[1] == 0)) {
                } else {
                  b_idx = xij->size[0] * (B->size[1] - 1);
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    i13 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i13; ic++) {
                      b_y->data[ic] = 0.0;
                    }

                    vstride += b_m;
                  }

                  br = 0;
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (B->data[ib] != 0.0) {
                        ia = ar;
                        i14 = vstride + b_m;
                        for (ic = vstride; ic + 1 <= i14; ic++) {
                          ia++;
                          b_y->data[ic] += B->data[ib] * xij->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    vstride += b_m;
                  }
                }
              }

              if (!((b_y->size[0] == 0) || (b_y->size[1] == 0))) {
                trueCount = b_y->size[0];
              } else if (!(n == 0)) {
                trueCount = n;
              } else {
                trueCount = b_y->size[0];
                if (!(trueCount > 0)) {
                  trueCount = 0;
                }
              }

              empty_non_axis_sizes = (trueCount == 0);
              if (empty_non_axis_sizes || (!((b_y->size[0] == 0) || (b_y->size[1]
                     == 0)))) {
                b_idx = b_y->size[1];
              } else {
                b_idx = 0;
              }

              if (empty_non_axis_sizes || (!(n == 0))) {
                vstride = 1;
              } else {
                vstride = 0;
              }

              i13 = c_reshapes[1].f1->size[0] * c_reshapes[1].f1->size[1];
              c_reshapes[1].f1->size[0] = trueCount;
              c_reshapes[1].f1->size[1] = vstride;
              emxEnsureCapacity((emxArray__common *)c_reshapes[1].f1, i13,
                                sizeof(double));
              loop_ub = trueCount * vstride;
              for (i13 = 0; i13 < loop_ub; i13++) {
                c_reshapes[1].f1->data[i13] = 1.0;
              }

              i13 = onexi->size[0] * onexi->size[1];
              onexi->size[0] = trueCount;
              onexi->size[1] = b_idx + c_reshapes[1].f1->size[1];
              emxEnsureCapacity((emxArray__common *)onexi, i13, sizeof(double));
              for (i13 = 0; i13 < b_idx; i13++) {
                for (i14 = 0; i14 < trueCount; i14++) {
                  onexi->data[i14 + onexi->size[0] * i13] = b_y->data[i14 +
                    trueCount * i13];
                }
              }

              loop_ub = c_reshapes[1].f1->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = c_reshapes[1].f1->size[0];
                for (i14 = 0; i14 < br; i14++) {
                  onexi->data[i14 + onexi->size[0] * (i13 + b_idx)] =
                    c_reshapes[1].f1->data[i14 + c_reshapes[1].f1->size[0] * i13];
                }
              }

              i13 = xk->size[0] * xk->size[1];
              xk->size[0] = onexi->size[1];
              xk->size[1] = onexi->size[0];
              emxEnsureCapacity((emxArray__common *)xk, i13, sizeof(double));
              loop_ub = onexi->size[0];
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = onexi->size[1];
                for (i14 = 0; i14 < br; i14++) {
                  xk->data[i14 + xk->size[0] * i13] = onexi->data[i13 +
                    onexi->size[0] * i14] * U->data[i13 + U->size[0] * i14];
                }
              }

              /* abi = inv(xk*onexi+eye(size(B,2)+1)/n)*(xk*ky1)*ky2; */
              if ((xk->size[1] == 1) || (onexi->size[0] == 1)) {
                i13 = a->size[0] * a->size[1];
                a->size[0] = xk->size[0];
                a->size[1] = onexi->size[1];
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                loop_ub = xk->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = onexi->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    a->data[i13 + a->size[0] * i14] = 0.0;
                    vstride = xk->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      a->data[i13 + a->size[0] * i14] += xk->data[i13 + xk->
                        size[0] * b_m] * onexi->data[b_m + onexi->size[0] * i14];
                    }
                  }
                }
              } else {
                k = xk->size[1];
                vstride = xk->size[0];
                b_idx = onexi->size[1];
                i13 = a->size[0] * a->size[1];
                a->size[0] = vstride;
                a->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                b_m = xk->size[0];
                i13 = a->size[0] * a->size[1];
                emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
                loop_ub = a->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = a->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    a->data[i14 + a->size[0] * i13] = 0.0;
                  }
                }

                if ((xk->size[0] == 0) || (onexi->size[1] == 0)) {
                } else {
                  b_idx = xk->size[0] * (onexi->size[1] - 1);
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    i13 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i13; ic++) {
                      a->data[ic] = 0.0;
                    }

                    vstride += b_m;
                  }

                  br = 0;
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (onexi->data[ib] != 0.0) {
                        ia = ar;
                        i14 = vstride + b_m;
                        for (ic = vstride; ic + 1 <= i14; ic++) {
                          ia++;
                          a->data[ic] += onexi->data[ib] * xk->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    vstride += b_m;
                  }
                }
              }

              eye((double)B->size[1] + 1.0, KY);
              if ((xk->size[1] == 1) || (ky1->size[0] == 1)) {
                i13 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = xk->size[0];
                b_y->size[1] = ky1->size[1];
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                loop_ub = xk->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = ky1->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    b_y->data[i13 + b_y->size[0] * i14] = 0.0;
                    vstride = xk->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      b_y->data[i13 + b_y->size[0] * i14] += xk->data[i13 +
                        xk->size[0] * b_m] * ky1->data[b_m + ky1->size[0] * i14];
                    }
                  }
                }
              } else {
                k = xk->size[1];
                vstride = xk->size[0];
                b_idx = ky1->size[1];
                i13 = b_y->size[0] * b_y->size[1];
                b_y->size[0] = vstride;
                b_y->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                b_m = xk->size[0];
                i13 = b_y->size[0] * b_y->size[1];
                emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
                loop_ub = b_y->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = b_y->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    b_y->data[i14 + b_y->size[0] * i13] = 0.0;
                  }
                }

                if ((xk->size[0] == 0) || (ky1->size[1] == 0)) {
                } else {
                  b_idx = xk->size[0] * (ky1->size[1] - 1);
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    i13 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i13; ic++) {
                      b_y->data[ic] = 0.0;
                    }

                    vstride += b_m;
                  }

                  br = 0;
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (ky1->data[ib] != 0.0) {
                        ia = ar;
                        i14 = vstride + b_m;
                        for (ic = vstride; ic + 1 <= i14; ic++) {
                          ia++;
                          b_y->data[ic] += ky1->data[ib] * xk->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    vstride += b_m;
                  }
                }
              }

              i13 = d_a->size[0] * d_a->size[1];
              d_a->size[0] = a->size[0];
              d_a->size[1] = a->size[1];
              emxEnsureCapacity((emxArray__common *)d_a, i13, sizeof(double));
              loop_ub = a->size[0] * a->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                d_a->data[i13] = a->data[i13] + KY->data[i13] / (double)n;
              }

              mldivide(d_a, b_y, a);
              if ((a->size[1] == 1) || (ky2->size[0] == 1)) {
                i13 = abi->size[0] * abi->size[1];
                abi->size[0] = a->size[0];
                abi->size[1] = ky2->size[1];
                emxEnsureCapacity((emxArray__common *)abi, i13, sizeof(double));
                loop_ub = a->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = ky2->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    abi->data[i13 + abi->size[0] * i14] = 0.0;
                    vstride = a->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      abi->data[i13 + abi->size[0] * i14] += a->data[i13 +
                        a->size[0] * b_m] * ky2->data[b_m + ky2->size[0] * i14];
                    }
                  }
                }
              } else {
                k = a->size[1];
                vstride = a->size[0];
                b_idx = ky2->size[1];
                i13 = abi->size[0] * abi->size[1];
                abi->size[0] = vstride;
                abi->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)abi, i13, sizeof(double));
                b_m = a->size[0];
                i13 = abi->size[0] * abi->size[1];
                emxEnsureCapacity((emxArray__common *)abi, i13, sizeof(double));
                loop_ub = abi->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = abi->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    abi->data[i14 + abi->size[0] * i13] = 0.0;
                  }
                }

                if ((a->size[0] == 0) || (ky2->size[1] == 0)) {
                } else {
                  b_idx = a->size[0] * (ky2->size[1] - 1);
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    i13 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i13; ic++) {
                      abi->data[ic] = 0.0;
                    }

                    vstride += b_m;
                  }

                  br = 0;
                  vstride = 0;
                  while ((b_m > 0) && (vstride <= b_idx)) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (ky2->data[ib] != 0.0) {
                        ia = ar;
                        i14 = vstride + b_m;
                        for (ic = vstride; ic + 1 <= i14; ic++) {
                          ia++;
                          abi->data[ic] += ky2->data[ib] * a->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    vstride += b_m;
                  }
                }
              }

              i13 = KY->size[0] * KY->size[1];
              KY->size[0] = xij->size[1];
              KY->size[1] = xij->size[0];
              emxEnsureCapacity((emxArray__common *)KY, i13, sizeof(double));
              loop_ub = xij->size[0];
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = xij->size[1];
                for (i14 = 0; i14 < br; i14++) {
                  KY->data[i14 + KY->size[0] * i13] = xij->data[i13 + xij->size
                    [0] * i14] * U->data[i13 + U->size[0] * i14];
                }
              }

              loop_ub = abi->size[1];
              i13 = b_abi->size[0] * b_abi->size[1];
              b_abi->size[0] = 1;
              b_abi->size[1] = loop_ub;
              emxEnsureCapacity((emxArray__common *)b_abi, i13, sizeof(double));
              for (i13 = 0; i13 < loop_ub; i13++) {
                b_abi->data[b_abi->size[0] * i13] = abi->data[((int)(m + 1.0) +
                  abi->size[0] * i13) - 1];
              }

              repmat(b_abi, (double)n, b);
              i13 = b->size[0] * b->size[1];
              b->size[0] = ky->size[0];
              b->size[1] = ky->size[1];
              emxEnsureCapacity((emxArray__common *)b, i13, sizeof(double));
              loop_ub = ky->size[0] * ky->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                b->data[i13] = ky->data[i13] - b->data[i13];
              }

              if ((KY->size[1] == 1) || (b->size[0] == 1)) {
                i13 = U->size[0] * U->size[1];
                U->size[0] = KY->size[0];
                U->size[1] = b->size[1];
                emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
                loop_ub = KY->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = b->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    U->data[i13 + U->size[0] * i14] = 0.0;
                    vstride = KY->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      U->data[i13 + U->size[0] * i14] += KY->data[i13 + KY->
                        size[0] * b_m] * b->data[b_m + b->size[0] * i14];
                    }
                  }
                }
              } else {
                k = KY->size[1];
                vstride = KY->size[0];
                b_idx = b->size[1];
                i13 = U->size[0] * U->size[1];
                U->size[0] = vstride;
                U->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
                b_m = KY->size[0];
                i13 = U->size[0] * U->size[1];
                emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
                loop_ub = U->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = U->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    U->data[i14 + U->size[0] * i13] = 0.0;
                  }
                }

                if (b->size[1] != 0) {
                  b_idx = KY->size[0] * (b->size[1] - 1);
                  for (vstride = 0; vstride <= b_idx; vstride += b_m) {
                    i13 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i13; ic++) {
                      U->data[ic] = 0.0;
                    }
                  }

                  br = 0;
                  for (vstride = 0; vstride <= b_idx; vstride += b_m) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (b->data[ib] != 0.0) {
                        ia = ar;
                        i14 = vstride + b_m;
                        for (ic = vstride; ic + 1 <= i14; ic++) {
                          ia++;
                          U->data[ic] += b->data[ib] * KY->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                  }
                }
              }

              if ((KY->size[1] == 1) || (xij->size[0] == 1)) {
                i13 = xk->size[0] * xk->size[1];
                xk->size[0] = KY->size[0];
                xk->size[1] = xij->size[1];
                emxEnsureCapacity((emxArray__common *)xk, i13, sizeof(double));
                loop_ub = KY->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = xij->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    xk->data[i13 + xk->size[0] * i14] = 0.0;
                    vstride = KY->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      xk->data[i13 + xk->size[0] * i14] += KY->data[i13 +
                        KY->size[0] * b_m] * xij->data[b_m + xij->size[0] * i14];
                    }
                  }
                }
              } else {
                k = KY->size[1];
                vstride = KY->size[0];
                b_idx = xij->size[1];
                i13 = xk->size[0] * xk->size[1];
                xk->size[0] = vstride;
                xk->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)xk, i13, sizeof(double));
                b_m = KY->size[0];
                i13 = xk->size[0] * xk->size[1];
                emxEnsureCapacity((emxArray__common *)xk, i13, sizeof(double));
                loop_ub = xk->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = xk->size[0];
                  for (i14 = 0; i14 < br; i14++) {
                    xk->data[i14 + xk->size[0] * i13] = 0.0;
                  }
                }

                if (xij->size[1] != 0) {
                  b_idx = KY->size[0] * (xij->size[1] - 1);
                  for (vstride = 0; vstride <= b_idx; vstride += b_m) {
                    i13 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i13; ic++) {
                      xk->data[ic] = 0.0;
                    }
                  }

                  br = 0;
                  for (vstride = 0; vstride <= b_idx; vstride += b_m) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (xij->data[ib] != 0.0) {
                        ia = ar;
                        i14 = vstride + b_m;
                        for (ic = vstride; ic + 1 <= i14; ic++) {
                          ia++;
                          xk->data[ic] += xij->data[ib] * KY->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                  }
                }
              }

              for (trueCount = 0; trueCount < (int)m; trueCount++) {
                pp = ((1.0 + (double)trueCount) - 1.0) * p + 1.0;
                e_x = (1.0 + (double)trueCount) * p;
                if (e_x < pp) {
                  i13 = meanky->size[0] * meanky->size[1];
                  meanky->size[0] = 1;
                  meanky->size[1] = 0;
                  emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof
                                    (double));
                } else if (pp == pp) {
                  i13 = meanky->size[0] * meanky->size[1];
                  meanky->size[0] = 1;
                  meanky->size[1] = (int)std::floor(e_x - pp) + 1;
                  emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof
                                    (double));
                  loop_ub = (int)std::floor(e_x - pp);
                  for (i13 = 0; i13 <= loop_ub; i13++) {
                    meanky->data[meanky->size[0] * i13] = pp + (double)i13;
                  }
                } else {
                  ndbl = std::floor((e_x - pp) + 0.5);
                  apnd = pp + ndbl;
                  cdiff = apnd - e_x;
                  absa = std::abs(pp);
                  absb = std::abs(e_x);
                  if ((absa > absb) || rtIsNaN(absb)) {
                    absb = absa;
                  }

                  if (std::abs(cdiff) < 4.4408920985006262E-16 * absb) {
                    ndbl++;
                    apnd = e_x;
                  } else if (cdiff > 0.0) {
                    apnd = pp + (ndbl - 1.0);
                  } else {
                    ndbl++;
                  }

                  if (ndbl >= 0.0) {
                    b_idx = (int)ndbl;
                  } else {
                    b_idx = 0;
                  }

                  i13 = meanky->size[0] * meanky->size[1];
                  meanky->size[0] = 1;
                  meanky->size[1] = b_idx;
                  emxEnsureCapacity((emxArray__common *)meanky, i13, sizeof
                                    (double));
                  if (b_idx > 0) {
                    meanky->data[0] = pp;
                    if (b_idx > 1) {
                      meanky->data[b_idx - 1] = apnd;
                      br = (b_idx - 1) / 2;
                      for (k = 1; k < br; k++) {
                        meanky->data[k] = pp + (double)k;
                        meanky->data[(b_idx - k) - 1] = apnd - (double)k;
                      }

                      if (br << 1 == b_idx - 1) {
                        meanky->data[br] = (pp + apnd) / 2.0;
                      } else {
                        meanky->data[br] = pp + (double)br;
                        meanky->data[br + 1] = apnd - (double)br;
                      }
                    }
                  }
                }

                loop_ub = abi->size[1];
                i13 = D->size[0];
                D->size[0] = loop_ub;
                emxEnsureCapacity((emxArray__common *)D, i13, sizeof(double));
                for (i13 = 0; i13 < loop_ub; i13++) {
                  D->data[i13] = abi->data[trueCount + abi->size[0] * i13];
                }

                if ((U->size[1] == 1) || (D->size[0] == 1)) {
                  i13 = dxij->size[0];
                  dxij->size[0] = U->size[0];
                  emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
                  loop_ub = U->size[0];
                  for (i13 = 0; i13 < loop_ub; i13++) {
                    dxij->data[i13] = 0.0;
                    br = U->size[1];
                    for (i14 = 0; i14 < br; i14++) {
                      dxij->data[i13] += U->data[i13 + U->size[0] * i14] *
                        D->data[i14];
                    }
                  }
                } else {
                  k = U->size[1];
                  U_idx_0 = (unsigned int)U->size[0];
                  i13 = dxij->size[0];
                  dxij->size[0] = (int)U_idx_0;
                  emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
                  b_m = U->size[0];
                  b_idx = dxij->size[0];
                  i13 = dxij->size[0];
                  dxij->size[0] = b_idx;
                  emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
                  for (i13 = 0; i13 < b_idx; i13++) {
                    dxij->data[i13] = 0.0;
                  }

                  vstride = 0;
                  while (vstride <= 0) {
                    for (ic = 1; ic <= b_m; ic++) {
                      dxij->data[ic - 1] = 0.0;
                    }

                    vstride = b_m;
                  }

                  br = 0;
                  vstride = 0;
                  while (vstride <= 0) {
                    ar = 0;
                    i13 = br + k;
                    for (ib = br; ib + 1 <= i13; ib++) {
                      if (D->data[ib] != 0.0) {
                        ia = ar;
                        for (ic = 0; ic + 1 <= b_m; ic++) {
                          ia++;
                          dxij->data[ic] += D->data[ib] * U->data[ia - 1];
                        }
                      }

                      ar += b_m;
                    }

                    br += k;
                    vstride = b_m;
                  }
                }

                i13 = r4->size[0] * r4->size[1];
                r4->size[0] = 1;
                r4->size[1] = meanky->size[1];
                emxEnsureCapacity((emxArray__common *)r4, i13, sizeof(int));
                loop_ub = meanky->size[0] * meanky->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  r4->data[i13] = (int)meanky->data[i13];
                }

                i13 = b_dc->size[0];
                b_dc->size[0] = meanky->size[1];
                emxEnsureCapacity((emxArray__common *)b_dc, i13, sizeof(double));
                loop_ub = meanky->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  b_dc->data[i13] = dc->data[(int)meanky->data[meanky->size[0] *
                    i13] - 1] + dxij->data[i13];
                }

                loop_ub = r4->size[1];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  dc->data[r4->data[r4->size[0] * i13] - 1] = b_dc->data[(*(int
                    (*)[2])r4->size)[0] * i13];
                }
              }

              if (1.0 > m) {
                loop_ub = 0;
                br = 0;
              } else {
                loop_ub = (int)m;
                br = (int)m;
              }

              vstride = abi->size[1];
              i13 = a->size[0] * a->size[1];
              a->size[0] = loop_ub;
              a->size[1] = vstride;
              emxEnsureCapacity((emxArray__common *)a, i13, sizeof(double));
              for (i13 = 0; i13 < vstride; i13++) {
                for (i14 = 0; i14 < loop_ub; i14++) {
                  a->data[i14 + a->size[0] * i13] = abi->data[i14 + abi->size[0]
                    * i13];
                }
              }

              vstride = abi->size[1];
              i13 = b->size[0] * b->size[1];
              b->size[0] = vstride;
              b->size[1] = br;
              emxEnsureCapacity((emxArray__common *)b, i13, sizeof(double));
              for (i13 = 0; i13 < br; i13++) {
                for (i14 = 0; i14 < vstride; i14++) {
                  b->data[i14 + b->size[0] * i13] = abi->data[i13 + abi->size[0]
                    * i14];
                }
              }

              i13 = abi->size[1];
              if ((i13 == 1) || (b->size[0] == 1)) {
                i13 = U->size[0] * U->size[1];
                U->size[0] = a->size[0];
                U->size[1] = b->size[1];
                emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
                loop_ub = a->size[0];
                for (i13 = 0; i13 < loop_ub; i13++) {
                  br = b->size[1];
                  for (i14 = 0; i14 < br; i14++) {
                    U->data[i13 + U->size[0] * i14] = 0.0;
                    vstride = a->size[1];
                    for (b_m = 0; b_m < vstride; b_m++) {
                      U->data[i13 + U->size[0] * i14] += a->data[i13 + a->size[0]
                        * b_m] * b->data[b_m + b->size[0] * i14];
                    }
                  }
                }
              } else {
                i13 = abi->size[1];
                b_idx = b->size[1];
                i14 = U->size[0] * U->size[1];
                U->size[0] = loop_ub;
                U->size[1] = b_idx;
                emxEnsureCapacity((emxArray__common *)U, i14, sizeof(double));
                i14 = U->size[0] * U->size[1];
                emxEnsureCapacity((emxArray__common *)U, i14, sizeof(double));
                br = U->size[1];
                for (i14 = 0; i14 < br; i14++) {
                  vstride = U->size[0];
                  for (b_m = 0; b_m < vstride; b_m++) {
                    U->data[b_m + U->size[0] * i14] = 0.0;
                  }
                }

                if ((loop_ub == 0) || (b->size[1] == 0)) {
                } else {
                  b_idx = loop_ub * (b->size[1] - 1);
                  vstride = 0;
                  while ((loop_ub > 0) && (vstride <= b_idx)) {
                    i14 = vstride + loop_ub;
                    for (ic = vstride; ic + 1 <= i14; ic++) {
                      U->data[ic] = 0.0;
                    }

                    vstride += loop_ub;
                  }

                  br = 0;
                  vstride = 0;
                  while ((loop_ub > 0) && (vstride <= b_idx)) {
                    ar = 0;
                    i14 = br + i13;
                    for (ib = br; ib + 1 <= i14; ib++) {
                      if (b->data[ib] != 0.0) {
                        ia = ar;
                        b_m = vstride + loop_ub;
                        for (ic = vstride; ic + 1 <= b_m; ic++) {
                          ia++;
                          U->data[ic] += b->data[ib] * a->data[ia - 1];
                        }
                      }

                      ar += loop_ub;
                    }

                    br += i13;
                    vstride += loop_ub;
                  }
                }
              }

              kron(U, xk, a);
              i13 = dd->size[0] * dd->size[1];
              emxEnsureCapacity((emxArray__common *)dd, i13, sizeof(double));
              b_idx = dd->size[0];
              vstride = dd->size[1];
              loop_ub = b_idx * vstride;
              for (i13 = 0; i13 < loop_ub; i13++) {
                dd->data[i13] += a->data[i13];
              }

              i13 = d_D->size[0] * d_D->size[1];
              emxEnsureCapacity((emxArray__common *)d_D, i13, sizeof(double));
              b_idx = d_D->size[0];
              vstride = d_D->size[1];
              loop_ub = b_idx * vstride;
              for (i13 = 0; i13 < loop_ub; i13++) {
                d_D->data[i13] += U->data[i13];
              }
            }

            eye((double)dc->size[0], KY);
            i13 = b_dd->size[0] * b_dd->size[1];
            b_dd->size[0] = dd->size[0];
            b_dd->size[1] = dd->size[1];
            emxEnsureCapacity((emxArray__common *)b_dd, i13, sizeof(double));
            loop_ub = dd->size[0] * dd->size[1];
            for (i13 = 0; i13 < loop_ub; i13++) {
              b_dd->data[i13] = dd->data[i13] + KY->data[i13] / (double)n;
            }

            b_mldivide(b_dd, dc);
            if (((int)m == 1) || (d_D->size[0] == 1)) {
              i13 = b_y->size[0] * b_y->size[1];
              b_y->size[0] = (int)p;
              b_y->size[1] = d_D->size[1];
              emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
              loop_ub = (int)p;
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = d_D->size[1];
                for (i14 = 0; i14 < br; i14++) {
                  b_y->data[i13 + b_y->size[0] * i14] = 0.0;
                  vstride = (int)m;
                  for (b_m = 0; b_m < vstride; b_m++) {
                    b_y->data[i13 + b_y->size[0] * i14] += dc->data[i13 + (int)p
                      * b_m] * d_D->data[b_m + d_D->size[0] * i14];
                  }
                }
              }
            } else {
              b_idx = d_D->size[1];
              i13 = b_y->size[0] * b_y->size[1];
              b_y->size[0] = (int)p;
              b_y->size[1] = b_idx;
              emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
              i13 = b_y->size[0] * b_y->size[1];
              emxEnsureCapacity((emxArray__common *)b_y, i13, sizeof(double));
              loop_ub = b_y->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = b_y->size[0];
                for (i14 = 0; i14 < br; i14++) {
                  b_y->data[i14 + b_y->size[0] * i13] = 0.0;
                }
              }

              if (d_D->size[1] != 0) {
                b_idx = (int)p * (d_D->size[1] - 1);
                for (vstride = 0; vstride <= b_idx; vstride += (int)p) {
                  i13 = vstride + (int)p;
                  for (ic = vstride; ic + 1 <= i13; ic++) {
                    b_y->data[ic] = 0.0;
                  }
                }

                br = 0;
                for (vstride = 0; vstride <= b_idx; vstride += (int)p) {
                  ar = 0;
                  i13 = br + (int)m;
                  for (ib = br; ib + 1 <= i13; ib++) {
                    if (d_D->data[ib] != 0.0) {
                      ia = ar;
                      i14 = vstride + (int)p;
                      for (ic = vstride; ic + 1 <= i14; ic++) {
                        ia++;
                        b_y->data[ic] += d_D->data[ib] * dc->data[ia - 1];
                      }
                    }

                    ar += (int)p;
                  }

                  br += (int)m;
                }
              }
            }

            i13 = b->size[0] * b->size[1];
            b->size[0] = (int)m;
            b->size[1] = (int)p;
            emxEnsureCapacity((emxArray__common *)b, i13, sizeof(double));
            loop_ub = (int)p;
            for (i13 = 0; i13 < loop_ub; i13++) {
              br = (int)m;
              for (i14 = 0; i14 < br; i14++) {
                b->data[i14 + b->size[0] * i13] = dc->data[i13 + (int)p * i14];
              }
            }

            if ((b_y->size[1] == 1) || (b->size[0] == 1)) {
              i13 = U->size[0] * U->size[1];
              U->size[0] = b_y->size[0];
              U->size[1] = b->size[1];
              emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
              loop_ub = b_y->size[0];
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = b->size[1];
                for (i14 = 0; i14 < br; i14++) {
                  U->data[i13 + U->size[0] * i14] = 0.0;
                  vstride = b_y->size[1];
                  for (b_m = 0; b_m < vstride; b_m++) {
                    U->data[i13 + U->size[0] * i14] += b_y->data[i13 + b_y->
                      size[0] * b_m] * b->data[b_m + b->size[0] * i14];
                  }
                }
              }
            } else {
              k = b_y->size[1];
              vstride = b_y->size[0];
              b_idx = b->size[1];
              i13 = U->size[0] * U->size[1];
              U->size[0] = vstride;
              U->size[1] = b_idx;
              emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
              b_m = b_y->size[0];
              i13 = U->size[0] * U->size[1];
              emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
              loop_ub = U->size[1];
              for (i13 = 0; i13 < loop_ub; i13++) {
                br = U->size[0];
                for (i14 = 0; i14 < br; i14++) {
                  U->data[i14 + U->size[0] * i13] = 0.0;
                }
              }

              b_idx = b_y->size[0] * (b->size[1] - 1);
              for (vstride = 0; vstride <= b_idx; vstride += b_m) {
                i13 = vstride + b_m;
                for (ic = vstride; ic + 1 <= i13; ic++) {
                  U->data[ic] = 0.0;
                }
              }

              br = 0;
              for (vstride = 0; vstride <= b_idx; vstride += b_m) {
                ar = 0;
                i13 = br + k;
                for (ib = br; ib + 1 <= i13; ib++) {
                  if (b->data[ib] != 0.0) {
                    ia = ar;
                    i14 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i14; ic++) {
                      ia++;
                      U->data[ic] += b->data[ib] * b_y->data[ia - 1];
                    }
                  }

                  ar += b_m;
                }

                br += k;
              }
            }

            eig(U, Vc, Dc);
            i13 = B->size[0] * B->size[1];
            B->size[0] = Vc->size[0];
            B->size[1] = Vc->size[1];
            emxEnsureCapacity((emxArray__common *)B, i13, sizeof(double));
            loop_ub = Vc->size[0] * Vc->size[1];
            for (i13 = 0; i13 < loop_ub; i13++) {
              B->data[i13] = Vc->data[i13].re;
            }

            /* Change for C */
            c_diag(Dc, R);
            i13 = dxij->size[0];
            dxij->size[0] = R->size[0];
            emxEnsureCapacity((emxArray__common *)dxij, i13, sizeof(double));
            loop_ub = R->size[0];
            for (i13 = 0; i13 < loop_ub; i13++) {
              dxij->data[i13] = R->data[i13].re;
            }

            sort(dxij, iidx);
            i13 = D->size[0];
            D->size[0] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)D, i13, sizeof(double));
            loop_ub = iidx->size[0];
            for (i13 = 0; i13 < loop_ub; i13++) {
              D->data[i13] = iidx->data[i13];
            }

            loop_ub = B->size[0];
            i13 = U->size[0] * U->size[1];
            U->size[0] = loop_ub;
            U->size[1] = iidx->size[0];
            emxEnsureCapacity((emxArray__common *)U, i13, sizeof(double));
            br = iidx->size[0];
            for (i13 = 0; i13 < br; i13++) {
              for (i14 = 0; i14 < loop_ub; i14++) {
                U->data[i14 + U->size[0] * i13] = B->data[i14 + B->size[0] *
                  (iidx->data[i13] - 1)];
              }
            }

            if (1.0 > ip) {
              loop_ub = 0;
            } else {
              loop_ub = (int)ip;
            }

            br = B->size[0];
            i13 = B->size[0] * B->size[1];
            B->size[0] = br;
            B->size[1] = loop_ub;
            emxEnsureCapacity((emxArray__common *)B, i13, sizeof(double));
            for (i13 = 0; i13 < loop_ub; i13++) {
              for (i14 = 0; i14 < br; i14++) {
                B->data[i14 + B->size[0] * i13] = U->data[i14 + U->size[0] * i13];
              }
            }
          }
        }

        b_idx = (int)SEQ->data[iter_ip];
        loop_ub = B->size[1];
        for (i13 = 0; i13 < loop_ub; i13++) {
          br = B->size[0];
          for (i14 = 0; i14 < br; i14++) {
            BB->data[(i14 + BB->size[0] * i13) + BB->size[0] * BB->size[1] *
              (b_idx - 1)] = B->data[i14 + B->size[0] * i13];
          }
        }

        if (d_strcmp(method)) {
          i13 = (int)((1.0 + (-1.0 - (p - 1.0))) / -1.0);
          for (b_idx = 0; b_idx < i13; b_idx++) {
            vstride = ((int)p - b_idx) - 1;
            if (1 > vstride) {
              loop_ub = -1;
            } else {
              loop_ub = vstride - 1;
            }

            br = U->size[0] - 1;
            for (i14 = 0; i14 <= loop_ub; i14++) {
              for (b_m = 0; b_m <= br; b_m++) {
                BB->data[(b_m + BB->size[0] * i14) + BB->size[0] * BB->size[1] *
                  (vstride - 1)] = U->data[b_m + U->size[0] * i14];
              }
            }
          }

          exitg2 = true;
        } else {
          /* B = B(:,1:(ip-1)); */
          iter_ip++;
        }
      }

      emxFree_int32_T(&e_Ifast);
      emxFree_int32_T(&d_Ifast);
      emxFree_real_T(&c_abi);
      emxFree_real_T(&b_dc);
      emxFree_int32_T(&c_Ifast);
      emxFree_int32_T(&b_Ifast);
      emxFree_real_T(&b_B);
      emxFree_real_T(&e_V);
      emxFree_real_T(&e_a);
      emxFree_real_T(&b_xfast);
      emxFree_real_T(&b_Vfast);
      emxFree_real_T(&d_V);
      emxFree_real_T(&d_a);
      emxFree_real_T(&b_abi);
      emxFree_real_T(&b_dd);
      emxFreeMatrix_cell_wrap_0(c_reshapes);
      emxFree_int32_T(&r4);
      emxFree_creal_T(&R);
      emxFree_real_T(&d_D);
      emxFree_real_T(&dc);
      emxFree_real_T(&dd);
      emxFree_real_T(&abi);
      emxFree_real_T(&xk);
      emxFree_real_T(&xij);
      emxFree_real_T(&xfast);
      emxFree_real_T(&Vfast);
      emxFree_real_T(&Ifast);
      emxFree_real_T(&SEQ);
    }

    emxFree_real_T(&onexi);

    /* for ip = 1:p */
    /*     B = ss*BB(:,:,ip); */
    /*     for i = 1:size(B,2) */
    /*         B(:,i) = B(:,i)/(norm(B(:,i))+1.0e-20); */
    /*     end */
    /*     BB(:,:,ip) = B; */
    /* end */
    if (1.0 > max_dim) {
      loop_ub = 0;
    } else {
      loop_ub = (int)max_dim;
    }

    if (1.0 > max_dim) {
      br = 0;
    } else {
      br = (int)max_dim;
    }

    emxInit_real_T2(&b_BB, 3);
    b_idx = BB->size[0];
    i13 = b_BB->size[0] * b_BB->size[1] * b_BB->size[2];
    b_BB->size[0] = b_idx;
    b_BB->size[1] = loop_ub;
    b_BB->size[2] = br;
    emxEnsureCapacity((emxArray__common *)b_BB, i13, sizeof(double));
    for (i13 = 0; i13 < br; i13++) {
      for (i14 = 0; i14 < loop_ub; i14++) {
        for (b_m = 0; b_m < b_idx; b_m++) {
          b_BB->data[(b_m + b_BB->size[0] * i14) + b_BB->size[0] * b_BB->size[1]
            * i13] = BB->data[(b_m + BB->size[0] * i14) + BB->size[0] * BB->
            size[1] * i13];
        }
      }
    }

    i13 = BB->size[0] * BB->size[1] * BB->size[2];
    BB->size[0] = b_BB->size[0];
    BB->size[1] = b_BB->size[1];
    BB->size[2] = b_BB->size[2];
    emxEnsureCapacity((emxArray__common *)BB, i13, sizeof(double));
    loop_ub = b_BB->size[2];
    for (i13 = 0; i13 < loop_ub; i13++) {
      br = b_BB->size[1];
      for (i14 = 0; i14 < br; i14++) {
        vstride = b_BB->size[0];
        for (b_m = 0; b_m < vstride; b_m++) {
          BB->data[(b_m + BB->size[0] * i14) + BB->size[0] * BB->size[1] * i13] =
            b_BB->data[(b_m + b_BB->size[0] * i14) + b_BB->size[0] * b_BB->size
            [1] * i13];
        }
      }
    }

    emxFree_real_T(&b_BB);

    /* 2018-05-16 */
    i13 = BBvs->size[0] * BBvs->size[1] * BBvs->size[2];
    BBvs->size[0] = BB->size[0];
    BBvs->size[1] = BB->size[1];
    BBvs->size[2] = BB->size[2];
    emxEnsureCapacity((emxArray__common *)BBvs, i13, sizeof(double));
    loop_ub = BB->size[0] * BB->size[1] * BB->size[2];
    for (i13 = 0; i13 < loop_ub; i13++) {
      BBvs->data[i13] = BB->data[i13];
    }

    emxInit_real_T2(&BB0, 3);

    /* 2018-05-16 */
    i13 = BB0->size[0] * BB0->size[1] * BB0->size[2];
    BB0->size[0] = p0;
    BB0->size[1] = (int)max_dim;
    BB0->size[2] = (int)max_dim;
    emxEnsureCapacity((emxArray__common *)BB0, i13, sizeof(double));
    loop_ub = p0 * (int)max_dim * (int)max_dim;
    for (i13 = 0; i13 < loop_ub; i13++) {
      BB0->data[i13] = 0.0;
    }

    /* 2018-05-16 */
    trueCount = 0;
    emxInit_real_T(&c_B, 1);
    emxInit_real_T(&d_B, 1);
    while (trueCount <= (int)max_dim - 1) {
      /* 2018-05-16 */
      loop_ub = BB->size[0];
      br = BB->size[1];
      i13 = b->size[0] * b->size[1];
      b->size[0] = loop_ub;
      b->size[1] = br;
      emxEnsureCapacity((emxArray__common *)b, i13, sizeof(double));
      for (i13 = 0; i13 < br; i13++) {
        for (i14 = 0; i14 < loop_ub; i14++) {
          b->data[i14 + b->size[0] * i13] = BB->data[(i14 + BB->size[0] * i13) +
            BB->size[0] * BB->size[1] * trueCount];
        }
      }

      guard1 = false;
      if (ss->size[1] == 1) {
        guard1 = true;
      } else {
        i13 = BB->size[0];
        if (i13 == 1) {
          guard1 = true;
        } else {
          k = ss->size[1];
          i13 = BB->size[1];
          vstride = ss->size[0];
          i14 = B->size[0] * B->size[1];
          B->size[0] = vstride;
          B->size[1] = i13;
          emxEnsureCapacity((emxArray__common *)B, i14, sizeof(double));
          b_m = ss->size[0];
          i13 = B->size[0] * B->size[1];
          emxEnsureCapacity((emxArray__common *)B, i13, sizeof(double));
          loop_ub = B->size[1];
          for (i13 = 0; i13 < loop_ub; i13++) {
            br = B->size[0];
            for (i14 = 0; i14 < br; i14++) {
              B->data[i14 + B->size[0] * i13] = 0.0;
            }
          }

          if (ss->size[0] == 0) {
          } else {
            i13 = BB->size[1];
            if (i13 == 0) {
            } else {
              i13 = BB->size[1] - 1;
              b_idx = ss->size[0] * i13;
              vstride = 0;
              while ((b_m > 0) && (vstride <= b_idx)) {
                i13 = vstride + b_m;
                for (ic = vstride; ic + 1 <= i13; ic++) {
                  B->data[ic] = 0.0;
                }

                vstride += b_m;
              }

              br = 0;
              vstride = 0;
              while ((b_m > 0) && (vstride <= b_idx)) {
                ar = 0;
                i13 = br + k;
                for (ib = br; ib + 1 <= i13; ib++) {
                  if (b->data[ib] != 0.0) {
                    ia = ar;
                    i14 = vstride + b_m;
                    for (ic = vstride; ic + 1 <= i14; ic++) {
                      ia++;
                      B->data[ic] += b->data[ib] * ss->data[ia - 1];
                    }
                  }

                  ar += b_m;
                }

                br += k;
                vstride += b_m;
              }
            }
          }
        }
      }

      if (guard1) {
        i13 = B->size[0] * B->size[1];
        B->size[0] = ss->size[0];
        B->size[1] = b->size[1];
        emxEnsureCapacity((emxArray__common *)B, i13, sizeof(double));
        loop_ub = ss->size[0];
        for (i13 = 0; i13 < loop_ub; i13++) {
          br = b->size[1];
          for (i14 = 0; i14 < br; i14++) {
            B->data[i13 + B->size[0] * i14] = 0.0;
            vstride = ss->size[1];
            for (b_m = 0; b_m < vstride; b_m++) {
              B->data[i13 + B->size[0] * i14] += ss->data[i13 + ss->size[0] *
                b_m] * b->data[b_m + b->size[0] * i14];
            }
          }
        }
      }

      i13 = B->size[1];
      for (ar = 0; ar < i13; ar++) {
        loop_ub = B->size[0];
        i14 = c_B->size[0];
        c_B->size[0] = loop_ub;
        emxEnsureCapacity((emxArray__common *)c_B, i14, sizeof(double));
        for (i14 = 0; i14 < loop_ub; i14++) {
          c_B->data[i14] = B->data[i14 + B->size[0] * ar];
        }

        pp = norm(c_B) + 1.0E-20;
        b_idx = B->size[0];
        i14 = d_B->size[0];
        d_B->size[0] = b_idx;
        emxEnsureCapacity((emxArray__common *)d_B, i14, sizeof(double));
        for (i14 = 0; i14 < b_idx; i14++) {
          d_B->data[i14] = B->data[i14 + B->size[0] * ar] / pp;
        }

        loop_ub = d_B->size[0];
        for (i14 = 0; i14 < loop_ub; i14++) {
          B->data[i14 + B->size[0] * ar] = d_B->data[i14];
        }
      }

      i13 = iidx->size[0];
      iidx->size[0] = (int)p;
      emxEnsureCapacity((emxArray__common *)iidx, i13, sizeof(int));
      loop_ub = (int)p;
      for (i13 = 0; i13 < loop_ub; i13++) {
        iidx->data[i13] = (int)sele->data[i13] - 1;
      }

      loop_ub = B->size[1];
      for (i13 = 0; i13 < loop_ub; i13++) {
        br = B->size[0];
        for (i14 = 0; i14 < br; i14++) {
          BB0->data[(iidx->data[i14] + BB0->size[0] * i13) + BB0->size[0] *
            BB0->size[1] * trueCount] = B->data[i14 + B->size[0] * i13];
        }
      }

      trueCount++;
    }

    emxFree_real_T(&d_B);
    emxFree_real_T(&c_B);
    emxFree_real_T(&B);
    i13 = BB->size[0] * BB->size[1] * BB->size[2];
    BB->size[0] = BB0->size[0];
    BB->size[1] = BB0->size[1];
    BB->size[2] = BB0->size[2];
    emxEnsureCapacity((emxArray__common *)BB, i13, sizeof(double));
    loop_ub = BB0->size[0] * BB0->size[1] * BB0->size[2];
    for (i13 = 0; i13 < loop_ub; i13++) {
      BB->data[i13] = BB0->data[i13];
    }

    emxFree_real_T(&BB0);

    /* x=x*inv(ss); */
    /*  for ip = 1:p */
    /*       cv(ip) = CVm(x*BB(:,1:ip, ip), ky); */
    /*  end */
    /* BB1D = reshape(BB0,p0*p*p,1); */
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&b_y);
  emxFree_creal_T(&Dc);
  emxFree_creal_T(&Vc);
  emxFree_real_T(&dxij);
  emxFree_real_T(&U);
  emxFree_real_T(&meanky);
  emxFree_real_T(&KY);
  emxFree_real_T(&DD);
  emxFree_real_T(&ky2);
  emxFree_real_T(&ky1);
  emxFree_real_T(&D);
  emxFree_real_T(&V);
  emxFree_real_T(&ss);
  emxFree_real_T(&sele);
}

/* End of code generation (MAVEfast.cpp) */
