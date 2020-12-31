/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * prctile.cpp
 *
 * Code generation for function 'prctile'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "prctile.h"
#include "MAVEfast_emxutil.h"

/* Function Declarations */
static double rt_roundd_snf(double u);

/* Function Definitions */
static double rt_roundd_snf(double u)
{
  double y;
  if (std::abs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = std::floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = std::ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

void prctile(const emxArray_real_T *x, const emxArray_real_T *p, emxArray_real_T
             *y)
{
  int i;
  emxArray_int32_T *idx;
  emxArray_int32_T *iwork;
  int n;
  unsigned int unnamed_idx_0;
  int i2;
  int k;
  boolean_T b_p;
  int j;
  int pEnd;
  int c_p;
  int q;
  int qEnd;
  double r;
  double b_i;
  int kEnd;
  i = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = p->size[1];
  emxEnsureCapacity((emxArray__common *)y, i, sizeof(double));
  emxInit_int32_T(&idx, 1);
  emxInit_int32_T(&iwork, 1);
  if ((x->size[0] == 0) || (p->size[1] == 0)) {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i, sizeof(double));
    i2 = y->size[1];
    for (i = 0; i < i2; i++) {
      y->data[y->size[0] * i] = rtNaN;
    }
  } else {
    n = x->size[0] + 1;
    unnamed_idx_0 = (unsigned int)x->size[0];
    i = idx->size[0];
    idx->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)idx, i, sizeof(int));
    i2 = (int)unnamed_idx_0;
    for (i = 0; i < i2; i++) {
      idx->data[i] = 0;
    }

    i = iwork->size[0];
    iwork->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)iwork, i, sizeof(int));
    for (k = 1; k <= n - 2; k += 2) {
      if ((x->data[k - 1] <= x->data[k]) || rtIsNaN(x->data[k])) {
        b_p = true;
      } else {
        b_p = false;
      }

      if (b_p) {
        idx->data[k - 1] = k;
        idx->data[k] = k + 1;
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    if ((x->size[0] & 1) != 0) {
      idx->data[x->size[0] - 1] = x->size[0];
    }

    i = 2;
    while (i < n - 1) {
      i2 = i << 1;
      j = 1;
      for (pEnd = 1 + i; pEnd < n; pEnd = qEnd + i) {
        c_p = j;
        q = pEnd - 1;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if ((x->data[idx->data[c_p - 1] - 1] <= x->data[idx->data[q] - 1]) ||
              rtIsNaN(x->data[idx->data[q] - 1])) {
            b_p = true;
          } else {
            b_p = false;
          }

          if (b_p) {
            iwork->data[k] = idx->data[c_p - 1];
            c_p++;
            if (c_p == pEnd) {
              while (q + 1 < qEnd) {
                k++;
                iwork->data[k] = idx->data[q];
                q++;
              }
            }
          } else {
            iwork->data[k] = idx->data[q];
            q++;
            if (q + 1 == qEnd) {
              while (c_p < pEnd) {
                k++;
                iwork->data[k] = idx->data[c_p - 1];
                c_p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= kEnd; k++) {
          idx->data[(j + k) - 1] = iwork->data[k];
        }

        j = qEnd;
      }

      i = i2;
    }

    i = x->size[0];
    while ((i > 0) && rtIsNaN(x->data[idx->data[i - 1] - 1])) {
      i--;
    }

    if (i < 1) {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)y, i, sizeof(double));
      i2 = y->size[1];
      for (i = 0; i < i2; i++) {
        y->data[y->size[0] * i] = rtNaN;
      }
    } else if (i == 1) {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)y, i, sizeof(double));
      i2 = y->size[1];
      for (i = 0; i < i2; i++) {
        y->data[y->size[0] * i] = x->data[idx->data[0] - 1];
      }
    } else {
      for (k = 0; k < p->size[1]; k++) {
        r = p->data[k] / 100.0 * (double)i;
        b_i = rt_roundd_snf(r);
        if (b_i < 1.0) {
          y->data[k] = x->data[idx->data[0] - 1];
        } else if (i <= b_i) {
          y->data[k] = x->data[idx->data[i - 1] - 1];
        } else {
          r -= b_i;
          y->data[k] = (0.5 - r) * x->data[idx->data[(int)b_i - 1] - 1] + (0.5 +
            r) * x->data[idx->data[(int)(b_i + 1.0) - 1] - 1];
        }
      }
    }
  }

  emxFree_int32_T(&iwork);
  emxFree_int32_T(&idx);
}

/* End of code generation (prctile.cpp) */
