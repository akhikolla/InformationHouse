/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * mldivide.cpp
 *
 * Code generation for function 'mldivide'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "mldivide.h"
#include "MAVEfast_emxutil.h"
#include "xgetrf.h"
#include "qrsolve.h"
#include "xgeqp3.h"

/* Function Definitions */
void b_mldivide(const emxArray_real_T *A, emxArray_real_T *B)
{
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  emxArray_real_T *b_B;
  unsigned int unnamed_idx_0;
  int m;
  int mn;
  int loop_ub;
  int rankA;
  double wj;
  int i;
  emxInit_real_T1(&b_A, 2);
  emxInit_real_T(&tau, 1);
  emxInit_int32_T1(&jpvt, 2);
  emxInit_real_T(&b_B, 1);
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0)) {
    unnamed_idx_0 = (unsigned int)A->size[1];
    m = B->size[0];
    B->size[0] = (int)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)B, m, sizeof(double));
    loop_ub = (int)unnamed_idx_0;
    for (m = 0; m < loop_ub; m++) {
      B->data[m] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    mn = A->size[1];
    m = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, m, sizeof(double));
    loop_ub = A->size[0] * A->size[1];
    for (m = 0; m < loop_ub; m++) {
      b_A->data[m] = A->data[m];
    }

    xgetrf(A->size[1], A->size[1], b_A, A->size[1], jpvt, &loop_ub);
    for (loop_ub = 0; loop_ub + 1 < mn; loop_ub++) {
      if (jpvt->data[loop_ub] != loop_ub + 1) {
        wj = B->data[loop_ub];
        B->data[loop_ub] = B->data[jpvt->data[loop_ub] - 1];
        B->data[jpvt->data[loop_ub] - 1] = wj;
      }
    }

    for (loop_ub = 0; loop_ub + 1 <= mn; loop_ub++) {
      m = mn * loop_ub;
      if (B->data[loop_ub] != 0.0) {
        for (i = loop_ub + 1; i + 1 <= mn; i++) {
          B->data[i] -= B->data[loop_ub] * b_A->data[i + m];
        }
      }
    }

    for (loop_ub = A->size[1] - 1; loop_ub + 1 > 0; loop_ub--) {
      m = mn * loop_ub;
      if (B->data[loop_ub] != 0.0) {
        wj = b_A->data[loop_ub + m];
        B->data[loop_ub] /= wj;
        for (i = 0; i + 1 <= loop_ub; i++) {
          B->data[i] -= B->data[loop_ub] * b_A->data[i + m];
        }
      }
    }
  } else {
    m = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity((emxArray__common *)b_A, m, sizeof(double));
    loop_ub = A->size[0] * A->size[1];
    for (m = 0; m < loop_ub; m++) {
      b_A->data[m] = A->data[m];
    }

    xgeqp3(b_A, tau, jpvt);
    rankA = rankFromQR(b_A);
    m = b_B->size[0];
    b_B->size[0] = B->size[0];
    emxEnsureCapacity((emxArray__common *)b_B, m, sizeof(double));
    loop_ub = B->size[0];
    for (m = 0; m < loop_ub; m++) {
      b_B->data[m] = B->data[m];
    }

    loop_ub = b_A->size[1];
    m = B->size[0];
    B->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)B, m, sizeof(double));
    for (m = 0; m < loop_ub; m++) {
      B->data[m] = 0.0;
    }

    m = b_A->size[0];
    loop_ub = b_A->size[0];
    mn = b_A->size[1];
    if (loop_ub < mn) {
      mn = loop_ub;
    }

    for (loop_ub = 0; loop_ub + 1 <= mn; loop_ub++) {
      if (tau->data[loop_ub] != 0.0) {
        wj = b_B->data[loop_ub];
        for (i = loop_ub + 1; i + 1 <= m; i++) {
          wj += b_A->data[i + b_A->size[0] * loop_ub] * b_B->data[i];
        }

        wj *= tau->data[loop_ub];
        if (wj != 0.0) {
          b_B->data[loop_ub] -= wj;
          for (i = loop_ub + 1; i + 1 <= m; i++) {
            b_B->data[i] -= b_A->data[i + b_A->size[0] * loop_ub] * wj;
          }
        }
      }
    }

    for (i = 0; i + 1 <= rankA; i++) {
      B->data[jpvt->data[i] - 1] = b_B->data[i];
    }

    for (loop_ub = rankA - 1; loop_ub + 1 > 0; loop_ub--) {
      B->data[jpvt->data[loop_ub] - 1] /= b_A->data[loop_ub + b_A->size[0] *
        loop_ub];
      for (i = 0; i + 1 <= loop_ub; i++) {
        B->data[jpvt->data[i] - 1] -= B->data[jpvt->data[loop_ub] - 1] *
          b_A->data[i + b_A->size[0] * loop_ub];
      }
    }
  }

  emxFree_real_T(&b_B);
  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
}

void mldivide(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *Y)
{
  emxArray_real_T *b_A;
  unsigned int unnamed_idx_0;
  emxArray_int32_T *jpvt;
  int mn;
  int n;
  int jBcol;
  emxArray_real_T *tau;
  emxArray_real_T *b_B;
  int nb;
  int m;
  int j;
  int k;
  double wj;
  int i;
  double y;
  if (B->size[1] == 0) {
    unnamed_idx_0 = (unsigned int)A->size[1];
    mn = Y->size[0] * Y->size[1];
    Y->size[0] = (int)unnamed_idx_0;
    Y->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)Y, mn, sizeof(double));
  } else {
    emxInit_real_T1(&b_A, 2);
    emxInit_int32_T1(&jpvt, 2);
    if (A->size[0] == A->size[1]) {
      n = A->size[1];
      mn = b_A->size[0] * b_A->size[1];
      b_A->size[0] = A->size[0];
      b_A->size[1] = A->size[1];
      emxEnsureCapacity((emxArray__common *)b_A, mn, sizeof(double));
      jBcol = A->size[0] * A->size[1];
      for (mn = 0; mn < jBcol; mn++) {
        b_A->data[mn] = A->data[mn];
      }

      xgetrf(A->size[1], A->size[1], b_A, A->size[1], jpvt, &jBcol);
      nb = B->size[1];
      mn = Y->size[0] * Y->size[1];
      Y->size[0] = B->size[0];
      Y->size[1] = B->size[1];
      emxEnsureCapacity((emxArray__common *)Y, mn, sizeof(double));
      jBcol = B->size[0] * B->size[1];
      for (mn = 0; mn < jBcol; mn++) {
        Y->data[mn] = B->data[mn];
      }

      for (m = 0; m + 1 < n; m++) {
        if (jpvt->data[m] != m + 1) {
          jBcol = jpvt->data[m] - 1;
          for (mn = 0; mn + 1 <= nb; mn++) {
            wj = Y->data[m + Y->size[0] * mn];
            Y->data[m + Y->size[0] * mn] = Y->data[jBcol + Y->size[0] * mn];
            Y->data[jBcol + Y->size[0] * mn] = wj;
          }
        }
      }

      for (j = 1; j <= nb; j++) {
        jBcol = n * (j - 1);
        for (k = 0; k + 1 <= n; k++) {
          m = n * k;
          if (Y->data[k + jBcol] != 0.0) {
            for (i = k + 1; i + 1 <= n; i++) {
              Y->data[i + jBcol] -= Y->data[k + jBcol] * b_A->data[i + m];
            }
          }
        }
      }

      for (j = 1; j <= nb; j++) {
        jBcol = n * (j - 1);
        for (k = n - 1; k + 1 > 0; k--) {
          m = n * k;
          if (Y->data[k + jBcol] != 0.0) {
            wj = Y->data[k + jBcol];
            y = b_A->data[k + m];
            Y->data[k + jBcol] = wj / y;
            for (i = 0; i + 1 <= k; i++) {
              Y->data[i + jBcol] -= Y->data[k + jBcol] * b_A->data[i + m];
            }
          }
        }
      }
    } else {
      mn = b_A->size[0] * b_A->size[1];
      b_A->size[0] = A->size[0];
      b_A->size[1] = A->size[1];
      emxEnsureCapacity((emxArray__common *)b_A, mn, sizeof(double));
      jBcol = A->size[0] * A->size[1];
      for (mn = 0; mn < jBcol; mn++) {
        b_A->data[mn] = A->data[mn];
      }

      emxInit_real_T(&tau, 1);
      emxInit_real_T1(&b_B, 2);
      xgeqp3(b_A, tau, jpvt);
      nb = rankFromQR(b_A);
      mn = b_B->size[0] * b_B->size[1];
      b_B->size[0] = B->size[0];
      b_B->size[1] = B->size[1];
      emxEnsureCapacity((emxArray__common *)b_B, mn, sizeof(double));
      jBcol = B->size[0] * B->size[1];
      for (mn = 0; mn < jBcol; mn++) {
        b_B->data[mn] = B->data[mn];
      }

      jBcol = b_A->size[1];
      m = B->size[1];
      mn = Y->size[0] * Y->size[1];
      Y->size[0] = jBcol;
      Y->size[1] = m;
      emxEnsureCapacity((emxArray__common *)Y, mn, sizeof(double));
      jBcol *= m;
      for (mn = 0; mn < jBcol; mn++) {
        Y->data[mn] = 0.0;
      }

      m = b_A->size[0];
      jBcol = b_A->size[0];
      mn = b_A->size[1];
      if (jBcol < mn) {
        mn = jBcol;
      }

      for (j = 0; j + 1 <= mn; j++) {
        if (tau->data[j] != 0.0) {
          for (k = 0; k + 1 <= B->size[1]; k++) {
            wj = b_B->data[j + b_B->size[0] * k];
            for (i = j + 1; i + 1 <= m; i++) {
              wj += b_A->data[i + b_A->size[0] * j] * b_B->data[i + b_B->size[0]
                * k];
            }

            wj *= tau->data[j];
            if (wj != 0.0) {
              b_B->data[j + b_B->size[0] * k] -= wj;
              for (i = j + 1; i + 1 <= m; i++) {
                b_B->data[i + b_B->size[0] * k] -= b_A->data[i + b_A->size[0] *
                  j] * wj;
              }
            }
          }
        }
      }

      emxFree_real_T(&tau);
      for (k = 0; k + 1 <= B->size[1]; k++) {
        for (i = 0; i + 1 <= nb; i++) {
          Y->data[(jpvt->data[i] + Y->size[0] * k) - 1] = b_B->data[i +
            b_B->size[0] * k];
        }

        for (j = nb - 1; j + 1 > 0; j--) {
          Y->data[(jpvt->data[j] + Y->size[0] * k) - 1] /= b_A->data[j +
            b_A->size[0] * j];
          for (i = 0; i + 1 <= j; i++) {
            Y->data[(jpvt->data[i] + Y->size[0] * k) - 1] -= Y->data[(jpvt->
              data[j] + Y->size[0] * k) - 1] * b_A->data[i + b_A->size[0] * j];
          }
        }
      }

      emxFree_real_T(&b_B);
    }

    emxFree_int32_T(&jpvt);
    emxFree_real_T(&b_A);
  }
}

/* End of code generation (mldivide.cpp) */
