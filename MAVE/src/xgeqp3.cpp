/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xgeqp3.cpp
 *
 * Code generation for function 'xgeqp3'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xgeqp3.h"
#include "xnrm2.h"
#include "xzlarfg.h"
#include "MAVEfast_emxutil.h"
#include "colon.h"

/* Function Definitions */
void xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, emxArray_int32_T *jpvt)
{
  int m;
  int n;
  int k;
  int mn;
  int i23;
  emxArray_real_T *work;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int nmi;
  int i;
  int i_i;
  int mmi;
  int pvt;
  int ix;
  double smax;
  int iy;
  double s;
  int i_ip1;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int exitg1;
  double absxk;
  double t;
  m = A->size[0];
  n = A->size[1];
  k = A->size[0];
  mn = A->size[1];
  if (k < mn) {
    mn = k;
  }

  i23 = tau->size[0];
  tau->size[0] = mn;
  emxEnsureCapacity((emxArray__common *)tau, i23, sizeof(double));
  eml_signed_integer_colon(A->size[1], jpvt);
  if (!((A->size[0] == 0) || (A->size[1] == 0))) {
    emxInit_real_T(&work, 1);
    k = A->size[1];
    i23 = work->size[0];
    work->size[0] = k;
    emxEnsureCapacity((emxArray__common *)work, i23, sizeof(double));
    for (i23 = 0; i23 < k; i23++) {
      work->data[i23] = 0.0;
    }

    emxInit_real_T(&vn1, 1);
    emxInit_real_T(&vn2, 1);
    k = A->size[1];
    i23 = vn1->size[0];
    vn1->size[0] = k;
    emxEnsureCapacity((emxArray__common *)vn1, i23, sizeof(double));
    i23 = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity((emxArray__common *)vn2, i23, sizeof(double));
    k = 1;
    for (nmi = 0; nmi + 1 <= n; nmi++) {
      vn1->data[nmi] = xnrm2(m, A, k);
      vn2->data[nmi] = vn1->data[nmi];
      k += m;
    }

    for (i = 0; i + 1 <= mn; i++) {
      i_i = i + i * m;
      nmi = n - i;
      mmi = (m - i) - 1;
      if (nmi < 1) {
        pvt = -1;
      } else {
        pvt = 0;
        if (nmi > 1) {
          ix = i;
          smax = std::abs(vn1->data[i]);
          for (k = 0; k + 2 <= nmi; k++) {
            ix++;
            s = std::abs(vn1->data[ix]);
            if (s > smax) {
              pvt = k + 1;
              smax = s;
            }
          }
        }
      }

      pvt += i;
      if (pvt + 1 != i + 1) {
        ix = m * pvt;
        iy = m * i;
        for (k = 1; k <= m; k++) {
          smax = A->data[ix];
          A->data[ix] = A->data[iy];
          A->data[iy] = smax;
          ix++;
          iy++;
        }

        k = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i];
        jpvt->data[i] = k;
        vn1->data[pvt] = vn1->data[i];
        vn2->data[pvt] = vn2->data[i];
      }

      if (i + 1 < m) {
        s = A->data[i_i];
        tau->data[i] = xzlarfg(1 + mmi, &s, A, i_i + 2);
        A->data[i_i] = s;
      } else {
        tau->data[i] = 0.0;
      }

      if (i + 1 < n) {
        s = A->data[i_i];
        A->data[i_i] = 1.0;
        i_ip1 = (i + (i + 1) * m) + 1;
        if (tau->data[i] != 0.0) {
          lastv = mmi;
          pvt = i_i + mmi;
          while ((lastv + 1 > 0) && (A->data[pvt] == 0.0)) {
            lastv--;
            pvt--;
          }

          lastc = nmi - 1;
          exitg2 = false;
          while ((!exitg2) && (lastc > 0)) {
            k = i_ip1 + (lastc - 1) * m;
            nmi = k;
            do {
              exitg1 = 0;
              if (nmi <= k + lastv) {
                if (A->data[nmi - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  nmi++;
                }
              } else {
                lastc--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = -1;
          lastc = 0;
        }

        if (lastv + 1 > 0) {
          if (lastc != 0) {
            for (iy = 1; iy <= lastc; iy++) {
              work->data[iy - 1] = 0.0;
            }

            iy = 0;
            i23 = i_ip1 + m * (lastc - 1);
            pvt = i_ip1;
            while ((m > 0) && (pvt <= i23)) {
              ix = i_i;
              smax = 0.0;
              k = pvt + lastv;
              for (nmi = pvt; nmi <= k; nmi++) {
                smax += A->data[nmi - 1] * A->data[ix];
                ix++;
              }

              work->data[iy] += smax;
              iy++;
              pvt += m;
            }
          }

          if (!(-tau->data[i] == 0.0)) {
            pvt = 0;
            for (nmi = 1; nmi <= lastc; nmi++) {
              if (work->data[pvt] != 0.0) {
                smax = work->data[pvt] * -tau->data[i];
                ix = i_i;
                i23 = lastv + i_ip1;
                for (k = i_ip1; k <= i23; k++) {
                  A->data[k - 1] += A->data[ix] * smax;
                  ix++;
                }
              }

              pvt++;
              i_ip1 += m;
            }
          }
        }

        A->data[i_i] = s;
      }

      for (nmi = i + 1; nmi + 1 <= n; nmi++) {
        k = (i + m * nmi) + 1;
        if (vn1->data[nmi] != 0.0) {
          smax = std::abs(A->data[i + A->size[0] * nmi]) / vn1->data[nmi];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          s = vn1->data[nmi] / vn2->data[nmi];
          s = smax * (s * s);
          if (s <= 1.4901161193847656E-8) {
            if (i + 1 < m) {
              smax = 0.0;
              if (!(mmi < 1)) {
                if (mmi == 1) {
                  smax = std::abs(A->data[k]);
                } else {
                  s = 2.2250738585072014E-308;
                  pvt = k + mmi;
                  while (k + 1 <= pvt) {
                    absxk = std::abs(A->data[k]);
                    if (absxk > s) {
                      t = s / absxk;
                      smax = 1.0 + smax * t * t;
                      s = absxk;
                    } else {
                      t = absxk / s;
                      smax += t * t;
                    }

                    k++;
                  }

                  smax = s * std::sqrt(smax);
                }
              }

              vn1->data[nmi] = smax;
              vn2->data[nmi] = vn1->data[nmi];
            } else {
              vn1->data[nmi] = 0.0;
              vn2->data[nmi] = 0.0;
            }
          } else {
            vn1->data[nmi] *= std::sqrt(smax);
          }
        }
      }
    }

    emxFree_real_T(&vn2);
    emxFree_real_T(&vn1);
    emxFree_real_T(&work);
  }
}

/* End of code generation (xgeqp3.cpp) */
