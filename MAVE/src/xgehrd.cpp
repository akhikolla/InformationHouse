/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xgehrd.cpp
 *
 * Code generation for function 'xgehrd'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xgehrd.h"
#include "xzlarf.h"
#include "xzlarfg.h"
#include "MAVEfast_emxutil.h"

/* Function Definitions */
void xgehrd(emxArray_real_T *a, emxArray_real_T *tau)
{
  int n;
  int ntau;
  emxArray_real_T *work;
  int i16;
  int i;
  int im1n;
  int in;
  double alpha1;
  int jy;
  int lastv;
  int lastc;
  boolean_T exitg2;
  int iac;
  int ix;
  int ia;
  int exitg1;
  double temp;
  int ijA;
  n = a->size[0];
  if (a->size[0] < 1) {
    ntau = 0;
  } else {
    ntau = a->size[0] - 1;
  }

  emxInit_real_T(&work, 1);
  i16 = tau->size[0];
  tau->size[0] = ntau;
  emxEnsureCapacity((emxArray__common *)tau, i16, sizeof(double));
  ntau = a->size[0];
  i16 = work->size[0];
  work->size[0] = ntau;
  emxEnsureCapacity((emxArray__common *)work, i16, sizeof(double));
  for (i16 = 0; i16 < ntau; i16++) {
    work->data[i16] = 0.0;
  }

  for (i = 0; i + 1 < n; i++) {
    im1n = i * n + 2;
    in = (i + 1) * n;
    alpha1 = a->data[(i + a->size[0] * i) + 1];
    ntau = i + 3;
    if (!(ntau < n)) {
      ntau = n;
    }

    tau->data[i] = xzlarfg((n - i) - 1, &alpha1, a, ntau + i * n);
    a->data[(i + a->size[0] * i) + 1] = 1.0;
    ntau = (n - i) - 3;
    jy = (i + im1n) - 1;
    if (tau->data[i] != 0.0) {
      lastv = ntau + 2;
      ntau += jy;
      while ((lastv > 0) && (a->data[ntau + 1] == 0.0)) {
        lastv--;
        ntau--;
      }

      lastc = n;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        ntau = in + lastc;
        ia = ntau;
        do {
          exitg1 = 0;
          if ((n > 0) && (ia <= ntau + (lastv - 1) * n)) {
            if (a->data[ia - 1] != 0.0) {
              exitg1 = 1;
            } else {
              ia += n;
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
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      if (lastc != 0) {
        for (ntau = 1; ntau <= lastc; ntau++) {
          work->data[ntau - 1] = 0.0;
        }

        ix = jy;
        i16 = (in + n * (lastv - 1)) + 1;
        iac = in + 1;
        while ((n > 0) && (iac <= i16)) {
          ntau = 0;
          ijA = (iac + lastc) - 1;
          for (ia = iac; ia <= ijA; ia++) {
            work->data[ntau] += a->data[ia - 1] * a->data[ix];
            ntau++;
          }

          ix++;
          iac += n;
        }
      }

      if (!(-tau->data[i] == 0.0)) {
        ntau = in;
        for (iac = 1; iac <= lastv; iac++) {
          if (a->data[jy] != 0.0) {
            temp = a->data[jy] * -tau->data[i];
            ix = 0;
            i16 = lastc + ntau;
            for (ijA = ntau; ijA + 1 <= i16; ijA++) {
              a->data[ijA] += work->data[ix] * temp;
              ix++;
            }
          }

          jy++;
          ntau += n;
        }
      }
    }

    xzlarf((n - i) - 1, (n - i) - 1, i + im1n, tau->data[i], a, (i + in) + 2, n,
           work);
    a->data[(i + a->size[0] * i) + 1] = alpha1;
  }

  emxFree_real_T(&work);
}

/* End of code generation (xgehrd.cpp) */
