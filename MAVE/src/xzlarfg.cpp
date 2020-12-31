/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzlarfg.cpp
 *
 * Code generation for function 'xzlarfg'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xzlarfg.h"
#include "xscal.h"
#include "xnrm2.h"
#include "MAVEfast_rtwutil.h"

/* Function Definitions */
double b_xzlarfg(int n, double *alpha1, double x[3])
{
  double tau;
  double xnorm;
  int knt;
  int k;
  tau = 0.0;
  if (!(n <= 0)) {
    xnorm = b_xnrm2(n - 1, x);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      if (std::abs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          for (k = 1; k + 1 <= n; k++) {
            x[k] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(xnorm) >= 1.0020841800044864E-292));

        xnorm = rt_hypotd_snf(*alpha1, b_xnrm2(n - 1, x));
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        tau = (xnorm - *alpha1) / xnorm;
        *alpha1 = 1.0 / (*alpha1 - xnorm);
        for (k = 1; k + 1 <= n; k++) {
          x[k] *= *alpha1;
        }

        for (k = 1; k <= knt; k++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        tau = (xnorm - *alpha1) / xnorm;
        *alpha1 = 1.0 / (*alpha1 - xnorm);
        for (k = 1; k + 1 <= n; k++) {
          x[k] *= *alpha1;
        }

        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

double xzlarfg(int n, double *alpha1, emxArray_real_T *x, int ix0)
{
  double tau;
  double xnorm;
  int knt;
  int k;
  tau = 0.0;
  if (!(n <= 0)) {
    xnorm = xnrm2(n - 1, x, ix0);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      if (std::abs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          knt++;
          xscal(n - 1, 9.9792015476736E+291, x, ix0);
          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(std::abs(xnorm) >= 1.0020841800044864E-292));

        xnorm = xnrm2(n - 1, x, ix0);
        xnorm = rt_hypotd_snf(*alpha1, xnorm);
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        tau = (xnorm - *alpha1) / xnorm;
        xscal(n - 1, 1.0 / (*alpha1 - xnorm), x, ix0);
        for (k = 1; k <= knt; k++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        tau = (xnorm - *alpha1) / xnorm;
        xscal(n - 1, 1.0 / (*alpha1 - xnorm), x, ix0);
        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

/* End of code generation (xzlarfg.cpp) */
