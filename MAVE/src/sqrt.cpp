/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sqrt.cpp
 *
 * Code generation for function 'sqrt'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "sqrt.h"
#include "xzlarfg.h"
#include "MAVEfast_rtwutil.h"

/* Function Definitions */
void b_sqrt(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[0] * x->size[1];
  for (k = 0; k + 1 <= nx; k++) {
    x->data[k] = std::sqrt(x->data[k]);
  }
}

void c_sqrt(creal_T *x)
{
  double xr;
  double xi;
  double absxi;
  double absxr;
  xr = x->re;
  xi = x->im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      absxi = 0.0;
      xr = std::sqrt(-xr);
    } else {
      absxi = std::sqrt(xr);
      xr = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      absxi = std::sqrt(-xi / 2.0);
      xr = -absxi;
    } else {
      absxi = std::sqrt(xi / 2.0);
      xr = absxi;
    }
  } else if (rtIsNaN(xr)) {
    absxi = xr;
  } else if (rtIsNaN(xi)) {
    absxi = xi;
    xr = xi;
  } else if (rtIsInf(xi)) {
    absxi = std::abs(xi);
    xr = xi;
  } else if (rtIsInf(xr)) {
    if (xr < 0.0) {
      absxi = 0.0;
      xr = xi * -xr;
    } else {
      absxi = xr;
      xr = 0.0;
    }
  } else {
    absxr = std::abs(xr);
    absxi = std::abs(xi);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi);
      if (absxi > absxr) {
        absxi = std::sqrt(absxi) * std::sqrt(1.0 + absxr / absxi);
      } else {
        absxi = std::sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxi = std::sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }

    if (xr > 0.0) {
      xr = 0.5 * (xi / absxi);
    } else {
      if (xi < 0.0) {
        xr = -absxi;
      } else {
        xr = absxi;
      }

      absxi = 0.5 * (xi / xr);
    }
  }

  x->re = absxi;
  x->im = xr;
}

void d_sqrt(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[0];
  for (k = 0; k + 1 <= nx; k++) {
    x->data[k] = std::sqrt(x->data[k]);
  }
}

/* End of code generation (sqrt.cpp) */
