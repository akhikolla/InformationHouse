/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xdlanv2.cpp
 *
 * Code generation for function 'xdlanv2'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MAVEfast.h"
#include "xdlanv2.h"
#include "xzlarfg.h"
#include "MAVEfast_rtwutil.h"

/* Function Definitions */
void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r, double
             *rt1i, double *rt2r, double *rt2i, double *cs, double *sn)
{
  double temp;
  double p;
  double z;
  double scale;
  double bcmax;
  double b_z;
  int b_b;
  int b_c;
  double bcmis;
  double tau;
  double b_p;
  int b_scale;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    temp = *d;
    *d = *a;
    *a = temp;
    *b = -*c;
    *c = 0.0;
  } else if ((*a - *d == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
    *cs = 1.0;
    *sn = 0.0;
  } else {
    temp = *a - *d;
    p = 0.5 * temp;
    z = std::abs(*b);
    scale = std::abs(*c);
    if ((z > scale) || rtIsNaN(scale)) {
      bcmax = z;
    } else {
      bcmax = scale;
    }

    z = std::abs(*b);
    scale = std::abs(*c);
    if ((z < scale) || rtIsNaN(scale)) {
      b_z = z;
    } else {
      b_z = scale;
    }

    if (!(*b < 0.0)) {
      b_b = 1;
    } else {
      b_b = -1;
    }

    if (!(*c < 0.0)) {
      b_c = 1;
    } else {
      b_c = -1;
    }

    bcmis = b_z * (double)b_b * (double)b_c;
    z = std::abs(p);
    if ((z > bcmax) || rtIsNaN(bcmax)) {
      scale = z;
    } else {
      scale = bcmax;
    }

    z = p / scale * p + bcmax / scale * bcmis;
    if (z >= 8.8817841970012523E-16) {
      *a = std::sqrt(scale) * std::sqrt(z);
      if (!(p < 0.0)) {
        b_p = *a;
      } else {
        b_p = -*a;
      }

      z = p + b_p;
      *a = *d + z;
      *d -= bcmax / z * bcmis;
      tau = rt_hypotd_snf(*c, z);
      *cs = z / tau;
      *sn = *c / tau;
      *b -= *c;
      *c = 0.0;
    } else {
      scale = *b + *c;
      tau = rt_hypotd_snf(scale, temp);
      *cs = std::sqrt(0.5 * (1.0 + std::abs(scale) / tau));
      if (!(scale < 0.0)) {
        b_scale = 1;
      } else {
        b_scale = -1;
      }

      *sn = -(p / (tau * *cs)) * (double)b_scale;
      bcmax = *a * *cs + *b * *sn;
      z = -*a * *sn + *b * *cs;
      bcmis = *c * *cs + *d * *sn;
      scale = -*c * *sn + *d * *cs;
      *b = z * *cs + scale * *sn;
      *c = -bcmax * *sn + bcmis * *cs;
      temp = 0.5 * ((bcmax * *cs + bcmis * *sn) + (-z * *sn + scale * *cs));
      *a = temp;
      *d = temp;
      if (*c != 0.0) {
        if (*b != 0.0) {
          if ((*b < 0.0) == (*c < 0.0)) {
            scale = std::sqrt(std::abs(*b));
            bcmis = std::sqrt(std::abs(*c));
            *a = scale * bcmis;
            if (!(*c < 0.0)) {
              p = *a;
            } else {
              p = -*a;
            }

            tau = 1.0 / std::sqrt(std::abs(*b + *c));
            *a = temp + p;
            *d = temp - p;
            *b -= *c;
            *c = 0.0;
            bcmax = scale * tau;
            scale = bcmis * tau;
            temp = *cs * bcmax - *sn * scale;
            *sn = *cs * scale + *sn * bcmax;
            *cs = temp;
          }
        } else {
          *b = -*c;
          *c = 0.0;
          temp = *cs;
          *cs = -*sn;
          *sn = temp;
        }
      }
    }
  }

  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = std::sqrt(std::abs(*b)) * std::sqrt(std::abs(*c));
    *rt2i = -*rt1i;
  }
}

/* End of code generation (xdlanv2.cpp) */
