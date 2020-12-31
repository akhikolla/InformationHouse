/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * mean.h
 *
 * Code generation for function 'mean'
 *
 */

#ifndef MEAN_H
#define MEAN_H

/* Include files */
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "MAVEfast_types.h"

/* Function Declarations */
extern void b_mean(const emxArray_real_T *x, emxArray_real_T *y);
extern double c_mean(const emxArray_real_T *x);
extern double d_mean(const emxArray_real_T *x);
extern void mean(const emxArray_real_T *x, emxArray_real_T *y);

#endif

/* End of code generation (mean.h) */
