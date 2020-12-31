/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * abs.h
 *
 * Code generation for function 'abs'
 *
 */

#ifndef ABS_H
#define ABS_H

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
extern void b_abs(const emxArray_boolean_T *x, emxArray_real_T *y);
extern void c_abs(const emxArray_real_T *x, emxArray_real_T *y);
extern void d_abs(const emxArray_real_T *x, emxArray_real_T *y);

#endif

/* End of code generation (abs.h) */
