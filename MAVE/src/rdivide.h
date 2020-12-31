/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * rdivide.h
 *
 * Code generation for function 'rdivide'
 *
 */

#ifndef RDIVIDE_H
#define RDIVIDE_H

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
extern void b_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z);
extern void rdivide(const emxArray_real_T *y, emxArray_real_T *z);

#endif

/* End of code generation (rdivide.h) */
