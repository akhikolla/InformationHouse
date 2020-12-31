/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sortIdx.h
 *
 * Code generation for function 'sortIdx'
 *
 */

#ifndef SORTIDX_H
#define SORTIDX_H

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
extern void b_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
extern void c_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
extern void sortIdx(const emxArray_real_T *x, emxArray_int32_T *idx);

#endif

/* End of code generation (sortIdx.h) */
