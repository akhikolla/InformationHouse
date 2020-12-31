/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * mldivide.h
 *
 * Code generation for function 'mldivide'
 *
 */

#ifndef MLDIVIDE_H
#define MLDIVIDE_H

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
extern void b_mldivide(const emxArray_real_T *A, emxArray_real_T *B);
extern void mldivide(const emxArray_real_T *A, const emxArray_real_T *B,
                     emxArray_real_T *Y);

#endif

/* End of code generation (mldivide.h) */
