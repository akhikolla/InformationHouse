/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * schur.h
 *
 * Code generation for function 'schur'
 *
 */

#ifndef SCHUR_H
#define SCHUR_H

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
extern void schur(const emxArray_real_T *A, emxArray_creal_T *V,
                  emxArray_creal_T *T);

#endif

/* End of code generation (schur.h) */
