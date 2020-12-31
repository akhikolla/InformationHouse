/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * diag.h
 *
 * Code generation for function 'diag'
 *
 */

#ifndef DIAG_H
#define DIAG_H

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
extern void b_diag(const emxArray_real_T *v, emxArray_real_T *d);
extern void c_diag(const emxArray_creal_T *v, emxArray_creal_T *d);
extern void diag(const emxArray_real_T *v, emxArray_real_T *d);

#endif

/* End of code generation (diag.h) */
