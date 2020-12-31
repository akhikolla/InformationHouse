/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * repmat.h
 *
 * Code generation for function 'repmat'
 *
 */

#ifndef REPMAT_H
#define REPMAT_H

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
extern void b_repmat(const emxArray_real_T *a, double varargin_2,
                     emxArray_real_T *b);
extern void repmat(const emxArray_real_T *a, double varargin_1, emxArray_real_T *
                   b);

#endif

/* End of code generation (repmat.h) */
