/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzhgeqz.h
 *
 * Code generation for function 'xzhgeqz'
 *
 */

#ifndef XZHGEQZ_H
#define XZHGEQZ_H

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
extern void xzhgeqz(emxArray_creal_T *A, int ilo, int ihi, emxArray_creal_T *Z,
                    int *info, emxArray_creal_T *alpha1, emxArray_creal_T *beta1);

#endif

/* End of code generation (xzhgeqz.h) */
