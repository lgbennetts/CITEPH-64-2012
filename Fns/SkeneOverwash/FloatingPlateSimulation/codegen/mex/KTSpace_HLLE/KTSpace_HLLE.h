/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * KTSpace_HLLE.h
 *
 * Code generation for function 'KTSpace_HLLE'
 *
 */

#pragma once

/* Include files */
#include "KTSpace_HLLE_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void KTSpace_HLLE(const emlrtStack *sp, const emxArray_real_T *u, real_T dx,
                  const emxArray_real_T *x, real_T theta, real_T mindepth,
                  emxArray_real_T *spacediscrete, real_T *amax);

/* End of code generation (KTSpace_HLLE.h) */
