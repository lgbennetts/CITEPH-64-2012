/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * minmod.h
 *
 * Code generation for function 'minmod'
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
void minmod(const emlrtStack *sp, const emxArray_real_T *u, real_T theta,
            emxArray_real_T *ux);

/* End of code generation (minmod.h) */
