/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * KTSpace_HLLE_initialize.c
 *
 * Code generation for function 'KTSpace_HLLE_initialize'
 *
 */

/* Include files */
#include "KTSpace_HLLE_initialize.h"
#include "KTSpace_HLLE_data.h"
#include "_coder_KTSpace_HLLE_mex.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void KTSpace_HLLE_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (KTSpace_HLLE_initialize.c) */
