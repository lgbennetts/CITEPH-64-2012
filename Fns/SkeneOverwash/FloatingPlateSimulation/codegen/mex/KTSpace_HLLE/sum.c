/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "KTSpace_HLLE_data.h"
#include "KTSpace_HLLE_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo sb_emlrtRSI = { 20, /* lineNo */
  "sum",                               /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

static emlrtRSInfo tb_emlrtRSI = { 194,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

/* Function Definitions */
real_T sum(const emlrtStack *sp, const emxArray_real_T *x)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  real_T y;
  int32_T k;
  int32_T vlen;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &sb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  b_st.site = &pb_emlrtRSI;
  vlen = x->size[1];
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    c_st.site = &qb_emlrtRSI;
    y = x->data[0];
    d_st.site = &tb_emlrtRSI;
    if ((2 <= x->size[1]) && (x->size[1] > 2147483646)) {
      e_st.site = &bb_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    for (k = 2; k <= vlen; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

/* End of code generation (sum.c) */
