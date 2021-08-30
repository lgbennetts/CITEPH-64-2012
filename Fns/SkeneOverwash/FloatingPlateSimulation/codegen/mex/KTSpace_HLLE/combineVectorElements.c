/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * combineVectorElements.c
 *
 * Code generation for function 'combineVectorElements'
 *
 */

/* Include files */
#include "combineVectorElements.h"
#include "KTSpace_HLLE_data.h"
#include "KTSpace_HLLE_emxutil.h"
#include "KTSpace_HLLE_types.h"
#include "eml_int_forloop_overflow_check.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo rb_emlrtRSI = { 173,/* lineNo */
  "colMajorFlatIter",                  /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pathName */
};

static emlrtRTEInfo u_emlrtRTEI = { 166,/* lineNo */
  24,                                  /* colNo */
  "combineVectorElements",             /* fName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/private/combineVectorElements.m"/* pName */
};

/* Function Definitions */
void combineVectorElements(const emlrtStack *sp, const emxArray_boolean_T *x,
  emxArray_int32_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  int32_T i;
  int32_T npages;
  int32_T xpageoffset;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  if (x->size[1] == 0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    st.site = &qb_emlrtRSI;
    npages = x->size[1];
    xpageoffset = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_int32_T(&st, y, xpageoffset, &u_emlrtRTEI);
    b_st.site = &rb_emlrtRSI;
    if (x->size[1] > 2147483646) {
      c_st.site = &bb_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (i = 0; i < npages; i++) {
      xpageoffset = i << 1;
      y->data[i] = x->data[xpageoffset] + x->data[xpageoffset + 1];
    }
  }
}

/* End of code generation (combineVectorElements.c) */
