/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * minmod.c
 *
 * Code generation for function 'minmod'
 *
 */

/* Include files */
#include "minmod.h"
#include "KTSpace_HLLE_data.h"
#include "KTSpace_HLLE_emxutil.h"
#include "KTSpace_HLLE_types.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtBCInfo jd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  14,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo kd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  19,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ld_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  8,                                   /* colNo */
  "ux",                                /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo md_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  7,                                   /* lineNo */
  16,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo nd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  7,                                   /* lineNo */
  21,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo od_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  8,                                   /* lineNo */
  16,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo pd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  8,                                   /* lineNo */
  23,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo qd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  9,                                   /* lineNo */
  16,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo rd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  9,                                   /* lineNo */
  23,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo sd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  22,                                  /* lineNo */
  12,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo td_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  22,                                  /* lineNo */
  17,                                  /* colNo */
  "u",                                 /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ud_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  22,                                  /* lineNo */
  5,                                   /* colNo */
  "ux",                                /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo vd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  12,                                  /* lineNo */
  13,                                  /* colNo */
  "ux",                                /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo wd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  16,                                  /* lineNo */
  13,                                  /* colNo */
  "ux",                                /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo xd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  14,                                  /* lineNo */
  13,                                  /* colNo */
  "ux",                                /* aName */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo s_emlrtRTEI = { 3, /* lineNo */
  5,                                   /* colNo */
  "minmod",                            /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/minmod.m"/* pName */
};

/* Function Definitions */
void minmod(const emlrtStack *sp, const emxArray_real_T *u, real_T theta,
            emxArray_real_T *ux)
{
  real_T term1;
  real_T term2;
  real_T term3;
  int32_T i;
  int32_T loop_ub;
  i = ux->size[0] * ux->size[1];
  ux->size[0] = 1;
  ux->size[1] = u->size[1];
  emxEnsureCapacity_real_T(sp, ux, i, &s_emlrtRTEI);
  loop_ub = u->size[0] * u->size[1];
  for (i = 0; i < loop_ub; i++) {
    ux->data[i] = u->data[i] * 0.0;
  }

  i = u->size[1];
  for (loop_ub = 0; loop_ub <= i - 3; loop_ub++) {
    if ((loop_ub + 2 < 1) || (loop_ub + 2 > u->size[1])) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 2, 1, u->size[1], &md_emlrtBCI, sp);
    }

    if ((loop_ub + 1 < 1) || (loop_ub + 1 > u->size[1])) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1, u->size[1], &nd_emlrtBCI, sp);
    }

    term1 = (u->data[loop_ub + 1] - u->data[loop_ub]) * theta;
    if ((loop_ub + 3 < 1) || (loop_ub + 3 > u->size[1])) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 3, 1, u->size[1], &od_emlrtBCI, sp);
    }

    if ((loop_ub + 1 < 1) || (loop_ub + 1 > u->size[1])) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1, u->size[1], &pd_emlrtBCI, sp);
    }

    term2 = (u->data[loop_ub + 2] - u->data[loop_ub]) / 2.0;
    if ((loop_ub + 3 < 1) || (loop_ub + 3 > u->size[1])) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 3, 1, u->size[1], &qd_emlrtBCI, sp);
    }

    if ((loop_ub + 2 < 1) || (loop_ub + 2 > u->size[1])) {
      emlrtDynamicBoundsCheckR2012b(loop_ub + 2, 1, u->size[1], &rd_emlrtBCI, sp);
    }

    term3 = (u->data[loop_ub + 2] - u->data[loop_ub + 1]) * theta;
    if ((term1 > 0.0) && (term2 > 0.0) && (term3 > 0.0)) {
      if ((loop_ub + 2 < 1) || (loop_ub + 2 > ux->size[1])) {
        emlrtDynamicBoundsCheckR2012b(loop_ub + 2, 1, ux->size[1], &vd_emlrtBCI,
          sp);
      }

      ux->data[loop_ub + 1] = muDoubleScalarMin(term1, muDoubleScalarMin(term2,
        term3));
    } else if ((term1 < 0.0) && (term2 < 0.0) && (term3 < 0.0)) {
      if ((loop_ub + 2 < 1) || (loop_ub + 2 > ux->size[1])) {
        emlrtDynamicBoundsCheckR2012b(loop_ub + 2, 1, ux->size[1], &xd_emlrtBCI,
          sp);
      }

      ux->data[loop_ub + 1] = muDoubleScalarMax(term1, muDoubleScalarMax(term2,
        term3));
    } else {
      if ((loop_ub + 2 < 1) || (loop_ub + 2 > ux->size[1])) {
        emlrtDynamicBoundsCheckR2012b(loop_ub + 2, 1, ux->size[1], &wd_emlrtBCI,
          sp);
      }

      ux->data[loop_ub + 1] = 0.0;
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  if ((ux->size[1] < 1) || (ux->size[1] > u->size[1])) {
    emlrtDynamicBoundsCheckR2012b(ux->size[1], 1, u->size[1], &sd_emlrtBCI, sp);
  }

  if ((ux->size[1] - 1 < 1) || (ux->size[1] - 1 > u->size[1])) {
    emlrtDynamicBoundsCheckR2012b(ux->size[1] - 1, 1, u->size[1], &td_emlrtBCI,
      sp);
  }

  if (ux->size[1] < 1) {
    emlrtDynamicBoundsCheckR2012b(ux->size[1], 1, ux->size[1], &ud_emlrtBCI, sp);
  }

  ux->data[ux->size[1] - 1] = (u->data[ux->size[1] - 1] - u->data[ux->size[1] -
    2]) * theta;
  if (1 > ux->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, ux->size[1], &ld_emlrtBCI, sp);
  }

  if (2 > u->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, u->size[1], &jd_emlrtBCI, sp);
  }

  if (1 > u->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, u->size[1], &kd_emlrtBCI, sp);
  }

  ux->data[0] = (u->data[1] - u->data[0]) * theta;
}

/* End of code generation (minmod.c) */
