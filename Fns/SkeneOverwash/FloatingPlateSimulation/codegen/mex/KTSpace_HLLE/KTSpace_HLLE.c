/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * KTSpace_HLLE.c
 *
 * Code generation for function 'KTSpace_HLLE'
 *
 */

/* Include files */
#include "KTSpace_HLLE.h"
#include "KTSpace_HLLE_data.h"
#include "KTSpace_HLLE_emxutil.h"
#include "KTSpace_HLLE_types.h"
#include "abs.h"
#include "combineVectorElements.h"
#include "eml_int_forloop_overflow_check.h"
#include "f.h"
#include "minmod.h"
#include "rt_nonfinite.h"
#include "spectralradius.h"
#include "sum.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 13,    /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 34,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 35,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 36,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 37,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 40,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 48,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 50,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 67,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 85,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 87,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 91,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 98,  /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 104, /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 107, /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 140, /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 145, /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 146, /* lineNo */
  "KTSpace_HLLE",                      /* fcnName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 14,  /* lineNo */
  "min",                               /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/min.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 46,  /* lineNo */
  "minOrMax",                          /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 92,  /* lineNo */
  "minimum",                           /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 157, /* lineNo */
  "unaryMinOrMax",                     /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 1017,/* lineNo */
  "minRealVectorOmitNaN",              /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 934, /* lineNo */
  "minOrMaxRealVector",                /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 994,/* lineNo */
  "minOrMaxRealVectorKernel",          /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo cb_emlrtRSI = { 926,/* lineNo */
  "minOrMaxRealVector",                /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo db_emlrtRSI = { 977,/* lineNo */
  "findFirst",                         /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo eb_emlrtRSI = { 14, /* lineNo */
  "max",                               /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/max.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 44, /* lineNo */
  "minOrMax",                          /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo gb_emlrtRSI = { 79, /* lineNo */
  "maximum",                           /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/minOrMax.m"/* pathName */
};

static emlrtRSInfo hb_emlrtRSI = { 147,/* lineNo */
  "unaryMinOrMax",                     /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 1021,/* lineNo */
  "maxRealVectorOmitNaN",              /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pathName */
};

static emlrtRSInfo ob_emlrtRSI = { 15, /* lineNo */
  "sum",                               /* fcnName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/sum.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 97,  /* lineNo */
  27,                                  /* colNo */
  "unaryMinOrMax",                     /* fName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  143,                                 /* lineNo */
  29,                                  /* colNo */
  "spacediscrete",                     /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  77,                                  /* lineNo */
  21,                                  /* colNo */
  "spacediscrete",                     /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  71,                                  /* lineNo */
  25,                                  /* colNo */
  "spacediscrete",                     /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  59,                                  /* lineNo */
  21,                                  /* colNo */
  "spacediscrete",                     /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  140,                                 /* lineNo */
  38,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  140,                                 /* lineNo */
  26,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  68,                                  /* lineNo */
  75,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  68,                                  /* lineNo */
  65,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  67,                                  /* lineNo */
  60,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  67,                                  /* lineNo */
  47,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo k_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  67,                                  /* lineNo */
  38,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  67,                                  /* lineNo */
  26,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo m_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  56,                                  /* lineNo */
  22,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  56,                                  /* lineNo */
  14,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo o_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  53,                                  /* lineNo */
  24,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo p_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  53,                                  /* lineNo */
  14,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo q_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  51,                                  /* lineNo */
  85,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo r_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  51,                                  /* lineNo */
  73,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo s_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  51,                                  /* lineNo */
  58,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo t_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  51,                                  /* lineNo */
  48,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo u_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  51,                                  /* lineNo */
  35,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo v_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  51,                                  /* lineNo */
  25,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo w_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  50,                                  /* lineNo */
  81,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo x_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  50,                                  /* lineNo */
  71,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo y_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  50,                                  /* lineNo */
  57,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ab_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  50,                                  /* lineNo */
  42,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo bb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  50,                                  /* lineNo */
  31,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo cb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  50,                                  /* lineNo */
  17,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo db_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  49,                                  /* lineNo */
  75,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo eb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  49,                                  /* lineNo */
  65,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo fb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  49,                                  /* lineNo */
  52,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo gb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  49,                                  /* lineNo */
  44,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo hb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  49,                                  /* lineNo */
  33,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ib_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  49,                                  /* lineNo */
  25,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo jb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  48,                                  /* lineNo */
  72,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo kb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  48,                                  /* lineNo */
  64,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo lb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  48,                                  /* lineNo */
  52,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo mb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  48,                                  /* lineNo */
  39,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo nb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  48,                                  /* lineNo */
  30,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ob_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  48,                                  /* lineNo */
  18,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo pb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  35,                                  /* lineNo */
  33,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo qb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  34,                                  /* lineNo */
  34,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo rb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  27,                                  /* lineNo */
  13,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo sb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  25,                                  /* lineNo */
  13,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo tb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  17,                                  /* colNo */
  "uplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ub_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  28,                                  /* lineNo */
  12,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo vb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  12,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo wb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  21,                                  /* lineNo */
  16,                                  /* colNo */
  "uneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo xb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  28,                                  /* lineNo */
  19,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo yb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  27,                                  /* lineNo */
  32,                                  /* colNo */
  "ux",                                /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ac_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  27,                                  /* lineNo */
  20,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo bc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  25,                                  /* lineNo */
  20,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo cc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  33,                                  /* colNo */
  "ux",                                /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo dc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  19,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ec_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  21,                                  /* lineNo */
  37,                                  /* colNo */
  "ux",                                /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo fc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  21,                                  /* lineNo */
  23,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo gc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  36,                                  /* colNo */
  "ux",                                /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo hc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  24,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ic_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  36,                                  /* lineNo */
  9,                                   /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo jc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  37,                                  /* lineNo */
  9,                                   /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo kc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  64,                                  /* lineNo */
  12,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo lc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  64,                                  /* lineNo */
  23,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo mc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  68,                                  /* lineNo */
  19,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo nc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  68,                                  /* lineNo */
  28,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo oc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  68,                                  /* lineNo */
  38,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo pc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  68,                                  /* lineNo */
  47,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo qc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  84,                                  /* lineNo */
  12,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo rc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  86,                                  /* lineNo */
  12,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo sc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  90,                                  /* lineNo */
  15,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo tc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  67,                                  /* lineNo */
  66,                                  /* colNo */
  "aplus",                             /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo uc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  67,                                  /* lineNo */
  75,                                  /* colNo */
  "aneg",                              /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo vc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  91,                                  /* lineNo */
  30,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo wc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  91,                                  /* lineNo */
  18,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo xc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  97,                                  /* lineNo */
  16,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo yc_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  98,                                  /* lineNo */
  33,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ad_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  98,                                  /* lineNo */
  19,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo bd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  116,                                 /* lineNo */
  17,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo cd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  116,                                 /* lineNo */
  26,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo dd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  116,                                 /* lineNo */
  47,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ed_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  116,                                 /* lineNo */
  56,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo fd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  117,                                 /* lineNo */
  17,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo gd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  117,                                 /* lineNo */
  26,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo hd_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  117,                                 /* lineNo */
  47,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo id_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  117,                                 /* lineNo */
  56,                                  /* colNo */
  "u",                                 /* aName */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo emlrtECI = { -1,    /* nDims */
  13,                                  /* lineNo */
  9,                                   /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo e_emlrtRTEI = { 13,/* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

static emlrtRTEInfo f_emlrtRTEI = { 10,/* lineNo */
  8,                                   /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo g_emlrtRTEI = { 10,/* lineNo */
  5,                                   /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 13,/* lineNo */
  24,                                  /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 16,/* lineNo */
  5,                                   /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo j_emlrtRTEI = { 17,/* lineNo */
  5,                                   /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo k_emlrtRTEI = { 30,/* lineNo */
  5,                                   /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo l_emlrtRTEI = { 31,/* lineNo */
  5,                                   /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo m_emlrtRTEI = { 40,/* lineNo */
  14,                                  /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

static emlrtRTEInfo n_emlrtRTEI = { 147,/* lineNo */
  38,                                  /* colNo */
  "unaryMinOrMax",                     /* fName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/eml/+coder/+internal/unaryMinOrMax.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 15,/* lineNo */
  13,                                  /* colNo */
  "isnan",                             /* fName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/elmat/isnan.m"/* pName */
};

static emlrtRTEInfo p_emlrtRTEI = { 16,/* lineNo */
  13,                                  /* colNo */
  "sum",                               /* fName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/sum.m"/* pName */
};

static emlrtRTEInfo q_emlrtRTEI = { 15,/* lineNo */
  9,                                   /* colNo */
  "sum",                               /* fName */
  "/Applications/MATLAB_R2020b.app/toolbox/eml/lib/matlab/datafun/sum.m"/* pName */
};

static emlrtRTEInfo r_emlrtRTEI = { 40,/* lineNo */
  25,                                  /* colNo */
  "KTSpace_HLLE",                      /* fName */
  "/Users/a1229158/Documents/GitHub/CITEPH-64-2012/Fns/SkeneOverwash/FloatingPlateSimulation/KTSpace_HLLE.m"/* pName */
};

/* Function Definitions */
void KTSpace_HLLE(const emlrtStack *sp, const emxArray_real_T *u, real_T dx,
                  const emxArray_real_T *x, real_T theta, real_T mindepth,
                  emxArray_real_T *spacediscrete, real_T *amax)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  emlrtStack st;
  emxArray_boolean_T *r1;
  emxArray_int32_T *nz;
  emxArray_real_T *aneg;
  emxArray_real_T *aplus;
  emxArray_real_T *b_x;
  emxArray_real_T *r;
  emxArray_real_T *uneg;
  emxArray_real_T *uplus;
  emxArray_real_T *ux;
  emxArray_real_T *varargin_1;
  real_T varargin_1_tmp[5];
  real_T b_varargin_1[3];
  real_T a[2];
  real_T b[2];
  real_T b_uplus[2];
  real_T a1_tmp;
  real_T a2_tmp;
  real_T alpha1;
  real_T alpha2;
  real_T cbar;
  real_T cl;
  real_T cr;
  real_T ex;
  real_T hl;
  real_T hneg_idx_0;
  real_T hneg_idx_1;
  real_T hr;
  real_T ubarl;
  real_T ubarr;
  real_T ul;
  real_T ur;
  int32_T iv[2];
  int32_T b_i;
  int32_T i;
  int32_T idx;
  int32_T k;
  int32_T loop_ub;
  boolean_T b_b;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  i_st.prev = &h_st;
  i_st.tls = h_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /*      u(1,:)=u(1,:)*pi*sqrt(2/24) */
  /*      u(2,:)=u(2,:)*sqrt(29/24)*sqrt(pi*sqrt(2/24)) */
  /* theta=1; */
  k = spacediscrete->size[0] * spacediscrete->size[1];
  spacediscrete->size[0] = 2;
  spacediscrete->size[1] = u->size[1];
  emxEnsureCapacity_real_T(sp, spacediscrete, k, &f_emlrtRTEI);
  loop_ub = u->size[0] * u->size[1];
  for (k = 0; k < loop_ub; k++) {
    spacediscrete->data[k] = u->data[k] * 0.0;
  }

  emxInit_real_T(sp, &ux, 2, &g_emlrtRTEI, true);
  k = ux->size[0] * ux->size[1];
  ux->size[0] = 2;
  ux->size[1] = spacediscrete->size[1];
  emxEnsureCapacity_real_T(sp, ux, k, &g_emlrtRTEI);
  loop_ub = spacediscrete->size[0] * spacediscrete->size[1];
  for (k = 0; k < loop_ub; k++) {
    ux->data[k] = spacediscrete->data[k];
  }

  emxInit_real_T(sp, &aneg, 2, &k_emlrtRTEI, true);
  loop_ub = u->size[1];
  iv[0] = 1;
  k = aneg->size[0] * aneg->size[1];
  aneg->size[0] = 1;
  aneg->size[1] = u->size[1];
  emxEnsureCapacity_real_T(sp, aneg, k, &h_emlrtRTEI);
  for (k = 0; k < loop_ub; k++) {
    aneg->data[k] = u->data[2 * k];
  }

  emxInit_real_T(sp, &b_x, 2, &n_emlrtRTEI, true);
  st.site = &emlrtRSI;
  minmod(&st, aneg, theta, b_x);
  iv[1] = ux->size[1];
  emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &b_x->size[0], 2, &emlrtECI, sp);
  idx = b_x->size[1];
  for (k = 0; k < idx; k++) {
    ux->data[2 * k] = b_x->data[k];
  }

  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(sp);
  }

  k = aneg->size[0] * aneg->size[1];
  aneg->size[0] = 1;
  aneg->size[1] = u->size[1];
  emxEnsureCapacity_real_T(sp, aneg, k, &h_emlrtRTEI);
  for (k = 0; k < loop_ub; k++) {
    aneg->data[k] = u->data[2 * k + 1];
  }

  st.site = &emlrtRSI;
  minmod(&st, aneg, theta, b_x);
  iv[1] = ux->size[1];
  emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &b_x->size[0], 2, &emlrtECI, sp);
  idx = b_x->size[1];
  for (k = 0; k < idx; k++) {
    ux->data[2 * k + 1] = b_x->data[k];
  }

  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(sp);
  }

  emxInit_real_T(sp, &uplus, 2, &i_emlrtRTEI, true);
  k = uplus->size[0] * uplus->size[1];
  uplus->size[0] = 2;
  uplus->size[1] = spacediscrete->size[1];
  emxEnsureCapacity_real_T(sp, uplus, k, &i_emlrtRTEI);
  loop_ub = spacediscrete->size[0] * spacediscrete->size[1];
  for (k = 0; k < loop_ub; k++) {
    uplus->data[k] = spacediscrete->data[k];
  }

  emxInit_real_T(sp, &uneg, 2, &j_emlrtRTEI, true);
  k = uneg->size[0] * uneg->size[1];
  uneg->size[0] = 2;
  uneg->size[1] = spacediscrete->size[1];
  emxEnsureCapacity_real_T(sp, uneg, k, &j_emlrtRTEI);
  loop_ub = spacediscrete->size[0] * spacediscrete->size[1];
  for (k = 0; k < loop_ub; k++) {
    uneg->data[k] = spacediscrete->data[k];
  }

  k = u->size[1] - 2;
  i = 2;
  for (b_i = 0; b_i <= k; b_i++) {
    i = b_i + 2;
    if ((b_i + 2 < 1) || (b_i + 2 > ux->size[1])) {
      emlrtDynamicBoundsCheckR2012b(b_i + 2, 1, ux->size[1], &gc_emlrtBCI, sp);
    }

    if ((b_i + 2 < 1) || (b_i + 2 > u->size[1])) {
      emlrtDynamicBoundsCheckR2012b(b_i + 2, 1, u->size[1], &hc_emlrtBCI, sp);
    }

    if ((b_i + 2 < 1) || (b_i + 2 > uplus->size[1])) {
      emlrtDynamicBoundsCheckR2012b(b_i + 2, 1, uplus->size[1], &tb_emlrtBCI, sp);
    }

    idx = 2 * (b_i + 1);
    uplus->data[idx] = u->data[idx] - 0.5 * ux->data[idx];
    uplus->data[idx + 1] = u->data[idx + 1] - 0.5 * ux->data[idx + 1];
    if (((int32_T)(((real_T)b_i + 2.0) - 1.0) < 1) || ((int32_T)(((real_T)b_i +
           2.0) - 1.0) > ux->size[1])) {
      emlrtDynamicBoundsCheckR2012b((int32_T)(((real_T)b_i + 2.0) - 1.0), 1,
        ux->size[1], &ec_emlrtBCI, sp);
    }

    if ((b_i + 1 < 1) || (b_i + 1 > u->size[1])) {
      emlrtDynamicBoundsCheckR2012b(b_i + 1, 1, u->size[1], &fc_emlrtBCI, sp);
    }

    if ((b_i + 2 < 1) || (b_i + 2 > uneg->size[1])) {
      emlrtDynamicBoundsCheckR2012b(b_i + 2, 1, uneg->size[1], &wb_emlrtBCI, sp);
    }

    uneg->data[idx] = u->data[2 * b_i] + 0.5 * ux->data[2 * b_i];
    loop_ub = 2 * b_i + 1;
    uneg->data[idx + 1] = u->data[loop_ub] + 0.5 * ux->data[loop_ub];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  if ((i - 1 < 1) || (i - 1 > ux->size[1])) {
    emlrtDynamicBoundsCheckR2012b(i - 1, 1, ux->size[1], &cc_emlrtBCI, sp);
  }

  if ((i - 1 < 1) || (i - 1 > u->size[1])) {
    emlrtDynamicBoundsCheckR2012b(i - 1, 1, u->size[1], &dc_emlrtBCI, sp);
  }

  if ((i < 1) || (i > uneg->size[1])) {
    emlrtDynamicBoundsCheckR2012b(i, 1, uneg->size[1], &vb_emlrtBCI, sp);
  }

  k = 2 * (i - 2);
  idx = 2 * (i - 1);
  uneg->data[idx] = u->data[k] + 0.5 * ux->data[k];
  uneg->data[idx + 1] = u->data[k + 1] + 0.5 * ux->data[k + 1];
  if (i > u->size[1]) {
    emlrtDynamicBoundsCheckR2012b(i, 1, u->size[1], &bc_emlrtBCI, sp);
  }

  if (i > uplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(i, 1, uplus->size[1], &sb_emlrtBCI, sp);
  }

  uplus->data[idx] = u->data[idx];
  uplus->data[idx + 1] = u->data[idx + 1];
  if (1 > uplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, uplus->size[1], &rb_emlrtBCI, sp);
  }

  if (1 > u->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, u->size[1], &ac_emlrtBCI, sp);
  }

  if (1 > ux->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, ux->size[1], &yb_emlrtBCI, sp);
  }

  uplus->data[0] = u->data[0] - 0.5 * ux->data[0];
  uplus->data[1] = u->data[1] - 0.5 * ux->data[1];
  emxFree_real_T(&ux);
  if (1 > uneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, uneg->size[1], &ub_emlrtBCI, sp);
  }

  if (1 > u->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, u->size[1], &xb_emlrtBCI, sp);
  }

  uneg->data[0] = u->data[0];
  uneg->data[1] = u->data[1];
  loop_ub = u->size[1];
  k = aneg->size[0] * aneg->size[1];
  aneg->size[0] = 1;
  aneg->size[1] = u->size[1];
  emxEnsureCapacity_real_T(sp, aneg, k, &k_emlrtRTEI);
  for (k = 0; k < loop_ub; k++) {
    aneg->data[k] = u->data[2 * k] * 0.0;
  }

  emxInit_real_T(sp, &aplus, 2, &l_emlrtRTEI, true);
  k = aplus->size[0] * aplus->size[1];
  aplus->size[0] = 1;
  aplus->size[1] = aneg->size[1];
  emxEnsureCapacity_real_T(sp, aplus, k, &l_emlrtRTEI);
  loop_ub = aneg->size[0] * aneg->size[1];
  for (k = 0; k < loop_ub; k++) {
    aplus->data[k] = aneg->data[k];
  }

  if (uneg->size[1] == 0) {
    b_i = 0;
  } else {
    b_i = muIntScalarMax_sint32(2, uneg->size[1]);
  }

  if (0 <= b_i - 1) {
    varargin_1_tmp[4] = 0.0;
  }

  for (i = 0; i < b_i; i++) {
    if (i + 1 > uplus->size[1]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, uplus->size[1], &qb_emlrtBCI, sp);
    }

    b_uplus[0] = uplus->data[2 * i];
    b_uplus[1] = uplus->data[1 + 2 * i];
    st.site = &b_emlrtRSI;
    spectralradius(b_uplus, a);
    if (i + 1 > uneg->size[1]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, uneg->size[1], &pb_emlrtBCI, sp);
    }

    b_uplus[0] = uneg->data[2 * i];
    b_uplus[1] = uneg->data[1 + 2 * i];
    st.site = &c_emlrtRSI;
    spectralradius(b_uplus, b);
    st.site = &d_emlrtRSI;
    varargin_1_tmp[0] = a[0];
    varargin_1_tmp[2] = b[0];
    varargin_1_tmp[1] = a[1];
    varargin_1_tmp[3] = b[1];
    b_st.site = &t_emlrtRSI;
    c_st.site = &u_emlrtRSI;
    d_st.site = &v_emlrtRSI;
    e_st.site = &w_emlrtRSI;
    f_st.site = &x_emlrtRSI;
    g_st.site = &cb_emlrtRSI;
    b_b = muDoubleScalarIsNaN(a[0]);
    if (!b_b) {
      idx = 1;
    } else {
      idx = 0;
      h_st.site = &db_emlrtRSI;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= 5)) {
        if (!muDoubleScalarIsNaN(varargin_1_tmp[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      if (i + 1 > aneg->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, aneg->size[1], &ic_emlrtBCI,
          &f_st);
      }

      aneg->data[i] = a[0];
    } else {
      g_st.site = &y_emlrtRSI;
      ex = varargin_1_tmp[idx - 1];
      loop_ub = idx + 1;
      h_st.site = &ab_emlrtRSI;
      for (k = loop_ub; k < 6; k++) {
        alpha2 = varargin_1_tmp[k - 1];
        if (ex > alpha2) {
          ex = alpha2;
        }
      }

      if (i + 1 > aneg->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, aneg->size[1], &ic_emlrtBCI,
          &f_st);
      }

      aneg->data[i] = ex;
    }

    st.site = &e_emlrtRSI;
    b_st.site = &eb_emlrtRSI;
    c_st.site = &fb_emlrtRSI;
    d_st.site = &gb_emlrtRSI;
    e_st.site = &hb_emlrtRSI;
    f_st.site = &ib_emlrtRSI;
    g_st.site = &cb_emlrtRSI;
    if (!b_b) {
      idx = 1;
    } else {
      idx = 0;
      h_st.site = &db_emlrtRSI;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= 5)) {
        if (!muDoubleScalarIsNaN(varargin_1_tmp[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      if (i + 1 > aplus->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, aplus->size[1], &jc_emlrtBCI,
          &f_st);
      }

      aplus->data[i] = a[0];
    } else {
      g_st.site = &y_emlrtRSI;
      ex = varargin_1_tmp[idx - 1];
      loop_ub = idx + 1;
      h_st.site = &ab_emlrtRSI;
      for (k = loop_ub; k < 6; k++) {
        alpha2 = varargin_1_tmp[k - 1];
        if (ex < alpha2) {
          ex = alpha2;
        }
      }

      if (i + 1 > aplus->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, aplus->size[1], &jc_emlrtBCI,
          &f_st);
      }

      aplus->data[i] = ex;
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxInit_real_T(sp, &varargin_1, 2, &m_emlrtRTEI, true);
  emxInit_real_T(sp, &r, 2, &r_emlrtRTEI, true);
  st.site = &f_emlrtRSI;
  b_st.site = &f_emlrtRSI;
  b_abs(&b_st, aneg, b_x);
  b_st.site = &f_emlrtRSI;
  b_abs(&b_st, aplus, r);
  k = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = 1;
  varargin_1->size[1] = b_x->size[1] + r->size[1];
  emxEnsureCapacity_real_T(&st, varargin_1, k, &m_emlrtRTEI);
  loop_ub = b_x->size[1];
  for (k = 0; k < loop_ub; k++) {
    varargin_1->data[k] = b_x->data[k];
  }

  loop_ub = r->size[1];
  for (k = 0; k < loop_ub; k++) {
    varargin_1->data[k + b_x->size[1]] = r->data[k];
  }

  emxFree_real_T(&r);
  b_st.site = &eb_emlrtRSI;
  c_st.site = &fb_emlrtRSI;
  d_st.site = &gb_emlrtRSI;
  if (varargin_1->size[1] < 1) {
    emlrtErrorWithMessageIdR2018a(&d_st, &emlrtRTEI,
      "Coder:toolbox:eml_min_or_max_varDimZero",
      "Coder:toolbox:eml_min_or_max_varDimZero", 0);
  }

  e_st.site = &hb_emlrtRSI;
  f_st.site = &ib_emlrtRSI;
  b_i = varargin_1->size[1];
  if (varargin_1->size[1] <= 2) {
    if (varargin_1->size[1] == 1) {
      ex = varargin_1->data[0];
    } else if ((varargin_1->data[0] < varargin_1->data[1]) ||
               (muDoubleScalarIsNaN(varargin_1->data[0]) &&
                (!muDoubleScalarIsNaN(varargin_1->data[1])))) {
      ex = varargin_1->data[1];
    } else {
      ex = varargin_1->data[0];
    }
  } else {
    g_st.site = &cb_emlrtRSI;
    if (!muDoubleScalarIsNaN(varargin_1->data[0])) {
      idx = 1;
    } else {
      idx = 0;
      h_st.site = &db_emlrtRSI;
      if (varargin_1->size[1] > 2147483646) {
        i_st.site = &bb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= varargin_1->size[1])) {
        if (!muDoubleScalarIsNaN(varargin_1->data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      ex = varargin_1->data[0];
    } else {
      g_st.site = &y_emlrtRSI;
      ex = varargin_1->data[idx - 1];
      loop_ub = idx + 1;
      h_st.site = &ab_emlrtRSI;
      if ((idx + 1 <= varargin_1->size[1]) && (varargin_1->size[1] > 2147483646))
      {
        i_st.site = &bb_emlrtRSI;
        check_forloop_overflow_error(&i_st);
      }

      for (k = loop_ub; k <= b_i; k++) {
        alpha2 = varargin_1->data[k - 1];
        if (ex < alpha2) {
          ex = alpha2;
        }
      }
    }
  }

  emxFree_real_T(&varargin_1);
  if (2 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aplus->size[1], &ob_emlrtBCI, sp);
  }

  if (2 > uneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, uneg->size[1], &nb_emlrtBCI, sp);
  }

  if (2 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aneg->size[1], &mb_emlrtBCI, sp);
  }

  if (2 > uplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, uplus->size[1], &lb_emlrtBCI, sp);
  }

  if (2 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aplus->size[1], &kb_emlrtBCI, sp);
  }

  if (2 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aneg->size[1], &jb_emlrtBCI, sp);
  }

  if (2 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aplus->size[1], &ib_emlrtBCI, sp);
  }

  if (2 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aneg->size[1], &hb_emlrtBCI, sp);
  }

  if (2 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aplus->size[1], &gb_emlrtBCI, sp);
  }

  if (2 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aneg->size[1], &fb_emlrtBCI, sp);
  }

  if (2 > uplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, uplus->size[1], &eb_emlrtBCI, sp);
  }

  if (2 > uneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, uneg->size[1], &db_emlrtBCI, sp);
  }

  st.site = &g_emlrtRSI;
  f(*(real_T (*)[2])&uneg->data[2], b_uplus);
  a[0] = b_uplus[0];
  a[1] = b_uplus[1];
  st.site = &g_emlrtRSI;
  f(*(real_T (*)[2])&uplus->data[2], b_uplus);
  ubarr = aplus->data[1] - aneg->data[1];
  cbar = aplus->data[1] * aneg->data[1] / ubarr;
  ubarl = aplus->data[1];
  hl = aneg->data[1];
  a[0] = (ubarl * a[0] - hl * b_uplus[0]) / ubarr + cbar * (uplus->data[2] -
    uneg->data[2]);
  a[1] = (ubarl * a[1] - hl * b_uplus[1]) / ubarr + cbar * (uplus->data[3] -
    uneg->data[3]);
  if (1 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aplus->size[1], &cb_emlrtBCI, sp);
  }

  if (1 > uneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, uneg->size[1], &bb_emlrtBCI, sp);
  }

  if (1 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aneg->size[1], &ab_emlrtBCI, sp);
  }

  if (1 > uplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, uplus->size[1], &y_emlrtBCI, sp);
  }

  if (1 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aplus->size[1], &x_emlrtBCI, sp);
  }

  if (1 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aneg->size[1], &w_emlrtBCI, sp);
  }

  if (1 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aplus->size[1], &v_emlrtBCI, sp);
  }

  if (1 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aneg->size[1], &u_emlrtBCI, sp);
  }

  if (1 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aplus->size[1], &t_emlrtBCI, sp);
  }

  if (1 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aneg->size[1], &s_emlrtBCI, sp);
  }

  if (1 > uplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, uplus->size[1], &r_emlrtBCI, sp);
  }

  if (1 > uneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, uneg->size[1], &q_emlrtBCI, sp);
  }

  st.site = &h_emlrtRSI;
  f(*(real_T (*)[2])&uneg->data[0], b_uplus);
  st.site = &h_emlrtRSI;
  f(*(real_T (*)[2])&uplus->data[0], b_uplus);
  if (1 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aplus->size[1], &p_emlrtBCI, sp);
  }

  if (1 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, aneg->size[1], &o_emlrtBCI, sp);
  }

  if (aplus->data[0] - aneg->data[0] == 0.0) {
    a[0] = 0.0;
    a[1] = 0.0;
  }

  if (2 > aplus->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aplus->size[1], &n_emlrtBCI, sp);
  }

  if (2 > aneg->size[1]) {
    emlrtDynamicBoundsCheckR2012b(2, 1, aneg->size[1], &m_emlrtBCI, sp);
  }

  if (1 > spacediscrete->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, spacediscrete->size[1], &d_emlrtBCI, sp);
  }

  if (u->size[1] == 0) {
    b_i = -3;
  } else if (2 > u->size[1]) {
    b_i = -1;
  } else {
    b_i = u->size[1] - 3;
  }

  for (i = 0; i <= b_i; i++) {
    hneg_idx_0 = a[0];
    hneg_idx_1 = a[1];
    if ((i + 3 < 1) || (i + 3 > aplus->size[1])) {
      emlrtDynamicBoundsCheckR2012b(i + 3, 1, aplus->size[1], &kc_emlrtBCI, sp);
    }

    if ((i + 3 < 1) || (i + 3 > aneg->size[1])) {
      emlrtDynamicBoundsCheckR2012b(i + 3, 1, aneg->size[1], &lc_emlrtBCI, sp);
    }

    if (aplus->data[i + 2] - aneg->data[i + 2] == 0.0) {
      alpha2 = 0.0;
      a[0] = 0.0;
      ubarr = 0.0;
      a[1] = 0.0;
    } else {
      if ((i + 3 < 1) || (i + 3 > aplus->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, aplus->size[1], &l_emlrtBCI, sp);
      }

      if ((i + 3 < 1) || (i + 3 > uneg->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, uneg->size[1], &k_emlrtBCI, sp);
      }

      st.site = &i_emlrtRSI;
      f(*(real_T (*)[2])&uneg->data[2 * (i + 2)], b_uplus);
      alpha2 = b_uplus[0];
      ubarr = b_uplus[1];
      if ((i + 3 < 1) || (i + 3 > aneg->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, aneg->size[1], &j_emlrtBCI, sp);
      }

      if ((i + 3 < 1) || (i + 3 > uplus->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, uplus->size[1], &i_emlrtBCI, sp);
      }

      st.site = &i_emlrtRSI;
      f(*(real_T (*)[2])&uplus->data[2 * (i + 2)], b_uplus);
      if ((i + 3 < 1) || (i + 3 > aplus->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, aplus->size[1], &mc_emlrtBCI, sp);
      }

      if ((i + 3 < 1) || (i + 3 > aneg->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, aneg->size[1], &nc_emlrtBCI, sp);
      }

      if ((i + 3 < 1) || (i + 3 > aplus->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, aplus->size[1], &oc_emlrtBCI, sp);
      }

      if ((i + 3 < 1) || (i + 3 > aneg->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, aneg->size[1], &pc_emlrtBCI, sp);
      }

      cbar = aplus->data[i + 2] * aneg->data[i + 2] / (aplus->data[i + 2] -
        aneg->data[i + 2]);
      if ((i + 3 < 1) || (i + 3 > uplus->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, uplus->size[1], &h_emlrtBCI, sp);
      }

      if ((i + 3 < 1) || (i + 3 > uneg->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, uneg->size[1], &g_emlrtBCI, sp);
      }

      ubarl = aplus->data[i + 2];
      hl = aneg->data[i + 2];
      if ((i + 3 < 1) || (i + 3 > aplus->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, aplus->size[1], &tc_emlrtBCI, sp);
      }

      if ((i + 3 < 1) || (i + 3 > aneg->size[1])) {
        emlrtDynamicBoundsCheckR2012b(i + 3, 1, aneg->size[1], &uc_emlrtBCI, sp);
      }

      hr = aplus->data[i + 2] - aneg->data[i + 2];
      k = 2 * (i + 2);
      alpha2 = (ubarl * alpha2 - hl * b_uplus[0]) / hr + cbar * (uplus->data[k]
        - uneg->data[k]);
      a[0] = alpha2;
      ubarr = (ubarl * ubarr - hl * b_uplus[1]) / hr + cbar * (uplus->data[k + 1]
        - uneg->data[k + 1]);
      a[1] = ubarr;
    }

    if ((i + 2 < 1) || (i + 2 > spacediscrete->size[1])) {
      emlrtDynamicBoundsCheckR2012b(i + 2, 1, spacediscrete->size[1],
        &c_emlrtBCI, sp);
    }

    k = 2 * (i + 1);
    spacediscrete->data[k] = -(alpha2 - hneg_idx_0) / dx;
    spacediscrete->data[k + 1] = -(ubarr - hneg_idx_1) / dx;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(&aplus);
  emxFree_real_T(&uneg);
  emxFree_real_T(&uplus);

  /* spacediscrete=spacediscrete*0; */
  if (1 > spacediscrete->size[1]) {
    emlrtDynamicBoundsCheckR2012b(1, 1, spacediscrete->size[1], &b_emlrtBCI, sp);
  }

  spacediscrete->data[0] = 0.0;
  a[0] = 0.0;
  spacediscrete->data[1] = 0.0;
  a[1] = 0.0;
  emxInit_int32_T(sp, &nz, 2, &q_emlrtRTEI, true);
  emxInit_boolean_T(sp, &r1, 2, &o_emlrtRTEI, true);
  for (i = 0; i < 2; i++) {
    /* ,length(x)-2,length(x)-1] */
    hneg_idx_0 = a[0];
    hneg_idx_1 = a[1];
    if (i + 1 > u->size[1]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &qc_emlrtBCI, sp);
    }

    hl = u->data[2 * i];
    cl = muDoubleScalarAbs(hl * 9.81);
    st.site = &j_emlrtRSI;
    if (cl < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    cl = muDoubleScalarSqrt(cl);
    if (i + 2 > u->size[1]) {
      emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &rc_emlrtBCI, sp);
    }

    b_i = 2 * (i + 1);
    hr = u->data[b_i];
    cr = muDoubleScalarAbs(hr * 9.81);
    st.site = &k_emlrtRSI;
    if (cr < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    cr = muDoubleScalarSqrt(cr);
    if (hl > mindepth) {
      if (i + 1 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &sc_emlrtBCI, sp);
      }

      idx = 2 * i + 1;
      ul = u->data[idx] / u->data[2 * i];
      if (i + 1 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &vc_emlrtBCI, sp);
      }

      alpha2 = u->data[2 * i];
      st.site = &l_emlrtRSI;
      if (alpha2 < 0.0) {
        emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      alpha2 = muDoubleScalarSqrt(alpha2);
      if (i + 1 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &wc_emlrtBCI, sp);
      }

      ubarl = u->data[idx] / alpha2;
    } else {
      ul = 0.0;
      ubarl = 0.0;
    }

    if (hr > mindepth) {
      if (i + 2 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &xc_emlrtBCI, sp);
      }

      ur = u->data[b_i + 1] / u->data[b_i];
      if (i + 2 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &yc_emlrtBCI, sp);
      }

      alpha2 = u->data[b_i];
      st.site = &m_emlrtRSI;
      if (alpha2 < 0.0) {
        emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      alpha2 = muDoubleScalarSqrt(alpha2);
      if (i + 2 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &ad_emlrtBCI, sp);
      }

      ubarr = u->data[b_i + 1] / alpha2;
    } else {
      ur = 0.0;
      ubarr = 0.0;
    }

    cbar = muDoubleScalarAbs(4.905 * (hl + hr));
    st.site = &n_emlrtRSI;
    if (cbar < 0.0) {
      emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    cbar = muDoubleScalarSqrt(cbar);
    if ((hl > 0.0) && (hr > 0.0)) {
      st.site = &o_emlrtRSI;
      if (hl < 0.0) {
        emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      st.site = &o_emlrtRSI;
      if (hr < 0.0) {
        emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
          "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3,
          4, 4, "sqrt");
      }

      hr = (ubarr + ubarl) / (muDoubleScalarSqrt(hl) + muDoubleScalarSqrt(hr));
    } else {
      hr = 0.0;
    }

    a1_tmp = hr + cbar;
    a2_tmp = hr - cbar;
    if (cbar != 0.0) {
      if (i + 2 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &bd_emlrtBCI, sp);
      }

      if (i + 1 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &cd_emlrtBCI, sp);
      }

      if (i + 2 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &dd_emlrtBCI, sp);
      }

      if (i + 1 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &ed_emlrtBCI, sp);
      }

      idx = 2 * (i + 1) + 1;
      loop_ub = 2 * i + 1;
      alpha1 = ((u->data[idx] - u->data[loop_ub]) + (-hr + cbar) * (u->data[b_i]
                 - u->data[2 * i])) / (2.0 * cbar);
      if (i + 2 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &fd_emlrtBCI, sp);
      }

      if (i + 1 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &gd_emlrtBCI, sp);
      }

      if (i + 2 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &hd_emlrtBCI, sp);
      }

      if (i + 1 > u->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &id_emlrtBCI, sp);
      }

      alpha2 = ((u->data[idx] - u->data[loop_ub]) + (-hr - cbar) * (u->data[b_i]
                 - u->data[2 * i])) / (-2.0 * cbar);
    } else {
      alpha1 = 0.0;
      alpha2 = 0.0;
    }

    ur = muDoubleScalarMax(ur + cr, a1_tmp);
    cr = muDoubleScalarMin(ul - cl, a2_tmp);
    ubarr = muDoubleScalarMax(ur, 0.0);
    hl = muDoubleScalarMin(cr, 0.0);
    if (ubarr > hl) {
      hr = ubarr + hl;
      ubarl = 2.0 * ubarr * hl;
      ubarr -= hl;
      hl = (hr * a1_tmp - ubarl) / ubarr;
      ubarl = (hr * a2_tmp - ubarl) / ubarr;
    } else {
      hl = 0.0;
      ubarl = 0.0;
    }

    cbar = hl * alpha1;
    ubarr = ubarl * alpha2;
    if (i + 2 > u->size[1]) {
      emlrtDynamicBoundsCheckR2012b(i + 2, 1, u->size[1], &f_emlrtBCI, sp);
    }

    st.site = &p_emlrtRSI;
    f(*(real_T (*)[2])&u->data[2 * (i + 1)], b_uplus);
    b[0] = b_uplus[0];
    b[1] = b_uplus[1];
    if (i + 1 > u->size[1]) {
      emlrtDynamicBoundsCheckR2012b(i + 1, 1, u->size[1], &e_emlrtBCI, sp);
    }

    st.site = &p_emlrtRSI;
    f(*(real_T (*)[2])&u->data[2 * i], b_uplus);
    alpha2 = 0.5 * (((b[0] + b_uplus[0]) - cbar) - ubarr);
    a[0] = alpha2;
    ubarr = 0.5 * (((b[1] + b_uplus[1]) - cbar * a1_tmp) - ubarr * a2_tmp);
    a[1] = ubarr;
    if ((i + 1 == 2) || (1 == x->size[1] - 1)) {
      if (i + 1 > spacediscrete->size[1]) {
        emlrtDynamicBoundsCheckR2012b(i + 1, 1, spacediscrete->size[1],
          &emlrtBCI, sp);
      }

      spacediscrete->data[2 * i] = -(alpha2 - hneg_idx_0) / dx;
      spacediscrete->data[2 * i + 1] = -(ubarr - hneg_idx_1) / dx;
    }

    st.site = &q_emlrtRSI;
    b_varargin_1[0] = ex;
    hr = muDoubleScalarAbs(ur);
    b_varargin_1[1] = hr;
    ubarr = muDoubleScalarAbs(cr);
    b_varargin_1[2] = ubarr;
    b_st.site = &eb_emlrtRSI;
    c_st.site = &fb_emlrtRSI;
    d_st.site = &gb_emlrtRSI;
    e_st.site = &hb_emlrtRSI;
    k = b_x->size[0] * b_x->size[1];
    b_x->size[0] = 1;
    b_x->size[1] = 3;
    emxEnsureCapacity_real_T(&e_st, b_x, k, &n_emlrtRTEI);
    b_x->data[0] = ex;
    b_x->data[1] = hr;
    b_x->data[2] = ubarr;
    f_st.site = &ib_emlrtRSI;
    g_st.site = &cb_emlrtRSI;
    if (!muDoubleScalarIsNaN(b_x->data[0])) {
      idx = 1;
    } else {
      idx = 0;
      h_st.site = &db_emlrtRSI;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= 3)) {
        if (!muDoubleScalarIsNaN(b_x->data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx != 0) {
      g_st.site = &y_emlrtRSI;
      ex = b_varargin_1[idx - 1];
      loop_ub = idx + 1;
      h_st.site = &ab_emlrtRSI;
      for (k = loop_ub; k < 4; k++) {
        alpha2 = b_varargin_1[k - 1];
        if (ex < alpha2) {
          ex = alpha2;
        }
      }
    }

    st.site = &r_emlrtRSI;
    b_st.site = &ob_emlrtRSI;
    k = r1->size[0] * r1->size[1];
    r1->size[0] = 2;
    r1->size[1] = spacediscrete->size[1];
    emxEnsureCapacity_boolean_T(&b_st, r1, k, &o_emlrtRTEI);
    loop_ub = spacediscrete->size[0] * spacediscrete->size[1];
    for (k = 0; k < loop_ub; k++) {
      r1->data[k] = muDoubleScalarIsNaN(spacediscrete->data[k]);
    }

    c_st.site = &pb_emlrtRSI;
    combineVectorElements(&c_st, r1, nz);
    k = aneg->size[0] * aneg->size[1];
    aneg->size[0] = 1;
    aneg->size[1] = nz->size[1];
    emxEnsureCapacity_real_T(sp, aneg, k, &p_emlrtRTEI);
    loop_ub = nz->size[0] * nz->size[1];
    for (k = 0; k < loop_ub; k++) {
      aneg->data[k] = nz->data[k];
    }

    st.site = &r_emlrtRSI;
    sum(&st, aneg);
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_boolean_T(&r1);
  emxFree_real_T(&b_x);
  emxFree_int32_T(&nz);
  emxFree_real_T(&aneg);
  *amax = ex;
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (KTSpace_HLLE.c) */
