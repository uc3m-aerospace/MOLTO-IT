/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_propagate_spirals_try_api.h
 *
 * Code generation for function '_coder_propagate_spirals_try_api'
 *
 */

#ifndef _CODER_PROPAGATE_SPIRALS_TRY_API_H
#define _CODER_PROPAGATE_SPIRALS_TRY_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_propagate_spirals_try_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void propagate_spirals_try(real_T v0, real_T r0, real_T theta0, real_T
  psi0, real_T thetaf, real_T ee, real_T *Time, real_T *v, real_T *r, real_T
  *theta, real_T *psi, real_T *flag);
extern void propagate_spirals_try_api(const mxArray * const prhs[6], const
  mxArray *plhs[6]);
extern void propagate_spirals_try_atexit(void);
extern void propagate_spirals_try_initialize(void);
extern void propagate_spirals_try_terminate(void);
extern void propagate_spirals_try_xil_terminate(void);

#endif

/* End of code generation (_coder_propagate_spirals_try_api.h) */
