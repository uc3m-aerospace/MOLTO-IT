//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: propagate_spirals.h
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 26-Jun-2019 16:58:51
//
#ifndef PROPAGATE_SPIRALS_H
#define PROPAGATE_SPIRALS_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "propagate_spirals_types.h"

// Function Declarations
extern void propagate_spirals(double v0, double r0, double theta0, double psi0,
  double thetaf, double ee, double *Time, double *v, double *r, double *theta,
  double *psi, double *flag);

#endif

//
// File trailer for propagate_spirals.h
//
// [EOF]
//
