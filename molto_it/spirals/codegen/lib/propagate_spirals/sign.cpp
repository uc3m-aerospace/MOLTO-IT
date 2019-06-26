//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sign.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 26-Jun-2019 16:58:51
//

// Include Files
#include "rt_nonfinite.h"
#include "propagate_spirals.h"
#include "sign.h"

// Function Definitions

//
// Arguments    : double *x
// Return Type  : void
//
void b_sign(double *x)
{
  if (*x < 0.0) {
    *x = -1.0;
  } else if (*x > 0.0) {
    *x = 1.0;
  } else {
    if (*x == 0.0) {
      *x = 0.0;
    }
  }
}

//
// File trailer for sign.cpp
//
// [EOF]
//
