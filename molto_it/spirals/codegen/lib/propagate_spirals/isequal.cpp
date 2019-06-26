//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: isequal.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 26-Jun-2019 16:58:51
//

// Include Files
#include "rt_nonfinite.h"
#include "propagate_spirals.h"
#include "isequal.h"

// Function Definitions

//
// Arguments    : double varargin_1
//                double varargin_2
// Return Type  : boolean_T
//
boolean_T isequal(double varargin_1, double varargin_2)
{
  boolean_T p;
  boolean_T b_p;
  p = false;
  b_p = true;
  if (!(varargin_1 == varargin_2)) {
    b_p = false;
  }

  if (b_p) {
    p = true;
  }

  return p;
}

//
// File trailer for isequal.cpp
//
// [EOF]
//
