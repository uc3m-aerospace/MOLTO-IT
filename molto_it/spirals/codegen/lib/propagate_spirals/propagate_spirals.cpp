//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: propagate_spirals.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 26-Jun-2019 16:58:51
//

// Include Files
#include <cmath>
#include <math.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "propagate_spirals.h"
#include "isequal.h"
#include "sqrt.h"
#include "asin.h"
#include "acosh.h"
#include "log1.h"
#include "asinh.h"
#include "sign.h"

// Function Declarations
static void ELIT(double HK, double PHI, double *FE, double *EE);
static double ELIT3(double PHI, double HK, double C);
static double rt_atan2d_snf(double u0, double u1);
static double rt_powd_snf(double u0, double u1);

// Function Definitions

//
// ==================================================
//        Purpose: Compute complete and incomplete elliptic
//                 integrals F(k,phi) and E(k,phi)
//        Input  : HK  --- Modulus k ( 0 � k � 1 )
//                 Phi --- Argument ( in degrees )
//        Output : FE  --- F(k,phi)
//                 EE  --- E(k,phi)
//        ==================================================
// Arguments    : double HK
//                double PHI
//                double *FE
//                double *EE
// Return Type  : void
//
static void ELIT(double HK, double PHI, double *FE, double *EE)
{
  double G;
  double A0;
  double R;
  double B0;
  double D0;
  double D;
  double FAC;
  int N;
  boolean_T exitg1;
  double A;
  double B;
  double C;
  G = 0.0;
  A0 = 1.0;
  R = HK * HK;
  B0 = std::sqrt(1.0 - R);
  D0 = 0.017453292519943278 * PHI;
  D = 0.0;
  if (isequal(HK, 1.0) && isequal(PHI, 90.0)) {
    *FE = rtNaN;
    *EE = 1.0;
  } else if (isequal(HK, 1.0)) {
    *EE = std::sin(D0);
    *FE = std::log((1.0 + *EE) / std::cos(D0));
  } else {
    FAC = 1.0;
    N = 0;
    exitg1 = false;
    while ((!exitg1) && (N < 40)) {
      A = (A0 + B0) / 2.0;
      B = std::sqrt(A0 * B0);
      C = (A0 - B0) / 2.0;
      FAC *= 2.0;
      R += FAC * C * C;
      if (!isequal(PHI, 90.0)) {
        D = D0 + std::atan(B0 / A0 * std::tan(D0));
        G += C * std::sin(D);
        D0 = D + 3.14159265358979 * std::floor(D / 3.14159265358979 + 0.5);
      }

      A0 = A;
      B0 = B;
      if (C < 1.0E-7) {
        exitg1 = true;
      } else {
        N++;
      }
    }

    A0 = 3.14159265358979 / (2.0 * A);
    *EE = 3.14159265358979 * (2.0 - R) / (4.0 * A);
    if (isequal(PHI, 90.0)) {
      *FE = A0;
    } else {
      *FE = D / (FAC * A);
      *EE = *FE * *EE / A0 + G;
    }
  }
}

//
// =========================================================
//        Purpose: Compute the elliptic integral of the third kind
//                 using Gauss-Legendre quadrature
//        Input :  Phi --- Argument ( in degrees )
//                  k  --- Modulus   ( 0 � k � 1.0 )
//                  c  --- Parameter ( 0 � c � 1.0 )
//        Output:  EL3 --- Value of the elliptic integral of the
//                         third kind
//        =========================================================
// Arguments    : double PHI
//                double HK
//                double C
// Return Type  : double
//
static double ELIT3(double PHI, double HK, double C)
{
  double EL3;
  boolean_T LB1;
  boolean_T LB2;
  double C1;
  int b_I;
  double C0;
  static const double dv0[10] = { 0.99312859918509488, 0.96397192727791381,
    0.912234428251326, 0.83911697182221878, 0.7463319064601508,
    0.636053680726515, 0.51086700195082713, 0.37370608871541949,
    0.2277858511416451, 0.076526521133497338 };

  double EL3_tmp;
  double b_EL3_tmp;
  static const double dv1[10] = { 0.017614007139152121, 0.040601429800386939,
    0.062672048334109068, 0.083276741576704755, 0.1019301198172404,
    0.1181945319615184, 0.13168863844917661, 0.142096109318382,
    0.14917298647260371, 0.15275338713072581 };

  //
  if (isequal(HK, 1.0) && (std::abs(PHI - 90.0) <= 1.0E-8)) {
    LB1 = true;
  } else {
    LB1 = false;
  }

  if (isequal(C, 1.0) && (std::abs(PHI - 90.0) <= 1.0E-8)) {
    LB2 = true;
  } else {
    LB2 = false;
  }

  //
  if (LB1 || LB2) {
    EL3 = rtNaN;
  } else {
    //
    C1 = 0.0087266462599716 * PHI;
    EL3 = 0.0;
    for (b_I = 0; b_I < 10; b_I++) {
      //
      C0 = C1 * dv0[b_I];

      //
      //
      //
      EL3_tmp = std::sin(C1 + C0);
      C0 = std::sin(C1 - C0);
      b_EL3_tmp = HK * HK;
      EL3 += dv1[b_I] * (1.0 / ((1.0 - C * EL3_tmp * EL3_tmp) * std::sqrt(1.0 -
        b_EL3_tmp * EL3_tmp * EL3_tmp)) + 1.0 / ((1.0 - C * C0 * C0) * std::sqrt
        (1.0 - b_EL3_tmp * C0 * C0)));

      //
    }

    EL3 *= C1;
  }

  return EL3;
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2((double)b_u0, (double)b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d6;
  double d7;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d6 = std::abs(u0);
    d7 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d6 == 1.0) {
        y = 1.0;
      } else if (d6 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d7 == 0.0) {
      y = 1.0;
    } else if (d7 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//               CONTROLLED GENERALIZED LOGARITHMIC SPIRALS
//                    Written by David Morante (UC3M)
//                         dmorante@ing.uc3m.es
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//  INPUTS ::
//         r0      :: Initial Radius
//         v0      :: Magnitude of the initial velocity
//         theta0  :: Initial Polar angle [radians]
//         psi0    :: Initial Flight Path angle
//         thetaf  :: Final value of the polar angle [radians]
//         ee      :: The spiral control parameter
//         Npoints :: Orientations at which the state is to be computed (scalar or vector)
//
//  OUTPUTS ::
//         theta  :: Npoints Polar angles between theta0 and thetaf
//         r      :: Radius at polar angles theta
//         v      :: Velocity at polar angles theta
//         psi    :: Flight Path angle at polar angles theta
//         Time   :: Time of flight at polar angles theta
//         Flag   :: Status of the propagation. Posible values:
//                1 --> Normal Propagation
//               -1 --> Asymptote exceeded. Reduce the value of thetaf
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  Equations are based on Roa et al. (Nov.15)
//  'Controlled generalized logarithmic spirals for low-thrust mission desing'
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  NOTE 1 ! ALL VARIBLES SHOULD NON-DIMENSIONAL GIVEN THAT mu=1
//  NOTE 2 ! ee = 1/2 correspond to the tangential steering case
//  NOTE 3 ! We consider that when phi0 = pi/2 we are in raising regime
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//  Initialized Status Variable
// Arguments    : double v0
//                double r0
//                double theta0
//                double psi0
//                double thetaf
//                double ee
//                double *Time
//                double *v
//                double *r
//                double *theta
//                double *psi
//                double *flag
// Return Type  : void
//
void propagate_spirals(double v0, double r0, double theta0, double psi0, double
  thetaf, double ee, double *Time, double *v, double *r, double *theta, double
  *psi, double *flag)
{
  double K1_tmp;
  double b_K1_tmp;
  double K1;
  double kp;
  double K2;
  double Hf;
  double regime0;
  double ell;
  double beta;
  double d0;
  double d1;
  double rmin;
  double thetam;
  int regime;
  double d2;
  double kk;
  double nn;
  double phi0;
  double nn_tmp;
  double phi;
  double E0;
  double E;
  double d3;
  double d4;
  double d5;
  double Time_tmp;
  *flag = 1.0;
  *Time = 0.0;

  //
  //  Get the constants of Motion
  //
  K1_tmp = v0 * v0;
  b_K1_tmp = 2.0 * (1.0 - ee);
  K1 = K1_tmp - b_K1_tmp / r0;
  kp = std::sin(psi0);
  K2 = K1_tmp * r0 * kp;
  *theta = thetaf;

  //
  //  if K2 < 0
  //      flag = -1;
  //  end
  //  %
  //  if K1 < -abs(eps) && K2 > 2 * ( 1 - ee )
  //      flag = -1;
  //  end
  //  %
  //  if K1 < -abs(eps) && K2 > K1 * r0 + 2 * ( 1 - ee )
  //      flag = -1;
  //  end
  //  %
  //  if K1 > abs(eps) && r0 < ( K2 - 2 * ( 1 - ee )) / K1
  //      flag = -1;
  //  end
  //
  //  if abs(K1)< abs(eps)&& K2 ==1 && ee == 0.5
  //      flag = -1;
  //  end
  //
  //  Get the initial Regime (+1--> Raising; -1 --> Lowering)
  //
  Hf = std::cos(psi0);
  regime0 = Hf;
  b_sign(&regime0);

  //
  if (K1 < 0.0) {
    //  Elliptical spiral
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //  THE PARTICLE NEVER SCAPE (rmax) when propagated backward or forward
    //  The particle reach the origin of the central body
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
    //  Maximum radius
    //
    beta = (b_K1_tmp - K2) / -K1;

    //
    ell = 4.0 * ((1.0 - ee) * (1.0 - ee)) - K2 * K2;
    b_sqrt(&ell);

    //
    //  Axis of symmetry
    //
    d0 = beta / r0;
    d0 -= 2.0 * (1.0 - ee) / K2 * (1.0 - d0);
    b_acosh(&d0);
    thetam = theta0 + regime0 * K2 / ell * std::abs(d0);

    //
    //  Regime at each theta
    //
    regime = (thetaf <= thetam) - (thetaf > thetam);

    //
    //  Compute Trajectory
    //
    *r = beta * (2.0 * (1.0 - ee) + K2) / (b_K1_tmp + K2 * std::cosh(ell / K2 *
      (thetaf - thetam)));
    *v = K1 + 2.0 / *r * (1.0 - ee);
    b_sqrt(v);

    //
    Hf = b_K1_tmp + K1 * *r;
    d0 = Hf * Hf - K2 * K2;
    b_sqrt(&d0);
    *psi = rt_atan2d_snf(K2 / Hf, (double)regime * (d0 / Hf));

    //
    //  Time of flight
    //
    thetam = K2 / beta;
    b_sqrt(&thetam);
    kk = -K1 * beta / (4.0 * (1.0 - ee));
    b_sqrt(&kk);
    kp = 1.0 - kk * kk;
    b_sqrt(&kp);
    nn = K1 * beta / (2.0 * K2);

    //
    //  Elliptic integrals
    //
    d0 = 2.0 / (1.0 + std::sin(psi0));
    b_sqrt(&d0);
    phi0 = thetam / v0 * d0;
    b_asin(&phi0);

    //
    ELIT(kk, 90.0, &Hf, &K1_tmp);
    ELIT(kk, phi0 * 180.0 / 3.14159265358979, &Hf, &nn_tmp);

    //
    ell = ELIT3(90.0, kk, nn);

    //
    //
    d0 = 2.0 / (1.0 + std::sin(*psi));
    b_sqrt(&d0);
    phi = thetam / *v * d0;
    b_asin(&phi);

    //
    ELIT(kk, phi * 180.0 / 3.14159265358979, &Hf, &beta);

    //
    //
    //
    d0 = (1.0 - std::sin(psi0)) / (1.0 + std::sin(psi0));
    b_sqrt(&d0);
    d1 = 1.0 - ee;
    b_sqrt(&d1);
    d2 = K2;
    b_sqrt(&d2);
    d3 = (1.0 - std::sin(*psi)) / (1.0 + std::sin(*psi));
    b_sqrt(&d3);
    d4 = 1.0 - ee;
    b_sqrt(&d4);
    d5 = K2;
    b_sqrt(&d5);
    Time_tmp = b_K1_tmp * (kp * kp);
    rmin = rt_powd_snf(-K1, 1.5);
    *Time = -regime0 * (r0 * v0 / K1 * d0 + 2.0 * (Time_tmp * (ELIT3(phi0 *
      180.0 / 3.14159265358979, kk, nn) - ell) - K2 * (nn_tmp - K1_tmp)) * d1 /
                        (rmin * d2)) + (double)regime * (*r * *v / K1 * d3 + 2.0
      * (Time_tmp * (ELIT3(phi * 180.0 / 3.14159265358979, kk, nn) - ell) - K2 *
         (beta - K1_tmp)) * d4 / (rmin * d5));
  } else if (K1 > 0.0) {
    //  Hiperbolic  Spiral
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //  THE PARTICLE REACHES INFINITY WITH A NON-ZERO VELOCITY
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
    if (K2 < b_K1_tmp) {
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      //  Type I Hiperbolic Spiral ( It has an asymthotic value )
      //  Raising regime[ reaches infinity ], Lowering regime[ reaches origin ]
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      //
      ell = 4.0 * ((1.0 - ee) * (1.0 - ee)) - K2 * K2;
      b_sqrt(&ell);

      //
      thetam = b_K1_tmp + ell;

      //
      //  Get the Position of the Asymptote
      //
      d0 = K2 * (((thetam - ell) - K2 * std::sin(psi0)) + ell * std::abs(std::
                  cos(psi0))) / (r0 * K1 * thetam * kp);
      b_log(&d0);
      Hf = theta0 + regime0 * K2 / ell * d0;

      //
      //  Check that the desired theta is lower than Asymtote
      //
      if ((theta0 < Hf) && (Hf < thetaf)) {
        *flag = -1.0;
      }

      //
      beta = regime0 * ell / K2 * (Hf - thetaf);

      //
      Hf = std::sinh(beta / 2.0);
      *r = thetam * (ell * ell) / (K1 * (Hf * (4.0 * thetam * (1.0 - ee) * Hf +
        (thetam * thetam - K2 * K2) * std::cosh(beta / 2.0))));
      *v = K1 + 2.0 / *r * (1.0 - ee);
      b_sqrt(v);

      //
      beta = K1 * *r;
      Hf = b_K1_tmp + beta;

      //
      d0 = Hf * Hf - K2 * K2;
      b_sqrt(&d0);
      *psi = rt_atan2d_snf(K2 / Hf, regime0 * (d0 / Hf));

      //
      //  Time of Flight
      //
      d0 = (2.0 * (1.0 - ee) + K2) / (1.0 - ee);
      b_sqrt(&d0);
      kk = 0.5 * d0;
      nn = (2.0 * (1.0 - ee) + K2) / (2.0 * K2);

      //  Ellliptic integrals
      Hf = nn * K2;
      phi0 = K1 * r0 * kp / (Hf * (1.0 - std::sin(psi0)));
      b_sqrt(&phi0);
      b_asin(&phi0);
      thetam = std::sin(*psi);
      phi = beta * thetam / (Hf * (1.0 - thetam));
      b_sqrt(&phi);
      b_asin(&phi);

      //
      //  Do only if any asymptote is reached
      //
      if (*flag == 1.0) {
        ELIT(kk, phi0 * 180.0 / 3.14159265358979, &Hf, &E0);
        ELIT(kk, phi * 180.0 / 3.14159265358979, &Hf, &E);
        d0 = (1.0 + kp) / (1.0 - std::sin(psi0));
        b_sqrt(&d0);
        d1 = K2 * (1.0 - ee);
        d2 = d1;
        b_sqrt(&d2);
        d3 = (1.0 + thetam) / (1.0 - std::sin(*psi));
        b_sqrt(&d3);
        b_sqrt(&d1);
        *Time = -regime0 * (r0 * v0 / K1 * d0 - 2.0 * (E0 - (1.0 - nn) * ELIT3
          (phi0 * 180.0 / 3.14159265358979, kk, nn)) * d2 / rt_powd_snf(K1, 1.5))
          + regime0 * (*r * *v / K1 * d3 - 2.0 * (E - (1.0 - nn) * ELIT3(phi *
          180.0 / 3.14159265358979, kk, nn)) * d1 / rt_powd_snf(K1, 1.5));
      }
    } else if (K2 > b_K1_tmp) {
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      //  Type II Hiperbolic Spiral
      //  There is a transition between raising and lowering regime [rmin]
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      //
      //  Minimum Radius
      //
      rmin = (K2 - b_K1_tmp) / K1;

      //
      ell = K2 * K2 - 4.0 * ((1.0 - ee) * (1.0 - ee));
      b_sqrt(&ell);

      //
      thetam = theta0 - regime0 * K2 / ell * (1.5707963267948966 + std::atan
        ((b_K1_tmp - K2 * kp) / (ell * std::abs(Hf))));

      //
      //  Get the position of the Asymptotes
      //
      Hf = K2 / ell * (1.5707963267948966 + std::atan(b_K1_tmp / ell));
      beta = thetam + Hf;
      Hf = thetam - Hf;

      //
      //  Check that the desired theta is lower than Asymtote
      //
      if (((thetaf > beta) && (beta > theta0)) || ((thetaf > Hf) && (Hf > theta0)))
      {
        *flag = -1.0;
      }

      //
      //
      *r = rmin * ((b_K1_tmp + K2) / (b_K1_tmp + K2 * std::cos(ell / K2 *
        (thetaf - thetam))));
      *v = K1 + 2.0 / *r * (1.0 - ee);
      b_sqrt(v);

      //
      regime = (thetaf >= thetam) - (thetaf < thetam);

      //
      Hf = b_K1_tmp + K1 * *r;

      //
      d0 = Hf * Hf - K2 * K2;
      b_sqrt(&d0);
      *psi = rt_atan2d_snf(K2 / Hf, (double)regime * (d0 / Hf));

      //
      //  Time of Flight
      //
      d0 = 1.0 - ee;
      b_sqrt(&d0);
      d1 = K2 + b_K1_tmp;
      d2 = d1;
      b_sqrt(&d2);
      kk = 2.0 * d0 / d2;
      nn_tmp = b_K1_tmp / K2;

      //
      //  Ellliptic integrals
      //
      d0 = b_K1_tmp * (1.0 - kp);
      b_sqrt(&d0);
      d2 = K2 - b_K1_tmp * kp;
      b_sqrt(&d2);
      phi0 = d0 / (kk * d2);
      b_asin(&phi0);
      d0 = b_K1_tmp * (1.0 - std::sin(*psi));
      b_sqrt(&d0);
      d2 = K2 - b_K1_tmp * std::sin(*psi);
      b_sqrt(&d2);
      phi = d0 / (kk * d2);
      b_asin(&phi);
      if (*flag == 1.0) {
        //
        ELIT(kk, phi0 * 180.0 / 3.14159265358979, &ell, &E0);

        //
        ELIT(kk, phi * 180.0 / 3.14159265358979, &kp, &E);

        //
        //
        d0 = K1 * K2 * d1;
        d2 = d0;
        b_sqrt(&d2);
        d3 = r0 * (v0 * v0) - K2;
        d4 = 2.0 * K1 * r0 * d3;
        b_sqrt(&d4);
        d3 = K2 * r0 * K1_tmp + d3 * (1.0 - ee);
        b_sqrt(&d3);
        d3 = d4 / (2.0 * d3);
        b_asinh(&d3);
        d4 = r0 * r0 * rt_powd_snf(v0, 4.0) - K2 * K2;
        b_sqrt(&d4);
        b_sqrt(&d0);
        d5 = *v * *v;
        thetam = *r * d5 - K2;
        Hf = 2.0 * K1 * *r * thetam;
        b_sqrt(&Hf);
        d5 = K2 * *r * d5 + thetam * (1.0 - ee);
        b_sqrt(&d5);
        d5 = Hf / (2.0 * d5);
        b_asinh(&d5);
        thetam = *r * *r * rt_powd_snf(*v, 4.0) - K2 * K2;
        b_sqrt(&thetam);
        Time_tmp = d1 * K2;
        rmin *= K1;
        Hf = b_K1_tmp / rt_powd_snf(K1, 1.5);
        beta = K1 * K1;
        *Time = regime0 * (((Time_tmp * E0 - rmin * (K2 * ell + b_K1_tmp * ELIT3
          (phi0 * 180.0 / 3.14159265358979, kk, nn_tmp))) / (K1 * d2) + Hf * d3)
                           + -v0 / beta * d4) - (double)regime * (((Time_tmp * E
          - rmin * (K2 * kp + b_K1_tmp * ELIT3(phi * 180.0 / 3.14159265358979,
          kk, nn_tmp))) / (K1 * d0) + Hf * d5) + -*v / beta * thetam);
      }
    } else {
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      //  Limit case Hiperbolyc spiral
      // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      //
      //  Get the asymptote
      //
      d0 = 4.0 * (1.0 - ee);
      d1 = 1.0 + d0 / (K1 * r0);
      b_sqrt(&d1);
      beta = thetaf - (theta0 - regime0 * (1.0 - d1));

      //
      *r = d0 / (K1 * beta * (beta - regime0 * 2.0));
      *v = K1 + 2.0 / *r * (1.0 - ee);
      b_sqrt(v);

      //
      Hf = b_K1_tmp + K1 * *r;

      //
      d0 = Hf * Hf - K2 * K2;
      b_sqrt(&d0);
      *psi = rt_atan2d_snf(K2 / Hf, regime0 * (d0 / Hf));

      //
      //  Time of flight
      //
      d0 = r0 * K1_tmp;
      d1 = r0 * (d0 + b_K1_tmp);
      b_sqrt(&d1);
      beta = v0 * d1;
      d1 = *r * (*v * *v);
      d2 = *r * (d1 + b_K1_tmp);
      b_sqrt(&d2);
      Hf = *v * d2;
      *Time = regime0 / rt_powd_snf(K1, 1.5) * ((Hf - beta) + (1.0 - ee) * std::
        log((((d0 + 1.0) - ee) + beta) / (((d1 + 1.0) - ee) + Hf)));
    }
  } else {
    //  Parabolic spiral
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //  THE FLIGHT PATH ANGLE REMAINS CONSTANT
    //  THE PARTICLE REACHES INFINITY WITH A ZERO VELOCITY
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //
    ell = std::sqrt(4.0 * ((1.0 - ee) * (1.0 - ee)) - K2 * K2);

    //
    *r = r0 * std::exp(regime0 * ell * (thetaf - theta0) / K2);
    *v = std::sqrt(K1 + 2.0 / *r * (1.0 - ee));

    //
    *psi = psi0;

    //
    //  Time of flight
    //
    *Time = regime0 * 2.0 * std::sqrt(b_K1_tmp) / 3.0 / ell * (rt_powd_snf(*r,
      1.5) - rt_powd_snf(r0, 1.5));

    //
  }

  //
  //
  if (rtIsNaN(*Time)) {
    *Time = -1.0;
    *flag = -1.0;
  }

  //
  if (*flag == -1.0) {
    *Time = -1.0;
    *v = -1.0;
    *r = -1.0;
    *theta = -1.0;
    *psi = -1.0;
  }

  //
}

//
// File trailer for propagate_spirals.cpp
//
// [EOF]
//
