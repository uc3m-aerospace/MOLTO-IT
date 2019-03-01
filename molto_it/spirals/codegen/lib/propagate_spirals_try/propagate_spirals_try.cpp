/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * propagate_spirals_try.cpp
 *
 * Code generation for function 'propagate_spirals_try'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "propagate_spirals_try.h"
#include "isequal.h"
#include "sqrt.h"
#include "asin.h"
#include "acosh.h"
#include "log1.h"
#include "asinh.h"
#include "sign.h"

/* Function Declarations */
static void ELIT(double HK, double PHI, double *FE, double *EE);
static double ELIT3(double PHI, double HK, double C);
static double rt_atan2d_snf(double u0, double u1);
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
static void ELIT(double HK, double PHI, double *FE, double *EE)
{
  double G;
  double A0;
  double B0;
  double D0;
  double R;
  double D;
  double FAC;
  int N;
  boolean_T exitg1;
  double A;
  double x;
  double C;

  /*  */
  /*        ================================================== */
  /*        Purpose: Compute complete and incomplete elliptic */
  /*                 integrals F(k,phi) and E(k,phi) */
  /*        Input  : HK  --- Modulus k ( 0 Û k Û 1 ) */
  /*                 Phi --- Argument ( in degrees ) */
  /*        Output : FE  --- F(k,phi) */
  /*                 EE  --- E(k,phi) */
  /*        ================================================== */
  /*  */
  G = 0.0;
  A0 = 1.0;
  B0 = std::sqrt(1.0 - HK * HK);
  D0 = 0.017453292519943278 * PHI;
  R = HK * HK;
  D = 0.0;
  if (isequal(HK, 1.0) && isequal(PHI, 90.0)) {
    *FE = rtNaN;
    *EE = 1.0;
  } else if (isequal(HK, 1.0)) {
    *FE = std::log((1.0 + std::sin(D0)) / std::cos(D0));
    *EE = std::sin(D0);
  } else {
    FAC = 1.0;
    N = 0;
    exitg1 = false;
    while ((!exitg1) && (N < 40)) {
      A = (A0 + B0) / 2.0;
      x = A0 * B0;
      C = (A0 - B0) / 2.0;
      FAC *= 2.0;
      R += FAC * C * C;
      if (!isequal(PHI, 90.0)) {
        D = D0 + std::atan(B0 / A0 * std::tan(D0));
        G += C * std::sin(D);
        D0 = D + 3.14159265358979 * std::floor(D / 3.14159265358979 + 0.5);
      }

      A0 = A;
      B0 = std::sqrt(x);
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

static double ELIT3(double PHI, double HK, double C)
{
  double EL3;
  boolean_T LB1;
  boolean_T LB2;
  double C1;
  int I;
  double C0;
  static const double dv0[10] = { 0.99312859918509488, 0.96397192727791381,
    0.912234428251326, 0.83911697182221878, 0.7463319064601508,
    0.636053680726515, 0.51086700195082713, 0.37370608871541949,
    0.2277858511416451, 0.076526521133497338 };

  double T1;
  static const double dv1[10] = { 0.017614007139152121, 0.040601429800386939,
    0.062672048334109068, 0.083276741576704755, 0.1019301198172404,
    0.1181945319615184, 0.13168863844917661, 0.142096109318382,
    0.14917298647260371, 0.15275338713072581 };

  /*  */
  /*        ========================================================= */
  /*        Purpose: Compute the elliptic integral of the third kind */
  /*                 using Gauss-Legendre quadrature */
  /*        Input :  Phi --- Argument ( in degrees ) */
  /*                  k  --- Modulus   ( 0 Û k Û 1.0 ) */
  /*                  c  --- Parameter ( 0 Û c Û 1.0 ) */
  /*        Output:  EL3 --- Value of the elliptic integral of the */
  /*                         third kind */
  /*        ========================================================= */
  /*  */
  /*  */
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

  /*  */
  if (LB1 || LB2) {
    EL3 = rtNaN;
  } else {
    /*  */
    C1 = 0.0087266462599716 * PHI;
    EL3 = 0.0;
    for (I = 0; I < 10; I++) {
      /*  */
      C0 = C1 * dv0[I];
      T1 = C1 + C0;
      C0 = C1 - C0;

      /*  */
      /*  */
      /*  */
      EL3 += dv1[I] * (1.0 / ((1.0 - C * std::sin(T1) * std::sin(T1)) * std::
        sqrt(1.0 - HK * HK * std::sin(T1) * std::sin(T1))) + 1.0 / ((1.0 - C *
        std::sin(C0) * std::sin(C0)) * std::sqrt(1.0 - HK * HK * std::sin(C0) *
        std::sin(C0))));

      /*  */
    }

    EL3 *= C1;
  }

  return EL3;
}

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

void propagate_spirals_try(double v0, double r0, double theta0, double psi0,
  double thetaf, double ee, double *Time, double *v, double *r, double *theta,
  double *psi, double *flag)
{
  double K1;
  double K2;
  double regime0;
  double ell;
  double H0;
  double d0;
  double rmin;
  double beta;
  double thetam;
  int regime;
  double kk;
  double kp;
  double nn;
  double phi0;
  double d1;
  double phi;
  double E01;
  double d2;
  double d3;
  double d4;
  double d5;

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*               CONTROLLED GENERALIZED LOGARITHMIC SPIRALS */
  /*                    Written by David Morante (UC3M) */
  /*                         dmorante@ing.uc3m.es */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  */
  /*  INPUTS :: */
  /*         r0      :: Initial Radius */
  /*         v0      :: Magnitude of the initial velocity */
  /*         theta0  :: Initial Polar angle [radians] */
  /*         psi0    :: Initial Flight Path angle */
  /*         thetaf  :: Final value of the polar angle [radians] */
  /*         ee      :: The spiral control parameter */
  /*         Npoints :: Orientations at which the state is to be computed (scalar or vector) */
  /*  */
  /*  OUTPUTS :: */
  /*         theta  :: Npoints Polar angles between theta0 and thetaf */
  /*         r      :: Radius at polar angles theta */
  /*         v      :: Velocity at polar angles theta */
  /*         psi    :: Flight Path angle at polar angles theta */
  /*         Time   :: Time of flight at polar angles theta */
  /*         Flag   :: Status of the propagation. Posible values: */
  /*                1 --> Normal Propagation */
  /*               -1 --> Asymptote exceeded. Reduce the value of thetaf */
  /*  */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  Equations are based on Roa et al. (Nov.15) */
  /*  'Controlled generalized logarithmic spirals for low-thrust mission desing' */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  NOTE 1 ! ALL VARIBLES SHOULD NON-DIMENSIONAL GIVEN THAT mu=1 */
  /*  NOTE 2 ! ee = 1/2 correspond to the tangential steering case */
  /*  NOTE 3 ! We consider that when phi0 = pi/2 we are in raising regime */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  */
  /*  Initialized Status Variable */
  /*  */
  *flag = 1.0;
  *Time = 0.0;

  /*  */
  /*  Get the constants of Motion */
  /*  */
  K1 = v0 * v0 - 2.0 * (1.0 - ee) / r0;
  K2 = v0 * v0 * r0 * std::sin(psi0);
  *theta = thetaf;

  /*  */
  /*  if K2 < 0 */
  /*      flag = -1; */
  /*  end */
  /*  % */
  /*  if K1 < -abs(eps) && K2 > 2 * ( 1 - ee ) */
  /*      flag = -1; */
  /*  end */
  /*  % */
  /*  if K1 < -abs(eps) && K2 > K1 * r0 + 2 * ( 1 - ee ) */
  /*      flag = -1; */
  /*  end */
  /*  % */
  /*  if K1 > abs(eps) && r0 < ( K2 - 2 * ( 1 - ee )) / K1 */
  /*      flag = -1; */
  /*  end */
  /*   */
  /*  if abs(K1)< abs(eps)&& K2 ==1 && ee == 0.5 */
  /*      flag = -1; */
  /*  end */
  /*  */
  /*  Get the initial Regime (+1--> Raising; -1 --> Lowering) */
  /*  */
  regime0 = std::cos(psi0);
  b_sign(&regime0);

  /*  */
  if (K1 < 0.0) {
    /*  Elliptical spiral */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  THE PARTICLE NEVER SCAPE (rmax) when propagated backward or forward */
    /*  The particle reach the origin of the central body */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  */
    /*  Maximum radius */
    /*  */
    H0 = (2.0 * (1.0 - ee) - K2) / -K1;

    /*  */
    ell = 4.0 * ((1.0 - ee) * (1.0 - ee)) - K2 * K2;
    b_sqrt(&ell);

    /*  */
    /*  Axis of symmetry */
    /*  */
    d0 = H0 / r0 - 2.0 * (1.0 - ee) / K2 * (1.0 - H0 / r0);
    b_acosh(&d0);
    thetam = theta0 + regime0 * K2 / ell * std::abs(d0);

    /*  */
    /*  Regime at each theta */
    /*  */
    regime = (thetaf <= thetam) - (thetaf > thetam);

    /*  */
    /*  Compute Trajectory */
    /*  */
    *r = H0 * (2.0 * (1.0 - ee) + K2) / (2.0 * (1.0 - ee) + K2 * std::cosh(ell /
      K2 * (thetaf - thetam)));
    *v = K1 + 2.0 / *r * (1.0 - ee);
    b_sqrt(v);

    /*  */
    beta = 2.0 * (1.0 - ee) + K1 * *r;
    d0 = beta * beta - K2 * K2;
    b_sqrt(&d0);
    *psi = rt_atan2d_snf(K2 / (2.0 * (1.0 - ee) + K1 * *r), (double)regime * (d0
      / (2.0 * (1.0 - ee) + K1 * *r)));

    /*  */
    /*  Time of flight */
    /*  */
    thetam = K2 / H0;
    b_sqrt(&thetam);
    kk = -K1 * H0 / (4.0 * (1.0 - ee));
    b_sqrt(&kk);
    kp = 1.0 - kk * kk;
    b_sqrt(&kp);
    nn = K1 * H0 / (2.0 * K2);

    /*  */
    /*  Elliptic integrals */
    /*  */
    d0 = 2.0 / (1.0 + std::sin(psi0));
    b_sqrt(&d0);
    phi0 = thetam / v0 * d0;
    b_asin(&phi0);

    /*  */
    ELIT(kk, 90.0, &beta, &E01);
    ELIT(kk, phi0 * 180.0 / 3.14159265358979, &beta, &rmin);

    /*  */
    ell = ELIT3(90.0, kk, nn);

    /*  */
    /*  */
    d0 = 2.0 / (1.0 + std::sin(*psi));
    b_sqrt(&d0);
    phi = thetam / *v * d0;
    b_asin(&phi);

    /*  */
    ELIT(kk, phi * 180.0 / 3.14159265358979, &beta, &H0);

    /*  */
    /*  */
    /*  */
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
    *Time = -regime0 * (r0 * v0 / K1 * d0 + 2.0 * (2.0 * (1.0 - ee) * (kp * kp) *
      (ELIT3(phi0 * 180.0 / 3.14159265358979, kk, nn) - ell) - K2 * (rmin - E01))
                        * d1 / (rt_powd_snf(-K1, 1.5) * d2)) + (double)regime *
      (*r * *v / K1 * d3 + 2.0 * (2.0 * (1.0 - ee) * (kp * kp) * (ELIT3(phi *
          180.0 / 3.14159265358979, kk, nn) - ell) - K2 * (H0 - E01)) * d4 /
       (rt_powd_snf(-K1, 1.5) * d5));
  } else if (K1 > 0.0) {
    /*  Hiperbolic  Spiral */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  THE PARTICLE REACHES INFINITY WITH A NON-ZERO VELOCITY */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  */
    if (K2 < 2.0 * (1.0 - ee)) {
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      /*  Type I Hiperbolic Spiral ( It has an asymthotic value ) */
      /*  Raising regime[ reaches infinity ], Lowering regime[ reaches origin ] */
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      /*  */
      ell = 4.0 * ((1.0 - ee) * (1.0 - ee)) - K2 * K2;
      b_sqrt(&ell);

      /*  */
      H0 = 2.0 * (1.0 - ee) + ell;

      /*  */
      /*  Get the Position of the Asymptote */
      /*  */
      d0 = K2 * (((H0 - ell) - K2 * std::sin(psi0)) + ell * std::abs(std::cos
                  (psi0))) / (r0 * K1 * H0 * std::sin(psi0));
      b_log(&d0);
      beta = theta0 + regime0 * K2 / ell * d0;

      /*  */
      /*  Check that the desired theta is lower than Asymtote */
      /*  */
      if ((theta0 < beta) && (beta < thetaf)) {
        *flag = -1.0;
      }

      /*  */
      beta = regime0 * ell / K2 * (beta - thetaf);

      /*  */
      *r = H0 * (ell * ell) / (K1 * (std::sinh(beta / 2.0) * (4.0 * H0 * (1.0 -
        ee) * std::sinh(beta / 2.0) + (H0 * H0 - K2 * K2) * std::cosh(beta / 2.0))));
      *v = K1 + 2.0 / *r * (1.0 - ee);
      b_sqrt(v);

      /*  */
      beta = 2.0 * (1.0 - ee) + K1 * *r;

      /*  */
      d0 = beta * beta - K2 * K2;
      b_sqrt(&d0);
      *psi = rt_atan2d_snf(K2 / (2.0 * (1.0 - ee) + K1 * *r), regime0 * (d0 /
        (2.0 * (1.0 - ee) + K1 * *r)));

      /*  */
      /*  Time of Flight */
      /*  */
      d0 = (2.0 * (1.0 - ee) + K2) / (1.0 - ee);
      b_sqrt(&d0);
      kk = 0.5 * d0;
      nn = (2.0 * (1.0 - ee) + K2) / (2.0 * K2);

      /*  Ellliptic integrals */
      phi0 = K1 * r0 * std::sin(psi0) / (nn * K2 * (1.0 - std::sin(psi0)));
      b_sqrt(&phi0);
      b_asin(&phi0);
      phi = K1 * *r * std::sin(*psi) / (nn * K2 * (1.0 - std::sin(*psi)));
      b_sqrt(&phi);
      b_asin(&phi);

      /*  */
      /*  Do only if any asymptote is reached */
      /*  */
      if (*flag == 1.0) {
        ELIT(kk, phi0 * 180.0 / 3.14159265358979, &beta, &ell);
        ELIT(kk, phi * 180.0 / 3.14159265358979, &beta, &thetam);
        d0 = (1.0 + std::sin(psi0)) / (1.0 - std::sin(psi0));
        b_sqrt(&d0);
        d1 = K2 * (1.0 - ee);
        b_sqrt(&d1);
        d2 = (1.0 + std::sin(*psi)) / (1.0 - std::sin(*psi));
        b_sqrt(&d2);
        d3 = K2 * (1.0 - ee);
        b_sqrt(&d3);
        *Time = -regime0 * (r0 * v0 / K1 * d0 - 2.0 * (ell - (1.0 - nn) * ELIT3
          (phi0 * 180.0 / 3.14159265358979, kk, nn)) * d1 / rt_powd_snf(K1, 1.5))
          + regime0 * (*r * *v / K1 * d2 - 2.0 * (thetam - (1.0 - nn) * ELIT3
          (phi * 180.0 / 3.14159265358979, kk, nn)) * d3 / rt_powd_snf(K1, 1.5));
      }
    } else if (K2 > 2.0 * (1.0 - ee)) {
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      /*  Type II Hiperbolic Spiral */
      /*  There is a transition between raising and lowering regime [rmin] */
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      /*  */
      /*  Minimum Radius */
      /*  */
      rmin = (K2 - 2.0 * (1.0 - ee)) / K1;

      /*  */
      ell = K2 * K2 - 4.0 * ((1.0 - ee) * (1.0 - ee));
      b_sqrt(&ell);

      /*  */
      thetam = theta0 - regime0 * K2 / ell * (1.5707963267948966 + std::atan
        ((2.0 * (1.0 - ee) - K2 * std::sin(psi0)) / (ell * std::abs(std::cos
        (psi0)))));

      /*  */
      /*  Get the position of the Asymptotes */
      /*  */
      beta = thetam + K2 / ell * (1.5707963267948966 + std::atan(2.0 * (1.0 - ee)
        / ell));
      H0 = thetam - K2 / ell * (1.5707963267948966 + std::atan(2.0 * (1.0 - ee) /
        ell));

      /*  */
      /*  Check that the desired theta is lower than Asymtote */
      /*  */
      if (((thetaf > beta) && (beta > theta0)) || ((thetaf > H0) && (H0 > theta0)))
      {
        *flag = -1.0;
      }

      /*  */
      /*  */
      *r = rmin * ((2.0 * (1.0 - ee) + K2) / (2.0 * (1.0 - ee) + K2 * std::cos
        (ell / K2 * (thetaf - thetam))));
      *v = K1 + 2.0 / *r * (1.0 - ee);
      b_sqrt(v);

      /*  */
      regime = (thetaf >= thetam) - (thetaf < thetam);

      /*  */
      beta = 2.0 * (1.0 - ee) + K1 * *r;

      /*  */
      d0 = beta * beta - K2 * K2;
      b_sqrt(&d0);
      *psi = rt_atan2d_snf(K2 / (2.0 * (1.0 - ee) + K1 * *r), (double)regime *
                           (d0 / (2.0 * (1.0 - ee) + K1 * *r)));

      /*  */
      /*  Time of Flight */
      /*  */
      d0 = 1.0 - ee;
      b_sqrt(&d0);
      d1 = K2 + 2.0 * (1.0 - ee);
      b_sqrt(&d1);
      kk = 2.0 * d0 / d1;
      nn = 2.0 * (1.0 - ee) / K2;

      /*  */
      /*  Ellliptic integrals */
      /*  */
      d0 = 2.0 * (1.0 - ee) * (1.0 - std::sin(psi0));
      b_sqrt(&d0);
      d1 = K2 - 2.0 * (1.0 - ee) * std::sin(psi0);
      b_sqrt(&d1);
      phi0 = d0 / (kk * d1);
      b_asin(&phi0);
      d0 = 2.0 * (1.0 - ee) * (1.0 - std::sin(*psi));
      b_sqrt(&d0);
      d1 = K2 - 2.0 * (1.0 - ee) * std::sin(*psi);
      b_sqrt(&d1);
      phi = d0 / (kk * d1);
      b_asin(&phi);
      if (*flag == 1.0) {
        /*  */
        ELIT(kk, phi0 * 180.0 / 3.14159265358979, &beta, &ell);

        /*  */
        ELIT(kk, phi * 180.0 / 3.14159265358979, &H0, &thetam);

        /*  */
        /*  */
        d0 = K1 * K2 * (K2 + 2.0 * (1.0 - ee));
        b_sqrt(&d0);
        d1 = 2.0 * K1 * r0 * (r0 * (v0 * v0) - K2);
        b_sqrt(&d1);
        d2 = K2 * r0 * (v0 * v0) + (r0 * (v0 * v0) - K2) * (1.0 - ee);
        b_sqrt(&d2);
        d1 /= 2.0 * d2;
        b_asinh(&d1);
        d2 = r0 * r0 * rt_powd_snf(v0, 4.0) - K2 * K2;
        b_sqrt(&d2);
        d3 = K1 * K2 * (K2 + 2.0 * (1.0 - ee));
        b_sqrt(&d3);
        d4 = 2.0 * K1 * *r * (*r * (*v * *v) - K2);
        b_sqrt(&d4);
        d5 = K2 * *r * (*v * *v) + (*r * (*v * *v) - K2) * (1.0 - ee);
        b_sqrt(&d5);
        d4 /= 2.0 * d5;
        b_asinh(&d4);
        d5 = *r * *r * rt_powd_snf(*v, 4.0) - K2 * K2;
        b_sqrt(&d5);
        *Time = regime0 * ((((K2 + 2.0 * (1.0 - ee)) * K2 * ell - K1 * rmin *
                             (K2 * beta + 2.0 * (1.0 - ee) * ELIT3(phi0 * 180.0 /
          3.14159265358979, kk, nn))) / (K1 * d0) + 2.0 * (1.0 - ee) /
                            rt_powd_snf(K1, 1.5) * d1) + -v0 / (K1 * K1) * d2) -
          (double)regime * ((((K2 + 2.0 * (1.0 - ee)) * K2 * thetam - K1 * rmin *
                              (K2 * H0 + 2.0 * (1.0 - ee) * ELIT3(phi * 180.0 /
          3.14159265358979, kk, nn))) / (K1 * d3) + 2.0 * (1.0 - ee) /
                             rt_powd_snf(K1, 1.5) * d4) + -*v / (K1 * K1) * d5);
      }
    } else {
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      /*  Limit case Hiperbolyc spiral */
      /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
      /*  */
      /*  Get the asymptote */
      /*  */
      d0 = 1.0 + 4.0 * (1.0 - ee) / (K1 * r0);
      b_sqrt(&d0);
      beta = thetaf - (theta0 - regime0 * (1.0 - d0));

      /*  */
      *r = 4.0 * (1.0 - ee) / (K1 * beta * (beta - regime0 * 2.0));
      *v = K1 + 2.0 / *r * (1.0 - ee);
      b_sqrt(v);

      /*  */
      beta = 2.0 * (1.0 - ee) + K1 * *r;

      /*  */
      d0 = beta * beta - K2 * K2;
      b_sqrt(&d0);
      *psi = rt_atan2d_snf(K2 / (2.0 * (1.0 - ee) + K1 * *r), regime0 * (d0 /
        (2.0 * (1.0 - ee) + K1 * *r)));

      /*  */
      /*  Time of flight */
      /*  */
      d0 = r0 * (r0 * (v0 * v0) + 2.0 * (1.0 - ee));
      b_sqrt(&d0);
      H0 = v0 * d0;
      d0 = *r * (*r * (*v * *v) + 2.0 * (1.0 - ee));
      b_sqrt(&d0);
      beta = *v * d0;
      *Time = regime0 / rt_powd_snf(K1, 1.5) * ((beta - H0) + (1.0 - ee) * std::
        log((((r0 * (v0 * v0) + 1.0) - ee) + H0) / (((*r * (*v * *v) + 1.0) - ee)
        + beta)));
    }
  } else {
    /*  Parabolic spiral */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  THE FLIGHT PATH ANGLE REMAINS CONSTANT */
    /*  THE PARTICLE REACHES INFINITY WITH A ZERO VELOCITY */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /*  */
    ell = std::sqrt(4.0 * ((1.0 - ee) * (1.0 - ee)) - K2 * K2);

    /*  */
    *r = r0 * std::exp(regime0 * ell * (thetaf - theta0) / K2);
    *v = std::sqrt(K1 + 2.0 / *r * (1.0 - ee));

    /*  */
    *psi = psi0;

    /*  */
    /*  Time of flight */
    /*  */
    *Time = regime0 * 2.0 * std::sqrt(2.0 * (1.0 - ee)) / 3.0 / ell *
      (rt_powd_snf(*r, 1.5) - rt_powd_snf(r0, 1.5));

    /*  */
  }

  /*  */
  /*  */
  if (rtIsNaN(*Time)) {
    *Time = -1.0;
    *flag = -1.0;
  }

  /*  */
  if (*flag == -1.0) {
    *Time = -1.0;
    *v = -1.0;
    *r = -1.0;
    *theta = -1.0;
    *psi = -1.0;
  }

  /*  */
}

/* End of code generation (propagate_spirals_try.cpp) */
