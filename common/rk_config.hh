/*!
\file rk_config.hh
\brief Runge-Kutta configurator (Butcher tables)
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_RK_CONFIG_HH
#define SPECTRUM_RK_CONFIG_HH

#include "config.h"

#ifdef RK_INTEGRATOR_TYPE

#include <string>

namespace Spectrum {

//! How many RK stages to expect
#define MAX_RK_STAGES 7

//! Absolute tolerance (used only to avoid divission by zero in computing the relative tolerance)
const double rk_tol_abs = 1.0E-9;

//! Relative tolerance
const double rk_tol_rel = 1.0E-9;

//! Safety factor: the step size cannot change by more than this factor during one iteration
const double rk_safety = 2.0;
 
template <uint8_t rk_stages> struct ButcherTable {

//! Readable name
   std::string name;

//! Number of stages
   uint8_t stages;

//! Formal order of convergence
   uint8_t order;

//! Adaptive or not
   bool adaptive;

//! Implicit or not
   bool implicit;

//! Time coefficients
   double a[rk_stages];

//! Slope coefficients
   double b[rk_stages][rk_stages];

//! Higher order weigts
   double v[rk_stages];

//! Lower order weigts
   double w[rk_stages];

//! Default constructor
   ButcherTable(void) = delete;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Order 1
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Euler method (1 E)
#if RK_INTEGRATOR_TYPE == 0

const uint8_t n_stages = 1;
const ButcherTable <n_stages> RK_Table = {"Euler first order explicit", n_stages, 1, false, false,
      {0.0},
      {0.0},
      {1.0},
      {1.0}};

//! Backward Euler method (1 I)
#elif RK_INTEGRATOR_TYPE == 1

const uint8_t n_stages = 1;
const ButcherTable <n_stages> RK_Table = {"Backward Euler first order implicit", n_stages, 1, false, true,
      {1.0},
      {1.0},
      {1.0},
      {1.0}};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Order 2
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Implicit midpoint method (2 I)
#elif RK_INTEGRATOR_TYPE == 2

const uint8_t n_stages = 1;
const ButcherTable <n_stages> RK_Table = {"Midpoint second order implicit", n_stages, 2, false, true,
      {1.0 / 2.0},
      {1.0 / 2.0},
      {1.0},
      {1.0}};

//! Midpoint method (2 E)
#elif RK_INTEGRATOR_TYPE == 3

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Midpoint second order explicit", n_stages, 2, false, false,
      {0.0, 1.0 / 2.0},
      {{0.0, 0.0},
       {1.0 / 2.0, 0.0}},
      {0.0, 1.0},
      {0.0, 1.0}};

//! Kraaijevanger-Spijker method (2 I), not convex
#elif RK_INTEGRATOR_TYPE == 4

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Kraaijevanger-Spijker second order implicit", n_stages, 2, false, true,
      {1.0 / 2.0, 3.0 / 2.0},
      {{1.0 / 2.0, 0.0},
       {-1.0 / 2.0, 2.0}},
      {-1.0 / 2.0, 3.0 / 2.0},
      {-1.0 / 2.0, 3.0 / 2.0}};

//! Qin-Zhang method (2 I)
#elif RK_INTEGRATOR_TYPE == 5

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Qin-Zhang second order implicit", n_stages, 2, false, true,
      {1.0 / 4.0, 3.0 / 4.0},
      {{1.0 / 4.0, 0.0},
       {1.0 / 2.0, 1.0 / 4.0}},
      {1.0 / 2.0, 1.0 / 2.0},
      {1.0 / 2.0, 1.0 / 2.0}};

//! Ralston's method (2 E)
#elif RK_INTEGRATOR_TYPE == 6

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Ralston second order explicit", n_stages, 2, false, false,
      {0.0, 2.0 / 3.0},
      {{0.0, 0.0},
       {2.0 / 3.0, 0.0}},
      {1.0 / 4.0, 3.0 / 4.0},
      {1.0 / 4.0, 3.0 / 4.0}};

//! Heun's method (2 E)
#elif RK_INTEGRATOR_TYPE == 7

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Heun second order explicit", n_stages, 2, false, false,
      {0.0, 1.0},
      {{0.0, 0.0},
       {1.0, 0.0}},
      {1.0 / 2.0, 1.0 / 2.0},
      {1.0 / 2.0, 1.0 / 2.0}};

//! Crank-Nicolson method (2 I)
#elif RK_INTEGRATOR_TYPE == 8

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Crank-Nicolson second order implicit", n_stages, 2, false, true,
      {0.0, 1.0},
      {{0.0, 0.0},
       {1.0 / 2.0, 1.0 / 2.0}},
      {1.0 / 2.0, 1.0 / 2.0},
      {1.0 / 2.0, 1.0 / 2.0}};

//! Gauss-Legendre (2 I)
#elif RK_INTEGRATOR_TYPE == 30

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Gauss-Legendre second order implicit", n_stages, 2, false, true,
      {1.0 / 2.0 - sqrtthr / 6.0, 1.0 / 2.0 + sqrtthr / 6.0},
      {{1.0 / 4.0, 1.0 / 4.0 - sqrtthr / 6.0},
       {1.0 / 4.0 + sqrtthr / 6.0, 1.0 / 4.0}},
      {1.0 / 2.0, 1.0 / 2.0},
      {1.0 / 2.0, 1.0 / 2.0}};

//! Heun-Euler method (2/1 E)
#elif RK_INTEGRATOR_TYPE == 9

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Heun-Euler second order adaptive explicit", n_stages, 2, true, false,
      {0.0, 1.0},
      {{0.0, 0.0},
       {1.0, 0.0}},
      {1.0, 0.0},
      {1.0 / 2.0, 1.0 / 2.0}};

//! Lobatto IIIA method (2/1 I)
#elif RK_INTEGRATOR_TYPE == 10

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Lobatto IIIA second order adaptive implicit", n_stages, 2, true, true,
      {0.0, 1.0},
      {{0.0, 0.0},
       {1.0 / 2.0, 1.0 / 2.0}},
      {1.0, 0.0},
      {1.0 / 2.0, 1.0 / 2.0}};

//! Lobatto IIIB method (2/1 I)
#elif RK_INTEGRATOR_TYPE == 11

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Lobatto IIIB second order adaptive implicit", n_stages, 2, true, true,
      {1.0 / 2.0, 1.0 / 2.0},
      {{1.0 / 2.0, 0.0},
       {1.0 / 2.0, 0.0}},
      {1.0, 0.0},
      {1.0 / 2.0, 1.0 / 2.0}};

//! Lobatto IIIC method (2/1 I)
#elif RK_INTEGRATOR_TYPE == 12

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Lobatto IIIC second order adaptive implicit", n_stages, 2, true, true,
      {0.0, 1.0},
      {{1.0 / 2.0, -1.0 / 2.0},
       {1.0 / 2.0, 1.0 / 2.0}},
      {1.0, 0.0},
      {1.0 / 2.0, 1.0 / 2.0}};

//! Fehlberg's method (2/1 E)
#elif RK_INTEGRATOR_TYPE == 13

const uint8_t n_stages = 3;
const ButcherTable <n_stages> RK_Table = {"Fehlberg second order adaptive explicit", n_stages, 2, true, false,
      {0.0, 1.0 / 2.0, 1.0},
      {{0.0, 0.0, 0.0},
       {1.0 / 2.0, 0.0, 0.0},
       {1.0 / 256.0, 255.0 / 256.0, 0.0}},
      {1.0 / 256.0, 255.0 / 256.0, 0.0},
      {1.0 / 512.0, 255.0 / 256.0, 1.0 / 512.0}};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Order 3
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Radau IA method (3 I)
#elif RK_INTEGRATOR_TYPE == 14

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Radau IA third order implicit", n_stages, 3, false, true,
      {0.0, 2.0 / 3.0},
      {{1.0 / 4.0, -1.0 / 4.0},
       {1.0 / 4.0, 5.0 / 12.0}},
      {1.0 / 4.0, 3.0 / 4.0},
      {1.0 / 4.0, 3.0 / 4.0}};

//! Radau IIA method (3 I)
#elif RK_INTEGRATOR_TYPE == 15

const uint8_t n_stages = 2;
const ButcherTable <n_stages> RK_Table = {"Radau IIA third order implicit", n_stages, 3, false, true,
      {1.0 / 3.0, 1.0},
      {{5.0 / 12.0, -1.0 / 12.0},
       {3.0 / 4.0, 1.0 / 4.0}},
      {3.0 / 4.0, 1.0 / 4.0},
      {3.0 / 4.0, 1.0 / 4.0}};

//! Kutta's third-order method (3 E)
#elif RK_INTEGRATOR_TYPE == 16

const uint8_t n_stages = 3;
const ButcherTable <n_stages> RK_Table = {"Kutta third order explicit", n_stages, 3, false, false,
      {0.0, 1.0 / 2.0, 1.0},
      {{0.0, 0.0, 0.0},
       {1.0 / 2.0, 0.0, 0.0},
       {-1.0, 2.0, 0.0}},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}};

//! Heun's third-order method (3 E)
#elif RK_INTEGRATOR_TYPE == 17

const uint8_t n_stages = 3;
const ButcherTable <n_stages> RK_Table = {"Heun third order explicit", n_stages, 3, false, false,
      {0.0, 1.0 / 3.0, 2.0 / 3.0},
      {{0.0, 0.0, 0.0},
       {1.0 / 3.0, 0.0, 0.0},
       {0.0, 2.0 / 3.0, 0.0}},
      {1.0 / 4.0, 0.0, 3.0 / 4.0},
      {1.0 / 4.0, 0.0, 3.0 / 4.0}};

//! Ralston's's third-order method (3 E)
#elif RK_INTEGRATOR_TYPE == 18

const uint8_t n_stages = 3;
const ButcherTable <n_stages> RK_Table = {"Ralston third order explicit", n_stages, 3, false, false,
      {0.0, 1.0 / 2.0, 3.0 / 4.0},
      {{0.0, 0.0, 0.0},
       {1.0 / 2.0, 0.0, 0.0},
       {0.0, 3.0 / 4.0, 0.0}},
      {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0},
      {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0}};

//! Strong Stability Preserving method (3 E)
#elif RK_INTEGRATOR_TYPE == 19

const uint8_t n_stages = 3;
const ButcherTable <n_stages> RK_Table = {"Strong Stability Preserving third order explicit", n_stages, 3, false, false,
      {0.0, 1.0, 1.0 / 2.0},
      {{0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {1.0 / 4.0, 1.0 / 4.0, 0.0}},
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},
      {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0}};

//! Implicit Runge-Kutta method (3 I)
#elif RK_INTEGRATOR_TYPE == 20

const uint8_t n_stages = 4;
const ButcherTable <n_stages> RK_Table = {"Runge-Kutta third order implicit", n_stages, 3, false, true,
      {1.0 / 2.0, 2.0 / 3.0, 1.0 / 2.0, 1.0},
      {{1.0 / 2.0, 0.0, 0.0, 0.0},
       {1.0 / 6.0, 1.0 / 2.0, 0.0, 0.0},
       {-1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0, 0.0},
       {3.0 / 2.0, -3.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0}},
      {3.0 / 2.0, -3.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0},
      {3.0 / 2.0, -3.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0}};


//! Bogacki–Shampine method (3/2 E)
#elif RK_INTEGRATOR_TYPE == 21

const uint8_t n_stages = 4;
const ButcherTable <n_stages> RK_Table = {"Bogacki–Shampine third order adaptive explicit", n_stages, 3, true, false,
      {0.0, 1.0 / 2.0, 3.0 / 4.0, 1.0},
      {{0.0, 0.0, 0.0, 0.0},
       {1.0 / 2.0, 0.0, 0.0, 0.0},
       {0.0, 3.0 / 4.0, 0.0, 0.0},
       {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0}},
      {7.0 / 24.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 8.0},
      {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0}};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Order 4
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Lobatto IIIA fourth-order method (4/3 I)
#elif RK_INTEGRATOR_TYPE == 22

const uint8_t n_stages = 3;
const ButcherTable <n_stages> RK_Table = {"Lobatto IIIA fourth order adaptive implicit", n_stages, 4, true, true,
      {0.0, 1.0 / 2.0, 1.0},
      {{0.0, 0.0, 0.0},
       {5.0 / 24.0, 1.0 / 3.0, -1.0 / 24.0},
       {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {-1.0 / 2.0, 2.0, -1.0 / 2.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}};

//! Lobatto IIIB fourth-order method (4/3 I)
#elif RK_INTEGRATOR_TYPE == 23

const uint8_t n_stages = 3;
const ButcherTable <n_stages> RK_Table = {"Lobatto IIIB fourth order adaptive implicit", n_stages, 4, true, true,
      {0.0, 1.0 / 2.0, 1.0},
      {{1.0 / 6.0, -1.0 / 6.0, 0.0},
       {1.0 / 6.0, 1.0 / 3.0, 0.0},
       {1.0 / 6.0, 5.0 / 6.0, 0.0}},
      {-1.0 / 2.0, 2.0, -1.0 / 2.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}};

//! Lobatto IIIC fourth-order method (4/3 I)
#elif RK_INTEGRATOR_TYPE == 24

const uint8_t n_stages = 3;
const ButcherTable <n_stages> RK_Table = {"Lobatto IIIC fourth order adaptive implicit", n_stages, 4, true, true,
      {0.0, 1.0 / 2.0, 1.0},
      {{1.0 / 6.0, -1.0 / 3.0, 1.0 / 6.0},
       {1.0 / 6.0, 5.0 / 12.0, -1.0 / 12.0},
       {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}},
      {-1.0 / 2.0, 2.0, -1.0 / 2.0},
      {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}};

//! Classic Runge-Kutta method (4 E)
#elif RK_INTEGRATOR_TYPE == 25

const uint8_t n_stages = 4;
const ButcherTable <n_stages> RK_Table = {"Runge-Kutta fourth order explicit", n_stages, 4, false, false,
      {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0},
      {{0.0, 0.0, 0.0, 0.0},
       {1.0 / 2.0, 0.0, 0.0, 0.0},
       {0.0, 1.0 / 2.0, 0.0, 0.0},
       {0.0, 0.0, 1.0, 0.0}},
      {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0},
      {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0}};

//! Kutta's 3/8 rule (4 E)
#elif RK_INTEGRATOR_TYPE == 26

const uint8_t n_stages = 4;
const ButcherTable <n_stages> RK_Table = {"Kutta 3/8 fourth order explicit", n_stages, 4, false, false,
      {0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0},
      {{0.0, 0.0, 0.0, 0.0},
       {1.0 / 3.0, 0.0, 0.0, 0.0},
       {-1.0 / 3.0, 1.0, 0.0, 0.0},
       {1.0, -1.0, 1.0, 0.0}},
      {1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0},
      {1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0}};

//----------------------------------------------------------------------------------------------------------------------------------------------------
// Order 5
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Runge-Kutta-Fehlberg method (5/4 E)
#elif RK_INTEGRATOR_TYPE == 27

const uint8_t n_stages = 6;
const ButcherTable <n_stages> RK_Table = {"Runge-Kutta-Fehlberg fifth order explicit", n_stages, 5, true, false,
      {0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0},
      {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0 / 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0, 0.0},
       {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0, 0.0},
       {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0, 0.0},
       {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0, 0.0}},
      {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0},
      {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0.0}};

//! Cash-Karp method (5/4 E)
#elif RK_INTEGRATOR_TYPE == 28

const uint8_t n_stages = 6;
const ButcherTable <n_stages> RK_Table = {"Cash-Karp fifth order explicit", n_stages, 5, true, false,
      {0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0},
      {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0},
       {3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0, 0.0, 0.0, 0.0},
       {-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0, 0.0, 0.0},
       {1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0, 0.0}},
      {37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0},
      {2825.0 / 27648.0, 0.0, 18575.0 / 48384.0, 13525.0 / 55296.0, 277.0 / 14336.0, 1.0 / 4.0}};

//! Dormand-Prince method (5/4 E)
#elif RK_INTEGRATOR_TYPE == 29

const uint8_t n_stages = 7;
const ButcherTable <n_stages> RK_Table = {"Dormand-Prince fifth order explicit", n_stages, 5, true, false,
      {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0},
      {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
       {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0},
       {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
       {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}},
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0},
      {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0}};

#else
#error Unsupported RK integrator type
#endif

};

#endif

#endif
