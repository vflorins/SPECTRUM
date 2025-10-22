//
// Created by Lucius Schoenbaum on 9/9/25.
//

#ifndef SPECTRUM_RK_LISTS_HH
#define SPECTRUM_RK_LISTS_HH

#include <array>
#include <string_view>

#include "common/definitions.hh"

namespace Spectrum {

enum class RKIntegrator {
//! Euler method (1 E)
   Euler_1E,
//! Backward Euler method (1 I)
   Euler_1I,
//! Implicit midpoint method (2 I)
   Midpoint_2I,
//! Midpoint method (2 E)
   Midpoint_2E,
//! Kraaijevanger-Spijker method (2 I), not convex
   KraaijevangerSpijker_2I,
//! Qin-Zhang method (2 I)
   QinZhang_2I,
//! Ralston's method (2 E)
   Ralston_2E,
//! Heun's method (2 E)
   Heun_2E,
//! Crank-Nicolson method (2 I)
   CrankNicolson_2I,
//! Gauss-Legendre (2 I)
   GaussLegendre_2I,
//! Heun-Euler method (2/1 E)
   HeunEuler_21E,
//! Lobatto IIIA method (2/1 I)
   LobattoIIIA_21I,
//! Lobatto IIIB method (2/1 I)
   LobattoIIIB_21I,
//! Lobatto IIIC method (2/1 I)
   LobattoIIIC_21I,
//! Fehlberg's method (2/1 E)
   Fehlberg_21E,
//! Radau IA method (3 I)
   RadauIA_3I,
//! Radau IIA method (3 I)
   RadauIIA_3I,
//! Kutta's third-order method (3 E)
   Kutta_3E,
//! Heun's third-order method (3 E)
   Heun_3E,
//! Ralston's's third-order method (3 E)
   Ralston_3E,
//! Strong Stability Preserving method (3 E)
   StrongStability_3E,
//! Implicit Runge-Kutta method (3 I)
   RungeKutta_3I,
//! Bogacki–Shampine method (3/2 E)
   BogackiShampine_32E,
//! Lobatto IIIA fourth-order method (4/3 I)
   LobattoIIIA_43I,
//! Lobatto IIIB fourth-order method (4/3 I)
   LobattoIIIB_43I,
//! Lobatto IIIC fourth-order method (4/3 I)
   LobattoIIIC_43I,
//! Classic Runge-Kutta method (4 E)
   RungeKutta_4E,
//! Kutta's 3/8 rule (4 E)
   Kutta38_4E,
//! Runge-Kutta-Fehlberg method (5/4 E)
   RungeKuttaFehlberg_54E,
//! Cash-Karp method (5/4 E)
   CashKarp_54E,
//! Dormand-Prince method (5/4 E)
   DormandPrince_54E,
//! Runge-Kutta-Fehlberg method (6/5 E)
   RungeKuttaFehlberg_65E,
//! Runge-Kutta-Fehlberg method (7/6 E)
   RungeKuttaFehlberg_76E,
//! Runge-Kutta-Fehlberg method (8/7 E)
   RungeKuttaFehlberg_87E
};

struct ButcherTableData {
//! Readable name
   std::string_view name;
//! Number of stages
   uint8_t rk_stages;
//! Formal order of convergence
   uint8_t order;
//! Adaptive or not
   bool adaptive;
//! Implicit or not
   bool implicit;
};

/*!
Number of stages of RK method (look-up is via the enum value).
\author Lucius Schoenbaum
\date 09/09/2025
*/
const constexpr std::array<ButcherTableData, 33> RKData = {
      /* First Order */
      ButcherTableData{"Euler first order explicit", 1, 1, false, false,},
      ButcherTableData{"Backward Euler first order implicit", 1, 1, false, true,},
      /* Second Order */
      ButcherTableData{"Midpoint second order implicit", 1, 2, false, true,},
      ButcherTableData{"Midpoint second order explicit", 2, 2, false, false,},
      ButcherTableData{"Kraaijevanger-Spijker second order implicit", 2, 2, false, true,},
      ButcherTableData{"Qin-Zhang second order implicit", 2, 2, false, true,},
      ButcherTableData{"Ralston second order explicit", 2, 2, false, false,},
      ButcherTableData{"Heun second order explicit", 2, 2, false, false,},
      ButcherTableData{"Crank-Nicolson second order implicit", 2, 2, false, true,},
      ButcherTableData{"Gauss-Legendre second order implicit", 2, 2, false, true,},
      ButcherTableData{"Heun-Euler second order adaptive explicit", 2, 2, true, false,},
      ButcherTableData{"Lobatto IIIA second order adaptive implicit", 2, 2, true, true,},
      ButcherTableData{"Lobatto IIIB second order adaptive implicit", 2, 2, true, true,},
      ButcherTableData{"Lobatto IIIC second order adaptive implicit", 2, 2, true, true,},
      ButcherTableData{"Fehlberg second order adaptive explicit", 3, 2, true, false,},
      /* Third Order */
      ButcherTableData{"Radau IA third order implicit", 2, 3, false, true,},
      ButcherTableData{"Radau IIA third order implicit", 2, 3, false, true,},
      ButcherTableData{"Kutta third order explicit", 3, 3, false, false,},
      ButcherTableData{"Heun third order explicit", 3, 3, false, false,},
      ButcherTableData{"Ralston third order explicit", 3, 3, false, false,},
      ButcherTableData{"Strong Stability Preserving third order explicit", 3, 3, false, false,},
      ButcherTableData{"Runge-Kutta third order implicit", 4, 3, false, true,},
      ButcherTableData{"Bogacki–Shampine third order adaptive explicit", 4, 3, true, false,},
      /* Fourth Order */
      ButcherTableData{"Lobatto IIIA fourth order adaptive implicit", 3, 4, true, true,},
      ButcherTableData{"Lobatto IIIB fourth order adaptive implicit", 3, 4, true, true,},
      ButcherTableData{"Lobatto IIIC fourth order adaptive implicit", 3, 4, true, true,},
      ButcherTableData{"Runge-Kutta fourth order explicit", 4, 4, false, false,},
      ButcherTableData{"Kutta 3/8 fourth order explicit", 4, 4, false, false,},
      /* Fifth or Higher Order */
      ButcherTableData{"Runge-Kutta-Fehlberg fifth order explicit", 6, 5, true, false,},
      ButcherTableData{"Cash-Karp fifth order explicit", 6, 5, true, false,},
      ButcherTableData{"Dormand-Prince fifth order explicit", 7, 5, true, false,},
      ButcherTableData{"Runge-Kutta-Fehlberg sixth order explicit", 8, 6, true, false,},
      ButcherTableData{"Runge-Kutta-Fehlberg seventh order explicit", 10, 7, true, false,},
};


template <RKIntegrator rk_integrator>
struct ButcherTable
{
//! Data about the RK scheme that is dependent on the Butcher Table
   static constexpr ButcherTableData data = RKData[static_cast<size_t>(rk_integrator)];

//! Absolute tolerance (used only to avoid divission by zero in computing the relative tolerance)
   static constexpr double rk_tol_abs = 1.0E-9;

//! Relative tolerance
   static constexpr double rk_tol_rel = 1.0E-9;

//! Safety factor: the step size cannot change by more than this factor during one iteration
   static constexpr double rk_safety = 2.0;

//! When computing the new time step, use this factor (should be very close to 1) - this fine adjustment can affect performance
   static constexpr double rk_adjust = 0.995;

//! Time coefficients
   double a[data.rk_stages];

//! Slope coefficients
   double b[data.rk_stages][data.rk_stages];

//! Higher order weights
   double v[data.rk_stages];

//! Lower order weights
   double w[data.rk_stages];

   ButcherTable() {
      if constexpr (rk_integrator == RKIntegrator::Euler_1E) {
         a = {0.0};
         b = {0.0};
         v = {1.0};
         w = {1.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Euler_1I) {
         a = {1.0};
         b = {1.0};
         v = {1.0};
         w = {1.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Midpoint_2I) {
         a = {1.0 / 2.0};
         b = {1.0 / 2.0};
         v = {1.0};
         w = {1.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Midpoint_2E) {
         a = {0.0, 1.0 / 2.0};
         b = {{0.0, 0.0},
              {1.0 / 2.0, 0.0}};
         v = {0.0, 1.0};
         w = {0.0, 1.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::KraaijevangerSpijker_2I) {
         a = {1.0 / 2.0, 3.0 / 2.0};
         b = {{1.0 / 2.0, 0.0},
               {-1.0 / 2.0, 2.0}};
         v = {-1.0 / 2.0, 3.0 / 2.0};
         w = {-1.0 / 2.0, 3.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::QinZhang_2I) {
         a = {1.0 / 4.0, 3.0 / 4.0};
         b = {{1.0 / 4.0, 0.0},
            {1.0 / 2.0, 1.0 / 4.0}};
         v = {1.0 / 2.0, 1.0 / 2.0};
         w = {1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Ralston_2E) {
         a = {0.0, 2.0 / 3.0};
         b = {{0.0, 0.0},
            {2.0 / 3.0, 0.0}};
         v = {1.0 / 4.0, 3.0 / 4.0};
         w = {1.0 / 4.0, 3.0 / 4.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Heun_2E) {
         a = {0.0, 1.0};
         b = {{0.0, 0.0},
               {1.0, 0.0}};
         v = {1.0 / 2.0, 1.0 / 2.0};
         w = {1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::CrankNicolson_2I) {
         a = {0.0, 1.0};
         b = {{0.0, 0.0},
            {1.0 / 2.0, 1.0 / 2.0}};
         v = {1.0 / 2.0, 1.0 / 2.0};
         w = {1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::GaussLegendre_2I) {
         a = {1.0 / 2.0 - M_SQRT3 / 6.0, 1.0 / 2.0 + M_SQRT3 / 6.0};
         b = {{1.0 / 4.0, 1.0 / 4.0 - M_SQRT3 / 6.0},
            {1.0 / 4.0 + M_SQRT3 / 6.0, 1.0 / 4.0}};
         v = {1.0 / 2.0, 1.0 / 2.0};
         w = {1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::HeunEuler_21E) {
         a = {0.0, 1.0};
         b = {{0.0, 0.0},
            {1.0, 0.0}};
         v = {1.0, 0.0};
         w = {1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::LobattoIIIA_21I) {
         a = {0.0, 1.0};
         b = {{0.0, 0.0},
            {1.0 / 2.0, 1.0 / 2.0}};
         v = {1.0, 0.0};
         w = {1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::LobattoIIIB_21I) {
         a = {1.0 / 2.0, 1.0 / 2.0};
         b = {{1.0 / 2.0, 0.0},
            {1.0 / 2.0, 0.0}};
         v = {1.0, 0.0};
         w = {1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::LobattoIIIC_21I) {
         a = {0.0, 1.0};
         b = {{1.0 / 2.0, -1.0 / 2.0},
            {1.0 / 2.0, 1.0 / 2.0}};
         v = {1.0, 0.0};
         w = {1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Fehlberg_21E) {
         a = {0.0, 1.0 / 2.0, 1.0};
         b = {{0.0, 0.0, 0.0},
            {1.0 / 2.0, 0.0, 0.0},
            {1.0 / 256.0, 255.0 / 256.0, 0.0}};
         v = {1.0 / 256.0, 255.0 / 256.0, 0.0};
         w = {1.0 / 512.0, 255.0 / 256.0, 1.0 / 512.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::RadauIA_3I) {
         a = {0.0, 2.0 / 3.0};
         b = {{1.0 / 4.0, -1.0 / 4.0},
            {1.0 / 4.0, 5.0 / 12.0}};
         v = {1.0 / 4.0, 3.0 / 4.0};
         w = {1.0 / 4.0, 3.0 / 4.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::RadauIIA_3I) {
         a = {1.0 / 3.0, 1.0};
         b = {{5.0 / 12.0, -1.0 / 12.0},
            {3.0 / 4.0, 1.0 / 4.0}};
         v = {3.0 / 4.0, 1.0 / 4.0};
         w = {3.0 / 4.0, 1.0 / 4.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Kutta_3E) {
         a = {0.0, 1.0 / 2.0, 1.0};
         b = {{0.0, 0.0, 0.0},
            {1.0 / 2.0, 0.0, 0.0},
            {-1.0, 2.0, 0.0}};
         v = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
         w = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Heun_3E) {
         a = {0.0, 1.0 / 3.0, 2.0 / 3.0};
         b = {{0.0, 0.0, 0.0},
            {1.0 / 3.0, 0.0, 0.0},
            {0.0, 2.0 / 3.0, 0.0}};
         v = {1.0 / 4.0, 0.0, 3.0 / 4.0};
         w = {1.0 / 4.0, 0.0, 3.0 / 4.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Ralston_3E) {
         a = {0.0, 1.0 / 2.0, 3.0 / 4.0};
         b = {{0.0, 0.0, 0.0},
            {1.0 / 2.0, 0.0, 0.0},
            {0.0, 3.0 / 4.0, 0.0}};
         v = {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0};
         w = {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::StrongStability_3E) {
         a = {0.0, 1.0, 1.0 / 2.0};
         b = {{0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {1.0 / 4.0, 1.0 / 4.0, 0.0}};
         v = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
         w = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::RungeKutta_3I) {
         a = {1.0 / 2.0, 2.0 / 3.0, 1.0 / 2.0, 1.0};
         b = {{1.0 / 2.0, 0.0, 0.0, 0.0},
            {1.0 / 6.0, 1.0 / 2.0, 0.0, 0.0},
            {-1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0, 0.0},
            {3.0 / 2.0, -3.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0}};
         v = {3.0 / 2.0, -3.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0};
         w = {3.0 / 2.0, -3.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::BogackiShampine_32E) {
         a = {0.0, 1.0 / 2.0, 3.0 / 4.0, 1.0};
         b = {{0.0, 0.0, 0.0, 0.0},
            {1.0 / 2.0, 0.0, 0.0, 0.0},
            {0.0, 3.0 / 4.0, 0.0, 0.0},
            {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0}};
         v = {7.0 / 24.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 8.0};
         w = {2.0 / 9.0, 1.0 / 3.0, 4.0 / 9.0, 0.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::LobattoIIIA_43I) {
         a = {0.0, 1.0 / 2.0, 1.0};
         b = {{0.0, 0.0, 0.0},
            {5.0 / 24.0, 1.0 / 3.0, -1.0 / 24.0},
            {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}};
         v = {-1.0 / 2.0, 2.0, -1.0 / 2.0};
         w = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::LobattoIIIB_43I) {
         a = {0.0, 1.0 / 2.0, 1.0};
         b = {{1.0 / 6.0, -1.0 / 6.0, 0.0},
            {1.0 / 6.0, 1.0 / 3.0, 0.0},
            {1.0 / 6.0, 5.0 / 6.0, 0.0}};
         v = {-1.0 / 2.0, 2.0, -1.0 / 2.0};
         w = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::LobattoIIIC_43I) {
         a = {0.0, 1.0 / 2.0, 1.0};
         b = {{1.0 / 6.0, -1.0 / 3.0, 1.0 / 6.0},
            {1.0 / 6.0, 5.0 / 12.0, -1.0 / 12.0},
            {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0}};
         v = {-1.0 / 2.0, 2.0, -1.0 / 2.0};
         w = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::RungeKutta_4E) {
         a = {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
         b = {{0.0, 0.0, 0.0, 0.0},
            {1.0 / 2.0, 0.0, 0.0, 0.0},
            {0.0, 1.0 / 2.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0}};
         v = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
         w = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::Kutta38_4E) {
         a = {0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0};
         b = {{0.0, 0.0, 0.0, 0.0},
            {1.0 / 3.0, 0.0, 0.0, 0.0},
            {-1.0 / 3.0, 1.0, 0.0, 0.0},
            {1.0, -1.0, 1.0, 0.0}};
         v = {1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0};
         w = {1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::RungeKuttaFehlberg_54E) {
         a = {0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0};
         b = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 4.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {3.0 / 32.0, 9.0 / 32.0, 0.0, 0.0, 0.0, 0.0},
            {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0, 0.0, 0.0, 0.0},
            {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0, 0.0, 0.0},
            {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0, 0.0}};
         v = {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};
         w = {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::CashKarp_54E) {
         a = {0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0};
         b = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0},
            {3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0, 0.0, 0.0, 0.0},
            {-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0, 0.0, 0.0},
            {1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0, 0.0}};
         v = {37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0};
         w = {2825.0 / 27648.0, 0.0, 18575.0 / 48384.0, 13525.0 / 55296.0, 277.0 / 14336.0, 1.0 / 4.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::DormandPrince_54E) {
         a = {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0};
         b = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0},
            {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0},
            {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0},
            {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}};
         v = {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
         w = {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::RungeKuttaFehlberg_65E) {
         a = {0.0, 1.0 / 6.0, 4.0 / 15.0, 2.0 / 3.0, 4.0 / 5.0, 1.0, 0.0, 1.0};
         b = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {4.0 / 75.0, 16.0 / 75.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {5.0 / 6.0, -8.0 / 3.0, 5.0 / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {-8.0 / 5.0, 144.0 / 25.0, -4.0, 16.0 / 25.0, 0.0, 0.0, 0.0, 0.0},
            {361.0 / 320.0, -18.0 / 5.0, 407.0 / 128.0, -11.0 / 80.0, 55.0 / 128.0, 0.0, 0.0, 0.0},
            {-11.0 / 640.0, 0.0, 11.0 / 256.0, -11.0 / 160.0, 11.0 / 256.0, 0.0, 0.0, 0.0},
            {93.0 / 640.0, -18.0 / 5.0, 803.0 / 256.0, -11.0 / 160.0, 99.0 / 256.0, 0.0, 1.0, 0.0}};
         v = {7.0 / 1408.0, 0.0, 1125.0 / 2816.0, 9.0 / 32.0, 125.0 / 768.0, 0.0, 5.0 / 66.0, 5.0 / 66.0};
         w = {31.0 / 384.0, 0.0, 1125.0 / 2816.0, 9.0 / 32.0, 125.0 / 768.0, 5.0 / 66.0, 0.0, 0.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::RungeKuttaFehlberg_76E) {
         a = {0.0, 2.0 / 33.0, 4.0 / 33.0, 2.0 / 11.0, 1.0 / 2.0, 2.0 / 3.0, 6.0 / 7.0, 1.0, 0.0, 1.0};
         b = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {2.0 / 33.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 4.0 / 33.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 22.0, 0.0, 3.0 / 22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {43.0 / 64.0, 0.0, -165.0 / 64.0, 77.0 / 32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {-2383.0 / 486.0, 0.0 / 5.0, 1067.0 / 54.0, -26312.0 / 1701.0, 2176.0 / 1701.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {10077.0 / 4802.0, 0.0, -5643.0 / 686.0, 116259.0 / 16807.0, -6420.0 / 16807.0, 1053.0 / 2401.0, 0.0, 0.0, 0.0, 0.0},
            {-733.0 / 176.0, 0.0, 141.0 / 8.0, -335763.0 / 23296.0, 216.0 / 77.0, -4617.0 / 2816.0, 7203.0 / 9152.0, 0.0, 0.0, 0.0},
            {15.0 / 352.0, 0.0, 0.0, -5445.0 / 46592.0, 18.0 / 77.0, -1215.0 / 5632.0, 1029.0 / 18304.0, 0.0, 0.0, 0.0},
            {-1833.0 / 352.0, 0.0, 141.0 / 8.0, -51237.0 / 3584.0, 18.0 / 7.0, -729.0 / 512.0, 1029.0 / 1408.0, 0.0, 1.0, 0.0}};
         v = {11.0 / 864.0, 0.0, 0.0, 1771561.0 / 6289920.0, 32.0 / 105.0, 243.0 / 2560.0, 16807.0 / 74880.0, 0.0, 11.0 / 270.0, 11.0 / 270.0};
         w = {77.0 / 1440.0, 0.0, 0.0, 1771561.0 / 6289920.0, 32.0 / 105.0, 243.0 / 2560.0, 16807.0 / 74880.0, 11.0 / 270.0, 0.0, 0.0};
      }
      else if constexpr (rk_integrator == RKIntegrator::RungeKuttaFehlberg_87E) {
         a = {0.0, 2.0 / 27.0, 1.0 / 9.0, 1.0 / 6.0, 5.0 / 12.0, 1.0 / 2.0, 5.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0, 1.0, 0.0, 1.0};
         b = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {2.0 / 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 36.0, 1.0 / 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 24.0, 0.0, 1.0 / 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {5.0 / 12.0, 0.0, -25.0 / 16.0, 25.0 / 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 20.0, 0.0, 0.0, 1.0 / 4.0, 1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {-25.0 / 108.0, 0.0, 0.0, 125.0 / 108.0, -65.0 / 27.0, 125.0 / 54.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {31.0 / 300.0, 0.0, 0.0, 0.0, 61.0 / 225.0, -2.0 / 9.0, 13.0 / 900.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {2.0, 0.0, 0.0, -53.0 / 6.0, 704.0 / 45.0, -107.0 / 9.0, 67.0 / 9.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {-91.0 / 108.0, 0.0, 0.0, 23.0 / 108.0, -976.0 / 135.0, 311.0 / 54.0, -19.0 / 60.0, 17.0 / 6.0, -1.0 / 12.0, 0.0, 0.0, 0.0, 0.0},
            {2383.0 / 4100.0, 0.0, 0.0, -341.0 / 164.0, 4496.0 / 1025.0, -301.0 / 82.0, 2133.0 / 4100.0, 45.0 / 82.0, 45.0 / 164.0, 18.0 / 41.0, 0.0, 0.0, 0.0},
            {3.0 / 205.0, 0.0, 0.0, 0.0, 0.0, -6.0 / 41.0, -3.0 / 205.0, -3.0 / 41.0, 3.0 / 41.0, 6.0 / 41.0, 0.0, 0.0, 0.0},
            {-1777.0 / 4100.0, 0.0, 0.0, -341.0 / 164.0, 3396.0 / 1025.0, -289.0 / 82.0, 2193.0 / 4100.0, 51.0 / 82.0, 33.0 / 164.0, 12.0 / 41.0, 0.0, 1.0, 0.0}};
         v = {0.0, 0.0, 0.0, 0.0, 0.0, 34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0, 9.0 / 280.0, 9.0 / 280.0, 0.0, 41.0 / 840.0, 41.0 / 840.0};
         w = {41.0 / 840.0, 0.0, 0.0, 0.0, 0.0, 34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0, 9.0 / 280.0, 9.0 / 280.0, 41.0 / 840.0, 0.0, 0.0};
      }
      else {
         static_assert("Invalid RK Method");
      };
   }

};

};

#endif
