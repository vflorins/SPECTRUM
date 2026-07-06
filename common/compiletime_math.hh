/*!
\file compiletime_lists.hh
\brief Defines global compile-time constant options that are available, and some associated constant lists.
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_COMPILETIME_MATH_HH
#define SPECTRUM_COMPILETIME_MATH_HH

#include <limits>
#include <cmath>
#include <cstdint>
#include <bit>

namespace Spectrum {


/*!
\brief Computes the fabs function
\author Lucius Schoenbaum
\date 11/07/2025
\param[in] x The argument
\return \f$\abs(x)\f$
*/
template <typename T>
constexpr T cfabs(T x) noexcept {
   return x < 0 ? -x : x;
}


/*!
\brief Computes the square root as constexpr,
while allowing constexpr-valid functions to be defined
that have the performance degradation of constexpr evaluation
when forced to execute at runtime.
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 11/07/2025
\param[in] x The argument
\return \f$\sqrt{x}\f$
*/
SPECTRUM_DEVICE_FUNC inline constexpr double csqrt(double x) {
   if (std::is_constant_evaluated()) {
      if (x < 0) return std::numeric_limits<double>::quiet_NaN();
      if (x == 0) return 0;
      const int maxiter = 10;
      int niter = 0;
      double y = x;

// Initial reduction based on a bitwise representation of a "double"
      uint64_t y_bit = std::bit_cast<uint64_t>(y);
      uint64_t processed = (((((y_bit & 0x7FF0000000000000ULL) - 0x3FF0000000000000ULL) >> 1) + 0x3FF0000000000000ULL) & 0x7FF0000000000000ULL) | (y_bit & 0x000FFFFFFFFFFFFFULL);
      y = std::bit_cast<double>(processed);

// Henon's method (Newton iterations)
      while ((cfabs((y * y - x) / x) > sp_tiny) && (niter < maxiter)) {
         y = 0.5 * (y + x / y);
         niter++;
      };
      return y;
   } else {
      return std::sqrt(x);
   }

}



};

#endif
