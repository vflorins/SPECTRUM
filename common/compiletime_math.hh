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

namespace Spectrum {

constexpr double csqrt(double x) {
   if (x < 0) return std::numeric_limits<double>::quiet_NaN();
   if (x == 0) return 0;
   double a = x;
   for (int i = 0; i < 50; ++i) {
      a = 0.5*(a + x/a);
   }
   return a;
}





};

#endif
