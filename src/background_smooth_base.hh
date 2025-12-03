/*!
\file background_smooth_base.hh
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/
#ifndef SPECTRUM_BACKGROUND_SMOOTH_BASE_HH
#define SPECTRUM_BACKGROUND_SMOOTH_BASE_HH

#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSmoothDiscontinuity class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Smooth discontinuity base class
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

*/
template <typename HConfig_>
class BackgroundSmoothBase {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundSmoothBase";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

   static constexpr int smooth_discontinuity_order = Config::smooth_discontinuity_order;

   static constexpr double tanh_width_factor = Config::tanh_width_factor;

protected:

/*!
\author Juan G Alonso Guzman
\date 05/30/2024
\param [in] x Relative transition region location
\return Relative value of discontinuous quantity
*/
   double Transition(double x)
   {
      double y = x + 0.5;
      if constexpr (smooth_discontinuity_order == 0) {
         // not continous
         if (x < -0.5) return 0.0;
         else if (x > 0.5) return 1.0;
         else return y;
      }
      else if constexpr (smooth_discontinuity_order == 1){
         // differentiable
         if (x < -0.5) return 0.0;
         else if (x > 0.5) return 1.0;
         else return Sqr(y) * (3.0 - 2.0 * y);
      }
      else if constexpr (smooth_discontinuity_order == 2) {
         // twice differentiable
         if (x < -0.5) return 0.0;
         else if (x > 0.5) return 1.0;
         else return Cube(y) * (10.0 - 15.0 * y + 6.0 * Sqr(y));
      }
      else if constexpr (smooth_discontinuity_order == 3) {
         // thrice differentiable
         if (x < -0.5) return 0.0;
         else if (x > 0.5) return 1.0;
         else return Sqr(Sqr(y)) * (35.0 - 84.0 * y + 70.0 * Sqr(y) - 20.0 * Cube(y));
      }
      else {
         // smooth
         return 0.5 * (1.0 + tanh(tanh_width_factor * x));
      }
   };

/*!
\author Juan G Alonso Guzman
\date 05/30/2024
\param [in] x Relative transition region location
\return Derivative of relative value of discontinuous quantity
*/
   static double TransitionDerivative(double x)
   {
      double y = x + 0.5;
      if constexpr (smooth_discontinuity_order == 0) {
         if (x < -0.5) return 0.0;
         if (x > 0.5) return 0.0;
         return 1.0;
      }
      else if constexpr (smooth_discontinuity_order == 1){
         if (x < -0.5) return 0.0;
         if (x > 0.5) return 0.0;
         return 6.0 * y * (1.0 - y);
      }
      else if constexpr (smooth_discontinuity_order == 2) {
         if (x < -0.5) return 0.0;
         if (x > 0.5) return 0.0;
         return 30.0 * Sqr(y) * (1.0 - 2.0 * y + Sqr(y));
      }
      else if constexpr (smooth_discontinuity_order == 3) {
         if (x < -0.5) return 0.0;
         if (x > 0.5) return 0.0;
         return 140.0 * Cube(y) * (1.0 - 3.0 * y + 3.0 * Sqr(y) - 1.0 * Cube(y));
      }
      else {
         return 0.5 * tanh_width_factor * (1.0 - Sqr(tanh(tanh_width_factor * x)));
      }
   };


};

};






#endif
