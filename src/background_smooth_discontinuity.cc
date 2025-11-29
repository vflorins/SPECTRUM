/*!
\file background_smooth_discontinuity.cc
\brief Implements a smooth discontinuity field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_smooth_discontinuity.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSmoothDiscontinuity methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\date 05/30/2024
\param [in] x Relative transition region location
\return Relative value of discontinuous quantity
*/
template <typename HConfig>
constexpr double BackgroundSmoothDiscontinuity<HConfig>::DiscontinuityTransition(double x)
{
   double y = x + 0.5;
   if constexpr (Config::smooth_discontinuity_order == 0) {
      // not continous
      if (x < -0.5) return 0.0;
      else if (x > 0.5) return 1.0;
      else return y;
   }
   else if constexpr (Config::smooth_discontinuity_order == 1){
      // differentiable
      if (x < -0.5) return 0.0;
      else if (x > 0.5) return 1.0;
      else return Sqr(y) * (3.0 - 2.0 * y);
   }
   else if constexpr (Config::smooth_discontinuity_order == 2) {
      // twice differentiable
      if (x < -0.5) return 0.0;
      else if (x > 0.5) return 1.0;
      else return Cube(y) * (10.0 - 15.0 * y + 6.0 * Sqr(y));
   }
   else if constexpr (Config::smooth_discontinuity_order == 3) {
      // thrice differentiable
      if (x < -0.5) return 0.0;
      else if (x > 0.5) return 1.0;
      else return Sqr(Sqr(y)) * (35.0 - 84.0 * y + 70.0 * Sqr(y) - 20.0 * Cube(y));
   }
   else {
      // smooth
      return 0.5 * (1.0 + tanh(Config::tanh_width_factor * x));
   }
};

/*!
\author Juan G Alonso Guzman
\date 05/30/2024
\param [in] x Relative transition region location 
\return Derivative of relative value of discontinuous quantity
*/
template <typename HConfig>
constexpr double BackgroundSmoothDiscontinuity<HConfig>::DiscontinuityTransitionDerivative(double x)
{
   double y = x + 0.5;
   if constexpr (Config::smooth_discontinuity_order == 0) {
      if (x < -0.5) return 0.0;
      if (x > 0.5) return 0.0;
      return 1.0;
   }
   else if constexpr (Config::smooth_discontinuity_order == 1){
      if (x < -0.5) return 0.0;
      if (x > 0.5) return 0.0;
      return 6.0 * y * (1.0 - y);
   }
   else if constexpr (Config::smooth_discontinuity_order == 2) {
      if (x < -0.5) return 0.0;
      if (x > 0.5) return 0.0;
      return 30.0 * Sqr(y) * (1.0 - 2.0 * y + Sqr(y));
   }
   else if constexpr (Config::smooth_discontinuity_order == 3) {
      if (x < -0.5) return 0.0;
      if (x > 0.5) return 0.0;
      return 140.0 * Cube(y) * (1.0 - 3.0 * y + 3.0 * Sqr(y) - 1.0 * Cube(y));
   }
   else {
      return 0.5 * Config::tanh_width_factor * (1.0 - Sqr(tanh(Config::tanh_width_factor * x)));
   }
};

/*!
\author Juan G Alonso Guzman
\date 05/14/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundSmoothDiscontinuity<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   double a1, a2;
   auto ds_discont = ((coords.Pos() - r0) * n_discont - v_discont * coords.Time()) / width_discont;

   a1 = DiscontinuityTransition(ds_discont);
   a2 = 1.0 - a1;

   if constexpr (RequestedFields::Fluv_found()) fields.Fluv('w') = u0 * a1 + u1 * a2;
   if constexpr (RequestedFields::Mag_found()) fields.Mag('w') = B0 * a1 + B1 * a2;
   if constexpr (RequestedFields::Elc_found()) fields.Elc('w') = -(fields.Fluv() ^ fields.Mag()) / c_code;
   if constexpr (RequestedFields::Iv0_found()) fields.Iv0('w') = 1.0 * a1 + 2.0 * a2; // same as 2.0 - a1

   return 0;
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundSmoothDiscontinuity<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   auto ds_discont = ((coords.Pos() - r0) * n_discont - v_discont * coords.Time()) / width_discont;
   if constexpr (RequestedFields::DelFluv_found()) {
      fields.DelFluv('w').Dyadic(n_discont, u0 - u1);
      fields.DelFluv('w') *= DiscontinuityTransitionDerivative(ds_discont) / width_discont;
   };
   if constexpr (RequestedFields::DelMag_found()) {
      fields.DelMag('w').Dyadic(n_discont, B0 - B1);
      fields.DelMag('w') *= DiscontinuityTransitionDerivative(ds_discont) / width_discont;
   }
   if constexpr (RequestedFields::DelAbsMag_found()) {
      fields.DelAbsMag('w') = fields.DelMag() * fields.HatMag();
   };
   if constexpr (RequestedFields::DelElc_found()) {
      fields.DelElc('w') = -((fields.DelFluv() ^ fields.Mag()) + (fields.Fluv() ^ fields.DelMag())) / c_code;
   };
   if constexpr (RequestedFields::DotFluv_found()) {
      fields.DotFluv('w') = (DiscontinuityTransitionDerivative(ds_discont) * v_discont / width_discont) * (u1 - u0);
   };
   if constexpr (RequestedFields::DotMag_found()) {
      fields.DotMag('w') = (DiscontinuityTransitionDerivative(ds_discont) * v_discont / width_discont) * (B1 - B0);
   };
   if constexpr (RequestedFields::DotElc_found()) {
      fields.DotElc('w') = -((fields.DotFluv() ^ fields.Mag()) + (fields.Fluv() ^ fields.DotMag())) / c_code;
   };
   return 0;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
*/
template <typename HConfig>
template <typename Coordinates>
status_t BackgroundSmoothDiscontinuity<HConfig>::EvaluateDmax(Coordinates& coords, double& dmax)
{
   auto ds_discont = ((coords.Pos() - r0) * n_discont - v_discont * coords.Time()) / width_discont;
   dmax = fmin(dmax_fraction * width_discont * fmax(1.0, fabs(ds_discont)), dmax0);
   return 0;
};

};
