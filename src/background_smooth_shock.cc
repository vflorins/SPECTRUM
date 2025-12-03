/*!
\file background_smooth_shock.cc
\brief Implements a smooth shock field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_smooth_shock.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSmoothShock methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Juan G Alonso Guzman
\date 05/14/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundSmoothShock<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   double a1, a2;
   auto ds_shock = ((coords.Pos() - r0) * n_shock - v_shock * coords.Time()) / width_shock;

   a1 = Transition(ds_shock);
   a2 = 1.0 - a1;

// TODO: Change the way Bvec transitions to depend directly on Uvec to keep some product constant
   if constexpr (RequestedFields::Fluv_found()) fields.Fluv('w') = u0 * a1 + u1 * a2;
   if constexpr (RequestedFields::Mag_found()) fields.Mag('w') = B0 * a1 + B1 * a2;
   if constexpr (RequestedFields::Ele_found()) fields.Ele('w') = -(fields.Fluv() ^ fields.Mag()) / c_code;
   if constexpr (RequestedFields::Iv0_found()) fields.Iv0('w') = 1.0 * a1 + 2.0 * a2; // same as 2.0 - a1

   return 0;
};

/*!
\author Juan G Alonso Guzman
\date 10/20/2023
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundSmoothShock<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   auto ds_shock = ((coords.Pos() - r0) * n_shock - v_shock * coords.Time()) / width_shock;
   if constexpr (RequestedFields::DelFluv_found()) {
      fields.DelFluv('w').Dyadic(n_shock, u0 - u1);
      fields.DelFluv('w') *= TransitionDerivative(ds_shock) / width_shock;
   };
   if constexpr (RequestedFields::DelMag_found()) {
      fields.DelMag('w').Dyadic(n_shock, B0 - B1);
      fields.DelMag('w') *= TransitionDerivative(ds_shock) / width_shock;
   };
   if constexpr (RequestedFields::DelAbsMag_found()) {
      fields.DelAbsMag('w') = fields.DelMag() * fields.HatMag();
   };
   if constexpr (RequestedFields::DelEle_found()) {
      fields.DelEle('w') = -((fields.DelFluv() ^ fields.Mag()) + (fields.Fluv() ^ fields.DelMag())) / c_code;
   };
   if constexpr (RequestedFields::DotFluv_found()) {
      fields.DotFluv('w') = (TransitionDerivative(ds_shock) * v_shock / width_shock) * (u1 - u0);
   };
   if constexpr (RequestedFields::DotMag_found()) {
      fields.DotMag('w') = (TransitionDerivative(ds_shock) * v_shock / width_shock) * (B1 - B0);
   };
   if constexpr (RequestedFields::DotEle_found()) {
      fields.DotEle('w') = -((fields.DotFluv() ^ fields.Mag()) + (fields.Fluv() ^ fields.DotMag())) / c_code;
   };
   return 0;
};

/*!
\author Juan G Alonso Guzman
\date 02/28/2025
*/
template <typename HConfig>
template <typename Coordinates>
status_t BackgroundSmoothShock<HConfig>::EvaluateDmax(Coordinates& coords, double* dmax)
{
   auto ds_shock = ((coords.Pos() - r0) * n_shock - v_shock * coords.Time()) / width_shock;
   *dmax = fmin(dmax_fraction * width_shock * fmax(1.0, fabs(ds_shock)), dmax0);
   return 0;
};

};
