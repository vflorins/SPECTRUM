/*!
\file background_shock.cc
\brief Implements a simple planar shock field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_shock.hh"
#include "common/print_warn.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundShock methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Juan G Alonso Guzman
\date 05/14/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundShock<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
// Upstream
   if ((coords.Pos() - r0) * n_shock - v_shock * coords.Time() > 0) {
      if constexpr (RequestedFields::Fluv_found()) fields.Fluv('w') = u0;
      if constexpr (RequestedFields::Mag_found()) fields.Mag('w') = B0;
      if constexpr (RequestedFields::Iv0_found()) fields.Iv0('w') = 1.0;
   }
// Downstream
   else {
      if constexpr (RequestedFields::Fluv_found()) fields.Fluv('w') = u1;
      if constexpr (RequestedFields::Mag_found()) fields.Mag('w') = B1;
      if constexpr (RequestedFields::Iv0_found()) fields.Iv0('w') = 2.0;
   };

   if constexpr (RequestedFields::Ele_found()) fields.Ele('w') = -(fields.Fluv() ^ fields.Mag()) / c_code;

   return 0;
};

/*!
\author Juan G Alonso Guzman
\date 10/14/2022
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundShock<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
// Spatial derivatives are zero
   if constexpr (RequestedFields::DelFluv_found()) fields.DelFluv('w') = gm_zeros;
   if constexpr (RequestedFields::DelMag_found()) fields.DelMag('w') = gm_zeros;
   if constexpr (RequestedFields::DelEle_found()) fields.DelEle('w') = gm_zeros;

// Time derivatives are zero
   if constexpr (RequestedFields::DotFluv_found()) fields.DotFluv('w') = gv_zeros;
   if constexpr (RequestedFields::DotMag_found()) fields.DotMag('w') = gv_zeros;
   if constexpr (RequestedFields::DotEle_found()) fields.DotEle('w') = gv_zeros;

   return 0;
};


/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename HConfig>
template <typename Coordinates>
status_t BackgroundShock<HConfig>::EvaluateDmax(Coordinates& coords, double* dmax)
{
   *dmax = dmax0;
   return 0;
};


};
