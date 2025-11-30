/*!
\file background_dipole.cc
\brief Implements a dipole magnetic field background without a flow
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_dipole.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDipole methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Vladimir Florinski
\date 03/25/2022
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundDipole<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   if constexpr (RequestedFields::Fluv_found()) {
      fields.Fluv('w') = gv_zeros;
   }
   if constexpr (RequestedFields::Mag_found()) {
      GeoVector posprime = coords.Pos() - r0;
      double r2 = posprime.Norm();
      double r5 = Cube(r2);
      r2 *= r2;
      r5 *= r2;
      fields.Mag('w') = (3.0 * (posprime * M) * posprime - r2 * M) / r5;
   }
   if constexpr (RequestedFields::Elc_found()) {
      fields.Elc('w') = gv_zeros;
   }
   if constexpr (RequestedFields::Iv0_found()) {
      fields.Iv0('w') = 1.0;
   }
   return 0;
};

/*!
\author Vladimir Florinski
\author Juan G Alonso Guzman
\date 10/14/2022
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundDipole<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   if constexpr (RequestedFields::DelFluv_found())
      fields.DelFluv('w') = gm_zeros;
   if constexpr (RequestedFields::DelMag_found() || RequestedFields::DelAbsMag_found()) {
      GeoVector posprime = coords.Pos() - r0;
      double r2 = posprime.Norm();
      double r5 = Cube(r2);
      r2 *= r2;
      r5 *= r2;
      GeoMatrix mr, rm, rr;
      mr.Dyadic(M,posprime);
      rm.Dyadic(posprime,M);
      rr.Dyadic(posprime);

      auto DelMag = 3.0 * (mr + rm + (M * posprime) * (gm_unit - 5.0 * rr / r2)) / r5;

      if constexpr (RequestedFields::DelMag_found())
         fields.DelMag('w') = DelMag;
      if constexpr (RequestedFields::DelAbsMag_found()) {
         fields.DelAbsMag('w') = DelMag * fields.HatMag();
      }
   };
   if constexpr (RequestedFields::DelElc_found())
      fields.DelElc('w') = gm_zeros;
   if constexpr (RequestedFields::DotFluv_found())
      fields.DotFluv('w') = gv_zeros;
   if constexpr (RequestedFields::DotMag_found())
      fields.DotMag('w') = gv_zeros;
   if constexpr (RequestedFields::DotAbsMag_found())
      fields.DotAbsMag('w') = 0.0;
   if constexpr (RequestedFields::DotElc_found())
      fields.DotElc('w') = gv_zeros;
   return 0;
};

/*!
\author Vladimir Florinski
\author Lucius Schoenbaum
\date 09/22/2025
*/
template <typename HConfig>
template <typename Coordinates>
status_t BackgroundDipole<HConfig>::EvaluateDmax(Coordinates& coords, double& dmax)
{
   dmax = fmin(dmax_fraction * (coords.Pos() - r0).Norm(), dmax0);
   return 0;
};

};
