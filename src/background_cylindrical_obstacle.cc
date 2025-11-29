/*!
\file background_cylindrical_obstacle.cc
\brief Implements a magnetic field background around a cylindrical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_cylindrical_obstacle.hh"
#include <common/matrix.hh>

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCylindricalObstacle methods
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
Compute the internal u, B, and E fields
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundCylindricalObstacle<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{

   GeoVector posprime = coords.Pos() - r0;
   posprime.SubtractParallel(axis);
   double posprimenorm = posprime.Norm();

   if constexpr (RequestedFields::Fluv_found())
      fields.Fluv('w') = gv_zeros;
   if constexpr (RequestedFields::Mag_found()) {
      if (posprimenorm < r_cylinder)
         fields.Mag('w') = gv_zeros;
      else {
         double s2 = posprime.Norm2();
         fields.Mag('w') = B0 - Sqr(r_cylinder) / s2 * (2.0 * (posprime * B0) / s2 * posprime - B0);
      };
   };
   if constexpr (RequestedFields::Elc_found())
      fields.Elc('w') = gv_zeros;
   if constexpr (RequestedFields::Iv1_found()) {
      if (posprimenorm < r_cylinder) fields.Iv1() = 0.0;
      else fields.Iv1('w') = 1.0;
   }

   return 0;
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 09/08/2025
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundCylindricalObstacle<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   GeoVector posprime = coords.Pos() - r0;
   posprime.SubtractParallel(axis);
   double posprimenorm = posprime.Norm();

   if constexpr (RequestedFields::DelFluv_found())
      fields.DelFluv('w') = gm_zeros;
   if constexpr (RequestedFields::DelMag_found()) {
      if (posprimenorm < r_cylinder)
         fields.DelMag('w') = gm_zeros;
      else {
// TODO: complete
      };
   };
   if constexpr (RequestedFields::DelElc_found()) fields.DelElc('w') = gm_zeros;
   if constexpr (RequestedFields::DotFluv_found()) fields.DotFluv('w') = gv_zeros;
   if constexpr (RequestedFields::DotMag_found()) fields.DotMag('w') = gv_zeros;
   if constexpr (RequestedFields::DotElc_found()) fields.DotElc('w') = gv_zeros;
   return 0;
};

/*!
\author Juan G Alonso Guzman
\author Lucius Schoenbaum
\date 08/06/2025
*/
template <typename HConfig>
template <typename Coordinates>
status_t BackgroundCylindricalObstacle<HConfig>::EvaluateDmax(Coordinates& coords, double& dmax)
{
   GeoVector posprime = coords.Pos() - r0;
   posprime.SubtractParallel(axis);
   dmax = fmin(dmax_fraction * posprime.Norm(), dmax0);
   return 0;
};


};

