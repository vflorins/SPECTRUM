/*!
\file background_magnetized_cylinder.cc
\brief Implements a background of an infinite cylinder magnetized perpendicular to its axis
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#include "background_magnetized_cylinder.hh"

namespace Spectrum {

using namespace BackgroundOptions;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundMagnetizedCylinder methods
//----------------------------------------------------------------------------------------------------------------------------------------------------


/*!
\author Vladimir Florinski
\date 02/13/2024
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundMagnetizedCylinder<HConfig>::EvaluateBackground(Coordinates& coords, Fields& fields)
{
   BackgroundCylindricalObstacle::template EvaluateBackground<Coordinates, Fields, RequestedFields>(coords, fields);
   if constexpr (RequestedFields::Mag_found())
      fields.Mag('w') = B0 - fields.Mag();
   return 0;
};

/*!
\author Juan G Alonso Guzman
\date 03/11/2024
*/
template <typename HConfig>
template <typename Coordinates, typename Fields, typename RequestedFields>
status_t BackgroundMagnetizedCylinder<HConfig>::EvaluateBackgroundDerivatives(Coordinates& coords, Fields& fields)
{
   BackgroundCylindricalObstacle::template EvaluateBackgroundDerivatives<Coordinates, Fields, RequestedFields>(coords, fields);
   if constexpr (RequestedFields::DelMag_found()) fields.DelMag('w') *= -1.0;
   return 0;
};

/*!
\author Lucius Schoenbaum
\date 11/24/2025
*/
template <typename HConfig>
template <typename Coordinates>
status_t BackgroundMagnetizedCylinder<HConfig>::EvaluateDmax(Coordinates& coords, double& dmax)
{
   return BackgroundCylindricalObstacle::template EvaluateBackgroundDerivatives<Coordinates>(coords, dmax);
};


};
