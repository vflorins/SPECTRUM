/*!
\file background_magnetized_cylinder.hh
\brief Declares a background of an infinite cylinder magnetized perpendicular to its axis
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_MAGNETIZED_CYLINDER_HH
#define SPECTRUM_BACKGROUND_MAGNETIZED_CYLINDER_HH

#include "background_cylindrical_obstacle.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundMagnetizedCylinder class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Magnetic field of a uniformly magnmetized cylinder, inside and outside
\author Vladimir Florinski
\author Juan G Alonso Guzman

Parameters: (BackgroundCylindricalObstacle)
*/
template <typename HConfig_>
class BackgroundMagnetizedCylinder: public BackgroundCylindricalObstacle<HConfig_> {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundMagnetizedCylinder";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;
   using BackgroundCylindricalObstacle = BackgroundCylindricalObstacle<HConfig>;
   using BackgroundCylindricalObstacle::B0;

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

public:

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

//! Compute the maximum distance per time step
   template <typename Coordinates>
   static status_t EvaluateDmax(Coordinates&, double*);

};

};

// Something like this is needed for templated classes
#include "background_magnetized_cylinder.cc"

#endif
