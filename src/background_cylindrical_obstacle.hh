/*!
\file background_cylindrical_obstacle.hh
\brief Declares a magnetic field background around a cylindrical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_CYLINDRICAL_OBSTACLE_HH
#define SPECTRUM_BACKGROUND_CYLINDRICAL_OBSTACLE_HH

#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCylindricalObstacle class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Magnetic field around a cylindrical obstacle
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

Parameters: (BackgroundBase), GeoVector axis, double r_obstacle, double dmax_fraction
*/
template <typename HConfig_>
class BackgroundCylindricalObstacle {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundCylindricalObstacle";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

protected:

   static constexpr GeoVector r0 = Config::r0;

   //! Axis of the cylinder (persistent)
   static constexpr GeoVector axis = Config::axis.Normalize();

   static constexpr GeoVector B0 = Config::B0.SubtractParallel(axis);

   static constexpr double dmax0 = Config::dmax0;

//! Maximum fraction of the radial distance per step (persistent)
   static constexpr double dmax_fraction = Config::dmax_fraction;

//! Radius of cylindrical obstacle (persistent)
   static constexpr double r_cylinder = Config::r_cylinder;

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
#include "background_cylindrical_obstacle.cc"

#endif
