/*!
\file background_spherical_obstacle.hh
\brief Declares a magnetic field background around a spherical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SPHERICAL_OBSTACLE_HH
#define SPECTRUM_BACKGROUND_SPHERICAL_OBSTACLE_HH

#include "common/vectors.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSphericalObstacle class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Magnetic field around a spherical obstacle
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), double r_sphere, double dmax_fraction
*/
template <typename HConfig_>
class BackgroundSphericalObstacle {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundSphericalObstacle";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

   static constexpr double dmax0 = Config::dmax0;

//! Maximum fraction of the radial distance per step (persistent)
   static constexpr double dmax_fraction = Config::dmax_fraction;

   static constexpr GeoVector r0 = Config::r0;

//! Radius of spherical obstacle (persistent)
   static constexpr double r_ref = Config::r_ref;

   static constexpr GeoVector B0 = Config::B0;

//! Dipole moment (persistent)
   static constexpr GeoVector M = B0*Cube(r_ref);

public:

   //! Compute the maximum distance per time step
   template <typename Coordinates>
   static status_t EvaluateDmax(Coordinates&, double*);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

};

};

// Something like this is needed for templated classes
#include "background_spherical_obstacle.cc"

#endif
