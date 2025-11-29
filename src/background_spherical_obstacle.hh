/*!
\file background_spherical_obstacle.hh
\brief Declares a magnetic field background around a spherical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SPHERICAL_OBSTACLE_HH
#define SPECTRUM_BACKGROUND_SPHERICAL_OBSTACLE_HH

#include "utils_numerical_derivatives.hh"

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
class BackgroundSphericalObstacle : public BackgroundBase<HConfig_>{
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundSphericalObstacle";

public:

   using HConfig = HConfig_;
   using BackgroundConfig = HConfig::BackgroundConfig;

//   using BackgroundBase = BackgroundBase<HConfig>;
//   using BackgroundBase::_status;
//   using BackgroundBase::container;
//   using BackgroundBase::_ddata;
//   using BackgroundBase::dmax0;
//   using BackgroundBase::r0;
//   using BackgroundBase::u0;
//   using BackgroundBase::B0;
//   // methods
//   using BackgroundBase::EvaluateDmax;
//   using BackgroundBase::GetDmax;
//   using BackgroundBase::StopServerFront;
//   using BackgroundBase::SetupBackground;

   using BackgroundConfig::derivative_method;

protected:

////! Radius of spherical obstacle (persistent)
//   static double r_sphere;
//
////! Dipole moment (persistent)
//   static GeoVector M;
//
////! Maximum fraction of the radial distance per step (persistent)
//   static double dmax_fraction;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct);

   //! Compute the maximum distance per time step
   template <typename Coordinates>
   static double EvaluateDmax(Coordinates&);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

//! Default constructor
   BackgroundSphericalObstacle(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundSphericalObstacle(const std::string_view& name_in, status_t status_in);

//! Copy constructor
   BackgroundSphericalObstacle(const BackgroundSphericalObstacle& other);

//! Destructor
   ~BackgroundSphericalObstacle() = default;

};

};

// Something like this is needed for templated classes
#include "background_spherical_obstacle.cc"

#endif
