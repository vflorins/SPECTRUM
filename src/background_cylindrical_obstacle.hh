/*!
\file background_cylindrical_obstacle.hh
\brief Declares a magnetic field background around a cylindrical object without a flow
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_CYLINDRICAL_OBSTACLE_HH
#define SPECTRUM_BACKGROUND_CYLINDRICAL_OBSTACLE_HH

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundCylindricalObstacle class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Magnetic field around a cylindrical obstacle
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), GeoVector axis, double r_obstacle, double dmax_fraction
*/
template <typename HConfig_>
class BackgroundCylindricalObstacle : public BackgroundBase<HConfig_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundCylindricalObstacle";

public:

   using HConfig = HConfig_;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundCylindricalObstacle<HConfig>>, typename HConfig::BackgroundConfig>;
   using BackgroundBase = BackgroundBase<HConfig>;
   using BackgroundBase::_status;
   using BackgroundBase::container;
   using BackgroundBase::_ddata;
   using BackgroundBase::dmax0;
   using BackgroundBase::r0;
   using BackgroundBase::u0;
   using BackgroundBase::B0;
   // methods
   using BackgroundBase::EvaluateAbsMag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;

   using BackgroundConfig::derivative_method;

protected:

//! Axis of the cylinder (persistent)
   GeoVector axis;

//! Radius of cylindrical obstacle (persistent)
   double r_cylinder;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct);

//! Compute the maximum distance per time step
   template <typename Coordinates>
   void EvaluateDmax(Coordinates&);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

//! Default constructor
   BackgroundCylindricalObstacle(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundCylindricalObstacle(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   BackgroundCylindricalObstacle(const BackgroundCylindricalObstacle& other);

//! Destructor
   ~BackgroundCylindricalObstacle() = default;

//! Clone function
   CloneFunctionBackground(BackgroundCylindricalObstacle);

};

};

// Something like this is needed for templated classes
#include "background_cylindrical_obstacle.cc"

#endif
