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
class BackgroundMagnetizedCylinder : public BackgroundCylindricalObstacle<HConfig_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundMagnetizedCylinder";

public:

   using HConfig = HConfig_;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundMagnetizedCylinder<HConfig>>, typename HConfig::BackgroundConfig>;
   using BackgroundCoordinates = BackgroundConfig::Coordinates;
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
   using BackgroundCylindricalObstacle = BackgroundCylindricalObstacle<HConfig>;

   using BackgroundConfig::derivative_method;

protected:

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

//! Default constructor
   BackgroundMagnetizedCylinder(void);

//! Copy constructor
   BackgroundMagnetizedCylinder(const BackgroundMagnetizedCylinder& other);

//! Destructor
   ~BackgroundMagnetizedCylinder() = default;

//! Clone function
   CloneFunctionBackground(BackgroundMagnetizedCylinder);

};

};

// Something like this is needed for templated classes
#include "background_magnetized_cylinder.cc"

#endif
