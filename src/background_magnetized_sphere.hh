/*!
\file background_magnetized_sphere.hh
\brief Declares a background of a magnetized sphere
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_MAGNETIZED_SPHERE_HH
#define SPECTRUM_BACKGROUND_MAGNETIZED_SPHERE_HH

#include "background_spherical_obstacle.hh"

namespace Spectrum {

//! Method for computing derivatives of B (0: analytical, 1: numerical)
#define MAGNETIZED_SPHERE_DERIVATIVE_METHOD 0

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundMagnetizedSphere class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_magnetized_sphere = "BackgroundMagnetizedSphere";

/*!
\brief Magnetic field of a uniformly magnmetized sphere, inside and outside
\author Juan G Alonso Guzman

Parameters: (BackgroundSphericalObstacle)
*/
template <typename Fields_>
class BackgroundMagnetizedSphere : public BackgroundSphericalObstacle<Fields_> {
public:

   using Fields = Fields_;
   using BackgroundBase = BackgroundBase<Fields>;
   using BackgroundSphericalObstacle = BackgroundSphericalObstacle<Fields>;
   using BackgroundSphericalObstacle::_status;
   using BackgroundSphericalObstacle::_fields;
   using BackgroundSphericalObstacle::_ddata;
   using BackgroundSphericalObstacle::_pos;
   using BackgroundSphericalObstacle::container;
   using BackgroundSphericalObstacle::r0;
   using BackgroundSphericalObstacle::B0;
   using BackgroundSphericalObstacle::dmax0;
   // methods
   using BackgroundSphericalObstacle::EvaluateBmag;
   using BackgroundSphericalObstacle::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundSphericalObstacle::StopServerFront;
   using BackgroundSphericalObstacle::SetupBackground;
   using BackgroundSphericalObstacle::EvaluateBackground;
   using BackgroundSphericalObstacle::EvaluateBackgroundDerivatives;
   using BackgroundSphericalObstacle::NumericalDerivatives;

protected:

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

public:

//! Default constructor
   BackgroundMagnetizedSphere(void);

//! Copy constructor
   BackgroundMagnetizedSphere(const BackgroundMagnetizedSphere& other);

//! Destructor
   ~BackgroundMagnetizedSphere() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundMagnetizedSphere);
};

};

// Something like this is needed for templated classes
#include "background_magnetized_sphere.cc"

#endif
