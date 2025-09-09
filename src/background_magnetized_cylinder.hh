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
template <typename HyperParams_>
class BackgroundMagnetizedCylinder : public BackgroundCylindricalObstacle<HyperParams_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundMagnetizedCylinder";

public:

   using HyperParams = HyperParams_;
   using BackgroundBase = BackgroundBase<HyperParams>;
   using BackgroundCylindricalObstacle = BackgroundCylindricalObstacle<HyperParams>;
   using BackgroundBase::_status;
   using BackgroundBase::_fields;
   using BackgroundBase::_ddata;
   using BackgroundBase::_pos;
   using BackgroundBase::container;
   using BackgroundBase::r0;
   using BackgroundBase::B0;
   using BackgroundBase::dmax0;
   // methods
   using BackgroundBase::EvaluateBmag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;
//   using BackgroundBase::EvaluateBackground;
//   using BackgroundBase::EvaluateBackgroundDerivatives;
   using BackgroundBase::NumericalDerivatives;

protected:

//! Compute the internal u, B, and E fields
   template <typename Fields>
   void EvaluateBackground(Fields&);

//! Compute the internal derivatives of the fields
   template <typename Fields>
   void EvaluateBackgroundDerivatives(Fields&);

public:

//! Default constructor
   BackgroundMagnetizedCylinder(void);

//! Copy constructor
   BackgroundMagnetizedCylinder(const BackgroundMagnetizedCylinder& other);

//! Destructor
   ~BackgroundMagnetizedCylinder() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundMagnetizedCylinder);

};

};

// Something like this is needed for templated classes
#include "background_magnetized_cylinder.cc"

#endif
