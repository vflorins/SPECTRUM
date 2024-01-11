/*!
\file background_solarwind.hh
\brief Declares a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef _BACKGROUND_SOLARWIND_HH
#define _BACKGROUND_SOLARWIND_HH

#include "background_base.hh"

namespace Spectrum {

//! Method for computing derivatives (0: analytical, 1: Numerical)
#define SOLARWIND_DERIVATIVE_METHOD 1

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_solarwind = "BackgroundSolarWind";

/*!
\brief Plasma background calculator for a radially expanding solar wind
\author Vladimir Florinski

Parameters: (BackgroundBase), GeoVector Omega, double r_ref, double dmax_fraction
*/
class BackgroundSolarWind : public BackgroundBase {

protected:

//! Angular velocity vector of a rotating star (persistent)
   GeoVector Omega;

//! Reference radius (persistent)
   double r_ref;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

//! Local coordinate system tied to the rotation axis (persistent)
   GeoVector eprime[3];

//! Velocity magnitude (persistent)
   double ur0;

//! Radial magnetic field at "r_ref" (persistent)
   double Br0;

//! Angular frequency magnitude (persistent)
   double w0;

// TODO implement magnetic axis tilt

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

//! Compute the maximum distance per time step
   void EvaluateDmax(void) override;

public:

//! Default constructor
   BackgroundSolarWind(void);

//! Copy constructor
   BackgroundSolarWind(const BackgroundSolarWind& other);

//! Destructor
   ~BackgroundSolarWind() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundSolarWind);
};

};

#endif
