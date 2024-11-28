/*!
\file background_solarwind_termshock.hh
\brief Declares a plasma background class for the constant speed supersonic wind of a rotating star
\author Juan G Alonso Guzman
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef _BACKGROUND_SOLARWIND_TERMSHOCK_HH
#define _BACKGROUND_SOLARWIND_TERMSHOCK_HH

#include "background_solarwind.hh"

namespace Spectrum {

//! Integer exponent of decrease of solar wind speed beyond the termination shock
#define SOLARWIND_TERMSHOCK_SPEED_EXPONENT 2

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWindTermShock class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_solarwind_termshock = "BackgroundSolarWindTermShock";

/*!
\brief Plasma background calculator for a radially expanding solar wind with a spherical termination shock
\author Juan G Alonso Guzman

Parameters: (BackgroundSolarWind), double r_TS, double w_TS, double s_TS
*/
class BackgroundSolarWindTermShock : public BackgroundSolarWind {

protected:

//! Radius of termination shock (persistent)
   double r_TS;

//! Width of termination shock (persistent)
   double w_TS;

//! Strength of termination shock (persistent)
   double s_TS;

//! Inverse of s_TS (persistent)
   double s_TS_inv;

//! Maximum displacement in the shock region (persistent)
   double dmax_TS;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Modify radial flow (if necessary)
   void ModifyUr(const double r, double &ur_mod) override;

//! Get time lag for time dependent current sheet (if necessary)
   double TimeLag(const double r) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

//! Compute the maximum distance per time step
   void EvaluateDmax(void) override;

public:

//! Default constructor
   BackgroundSolarWindTermShock(void);

//! Copy constructor
   BackgroundSolarWindTermShock(const BackgroundSolarWindTermShock& other);

//! Destructor
   ~BackgroundSolarWindTermShock() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundSolarWindTermShock);
};

};

#endif
