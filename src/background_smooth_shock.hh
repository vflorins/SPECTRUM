/*!
\file background_smooth_shock.hh
\brief Declares a smooth MHD shock field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SMOOTH_SHOCK_HH
#define SPECTRUM_BACKGROUND_SMOOTH_SHOCK_HH

#include "background_shock.hh"

namespace Spectrum {

//! Flag to control smoothness of shock
#define SMOOTH_SHOCK_ORDER 4

//! Method for computing derivatives (0: analytical, 1: numerical)
#define SMOOTHSHOCK_DERIVATIVE_METHOD 0

//! Scaling factor to better match shock width when using smooth shock (tanh)
const double tanh_width_factor = 4.0;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSmoothShock class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BackgroundSmoothShock class
const std::string bg_name_smooth_shock = "BackgroundSmoothShock";

/*!
\brief Planar MHD shock with a smooth transition region
\author Juan G Alonso Guzman

Parameters: (BackgroundShock), double width_shock, double dmax_fraction
*/
class BackgroundSmoothShock : public BackgroundShock {

protected:

//! Width of shock transition region (persistent)
   double width_shock;

//! Fraction of the shock width to assign to dmax near shock (persistent)
   double dmax_fraction;

//! Distance from shock within which to modify dmax
   double dmax_limit;

//! Relative distance to shock (transient)
   double ds_shock;

//! Shock transition region function
   double ShockTransition(double x);

//! Derivative of shock transition region function
   double ShockTransitionDerivative(double x);

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
   BackgroundSmoothShock(void);

//! Copy constructor
   BackgroundSmoothShock(const BackgroundSmoothShock& other);

//! Destructor
   ~BackgroundSmoothShock() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundSmoothShock);
};

};

#endif
