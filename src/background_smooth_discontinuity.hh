/*!
\file background_smooth_discontinuity.hh
\brief Declares a smooth MHD discontinuity field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SMOOTH_SHOCK_HH
#define SPECTRUM_BACKGROUND_SMOOTH_SHOCK_HH

#include "background_discontinuity.hh"

namespace Spectrum {

//! Flag to control smoothness of discontinuity
#define SMOOTH_DISCONT_ORDER 4

//! Method for computing derivatives (0: analytical, 1: numerical)
#define SMOOTHDISCONT_DERIVATIVE_METHOD 0

//! Scaling factor to better match discontinuity width when using smooth discontinuity (tanh)
const double tanh_width_factor = 4.0;

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSmoothDiscontinuity class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BackgroundSmoothDiscontinuity class
const std::string bg_name_smooth_discontinuity = "BackgroundSmoothDiscontinuity";

/*!
\brief Planar MHD discontinuity with a smooth transition region
\author Juan G Alonso Guzman

Parameters: (BackgroundShock), double width_discont, double dmax_fraction
*/
class BackgroundSmoothDiscontinuity : public BackgroundShock {

protected:

//! Width of discontinuity transition region (persistent)
   double width_discont;

//! Fraction of the discontinuity width to assign to dmax near discontinuity (persistent)
   double dmax_fraction;

//! Relative distance to discontinuity (transient)
   double ds_discont;

//! Shock transition region function
   double DiscontinuityTransition(double x);

//! Derivative of discontinuity transition region function
   double DiscontinuityTransitionDerivative(double x);

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
   BackgroundSmoothDiscontinuity(void);

//! Copy constructor
   BackgroundSmoothDiscontinuity(const BackgroundSmoothDiscontinuity& other);

//! Destructor
   ~BackgroundSmoothDiscontinuity() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundSmoothDiscontinuity);
};

};

#endif
