/*!
\file background_dipole.hh
\brief Declares a dipole magnetic field background without a flow
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_DIPOLE_HH
#define SPECTRUM_BACKGROUND_DIPOLE_HH

#include "background_base.hh"

namespace Spectrum {

//! Method for computing derivatives of B (0: analytical, 1: Numerical)
#define DIPOLE_DERIVATIVE_METHOD 0

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundDipole class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the class
const std::string bg_name_dipole = "BackgroundDipole";

/*!
\brief Magnetic field of a dipole
\author Vladimir Florinski

Parameters: (BackgroundBase), double r_ref, double dmax_fraction
*/
class BackgroundDipole : public BackgroundBase {

protected:

//! Dipole moment (persistent)
   GeoVector M;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

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
   BackgroundDipole(void);

//! Copy constructor
   BackgroundDipole(const BackgroundDipole& other);

//! Destructor
   ~BackgroundDipole() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundDipole);
};

};

#endif
