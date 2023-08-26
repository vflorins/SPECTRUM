/*!
\file background_shock.hh
\brief Declares a simple shock field background
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SHOCK_HH
#define SPECTRUM_BACKGROUND_SHOCK_HH

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundShock class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BackgroundShock class
const std::string bg_name_shock = "BackgroundShock";

/*!
\brief Constant EM field, mainly for testing
\author Juan G Alonso Guzman

Parameters: (BackgroundBase)
*/
class BackgroundShock : public BackgroundBase {

protected:

//! Shock starting position (persistent)
   GeoVector r0_shock;

//! Shock normal (persistent)
   GeoVector n_shock;

//! Shock velocity (persistent)
   double v_shock;

//! Downstream flow vector (persistent), "u0" is upstream flow vector
   GeoVector u1;

//! Downstream magnetic field (persistent), "B0" is upstream magnetic field
   GeoVector B1;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

public:

//! Default constructor
   BackgroundShock(void);

//! Copy constructor
   BackgroundShock(const BackgroundShock& other);

//! Clone function
   CloneFunctionBackground(BackgroundShock);
};

};

#endif
