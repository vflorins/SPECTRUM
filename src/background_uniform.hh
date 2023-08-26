/*!
\file background_uniform.hh
\brief Declares a simple uniform field background, mainly for testing
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_UNIFORM_HH
#define SPECTRUM_BACKGROUND_UNIFORM_HH

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

//! Readable name of the BackgroundUniform class
const std::string bg_name_uniform = "BackgroundUniform";

/*!
\brief Constant EM field, mainly for testing
\author Vladimir Florinski

Parameters: (BackgroundBase)
*/
class BackgroundUniform : public BackgroundBase {

protected:

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   void EvaluateBackground(void) override;

//! Compute the internal u, B, and E derivatives
   void EvaluateBackgroundDerivatives(void) override;

public:

//! Default constructor
   BackgroundUniform(void);

//! Copy constructor
   BackgroundUniform(const BackgroundUniform& other);

//! Destructor
   ~BackgroundUniform() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundUniform);
};

};

#endif
