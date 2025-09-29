/*!
\file background_uniform.hh
\brief Declares a simple uniform field background, mainly for testing
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_UNIFORM_HH
#define SPECTRUM_BACKGROUND_UNIFORM_HH

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Constant EM field, mainly for testing
\author Vladimir Florinski

Parameters: (BackgroundBase)
*/
template <typename HConfig_>
class BackgroundUniform : public BackgroundBase<HConfig_> {
private:

   //! Readable name of the BackgroundUniform class
   static constexpr std::string bg_name = "BackgroundUniform";

public:

   using HConfig = HConfig_;
   using BackgroundCoordinates = HConfig::BackgroundCoordinates;
   using BackgroundBase = BackgroundBase<HConfig>;
   using BackgroundBase::_status;
   using BackgroundBase::container;
   using BackgroundBase::_ddata;
   using BackgroundBase::dmax0;
   using BackgroundBase::r0;
   using BackgroundBase::u0;
   using BackgroundBase::B0;
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

//! Electric field (persistent)
   GeoVector E0;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the internal u, B, and E fields
   template <typename Fields>
   void EvaluateBackground(BackgroundCoordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Fields>
   void EvaluateBackgroundDerivatives(BackgroundCoordinates&, Fields&);

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

// Something like this is needed for templated classes
#include "background_uniform.cc"

#endif
