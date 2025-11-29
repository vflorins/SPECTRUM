/*!
\file background_uniform.hh
\brief Declares a simple uniform field background, mainly for testing
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_UNIFORM_HH
#define SPECTRUM_BACKGROUND_UNIFORM_HH

#include "utils_numerical_derivatives.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundUniform class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Constant EM field, mainly for testing
\author Vladimir Florinski

*/
template <typename HConfig_>
class BackgroundUniform {
public:

   //! Readable name of the class
   static constexpr std::string_view name = "BackgroundUniform";

public:

   using HConfig = HConfig_;
   using Cfg = HConfig::BackgroundConfig;

protected:

//! Electric field (persistent)
   static constexpr GeoVector E0 = -(Cfg::u0 ^ Cfg::B0) / c_code;

////! Set up the field evaluator based on "params"
//   void SetupBackground(bool construct);

public:

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

   template <typename Coordinates>
   static double EvaluateDmax(Coordinates& coords);

};

};

// Something like this is needed for templated classes
#include "background_uniform.cc"

#endif
