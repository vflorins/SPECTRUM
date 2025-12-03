/*!
\file background_uniform.hh
\brief Declares a simple uniform field background, mainly for testing
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_UNIFORM_HH
#define SPECTRUM_BACKGROUND_UNIFORM_HH

#include "common/vectors.hh"

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
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

   static constexpr double dmax0 = Config::dmax0;

   static constexpr GeoVector u0 = Config::u0;

   static constexpr GeoVector B0 = Config::B0;

//! Electric field (persistent)
   static constexpr GeoVector E0 = -(u0 ^ B0) / c_code;

public:

   template <typename Coordinates>
   static status_t EvaluateDmax(Coordinates& coords, double*);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

};

};

// Something like this is needed for templated classes
#include "background_uniform.cc"

#endif
