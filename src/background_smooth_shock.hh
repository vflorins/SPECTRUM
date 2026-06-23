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

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSmoothShock class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Planar MHD shock with a smooth transition region
\author Juan G Alonso Guzman

Parameters: (BackgroundShock), double width_shock, double dmax_fraction
*/
template <typename HConfig_>
class BackgroundSmoothShock: public BackgroundShock<HConfig_>, BackgroundSmoothBase<HConfig_> {
public:

//! Readable name of the class
   const std::string name = "BackgroundSmoothShock";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

   using BackgroundSmoothBase = BackgroundSmoothBase<HConfig>;
   using BackgroundSmoothBase::Transition;
   using BackgroundSmoothBase::TransitionDerivative;

   using BackgroundShock = BackgroundShock<HConfig>;
   using BackgroundShock::dmax0;
   using BackgroundShock::r0;
   using BackgroundShock::u0;
   using BackgroundShock::B0;
   using BackgroundShock::u1;
   using BackgroundShock::B1;
   using BackgroundShock::n_shock;
   using BackgroundShock::v_shock;

   //! Width of shock transition region (persistent)
   static constexpr double width_shock = Config::width_shock;

//! Fraction of the shock width to assign to dmax near shock (persistent)
   static constexpr double dmax_fraction = Config::dmax_fraction;

public:

//! Compute the maximum distance per time step
   template <typename Coordinates>
   static status_t EvaluateDmax(Coordinates&, double*);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

};

};

// Something like this is needed for templated classes
#include "background_smooth_shock.cc"

#endif
