/*!
\file background_smooth_discontinuity.hh
\brief Declares a smooth MHD discontinuity field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SMOOTH_DISCONTINUITY_HH
#define SPECTRUM_BACKGROUND_SMOOTH_DISCONTINUITY_HH

#include "background_discontinuity.hh"
#include "background_smooth_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSmoothDiscontinuity class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Planar MHD discontinuity with a smooth transition region
\author Juan G Alonso Guzman

Parameters: (BackgroundShock), double width_discont, double dmax_fraction
*/
template <typename HConfig_>
class BackgroundSmoothDiscontinuity: public BackgroundDiscontinuity<HConfig_>, BackgroundSmoothBase<HConfig_> {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundSmoothDiscontinuity";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = false;
   static constexpr bool stochastic = false;

   using BackgroundSmoothBase = BackgroundSmoothBase<HConfig>;
   using BackgroundSmoothBase::Transition;
   using BackgroundSmoothBase::TransitionDerivative;

   using BackgroundDiscontinuity = BackgroundDiscontinuity<HConfig>;
   using BackgroundDiscontinuity::n_discont;
   using BackgroundDiscontinuity::v_discont;
   using BackgroundDiscontinuity::dmax0;
   using BackgroundDiscontinuity::r0;
   using BackgroundDiscontinuity::u0;
   using BackgroundDiscontinuity::u1;
   using BackgroundDiscontinuity::B0;
   using BackgroundDiscontinuity::B1;

public:

//! Width of discontinuity transition region (persistent)
   static constexpr double width_discont = Config::width_discont;

//! Fraction of the discontinuity width to assign to dmax near discontinuity (persistent)
   static constexpr double dmax_fraction = Config::dmax_fraction;

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

//! Compute the maximum distance per time step
   template <typename Coordinates>
   static status_t EvaluateDmax(Coordinates&, double*);

};

};

// Something like this is needed for templated classes
#include "background_smooth_discontinuity.cc"

#endif
