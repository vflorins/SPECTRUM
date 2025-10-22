/*!
\file background_smooth_discontinuity.hh
\brief Declares a smooth MHD discontinuity field background
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a SHOCK coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SMOOTH_DISCONTINUITY_HH
#define SPECTRUM_BACKGROUND_SMOOTH_DISCONTINUITY_HH

#include "background_discontinuity.hh"

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
class BackgroundSmoothDiscontinuity : public BackgroundDiscontinuity<HConfig_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundSmoothDiscontinuity";

public:

   using HConfig = HConfig_;
   using BackgroundDiscontinuity = BackgroundDiscontinuity<HConfig>;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundSmoothDiscontinuity<HConfig>>, typename HConfig::BackgroundConfig>;

   using BackgroundCoordinates = BackgroundConfig::Coordinates;
   using BackgroundBase = BackgroundBase<HConfig>;
   using BackgroundBase::_status;
   using BackgroundBase::container;
   using BackgroundBase::_ddata;
   using BackgroundBase::dmax0;
   using BackgroundBase::r0;
   using BackgroundBase::u0;
   using BackgroundBase::B0;
   // methods
   using BackgroundBase::EvaluateAbsMag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;

   using BackgroundDiscontinuity::n_discont;
   using BackgroundDiscontinuity::v_discont;
   using BackgroundDiscontinuity::u1;
   using BackgroundDiscontinuity::B1;

   using BackgroundConfig::derivative_method;
   using BackgroundConfig::smooth_discontinuity_order;

//! Scaling factor to better match discontinuity width when using smooth discontinuity (tanh)
   const double tanh_width_factor = 4.0;

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
   void SetupBackground(bool construct);

//! Compute the maximum distance per time step
   template <typename Coordinates>
   void EvaluateDmax(Coordinates&);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

//! Default constructor
   BackgroundSmoothDiscontinuity(void);

//! Copy constructor
   BackgroundSmoothDiscontinuity(const BackgroundSmoothDiscontinuity& other);

//! Destructor
   ~BackgroundSmoothDiscontinuity() = default;

//! Clone function
   CloneFunctionBackground(BackgroundSmoothDiscontinuity);

};

};

// Something like this is needed for templated classes
#include "background_smooth_discontinuity.cc"

#endif
