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
class BackgroundSmoothShock {//: public BackgroundShock<HConfig_> {
public:

//! Readable name of the class
   const std::string name = "BackgroundSmoothShock";

public:

   using HConfig = HConfig_;
   using Config = HConfig::BackgroundConfig;
   using BackgroundShock = BackgroundShock<HConfig>;

//   using BackgroundBase = BackgroundBase<HConfig>;
//   using BackgroundBase::_status;
//   using BackgroundBase::container;
//   using BackgroundBase::_ddata;
//   using BackgroundBase::dmax0;
//   using BackgroundBase::r0;
//   using BackgroundBase::u0;
//   using BackgroundBase::B0;
//   // methods
//   using BackgroundBase::EvaluateDmax;
//   using BackgroundBase::GetDmax;

   using BackgroundShock::u1;
   using BackgroundShock::B1;
   using BackgroundShock::n_shock;
   using BackgroundShock::v_shock;

protected:

////! Width of shock transition region (persistent)
//   static double width_shock;
//
////! Fraction of the shock width to assign to dmax near shock (persistent)
//   static double dmax_fraction;

//! Shock transition region function
   static double ShockTransition(double x);

//! Derivative of shock transition region function
   static double ShockTransitionDerivative(double x);

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct);

   //! Compute the maximum distance per time step
   template <typename Coordinates>
   static double EvaluateDmax(Coordinates&);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static void EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   static void EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

//! Default constructor
   BackgroundSmoothShock(void);

//! Copy constructor
   BackgroundSmoothShock(const BackgroundSmoothShock& other);

//! Destructor
   ~BackgroundSmoothShock() = default;

};

};

// Something like this is needed for templated classes
#include "background_smooth_shock.cc"

#endif
