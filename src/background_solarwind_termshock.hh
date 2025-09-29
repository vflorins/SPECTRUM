/*!
\file background_solarwind_termshock.hh
\brief Declares a plasma background class for the constant speed supersonic wind of a rotating star
\author Juan G Alonso Guzman
\author Vladimir Florinski

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SOLARWIND_TERMSHOCK_HH
#define SPECTRUM_BACKGROUND_SOLARWIND_TERMSHOCK_HH

#include "background_solarwind.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWindTermShock class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background calculator for a radially expanding solar wind with a spherical termination shock
\author Juan G Alonso Guzman

Parameters: (BackgroundSolarWind), double r_TS, double w_TS, double s_TS
*/
template <typename HConfig_>
class BackgroundSolarWindTermShock : public BackgroundSolarWind<HConfig_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundSolarWindTermShock";

public:

   using HConfig = HConfig_;
   using BackgroundSolarWind = BackgroundSolarWind<HConfig>;
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
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;
   using BackgroundBase::NumericalDerivatives;

   using BackgroundSolarWind::dmax_fraction;
   using BackgroundSolarWind::ur0;

protected:

//! Radius of termination shock (persistent)
   double r_TS;

//! Width of termination shock (persistent)
   double w_TS;

//! Strength of termination shock (persistent)
   double s_TS;

//! Inverse of s_TS (persistent)
   double s_TS_inv;

//! Maximum displacement in the shock region (persistent)
   double dmax_TS;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct) override;

//! Compute the maximum distance per time step
   void EvaluateDmax(BackgroundCoordinates&) override;

//! Modify radial flow (if necessary)
   void ModifyUr(double r, double &ur_mod) override;

//! Radial derivative of radial flow
   double dUrdr(double r, double v_norm);

//! Get time lag for time dependent current sheet (if necessary)
   double TimeLag(double r) override;

//! Compute the internal u, B, and E fields
   template <typename Fields>
   void EvaluateBackground(BackgroundCoordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Fields>
   void EvaluateBackgroundDerivatives(BackgroundCoordinates&, Fields&);

public:

//! Default constructor
   BackgroundSolarWindTermShock(void);

//! Copy constructor
   BackgroundSolarWindTermShock(const BackgroundSolarWindTermShock& other);

//! Destructor
   ~BackgroundSolarWindTermShock() override = default;

//! Clone function
   CloneFunctionBackground(BackgroundSolarWindTermShock);

};

};

// Something like this is needed for templated classes
#include "background_solarwind_termshock.cc"

#endif
