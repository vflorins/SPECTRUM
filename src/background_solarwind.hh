/*!
\file background_solarwind.hh
\brief Declares a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SOLARWIND_HH
#define SPECTRUM_BACKGROUND_SOLARWIND_HH

#include "background_base.hh"

namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background calculator for a radially expanding solar wind
\author Vladimir Florinski
\author Juan G Alonso Guzman

Parameters: (BackgroundBase), GeoVector Omega, double r_ref, double dmax_fraction
*/
template <typename HConfig_>
class BackgroundSolarWind : public BackgroundBase<HConfig_> {
private:

//! Readable name of the class
   static constexpr std::string_view bg_name = "BackgroundSolarWind";

public:

   using HConfig = HConfig_;
   using BackgroundConfig = Cond<std::same_as<typename HConfig::BackgroundConfig, Default>, BackgroundDefault<BackgroundSolarWind<HConfig>>, typename HConfig::BackgroundConfig>;
   using BackgroundCoordinates = BackgroundConfig::Coordinates;
   using BackgroundBase = BackgroundBase<HConfig>;
   using BackgroundBase::_status;
   using BackgroundBase::container;
   using BackgroundBase::_ddata;
   using BackgroundBase::dmax0;
   using BackgroundBase::r0;
   using BackgroundBase::u0;
   using BackgroundBase::B0;
   //
   // todo review (solar wind is the only background using this)
   using BackgroundBase::t0;
   // methods
   using BackgroundBase::EvaluateAbsMag;
   using BackgroundBase::EvaluateDmax;
   using BackgroundBase::GetDmax;
   using BackgroundBase::StopServerFront;
   using BackgroundBase::SetupBackground;

   using BackgroundConfig::derivative_method;
   using BackgroundConfig::solarwind_speed_latitude_profile;
   using BackgroundConfig::solarwind_current_sheet;
   using BackgroundConfig::solarwind_sectored_region;
   using BackgroundConfig::solarwind_polar_correction;

protected:

//! Heliopause radius
#ifdef USE_GSL
   static constexpr double hp_rad_sw = 117.0 * GSL_CONST_CGSM_ASTRONOMICAL_UNIT / unit_length_fluid;
#else
   static constexpr double hp_rad_sw = -1.0;
#endif

//! Magnetic axis tilt angle relative to the solar rotation axis
   static constexpr double tilt_ang_sw = 40.0 * M_PI / 180.0;

//! [solarwind_current_sheet == 3] Amplitude of variation to magnetic axis tilt angle
   static constexpr double dtilt_ang_sw = 35.0 * M_PI / 180.0;

//! [solarwind_current_sheet == 3] Solar cycle frequency
   static constexpr double W0_sw = M_2PI / (60.0 * 60.0 * 24.0 * 365.0 * 22.0) / unit_frequency_fluid;

//! [solarwind_current_sheet == 3] Factor to thin peaks and widen troughs (0.0: largest modification, 1.0: no modification)
   static constexpr double stilt_ang_sw = 1.0;

//! [solarwind_polar_correction > 0] Differential rotation factor
   static constexpr double delta_omega_sw = 0.05;

//! [solarwind_polar_correction == 2] Polar correction angle
   static constexpr double polar_offset_sw = 30.0 * M_PI / 180.0;

//! [solarwind_speed_latitude_profile > 0] Ratio of fast to slow wind speed
   static constexpr double fast_slow_ratio_sw = 2.0;

//! [solarwind_speed_latitude_profile > 0] Angle to transition between fast and slow
   static constexpr double fast_slow_lat_sw = M_PI_2 - 30.0 * M_PI / 180.0 - tilt_ang_sw;

//! [solarwind_speed_latitude_profile > 0] Transition speed coefficient
   static constexpr double fast_slow_dlat_sw = 20.0;


/*!
\brief Function to compress peaks and stretch troughs
\author Juan G Alonso Guzman
\date 07/16/2024
\param[in] t periodic time to stretch between 0 and M_2PI
\return stretched time
*/
   constexpr double CubicStretch(double t)
   {
      double t_pi = t / M_PI;
      return t * ((1.0 - stilt_ang_sw) * t_pi * (t_pi - 3.0) + 3.0 - 2.0 * stilt_ang_sw);
   };

protected:

   //! [solarwind_polar_correction == 2] Ratio of polar differential rotation to angular frequency of rotation
   const double dwt_sw = delta_omega_sw * sin(polar_offset_sw);

//! [solarwind_polar_correction == 2] Ratio of azimuthal differential rotation to angular frequency of rotation
   const double dwp_sw = delta_omega_sw * cos(polar_offset_sw);

   //! Angular velocity vector of a rotating star (persistent)
   GeoVector Omega;

//! Reference radius (persistent)
   double r_ref;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

//! Local coordinate system tied to the rotation axis (persistent)
   GeoVector eprime[3];

//! Velocity magnitude for slow wind (persistent)
   double ur0;

//! Radial magnetic field at "r_ref" (persistent)
   double Br0;

//! Angular frequency magnitude (persistent)
   double w0;

//! Position relative to origin (transient)
   GeoVector posprime;

//! [SOLARWIND_SPEED_LATITUDE_PROFILE == 1] Latitude separating transition region from slow wind (persistent)
//! [SOLARWIND_SPEED_LATITUDE_PROFILE == 2] Half of fast-slow ratio plus 1 (persistent)
   double fsl_pls;

//! [SOLARWIND_SPEED_LATITUDE_PROFILE == 1] Latitude separating transition region from fast wind (persistent)
//! SOLARWIND_SPEED_LATITUDE_PROFILE == 2] Half of fast-slow ratio minus 1 (persistent)
   double fsl_mns;

//! Set up the field evaluator based on "params"
   void SetupBackground(bool construct);

//! Modify radial flow (if necessary)
   virtual void ModifyUr(const double r, double &ur_mod);

//! Get time lag for time dependent current sheet (if necessary)
   virtual double TimeLag(const double r);

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
   BackgroundSolarWind(void);

//! Constructor with arguments (to speed up construction of derived classes)
   BackgroundSolarWind(const std::string& name_in, uint16_t status_in);

//! Copy constructor
   BackgroundSolarWind(const BackgroundSolarWind& other);

//! Destructor
   ~BackgroundSolarWind() = default;

//! Clone function
   CloneFunctionBackground(BackgroundSolarWind);

};

};

// Something like this is needed for templated classes
#include "background_solarwind.cc"

#endif
