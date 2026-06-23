/*!
\file background_solarwind.hh
\brief Declares a plasma background class for the constant speed supersonic wind of a rotating star
\author Vladimir Florinski
\author Juan G Alonso Guzman

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_SOLARWIND_HH
#define SPECTRUM_BACKGROUND_SOLARWIND_HH

#include "common/vectors.hh"

// todo replace later - still early stages testing
constexpr double constexpr_sin_est(double x) {
   return x - x/3.0*x/2.0*x + x/5.0*x/4.0*x/3.0*x/2.0*x;
}
constexpr double constexpr_cos_est(double x) {
   return 1.0 - x/2.0*x + x/4.0*x/3.0*x/2.0*x;
}


template <bool with_terminations_shock>
struct TerminationShock;


template<>
struct TerminationShock<true> {
//! Radius of termination shock (persistent)
   double r;
//! Width of termination shock (persistent)
   double w;
//! Strength of termination shock (persistent)
   double s;
//! Inverse of s_TS (persistent)
   double s_inv;
//! Maximum displacement in the shock region (persistent)
   double dmax;
};

template<>
struct TerminationShock<false> {};




namespace Spectrum {

//----------------------------------------------------------------------------------------------------------------------------------------------------
// BackgroundSolarWind class declaration
//----------------------------------------------------------------------------------------------------------------------------------------------------

/*!
\brief Plasma background calculator for a radially expanding solar wind with (optional) spherical termination shock
\author Vladimir Florinski
\author Juan G Alonso Guzman
\author Lucius Schoenbaum

Parameters: (BackgroundBase), GeoVector Omega, double r_ref, double dmax_fraction
*/
template <typename HConfig_>
class BackgroundSolarWind {
public:

//! Readable name of the class
   static constexpr std::string_view name = "BackgroundSolarWind";

public:

   using HConfig = HConfig_;
   using BackgroundConfig = HConfig::BackgroundConfig;

// secular config:
   static constexpr bool requires_setup = true;
   static constexpr bool stochastic = false;

   static constexpr auto derivative_method = BackgroundConfig::derivative_method;
   static constexpr auto solarwind_speed_latitude_profile = BackgroundConfig::solarwind_speed_latitude_profile;
   static constexpr auto solarwind_current_sheet = BackgroundConfig::solarwind_current_sheet;
   static constexpr auto solarwind_sectored_region = BackgroundConfig::solarwind_sectored_region;
   static constexpr auto solarwind_polar_correction = BackgroundConfig::solarwind_polar_correction;
   static constexpr auto with_termination_shock = BackgroundConfig::with_termination_shock;
   static constexpr auto termshock_speed_exponent = BackgroundConfig::termshock_speed_exponent;

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

//! [solarwind_polar_correction == 2] Ratio of polar differential rotation to angular frequency of rotation
   static constexpr double dwt_sw = delta_omega_sw * constexpr_sin_est(polar_offset_sw);

//! [solarwind_polar_correction == 2] Ratio of azimuthal differential rotation to angular frequency of rotation
   static constexpr double dwp_sw = delta_omega_sw * constexpr_cos_est(polar_offset_sw);

//! [solarwind_speed_latitude_profile > 0] Ratio of fast to slow wind speed
   static constexpr double fast_slow_ratio_sw = 2.0;

//! [solarwind_speed_latitude_profile > 0] Angle to transition between fast and slow
   static constexpr double fast_slow_lat_sw = M_PI_2 - 30.0 * M_PI / 180.0 - tilt_ang_sw;

//! [solarwind_speed_latitude_profile > 0] Transition speed coefficient
   static constexpr double fast_slow_dlat_sw = 20.0;

// Compute auxiliary quantities for fast-slow wind calculation
   static constexpr double construct_fsl_pls() {
      if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::constant) {
         /* unused*/
         return 0;
      }
      if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::linear_step) {
//! Latitude separating transition region from slow wind
         return fast_slow_lat_sw + 2.0 / fast_slow_dlat_sw;
      }
      else if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::smooth_step) {
//! Half of fast-slow ratio plus 1
         return 0.5 * (fast_slow_ratio_sw + 1.0);
      }
      return 0;
   }

// Compute auxiliary quantities for fast-slow wind calculation
   static constexpr double construct_fsl_mns() {
      if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::constant) {
         /* unused*/
         return 0;
      }
      if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::linear_step) {
//! Latitude separating transition region from fast wind
         return fast_slow_lat_sw - 2.0 / fast_slow_dlat_sw;
      }
      else if constexpr (solarwind_speed_latitude_profile == SpeedLatitudeProfile::smooth_step) {
//! Half of fast-slow ratio minus 1
         return 0.5 * (fast_slow_ratio_sw - 1.0);
      }
      return 0;
   }

/*!
\brief Function to compress peaks and stretch troughs
\author Juan G Alonso Guzman
\date 07/16/2024
\param[in] t periodic time to stretch between 0 and M_2PI
\return stretched time
*/
   static constexpr double CubicStretch(double t)
   {
      double t_pi = t / M_PI;
      return t * ((1.0 - stilt_ang_sw) * t_pi * (t_pi - 3.0) + 3.0 - 2.0 * stilt_ang_sw);
   };

   static constexpr double fsl_pls = construct_fsl_pls();

   static constexpr double fsl_mns = construct_fsl_mns();

protected:

   double dmax0;

   double t0;

   GeoVector r0;

   GeoVector u0;

   GeoVector B0;

   //! Velocity magnitude for slow wind (persistent)
   double ur0;

   //! Angular velocity vector of a rotating star (persistent)
   GeoVector Omega;

//! Reference radius (persistent)
   double r_ref;

//! Maximum fraction of the radial distance per step (persistent)
   double dmax_fraction;

//! Local coordinate system tied to the rotation axis (persistent)
   GeoVector eprime[3];

//! Angular frequency magnitude (persistent)
   double w0;

   TerminationShock<with_termination_shock> TS;

public:

//! Set up the field evaluator based on "params"
   void SetupBackground(DataContainer& container);

//! Modify radial flow (if necessary)
   virtual void ModifyUr(const double r, double &ur_mod);

//! Radial derivative of radial flow
   double dUrdr(double r, double v_norm);

   //! Get time lag for time dependent current sheet (if necessary)
   virtual double TimeLag(const double r);

//! Compute the maximum distance per time step
   template <typename Coordinates>
   status_t EvaluateDmax(Coordinates&, double*);

//! Compute the internal u, B, and E fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackground(Coordinates&, Fields&);

//! Compute the internal derivatives of the fields
   template <typename Coordinates, typename Fields, typename RequestedFields>
   status_t EvaluateBackgroundDerivatives(Coordinates&, Fields&);

public:

////! Default constructor
//   BackgroundSolarWind(void);
//
////! Constructor with arguments (to speed up construction of derived classes)
//   BackgroundSolarWind(const std::string_view& name_in, status_t status_in);
//
////! Copy constructor
//   BackgroundSolarWind(const BackgroundSolarWind& other);
//
////! Destructor
//   ~BackgroundSolarWind() = default;

};

};

// Something like this is needed for templated classes
#include "background_solarwind.cc"

#endif
