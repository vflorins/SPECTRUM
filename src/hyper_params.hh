/*!
\file hyper_params.hh
\brief (Hyper)parameters for a
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_BACKGROUND_PARAMS_HH
#define SPECTRUM_BACKGROUND_PARAMS_HH

#include "common/definitions.hh"

namespace Spectrum {

/*!
\brief (Hyper)parameters for a SPECTRUM trajectory simulation
\author Lucius Schoenbaum
\date 09/05/2025
Note that template arguments cannot be floating point.
*/
template <
      typename Coordinates_,
      BuildMode build_mode_ = BuildMode::debug,
      DerivativeMethod derivative_method_ = DerivativeMethod::numeric,
      TimeFlow time_flow_ = TimeFlow::forward,
      RKIntegrator rk_integrator_ = RKIntegrator::DormandPrince_54E,
      int num_trajectories = 1,
      int batch_size = 1,
      int num_numeric_grad_evals_ = 1,
      int incr_dmax_ratio_denom = 10000,
      int server_interpolation_order_ = 1,
      int smooth_discontinuity_order_ = 4,
      int server_num_ghost_cells_ = 2,
      SolarWindCurrentSheet solarwind_current_sheet_ = SolarWindCurrentSheet::disabled,
      SolarWindSectoredRegion solarwind_sectored_region_ = SolarWindSectoredRegion::nowhere,
      SolarWindPolarCorrection solarwind_polar_correction_ = SolarWindPolarCorrection::none,
      SolarWindSpeedLatitudeProfile solarwind_speed_latitude_profile_ = SolarWindSpeedLatitudeProfile::constant,
      SolarWindTermShockSpeedExponent solarwind_termshock_speed_exponent_ = SolarWindTermShockSpeedExponent::square,
      VLISMBochumModType VLISM_bochum_mod_type_ = VLISMBochumModType::scaled,
      VLISMBochumModRPos VLISM_bochum_mod_rpos_ = VLISMBochumModRPos::scale_rel_zero
>
struct HyperParams {
// todo line comments

   using Coordinates = Coordinates_;

   static constexpr BuildMode build_mode = build_mode_;

   static constexpr DerivativeMethod derivative_method = derivative_method_;

   static constexpr TimeFlow time_flow = time_flow_;

   static constexpr int num_numeric_grad_evals = num_numeric_grad_evals_;

   static constexpr double incr_dmax_ratio = 1.0/static_cast<double>(incr_dmax_ratio_denom);

// todo impl: SERVER_INTERP_ORDER
   static constexpr int server_interpolation_order = server_interpolation_order_;

//! Parameter controlling smoothness of discontinuity
// todo impl: SMOOTH_DISCONT_ORDER as well as SMOOTH_SHOCK_ORDER (review for merge of these)
   static constexpr int smooth_discontinuity_order = smooth_discontinuity_order_;


   static constexpr int server_num_ghost_cells = server_num_ghost_cells_;


//! Heliospheric current sheet (0: disabled, 1: flat, 2: wavy (Jokipii-Thomas 1981) and static, 3: wavy and time-dependent).
//#define SOLARWIND_CURRENT_SHEET 0
//   [test_parker_spiral.cc] Make sure that "SOLARWIND_CURRENT_SHEET" and SOLARWIND_POLAR_CORRECTION are (#)defined as 0 in src/background_solarwind.hh.
   static constexpr SolarWindCurrentSheet solarwind_current_sheet = solarwind_current_sheet_;



//! Magnetic topology region (0: nowhere, 1: same as HCS)
//#define SOLARWIND_SECTORED_REGION 0
   static constexpr SolarWindSectoredRegion solarwind_sectored_region = solarwind_sectored_region_;


//! Correction to Parker Spiral, mainly for polar regions (0: none, 1: Smith-Bieber 1991, 2: Zurbuchen et al. 1997, 3: Schwadron-McComas 2003)
//#define SOLARWIND_POLAR_CORRECTION 0
   static constexpr SolarWindPolarCorrection solarwind_polar_correction = solarwind_polar_correction_;


//! Latitudinal profile for bulk speed (0: constant, 1: linear step, 2: smooth step)
//#define SOLARWIND_SPEED_LATITUDE_PROFILE 0
   static constexpr SolarWindSpeedLatitudeProfile solarwind_speed_latitude_profile = solarwind_speed_latitude_profile_;



//! Integer exponent of decrease of solar wind speed beyond the termination shock
//#define SOLARWIND_TERMSHOCK_SPEED_EXPONENT 2
   static constexpr SolarWindTermShockSpeedExponent solarwind_termshock_speed_exponent = solarwind_termshock_speed_exponent_;


//! What function to use within 'get_ampfactor' (0 = none, 1 = zero, 2 = constant, 3 = scaled)
//#define MOD_TYPE 3
   static constexpr VLISMBochumModType VLISM_bochum_mod_type = VLISM_bochum_mod_type_;

//! Whether to scale relative to s=0 (0) or s=+inf (1)
//#define MOD_RPOS 0
   static constexpr VLISMBochumModRPos VLISM_bochum_mod_rpos = VLISM_bochum_mod_rpos_;


};

};

// TODO cf. background_waves.cc, we must locate the point where Backgrounds are being "RNGConnect"ed, and modify that init logic (so that it is more a-la-carte)

#endif
