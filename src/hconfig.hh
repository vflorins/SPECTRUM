/*!
\file hconfig.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM trajectory simulation
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_HCONFIG_HH
#define SPECTRUM_HCONFIG_HH

#include "common/definitions.hh"

namespace Spectrum {

/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM trajectory simulation
\author Lucius Schoenbaum
\date 09/05/2025
Note that template arguments cannot be floating point (cf. dmax ratio).
*/
template <
      typename Coordinates_, // live, including vel/mom typ.
      typename StoreCoordinates_, // only including mom typ. todo
      typename TrajectoryFields_,
      typename DiffusionFields_,
      BuildMode build_mode_ = BuildMode::debug,
      Specie specie_ = Specie::proton_core,
      DerivativeMethod derivative_method_ = DerivativeMethod::numeric,
      TimeFlow time_flow_ = TimeFlow::forward,
      RKIntegrator rk_integrator_ = RKIntegrator::DormandPrince_54E,
      int num_trajectories = 1,
      int batch_size = 1,
      // ----- Background ----- //
      int num_numeric_grad_evals_ = 1,
      std::ratio incr_dmax_ratio_ = {1,10000},
      int server_interpolation_order_ = 1,
      int smooth_discontinuity_order_ = 4,
      int server_num_ghost_cells_ = 2,
      SolarWindCurrentSheet solarwind_current_sheet_ = SolarWindCurrentSheet::disabled,
      SolarWindSectoredRegion solarwind_sectored_region_ = SolarWindSectoredRegion::nowhere,
      SolarWindPolarCorrection solarwind_polar_correction_ = SolarWindPolarCorrection::none,
      SolarWindSpeedLatitudeProfile solarwind_speed_latitude_profile_ = SolarWindSpeedLatitudeProfile::constant,
      SolarWindTermShockSpeedExponent solarwind_termshock_speed_exponent_ = SolarWindTermShockSpeedExponent::square,
      VLISMBochumModType VLISM_bochum_mod_type_ = VLISMBochumModType::scaled,
      VLISMBochumModRPos VLISM_bochum_mod_rpos_ = VLISMBochumModRPos::scale_rel_zero,
      // ----- Trajectory ----- //
      bool record_mag_extrema_ = false,
      bool record_trajectory_ = false,
      size_t record_trajectory_segment_presize_ = 0,
      int trajectory_adv_safety_level_ = 2
>
struct HConfig {
// todo line comments

   using Coordinates = Coordinates_;

   using TrajectoryFields = TrajectoryFields_;

   using DiffusionFields = DiffusionFields_;

   static constexpr BuildMode build_mode = build_mode_;

   static constexpr Specie specie = specie_;

   static constexpr DerivativeMethod derivative_method = derivative_method_;

   static constexpr TimeFlow time_flow = time_flow_;

   static constexpr int num_numeric_grad_evals = num_numeric_grad_evals_;

//! What fraction of "_dmax" to use to calculate the field increment
   static constexpr double incr_dmax_ratio = static_cast<double>(incr_dmax_ratio_);

// todo impl: SERVER_INTERP_ORDER
   static constexpr int server_interpolation_order = server_interpolation_order_;

//! Parameter controlling smoothness of discontinuity/shock
// 0: not continuous
// 1: differentiable
// 2: twice differentiable
// 3: thrice differentiable
// 4, 5, ... (any other value): smooth
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

// ----- Trajectory ----- //

   static constexpr bool record_mag_extrema = record_mag_extrema_;

   static constexpr bool record_trajectory = record_trajectory_;

   static constexpr size_t record_trajectory_segment_presize = record_trajectory_segment_presize_;

//! Trajectory advance safety level: 0 means no checks, 1 means check dt only, 2 means check dt, number of segments, and time adaptations per step.
   static constexpr int trajectory_adv_safety_level = trajectory_adv_safety_level_;

//#if TRAJ_ADV_SAFETY_LEVEL == 2 todo ??
//! Largest length for single trajectory
   static constexpr int max_trajectory_steps = 10000000;

//! Largest number of time step adaptations for a single time step
   static constexpr int max_time_adaptations = 100;
//#endif

//! Upper limit on the number of steps in debug mode
//! use -1 for unlimited
   static constexpr unsigned int n_max_calls = -1;


   // todo the charge to mass ratio
   // The factor multiplying "SpeciesCharges[]" is applied in order to marry particle and fluid scales. See "LarmorRadius()" and "CyclotronFrequency()" functions in physics.hh for reference.
   static constexpr double q = charge_mass_particle * SpeciesCharges[static_cast<size_t>(specie)];




};

};

#endif
