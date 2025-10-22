// config header file for main_test_dipole_periods.cc
// Created by Lucius Schoenbaum September 28, 2025

#pragma once

#include "common/fields.hh"
#include "src/hconfig.hh"

using namespace Spectrum;

using BackgroundCoordinates = Fields<
      FConfig<
            CoordinateSystem::cartesian,
            CoordinateSystem::p_mu_phi
      >,
      Time_t,
      Pos_t,
      Mom_t
>;

using BackgroundFields = Fields<
      FConfig<>,
      // todo
      Mag_t
>;

using TrajectoryCoordinates = Fields<
      FConfig<
            CoordinateSystem::cartesian,
            // todo review
            CoordinateSystem::p_mu_phi
      >,
      Time_t,
      Pos_t,
      Vel_t,
      Mom_t
>;

using TrajectoryFields = Fields<
      FConfig<>,
      Mag_t
>;

using RecordCoordinates = Fields<
      FConfig<>,
      Time_t,
      Vel_t
>;

using DiffusionCoordinates = Fields<
      FConfig<
            CoordinateSystem::cartesian,
            // todo review
            CoordinateSystem::p_mu_phi
      >,
      Time_t,
      Pos_t,
      Vel_t,
      Mom_t
>;

using DiffusionFields = Fields<
      FConfig<>,
      Mag_t
>;



using HConfig = HConfig<
BackgroundCoordinates,
BackgroundFields,
TrajectoryCoordinates,
TrajectoryFields,
RecordCoordinates,
DiffusionCoordinates,
DiffusionFields,
BuildMode::debug,
SpecieId::proton_core,
DerivativeMethod::numeric,
TimeFlow::forward,
RKIntegrator::DormandPrince_54E,
//int num_trajectories,
1,
//int batch_size,
1,
//// ----- Background ----- //
//int num_numeric_grad_evals_,
1,
////! What fraction of "_dmax" to use to calculate the field increment
//std::ratio incr_dmax_ratio_,
std::ratio<1,1>(),
//// todo impl: SERVER_INTERP_ORDER
//int server_interpolation_order_,
1,
////! Parameter controlling smoothness of discontinuity/shock
//// 0: not continuous
//// 1: differentiable
//// 2: twice differentiable
//// 3: thrice differentiable
//// 4, 5, ... (any other value): smooth
//// todo impl: SMOOTH_DISCONT_ORDER as well as SMOOTH_SHOCK_ORDER (review for merge of these)
//int smooth_discontinuity_order_,
1,
//////! number of ghost cells (server parameter)
//int server_num_ghost_cells_,
2,
////! Heliospheric current sheet (0: disabled, 1: flat, 2: wavy (Jokipii-Thomas 1981) and static, 3: wavy and time-dependent).
////#define SOLARWIND_CURRENT_SHEET 0
////   [test_parker_spiral.cc] Make sure that "SOLARWIND_CURRENT_SHEET" and SOLARWIND_POLAR_CORRECTION are (#)defined as 0 in src/background_solarwind.hh.
//CurrentSheet solarwind_current_sheet_,
CurrentSheet::disabled,
////! Magnetic topology region (0: nowhere, 1: same as HCS)
////#define SOLARWIND_SECTORED_REGION 0
//SectoredRegion solarwind_sectored_region_,
SectoredRegion::nowhere,
////! Correction to Parker Spiral, mainly for polar regions (0: none, 1: Smith-Bieber 1991, 2: Zurbuchen et al. 1997, 3: Schwadron-McComas 2003)
////#define SOLARWIND_POLAR_CORRECTION 0
//PolarCorrection solarwind_polar_correction_,
PolarCorrection::none,
////! Latitudinal profile for bulk speed (0: constant, 1: linear step, 2: smooth step)
////#define SOLARWIND_SPEED_LATITUDE_PROFILE 0
//SpeedLatitudeProfile solarwind_speed_latitude_profile_,
SpeedLatitudeProfile::constant,
////! Integer exponent of decrease of solar wind speed beyond the termination shock
////#define SOLARWIND_TERMSHOCK_SPEED_EXPONENT 2
//TermShockSpeedExponent solarwind_termshock_speed_exponent_,
TermShockSpeedExponent::one,
////! What function to use within 'get_ampfactor' (0 = none, 1 = zero, 2 = constant, 3 = scaled)
////#define MOD_TYPE 3
//ModType mod_type_,
ModType::constant,
////! Whether to scale relative to s=0 (0) or s=+inf (1)
////#define MOD_RPOS 0
//ModRPos mod_rpos_,
ModRPos::scale_rel_inf,
//// ----- Trajectory ----- //
//////! Whether to record magnetic field extrema
//bool record_mag_extrema_,
true,
//////! Whether to record segments
//bool record_trajectory_,
true,
//////! Default initial size of trajectory segment record
//size_t record_trajectory_segment_presize_,
10000,
////! Trajectory advance safety level: 0 means no checks, 1 means check dt only, 2 means check dt, number of segments, and time adaptations per step.
//TrajectoryOptions::SafetyLevel trajectory_adv_safety_level_,
TrajectoryOptions::SafetyLevel::high,
////! Largest length for single trajectory
//int max_trajectory_steps_,
100000,
////! Largest number of time step adaptations for a single time step
//int max_time_adaptations_,
10,
////! Upper limit on the number of steps in debug mode
////! use -1 for unlimited
//int n_max_calls_,
-1,
//////! CFL condition for advection
//std::ratio cfl_advection_,
std::ratio<1,10>(),
//////! CFL condition for diffusion
//std::ratio cfl_diffusion_,
std::ratio<1,10>(),
//////! CFL condition for acceleration
//std::ratio cfl_acceleration_,
std::ratio<1,10>(),
//////! CFL condition for pitch angle scattering
//std::ratio cfl_pitchangle_,
std::ratio<1,10>(),
////! Safety factor for drift-based time step (to modify "drift_vel" with a small fraction of the particle's velocity)
//std::ratio drift_safety_,
std::ratio<1,10>(),
////! How many time steps to allow before recording a mirror event
//std::ratio mirror_threshold_,
10,
////// Switch controlling how to calculate mu. "0" means computing it at the end of the step from magnetic moment conservation. "1" means advancing it in time according to the scheme (does not guarantee conservation of MM, but can be used with non-adiabatic terms).
//TrajectoryOptions::PPerpMethod pperp_method_,
TrajectoryOptions::PPerpMethod::mag_moment_conservation,
////! Flag to use gradient and curvature drifts in drift velocity calculation
//TrajectoryOptions::UseBDrifts use_B_drifts_,
TrajectoryOptions::UseBDrifts::none,
////! Which stochastic method to use for PA scattering, 0 = Euler, 1 = Milstein, 2 = RK2
//////! Which stochastic method to use for perpendicular diffusion, 0 = Euler, 1 = Milstein, 2 = RK2
//TrajectoryOptions::StochasticMethod stochastic_method_,
TrajectoryOptions::StochasticMethod::Euler,
//TrajectoryOptions::StochasticMethod stochastic_method_mu_,
TrajectoryOptions::StochasticMethod::Euler,
//TrajectoryOptions::StochasticMethod stochastic_method_perp_,
TrajectoryOptions::StochasticMethod::Euler,
//////! Whether to split the diffusive advance into two (one before and one after the advection).
////// #define SPLIT_SCATT
//bool split_scatt_,
true,
////#ifdef SPLIT_SCATT
//////! Fraction of stochastic step to take before deterministic step
////const double alpha = 0.5;
////#endif
//std::ratio split_scatt_fraction_,
std::ratio<1,10>(),
////#if CONST_DMUMAX == 1
//////! Desired accuracy in pitch angle cosine
////const double dmumax = 0.02;
////#else
//////! Desired accuracy in pitch angle (deg x [deg to rad conversion factor])
////const double dthetamax = 2.0 * M_PI / 180.0;
////#endif
//TrajectoryOptions::ConstDmumax const_dmumax_,
TrajectoryOptions::ConstDmumax::constant_dmu_max,
//////! Number of time steps per one orbit
//int steps_per_orbit_,
1,
//////! Which method of computation to use for divK: 0 = using direct central FD, 1 = using _spdata.grad quantities
//TrajectoryOptions::DivkMethod divk_method_,
TrajectoryOptions::DivkMethod::direct,
//////! Maximum allowed fraction of momentum change per step
//std::ratio dlnp_max_,
std::ratio<1,2>(),
//// ----- Diffusion ----- //
////! Flag to use QLT pitch angle scattering with WLNT perpendicular diffusion
//bool use_qlt_scatt_,
1,
////! diffusion coefficient (if constant)
////! Scattering frequency (persistent)
//std::ratio D0_,
std::ratio<1,1>(),
//std::ratio Dperp_,
      std::ratio<1,1>(),
//std::ratio Dpara_,
      std::ratio<1,1>(),
//std::ratio A2A_,
      std::ratio<1,1>(),
//std::ratio l_max_,
      std::ratio<1,1>(),
//std::ratio k_min_,
      std::ratio<1,1>(),
//std::ratio ps_index_,
      std::ratio<1,1>(),
//std::ratio ps_minus_,
      std::ratio<1,1>(),
//std::ratio A2T_,
      std::ratio<1,1>(),
//std::ratio A2L_,
      std::ratio<1,1>(),
//std::ratio ps_plus_,
      std::ratio<1,1>(),
//std::ratio l_max_HP_,
      std::ratio<1,1>(),
//std::ratio dl_max_,
      std::ratio<1,1>(),
//std::ratio nose_z_nose_,
      std::ratio<1,1>(),
//std::ratio nose_z_sheath_,
      std::ratio<1,1>(),
//std::ratio nose_dz_,
      std::ratio<1,1>(),
//std::ratio kappa0_,
      std::ratio<1,1>(),
//std::ratio kappa_ratio_,
      std::ratio<1,1>(),
//std::ratio U0_,
      std::ratio<1,1>(),
//std::ratio p0_,
      std::ratio<1,1>(),
//std::ratio T0_,
      std::ratio<1,1>(),
//std::ratio r0_,
      std::ratio<1,1>(),
//std::ratio lam0_,
      std::ratio<1,1>(),
//std::ratio R0_,
      std::ratio<1,1>(),
//std::ratio B0_,
      std::ratio<1,1>(),
//std::ratio pow_law_U_,
      std::ratio<1,1>(),
//std::ratio pow_law_p_,
      std::ratio<1,1>(),
//std::ratio pow_law_T_,
      std::ratio<1,1>(),
//std::ratio pow_law_r_,
      std::ratio<1,1>(),
//std::ratio pow_law_R_,
      std::ratio<1,1>(),
//std::ratio pow_law_B_,
      std::ratio<1,1>(),
//int LISM_idx_,
   1,
//std::ratio LISM_ind_,
      std::ratio<1,1>(),
//std::ratio kappa_ratio_inner_,
      std::ratio<1,1>(),
//std::ratio kappa_ratio_outer_,
      std::ratio<1,1>(),
//std::ratio lam_inner_,
      std::ratio<1,1>(),
//std::ratio lam_outer_,
      std::ratio<1,1>(),
//std::ratio kappa_inner_,
      std::ratio<1,1>(),
//std::ratio kappa_outer_,
      std::ratio<1,1>(),
//std::ratio lam_para_,
      std::ratio<1,1>(),
//std::ratio lam_perp_,
      std::ratio<1,1>(),
//std::ratio kappa_ratio_red_,
      std::ratio<1,1>(),
//std::ratio radial_limit_perp_red_,
      std::ratio<1,1>(),
//std::ratio solar_cycle_idx_,
      std::ratio<1,1>(),
//std::ratio solar_cycle_effect_,
      std::ratio<1,1>(),
//int Bmix_idx_,
1,
//std::ratio Bmix_ind_
      std::ratio<1,1>()
>;




