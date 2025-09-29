/*!
\file hconfig.hh
\brief (Hyper)parameters and config(uration) options for a SPECTRUM trajectory simulation
\author Lucius Schoenbaum

This file is part of the SPECTRUM suite of scientific numerical simulation codes. SPECTRUM stands for Space Plasma and Energetic Charged particle TRansport on Unstructured Meshes. The code simulates plasma or neutral particle flows using MHD equations on a grid, transport of cosmic rays using stochastic or grid based methods. The "unstructured" part refers to the use of a geodesic mesh providing a uniform coverage of the surface of a sphere.
*/

#ifndef SPECTRUM_HCONFIG_HH
#define SPECTRUM_HCONFIG_HH

#include "common/definitions.hh"
#include <ratio>

namespace Spectrum {

/*!
\brief (Hyper)parameters and config(uration) options for a SPECTRUM trajectory simulation
\author Lucius Schoenbaum
\date 09/05/2025
Default values are set in the config.py script.
Because template arguments cannot be floating point type, std::ratio is used for these.
*/
template <
typename BackgroundCoordinates_, // while live, including both vel and mom.
typename BackgroundFields_,
typename TrajectoryCoordinates_,
typename TrajectoryFields_,
typename RecordCoordinates_, // only including mom. todo
typename DiffusionCoordinates_,
typename DiffusionFields_,
BuildMode build_mode_,
SpecieId specieid,
DerivativeMethod derivative_method_,
TimeFlow time_flow_,
RKIntegrator rk_integrator_,
int num_trajectories,
int batch_size,
// ----- Background ----- //
int num_numeric_grad_evals_,
//! What fraction of "_dmax" to use to calculate the field increment
std::ratio incr_dmax_ratio_,
// todo impl: SERVER_INTERP_ORDER
int server_interpolation_order_,
//! Parameter controlling smoothness of discontinuity/shock
// 0: not continuous
// 1: differentiable
// 2: twice differentiable
// 3: thrice differentiable
// 4, 5, ... (any other value): smooth
// todo impl: SMOOTH_DISCONT_ORDER as well as SMOOTH_SHOCK_ORDER (review for merge of these)
int smooth_discontinuity_order_,
////! number of ghost cells (server parameter)
int server_num_ghost_cells_,
//! Heliospheric current sheet (0: disabled, 1: flat, 2: wavy (Jokipii-Thomas 1981) and static, 3: wavy and time-dependent).
//#define SOLARWIND_CURRENT_SHEET 0
//   [test_parker_spiral.cc] Make sure that "SOLARWIND_CURRENT_SHEET" and SOLARWIND_POLAR_CORRECTION are (#)defined as 0 in src/background_solarwind.hh.
SolarWindOptions::CurrentSheet solarwind_current_sheet_,
//! Magnetic topology region (0: nowhere, 1: same as HCS)
//#define SOLARWIND_SECTORED_REGION 0
SolarWindOptions::SectoredRegion solarwind_sectored_region_,
//! Correction to Parker Spiral, mainly for polar regions (0: none, 1: Smith-Bieber 1991, 2: Zurbuchen et al. 1997, 3: Schwadron-McComas 2003)
//#define SOLARWIND_POLAR_CORRECTION 0
SolarWindOptions::PolarCorrection solarwind_polar_correction_,
//! Latitudinal profile for bulk speed (0: constant, 1: linear step, 2: smooth step)
//#define SOLARWIND_SPEED_LATITUDE_PROFILE 0
SolarWindOptions::SpeedLatitudeProfile solarwind_speed_latitude_profile_,
//! Integer exponent of decrease of solar wind speed beyond the termination shock
//#define SOLARWIND_TERMSHOCK_SPEED_EXPONENT 2
SolarWindOptions::TermShockSpeedExponent solarwind_termshock_speed_exponent_,
//! What function to use within 'get_ampfactor' (0 = none, 1 = zero, 2 = constant, 3 = scaled)
//#define MOD_TYPE 3
VLISMBochumOptions::ModType VLISMBochum_mod_type_,
//! Whether to scale relative to s=0 (0) or s=+inf (1)
//#define MOD_RPOS 0
VLISMBochumOptions::ModRPos VLISMBochum_mod_rpos_,
// ----- Trajectory ----- //
////! Whether to record magnetic field extrema
bool record_mag_extrema_,
////! Whether to record segments
bool record_trajectory_,
////! Default initial size of trajectory segment record
size_t record_trajectory_segment_presize_,
//! Trajectory advance safety level: 0 means no checks, 1 means check dt only, 2 means check dt, number of segments, and time adaptations per step.
TrajectoryOptions::SafetyLevel trajectory_adv_safety_level_,
//! Largest length for single trajectory
int max_trajectory_steps_,
//! Largest number of time step adaptations for a single time step
int max_time_adaptations_,
//! Upper limit on the number of steps in debug mode
//! use -1 for unlimited
int n_max_calls_,
////! CFL condition for advection
std::ratio cfl_advection_,
////! CFL condition for diffusion
std::ratio cfl_diffusion_,
////! CFL condition for acceleration
std::ratio cfl_acceleration_,
////! CFL condition for pitch angle scattering
std::ratio cfl_pitchangle_,
//! Safety factor for drift-based time step (to modify "drift_vel" with a small fraction of the particle's velocity)
std::ratio drift_safety_,
//! How many time steps to allow before recording a mirror event
int mirror_threshold_,
//// Switch controlling how to calculate mu. "0" means computing it at the end of the step from magnetic moment conservation. "1" means advancing it in time according to the scheme (does not guarantee conservation of MM, but can be used with non-adiabatic terms).
TrajectoryOptions::PPerpMethod pperp_method_,
//! Flag to use gradient and curvature drifts in drift velocity calculation
TrajectoryOptions::UseBDrifts use_B_drifts_,
//! Which stochastic method to use for PA scattering, 0 = Euler, 1 = Milstein, 2 = RK2
////! Which stochastic method to use for perpendicular diffusion, 0 = Euler, 1 = Milstein, 2 = RK2
TrajectoryOptions::StochasticMethod stochastic_method_,
TrajectoryOptions::StochasticMethod stochastic_method_mu_,
TrajectoryOptions::StochasticMethod stochastic_method_perp_,
////! Whether to split the diffusive advance into two (one before and one after the advection).
//// #define SPLIT_SCATT
bool split_scatt_,
//#ifdef SPLIT_SCATT
////! Fraction of stochastic step to take before deterministic step
//const double alpha = 0.5;
//#endif
std::ratio split_scatt_fraction_,
//#if CONST_DMUMAX == 1
////! Desired accuracy in pitch angle cosine
//const double dmumax = 0.02;
//#else
////! Desired accuracy in pitch angle (deg x [deg to rad conversion factor])
//const double dthetamax = 2.0 * M_PI / 180.0;
//#endif
TrajectoryOptions::ConstDmumax const_dmumax_,
////! Number of time steps per one orbit
int steps_per_orbit_,
////! Which method of computation to use for divK: 0 = using direct central FD, 1 = using _spdata.grad quantities
TrajectoryOptions::DivkMethod divk_method_,
////! Maximum allowed fraction of momentum change per step
std::ratio dlnp_max_,
// ----- Diffusion ----- //
//! Flag to use QLT pitch angle scattering with WLNT perpendicular diffusion
bool use_qlt_scatt_,
//! diffusion coefficient (if constant)
//! Scattering frequency (persistent)
std::ratio D0_,
std::ratio Dperp_,
std::ratio Dpara_,
std::ratio A2A_,
std::ratio l_max_,
std::ratio k_min_,
std::ratio ps_index_,
std::ratio ps_minus_,
std::ratio A2T_,
std::ratio A2L_,
std::ratio ps_plus_,
std::ratio l_max_HP_,
std::ratio dl_max_,
std::ratio nose_z_nose_,
std::ratio nose_z_sheath_,
std::ratio nose_dz_,
std::ratio kappa0_,
std::ratio kappa_ratio_,
std::ratio U0_,
std::ratio p0_,
std::ratio T0_,
std::ratio r0_,
std::ratio lam0_,
std::ratio R0_,
std::ratio B0_,
std::ratio pow_law_U_,
std::ratio pow_law_p_,
std::ratio pow_law_T_,
std::ratio pow_law_r_,
std::ratio pow_law_R_,
std::ratio pow_law_B_,
int LISM_idx_,
std::ratio LISM_ind_,
std::ratio kappa_ratio_inner_,
std::ratio kappa_ratio_outer_,
std::ratio lam_inner_,
std::ratio lam_outer_,
std::ratio kappa_inner_,
std::ratio kappa_outer_,
std::ratio lam_para_,
std::ratio lam_perp_,
std::ratio kappa_ratio_red_,
std::ratio radial_limit_perp_red_,
std::ratio solar_cycle_idx_,
std::ratio solar_cycle_effect_,
int Bmix_idx_,
std::ratio Bmix_ind_
>
struct HConfig {

   using BackgroundCoordinates = BackgroundCoordinates_;
   using BackgroundFields = BackgroundFields_;
   using TrajectoryCoordinates = TrajectoryCoordinates_;
   using TrajectoryFields = TrajectoryFields_;
   using RecordCoordinates = RecordCoordinates_;
   using DiffusionCoordinates = DiffusionCoordinates_;
   using DiffusionFields = DiffusionFields_;

   static constexpr auto build_mode = build_mode_;
   static constexpr auto specie = Specie<specieid>();
   static constexpr auto derivative_method = derivative_method_;
   static constexpr auto time_flow = time_flow_;
   static constexpr auto num_numeric_grad_evals = num_numeric_grad_evals_;
   static constexpr auto incr_dmax_ratio = static_cast<double>(incr_dmax_ratio_);
   static constexpr auto server_interpolation_order = server_interpolation_order_;
   static constexpr auto smooth_discontinuity_order = smooth_discontinuity_order_;
   static constexpr auto server_num_ghost_cells = server_num_ghost_cells_;
   static constexpr auto solarwind_current_sheet = solarwind_current_sheet_;
   static constexpr auto solarwind_sectored_region = solarwind_sectored_region_;
   static constexpr auto solarwind_polar_correction = solarwind_polar_correction_;
   static constexpr auto solarwind_speed_latitude_profile = solarwind_speed_latitude_profile_;
   static constexpr auto solarwind_termshock_speed_exponent = solarwind_termshock_speed_exponent_;
   static constexpr auto VLISMBochum_mod_type = VLISMBochum_mod_type_;
   static constexpr auto VLISMBochum_mod_rpos = VLISMBochum_mod_rpos_;
   static constexpr auto record_mag_extrema = record_mag_extrema_;
   static constexpr auto record_trajectory = record_trajectory_;
   static constexpr auto record_trajectory_segment_presize = record_trajectory_segment_presize_;
   static constexpr auto trajectory_adv_safety_level = trajectory_adv_safety_level_;
   static constexpr auto max_trajectory_steps = max_trajectory_steps_;
   static constexpr auto max_time_adaptations = max_time_adaptations_;
   static constexpr auto n_max_calls = n_max_calls_;
   static constexpr auto cfl_advection = static_cast<double>(cfl_advection_);
   static constexpr auto cfl_diffusion = static_cast<double>(cfl_diffusion_);
   static constexpr auto cfl_acceleration = static_cast<double>(cfl_acceleration_);
   static constexpr auto cfl_pitchangle = static_cast<double>(cfl_pitchangle_);
   static constexpr auto drift_safety = static_cast<double>(drift_safety_);
   static constexpr auto mirror_threshold = static_cast<double>(mirror_threshold_);
   static constexpr auto pperp_method = pperp_method_;
   static constexpr auto use_B_drifts = use_B_drifts_;
   static constexpr auto stochastic_method = stochastic_method_;
   static constexpr auto stochastic_method_mu = stochastic_method_mu_;
   static constexpr auto stochastic_method_perp = stochastic_method_perp_;
   static constexpr auto split_scatt = split_scatt_;
   static constexpr auto split_scatt_fraction = static_cast<double>(split_scatt_fraction_);
   static constexpr auto const_dmumax = const_dmumax_;
   static constexpr auto steps_per_orbit = steps_per_orbit_;
   static constexpr auto divk_method = divk_method_;
   static constexpr auto dlnp_max = static_cast<double>(dlnp_max_);
   static constexpr auto use_qlt_scatt = static_cast<double>(use_qlt_scatt_);
   static constexpr auto D0 = static_cast<double>(D0_);
   static constexpr auto Dperp = static_cast<double>(Dperp_);
   static constexpr auto Dpara = static_cast<double>(Dpara_);
   static constexpr auto A2A = static_cast<double>(A2A_);
   static constexpr auto l_max = static_cast<double>(l_max_);
   static constexpr auto k_min = static_cast<double>(k_min_);
   static constexpr auto ps_index = static_cast<double>(ps_index_);
   static constexpr auto ps_minus = static_cast<double>(ps_minus_);
   static constexpr auto A2T = static_cast<double>(A2T_);
   static constexpr auto A2L = static_cast<double>(A2L_);
   static constexpr auto ps_plus = static_cast<double>(ps_plus_);
   static constexpr auto l_max_HP = static_cast<double>(l_max_HP_);
   static constexpr auto dl_max = static_cast<double>(dl_max_);
   static constexpr auto nose_z_nose = static_cast<double>(nose_z_nose_);
   static constexpr auto nose_z_sheath = static_cast<double>(nose_z_sheath_);
   static constexpr auto nose_dz = static_cast<double>(nose_dz_);
   static constexpr auto kappa0 = static_cast<double>(kappa0_);
   static constexpr auto kappa_ratio = static_cast<double>(kappa_ratio_);
   static constexpr auto U0 = static_cast<double>(U0_);
   static constexpr auto p0 = static_cast<double>(p0_);
   static constexpr auto T0 = static_cast<double>(T0_);
   static constexpr auto r0 = static_cast<double>(r0_);
   static constexpr auto lam0 = static_cast<double>(lam0_);
   static constexpr auto R0 = static_cast<double>(R0_);
   static constexpr auto B0 = static_cast<double>(B0_);
   static constexpr auto pow_law_U = static_cast<double>(pow_law_U_);
   static constexpr auto pow_law_p = static_cast<double>(pow_law_p_);
   static constexpr auto pow_law_T = static_cast<double>(pow_law_T_);
   static constexpr auto pow_law_r = static_cast<double>(pow_law_r_);
   static constexpr auto poaw_law_R = static_cast<double>(pow_law_R_);
   static constexpr auto pow_law_B = static_cast<double>(pow_law_B_);
   static constexpr auto LISM_idx = LISM_idx_;
   static constexpr auto LISM_ind = static_cast<double>(LISM_ind_);
   static constexpr auto kappa_ratio_inner = static_cast<double>(kappa_ratio_inner_);
   static constexpr auto kappa_ratio_outer = static_cast<double>(kappa_ratio_outer_);
   static constexpr auto lam_inner = static_cast<double>(lam_inner_);
   static constexpr auto lam_outer = static_cast<double>(lam_outer_);
   static constexpr auto kappa_inner = static_cast<double>(kappa_inner_);
   static constexpr auto kappa_outer = static_cast<double>(kappa_outer_);
   static constexpr auto lam_para = static_cast<double>(lam_para_);
   static constexpr auto lam_perp = static_cast<double>(lam_perp_);
   static constexpr auto kappa_ratio_red = static_cast<double>(kappa_ratio_red_);
   static constexpr auto radial_limit_perp_red = static_cast<double>(radial_limit_perp_red_);
   static constexpr auto solar_cycle_idx = static_cast<double>(solar_cycle_idx_);
   static constexpr auto solar_cycle_effect = static_cast<double>(solar_cycle_effect_);
   static constexpr auto Bmix_idx = Bmix_idx_;
   static constexpr auto Bmix_ind = static_cast<double>(Bmix_ind_);

   // ----- Derived constants ----- //

   // The factor multiplying "SpeciesCharges[]" is applied in order to marry particle and fluid scales. See "LarmorRadius()" and "CyclotronFrequency()" functions in physics.hh for reference.
   static constexpr auto q = charge_mass_particle * SpeciesCharges[static_cast<size_t>(specie)];


};

};

#endif
